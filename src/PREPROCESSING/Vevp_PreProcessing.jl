using Ferrite, Tensors, SparseArrays, LinearAlgebra, Printf

function doassemble_neumann!(r, dh, facetset, facetvalues, t)
    n_basefuncs = getnbasefunctions(facetvalues)
    re = zeros(n_basefuncs)
    for fc in FacetIterator(dh, facetset)
        reinit!(facetvalues, fc)
        fill!(re, 0)
        for q_point in 1:getnquadpoints(facetvalues)
            dΓ = getdetJdV(facetvalues, q_point)
            for i in 1:n_basefuncs
                δu = shape_value(facetvalues, q_point, i)
                re[i] -= (δu ⋅ t) * dΓ
            end
        end
        assemble!(r, celldofs(fc), re)
    end
    return r
end

function compute_stress_tangent(ϵ::Vector{Float64}, dϵ::Vector{Float64}, statev::Vector{Float64}, PROPS::Vector{Float64}, nprops::Int64, t::Int64)
    stress = zeros(6)
    ddsdde = zeros(6, 6)
    sse, spd, scd, rpl, ddsddt, drplde, drpldt = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    stran = ϵ
    dstran = dϵ
    time = [0.0, t]
    dtime = 1.0
    temp, dtemp, predef, dpred = 0.0, 0.0, 0.0, 0.0
    cmname = "VEVP"
    ndi, nshr, ntens, nstatv = 3, 3, 6, length(statev)
    coords = zeros(3)
    drot = zeros(3, 3)
    pnewdt = 0.0
    celent = 1.0
    dfgrd0 = Matrix{Float64}(I, 3, 3)
    dfgrd1 = Matrix{Float64}(I, 3, 3)
    noel, npt, layer, kspt, kstep, kinc = 1, 1, 1, 1, 1, 1

    call_umat(
        stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
        stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
        ndi, nshr, ntens, nstatv, PROPS, nprops, coords, drot, pnewdt,
        celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc
    )

    return stress, ddsdde, statev
end

function create_values(interpolation)
    # Setup quadrature rules
    qr = QuadratureRule{RefHexahedron}(2)
    facet_qr = FacetQuadratureRule{RefHexahedron}(3)

    # Cell and facet values for u
    cellvalues_u = CellValues(qr, interpolation)
    facetvalues_u = FacetValues(facet_qr, interpolation)

    return cellvalues_u, facetvalues_u
end

function create_dofhandler(grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation) # Add a displacement field with 3 components
    close!(dh)
    return dh
end

function create_bc(dh, grid)
    dbcs = ConstraintHandler(dh)
    # Clamped on the left side
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], dofs)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end


function mat2voit(tensor::SymmetricTensor{2, 3, Float64, 6})
    return [
        tensor[1, 1],  # xx
        tensor[2, 2],  # yy
        tensor[3, 3],  # zz
        tensor[1, 2],  # xy
        tensor[1, 3],  # xz
        tensor[2, 3]   # yz
    ]
end

function assemble_cell!(Ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    for q_point in 1:getnquadpoints(cellvalues)
        # Compute total strain at this quadrature point
        ϵ_ = function_symmetric_gradient(cellvalues, q_point, ue)
        ϵ = mat2voit(ϵ_)

        # Compute previous strain if you want to use strain increment
        ϵ_old = state_old[q_point][1:6]  # If you store previous strain, else zeros(6)
        dϵ = ϵ - ϵ_old

        # Call UMAT-based stress/tangent
        σ, D, state[q_point] = compute_stress_tangent(ϵ, dϵ, state_old[q_point], PROPS, nprops, t)
        state[q_point][1:6] .= σ # Store current stress in state

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ_ = shape_symmetric_gradient(cellvalues, q_point, i)
            δϵ = mat2voit(δϵ_)
            re[i] += (δϵ ⋅ σ) * dΩ
            for j in 1:i
                Δϵ_ = shape_symmetric_gradient(cellvalues, q_point, j)
                Δϵ = mat2voit(Δϵ_)
                Ke[i, j] += δϵ' * D * Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function doassemble!(K::SparseMatrixCSC, r::Vector, cellvalues::CellValues, dh::DofHandler,
                     PROPS::Vector{Float64}, u, states, states_old, nprops, t)
    assembler = start_assemble(K, r)
    nu = getnbasefunctions(cellvalues)
    re = zeros(nu)     # element residual vector
    ke = zeros(nu, nu) # element tangent matrix

    for (i, cell) in enumerate(CellIterator(dh))
        fill!(ke, 0)
        fill!(re, 0)
        eldofs = celldofs(cell)
        ue = u[eldofs]

        state = @view states[:, i]
        state_old = @view states_old[:, i]
        # Call UMAT-based assembly
        assemble_cell!(ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t)
        assemble!(assembler, eldofs, ke, re)
    end
    return K, r
end