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

function compute_stress_tangent(ϵ::Vector{Float64}, dϵ::Vector{Float64}, statev::Vector{Float64}, PROPS::Vector{Float64}, nprops::Int, t::Int64, F::Matrix{Float64})
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
    dfgrd1 = F
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
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end

function assemble_cell!(Ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, Xnode, dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    node_ids = getnodes(cell)  # global node indices for this cell

    X_nodes = [dh.grid.nodes[i].x for i in node_ids]  # reference coordinates as Vecs
    X_nodes_mat = hcat(X_nodes...)'  # (n_nodes × 3) matrix, each row is a node
    u_mat = reshape(ue, 3, :)'  # (n_nodes × 3) for 3D
    x_nodes = X_nodes_mat .+ u_mat  # (n_nodes × 3) current coordinates#ä


    for q_point in 1:getnquadpoints(cellvalues)
        # Compute total strain at this quadrature point
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue)
        # Compute previous strain if you want to use strain increment
        dϵ = zeros(6)

        dNdξ = hcat([shape_gradient(cellvalues, q_point, i)[:,1] for i in 1:length(node_ids)]...)'

        J_xξ = x_nodes' * dNdξ

        J_Xξ = X_nodes_mat' * dNdξ
        detJ = det(J_Xξ)
        F = J_xξ * inv(J_Xξ)


        # Call UMAT-based stress/tangent
        σ, D, state[q_point] = compute_stress_tangent(vec(collect(ϵ)), dϵ, state_old[q_point], PROPS, nprops, t, F)
        
        #if t == 1 && q_point == 1
        #    println("Tangent matrix D at first step:\n", D)
        #end

        if any(isnan, D) || any(isinf, D)
            error("NaN or Inf detected in tangent matrix D at step $t, q_point $q_point")
        end

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ_ = shape_symmetric_gradient(cellvalues, q_point, i)
            δϵ = tensor_to_voigt6(δϵ_)
            re[i] += dot(δϵ, σ) * dΩ
            for j in 1:i
                Δϵ_ = shape_symmetric_gradient(cellvalues, q_point, j)
                Δϵ = tensor_to_voigt6(Δϵ_)
                Ke[i, j] += δϵ' * D * Δϵ * dΩ
            end
        end


    end
    symmetrize_lower!(Ke)
end


# 3x3 matrix to 9-component Voigt vector (Fortran order)
function tensor_to_voigt(mat::AbstractMatrix)
    return [
        mat[1,1], mat[2,2], mat[3,3],
        mat[1,2], mat[1,3], mat[2,3],
        mat[2,1], mat[3,1], mat[3,2]
    ]
end

# 9-component Voigt vector to 3x3 matrix (Fortran order)
function voigt_to_tensor(voit::AbstractVector)
    return [
        voit[1] voit[4] voit[5];
        voit[4] voit[2] voit[6];
        voit[5] voit[6] voit[3]
    ]
end

function tensor_to_voigt6(mat::AbstractMatrix)
    # For symmetric 3x3 tensor: [11, 22, 33, 12, 13, 23]
    return [mat[1,1], mat[2,2], mat[3,3], mat[1,2], mat[1,3], mat[2,3]]
end


function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function doassemble!(K::SparseMatrixCSC, r::Vector, cellvalues::CellValues, dh::DofHandler,
                     PROPS::Vector{Float64}, u, states, states_old, nprops, t, Xnode)
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
        assemble_cell!(ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, Xnode, dh)
        assemble!(assembler, eldofs, ke, re)
    end
    return K, r
end

"""
function debug_compare_models(PROPS)
    ϵ = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]
    dϵ = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]
    statev = zeros(108)
    F = Matrix{Float64}(I, 3, 3)
    F[1,1] = 1.01
    t = 1
    nprops = length(PROPS)

    σ_umat, D_umat, statev_new = compute_stress_tangent(ϵ, dϵ, statev, PROPS, nprops, t, F)
    println("UMAT stress: ", σ_umat)
    println("UMAT tangent: ", D_umat)
    println("UMAT statev (first 10): ", statev_new[1:10])
end

"""