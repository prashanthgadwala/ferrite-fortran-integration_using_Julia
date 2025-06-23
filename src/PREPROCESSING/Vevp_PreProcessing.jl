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

    """println("==== UMAT INPUT CHECK ====")
    println("stran: ", stran)
    println("dstran: ", dstran)
    println("F (dfgrd1): ", dfgrd1)
    println("statev (first 10): ", statev[1:10])
    println("props (first 10): ", PROPS[1:10])
    println("dtime: ", dtime)
    println("noel: ", noel, " npt: ", npt)
    println("==========================")"""

    if any(isnan, stran) || any(isinf, stran)
        error("NaN or Inf in strain input to UMAT")
    end
    if any(isnan, dfgrd1) || any(isinf, dfgrd1)
        error("NaN or Inf in deformation gradient input to UMAT")
    end
    if any(isnan, statev) || any(isinf, statev)
        error("NaN or Inf in statev input to UMAT")
    end

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

function assemble_cell!(Ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    node_ids = getnodes(cell)  # global node indices for this cell
    X_nodes = [dh.grid.nodes[i].x for i in node_ids]  # reference coordinates as Vecs
    X_nodes_mat = hcat(X_nodes...)'  # (n_nodes × 3) matrix, each row is a node
    u_mat = reshape(ue, 3, :)'  # (n_nodes × 3) for 3D
    x_nodes = X_nodes_mat .+ u_mat  # (n_nodes × 3) current coordinates#


    for q_point in 1:getnquadpoints(cellvalues)
        # Compute total strain at this quadrature point
        ϵ_ = function_symmetric_gradient(cellvalues, q_point, ue)
        ϵ = tensor_to_voigt6(ϵ_)
        # Compute previous strain if you want to use strain increment

        dϵ = zeros(6)
      
        n_nodes = length(node_ids)
        qr = QuadratureRule{RefHexahedron}(2)
        interpolation = Lagrange{RefHexahedron, 1}()

        ξ = Ferrite.getpoints(qr)[q_point]  # Reference coordinates of this quadrature point
        dNdξ = zeros(8, 3)
        for i in 1:8
            dNdξ[i,:] =Ferrite.reference_shape_gradient(interpolation, ξ, i)
        end
        J_xξ = x_nodes' * dNdξ
        J_Xξ = X_nodes_mat' * dNdξ

  

        detJ = det(J_Xξ)
        F = J_xξ * inv(J_Xξ)
        # println("F = ", F)

        # Call UMAT-based stress/tangent state[q_point]
        σ, D, state[q_point] = compute_stress_tangent(ϵ, dϵ, state_old[q_point], PROPS, nprops, t, F)
        
        if any(isnan, D) || any(isinf, D)
            @show ue, ϵ, F, state_old[q_point], σ, D
            error("NaN or Inf detected in tangent matrix D at step $t, q_point $q_point")
        end
        if any(isnan, ue) || any(isinf, ue)
            error("NaN or Inf detected in element displacement ue")
        end
        if any(isnan, F) || any(isinf, F)
            error("NaN or Inf detected in deformation gradient F")
        end 

        w_q = Ferrite.getweights(cellvalues.qr)[q_point]
        dΩ = detJ * w_q
        # println("q_point = ", q_point, ", detJ = ", detJ, ", w_q = ", w_q, ", dΩ = ", dΩ)

        for i in 1:n_basefuncs
            δϵ_ = shape_symmetric_gradient(cellvalues, q_point, i)
            δϵ = tensor_to_voigt6(δϵ_)
            re[i] += dot(δϵ, σ) * dΩ
            #println("re[$i] = ", re[i])
            for j in 1:i
                Δϵ_ = shape_symmetric_gradient(cellvalues, q_point, j)
                Δϵ = tensor_to_voigt6(Δϵ_)
                Ke[i, j] += δϵ' * D * Δϵ * dΩ
                #println("Ke[$i, $j] = ", Ke[i, j])
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
        assemble_cell!(ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, dh)
        assemble!(assembler, eldofs, ke, re)
    end
    return K, r
end

