using Ferrite, Tensors, SparseArrays, LinearAlgebra, Printf

const qr = QuadratureRule{RefHexahedron}(2)
const facet_qr = FacetQuadratureRule{RefHexahedron}(3)
const interpolation = Lagrange{RefHexahedron, 1}()

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

    if any(isnan, stress) || any(isinf, stress) || any(isnan, ddsdde) || any(isinf, ddsdde)
        error("NaN or Inf in stress or tangent matrix returned from UMAT")
    end

    return stress, ddsdde, statev
end

function create_values(interpolation)
    cellvalues_u = CellValues(qr, interpolation)
    facetvalues_u = FacetValues(facet_qr, interpolation)
    return cellvalues_u, facetvalues_u
end

function create_dofhandler(grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

function create_bc(dh, grid)
    dbcs = ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end

function assemble_cell!(Ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    node_ids = getnodes(cell)
    X_nodes = [dh.grid.nodes[i].x for i in node_ids]
    X_nodes_mat = hcat(X_nodes...)'
    u_mat = reshape(ue, 3, :)'
    x_nodes = X_nodes_mat .+ u_mat

    for q_point in 1:getnquadpoints(cellvalues)
        ϵ_ = function_symmetric_gradient(cellvalues, q_point, ue)
        ϵ = tensor_to_voigt6(ϵ_)
        ϵ_old = tensor_to_voigt6(function_symmetric_gradient(cellvalues, q_point, ue .- ue)) # Placeholder for actual old strain
        dϵ = ϵ .- ϵ_old

        ξ = Ferrite.getpoints(qr)[q_point]
        dNdξ = zeros(8, 3)
        for i in 1:8
            dNdξ[i,:] = Ferrite.reference_shape_gradient(interpolation, ξ, i)
        end
        J_xξ = x_nodes' * dNdξ
        J_Xξ = X_nodes_mat' * dNdξ

        F = J_Xξ \ J_xξ

        σ, D, state[q_point] = compute_stress_tangent(ϵ, dϵ, state_old[q_point], PROPS, nprops, t, F)

        w_q = Ferrite.getweights(cellvalues.qr)[q_point]
        dΩ = det(J_Xξ) * w_q

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

function tensor_to_voigt(mat::AbstractMatrix)
    return [mat[1,1], mat[2,2], mat[3,3], mat[1,2], mat[1,3], mat[2,3], mat[2,1], mat[3,1], mat[3,2]]
end

function voigt_to_tensor(voit::AbstractVector)
    return [voit[1] voit[4] voit[5];
            voit[4] voit[2] voit[6];
            voit[5] voit[6] voit[3]]
end

function tensor_to_voigt6(mat::AbstractMatrix)
    return [mat[1,1], mat[2,2], mat[3,3], mat[1,2], mat[1,3], mat[2,3]]
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end

function doassemble!(K::SparseMatrixCSC, r::Vector, cellvalues::CellValues, dh::DofHandler,
                     PROPS::Vector{Float64}, u, states, states_old, nprops, t)
    assembler = start_assemble(K, r)
    nu = getnbasefunctions(cellvalues)
    re = zeros(nu)
    ke = zeros(nu, nu)

    for (i, cell) in enumerate(CellIterator(dh))
        fill!(ke, 0)
        fill!(re, 0)
        eldofs = celldofs(cell)
        ue = u[eldofs]

        state = @view states[:, i]
        state_old = @view states_old[:, i]

        assemble_cell!(ke, re, cell, cellvalues, PROPS, nprops, ue, state, state_old, t, dh)
        assemble!(assembler, eldofs, ke, re)
    end
    return K, r
end
