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
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], dofs)
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
    @show node_ids
    @show X_nodes_mat


    # node_ids = Ferrite.get_cell_node_ids(cell)  # or similar, depending on Ferrite version
    # X_nodes_mat   = hcat(node_ids...)' 
    # X_nodes_mat = reshape(X_nodes_mat, 3, :)'                  # (n_nodes × dim) reference coordinates
    # @show X_nodes_mat
    # @show size(X_nodes_mat)
    # @show ue
    # @show size(ue)
    # u_mat    = reshape(ue, 3, :)'            # (n_nodes × 3) for 3D
    # x_nodes  = X_nodes_mat .+ u_mat                     # (n_nodes × 3) current coordinates


    # # eldofs = celldofs(cell) gives the global DOF indices for the element
    # eldofs = celldofs(dh, cell)
    # # For a 3D vector field, every node has 3 DOFs, so:
    # node_ids = Int[]
    # for i in 1:3:length(eldofs)
    #     push!(node_ids, (eldofs[i]-1) ÷ 3 + 1)
    # end


    # node_ids = celldofs(cell)
    # # Get reference coordinates for the element's nodes
    # X_nodes_mat = hcat(Xnode[node_ids]...)'  # N_nodes x 3

    # for i in 1:10
    #     println(X_nodes_mat[i, :])
    # end


    # # node_ids = Ferrite.get_node_ids(cell.reference)  # or Ferrite.cellnodes(cell.ref)
    # # X_nodes_mat = hcat(Xnode[node_ids]...)'    # (27, 3)
    # X_nodes_mat = Ferrite.getcoordinates(cell)  # (27, 3) for Hex27
    # u_nodes = reshape(ue, 3, :)'               # (27, 3)
    # x_nodes = X_nodes_mat .+ u_nodes           # (27, 3)
    # @show X_nodes_mat
    # @show ue

    # u_nodes = reshape(ue, 3, :)'
    # x_nodes = X_nodes_mat .+ u_nodes

    # x_nodes = X_nodes_mat .+ ue

    # @show node_ids
    # @show X_nodes_mat
    # @show ue
    # @show x_nodes





    for q_point in 1:getnquadpoints(cellvalues)
        # Compute total strain at this quadrature point
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue)
        # Compute previous strain if you want to use strain increment
        dϵ = zeros(6)

        for i in 1:length(node_ids)
            grad = shape_gradient(cellvalues, q_point, i)
            println("q_point=$q_point, i=$i, grad=", grad, ", size=", size(grad))
        end

        @show size(shape_gradient(cellvalues, q_point, 1))
        n_basefuncs = getnbasefunctions(cellvalues)
        dNdξ = hcat([shape_gradient(cellvalues, q_point, i)[:,1] for i in 1:length(node_ids)]...)'
        @show size(dNdξ)
        J_xξ = x_nodes' * dNdξ
        @show J_xξ
        J_Xξ = X_nodes_mat' * dNdξ
        @show J_Xξ
        @show det(J_Xξ)
        F = J_xξ * inv(J_Xξ)
        @show size(dNdξ)
        @show size(X_nodes_mat)
        @show size(J_Xξ)     
        @show J_Xξ


        # Call UMAT-based stress/tangent
        σ, D, updated_state = compute_stress_tangent(collect(ϵ), dϵ, state_old[q_point], PROPS, nprops, t, F)
        updated_state[101:106] = σ
        state[q_point] = updated_state

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⋅ σ) * dΩ
            for j in 1:i
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
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