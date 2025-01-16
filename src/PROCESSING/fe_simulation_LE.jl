using Ferrite, FerriteGmsh, SparseArrays

function run_simulation(grid, C)
    dim = 2
    order = 1 # linear interpolation
    ip = Lagrange{RefTriangle, order}()^dim # vector valued interpolation

    qr = QuadratureRule{RefTriangle}(1) # 1 quadrature point
    qr_face = FacetQuadratureRule{RefTriangle}(1)

    cellvalues = CellValues(qr, ip)
    facetvalues = FacetValues(qr_face, ip)

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0, 2))
    add!(ch, Dirichlet(:u, getfacetset(grid, "left"),   (x, t) -> 0.0, 1))
    close!(ch)

    traction(x) = Vec(0.0, 20e3 * x[1])

    function assemble_external_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
        # Create a temporary array for the facet's local contributions to the external force vector
        fe_ext = zeros(getnbasefunctions(facetvalues))
        for facet in FacetIterator(dh, facetset)
            # Update the facetvalues to the correct facet number
            reinit!(facetvalues, facet)
            # Reset the temporary array for the next facet
            fill!(fe_ext, 0.0)
            # Access the cell's coordinates
            cell_coordinates = getcoordinates(facet)
            for qp in 1:getnquadpoints(facetvalues)
                # Calculate the global coordinate of the quadrature point.
                x = spatial_coordinate(facetvalues, qp, cell_coordinates)
                tₚ = prescribed_traction(x)
                # Get the integration weight for the current quadrature point.
                dΓ = getdetJdV(facetvalues, qp)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, qp, i)
                    fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
                end
            end
            # Add the local contributions to the correct indices in the global external force vector
            assemble!(f_ext, celldofs(facet), fe_ext)
        end
        return f_ext
    end

    function assemble_cell!(ke, cellvalues, C)
        for q_point in 1:getnquadpoints(cellvalues)
            # Get the integration weight for the quadrature point
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:getnbasefunctions(cellvalues)
                # Gradient of the test function
                ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
                for j in 1:getnbasefunctions(cellvalues)
                    # Symmetric gradient of the trial function
                    ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                    ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
                end
            end
        end
        return ke
    end

    function assemble_global!(K, dh, cellvalues, C)
        # Allocate the element stiffness matrix
        n_basefuncs = getnbasefunctions(cellvalues)
        ke = zeros(n_basefuncs, n_basefuncs)
        # Create an assembler
        assembler = start_assemble(K)
        # Loop over all cells
        for cell in CellIterator(dh)
            # Update the shape function gradients based on the cell coordinates
            reinit!(cellvalues, cell)
            # Reset the element stiffness matrix
            fill!(ke, 0.0)
            # Compute element contribution
            assemble_cell!(ke, cellvalues, C)
            # Assemble ke into K
            assemble!(assembler, celldofs(cell), ke)
        end
        return K
    end

    K = allocate_matrix(dh)
    assemble_global!(K, dh, cellvalues, C)

    f_ext = zeros(ndofs(dh))
    assemble_external_forces!(f_ext, dh, getfacetset(grid, "top"), facetvalues, traction)

    apply!(K, f_ext, ch)
    u = K \ f_ext

    return dh, cellvalues, u, qr
end