using Ferrite, FerriteGmsh, SparseArrays

function postprocess(grid, dh, cellvalues, u, C, qr)
    function calculate_stresses(grid, dh, cv, u, C)
        qp_stresses = [
            [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
            for _ in 1:getncells(grid)]
        avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:3)...)
        for cell in CellIterator(dh)
            reinit!(cv, cell)
            cell_stresses = qp_stresses[cellid(cell)]
            for q_point in 1:getnquadpoints(cv)
                ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
                cell_stresses[q_point] = C ⊡ ε
            end
            σ_avg = sum(cell_stresses) / getnquadpoints(cv)
            avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
            avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
            avg_cell_stresses[3][cellid(cell)] = σ_avg[1, 2]
        end
        return qp_stresses, avg_cell_stresses
    end

    qp_stresses, avg_cell_stresses = calculate_stresses(grid, dh, cellvalues, u, C)

    proj = L2Projector(Lagrange{RefTriangle, 1}(), grid)
    stress_field = project(proj, qp_stresses, qr)

    output_dir = joinpath(homedir(), "ferrite_output")
    mkpath(output_dir) 

    VTKGridFile(joinpath(output_dir, "linear_elasticity"), dh) do vtk
        write_solution(vtk, dh, u)
        for (i, key) in enumerate(("11", "22", "12"))
            write_cell_data(vtk, avg_cell_stresses[i], "sigma_" * key)
        end
        write_projection(vtk, proj, stress_field, "stress field")
        Ferrite.write_cellset(vtk, grid)
    end
end