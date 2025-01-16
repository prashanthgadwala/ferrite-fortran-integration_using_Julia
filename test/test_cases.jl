using Test

# Include source files
include("../src/PREPROCESSING/Mesh_Files/logo_mesh.jl")
include("../src/PREPROCESSING/Material_Models/Linear_Elasticity.jl")
include("../src/PROCESSING/fe_simulation_LE.jl")
include("../src/POSTPROCESSING/Linear_Elastic.jl")

# Test mesh generation
@testset "Mesh Generation" begin
    grid = generate_mesh()
    @test length(grid.nodes) > 0
    @test getncells(grid) > 0
end

# Test material model definition
@testset "Material Model Definition" begin
    C = define_material_model()
    @test size(C) == (2, 2, 2, 2)  # Update the expected size based on the actual dimensions of C
end

# Test FE simulation
@testset "FE Simulation" begin
    grid = generate_mesh()
    C = define_material_model()
    dh, cellvalues, u, qr = run_simulation(grid, C)
    @test length(u) == ndofs(dh)
end

# Test postprocessing
@testset "Postprocessing" begin
    grid = generate_mesh()
    C = define_material_model()
    dh, cellvalues, u, qr = run_simulation(grid, C)
    @testset "Calculate Stresses" begin
        qp_stresses, avg_cell_stresses = calculate_stresses(grid, dh, cellvalues, u, C)
        @test length(qp_stresses) == getncells(grid)
        @test length(avg_cell_stresses) == 3
    end
    @testset "VTK Output" begin
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
        @test isfile(joinpath(output_dir, "linear_elasticity.vtu"))
    end
end