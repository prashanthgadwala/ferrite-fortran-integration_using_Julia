include("../src/fe_simulation.jl")

using .FESimulation, Test

# Test simple mesh generation
@testset "Mesh Generation" begin
    mesh = generate_simple_mesh()
    @test length(mesh.nodes) > 0
    @test length(mesh.elements) > 0
end
