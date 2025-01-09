# src/main.jl

include("fe_simulation.jl")

using .FESimulation

function main()
    # Generate the mesh
    mesh = generate_simple_mesh()

    # Display the generated mesh
    println("Generated mesh:")
    println(mesh)

    # Optionally inspect nodes and elements
    println("Nodes:")
    for node in mesh.nodes
        println(node)
    end

    println("Elements:")
    for elem in mesh.elements
        println(elem)
    end
end

main()
