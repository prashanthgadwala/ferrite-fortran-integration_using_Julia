using Ferrite

# Function to generate a simple hexahedral mesh with one element
function generate_simple_mesh_sta()
    mesh = generate_grid(SerendipityQuadraticHexahedron, (1, 1, 1))
    return mesh
end

# Generate the mesh
mesh = generate_simple_mesh_sta()

# Display basic information about the mesh
println("Generated mesh:")
println(mesh)

# Optionally, inspect nodes and elements
println("Nodes:")
for node in mesh.nodes
    println(node)
end

println("Elements:")
for elem in mesh.cells
    println(elem)
end