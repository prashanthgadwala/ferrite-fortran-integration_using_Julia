using Ferrite

# Generate a simple hexahedral mesh with one element
mesh = generate_simple_mesh(QuadraticHexahedron, (1, 1, 1))

# Display basic information about the mesh
println("Generated mesh:")
println(mesh)

# Optionally, inspect nodes and elements
println("Nodes:")
for node in mesh.nodes
    println(node)
end

println("Elements:")
for elem in mesh.elements
    println(elem)
end