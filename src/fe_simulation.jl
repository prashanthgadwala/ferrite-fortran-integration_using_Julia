module FESimulation

using Ferrite

export generate_simple_mesh

function generate_simple_mesh()
    # Generate a simple hexahedral mesh with one element
    mesh = generate_grid(QuadraticHexahedron, (1, 1, 1))
    return mesh
end

end
