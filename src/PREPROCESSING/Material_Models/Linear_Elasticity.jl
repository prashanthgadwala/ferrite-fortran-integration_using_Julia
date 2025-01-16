using Ferrite, FerriteGmsh, SparseArrays

function define_material_model()
    Emod = 200e3 # Young's modulus [MPa]
    ν = 0.3      # Poisson's ratio [-]

    Gmod = Emod / (2(1 + ν))  # Shear modulus
    Kmod = Emod / (3(1 - 2ν)) # Bulk modulus

    C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2,2}))
    return C
end