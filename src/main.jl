include("PREPROCESSING/Mesh_Files/logo_mesh.jl")
include("PREPROCESSING/Material_Models/Linear_Elasticity.jl")

include("PROCESSING/fe_simulation_LE.jl")

include("POSTPROCESSING/Linear_Elastic.jl")


function main()
    # Preprocessing
    grid = generate_mesh()
    C = define_material_model()

    # Processing
    dh, cellvalues, u, qr = run_simulation(grid, C)

    # Postprocessing
    postprocess(grid, dh, cellvalues, u, C, qr)
end

main()
