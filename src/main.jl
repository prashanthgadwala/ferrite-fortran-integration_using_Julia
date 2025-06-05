include("./PREPROCESSING/Vevp_PreProcessing.jl")
include("./POSTPROCESSING/Vevp_PostProcess.jl")

# Wrapper function to call the Fortran UMAT subroutine
function call_umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
                   stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname_str,
                   ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
                   celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

    
                   # Prepare the CHARACTER*80 buffer
    cmname = Vector{UInt8}(undef, 80) # Allocate 80 bytes
    fill!(cmname, 0)                 # Fill with null bytes
    copyto!(cmname, cmname_str)      # Copy the string into the buffer

    ccall((:umat_, joinpath(@__DIR__, "PREPROCESSING", "Material_Models", "libumat.so")), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{UInt8}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Int32},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32},
           Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
          stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
          stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
          ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
          celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
end

function solve()
    # Define material properties in PROPS array

        PROPS = [
        2,                        # 1: Power Series Strain Approximation (dimensionless)
        1470.588416,            # 2: Bulk Modulus [Pa]
        563.9098439,            # 3: Shear Modulus [Pa]
        5.900948586,              # 4: Yield Exponent (dimensionless)
        0.33,                     # 5: Plastic Poisson Ratio (dimensionless)
        0.001,                    # 6: Viscoplastic Coefficient (dimensionless)
        10,                       # 7: Viscoplastic Exponent (dimensionless)
        2.086229688,            # 8: Initial Yield Limit - Compression [Pa]
        2164.115496,            # 9: Isotropic Hardening Parameter - Compression [Pa]
        4.450073598,            # 10: Isotropic Hardening Parameter - Compression [Pa]
        5401.554318,            # 11: Isotropic Hardening Parameter - Compression [Pa]
        1.66898375,             # 12: Initial Yield Limit - Tension [Pa]
        1731.292397,            # 13: Isotropic Hardening Parameter – Tension [Pa]
        3.560058879,            # 14: Isotropic Hardening Parameter – Tension [Pa]
        5401.554318,            # 15: Isotropic Hardening Parameter – Tension [Pa]
        0,                        # 16: Kinematic Hardening Parameter (dimensionless)
        0,                        # 17: Kinematic Hardening Parameter (dimensionless)
        0,                        # 18: Kinematic Hardening Parameter (dimensionless)
        298.5381412,            # 19: Bulk Modulus – Maxwell branch 1 [Pa]
        128.572959,             # 20: Bulk Modulus – Maxwell branch 2 [Pa]
        126.0116202,            # 21: Bulk Modulus – Maxwell branch 3 [Pa]
        56.16237131,            # 22: Bulk Modulus – Maxwell branch 4 [Pa]
        73.29632907,            # 23: Bulk Modulus – Maxwell branch 5 [Pa]
        34.58797299,            # 24: Bulk Modulus – Maxwell branch 6 [Pa]
        25.90272535,            # 25: Bulk Modulus – Maxwell branch 7 [Pa]
        117.6700904,            # 26: Bulk Modulus – Maxwell branch 8 [Pa]
        1000.0,                   # 27: Volumetric Relaxation Time – Maxwell Branch 1 [s]
        100.0,                    # 28: Volumetric Relaxation Time – Maxwell Branch 2 [s]
        10.0,                     # 29: Volumetric Relaxation Time – Maxwell Branch 3 [s]
        1.0,                      # 30: Volumetric Relaxation Time – Maxwell Branch 4 [s]
        0.1,                      # 31: Volumetric Relaxation Time – Maxwell Branch 5 [s]
        0.01,                     # 32: Volumetric Relaxation Time – Maxwell Branch 6 [s]
        0.001,                    # 33: Volumetric Relaxation Time – Maxwell Branch 7 [s]
        0.0001,                   # 34: Volumetric Relaxation Time – Maxwell Branch 8 [s]
        114.4770316,            # 35: Shear Modulus – Maxwell branch 1 [Pa]
        49.30241285,            # 36: Shear Modulus – Maxwell branch 2 [Pa]
        48.32024532,            # 37: Shear Modulus – Maxwell branch 3 [Pa]
        21.53594689,            # 38: Shear Modulus – Maxwell branch 4 [Pa]
        28.10611115,            # 39: Shear Modulus – Maxwell branch 5 [Pa]
        13.26305731,            # 40: Shear Modulus – Maxwell branch 6 [Pa]
        9.932624005,            # 41: Shear Modulus – Maxwell branch 7 [Pa]
        45.1216136,             # 42: Shear Modulus – Maxwell branch 8 [Pa]
        1000.0,                   # 43: Deviatoric Relaxation Time – Maxwell Branch 1 [s]
        100.0,                    # 44: Deviatoric Relaxation Time – Maxwell Branch 2 [s]
        10.0,                     # 45: Deviatoric Relaxation Time – Maxwell Branch 3 [s]
        1.0,                      # 46: Deviatoric Relaxation Time – Maxwell Branch 4 [s]
        0.1,                      # 47: Deviatoric Relaxation Time – Maxwell Branch 5 [s]
        0.01,                     # 48: Deviatoric Relaxation Time – Maxwell Branch 6 [s]
        0.001,                    # 49: Deviatoric Relaxation Time – Maxwell Branch 7 [s]
        0.0001                    # 50: Deviatoric Relaxation Time – Maxwell Branch 8 [s]
    ]
    
    nprops = length(PROPS) # Number of material properties


    # Geometry and mesh
    L = 10.0 # beam length [m]
    w = 1.0  # beam width [m]
    h = 1.0  # beam height [m]
    n = 4
    nels = (10*n, n, 2*n) # Number of elements in each spatial direction
    P1 = Vec((0.0, 0.0, 0.0))  # Start point for geometry
    P2 = Vec((L, w, h))        # End point for geometry
    grid = generate_grid(Hexahedron, nels, P1, P2)
    interpolation = Lagrange{RefHexahedron, 1}()^3
    println("Number of elements: ", getncells(grid))
    println("Number of nodes: ", getnnodes(grid))
    println("Number of right facets: ", length(getfacetset(grid, "right")))

    Xnode = [Ferrite.get_node_coordinate(grid, i) for i in 1:getnnodes(grid)]

    # Preprocessing
    dh = create_dofhandler(grid, interpolation)
    dbcs = create_bc(dh, grid)
    cellvalues, facetvalues = create_values(interpolation)

    # Pre-allocate solution vectors
    n_dofs = ndofs(dh)
    u = zeros(n_dofs)
    Δu = zeros(n_dofs)
    r = zeros(n_dofs)
    K = allocate_matrix(dh)

    # Material states
    nqp = getnquadpoints(cellvalues)
    nstatv = 108 # Number of state variables
    states = [zeros(nstatv) for _ in 1:nqp, _ in 1:getncells(grid)]
    states_old = [zeros(nstatv) for _ in 1:nqp, _ in 1:getncells(grid)]
    #println("First 10 statev after Newton step: ", states[1:min(10,end), 1])
    #println("First 10 statev after Newton step: ", states_old[1:min(10,end), 1])


    # Newton-Raphson loop
    n_timesteps = 10
    u_max = zeros(n_timesteps)
    traction_magnitude = 1.e-3 * range(0.1, 1.0, length=n_timesteps)
    NEWTON_TOL = 1e-6

    for timestep in 1:n_timesteps
        t = timestep
        traction = Vec((0.0, 0.0, traction_magnitude[timestep]))
        newton_itr = -1
        print("\n Time step @time = $timestep:\n")
        update!(dbcs, t)
        apply!(u, dbcs)
    
        while true; newton_itr += 1
            doassemble!(K, r, cellvalues, dh, PROPS, u, states, states_old, nprops, t)
            doassemble_neumann!(r, dh, getfacetset(grid, "right"), facetvalues, traction)
            #println("Assembled residual vector: ", r)
            #println("Assembled stiffness matrix: ", K)
            norm_r = norm(r[Ferrite.free_dofs(dbcs)])
            println("r before apply_zero! at constrained DOFs: ", count(iszero, r))
            apply_zero!(K, r, dbcs)
            println("r after apply_zero! at constrained DOFs: ", count(iszero, r))
            #println("Assembled stiffness matrix after applying BC: ", K)


            print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
            if norm_r < NEWTON_TOL
                break
            end
            Δu = Symmetric(K) \ r
            apply_zero!(Δu, dbcs)
            # if maximum(abs, Δu) > 0.1   # or some reasonable threshold
            #     println("Warning: Large displacement increment detected, scaling down Δu")
            #     Δu .*= 1.e-10 / maximum(abs, Δu)
            # end
            u -= Δu
            println("Max displacement: ", maximum(abs, u))
            println("Residual norm: ", norm_r)
        end

        for j in eachindex(states)
            for i in eachindex(states[j])
                states_old[j][i] = states[j][i]
            end
        end

        # states_old .= states

        # println("First 10 statev after Newton step: ", states[1:min(10,end), 1])
        
        u_max[timestep] = maximum(abs, u)
    end

    # Postprocessing
    # postprocess(grid, dh, states, PROPS, u)
    # plot_traction_displacement(u_max, traction_magnitude)

end

solve()