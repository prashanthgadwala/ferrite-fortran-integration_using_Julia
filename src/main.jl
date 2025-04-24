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

    ccall((:UMAT, joinpath(@__DIR__, "PREPROCESSING", "Material_Models", "libumat.so")), Cvoid,
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
        8.0,          # Number of Maxwell branches
        200.0e9,      # Young's modulus (E)
        0.3,          # Poisson's ratio (ν)
        200e6,        # Initial yield stress (σ₀)
        1.0e9,        # Hardening modulus (H)
        0.01,         # Viscosity parameter
        0.5,          # Exponent for viscoplasticity
        # Add other material properties as required by UMAT
    ]
    nprops = length(PROPS)


    # Geometry and mesh
    L = 10.0 # beam length [m]
    w = 1.0  # beam width [m]
    h = 1.0  # beam height [m]
    n = 2
    nels = (10n, n, 2n) # Number of elements in each spatial direction
    P1 = Vec((0.0, 0.0, 0.0))  # Start point for geometry
    P2 = Vec((L, w, h))        # End point for geometry
    grid = generate_grid(Hexahedron, nels, P1, P2)
    interpolation = Lagrange{RefHexahedron, 2}()^3

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
    states = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]
    states_old = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]

    # Newton-Raphson loop
    n_timesteps = 10
    u_max = zeros(n_timesteps)
    traction_magnitude = 1.e7 * range(0.5, 1.0, length=n_timesteps)
    NEWTON_TOL = 1e-6

    for timestep in 1:n_timesteps
        t = timestep
        traction = Vec((0.0, 0.0, traction_magnitude[timestep]))
        update!(dbcs, t)
        apply!(u, dbcs)

        for newton_itr in 1:10
            # Assemble stiffness matrix and residual vector
            doassemble!(K, r, cellvalues, dh, PROPS, u, states, states_old)

            # Call Fortran UMAT for material computations
            for cell in 1:getncells(grid)
                for qp in 1:nqp
                    # Extract required variables for UMAT
                    stress = states[qp, cell].σ
                    statev = zeros(108) # Example: adjust based on UMAT requirements
                    ddsdde = zeros(6, 6)
                    sse, spd, scd, rpl, ddsddt, drplde, drpldt = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    stran = zeros(6)
                    dstran = zeros(6)
                    time = [0.0, t]
                    dtime = 1.0
                    temp, dtemp, predef, dpred = 0.0, 0.0, 0.0, 0.0
                    cmname = "MaterialName"
                    ndi, nshr, ntens, nstatv = 3, 3, 6, 108
                    coords = zeros(3)
                    drot = zeros(3, 3)
                    pnewdt = 0.0
                    celent = 1.0
                    dfgrd0 = zeros(3, 3)
                    dfgrd1 = zeros(3, 3)
                    noel, npt, layer, kspt, kstep, kinc = 1, 1, 1, 1, 1, timestep

                    # Call UMAT
                    raw_stress = collect(states[qp, cell].σ)
                    raw_ddsdde = collect(ddsdde)
                    raw_stran = collect(stran)
                    raw_dstran = collect(dstran)
                    call_umat(raw_stress, statev, raw_ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
                    raw_stran, raw_dstran, time, dtime, temp, dtemp, predef, dpred, "MaterialName",
                    ndi, nshr, ntens, nstatv, PROPS, nprops, coords, drot, pnewdt,
                    celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
                end
            end

            # Neumann boundary conditions
            doassemble_neumann!(r, dh, getfacetset(grid, "right"), facetvalues, traction)
            norm_r = norm(r[Ferrite.free_dofs(dbcs)])

            if norm_r < NEWTON_TOL
                break
            end

            apply_zero!(K, r, dbcs)
            Δu = Symmetric(K) \ r
            u -= Δu
        end

        states_old .= states
        u_max[timestep] = maximum(abs, u)
    end

    # Postprocessing
    postprocess(grid, dh, states, nothing, u)
    plot_traction_displacement(u_max, traction_magnitude)
end

solve()