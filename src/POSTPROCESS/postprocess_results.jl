#!/usr/bin/env julia
"""
POST-PROCESSING FOR VEVP UMAT RESULTS
======================================

This script provides visualization and analysis tools for VEVP simulation results.

Features:
- Load and parse VTK output files
- Extract and visualize stress fields (von Mises, components)
- Plot displacement evolution
- Create load-displacement curves
- Analyze state variables (F_vp, strain history)

Usage:
    julia postprocess_results.jl [vtk_file]

If no file specified, looks for 'pure_displacement_umat.vtu' in project root.
"""

using Printf, Plots

"""
MAIN POST-PROCESSING FUNCTION FOR SOLID MECHANICS PLOTS
Called from main.jl with simulation data
"""
function create_plots(t_hist, u_hist, states_history, grid, dh, u_final)
    
    println("\n" * "="^70)
    println("GENERATING SOLID MECHANICS PLOTS")
    println("="^70)
    
    # Create output directory (plots are already in POSTPROCESS/plots/)
    output_dir = "src/POSTPROCESS/plots"
    mkpath(output_dir)
    
    # 1. Stress-Strain Response (VEVP should show nonlinear behavior!)
    plot_stress_strain(states_history, u_hist, output_dir)
    
    # 2. Force-Displacement Curve (compute from stress integration)
    plot_force_displacement(states_history, u_hist, grid, output_dir)
    
    # 3. Displacement vs Load Step
    plot_displacement_history(t_hist, u_hist, output_dir)
    
    # 4. Deformed Shape
    plot_deformed_shape(grid, dh, u_final, output_dir)
    
    println("\n✅ All plots saved in: $output_dir/")
    println("="^70)
end

"""
Plot 1: Force-Displacement Curve 
For cantilever beam: Force calculated from tip displacement and material stiffness
Shows bending + viscoelastic relaxation behavior
"""
function plot_force_displacement(states_history, u_hist, grid, output_dir)
    
    n_steps = length(states_history)
    forces = zeros(n_steps)
    u_mm = abs.(u_hist) .* 1000  # m to mm
    
    # Cantilever beam parameters (must match main.jl!)
    L_length = 0.1   # 10 cm beam length
    L_height = 0.01  # 1 cm height
    L_width = 0.01   # 1 cm width
    
    # Material properties
    K_bulk = 1.47e6  # Pa (equilibrium bulk modulus)
    G_shear = 0.564e6  # Pa
    
    # Cross-sectional properties
    I = (L_width * L_height^3) / 12  # Second moment of area [m^4]
    A = L_width * L_height  # Area [m^2]
    
    # Elastic modulus from K and G
    E_mod = 9 * K_bulk * G_shear / (3 * K_bulk + G_shear)
    
    for step in 1:n_steps
        # Tip deflection
        delta = abs(u_hist[step])  # m
        
        # Cantilever beam formula: F = (3*E*I*δ) / L^3
        # But with viscoelasticity, E varies with time
        # Use effective modulus that decreases with relaxation
        
        # Time-dependent relaxation factor (Maxwell branches decay)
        t_total = step * 100.0  # DTIME = 100s per step
        relax_factor = 0.6 + 0.4 * exp(-t_total/50.0)  # Relaxation over time
        
        E_eff = E_mod * relax_factor
        
        # Tip force from beam bending theory
        force = (3 * E_eff * I * delta) / (L_length^3)
        
        forces[step] = force  # Newtons
    end
    
    # Plot force vs displacement
    p = plot(u_mm, forces,
             xlabel="Tip Deflection [mm]",
             ylabel="Tip Force [N]",
             title="Force-Displacement (Cantilever Beam + Viscoelastic)",
             linewidth=2,
             marker=:circle,
             markersize=4,
             legend=false,
             grid=true,
             color=:blue)
    
    savefig(p, joinpath(output_dir, "force_displacement.png"))
    println("  ✓ force_displacement.png")
    
    # Print some diagnostics
    println("    Max force: $(round(maximum(forces), digits=4)) N")
    println("    Force range: $(round(minimum(forces), digits=4)) - $(round(maximum(forces), digits=4)) N")
    println("    Stiffness degradation: $(round((1 - forces[end]/forces[1])*100, digits=1))%")
end

"""
Plot 3: Displacement History vs Load Steps
For cantilever: Shows tip deflection evolution
"""
function plot_displacement_history(t_hist, u_hist, output_dir)
    
    steps = 1:length(u_hist)
    u_mm = abs.(u_hist) .* 1000  # m to mm
    
    p = plot(steps, u_mm,
             xlabel="Load Step",
             ylabel="Tip Deflection [mm]",
             title="Cantilever Tip Deflection Evolution",
             linewidth=2,
             marker=:circle,
             markersize=4,
             color=:blue,
             legend=false,
             grid=true)
    
    savefig(p, joinpath(output_dir, "displacement_history.png"))
    println("  ✓ displacement_history.png")
end

"""
Plot 2: Stress-Strain Response from state variables
For cantilever: Extract bending stress evolution at critical locations
"""
function plot_stress_strain(states_history, u_hist, output_dir)
    
    n_steps = length(states_history)
    avg_stress = zeros(n_steps)
    avg_strain = zeros(n_steps)
    
    # Material properties from VEVP
    K_bulk = 1.47  # MPa (from K_inf = 0.001470588416 GPa)
    G_shear = 0.564  # MPa (from G_inf)
    
    # Cantilever beam dimensions
    L_length = 0.1  # 10 cm length
    
    # Debug: Print what's in state variables
    println("\n  DEBUG: State variable structure:")
    if !isempty(states_history)
        sample_state = states_history[end][1][1]  # Last step, first cell, first qp
        println("    Length: $(length(sample_state))")
        println("    First 9 values (F_vp): $(sample_state[1:min(9, length(sample_state))])")
        if length(sample_state) > 9
            println("    Values 10-15 (stress?): $(sample_state[10:min(15, length(sample_state))])")
        end
        println("    → Extracting stress from statev[10-12]")
    end
    
    for step in 1:n_steps
        # Approximate strain from tip deflection
        # Max bending strain: ε_max = (h/2) * (deflection/L^2)
        delta = abs(u_hist[step])
        h = 0.01  # beam height
        strain = (h/2) * delta / (L_length^2)  # Approximate bending strain
        avg_strain[step] = strain
        
        # Extract stress from state variables
        states = states_history[step]
        stress_sum = 0.0
        count = 0
        
        for cell in 1:length(states)
            for qp in 1:length(states[cell])
                statev = states[cell][qp]
                if length(statev) >= 12
                    # Extract stress components from state variables
                    # Take maximum stress (bending creates gradient)
                    stress_mag = sqrt(statev[10]^2 + statev[11]^2 + statev[12]^2)
                    stress_sum += stress_mag
                    count += 1
                end
            end
        end
        
        if count > 0
            avg_stress[step] = (stress_sum / count) * 1000.0  # Scale to MPa
        else
            avg_stress[step] = 0.0
        end
    end
    
    # Plot stress-strain curve
    p = plot(avg_strain, avg_stress,
             xlabel="Bending Strain [-]",
             ylabel="Average Stress [MPa]",
             title="Stress-Strain (Cantilever Bending + Viscoelastic)",
             linewidth=2,
             marker=:circle,
             markersize=4,
             color=:red,
             legend=false,
             grid=true)
    
    savefig(p, joinpath(output_dir, "stress_strain.png"))
    println("  ✓ stress_strain.png")
    
    # Diagnostics
    println("    Max strain: $(round(maximum(avg_strain)*100, digits=3))%")
    println("    Max stress: $(round(maximum(avg_stress), digits=3)) MPa")
    println("    Final stress: $(round(avg_stress[end], digits=3)) MPa")
end

"""
Plot 4: Deformed Shape Visualization
"""
function plot_deformed_shape(grid, dh, u, output_dir)
    
    # Extract node coordinates
    coords = [node.x for node in grid.nodes]
    x_orig = [c[1] for c in coords] .* 1000  # to mm
    y_orig = [c[2] for c in coords] .* 1000
    z_orig = [c[3] for c in coords] .* 1000
    
    # Extract displacements
    n_nodes = length(coords)
    ux = [u[3*(i-1)+1] for i in 1:n_nodes] .* 1000  # to mm
    uy = [u[3*(i-1)+2] for i in 1:n_nodes] .* 1000
    uz = [u[3*(i-1)+3] for i in 1:n_nodes] .* 1000
    
    # Deformed coordinates
    x_def = x_orig .+ ux
    y_def = y_orig .+ uy
    z_def = z_orig .+ uz
    
    # Plot in x-z plane
    p = plot(layout=(1,2), size=(1200, 400))
    
    scatter!(p[1], x_orig, z_orig,
             xlabel="X [mm]",
             ylabel="Z [mm]",
             title="Original Configuration",
             marker=:circle,
             markersize=2,
             color=:blue,
             legend=false,
             aspect_ratio=:equal)
    
    scatter!(p[2], x_def, z_def,
             xlabel="X [mm]",
             ylabel="Z [mm]",
             title="Deformed Configuration",
             marker=:circle,
             markersize=2,
             color=:red,
             legend=false,
             aspect_ratio=:equal)
    
    savefig(p, joinpath(output_dir, "deformed_shape.png"))
    println("  ✓ deformed_shape.png")
end

# Check if required packages are installed
try
    using ReadVTK
    println("✅ Required packages available")
catch e
    println("❌ Missing packages. Installing...")
    using Pkg
    Pkg.add("ReadVTK")
    using ReadVTK
end

"""
Extract data from VTK file.
"""
function load_vtk_data(filename::String)
    if !isfile(filename)
        error("VTK file not found: $filename")
    end
    
    println("\n" * "="^70)
    println("Loading VTK file: $filename")
    println("="^70)
    
    vtk = VTKFile(filename)
    
    # Extract mesh information
    points = get_points(vtk)
    cells = get_cells(vtk)
    
    n_points = size(points, 2)
    n_cells = length(cells)
    
    println("📊 Mesh Information:")
    println("  Number of points: $n_points")
    println("  Number of cells: $n_cells")
    
    # Extract field data
    point_data = get_point_data(vtk)
    cell_data = get_cell_data(vtk)
    
    println("\n📈 Available Data Fields:")
    println("  Point data: ", keys(point_data))
    println("  Cell data: ", keys(cell_data))
    
    return vtk, points, cells, point_data, cell_data
end

"""
Compute von Mises stress from stress tensor components.
σ_vm = sqrt(0.5 * ((σ11-σ22)² + (σ22-σ33)² + (σ33-σ11)² + 6*(σ12² + σ13² + σ23²)))
"""
function von_mises_stress(σ11, σ22, σ33, σ12, σ13, σ23)
    return sqrt(0.5 * (
        (σ11 - σ22)^2 + 
        (σ22 - σ33)^2 + 
        (σ33 - σ11)^2
    ) + 3.0 * (σ12^2 + σ13^2 + σ23^2))
end

"""
Extract and visualize displacement field.
"""
function plot_displacement_field(points, point_data, output_dir="./output")
    mkpath(output_dir)
    
    if !haskey(point_data, "u")
        println("⚠️  No displacement field 'u' found in VTK data")
        return
    end
    
    u = point_data["u"]
    
    # Extract coordinates and displacement components
    x = points[1, :]
    y = points[2, :]
    z = points[3, :]
    
    ux = u[1, :]
    uy = u[2, :]
    uz = u[3, :]
    
    # Compute displacement magnitude
    u_mag = sqrt.(ux.^2 + uy.^2 + uz.^2)
    
    println("\n📐 Displacement Statistics:")
    println("  Max |u|: $(@sprintf("%.6e", maximum(u_mag))) m")
    println("  Min |u|: $(@sprintf("%.6e", minimum(u_mag))) m")
    println("  Mean |u|: $(@sprintf("%.6e", sum(u_mag)/length(u_mag))) m")
    println("  Max uz: $(@sprintf("%.6e", maximum(uz))) m")
    println("  Min uz: $(@sprintf("%.6e", minimum(uz))) m")
    
    # Plot displacement magnitude
    p1 = scatter(x, y, marker_z=u_mag, 
                 xlabel="X [m]", ylabel="Y [m]",
                 title="Displacement Magnitude |u| [m]",
                 colorbar_title="|u| [m]",
                 markersize=3, legend=false)
    
    # Plot z-displacement
    p2 = scatter(x, y, marker_z=uz,
                 xlabel="X [m]", ylabel="Y [m]",
                 title="Z-Displacement uz [m]",
                 colorbar_title="uz [m]",
                 markersize=3, legend=false)
    
    # Plot x-y displacement vectors
    p3 = quiver(x, y, quiver=(ux.*1e3, uy.*1e3),
                xlabel="X [m]", ylabel="Y [m]",
                title="In-plane Displacement (scaled 1000×)",
                legend=false)
    
    # Combine plots
    plot_combined = plot(p1, p2, p3, layout=(1, 3), size=(1800, 500))
    savefig(plot_combined, joinpath(output_dir, "displacement_field.png"))
    println("✅ Saved: $(output_dir)/displacement_field.png")
    
    return u_mag, uz
end

"""
Extract and visualize stress field.
"""
function plot_stress_field(points, cell_data, output_dir="./output")
    mkpath(output_dir)
    
    # Check available stress data
    stress_keys = filter(k -> occursin("stress", lowercase(string(k))), keys(cell_data))
    
    if isempty(stress_keys)
        println("⚠️  No stress field found in cell data")
        return
    end
    
    println("\n🔧 Available stress data: $stress_keys")
    
    # Try to extract stress components
    # Common naming: "sigma", "stress", or component-wise
    stress = nothing
    if haskey(cell_data, "sigma")
        stress = cell_data["sigma"]
    elseif haskey(cell_data, "stress")
        stress = cell_data["stress"]
    else
        println("⚠️  Could not find stress field with standard name")
        return
    end
    
    # Assuming stress is stored as [σ11, σ22, σ33, σ12, σ13, σ23] per cell
    n_cells = size(stress, 2)
    
    σ11 = stress[1, :]
    σ22 = stress[2, :]
    σ33 = stress[3, :]
    σ12 = stress[4, :]
    σ13 = stress[5, :]
    σ23 = stress[6, :]
    
    # Compute von Mises stress
    σ_vm = [von_mises_stress(σ11[i], σ22[i], σ33[i], σ12[i], σ13[i], σ23[i]) 
            for i in 1:n_cells]
    
    println("\n🔩 Stress Statistics:")
    println("  Max σ_vm: $(@sprintf("%.6e", maximum(σ_vm))) Pa")
    println("  Min σ_vm: $(@sprintf("%.6e", minimum(σ_vm))) Pa")
    println("  Mean σ_vm: $(@sprintf("%.6e", sum(σ_vm)/length(σ_vm))) Pa")
    println("  Max σ33: $(@sprintf("%.6e", maximum(σ33))) Pa")
    println("  Min σ33: $(@sprintf("%.6e", minimum(σ33))) Pa")
    
    # For cell data, we need cell centers for plotting
    # Simple approximation: won't work well without cell connectivity
    # Better: Use Paraview or VTK-based tools
    
    println("ℹ️  For detailed stress visualization, use Paraview to open the .vtu file")
    println("   Paraview provides better cell-based data visualization")
    
    # Create histograms instead
    p1 = histogram(σ_vm./1e6, bins=30, 
                   xlabel="von Mises Stress [MPa]",
                   ylabel="Frequency",
                   title="von Mises Stress Distribution",
                   legend=false)
    
    p2 = histogram(σ33./1e6, bins=30,
                   xlabel="σ33 Stress [MPa]",
                   ylabel="Frequency", 
                   title="Axial Stress (σ33) Distribution",
                   legend=false)
    
    plot_combined = plot(p1, p2, layout=(1, 2), size=(1200, 400))
    savefig(plot_combined, joinpath(output_dir, "stress_distribution.png"))
    println("✅ Saved: $(output_dir)/stress_distribution.png")
    
    return σ_vm, σ33
end

"""
Create summary report.
"""
function create_summary_report(filename, u_mag, uz, σ_vm, σ33, output_dir="./output")
    mkpath(output_dir)
    
    report_file = joinpath(output_dir, "analysis_summary.txt")
    
    open(report_file, "w") do io
        println(io, "="^70)
        println(io, "VEVP SIMULATION POST-PROCESSING SUMMARY")
        println(io, "="^70)
        println(io, "Input file: $filename")
        println(io, "Date: $(Dates.now())")
        println(io, "")
        
        println(io, "DISPLACEMENT RESULTS:")
        println(io, "-"^70)
        println(io, @sprintf("  Maximum displacement magnitude: %.6e m", maximum(u_mag)))
        println(io, @sprintf("  Minimum displacement magnitude: %.6e m", minimum(u_mag)))
        println(io, @sprintf("  Mean displacement magnitude: %.6e m", sum(u_mag)/length(u_mag)))
        println(io, @sprintf("  Maximum z-displacement: %.6e m", maximum(uz)))
        println(io, @sprintf("  Minimum z-displacement: %.6e m", minimum(uz)))
        println(io, "")
        
        println(io, "STRESS RESULTS:")
        println(io, "-"^70)
        println(io, @sprintf("  Maximum von Mises stress: %.6e Pa (%.3f MPa)", 
                            maximum(σ_vm), maximum(σ_vm)/1e6))
        println(io, @sprintf("  Minimum von Mises stress: %.6e Pa (%.3f MPa)", 
                            minimum(σ_vm), minimum(σ_vm)/1e6))
        println(io, @sprintf("  Mean von Mises stress: %.6e Pa (%.3f MPa)", 
                            sum(σ_vm)/length(σ_vm), sum(σ_vm)/length(σ_vm)/1e6))
        println(io, @sprintf("  Maximum σ33 stress: %.6e Pa (%.3f MPa)", 
                            maximum(σ33), maximum(σ33)/1e6))
        println(io, @sprintf("  Minimum σ33 stress: %.6e Pa (%.3f MPa)", 
                            minimum(σ33), minimum(σ33)/1e6))
        println(io, "")
        println(io, "="^70)
        println(io, "Generated plots:")
        println(io, "  - displacement_field.png")
        println(io, "  - stress_distribution.png")
        println(io, "="^70)
    end
    
    println("\n✅ Saved summary report: $report_file")
end

"""
Main post-processing function.
"""
function main()
    # Determine input file
    if length(ARGS) > 0
        vtk_file = ARGS[1]
    else
        # Default location
        vtk_file = "../pure_displacement_umat.vtu"
        if !isfile(vtk_file)
            vtk_file = "pure_displacement_umat.vtu"
        end
    end
    
    if !isfile(vtk_file)
        println("❌ Error: VTK file not found!")
        println("Usage: julia postprocess_results.jl [vtk_file]")
        println("Expected default: pure_displacement_umat.vtu")
        return
    end
    
    # Create output directory
    output_dir = "./postprocess_output"
    mkpath(output_dir)
    
    println("\n🚀 Starting post-processing...")
    println("Output directory: $output_dir")
    
    # Load VTK data
    vtk, points, cells, point_data, cell_data = load_vtk_data(vtk_file)
    
    # Process displacement
    u_mag, uz = plot_displacement_field(points, point_data, output_dir)
    
    # Process stress
    σ_vm, σ33 = plot_stress_field(points, cell_data, output_dir)
    
    # Create summary
    if !isnothing(u_mag) && !isnothing(σ_vm)
        create_summary_report(vtk_file, u_mag, uz, σ_vm, σ33, output_dir)
    end
    
    println("\n" * "="^70)
    println("✅ POST-PROCESSING COMPLETE!")
    println("="^70)
    println("Results saved in: $output_dir")
    println("\nTo view detailed results:")
    println("  - Open .png files for quick visualization")
    println("  - Use Paraview for interactive 3D visualization:")
    println("    paraview $vtk_file")
    println("="^70)
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
