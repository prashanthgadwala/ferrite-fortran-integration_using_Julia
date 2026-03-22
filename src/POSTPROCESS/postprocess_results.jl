#!/usr/bin/env julia
"""
Post-processing utilities for the unit-cube VEVP benchmark.

This module generates publication figures and VTK exports from solver outputs.
It is typically called from `main.jl` after convergence.
"""

using Printf, Plots
using WriteVTK
using Ferrite

"""
    create_plots(t_hist, u_hist, tip_disp_hist, states_history, stress_history, grid, dh, u_final)

Run the post-processing pipeline for a converged simulation.

# Arguments
- t_hist: Time history vector [n_steps]
- u_hist: Full displacement vectors for each converged step [n_steps]
- tip_disp_hist: Tip displacement in loading direction per step [n_steps]
- states_history: State variables at each time step [n_steps][n_cells][n_qpoints]
- stress_history: Cauchy stress at each time step [n_steps][n_cells][n_qpoints]
- grid: Ferrite grid object (unit cube mesh)
- dh: DofHandler with displacement field
- u_final: Final displacement vector [ndofs]

# Returns
No explicit return value. Plot files and VTK outputs are written to disk.
"""
function create_plots(t_hist, u_hist, tip_disp_hist, states_history, stress_history, grid, dh, u_final)
    
    println("\n" * "="^70)
    println("GENERATING VEVP CUBE TENSION TEST VISUALIZATIONS")
    println("="^70)
    
    # Create output directory
    output_dir = "src/POSTPROCESS/plots"
    mkpath(output_dir)
    
    # Viscoelastic strain evolution (E_ve components)
    plot_force_displacement(states_history, grid, output_dir)
    
    # Viscoplastic deformation (F_vp) and plastic strain (gamma)
    plot_stress_strain(states_history, output_dir)
    
    # Prescribed displacement history
    plot_displacement_history(t_hist, tip_disp_hist, output_dir)
    
    # 3D deformed-shape visualization
    plot_deformed_shape(grid, dh, u_final, output_dir)
    
    # VTK time-series export for ParaView
    export_vtk_with_results(t_hist, u_hist, states_history, stress_history, grid, dh, output_dir)
    
    println("\n✅ All plots saved in: $output_dir/")
    println("="^70)
end

"""
    export_vtk_with_results(t_hist, u_hist, states_history, stress_history, grid, dh, output_dir)

Export per-step VTK files and a PVD collection for ParaView animation.

The export includes displacement, selected state-variable summaries, and stress
components averaged per cell.

# Returns
No explicit return value. Writes `.vtu` and `.pvd` files.
"""
function export_vtk_with_results(t_hist, u_hist, states_history, stress_history, grid, dh, output_dir)
    
    println("\n📊 Exporting VTK files for ParaView...")
    
    # Create VTK subdirectory
    vtk_dir = joinpath(output_dir, "vtk_timesteps")
    mkpath(vtk_dir)
    
    n_steps = length(states_history)
    n_cells = getncells(grid)
    
    # Export each time step
    pvd = paraview_collection("$vtk_dir/results")
    
    for step in 1:n_steps
        # Get displacement at this step
        u_step = u_hist[step]
        
        # Create VTK file for this timestep
        vtk_file = joinpath(vtk_dir, "step_$(lpad(step, 4, '0'))")
        
        # Use VTKGridFile (Ferrite v1.0 API)
        VTKGridFile(vtk_file, grid) do vtk
            # Write displacement field
            write_solution(vtk, dh, u_step)
            
            # Extract state variables
            # STATEV layout (108 variables):
            # 1-9:   F_vp (viscoplastic def gradient, Voigt)
            # 10-18: E_ve (viscoelastic strain, Voigt)
            # 19:    gamma (accumulated plastic strain)
            # 20-28: b (kinematic hardening, Voigt)
            # 29-100: AA (8 VE branches, 9 components each)
            # 101-108: BB (8 VE branch scalars)
            
            states = states_history[step]
            stresses = stress_history[step]
            
            # Preallocate arrays for cell data
            gamma_plastic = zeros(n_cells)
            F_vp_11 = zeros(n_cells)
            F_vp_22 = zeros(n_cells)
            F_vp_33 = zeros(n_cells)
            E_ve_11 = zeros(n_cells)
            E_ve_22 = zeros(n_cells)
            E_ve_33 = zeros(n_cells)
            E_ve_eqv = zeros(n_cells)
            sigma_11 = zeros(n_cells)
            sigma_22 = zeros(n_cells)
            sigma_33 = zeros(n_cells)
            sigma_12 = zeros(n_cells)
            sigma_13 = zeros(n_cells)
            sigma_23 = zeros(n_cells)
            sigma_vm = zeros(n_cells)
            
            # Extract from state variables
            for cell_idx in 1:n_cells
                cell_states = states[cell_idx]
                cell_stresses = stresses[cell_idx]
                n_qp = length(cell_states)
                
                # Average over quadrature points
                F_vp_avg = zeros(3)
                E_ve_avg = zeros(3)
                gamma_avg = 0.0
                sigma_avg = zeros(6)
                
                for qp in 1:n_qp
                    statev = cell_states[qp]
                    sigma_qp = cell_stresses[qp]
                    
                    # F_vp diagonal components (1-3 in Voigt)
                    F_vp_avg[1] += statev[1]  # F_vp_11
                    F_vp_avg[2] += statev[2]  # F_vp_22
                    F_vp_avg[3] += statev[3]  # F_vp_33
                    
                    # E_ve diagonal components (10-12 in Voigt)
                    E_ve_avg[1] += statev[10]  # E_ve_11
                    E_ve_avg[2] += statev[11]  # E_ve_22
                    E_ve_avg[3] += statev[12]  # E_ve_33
                    
                    # Accumulated plastic strain (19)
                    gamma_avg += statev[19]

                    # Cauchy stress components (Voigt order: 11,22,33,12,13,23)
                    sigma_avg .+= sigma_qp
                end
                
                F_vp_avg ./= n_qp
                E_ve_avg ./= n_qp
                gamma_avg /= n_qp
                sigma_avg ./= n_qp
                
                # Store in arrays
                F_vp_11[cell_idx] = F_vp_avg[1]
                F_vp_22[cell_idx] = F_vp_avg[2]
                F_vp_33[cell_idx] = F_vp_avg[3]
                
                E_ve_11[cell_idx] = E_ve_avg[1]
                E_ve_22[cell_idx] = E_ve_avg[2]
                E_ve_33[cell_idx] = E_ve_avg[3]
                
                # Equivalent viscoelastic strain
                E_ve_eqv[cell_idx] = sqrt(E_ve_avg[1]^2 + E_ve_avg[2]^2 + E_ve_avg[3]^2)
                
                gamma_plastic[cell_idx] = gamma_avg
                sigma_11[cell_idx] = sigma_avg[1]
                sigma_22[cell_idx] = sigma_avg[2]
                sigma_33[cell_idx] = sigma_avg[3]
                sigma_12[cell_idx] = sigma_avg[4]
                sigma_13[cell_idx] = sigma_avg[5]
                sigma_23[cell_idx] = sigma_avg[6]
                sigma_vm[cell_idx] = von_mises_stress(
                    sigma_avg[1], sigma_avg[2], sigma_avg[3], sigma_avg[4], sigma_avg[5], sigma_avg[6]
                )
            end
            
            # Write cell data to VTK using Ferrite v1.0 API
            write_cell_data(vtk, F_vp_11, "F_vp_11")
            write_cell_data(vtk, F_vp_22, "F_vp_22")
            write_cell_data(vtk, F_vp_33, "F_vp_33")
            write_cell_data(vtk, E_ve_11, "E_ve_11")
            write_cell_data(vtk, E_ve_22, "E_ve_22")
            write_cell_data(vtk, E_ve_33, "E_ve_33")
            write_cell_data(vtk, E_ve_eqv, "E_ve_equivalent")
            write_cell_data(vtk, gamma_plastic, "gamma_plastic")
            write_cell_data(vtk, sigma_11, "sigma_11")
            write_cell_data(vtk, sigma_22, "sigma_22")
            write_cell_data(vtk, sigma_33, "sigma_33")
            write_cell_data(vtk, sigma_12, "sigma_12")
            write_cell_data(vtk, sigma_13, "sigma_13")
            write_cell_data(vtk, sigma_23, "sigma_23")
            write_cell_data(vtk, sigma_vm, "sigma_von_mises")
            
            # Add time information
            pvd[t_hist[step]] = vtk
        end
    end
    
    # Save the PVD collection file
    vtk_save(pvd)
    
    println("  ✅ Exported $n_steps VTK files to: $vtk_dir/")
    println("  ✅ ParaView collection file: $vtk_dir/results.pvd")
    println("\n  📌 To visualize in ParaView:")
    println("     1. Open ParaView")
    println("     2. File → Open → Select 'results.pvd'")
    println("     3. Apply → View time series animation")
    println("     4. Color by: E_ve_equivalent, gamma_plastic, F_vp_11, sigma_11, sigma_von_mises, etc.")
    println("     5. Add filters: Warp By Vector (u) to see deformation")
end

"""
    plot_force_displacement(states_history, grid, output_dir)

Plot viscoelastic strain component evolution over load steps.
"""
function plot_force_displacement(states_history, grid, output_dir)
    
    n_steps = length(states_history)
    E_ve_11_hist = zeros(n_steps)
    E_ve_22_hist = zeros(n_steps)
    E_ve_33_hist = zeros(n_steps)
    E_ve_mag_hist = zeros(n_steps)
    
    # Extract viscoelastic strain from state variables
    for step in 1:n_steps
        states = states_history[step]
        
        # Average over all cells and quadrature points
        E_ve_sum = zeros(3)
        count = 0
        
        for cell in 1:length(states)
            for qp in 1:length(states[cell])
                statev = states[cell][qp]
                
                # E_ve is stored in statev[10:18] (Voigt notation)
                E_ve_sum[1] += statev[10]  # E_ve_11
                E_ve_sum[2] += statev[11]  # E_ve_22
                E_ve_sum[3] += statev[12]  # E_ve_33
                count += 1
            end
        end
        
        E_ve_avg = E_ve_sum / count
        E_ve_11_hist[step] = E_ve_avg[1]
        E_ve_22_hist[step] = E_ve_avg[2]
        E_ve_33_hist[step] = E_ve_avg[3]
        E_ve_mag_hist[step] = sqrt(E_ve_avg[1]^2 + E_ve_avg[2]^2 + E_ve_avg[3]^2)
    end
    
    steps = 1:n_steps
    
    # Plot viscoelastic strain components
    p = plot(steps, [E_ve_11_hist, E_ve_22_hist, E_ve_33_hist, E_ve_mag_hist],
             xlabel="Load Step",
             ylabel="Viscoelastic Strain E_ve",
             title="VEVP: Viscoelastic Strain Evolution (Nonlinear)",
             label=["E_ve_11" "E_ve_22" "E_ve_33" "||E_ve||"],
             linewidth=3,
             marker=[:circle :square :diamond :star],
             markersize=5,
             grid=true,
             minorgrid=true,
             legend=:topleft)
    
    savefig(p, joinpath(output_dir, "force_displacement.png"))
    println("  ✓ force_displacement.png (viscoelastic strain evolution)")
    
    # Print diagnostics
    println("    Max E_ve magnitude: $(round(maximum(E_ve_mag_hist), digits=6))") 
    println("    Final E_ve_33 (loading dir): $(round(E_ve_33_hist[end], digits=6))")
end

"""
    plot_displacement_history(t_hist, tip_disp_hist, output_dir)

Plot prescribed displacement history for the displacement-controlled test.
"""
function plot_displacement_history(t_hist, tip_disp_hist, output_dir)
    
    steps = 1:length(tip_disp_hist)
    u_mm = tip_disp_hist .* 1000  # m to mm (keep sign)
    
    p = plot(steps, u_mm,
             xlabel="Load Step",
             ylabel="Applied Displacement [mm]",
              title="Displacement Loading History (Cube Tension Test)",
             linewidth=3,
             marker=:circle,
             markersize=5,
             color=:darkblue,
             legend=false,
             grid=true,
             minorgrid=true)
    
    savefig(p, joinpath(output_dir, "displacement_history.png"))
    println("  ✓ displacement_history.png")
    
    println("    Total displacement: $(round(u_mm[end], digits=3)) mm")
    println("    Loading rate: $(round(u_mm[end]/length(steps), digits=4)) mm/step")
end

"""
    plot_stress_strain(states_history, output_dir)

Plot viscoplastic deformation-gradient trends and accumulated plastic strain.
"""
function plot_stress_strain(states_history, output_dir)
    
    n_steps = length(states_history)
    F_vp_11_hist = zeros(n_steps)
    F_vp_22_hist = zeros(n_steps)
    F_vp_33_hist = zeros(n_steps)
    gamma_hist = zeros(n_steps)
    
    # Extract F_vp and gamma from state variables
    for step in 1:n_steps
        states = states_history[step]
        
        # Average over all cells and quadrature points
        F_vp_sum = zeros(3)
        gamma_sum = 0.0
        count = 0
        
        for cell in 1:length(states)
            for qp in 1:length(states[cell])
                statev = states[cell][qp]
                
                # F_vp diagonal components (1-3 in Voigt)
                F_vp_sum[1] += statev[1]  # F_vp_11
                F_vp_sum[2] += statev[2]  # F_vp_22
                F_vp_sum[3] += statev[3]  # F_vp_33
                
                # Accumulated plastic strain
                gamma_sum += statev[19]
                count += 1
            end
        end
        
        F_vp_avg = F_vp_sum / count
        F_vp_11_hist[step] = F_vp_avg[1]
        F_vp_22_hist[step] = F_vp_avg[2]
        F_vp_33_hist[step] = F_vp_avg[3]
        gamma_hist[step] = gamma_sum / count
    end
    
    steps = 1:n_steps
    
    # Create subplot layout
    p = plot(layout=(2,1), size=(900, 800))
    
    # Plot 1: F_vp components
    plot!(p[1], steps, [F_vp_11_hist, F_vp_22_hist, F_vp_33_hist],
          xlabel="Load Step",
          ylabel="F_vp (Viscoplastic Def. Gradient)",
          title="VEVP: Viscoplastic Deformation Evolution",
          label=["F_vp_11" "F_vp_22" "F_vp_33"],
          linewidth=3,
          marker=[:circle :square :diamond],
          markersize=5,
          grid=true,
          minorgrid=true,
          legend=:best)
    
    # Plot 2: Accumulated plastic strain
    plot!(p[2], steps, gamma_hist,
          xlabel="Load Step",
          ylabel="γ (Accumulated Plastic Strain)",
          title="Plastic Strain Accumulation (Nonlinear Growth)",
          linewidth=3,
          marker=:circle,
          markersize=5,
          color=:red,
          legend=false,
          grid=true,
          minorgrid=true)
    
    savefig(p, joinpath(output_dir, "stress_strain.png"))
    println("  ✓ stress_strain.png (F_vp and γ evolution)")
    
    # Diagnostics
    println("    Final γ (plastic strain): $(round(gamma_hist[end], digits=6))")
    println("    Final F_vp_33: $(round(F_vp_33_hist[end], digits=6))")
    if F_vp_33_hist[1] != 0
        println("    Plastic deformation: $((F_vp_33_hist[end] - F_vp_33_hist[1])/F_vp_33_hist[1]*100)%")
    end
end

"""
    plot_deformed_shape(grid, dh, u, output_dir)

Create side-by-side 3D plots of original and deformed nodal coordinates.
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
    
    # Displacement magnitude for coloring
    u_mag = sqrt.(ux.^2 + uy.^2 + uz.^2)
    
    # Create 3D plot with both configurations
    p = plot(layout=(1,2), size=(1400, 600), 
             camera=(30, 30),
             dpi=150)
    
    # Original configuration
    scatter!(p[1], x_orig, y_orig, z_orig,
             xlabel="X [mm]",
             ylabel="Y [mm]",
             zlabel="Z [mm]",
             title="Original Cube",
             marker=:circle,
             markersize=8,
             color=:blue,
             legend=false,
             camera=(30, 30),
             grid=true)
    
    # Deformed configuration with displacement coloring
    scatter!(p[2], x_def, y_def, z_def,
             xlabel="X [mm]",
             ylabel="Y [mm]",
             zlabel="Z [mm]",
             title="Deformed Cube (colored by |u|)",
             marker=:circle,
             markersize=8,
             marker_z=u_mag,
             color=:viridis,
             colorbar=true,
             colorbar_title="|u| [mm]",
             legend=false,
             camera=(30, 30),
             grid=true)
    
    savefig(p, joinpath(output_dir, "deformed_shape.png"))
    println("  ✓ deformed_shape.png (3D visualization)")
    
    println("    Max displacement: $(round(maximum(u_mag), digits=4)) mm")
    println("    Z-displacement range: [$(round(minimum(uz), digits=4)), $(round(maximum(uz), digits=4))] mm")
end

# Ensure optional VTK reader dependency is available for script-mode analysis.
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
    load_vtk_data(filename)

Load mesh points, cells, point-data, and cell-data fields from a VTK file.
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
    plot_displacement_field(points, point_data; output_dir="./output")

Generate displacement-field diagnostics and figures from VTK point data.
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
    plot_stress_field(points, cell_data; output_dir="./output")

Generate stress-distribution diagnostics and figures from VTK cell data.
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
    create_summary_report(filename, u_mag, uz, σ_vm, σ33; output_dir="./output")

Write a plain-text summary report for displacement and stress statistics.
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
    main()

CLI entry point for standalone VTK-based post-processing.
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
