export extract_vonmises, extract_sigma11, extract_disp_mag, postprocess, plot_traction_displacement

function extract_vonmises(states)
    n_cells = size(states, 2)
    mises_values = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ = state[101:106]
            mises_values[el] += vonMises(σ)
        end
        mises_values[el] /= length(cell_states)
    end
    return mises_values
end

function extract_sigma11(states)
    n_cells = size(states, 2)
    sigma_11 = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ = state[101:106]
            sigma_11[el] += σ[1]
        end
        sigma_11[el] /= length(cell_states)
    end
    return sigma_11
end

function extract_disp_mag(u, grid)
    disp_mag = zeros(getnnodes(grid))
    for i in 1:getnnodes(grid)
        u_node = u[(3i-2):(3i)]
        disp_mag[i] = norm(u_node)
    end
    return disp_mag
end

function extract_true_strain(states)
    n_cells = size(states, 2)
    strains = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            # Assuming state[94] is the axial logarithmic (true) strain (adjust index as needed)
            strains[el] += state[94]
        end
        strains[el] /= length(cell_states)
    end
    return strains
end

function plot_true_strain_vs_time(strain_history, time)
    using Plots
    plot(time, strain_history, xlabel="Time [s]", ylabel="True Strain", title="True Strain vs Time", legend=false)
end


function extract_avg_sigma11(states)
    n_cells = size(states, 2)
    sigma_11 = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ = state[101:106]
            sigma_11[el] += σ[1]
        end
        sigma_11[el] /= length(cell_states)
    end
    return mean(sigma_11)
end

function extract_avg_true_strain(states)
    n_cells = size(states, 2)
    strains = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            strains[el] += state[10]
        end
        strains[el] /= length(cell_states)
    end
    return mean(strains)
end

function plot_stress_vs_strain(stress_history, strain_history)
    using Plots
    plot(strain_history, stress_history, xlabel="True Strain", ylabel="Sigma_11 [Pa]", title="Stress vs Strain", legend=false)
end


function plot_traction_vs_displacement(u_max, traction_magnitude)
    using Plots
    plot(u_max, traction_magnitude, xlabel="Max Displacement [m]", ylabel="Traction [Pa]", title="Traction vs Max Displacement", legend=false)
end

function postprocess(grid, dh, states, props, u, filename="plasticity")
    mises_values = extract_vonmises(states)
    sigma_11 = extract_sigma11(states)
    disp_mag = extract_disp_mag(u, grid)

    # Write to VTK for ParaView/VTK visualization
    VTKGridFile(filename, dh) do vtk
        write_solution(vtk, dh, u)
        write_cell_data(vtk, mises_values, "von Mises [Pa]")
        write_cell_data(vtk, sigma_11, "Sigma_11 [Pa]")
        write_node_data(vtk, disp_mag, "Displacement magnitude [m]")
    end

    # Quick Julia plots
    plot(mises_values, title="Von Mises Stress per Element", xlabel="Element", ylabel="Von Mises [Pa]", legend=false)
    plot!(sigma_11, title="Sigma_11 per Element", xlabel="Element", ylabel="Sigma_11 [Pa]", legend=false)
    scatter(disp_mag, title="Displacement Magnitude at Nodes", xlabel="Node", ylabel="|u| [m]", legend=false)
end