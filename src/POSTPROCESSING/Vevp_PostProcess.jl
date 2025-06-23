using Plots

export extract_avg_true_strain, extract_avg_sigma11, extract_true_strain, extract_sigma11, extract_vonmises, extract_disp_mag, plot_true_strain_vs_time, plot_stress_vs_strain, plot_traction_vs_displacement

function extract_true_strain(states)
    n_cells = size(states, 2)
    strains = zeros(n_cells)
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            strains[el] += state[10]  # E_ve[1,1]
        end
        strains[el] /= length(cell_states)
    end
    return strains
end

function extract_avg_true_strain(states)
    mean(extract_true_strain(states))
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

function extract_avg_sigma11(states)
    mean(extract_sigma11(states))
end

function vonMises(σ::AbstractVector{<:Real})
    return sqrt(
        0.5 * (
            (σ[1] - σ[2])^2 +
            (σ[2] - σ[3])^2 +
            (σ[3] - σ[1])^2
        ) + 3 * (σ[4]^2 + σ[5]^2 + σ[6]^2)
    )
end

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

function extract_disp_mag(u, grid)
    disp_mag = zeros(getnnodes(grid))
    for i in 1:getnnodes(grid)
        u_node = u[(3i-2):(3i)]
        disp_mag[i] = norm(u_node)
    end
    return disp_mag
end

function plot_true_strain_vs_time(strain_history, time)
    plot(time, strain_history, xlabel="Time [s]", ylabel="True Strain", title="True Strain vs Time", legend=false)
end

function plot_stress_vs_strain(stress_history, strain_history)
    plot(strain_history, stress_history, xlabel="True Strain", ylabel="Sigma_11 [Pa]", title="Stress vs Strain", legend=false)
end

function plot_traction_vs_displacement(u_max, traction_magnitude)
    plot(u_max, traction_magnitude, xlabel="Max Displacement [m]", ylabel="Traction [Pa]", title="Traction vs Max Displacement", legend=false)
end

