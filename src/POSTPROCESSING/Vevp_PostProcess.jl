using Ferrite, Plots, Tensors

function deviatoric(σ::SymmetricTensor{2, 3})
    I = one(σ) # Identity tensor
    return σ - (1 / 3) * tr(σ) * I
end

function vonMises(σ::SymmetricTensor{2, 3})
    s = deviatoric(σ) # Deviatoric part of the stress tensor
    result = sqrt(3 / 2 * tr(s ⊡ s)) # von Mises stress
    return result
end

function postprocess(grid, dh, states, states_old, props, u, filename)
    mises_values = zeros(getncells(grid))
    κ_values = zeros(getncells(grid))

    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ_voigt = state[1:6]           # Stress (now stored during assembly)
            σ = voigt_to_tensor(σ_voigt)
            mises_values[el] += vonMises(σ)

            κ_values[el] += state[19]      # Hardening variable (gma_n)
        end
        mises_values[el] /= length(cell_states)
        κ_values[el] /= length(cell_states)
    end

    VTKGridFile(filename, dh) do vtk
        write_solution(vtk, dh, u)
        write_cell_data(vtk, mises_values, "von Mises [Pa]")
        write_cell_data(vtk, κ_values, "Hardening Variable")
    end

    println("Postprocessing completed.")
end

function voigt_to_tensor(voigt::Vector{Float64})
    result = SymmetricTensor{2, 3}([
        voigt[1], voigt[4], voigt[5],
        voigt[4], voigt[2], voigt[6],
        voigt[5], voigt[6], voigt[3]
    ])
    return result
end

function plot_traction_displacement(u_max, traction_magnitude)
    plot(
        vcat(0.0, u_max),                # Add the origin as a point
        vcat(0.0, traction_magnitude),
        linewidth=2,
        title="Traction-displacement",
        label=nothing,
        markershape=:auto
    )
    ylabel!("Traction [Pa]")
    xlabel!("Maximum deflection [m]")
end

function plot_stress_strain_hist(strain_hist, stress_hist)
    plot(
        strain_hist, stress_hist,
        linewidth=2,
        title="Stress vs. Strain (at one point)",
        xlabel="Equivalent Strain",
        ylabel="Equivalent Stress [Pa]",
        label=nothing,
        markershape=:auto
    )
end

function plot_hardening_hist(hardening_hist)
    plot(
        1:length(hardening_hist), hardening_hist,
        linewidth=2,
        title="Hardening Variable (at one point)",
        xlabel="Timestep",
        ylabel="Hardening Variable",
        label=nothing,
        markershape=:auto
    )
end

function plot_strain_rate_hist(strain_hist, Δt = 0.1)
    strain_rate_hist = Float64[]
    for i in 2:length(strain_hist)
        rate = (strain_hist[i] - strain_hist[i-1]) / Δt
        push!(strain_rate_hist, rate)
    end
    plot(
        2:length(strain_hist), strain_rate_hist,
        linewidth=2,
        title="Strain Rate (at one point)",
        xlabel="Timestep",
        ylabel="Equivalent Strain Rate [1/s]",
        label=nothing,
        markershape=:auto
    )
end