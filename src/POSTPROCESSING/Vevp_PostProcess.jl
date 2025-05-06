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

function plot_stress_strain(states)
    strains = Float64[]
    stresses = Float64[]
    for cell_states in eachcol(states)
        for state in cell_states
            ϵ_voigt = state[10:15]
            ϵ = voigt_to_tensor(ϵ_voigt)
            σ_voigt = state[1:6]
            σ = voigt_to_tensor(σ_voigt)
            eq_strain = sqrt(2.0 / 3.0 * tr(ϵ ⊡ ϵ))
            eq_stress = vonMises(σ)
            push!(strains, eq_strain)
            push!(stresses, eq_stress)
        end
    end
    plot(
        strains, stresses,
        linewidth=2,
        title="Stress vs. Strain",
        xlabel="Equivalent Strain",
        ylabel="Equivalent Stress [Pa]",
        label=nothing,
        markershape=:auto
    )
end

function plot_strain_rate(states, states_old, Δt)
    strain_rates = Float64[]
    for (cell_states, cell_states_old) in zip(eachcol(states), eachcol(states_old))
        for (state, state_old) in zip(cell_states, cell_states_old)
            ϵ_voigt = state[10:15]
            ϵ_old_voigt = state_old[10:15]
            dϵ_voigt = (ϵ_voigt - ϵ_old_voigt) / Δt
            dϵ = voigt_to_tensor(dϵ_voigt)
            eq_strain_rate = sqrt(2.0 / 3.0 * tr(dϵ ⊡ dϵ))
            push!(strain_rates, eq_strain_rate)
        end
    end
    plot(
        1:length(strain_rates), strain_rates,
        linewidth=2,
        title="Strain Rate",
        xlabel="Integration Point Index",
        ylabel="Equivalent Strain Rate [1/s]",
        label=nothing,
        markershape=:auto
    )
end

function plot_hardening(states)
    hardening_values = Float64[]
    for cell_states in eachcol(states)
        for state in cell_states
            k = state[19]
            push!(hardening_values, k)
        end
    end
    plot(
        1:length(hardening_values), hardening_values,
        linewidth=2,
        title="Hardening Variable",
        xlabel="Integration Point Index",
        ylabel="Hardening Variable",
        label=nothing,
        markershape=:auto
    )
end