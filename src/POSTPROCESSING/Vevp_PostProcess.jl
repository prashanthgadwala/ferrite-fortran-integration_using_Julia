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
    # Extract material properties
    H = props[5] # Hardening modulus

    mises_values = zeros(getncells(grid))
    κ_values = zeros(getncells(grid))

    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ_voigt = state[1:6]
            σ = voigt_to_tensor(σ_voigt) # Convert Voigt to tensor

            mises_values[el] += vonMises(σ) # von Mises stress (scalar)
            κ_values[el] += state[7] * H[1] # Ensure state[7] is a scalar
        end
        mises_values[el] /= length(cell_states) # Average von Mises stress
        κ_values[el] /= length(cell_states)     # Average drag stress
    end

    # Export results to VTK
    VTKGridFile(filename, dh) do vtk
        write_solution(vtk, dh, u) # Displacement field
        write_cell_data(vtk, mises_values, "von Mises [Pa]")
        write_cell_data(vtk, κ_values, "Drag stress [Pa]")
    end

    # Additional plots
    plot_stress_strain(states, filename="Stress_Strain")
    plot_strain_rate(states, states_old, 0.1, filename="Strain_Rate")
    plot_plastic_strain(states, filename="Plastic_Strain")
    plot_hardening(states, filename="Hardening")
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

function plot_stress_strain(states; filename="Stress_Strain")
    strains = []
    stresses = []

    for cell_states in eachcol(states)
        for state in cell_states
            # Extract strain tensor (e.g., indices 8–13 in Voigt notation)
            ϵ_voigt = state[8:13]
            ϵ = voigt_to_tensor(ϵ_voigt) # Convert Voigt to tensor

            # Extract stress tensor (first 6 elements in Voigt notation)
            σ_voigt = state[1:6]
            σ = voigt_to_tensor(σ_voigt) # Convert Voigt to tensor

            # Compute equivalent strain and stress
            eq_strain = sqrt(2.0 / 3.0 * tr(ϵ ⊡ ϵ)) # Equivalent strain
            eq_stress = vonMises(σ)                 # von Mises stress

            push!(strains, eq_strain)
            push!(stresses, eq_stress)
        end
    end

    # Plot stress vs. strain
    plot(
        strains, stresses,
        linewidth=2,
        title="Stress vs. Strain",
        xlabel="Equivalent Strain",
        ylabel="Equivalent Stress [Pa]",
        label=nothing,
        markershape=:auto
    )
    savefig("$filename.png")
end

function plot_strain_rate(states, states_old, Δt; filename="Strain_Rate")
    strain_rates = []

    for (cell_states, cell_states_old) in zip(eachcol(states), eachcol(states_old))
        for (state, state_old) in zip(cell_states, cell_states_old)
            # Extract strain tensors (e.g., indices 8–13 in Voigt notation)
            ϵ_voigt = state[8:13]
            ϵ_old_voigt = state_old[8:13]

            # Compute strain rate (difference in strain divided by Δt)
            dϵ_voigt = (ϵ_voigt - ϵ_old_voigt) / Δt
            dϵ = voigt_to_tensor(dϵ_voigt) # Convert Voigt to tensor

            # Compute equivalent strain rate
            eq_strain_rate = sqrt(2.0 / 3.0 * tr(dϵ ⊡ dϵ))
            push!(strain_rates, eq_strain_rate)
        end
    end

    # Plot strain rate
    plot(
        1:length(strain_rates), strain_rates,
        linewidth=2,
        title="Strain Rate",
        xlabel="Integration Point Index",
        ylabel="Equivalent Strain Rate [1/s]",
        label=nothing,
        markershape=:auto
    )
    savefig("$filename.png")
end

function plot_plastic_strain(states; filename="Plastic_Strain")
    plastic_strains = []

    for cell_states in eachcol(states)
        for state in cell_states
            # Extract plastic strain tensor (e.g., indices 14–19 in Voigt notation)
            ϵp_voigt = state[14:19]
            ϵp = voigt_to_tensor(ϵp_voigt) # Convert Voigt to tensor

            # Compute equivalent plastic strain
            eq_plastic_strain = sqrt(2.0 / 3.0 * tr(ϵp ⊡ ϵp))
            push!(plastic_strains, eq_plastic_strain)
        end
    end

    # Plot plastic strain
    plot(
        1:length(plastic_strains), plastic_strains,
        linewidth=2,
        title="Plastic Strain",
        xlabel="Integration Point Index",
        ylabel="Equivalent Plastic Strain",
        label=nothing,
        markershape=:auto
    )
    savefig("$filename.png")
end

function plot_hardening(states; filename="Hardening")
    hardening_values = []

    for cell_states in eachcol(states)
        for state in cell_states
            # Extract hardening variable (e.g., index 7)
            k = state[7]
            push!(hardening_values, k)
        end
    end

    # Plot hardening variable
    plot(
        1:length(hardening_values), hardening_values,
        linewidth=2,
        title="Hardening Variable",
        xlabel="Integration Point Index",
        ylabel="Hardening Variable",
        label=nothing,
        markershape=:auto
    )
    savefig("$filename.png")
end

