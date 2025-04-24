using Ferrite, Plots

function vonMises(σ)
    s = dev(σ)
    return sqrt(3.0 / 2.0 * s ⊡ s)
end

function postprocess(grid, dh, states, props, u, filename="plasticity")
    # Extract material properties from PROPS
    H = props[5] # Hardening modulus

    # Compute von Mises stress and drag stress
    mises_values = zeros(getncells(grid))
    κ_values = zeros(getncells(grid))
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            mises_values[el] += vonMises(state.σ)
            κ_values[el] += state.k * H
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