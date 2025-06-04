using Ferrite, Plots

function vonMises(σ::Vector{Float64})
    # Voigt: [σ11, σ22, σ33, σ12, σ13, σ23]
    s = σ .- sum(σ[1:3])/3  # deviatoric part (only for principal stresses)
    return sqrt(1.5 * (s[1]^2 + s[2]^2 + s[3]^2 + 2*(s[4]^2 + s[5]^2 + s[6]^2)))
end

function postprocess(grid, dh, states, props, u, filename="plasticity")
    n_cells = getncells(grid)
    n_qp = size(states, 1)

    mises_values = zeros(n_cells)
    sigma_11 = zeros(n_cells)
    disp_mag = zeros(getnnodes(grid))

    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            σ = state[101:106]  # Voigt stress from UMAT
            mises_values[el] += vonMises(σ)
            sigma_11[el] += σ[1]  # σ₁₁
            # If you want to extract other state variables, use their indices here
        end
        mises_values[el] /= length(cell_states)
        sigma_11[el] /= length(cell_states)
    end

    # Displacement magnitude at nodes
    for i in 1:getnnodes(grid)
        u_node = u[(3i-2):(3i)]
        disp_mag[i] = norm(u_node)
    end

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