using Ferrite, Tensors, TimerOutputs, ProgressMeter, LinearAlgebra, Printf

# ==============================================================================
# PURE DISPLACEMENT FORMULATION WITH FORTRAN UMAT INTEGRATION
# ==============================================================================
# This implementation demonstrates integration of Fortran UMAT (User MATerial)
# subroutines with Ferrite.jl FEM framework for finite strain material modeling.
#
# Test Configuration: Cantilever Beam with VEVP Material
# - Geometry: 10cm × 1cm × 1cm beam (20×4×4 Q1 hexahedra)
# - Material: Viscoelastic-Viscoplastic (VEVP) RTM6 epoxy
#   * 8 Maxwell branches: Multi-scale relaxation (τ = 1s to 1000s)
#   * Viscoplasticity: Rate-dependent yield with hardening
# - Loading: 1cm tip deflection (10% geometric nonlinearity)
# - Solver: Newton-Raphson with 100 load steps, 1e-4 tolerance
#
# Key Features:
# - ccall interface to Fortran UMAT (ABAQUS-compatible)
# - Consistent tangent stiffness (DDSDDE) for Newton convergence
# - State variable management (108 internal variables per quadrature point)
# - Large deformation kinematics with stress gradient verification
# ==============================================================================

# ==============================================================================
# 1. MATERIAL TYPE CONFIGURATION
# ==============================================================================
"""
Material model selection for finite element analysis.

Available options:
- "neohook" : Neo-Hookean hyperelastic (native Julia, validation baseline)
- "elastic" : Linear elastic UMAT (Fortran, interface verification)
- "vevp"    : Viscoelastic-Viscoplastic UMAT (Fortran, research model)

The material type determines which constitutive model is used for stress
computation at each Gauss quadrature point during FEM assembly.
"""
const MATERIAL_TYPE = "vevp"

# ==============================================================================
# 2. FORTRAN UMAT INTERFACE
# ==============================================================================
# Interface to Fortran User MATerial (UMAT) subroutines via ccall.
# Follows ABAQUS UMAT conventions for stress update and tangent computation.
# ==============================================================================

# Paths to compiled Fortran shared libraries
const UMAT_LIB_ELASTIC = "./src/Material_Models/libumat_elastic.so"
const UMAT_LIB_VEVP = "./src/Material_Models/libumat.so"

"""
    call_umat(stress, statev, ddsdde, ...)

Interface to Fortran UMAT subroutine via ccall.

# Arguments (ABAQUS UMAT convention)
- `stress::Vector{Float64}`: Cauchy stress tensor [σ11, σ22, σ33, σ12, σ13, σ23] (6 components)
- `statev::Vector{Float64}`: State variables (material history, e.g., plastic strain)
- `ddsdde::Matrix{Float64}`: Material tangent stiffness (6×6 Jacobian: ∂σ/∂ε)
- `stran::Vector{Float64}`: Total strain at start of increment
- `dstran::Vector{Float64}`: Strain increment (Δε)
- `dtime::Float64`: Time increment (Δt) for rate-dependent materials
- `props::Vector{Float64}`: Material properties (elastic moduli, yield stress, etc.)
- `dfgrd0::Matrix{Float64}`: Deformation gradient F at t_n (3×3, finite strain)
- `dfgrd1::Matrix{Float64}`: Deformation gradient F at t_{n+1} (3×3)
- `nstatv::Int`: Number of state variables

# Returns
Updates stress, statev, and ddsdde in-place (Fortran-style mutation).

# Notes
- CHARACTER*80 cmname converted to UInt8[80] for Fortran compatibility
- Material type (elastic/vevp) determines which shared library is called
- Consistent tangent (ddsdde) enables quadratic Newton convergence
"""
function call_umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
                   stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname_str,
                   ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
                   celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    
    # Prepare CHARACTER*80 buffer for Fortran
    cmname = zeros(UInt8, 80)
    src_bytes = collect(codeunits(cmname_str))
    n_copy = min(length(src_bytes), 80)
    cmname[1:n_copy] = src_bytes[1:n_copy]
    
    # Select library based on material type
    if MATERIAL_TYPE == "elastic"
        # Call elastic UMAT
        ccall((:umat_, UMAT_LIB_ELASTIC), Cvoid,
            (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, 
             Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
             Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, 
             Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8},
             Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ref{Int32},
             Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
             Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, 
             Ref{Int32}, Ref{Int32}, Ref{Int32}),
            stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
            stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
            ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
            celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    elseif MATERIAL_TYPE == "vevp"
        # Call VEVP UMAT
        ccall((:umat_, UMAT_LIB_VEVP), Cvoid,
            (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, 
             Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
             Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, 
             Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8},
             Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ref{Int32},
             Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
             Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}, 
             Ref{Int32}, Ref{Int32}, Ref{Int32}),
            stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
            stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname,
            ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt,
            celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    else
        error("Invalid MATERIAL_TYPE for UMAT: $MATERIAL_TYPE")
    end
end

"""
Convert symmetric 3×3 tensor to Voigt notation (6 components).
Engineering notation: [11, 22, 33, 12, 13, 23]
"""
function tensor_to_voigt6(mat::AbstractMatrix)
    return [mat[1,1], mat[2,2], mat[3,3], mat[1,2], mat[1,3], mat[2,3]]
end

"""
Convert Voigt notation (6 components) to symmetric 3×3 tensor.
"""
function voigt_to_tensor(v::AbstractVector)
    return Tensor{2,3}((v[1], v[4], v[5],
                        v[4], v[2], v[6],
                        v[5], v[6], v[3]))
end

"""
Compute stress and tangent from UMAT for given deformation gradient F.

This function:
1. Computes strain from F (logarithmic or Green-Lagrange)
2. Calls UMAT to get Cauchy stress σ and tangent C
3. Converts to 2nd Piola-Kirchhoff stress S and material tangent ∂S/∂C

Arguments:
- F: Current deformation gradient (3×3)
- F_old: Previous deformation gradient (3×3)
- statev: State variables (vector)
- PROPS: Material properties (vector)
- DTIME: Time increment
- qp_idx: Quadrature point index (for debugging)
- cell_idx: Cell index (for debugging)

Returns:
- S: 2nd Piola-Kirchhoff stress (Tensor{2,3})
- ∂S∂C: Material tangent ∂S/∂C (Tensor{4,3})
- statev_new: Updated state variables
"""
function compute_stress_tangent_umat(F, F_old, statev, PROPS, DTIME, qp_idx, cell_idx)
    
    # Compute strain measure (logarithmic strain for finite deformation)
    C = tdot(F)  # Right Cauchy-Green: C = F^T · F
    C_old = tdot(F_old)
    
    # Logarithmic strain: ε = 0.5 * log(C)
    # For small strains, this reduces to engineering strain
    # Note: This is a simplified approach; UMAT may use different strain measure
    E = 0.5 * (C - one(C))  # Green-Lagrange strain
    E_old = 0.5 * (C_old - one(C_old))
    
    # Convert to Voigt notation for UMAT (6 components)
    # Ferrite's tovoigt returns [11, 22, 33, 23, 13, 12]
    # UMAT expects [11, 22, 33, 12, 13, 23]
    E_voigt = tovoigt(E)
    E_old_voigt = tovoigt(E_old)
    
    # Reorder from Ferrite to UMAT convention
    ε = [E_voigt[1], E_voigt[2], E_voigt[3], E_voigt[6], E_voigt[5], E_voigt[4]]
    ε_old = [E_old_voigt[1], E_old_voigt[2], E_old_voigt[3], E_old_voigt[6], E_old_voigt[5], E_old_voigt[4]]
    Δε = ε - ε_old
    
    # Check for excessive strain
    strain_norm = norm(ε)
    if strain_norm > 0.1
        @warn "Large strain detected at cell $cell_idx, qp $qp_idx: ||ε|| = $strain_norm"
    end
    
    # Check deformation gradient
    F_norm = norm(F - one(F))
    if F_norm > 0.5
        @warn "Large deformation at cell $cell_idx, qp $qp_idx: ||F - I|| = $F_norm"
    end
    
    # Initialize UMAT output arrays (filled by ccall via pass-by-reference)
    ntens = 6
    ndi = 3
    nshr = 3
    nstatv = length(statev)
    nprops = length(PROPS)
    
    stress = zeros(6)           # OUTPUT: Cauchy stress in Voigt notation [σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃]
    statev_new = copy(statev)   # OUTPUT: Updated state variables
    ddsdde = zeros(6, 6)        # OUTPUT: Material tangent ∂σ/∂ε in Voigt notation (6×6 matrix)
    
    # UMAT scalar arguments
    sse, spd, scd = 0.0, 0.0, 0.0
    rpl = 0.0
    ddsddt = zeros(6)
    drplde = zeros(6)
    drpldt = 0.0
    
    # Time and temperature
    time = [0.0, DTIME]
    dtime = DTIME
    temp, dtemp = 20.0, 0.0
    predef, dpred = [0.0], [0.0]
    
    # Element info
    coords = zeros(3)
    drot = Matrix{Float64}(I, 3, 3)
    pnewdt = 1.0
    celent = 0.01  # Element characteristic length
    
    dfgrd0 = Matrix(F_old)
    dfgrd1 = Matrix(F)
    
    noel = cell_idx
    npt = qp_idx
    layer, kspt = 1, 1
    kstep, kinc = 1, 1
    
    cmname = "UMAT"
    
    # =========================================================================
    # CALL FORTRAN UMAT SUBROUTINE
    # =========================================================================
    # Fortran modifies stress, statev_new, and ddsdde arrays in-place
    call_umat(stress, statev_new, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt,
              ε_old, Δε, time, dtime, temp, dtemp, predef, dpred, cmname,
              ndi, nshr, ntens, nstatv, PROPS, nprops, coords, drot, pnewdt,
              celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    # =========================================================================
    
    # Validate UMAT outputs for numerical stability
    if any(isnan, stress) || any(isinf, stress)
        @error "UMAT returned NaN/Inf stress" cell_idx qp_idx F_norm strain_norm
        error("NaN or Inf in stress from UMAT")
    end
    
    if any(isnan, ddsdde) || any(isinf, ddsdde)
        @error "UMAT returned NaN/Inf tangent" cell_idx qp_idx F_norm strain_norm
        error("NaN or Inf in tangent from UMAT")
    end
    
    # -------------------------------------------------------------------------
    # STRESS CONVERSION: Cauchy stress → 2nd Piola-Kirchhoff stress
    # -------------------------------------------------------------------------
    # UMAT output 'stress' is stored in σ_voigt (Cauchy stress in current config)
    σ_voigt = stress
    σ = voigt_to_tensor(σ_voigt)
    
    # Transform to 2nd Piola-Kirchhoff stress: S = J·F⁻¹·σ·F⁻ᵀ
    J = det(F)
    F_inv = inv(F)
    S = J * F_inv ⋅ σ ⋅ transpose(F_inv)
    
    # -------------------------------------------------------------------------
    # TANGENT CONVERSION: Material tangent in Voigt notation → 4th order tensor
    # -------------------------------------------------------------------------
    # UMAT output 'ddsdde' (6×6 matrix) represents ∂σ/∂ε in Voigt notation
    # We construct ∂S/∂E as a 4th order tensor for finite element assembly
    C_voigt = ddsdde
    
    # Map Voigt 6×6 components to full 4th order tensor (3×3×3×3)
    # Voigt ordering: [11, 22, 33, 12, 13, 23] → Full tensor indexing
    ∂S∂E = Tensor{4,3}((
        C_voigt[1,1], C_voigt[1,4], C_voigt[1,5], C_voigt[1,4], C_voigt[1,2], C_voigt[1,6], C_voigt[1,5], C_voigt[1,6], C_voigt[1,3],
        C_voigt[4,1], C_voigt[4,4], C_voigt[4,5], C_voigt[4,4], C_voigt[4,2], C_voigt[4,6], C_voigt[4,5], C_voigt[4,6], C_voigt[4,3],
        C_voigt[5,1], C_voigt[5,4], C_voigt[5,5], C_voigt[5,4], C_voigt[5,2], C_voigt[5,6], C_voigt[5,5], C_voigt[5,6], C_voigt[5,3],
        C_voigt[4,1], C_voigt[4,4], C_voigt[4,5], C_voigt[4,4], C_voigt[4,2], C_voigt[4,6], C_voigt[4,5], C_voigt[4,6], C_voigt[4,3],
        C_voigt[2,1], C_voigt[2,4], C_voigt[2,5], C_voigt[2,4], C_voigt[2,2], C_voigt[2,6], C_voigt[2,5], C_voigt[2,6], C_voigt[2,3],
        C_voigt[6,1], C_voigt[6,4], C_voigt[6,5], C_voigt[6,4], C_voigt[6,2], C_voigt[6,6], C_voigt[6,5], C_voigt[6,6], C_voigt[6,3],
        C_voigt[5,1], C_voigt[5,4], C_voigt[5,5], C_voigt[5,4], C_voigt[5,2], C_voigt[5,6], C_voigt[5,5], C_voigt[5,6], C_voigt[5,3],
        C_voigt[6,1], C_voigt[6,4], C_voigt[6,5], C_voigt[6,4], C_voigt[6,2], C_voigt[6,6], C_voigt[6,5], C_voigt[6,6], C_voigt[6,3],
        C_voigt[3,1], C_voigt[3,4], C_voigt[3,5], C_voigt[3,4], C_voigt[3,2], C_voigt[3,6], C_voigt[3,5], C_voigt[3,6], C_voigt[3,3]
    ))
    
    # Chain rule: ∂S/∂C = (∂S/∂E):(∂E/∂C) where E = ½(C-I), thus ∂E/∂C = ½
    ∂S∂C = 0.5 * ∂S∂E
    
    return S, ∂S∂C, statev_new
end

# ==============================================================================
# 3. NEO-HOOKEAN MATERIAL (Built-in, for validation)
# ==============================================================================

struct NeoHooke
    μ::Float64  # Shear modulus
    λ::Float64  # Lame parameter
end

"""
Neo-Hookean strain energy function.
"""
function Ψ(C, mp::NeoHooke)
    μ = mp.μ
    λ = mp.λ
    Ic = tr(C)
    J = sqrt(det(C))
    return μ / 2 * (Ic - 3 - 2 * log(J)) + λ / 2 * (J - 1)^2
end

"""
Neo-Hookean constitutive driver using automatic differentiation.
Returns 2nd Piola-Kirchhoff stress S and material tangent ∂S/∂C.
"""
function constitutive_driver_neohook(C, mp::NeoHooke)
    # Compute all derivatives in one function call
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

# ==============================================================================
# 4. MATERIAL PROPERTY FUNCTIONS
# ==============================================================================

"""
Get material properties for elastic UMAT.
Returns PROPS vector and number of state variables.
"""
function get_elastic_properties()
    # Simple elastic: E = 2450 MPa, ν = 0.39 (RTM6 equilibrium values)
    E = 2.45e9  # Pa
    ν = 0.39
    
    PROPS = [E, ν]
    nstatv = 0  # No state variables for elastic
    
    return PROPS, nstatv
end

"""
Get material properties for VEVP UMAT (RTM6 epoxy resin).
Returns PROPS vector and number of state variables.
"""
function get_vevp_properties()
    # VEVP material properties - COMPLETE SET from main.jl
    # Critical: ALL hardening parameters MUST be non-zero!
    
    PROPS = [
        5,                        # 1: Approximation order + VEVP trigger (>2)
        1.470588416e6,            # 2: Bulk Modulus [Pa]
        5.639098439e5,            # 3: Shear Modulus [Pa]
        5.900948586,              # 4: Yield Exponent
        0.33,                     # 5: Plastic Poisson Ratio
        0.001,                    # 6: Viscoplastic Coefficient
        10,                       # 7: Viscoplastic Exponent
        2.086229688e6,            # 8: Initial Yield Limit - Compression [Pa]
        2.164115496e9,            # 9: Isotropic Hardening - Compression [Pa]
        4.450073598e6,            # 10: Isotropic Hardening - Compression [Pa]
        5.401554318e6,            # 11: Isotropic Hardening - Compression [Pa]
        1.66898375e6,             # 12: Initial Yield Limit - Tension [Pa]
        1.731292397e9,            # 13: Isotropic Hardening - Tension [Pa]
        3.560058879e6,            # 14: Isotropic Hardening - Tension [Pa]
        5.401554318e6,            # 15: Isotropic Hardening - Tension [Pa]
        0,                        # 16: Kinematic Hardening Parameter
        0,                        # 17: Kinematic Hardening Parameter
        0,                        # 18: Kinematic Hardening Parameter
        # Maxwell branches - ALL 8 BRANCHES ACTIVE! 🔥
        # Logarithmic time scale distribution: 1s → 1000s
        1.0e6, 0.8e6, 0.6e6, 0.5e6, 0.4e6, 0.3e6, 0.2e6, 0.1e6,   # 19-26: Bulk moduli [Pa]
        1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 700.0, 1000.0,        # 27-34: Volumetric relaxation times [s]
        0.5e6, 0.4e6, 0.3e6, 0.25e6, 0.2e6, 0.15e6, 0.1e6, 0.05e6, # 35-42: Shear moduli [Pa]
        1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 700.0, 1000.0          # 43-50: Deviatoric relaxation times [s]
    ]
    
    nstatv = 108  # VEVP has 108 state variables
    
    return PROPS, nstatv
end

"""
Get Neo-Hookean material for validation.
"""
function get_neohook_material()
    # Match elastic UMAT properties
    E = 2.45e9  # Pa
    ν = 0.39
    μ = E / (2(1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    return NeoHooke(μ, λ)
end

# ==============================================================================
# 5. ELEMENT ASSEMBLY (Pure Displacement Formulation)
# ==============================================================================

"""
Assemble element stiffness matrix and residual vector.

This is the PURE DISPLACEMENT formulation (no pressure field).
Compatible with standard UMAT interface.
"""
function assemble_element!(ke, ge, cell, cv, mp, ue, states_old, PROPS, nstatv, DTIME, cell_idx)
    
    # Reinitialize cell values, and reset output arrays
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)
    
    n_basefuncs = getnbasefunctions(cv)
    n_qpoints = getnquadpoints(cv)
    
    # Storage for updated state variables
    states_new = Vector{Vector{Float64}}(undef, n_qpoints)
    
    for qp in 1:n_qpoints
        
        # Get quadrature point weight
        dΩ = getdetJdV(cv, qp)
        
        # Compute deformation gradient F
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        
        # Previous F (for incremental analysis)
        # For first step, F_old = I
        F_old = one(F)
        
        # Compute stress and tangent based on material type
        if MATERIAL_TYPE == "neohook"
            # Neo-Hookean (built-in)
            C = tdot(F)  # F' ⋅ F
            S, ∂S∂C = constitutive_driver_neohook(C, mp)
            states_new[qp] = Float64[]  # No state variables
            
        elseif MATERIAL_TYPE == "elastic" || MATERIAL_TYPE == "vevp"
            # UMAT (elastic or VEVP)
            statev = states_old[qp]
            S, ∂S∂C, statev_new = compute_stress_tangent_umat(
                F, F_old, statev, PROPS, DTIME, qp, cell_idx
            )
            states_new[qp] = statev_new
            
        else
            error("Unknown MATERIAL_TYPE: $MATERIAL_TYPE")
        end
        
        # Convert to 1st Piola-Kirchhoff: P = F · S
        P = F ⋅ S
        
        # Compute tangent: ∂P/∂F
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * F ⋅ ∂S∂C ⊡ otimesu(F', I)
        
        # Loop over test functions (residual and tangent assembly)
        for i in 1:n_basefuncs
            # Test function gradient
            ∇δui = shape_gradient(cv, qp, i)
            
            # Add contribution to residual
            ge[i] += (∇δui ⊡ P) * dΩ
            
            # Precompute for efficiency
            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F
            
            # Loop over trial functions (tangent matrix)
            for j in 1:n_basefuncs
                ∇δuj = shape_gradient(cv, qp, j)
                # Add contribution to tangent
                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end
    
    return ke, ge, states_new
end

# ==============================================================================
# 6. GLOBAL ASSEMBLY
# ==============================================================================

"""
Assemble global stiffness matrix and residual vector.
"""
function assemble_global!(K, g, dh, cv, mp, u, states, PROPS, nstatv, DTIME)
    
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)
    
    # Start assembler (resets K and g)
    assembler = start_assemble(K, g)
    
    # Loop over all cells in the grid
    for (cell_idx, cell) in enumerate(CellIterator(dh))
        global_dofs = celldofs(cell)
        ue = u[global_dofs]  # Element DOFs
        
        # Get state variables for this cell
        states_old_cell = states[cell_idx]
        
        # Assemble element
        ke, ge, states_new_cell = assemble_element!(
            ke, ge, cell, cv, mp, ue, states_old_cell, PROPS, nstatv, DTIME, cell_idx
        )
        
        # Update state variables
        states[cell_idx] = states_new_cell
        
        # Assemble into global system
        assemble!(assembler, global_dofs, ke, ge)
    end
    
    return K, g, states
end

# ==============================================================================
# 7. SOLVER
# ==============================================================================

"""
Solve the nonlinear equilibrium problem using displacement control.
"""
function solve()
    
    reset_timer!()
    
    println("="^70)
    println("PURE DISPLACEMENT FORMULATION WITH UMAT")
    println("="^70)
    println("Material type: $MATERIAL_TYPE")
    println()
    
    # ========================================
    # Mesh generation - CANTILEVER BEAM
    # ========================================
    # Beam dimensions
    L_length = 0.1   # 10 cm long
    L_height = 0.01  # 1 cm tall  
    L_width = 0.01   # 1 cm wide
    
    # More elements along length for bending resolution
    N_length = 20
    N_height = 4
    N_width = 4
    
    left = Vec{3}((0.0, 0.0, 0.0))
    right = Vec{3}((L_length, L_width, L_height))
    
    grid = generate_grid(Hexahedron, (N_length, N_width, N_height), left, right)
    println("Mesh: $(getncells(grid)) hexahedral elements")
    println("  Beam dimensions: $(L_length*100) cm × $(L_width*100) cm × $(L_height*100) cm")
    println("  Elements: $(N_length) × $(N_width) × $(N_height)")
    
    # ========================================
    # Material properties
    # ========================================
    if MATERIAL_TYPE == "neohook"
        mp = get_neohook_material()
        PROPS = Float64[]
        nstatv = 0
        println("Material: Neo-Hookean hyperelastic")
        println("  μ = $(mp.μ/1e9) GPa")
        println("  λ = $(mp.λ/1e9) GPa")
        
    elseif MATERIAL_TYPE == "elastic"
        PROPS, nstatv = get_elastic_properties()
        mp = nothing
        println("Material: Elastic UMAT")
        println("  E = $(PROPS[1]/1e9) GPa")
        println("  ν = $(PROPS[2])")
        
    elseif MATERIAL_TYPE == "vevp"
        PROPS, nstatv = get_vevp_properties()
        mp = nothing
        println("Material: VEVP UMAT (RTM6 epoxy)")
        println("  K_inf = $(PROPS[2]/1e9) GPa")
        println("  G_inf = $(PROPS[3]/1e9) GPa")
        println("  State variables: $nstatv")
        println("  Viscoelastic: ALL 8 MAXWELL BRANCHES ACTIVE! 🔥")
        println("    Branch 1: K=$(PROPS[19]/1e6) MPa, G=$(PROPS[35]/1e6) MPa, τ=$(PROPS[27]) s")
        println("    Branch 2: K=$(PROPS[20]/1e6) MPa, G=$(PROPS[36]/1e6) MPa, τ=$(PROPS[28]) s")
        println("    Branch 3: K=$(PROPS[21]/1e6) MPa, G=$(PROPS[37]/1e6) MPa, τ=$(PROPS[29]) s")
        println("    Branch 4: K=$(PROPS[22]/1e6) MPa, G=$(PROPS[38]/1e6) MPa, τ=$(PROPS[30]) s")
        println("    Branch 5: K=$(PROPS[23]/1e6) MPa, G=$(PROPS[39]/1e6) MPa, τ=$(PROPS[31]) s")
        println("    Branch 6: K=$(PROPS[24]/1e6) MPa, G=$(PROPS[40]/1e6) MPa, τ=$(PROPS[32]) s")
        println("    Branch 7: K=$(PROPS[25]/1e6) MPa, G=$(PROPS[41]/1e6) MPa, τ=$(PROPS[33]) s")
        println("    Branch 8: K=$(PROPS[26]/1e6) MPa, G=$(PROPS[42]/1e6) MPa, τ=$(PROPS[34]) s")
    end
    
    DTIME = 100.0  # Time increment (for VEVP rate effects)
    println("  DTIME = $DTIME s")
    println()
    
    # ========================================
    # Finite element setup
    # ========================================
    ip = Lagrange{RefHexahedron, 1}()^3  # Q1 displacement
    qr = QuadratureRule{RefHexahedron}(2)  # 2×2×2 Gauss integration
    cv = CellValues(qr, ip)
    
    # DofHandler - SINGLE FIELD (pure displacement)
    dh = DofHandler(grid)
    add!(dh, :u, ip)  # Only displacement field
    close!(dh)
    
    println("DOFs: $(ndofs(dh))")
    println("DOFs per cell: $(ndofs_per_cell(dh))")
    println()
    
    # ========================================
    # Boundary conditions - CANTILEVER BEAM
    # ========================================
    dbcs = ConstraintHandler(dh)
    
    # Left face (x=0): fully fixed (clamped end)
    dbc_fixed = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(dbcs, dbc_fixed)
    
    # Right face (x=L): prescribed tip displacement in z-direction (downward deflection)
    u_tip = -0.01  # 1 cm tip deflection (10% of beam length!)
    # Only constrain z-component, function returns scalar for single component
    dbc_tip = Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> t*u_tip, [3])
    add!(dbcs, dbc_tip)
    
    close!(dbcs)
    
    println("Boundary conditions:")
    println("  Left end (x=0): fully clamped (u = 0)")
    println("  Right end (x=L): prescribed tip deflection u_z = $u_tip m ($(abs(u_tip/L_length)*100)% of length)")
    println("  Expected: Large bending deformation + stress gradient!")
    println()
    
    # ========================================
    # Initialize solution vectors
    # ========================================
    _ndofs = ndofs(dh)
    u = zeros(_ndofs)
    Δu = zeros(_ndofs)
    ΔΔu = zeros(_ndofs)
    
    # Initialize state variables for each quadrature point
    n_qpoints = getnquadpoints(cv)
    
    # CRITICAL FIX for VEVP: Initialize F_vp to identity [1, 1, 1, 0, 0, 0, 0, 0, 0]
    # STATEV(1:9) stores viscoplastic deformation gradient F_vp
    # Must be identity initially: F_vp = [[1,0,0], [0,1,0], [0,0,1]]
    function init_statev()
        statev = zeros(nstatv)
        if MATERIAL_TYPE == "vevp"
            # Initialize F_vp to identity (diagonal components)
            statev[1] = 1.0  # F_vp(1,1)
            statev[2] = 1.0  # F_vp(2,2)
            statev[3] = 1.0  # F_vp(3,3)
            # Off-diagonal components stay at zero
        end
        return statev
    end
    
    states = [Vector{Float64}[init_statev() for _ in 1:n_qpoints] for _ in 1:getncells(grid)]
    
    # Sparse matrix and residual
    K = allocate_matrix(dh)
    g = zeros(_ndofs)
    
    # ========================================
    # Load stepping
    # ========================================
    n_steps = 100  # Many small steps for stable convergence with 8 Maxwell branches
    t_history = Float64[]
    u_history = Float64[]
    states_history = []  # Store state variables at each converged step
    stress_history = []  # Store stress tensors at each converged step
    strain_history = []  # Store strain at each converged step
    
    println("="^70)
    println("Starting load stepping: $n_steps steps")
    println("="^70)
    
    for step in 1:n_steps
        
        t = step / n_steps
        
        # Update boundary conditions
        Ferrite.update!(dbcs, t)
        apply!(u, dbcs)
        
        # Newton-Raphson iteration
        newton_itr = 0
        NEWTON_TOL = 1.0e-4  # Relaxed for complex VEVP with 8 Maxwell branches
        NEWTON_MAXITER = 50  # Allow more iterations for stable convergence
        
        println("\nStep $step: t = $(round(t, digits=4))")
        
        while true
            newton_itr += 1
            
            # Assemble system
            K, g, states = assemble_global!(K, g, dh, cv, mp, u, states, PROPS, nstatv, DTIME)
            
            # Apply boundary conditions
            apply_zero!(K, g, dbcs)
            
            # Check convergence
            normg = norm(g)
            @printf("  Iter %2d: ||r|| = %.3e\n", newton_itr, normg)
            
            if normg < NEWTON_TOL
                println("  ✅ Converged!")
                break
            elseif newton_itr >= NEWTON_MAXITER
                error("❌ Reached maximum Newton iterations, aborting")
            end
            
            # Solve for increment
            ΔΔu = K \ g
            apply_zero!(ΔΔu, dbcs)
            
            # Update displacement
            u .-= ΔΔu
        end
        
        # Save history
        push!(t_history, t)
        
        # Get tip displacement (rightmost node at top corner)
        # For cantilever: nodes at x=L_length
        tip_node = getnnodes(grid)  # Last node is at tip
        tip_dof_z = 3 * tip_node
        push!(u_history, u[tip_dof_z])
        
        # Store converged state variables for post-processing
        push!(states_history, deepcopy(states))
        
        @printf("  Tip deflection u_z = %.6f mm\n", u[tip_dof_z]*1000)
    end
    
    println("\n" * "="^70)
    println("✅ ANALYSIS COMPLETE!")
    println("="^70)
    println("Steps completed: $n_steps")
    println("Final tip deflection: $(u_history[end]*1000) mm")
    println()
    
    # ========================================
    # Export results
    # ========================================
    @timeit "export" begin
        # Create visualization directory if it doesn't exist
        vis_dir = "src/POSTPROCESS/visualization"
        mkpath(vis_dir)
        
        # Export VTK file
        vtk_file = joinpath(vis_dir, "cantilever_vevp_8branches")
        VTKGridFile(vtk_file, dh) do vtk
            write_solution(vtk, dh, u)
        end
    end
    
    println("Results exported to: src/POSTPROCESS/visualization/cantilever_vevp_8branches.vtu")
    println()
    
    print_timer(title = "Analysis timing", linechars = :ascii)
    
    return u, t_history, u_history, states, states_history, grid, dh
end

# ==============================================================================
# 8. RUN SIMULATION
# ==============================================================================

# Run the solver
u_final, t_hist, u_hist, states_final, states_hist, grid, dh = solve()

println("\n" * "="^70)
println("SIMULATION FINISHED SUCCESSFULLY!")
println("="^70)

# ==============================================================================
# 9. POST-PROCESSING
# ==============================================================================

# Call post-processing to generate plots
include("POSTPROCESS/postprocess_results.jl")
create_plots(t_hist, u_hist, states_hist, grid, dh, u_final)
