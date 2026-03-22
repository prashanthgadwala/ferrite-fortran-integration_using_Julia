# Technical Reference: VEVP Implementation

Technical documentation for developers and researchers working with the VEVP material model integration.

---

## Table of Contents
1. [System Architecture](#system-architecture)
2. [VEVP Material Model](#vevp-material-model)
3. [Implementation Details](#implementation-details)
4. [Current Test Configuration](#current-test-configuration)
5. [Convergence Analysis](#convergence-analysis)
6. [Validation Strategy](#validation-strategy)

---

## System Architecture

### Integration Layers

```
Julia (Ferrite.jl) ←→ ccall ←→ Fortran UMAT ←→ Material physics
     ↓                                              ↓
  FEM solver                                  Stress integration
  Newton iteration                            State evolution
```

**Key responsibilities:**
- **Julia**: Mesh, assembly, linear solver, convergence
- **Fortran**: Material constitutive law, stress update, tangent

---

## VEVP Material Model

### Mathematical Formulation

**Total Response:**
```
σ = σ_equilibrium + σ_viscoelastic + σ_viscoplastic
```

### 1. Equilibrium (Long-Term Response)

```
K_∞ = 1.47 MPa    (bulk modulus)
G_∞ = 0.564 MPa   (shear modulus)
E_∞ ≈ 1.5 MPa     (Young's modulus)
ν ≈ 0.33          (Poisson's ratio)
```

### 2. Viscoelasticity (8 Maxwell Branches)

**Stress evolution:**
```
K(t) = K_∞ + Σᵢ Kᵢ exp(-t/τᵢᵛᵒˡ)
G(t) = G_∞ + Σᵢ Gᵢ exp(-t/τᵢᵈᵉᵛ)
```

**Branch Parameters (RTM6 Epoxy):**

| Branch | Kᵢ (MPa) | τᵛᵒˡ (s) | Gᵢ (MPa) | τᵈᵉᵛ (s) | Physics |
|--------|----------|----------|----------|----------|---------|
| 1      | 1.0      | 1        | 0.5      | 1        | Fast molecular rearrangement |
| 2      | 0.8      | 3        | 0.4      | 3        | Chain segment motion |
| 3      | 0.6      | 10       | 0.3      | 10       | Local polymer relaxation |
| 4      | 0.5      | 30       | 0.25     | 30       | Cooperative motion |
| 5      | 0.4      | 100      | 0.2      | 100      | Network rearrangement |
| 6      | 0.3      | 300      | 0.15     | 300      | Long-range relaxation |
| 7      | 0.2      | 700      | 0.1      | 700      | Structural adaptation |
| 8      | 0.1      | 1000     | 0.05     | 1000     | Matrix creep |

**Total Initial Stiffness:**
```
K₀ = K_∞ + ΣKᵢ = 1.47 + 3.8 = 5.27 MPa  (3.6× equilibrium)
G₀ = G_∞ + ΣGᵢ = 0.564 + 1.8 = 2.364 MPa (4.2× equilibrium)
```

### 3. Viscoplasticity (Rate-Dependent Yield)

**Flow rule:**
```
ε̇ᵛᵖ = γ₀ (σₑq/σy)ⁿ
```

**Hardening:**
```
σy = σy0 + H εpᵐ
```

**Parameters:**
```
σy0 = 2.08 MPa (compression) / 1.67 MPa (tension)
γ₀  = 0.001 s⁻¹
n   = 10 (rate sensitivity)
H   = 2.16 GPa (hardening modulus)
m   = 5.4 (hardening exponent)
```

### State Variables (108 total per Gauss point)

```
statev[1-9]:     F_vp (viscoplastic deformation gradient)
statev[10-18]:   Stress components
statev[19-26]:   Maxwell branch 1 internal variables
statev[27-34]:   Maxwell branch 2 internal variables
...
statev[91-98]:   Maxwell branch 8 internal variables
statev[99-108]:  Additional history variables
```

---

## Implementation Details

### 1. UMAT Interface (Julia ↔ Fortran)

**Location:** `src/main.jl` lines 47-170

**Key function:** `call_umat()`

**Conversion steps:**
```julia
# 1. Compute deformation gradient
F = I + ∇u

# 2. Strain (Green-Lagrange)
E = 0.5(F'F - I)

# 3. Voigt notation (6-component)
ε = [E11, E22, E33, 2E12, 2E13, 2E23]  # Factor of 2 for shear!

# 4. Call Fortran
ccall((:umat_, UMAT_LIB), Cvoid, (...))

# 5. Convert stress back to tensor
σ = [σ11 σ12 σ13]
    [σ12 σ22 σ23]
    [σ13 σ23 σ33]
```

**Critical details:**
- **CHARACTER*80**: Convert Julia string to `UInt8[80]` buffer
- **Column-major**: Fortran matrices transposed from Julia
- **Voigt convention**: Factor of 2 for shear components (engineering strain)
- **DDSDDE**: 6×6 tangent matrix (∂σ/∂ε) in Voigt notation

### 2. FEM Assembly

**Element loop** → Gauss point loop → UMAT call → Accumulate K, g

```julia
FOR each element:
  FOR each Gauss point:
    ∇u = shape_function_gradients · u_elem
    F = I + ∇u
    σ, DDSDDE = call_umat(F, Δt, props, statev)
    
    Ke += B' · DDSDDE · B · detJ · w
    ge += B' · σ · detJ · w
  END
END
```

**B matrix**: Strain-displacement operator (6×24 for Q1 hex)

### 3. Newton-Raphson

**Per load step:**
```julia
WHILE ||residual|| > tolerance:
  K, g = assemble_global()
  Δu = K \ g
  u -= Δu
END
```

**Settings**: tolerance = 1e-4, max_iter = 50

---

## Current Test Configuration

### Geometry
- **Cantilever beam**: 10 cm × 1 cm × 1 cm
- **Mesh**: 20×4×4 Q1 hexahedra (320 elements, 1575 DOFs)
- **Quadrature**: 2×2×2 Gauss points (2560 total, 108 state vars each)

### Loading
- **Displacement control**: 1 cm tip deflection (10% strain)
- **Load steps**: 100 (Δu = 0.1 mm per step)
- **Time**: Δt = 100 s per step (10,000 s total)

### Boundary Conditions
- Left face: Fully clamped (u = 0)
- Right face: uz = -10 mm prescribed

---

## Convergence Analysis

### Performance
- **Per step**: ~3-5 seconds (27 iterations × 2560 UMAT calls)
- **Full simulation**: ~6 minutes (100 steps)
- **Memory**: < 5 MB total

### Convergence Rate
- **Typical**: 20-30 iterations per step
- **Rate**: Linear (not quadratic)
- **Reason**: Approximate UMAT tangent (numerical, not analytical)

**Comparison**: ABAQUS/ANSYS show similar behavior (20-40 iterations) for complex materials

---

## Validation Strategy

1. **Analytical**: Compare tip deflection with beam theory (δ = FL³/3EI)
2. **Stress gradient**: Verify linear through-thickness distribution
3. **Relaxation**: Check multi-exponential force decay
4. **State variables**: Validate F_vp evolution and incompressibility

---

## Known Limitations

1. **Linear Newton convergence** (not quadratic) - approximate UMAT tangent
2. **Many load steps required** (100) - high initial stiffness from 8 branches
3. **Tolerance compromise** (1e-4) - cannot reach 1e-8 without 100+ iterations

**Future work**: Analytical tangent implementation for quadratic convergence

---

## References

- Zaïri et al. (2008): VEVP constitutive equations for rubber-modified polymers
- Belytschko et al. (2000): Nonlinear Finite Elements
- Ferrite.jl: https://ferrite-fem.github.io/Ferrite.jl/
- ABAQUS UMAT: User Subroutines Reference Guide

---

**For usage instructions**, see **[USER_GUIDE.md](USER_GUIDE.md)**
