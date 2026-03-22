# Technical Reference: VEVP Julia-Fortran Integration

Technical reference for developers and researchers working on the UMAT bridge and nonlinear FEM workflow.

## System Architecture

```text
Julia (Ferrite.jl solver/assembly)
        <-> ccall
Fortran UMAT (stress, tangent, state update)
```

Responsibilities:

- Julia (`src/main.jl`): grid, DOFs, boundary conditions, Newton loop, state management
- Fortran (`src/Material_Models/umat.f`): constitutive integration and material tangent
- Post-processing (`src/POSTPROCESS/postprocess_results.jl`): plots + VTK time series

## Benchmark and Model

Current implementation target:

- Unit-cube uniaxial tension benchmark
- VEVP material with 8 Maxwell branches
- 108 state variables per quadrature point

Representative parameter blocks are defined in `get_vevp_properties()` in `src/main.jl`:

- Equilibrium moduli: `K_inf`, `G_inf`
- Viscoplastic parameters: yield/hardening/rate terms
- Branch moduli and relaxation times: `1, 3, 10, 30, 100, 300, 700, 1000` seconds

## UMAT Interface Details

The core bridge is `call_umat(...)` in `src/main.jl`.

Important implementation points:

- Fortran symbol call via `ccall((:umat_, UMAT_LIB_VEVP), ...)`
- Character argument handling via explicit 80-byte buffer (`UInt8[80]`)
- In-place updates for `stress`, `statev`, and `ddsdde`
- Voigt mapping and tensor conversions handled in Julia before/after UMAT call

## FEM Setup in Source Code

From `solve()` in `src/main.jl`:

- Geometry: `L_length = L_width = L_height = 1.0`
- Mesh: `N_length = N_width = N_height = 2` (8 hexahedral cells)
- Interpolation: `Q2` displacement (`Lagrange{RefHexahedron, 2}()^3`)
- Quadrature: `3 x 3 x 3`
- BCs:
  - left face: `u_x = 0`
  - right face: displacement ramp in `x`
  - one corner pin for `u_y = 0`, `u_z = 0` (RBM suppression)

Default runtime controls in source:

- `u_load = 0.3`
- `n_steps = 100`
- `DTIME = 0.1`
- Newton loop with backtracking line search (`NEWTON_MAXITER = 150`)

## State Variables (108)

The code and report use the following practical layout convention:

- `1:9` viscoplastic deformation gradient components (`F_vp` in Voigt-style packing)
- `10:18` viscoelastic strain components (`E_ve`)
- `19` accumulated plastic measure (`gamma`)
- remaining entries store branch/history variables used by the UMAT

Initialization in `src/main.jl` sets `F_vp` diagonal entries to identity at start.

## Post-Processing Pipeline

`create_plots(...)` in `src/POSTPROCESS/postprocess_results.jl` writes:

- `src/POSTPROCESS/plots/force_displacement.png`
- `src/POSTPROCESS/plots/stress_strain.png`
- `src/POSTPROCESS/plots/displacement_history.png`
- `src/POSTPROCESS/plots/deformed_shape.png`
- `src/POSTPROCESS/plots/vtk_timesteps/results.pvd` and stepwise `.vtu`

The VTK export includes displacement and selected cell-averaged internal variables/stresses (for ParaView animation and inspection).

## Validation Notes

Project-report validation (same UMAT family against ABAQUS) reports:

- Mean normalized stress-shape error: `3.71%`
- Maximum normalized stress-shape error: `5.38%`

The report benchmark setup uses 50 increments with physical time scaling (`dt = 100 s`, total `5000 s`).

When reproducing report figures, ensure your runtime controls match report settings instead of development defaults.

## Known Practical Limitations

- Approximate/numerical tangent behavior can reduce Newton rate from ideal quadratic convergence.
- Convergence sensitivity increases at larger load increments.
- Reproducibility depends on matching both UMAT properties and time/load stepping choices.

## Related Documents

- [USER_GUIDE.md](USER_GUIDE.md)
- [PROJECT_REPORT.pdf](PROJECT_REPORT.pdf)
- [../README.md](../README.md)
