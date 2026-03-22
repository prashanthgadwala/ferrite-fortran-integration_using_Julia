# User Guide: Ferrite-Fortran VEVP Integration

Practical guide for installing, running, and modifying the unit-cube uniaxial VEVP benchmark.

## Quick Start

```bash
git clone https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia.git
cd ferrite-fortran-integration_using_Julia

julia -e 'using Pkg; Pkg.add(["Ferrite", "Tensors", "LinearAlgebra", "Printf", "Plots", "WriteVTK"])'

cd src/Material_Models
gfortran -shared -fPIC -O2 umat.f -o libumat.so
cd ../..

julia -e 'include("src/main.jl")'
```

## Prerequisites

- Julia 1.8+
- gfortran compiler
- ParaView (optional, for VTK visualization)

Install checks:

```bash
julia --version
gfortran --version
```

## What the Current Code Runs

The default setup in `src/main.jl` is:

- Geometry: unit cube (`1 x 1 x 1`)
- Mesh: `2 x 2 x 2` Q2 hexahedral elements
- Loading: uniaxial tension in x-direction up to `u_x = 0.3`
- Material model: VEVP UMAT, 8 Maxwell branches, 108 state variables
- Runtime controls: `n_steps = 100`, `DTIME = 0.1`

The project report documents a validation configuration with `50` increments and physical time scaling (`dt = 100 s`, total `5000 s`).

## Running the Simulation

```bash
julia -e 'include("src/main.jl")'
```

Expected console markers:

- `PURE DISPLACEMENT FORMULATION WITH UMAT`
- `Mesh: 8 hexahedral elements (refined cube)`
- `Material: VEVP UMAT`
- `Starting load stepping:`
- `ANALYSIS COMPLETE!`

## Output Files

Main output directory:

- `src/POSTPROCESS/plots/`

Generated figures:

- `force_displacement.png`: viscoelastic strain evolution (`E_ve` components)
- `stress_strain.png`: viscoplastic internal variable evolution (`F_vp`, `gamma`)
- `displacement_history.png`: prescribed displacement history
- `deformed_shape.png`: deformed configuration view

ParaView time series:

- `src/POSTPROCESS/plots/vtk_timesteps/results.pvd`
- `src/POSTPROCESS/plots/vtk_timesteps/step_*.vtu`

## Open Results Quickly

macOS:

```bash
open src/POSTPROCESS/plots/deformed_shape.png
open src/POSTPROCESS/plots/vtk_timesteps/results.pvd
```

Linux:

```bash
xdg-open src/POSTPROCESS/plots/deformed_shape.png
```

## Customization

All key runtime settings are in `src/main.jl` inside `solve()`.

Common edits:

- Change load steps: `n_steps = ...`
- Change time increment: `DTIME = ...`
- Change mesh density: `N_length`, `N_width`, `N_height`
- Change loading amplitude: `u_load = ...`

Material constants are defined in `get_vevp_properties()` in `src/main.jl`.

## Troubleshooting

1. `could not load library libumat.so`

```bash
cd src/Material_Models
gfortran -shared -fPIC -O2 umat.f -o libumat.so
cd ../..
```

2. Newton iterations fail to converge

- Increase `n_steps`
- Reduce load amplitude (`u_load`)
- Check UMAT build and path

3. No plots or VTK output

- Confirm `create_plots(...)` is executed (end of `src/main.jl`)
- Ensure `WriteVTK` is installed

## Recommended Reading

- [../README.md](../README.md)
- [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)
- [PROJECT_REPORT.pdf](PROJECT_REPORT.pdf)
