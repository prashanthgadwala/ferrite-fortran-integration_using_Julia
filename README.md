# Ferrite-Fortran Integration: VEVP Material Model Implementation

A reproducible scientific workflow for coupling legacy Fortran UMAT constitutive models with the Julia-based Ferrite.jl finite element framework.

## Abstract

This repository demonstrates a practical Julia-Fortran bridge for nonlinear finite-strain material simulation. A Fortran UMAT (ABAQUS-style interface) is compiled as a shared library and called from Julia via `ccall`, allowing reuse of legacy constitutive code while keeping finite element assembly, solver control, and postprocessing in Ferrite.jl.

The current benchmark is a unit-cube uniaxial tension test using a finite-strain viscoelastic-viscoplastic (VEVP) epoxy model with 8 Maxwell branches and 108 state variables per integration point.

## Why This Repository

- Reuses validated legacy UMAT implementations without full re-coding in Julia.
- Keeps constitutive updates compiled (Fortran), while preserving Julia-level workflow flexibility.
- Provides an end-to-end benchmark setup suitable for research reporting and reproducibility.

## Technical Highlights

### Constitutive Model (VEVP)

- Finite-strain viscoelastic-viscoplastic coupling.
- 8 Maxwell branches spanning relaxation times from 1 s to 1000 s.
- 108 internal state variables per integration point for history-dependent response.

### Julia-Fortran Integration

- Fortran UMAT compiled as a shared library (`libumat.so`).
- Called from Julia with `ccall` using in-place memory updates.
- ABAQUS-compatible stress/tangent/state variable interface (`STRESS`, `DDSDDE`, `STATEV`).

### Nonlinear FEM Workflow

- Ferrite.jl driver with Newton-Raphson iterations and line-search stabilization.
- Benchmark mesh: Q2 hexahedral unit cube (125 nodes, 375 displacement DOFs).
- Loading protocol: 50 increments, dt = 100 s, total time = 5000 s.

### Validation Snapshot

- Ferrite.jl and ABAQUS run with the same UMAT and loading definition.
- Normalized stress-shape agreement: mean error 3.71%, maximum error 5.38%.
- Typical Newton convergence for the VEVP case: 15-25 iterations per increment.

## Documentation

- [docs/USER_GUIDE.md](docs/USER_GUIDE.md): installation, execution, customization
- [docs/TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md): formulation and implementation details
- [docs/thesis.tex](docs/thesis.tex): most up-to-date report-level scientific documentation

## Getting Started (First-Time Setup)

Use this section if you are on a new machine.

### 1. Install Required Software

You need all of the following:

- Git
- Julia (recommended: 1.10 or newer)
- GFortran compiler

#### macOS

1. Install Homebrew (if needed):

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

2. Install dependencies:

```bash
brew install git julia gcc
```

Notes:
- `gcc` includes `gfortran`.
- On Apple Silicon, Homebrew binaries are typically in `/opt/homebrew/bin`.

#### Ubuntu/Debian

```bash
sudo apt update
sudo apt install -y git julia gfortran
```

If the distro Julia is too old, install from https://julialang.org/downloads/.

#### Windows

1. Install Git: https://git-scm.com/download/win
2. Install Julia: https://julialang.org/downloads/
3. Install MSYS2: https://www.msys2.org/
4. In MSYS2 shell, install gfortran:

```bash
pacman -S --needed mingw-w64-x86_64-gcc-fortran
```

5. Ensure `git`, `julia`, and `gfortran` are on your PATH.

### 2. Verify Installation

```bash
git --version
julia --version
gfortran --version
```

### 3. Clone Repository

```bash
git clone https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia.git
cd ferrite-fortran-integration_using_Julia
```

### 4. Install Julia Packages

```bash
julia -e 'using Pkg; Pkg.add(["Ferrite", "Tensors", "LinearAlgebra", "Printf", "Plots"])'
```

### 5. Compile UMAT Shared Library

```bash
cd src/Material_Models
gfortran -shared -fPIC -O2 umat.f -o libumat.so
cd ../..
```

Windows fallback:

```bash
gfortran -shared -O2 umat.f -o umat.dll
```

### 6. Run Simulation

```bash
julia -e 'include("src/main.jl")'
```

### 7. Check Outputs

Primary outputs are written under:

- `src/POSTPROCESS/plots/`
- `src/POSTPROCESS/visualization/`

## Quick Start (Configured Machine)

```bash
cd src/Material_Models && gfortran -shared -fPIC -O2 umat.f -o libumat.so && cd ../..
julia -e 'include("src/main.jl")'
```

## Reproducibility Checklist

For report-quality reproducibility, record:

- OS and architecture
- Julia version
- GFortran version
- Git commit hash
- Any modified material parameters (`PROPS` / UMAT inputs)

## Repository Layout

```text
src/
  main.jl                      # FEM driver and solver loop
  Material_Models/
    umat.f                     # Fortran UMAT implementation
    ABA_PARAM.INC              # UMAT include definitions
  POSTPROCESS/
    postprocess_results.jl     # plotting/postprocessing script

docs/
  USER_GUIDE.md
  TECHNICAL_REFERENCE.md
  thesis.tex                   # latest report-grade reference

test/
  square.f
  square.jl
```

## Common Issues

- `julia: command not found`
  Install Julia and restart terminal, or fix PATH.

- `gfortran: command not found`
  Install GCC/GFortran and restart terminal.

- Shared library compile errors
  Build from `src/Material_Models/` and verify compiler availability.

- First run is slow
  Expected behavior due to Julia package precompilation.

## Contributing

Contributions are welcome, especially in:

- solver scalability and parallel assembly
- adaptive time stepping
- additional constitutive models
- exact/analytic consistent tangent implementations

Please keep changes reproducible and document parameter/configuration updates.

## Citation

If you use this code in research, please cite:

```bibtex
@software{gadwala2025ferrite_vevp,
  author = {Gadwala, Prashanth},
  title = {Ferrite-Fortran Integration: VEVP Material Model Implementation},
  year = {2025},
  url = {https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Contact

Prashanth Gadwala  
Friedrich-Alexander University of Erlangen-Nuremberg  
Project: https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia
