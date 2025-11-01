# Ferrite-Fortran Integration: VEVP Material Model Implementation

This repository demonstrates the integration of advanced Fortran-based UMAT (User MATerial) subroutines with the Ferrite.jl finite element framework in Julia. The implementation focuses on finite strain viscoelastic-viscoplastic (VEVP) material modeling for polymer composites.

## Overview

Material constitutive models are fundamental to finite element analysis of complex material behavior. This project bridges high-performance Fortran material subroutines (commonly used in ABAQUS) with Julia's modern FEM ecosystem through Ferrite.jl, enabling:

- **Performance**: Native Fortran computation speed for material point integration
- **Flexibility**: Julia's high-level abstractions for FEM infrastructure
- **Extensibility**: Easy integration of legacy UMAT subroutines
- **Research**: Advanced material models (VEVP with 8 Maxwell branches)

## Key Features

### Material Model
- **Finite strain VEVP formulation** for RTM6 epoxy resin
- **Viscoelastic response**: 8 Maxwell branches with logarithmic relaxation time distribution (1s - 1000s)
- **Viscoplastic behavior**: Rate-dependent yield with isotropic/kinematic hardening
- **State variables**: 108 internal variables tracking deformation history

### Numerical Implementation
- **Pure displacement formulation** with Lagrange multipliers
- **Newton-Raphson solver** with automatic load stepping
- **Consistent tangent stiffness** from UMAT (DDSDDE matrix)
- **Adaptive convergence** tolerances for large deformation problems

### Validation
- **Cantilever beam test case**: 10% deflection with geometric nonlinearity
- **Stress gradient verification**: Bending-induced tension/compression zones
- **Multi-scale relaxation**: Time-dependent response across 3 orders of magnitude

## Documentation

- **[USER_GUIDE.md](docs/USER_GUIDE.md)** - Installation, running simulations, customization
- **[TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md)** - Mathematical formulation, implementation details

## Quick Start

**📖 Detailed instructions in [USER_GUIDE.md](docs/USER_GUIDE.md)**

```bash
# 1. Install Julia packages
julia -e 'using Pkg; Pkg.add(["Ferrite", "Tensors", "LinearAlgebra", "Printf", "Plots"])'

# 2. Compile Fortran UMAT
cd src/Material_Models && gfortran -shared -fPIC -O2 umat.f -o libumat.so && cd ../..

# 3. Run simulation
julia -e 'include("src/main.jl")'
```

Results saved to `src/POSTPROCESS/plots/` and `src/POSTPROCESS/visualization/`

## Technical Highlights

**📖 Full mathematical details in [TECHNICAL_REFERENCE.md](docs/TECHNICAL_REFERENCE.md)**

### VEVP Material Model
- **8 Maxwell branches**: Multi-scale viscoelastic relaxation (1s - 1000s)
- **Rate-dependent plasticity**: Nonlinear viscoplastic flow
- **108 state variables**: Full deformation history tracking

### Numerical Implementation
- **Finite strain formulation**: Large deformation capability
- **Newton-Raphson solver**: 100 load steps with 1e-4 tolerance
- **UMAT interface**: Fortran-Julia integration via ccall

## Results & Validation

- ✅ **Cantilever beam test**: 10% tip deflection with large deformation
- ✅ **Multi-scale relaxation**: Force decay across 3 orders of magnitude time scales
- ✅ **Stress gradients**: Bending theory validated (tension/compression zones)
- ✅ **Convergence**: 20-30 Newton iterations per step (typical for nonlinear FEM)

## Citation

If you use this code in your research, please cite:

```bibtex
@software{gadwala2025ferrite_vevp,
  author = {Gadwala, Prashanth},
  title = {Ferrite-Fortran Integration: VEVP Material Model Implementation},
  year = {2025},
  url = {https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia}
}
```

## Contributing

Areas for improvement:
- Parallel assembly for large meshes
- Adaptive time stepping
- Additional material models
- Automatic differentiation for exact tangent

## License

MIT License - see LICENSE file for details

## Contact

Prashanth Gadwala - Friedrich-Alexander University of Erlangen-Nuremberg
Project Link: https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia