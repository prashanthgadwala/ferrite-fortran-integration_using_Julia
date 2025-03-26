# ferrite-fortran-integration_using_Julia

This repository contains the implementation for a student thesis/project exploring the integration of advanced Fortran-based material models into the Ferrite.jl finite element (FE) toolbox in Julia. The project bridges the robust material modeling capabilities of Fortran with the modern features of Julia.

## Overview

Material models are essential for simulating how materials react under loading and other conditions in finite element analysis. By integrating Fortran subroutines with Julia's Ferrite.jl library, this project aims to enhance the capability of FE simulations while demonstrating an efficient workflow for combining traditional and modern programming approaches.

### Objectives
1. Understand Julia and Ferrite.jl.
2. Set up a simple load case for an FE mesh using Ferrite.jl's native material model.
3. Extend the implementation to call a Fortran-based material model for more advanced load cases.
4. Implement and test the `UMAT` subroutine for finite strain viscoelastic-viscoplastic (VEVP) material modeling.
5. Document the entire process and results.

## Repository Structure

- `src/`: Source code for Julia and Fortran implementations.
  - `preprocessing/`: Directory for preprocessing tasks and material model implementations.
    - `umat.f90`: Fortran UMAT subroutine for finite strain VEVP material modeling.
  - `processing/`: Directory for processing tasks.
    - `fe_simulation.jl`: FE simulation setup.
  - `postprocessing/`: Directory for postprocessing tasks.
  - `main.jl`: Main script for running simulations.
  - `fortran_wrapper.jl`: Julia-Fortran integration code.
- `test/`: Unit tests and example load case results.
- `docs/`: Documentation and results, including plots or images from simulations.
- `examples/`: Sample scripts to demonstrate basic and advanced use cases.

## Progress

### UMAT Subroutine
- Implemented the `UMAT` subroutine in Fortran for finite strain viscoelastic-viscoplastic (VEVP) material modeling.
- The subroutine includes:
  - A Newton-Raphson solver for nonlinear equations.
  - Numerical tangent computation.
  - Support for 8 Maxwell viscoelastic branches.
  - Hardcoded tolerances and maximum iterations for the solver.
- Helper subroutines for matrix operations (e.g., determinant, inverse, dot products) and material property calculations have been added.
- The `ABA_PARAM.INC` file has been created to define parameters such as the number of stress components, state variables, and material properties.

### Integration with Julia
- The Fortran `UMAT` subroutine will be compiled into a shared library (`.so`) for integration with Julia.
- A wrapper script (`fortran_wrapper.jl`) is being developed to call the Fortran subroutine from Julia.
- Initial tests with simple load cases are planned to validate the integration.

### Compilation
- The `UMAT` subroutine can be compiled using the following command:
  ```bash
  gfortran -o umat.so -shared -fPIC src/umat.f90


## Getting Started

### Prerequisites
- Julia (>= 1.8)
- Ferrite.jl library
- Fortran compiler (e.g., `gfortran`)
- BLAS/LAPACK libraries (if needed for Fortran integration)

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia.git
   cd ferrite-fortran-integration_using_Julia
2. Download and install Julia
   '''bash
   curl -fsSL https://install.julialang.org | sh
3. install Julia dependencies:
   ```bash
   import Pkg; Pkg.add("Ferrite")
3. Compile Fortran material model subroutine:
   ```bash
   gfortran -o material_model.so -shared -fPIC src/PREPROCESSING/Material_Models/umat.f90

### Usage
1. Run a simple load case:
   ```bash
   julia src/main.jl
2. Modify and test more advanced load cases in the examples/ directory.

### Contribution

Contributions are welcome! Please create a fork, make changes, and open a pull request. For major changes, discuss them first via an issue.

### License

This project is licensed under the MIT License.