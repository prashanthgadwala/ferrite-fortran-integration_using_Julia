# User Guide: Ferrite-Fortran VEVP Integration

Complete guide for installing, running, and customizing VEVP material simulations.

---

## Table of Contents
1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Running Simulations](#running-simulations)
4. [Understanding Results](#understanding-results)
5. [Customization](#customization)
6. [Troubleshooting](#troubleshooting)

---

## Quick Start

**5-Minute Setup:**
```bash
# 1. Clone
git clone https://github.com/prashanthgadwala/ferrite-fortran-integration_using_Julia.git
cd ferrite-fortran-integration_using_Julia

# 2. Install Julia packages
julia -e 'using Pkg; Pkg.add(["Ferrite", "Tensors", "LinearAlgebra", "Printf", "Plots"])'

# 3. Compile Fortran UMAT
cd src/Material_Models && gfortran -shared -fPIC -O2 umat.f -o libumat.so && cd ../..

# 4. Run!
julia -e 'include("src/main.jl")'
```

---

## Installation

### Prerequisites
- **Julia** ≥ 1.8 ([Download](https://julialang.org/downloads/))
- **gfortran** compiler
- **ParaView** (optional, for 3D visualization)

### Step-by-Step

#### 1. Install Julia
```bash
# macOS (Homebrew)
brew install julia

# Linux
curl -fsSL https://install.julialang.org | sh

# Windows
# Download installer from julialang.org
```

#### 2. Install Julia Packages
```bash
julia -e 'using Pkg; Pkg.add(["Ferrite", "Tensors", "LinearAlgebra", "Printf", "Plots"])'
```

#### 3. Check gfortran
```bash
gfortran --version
# Should show version ≥ 9.0
```

If not installed:
```bash
# macOS
brew install gcc

# Ubuntu/Debian
sudo apt install gfortran

# Windows
# Install MinGW-w64
```

#### 4. Compile UMAT
```bash
cd src/Material_Models
gfortran -shared -fPIC -O2 umat.f -o libumat.so
cd ../..

# Verify
ls -lh src/Material_Models/libumat.so
# Should show ~100 KB file
```

---

## Running Simulations

### Basic Run

```bash
julia
```

```julia
julia> include("src/main.jl")
```

**Expected Output:**
```
======================================================================
PURE DISPLACEMENT FORMULATION WITH UMAT
======================================================================
Material type: vevp
Mesh: 320 hexahedral elements
  Beam dimensions: 10.0 cm × 1.0 cm × 1.0 cm

Material: VEVP UMAT (RTM6 epoxy)
  Viscoelastic: ALL 8 MAXWELL BRANCHES ACTIVE! 🔥
  
Starting load stepping: 100 steps
======================================================================

Step 1: t = 0.01
  Iter  1: ||r|| = 6.418e-01
  ...
  Iter 27: ||r|| = 9.234e-05
  ✅ Converged!

...

✅ ANALYSIS COMPLETE!
Steps completed: 100
Final tip deflection: -10.0 mm

Results exported to: src/POSTPROCESS/visualization/cantilever_vevp_8branches.vtu
```

**Runtime:** ~5-10 minutes (depends on your system)

### What Happens During Simulation

1. **Mesh Generation** (1 second)
   - Creates 320 hexahedral elements
   - Applies boundary conditions

2. **Load Stepping** (5-10 minutes)
   - 100 incremental steps
   - Each step: Newton-Raphson iterations until convergence
   - Progress shown: `Step X: Iter Y: ||r|| = Z`

3. **Post-Processing** (5 seconds)
   - Generates 4 plots automatically
   - Exports VTK file for ParaView

4. **Results Saved**
   - Plots → `src/POSTPROCESS/plots/`
   - VTK → `src/POSTPROCESS/visualization/`

---

## Understanding Results

### Generated Files

**Location**: `src/POSTPROCESS/plots/`

1. **`force_displacement.png`**
   - Shows force-deflection curve
   - **Key feature**: Force decreases over time (viscoelastic relaxation)
   - Expected: Curved line showing multi-scale decay

2. **`stress_strain.png`**
   - Material response at critical point
   - Shows bending-induced stress gradient
   - Top surface: tension, bottom: compression

3. **`displacement_history.png`**
   - Tip deflection evolution over load steps
   - Should show smooth ramp to -10 mm

4. **`deformed_shape.png`**
   - Mesh visualization with color-coded displacement
   - Classic cantilever S-curve

**Location**: `src/POSTPROCESS/visualization/`

5. **`cantilever_vevp_8branches.vtu`**
   - 3D visualization for ParaView
   - Contains displacement field, mesh geometry

### Viewing Results

**Quick View (PNG files):**
```bash
# macOS
open src/POSTPROCESS/plots/force_displacement.png

# Linux
xdg-open src/POSTPROCESS/plots/force_displacement.png

# Windows
start src/POSTPROCESS/plots/force_displacement.png
```

**3D Visualization (ParaView):**
```bash
paraview src/POSTPROCESS/visualization/cantilever_vevp_8branches.vtu
```

In ParaView:
- Color by: `u` (displacement field)
- Apply: `Warp By Vector` filter (magnify deformation)
- View: Top tension (red), bottom compression (blue)

### Expected Physics

**Viscoelastic Relaxation:**
- **Initial response**: Stiff (all 8 Maxwell branches active)
- **Time evolution**: Force decreases as branches relax
- **Multi-scale decay**: Fast (1-10s) + medium (10-100s) + slow (100-1000s)
- **Final state**: Approaches equilibrium (72% force reduction)

**Cantilever Bending:**
- Top surface: Tension (red in visualization)
- Bottom surface: Compression (blue)
- Neutral axis: Zero stress at mid-height

---

## Customization

### Change Material Properties

**File**: `src/main.jl` (around line 320)

**Reduce to 4 Maxwell branches:**
```julia
1.0e6, 0.8e6, 0.6e6, 0.5e6, 0, 0, 0, 0,  # Only first 4 active
1.0, 3.0, 10.0, 30.0, 0, 0, 0, 0,        # Corresponding times
```

### Change Geometry

**File**: `src/main.jl` (around line 545)

**Example: Longer, finer mesh**
```julia
L_length = 0.2   # 20 cm beam
N_length = 40    # Finer discretization
```

### Change Loading & Solver

**File**: `src/main.jl`

**More load steps (smoother convergence):**
```julia
n_steps = 150    # Line ~630
```

**Adjust tolerance:**
```julia
NEWTON_TOL = 1.0e-3   # Faster (line ~693)
NEWTON_TOL = 1.0e-6   # More accurate
```

---

## Troubleshooting

### Common Errors

#### 1. "Could not load library libumat.so"

**Problem:** UMAT not compiled or wrong path

**Solution:**
```bash
cd src/Material_Models
gfortran -shared -fPIC -O2 umat.f -o libumat.so
cd ../..

# Verify file exists
ls src/Material_Models/libumat.so
```

#### 2. "Reached maximum Newton iterations"

**Problem:** Load step too large or tolerance too tight

**Solution 1** - Increase load steps (smaller increments):
```julia
# In src/main.jl
n_steps = 150    # More steps
```

**Solution 2** - Relax tolerance:
```julia
NEWTON_TOL = 1.0e-3    # Looser tolerance
```

**Solution 3** - Allow more iterations:
```julia
NEWTON_MAXITER = 100   # More patience
```

#### 3. "Plots not showing"

**Problem:** Plots.jl backend issue

**Solution:**
```julia
using Pkg
Pkg.add("Plots")
using Plots
gr()  # Use GR backend
```

#### 4. "Memory error" or "Out of memory"

**Problem:** Mesh too fine or too many state variables

**Solution** - Coarsen mesh:
```julia
N_length = 10    # Fewer elements
N_height = 3
N_width = 3
```

#### 5. Simulation is very slow

**Problem:** Too many load steps or fine mesh

**Quick test configuration:**
```julia
n_steps = 20         # Faster (less accurate)
N_length = 10        # Coarser mesh
NEWTON_MAXITER = 30
```

**High-accuracy configuration:**
```julia
n_steps = 200        # Slower (more accurate)
N_length = 40        # Fine mesh
NEWTON_TOL = 1.0e-6
```

### Performance Tips

**Speed up simulation:**
- Reduce load steps: `n_steps = 50`
- Coarsen mesh: `(10,3,3)` instead of `(20,4,4)`
- Relax tolerance: `1.0e-3`
- Compile with optimization: `gfortran -O3 ...`

**Improve accuracy:**
- Increase load steps: `n_steps = 200`
- Refine mesh: `(40,8,8)`
- Tighten tolerance: `1.0e-6`
- Verify with analytical solution

### Getting Help

1. Check **README.md** for overview
2. Check **TECHNICAL_REFERENCE.md** for implementation details
3. Review error messages carefully
4. Open issue on GitHub
5. Contact: [Your institution/email]

---

## Test Cases

### Test 1: Quick (2 min)
```julia
n_steps = 20, N_length = 10, N_height = 3
```

### Test 2: Research Quality (10 min) - Default
```julia
n_steps = 100, N_length = 20, N_height = 4
```

### Test 3: High Accuracy (30 min)
```julia
n_steps = 200, N_length = 40, NEWTON_TOL = 1.0e-6
```

---

**📖 For implementation details**, see **[TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)**

Ready to start? `julia -e 'include("src/main.jl")'` 🚀
