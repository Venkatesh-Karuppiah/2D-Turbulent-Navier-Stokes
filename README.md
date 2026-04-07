# Turbulent Navier-Stokes Solver

A Fortran 90 finite-volume solver for the Reynolds-Averaged Navier-Stokes (RANS) equations on structured quadrilateral grids, using Jameson's 4-stage Runge-Kutta explicit time-marching scheme with artificial dissipation and an algebraic Baldwin-Lomax turbulence model.

## Solvers

| Binary | Description |
|--------|-------------|
| **`twodturb`** | 2D Cartesian compressible turbulent flow solver |
| **`axiturb`** | Axisymmetric compressible turbulent flow solver |

Both solvers share a common codebase (~97% identical physics) through a modular architecture. The `IAXIS` flag in `flow.in` selects the geometric formulation at runtime.

## Quick Start

```bash
# Build both solvers
make all

# Run a case (from the case directory)
cd examples/case1
../../bin/twodturb     # 2D Cartesian
../../bin/axiturb      # Axisymmetric
```

## Project Structure

```
├── Makefile                          # Build system
├── src/
│   ├── solvers/
│   │   ├── constants.f90             # Physical & numerical constants (PARAMETER)
│   │   ├── solver_state.f90          # Shared arrays, DISSIPATION, POSTPROCESS
│   │   ├── geometry.f90              # Mesh I/O, surface vectors, volumes (IAXIS-aware)
│   │   ├── viscous.f90               # Viscous flux computation (IAXIS-aware)
│   │   ├── source_terms.f90          # Turbulence source terms (IAXIS-aware)
│   │   ├── boundary.f90              # Ghost cell boundary conditions (IAXIS-aware)
│   │   ├── residual.f90              # Inviscid flux residual & RK update
│   │   ├── twodturb.f90              # 2D Cartesian solver (main program)
│   │   └── axiturb.f90              # Axisymmetric solver (main program)
│   └── gridgen/
│       ├── generateInputGeometryFiles.f90
│       ├── jpl.f90
│       ├── jplwithprop.f90
│       ├── nasab2.f90
│       └── ramp.f90
├── examples/
│   └── case1/
│       ├── grid.in                   # Mesh: node coordinates + element connectivity
│       ├── bc.in                     # Boundary: ghost cell definitions
│       └── flow.in                   # Flow: physical parameters + solver settings
└── utils/
    └── displayMesh.py                # Mesh visualization utility
```

## Input Files

### `flow.in` — Case Configuration

Single line with 17 space/tab-separated values followed by a header line:

```
grid.in  bc.in  1.4  28.96  1  100000  0  0.25  0.72  0.9  873500.0  298.0  -1.0e5  0.01  0.00001  1  0
FNAMEGRID  FNAMEBC  GAMA  AMOL  LOCAL  NITER  ISTART  CFL  PRANDTL  PRANDTLT  P0  T0  PAMB  AMACHINIT  RTRED  NPSEUDO  IAXIS
```

| Parameter | Description |
|-----------|-------------|
| `FNAMEGRID` | Grid input filename |
| `FNAMEBC` | Boundary condition input filename |
| `GAMA` | Ratio of specific heats (γ) |
| `AMOL` | Molecular weight [kg/kmol] |
| `LOCAL` | Time stepping: 0 = global, 1 = local |
| `NITER` | Number of iterations |
| `ISTART` | 0 = cold start, 1 = restart from `restart.in` |
| `CFL` | CFL number |
| `PRANDTL` | Laminar Prandtl number |
| `PRANDTLT` | Turbulent Prandtl number |
| `P0` | Stagnation pressure [Pa] |
| `T0` | Stagnation temperature [K] |
| `PAMB` | Ambient pressure [Pa] (negative = extrapolated) |
| `AMACHINIT` | Initial Mach number |
| `RTRED` | Tolerance reduction factor |
| `NPSEUDO` | Pseudo time exponent for turbulence source |
| **`IAXIS`** | **0 = 2D Cartesian, 1 = Axisymmetric** |

### `grid.in` — Mesh Definition

```
  <NNODES>  <NELEMENTS>
  <NODE_ID>  <X[mm]>  <Y[mm]>
  ...
  <ELEM_ID>  <N1>  <N2>  <N3>  <N4>  <NC1>  <NC2>  <NC3>  <NC4>
  ...
```

Coordinates are specified in **millimeters** and converted to meters internally. `NC1..NC4` are the neighbor element indices for each face (-1 if boundary).

### `bc.in` — Boundary Conditions

```
  <NGHOSTS>
  <PARENT>  <GHOST>  <TYPE>  <SIDE>
  ...
```

| Type | Boundary Condition |
|------|-------------------|
| 2 | Subsonic inlet (isentropic) |
| 3 | Farfield (Riemann extrapolation) |
| 4 | Symmetry / slip wall |
| 5 | Wall / axis (no-slip reflection) |
| 6 | Combustion / mass injection (with G-threshold correction) |
| 7 | Boundary layer inlet profile (1/7th power law) |

## Numerical Method

### Governing Equations
Compressible RANS in integral form on unstructured quadrilateral grids, solved in conservative variables:

```
Q = [ρ, ρu, ρv, ρE, ρν̃]ᵀ
```

### Spatial Discretization
- **Finite volume** on vertex-centered dual cells
- Central differencing for inviscid fluxes (2nd order)
- Viscous fluxes via thin-layer Navier-Stokes approximation with gradient reconstruction at face midpoints

### Artificial Dissipation
Jameson's blended 2nd/4th-order scheme:
- 4th-order background dissipation (smooth regions)
- 2nd-order near shocks/discontinuities (pressure-switch sensor)
- Separate treatment for turbulence transport equation

### Time Integration
Jameson's 4-stage Runge-Kutta with coefficients:
- α₁ = 0.25, α₂ = 1/3, α₃ = 0.5, α₄ = 1.0
- Local or global time stepping
- Residual smoothing not implemented

### Turbulence Model
Algebraic Baldwin-Lomax model with:
- Sutherland's law for molecular viscosity: μ = 1.1848×10⁻⁷ √(M) T⁰·⁶
- tanh-based damping functions
- Turbulent production via velocity magnitude gradients

### Geometry: 2D vs Axisymmetric

The key differences controlled by `IAXIS`:

| Quantity | 2D Cartesian | Axisymmetric |
|----------|-------------|--------------|
| Surface vector | S = (Δy, -Δx) | S = (Δy·ȳ, -Δx·ȳ) |
| Cell volume | Area | Area × 2π·ȳ |
| Velocity divergence | ∂u/∂x + ∂v/∂y | ∂u/∂x + ∂v/∂y + v/y |
| Turbulence production | Standard | + radial strain terms |
| Radial source | None | S₃ = (P - τᵧᵧ)/y |
| Combustion GPORT | GPORT/DHYD | 2·GPORT/DHYD² |

## Output Files

| File | Description |
|------|-------------|
| `tec_plot.dat` | Tecplot FE format: nodal X, Y, U, V, P[bar], ρ, M, T |
| `restart.in` | Final solution + solver state for restart |
| `restartNNNNNN.in` | Intermediate restart files every 20,000 iterations |
| Convergence printed to stdout every 1,000 iterations |

## Grid Generation

The `src/gridgen/` directory contains standalone programs for generating orthogonal body-fitted grids. Each reads a `.neu` (GAMBIT-neutral) file and produces `grid.in` + `bc.in`:

```bash
cd src/gridgen
gfortran -O3 -o generateInputGeometryFiles generateInputGeometryFiles.f90
./generateInputGeometryFiles   # reads tape7.neu, produces grid.in + bc.in
```

## Build

```bash
make all          # Build twodturb and axiturb
make clean        # Remove build artifacts
```

Compiler defaults to `gfortran` with `-O3 -fdefault-real-8 -Wall`. Adjust `FC` and `FFLAGS` in the Makefile as needed.

## Architecture Decisions

### Module Dependency Graph

```
constants        →  PARAMETER values only (no mutable state)
    ↑
solver_state     →  All mutable arrays + DISSIPATION + POSTPROCESS
    ↑
geometry  viscous  source_terms  boundary  residual  →  Physics modules with IAXIS branching
    ↑
twodturb  axiturb  →  Main programs (~230 lines each)
```

### Why This Structure

The original codebase contained two ~2800-line files that were **97% identical**. The refactoring extracts:

1. **`constants`** — Immutable physical and numerical parameters, eliminating magic numbers scattered through the code
2. **`solver_state`** — Single source of truth for all field arrays, plus the two subroutines that are byte-for-byte identical (dissipation, postprocessing)
3. **Physics modules** — Each subroutine with 2D/axisymmetric differences gets its own module with an `IAXIS` flag. The differences are additive corrections (extra terms in divergence, source terms, surface weighting) rather than completely different algorithms
4. **Thin main programs** — Each solver is now just the RK4 loop, I/O, and initialization, delegating all physics to shared modules

This means a bug fix in dissipation, boundary conditions, viscous fluxes, or the RK4 scheme automatically applies to both solvers.
