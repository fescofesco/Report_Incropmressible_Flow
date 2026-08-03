# LaTeX Report - Incompressible Flow Simulation

This directory contains the LaTeX source files for the project report.

## Files

- **main.tex**: Main LaTeX document containing:
  - Mathematical formulation in cylindrical coordinates
  - Non-dimensionalization of governing equations
  - Numerical discretization (Finite Volume Method, 2nd order)
  - Solution procedure (pressure correction, Poisson solver)
  - Analytical solutions for fully developed flow
  - Boundary conditions

## Compiling the Report

### Prerequisites

Install a LaTeX distribution:
- **Windows**: MiKTeX or TeX Live
- **macOS**: MacTeX
- **Linux**: TeX Live

### Compilation

#### Option 1: Using pdflatex (recommended)

```bash
cd Report
pdflatex main.tex
pdflatex main.tex  # Run twice for references
```

#### Option 2: Using latexmk

```bash
cd Report
latexmk -pdf main.tex
```

#### Option 3: Using an IDE

Open `main.tex` in:
- **TeXstudio**
- **Overleaf** (upload all files)
- **VS Code** with LaTeX Workshop extension

## Structure

The report follows the guidelines specified in `Guidel_lines_for_report.pdf`:

1. **Title Page**: Name, course number, date
2. **Assignment**: Copy of the original assignment (included via `\includepdf`)
3. **Mathematical Formulation**:
   - Governing equations in cylindrical coordinates
   - Non-dimensionalization
   - Non-dimensional parameters (Re, Pr)
4. **Numerical Solution**:
   - Computational domain
   - Staggered grid structure
   - Spatial discretization (FVM, 2nd order)
   - Temporal discretization (explicit, 2nd order)
   - Boundary conditions
   - Pressure correction method
   - Poisson solver
5. **Analytical Solutions**:
   - Fully developed velocity profile
   - Fully developed temperature profile
6. **Results** (to be added after simulation):
   - Velocity and temperature contours
   - Profiles at selected positions
   - Comparison with analytical solutions
   - Convergence history

## Figures

The geometry diagram is automatically included from:
```
../Plots_python/geometry_diagram.png
```

Additional figures (simulation results) will be added to `../Plots_python/` after running the simulations.

## Key Equations

### Non-Dimensional Numbers

- Reynolds number: Re = ρW<sub>in</sub>D/μ = 100
- Prandtl number: Pr = μc<sub>v</sub>/λ = 5

### Analytical Solutions (Fully Developed)

**Velocity:**
```
w/W_in = 2[1 - (2r/D)²]
```

**Temperature:**
```
θ = 4(z/D) + Re·Pr[½(2r/D)² - ⅛(2r/D)⁴ - 7/48]
```

where θ = (T - T<sub>in</sub>)ρW<sub>in</sub>c<sub>v</sub>/q<sub>w</sub>

Note: the assignment sheet states this with the radial-term coefficients doubled
((2r/D)² - ¼(2r/D)⁴ - ¾). The halved coefficients above are used in the code, based
on independent re-derivation and empirical validation against the CFD solution —
see `temperature_formula_review.md`.

## Notes

- All variable names follow the conventions in [../claude.md](../claude.md)
- Consistent terminology is used throughout (e.g., "axial", "radial")
- The report can be written in English or German (currently in English)
