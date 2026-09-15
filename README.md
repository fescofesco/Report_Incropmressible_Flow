# Incompressible Flow Simulation - Cylindrical Pipe

2D laminar flow simulation in a heated cylindrical pipe with constant wall heat flux.

## Assignment

- **Course**: Numerical Simulation and Modelling of Incompressible Flow (LV-Nr.: 321.053 EX)
- **Semester**: WS 2025/26
- **Topic**: Development of flow at the entrance of a cylindrical pipe

## Project Structure

```
Report_Incropmressible_Flow/
├── Assignemnt/                  # Assignment PDFs
│   ├── Assignment_WS25.pdf
│   └── Guidel_lines_for_report.pdf
├── claude.md                    # IMPORTANT: Variable definitions and conventions
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies
├── Heun/                        # Projected Heun/RK2 implementation and plots
│   ├── src_python/
│   ├── src_matlab/
│   ├── Plots_python/
│   └── Plots_matlab/
├── AdamsBashfort2/              # AB2 incremental pressure-correction variant
│   ├── src_python/
│   ├── src_matlab/
│   ├── Plots_python/
│   └── Plots_matlab/
└── Report/                      # LaTeX report files
```

## Setup Instructions

### 1. Create Conda Environment

Open **Anaconda Prompt** or **Miniconda Prompt** and run:

```bash
cd C:\Users\Admin\Documents\Repos\Report_Incropmressible_Flow
conda env create -f environment.yml
```

Or manually create the environment:

```bash
conda create -n flow python=3.11 numpy matplotlib scipy pandas jupyter -y
```

### 2. Activate Environment

```bash
conda activate flow
```

### 3. Verify Installation

```bash
python -c "import numpy; import matplotlib; import scipy; print('All packages installed successfully!')"
```

## Important Files

### **claude.md** - MUST READ!
This file defines:
- All variable names and their meanings
- Naming conventions (always use these!)
- Code structure (calculation functions vs. plot functions)
- Non-dimensional variables
- File organization

**Always refer to claude.md when writing code to maintain consistency!**

## Simulation Parameters

- **Reynolds number**: Re = 100
- **Prandtl number**: Pr = 5
- **Pipe geometry**: L/D = 50
- **Boundary conditions**:
  - Uniform inlet: w|z=0 = W_in, T|z=0 = T_in
  - Constant wall heat flux: q_w = const

## Analytical Solutions (Fully Developed)

**Velocity profile:**
```
w/W_in = 2[1 - (2r/D)²]
```

**Temperature profile:**
```
θ = 4(z/D) + Re·Pr[½(2r/D)² - ⅛(2r/D)⁴ - 7/48]
```

where θ = (T - T_in)·ρ·W_in·c_v/q_w

This is the approved form (direct integration of the non-dimensional energy
equation; constant −7/48 from the global energy balance, so θ_bulk = 4z* for all
z ≥ 0; gives Nu = 48/11). An earlier transcription of the assignment omitted the
½ factor on the radial terms and used −¾; that form is not used anywhere in the
code. Full derivation and sources: Appendix B of `Report/main.pdf`.

## Code Philosophy

### Calculation Functions
- Perform all computations
- Return numerical results
- No plotting inside
- Example: `calculate_velocity_field()`, `solve_poisson_equation()`

### Plot Functions
- Take pre-computed data as input
- Return matplotlib axes
- No calculations inside
- Example: `plot_velocity_contour()`, `plot_temperature_profile()`

## Next Steps

1. ✅ Project structure created
2. ✅ Variable definitions in claude.md
3. ✅ Constants and plotting modules created
4. ⏳ Implement grid generation
5. ⏳ Implement finite volume discretization
6. ⏳ Implement time integration (2nd order explicit)
7. ⏳ Implement Poisson solver for pressure correction
8. ⏳ Implement main simulation loop
9. ⏳ Run simulation and validate against analytical solutions
10. ⏳ Generate plots for report
11. ⏳ Port to MATLAB

## Development Notes

- **Python first**, then port to MATLAB
- Follow naming conventions in [claude.md](claude.md)
- Keep calculation and plotting functions separate
- Use consistent terminology (see claude.md)
- Document all functions with docstrings

## Contact

Felix Scope
felix_scope@yahoo.de
