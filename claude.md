# Claude Code Reference for Incompressible Flow Simulation

## Project Overview
2D laminar flow simulation in a heated cylindrical pipe using finite volume method.
- **Python first**, then port to MATLAB
- **Calculation functions**: Perform all computations and return results
- **Plot functions**: Take computed data, return matplotlib axes (no calculations)

---

## Physical Properties (Dimensional)

### Fluid Properties
| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `rho` | Density | kg/m³ | Constant fluid density (ρ) |
| `mu` | Dynamic viscosity | Pa·s | Constant dynamic viscosity (µ) |
| `lambda_thermal` | Thermal conductivity | W/(m·K) | Constant thermal conductivity (λ) |
| `c_v` | Specific heat capacity | J/(kg·K) | Constant specific heat at constant volume |

### Heat Transfer
| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `q_w` | Wall heat flux | W/m² | Constant wall heat flux |

---

## Geometry (Dimensional)

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `D` | Pipe diameter | m | Diameter of cylindrical pipe |
| `L` | Pipe length | m | Total axial length (L = 50·D) |
| `r` | Radial coordinate | m | Radial position (0 ≤ r ≤ D/2) |
| `z` | Axial coordinate | m | Axial position (0 ≤ z ≤ L) |

---

## Inlet/Boundary Conditions (Dimensional)

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `W_in` | Inlet velocity | m/s | Uniform inlet axial velocity |
| `T_in` | Inlet temperature | K | Uniform inlet temperature |

---

## Flow Variables (Dimensional)

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `w` | Axial velocity | m/s | Axial (z-direction) velocity component |
| `u` | Radial velocity | m/s | Radial (r-direction) velocity component |
| `p` | Pressure | Pa | Static pressure field |
| `T` | Temperature | K | Temperature field |

---

## Non-Dimensional Numbers

| Variable | Name | Formula | Description |
|----------|------|---------|-------------|
| `Re` | Reynolds number | Re = ρ·W_in·D/µ | Ratio of inertial to viscous forces (Re = 100) |
| `Pr` | Prandtl number | Pr = µ·c_v/λ | Ratio of momentum to thermal diffusivity (Pr = 5) |

---

## Non-Dimensional Variables

### Reference Scales
- **Length scale**: D
- **Velocity scale**: W_in
- **Pressure scale**: ρ·W_in²
- **Temperature scale**: q_w/(ρ·c_v·W_in)

### Non-Dimensional Coordinates
| Variable | Name | Formula | Description |
|----------|------|---------|-------------|
| `r_star` | Non-dim radial coord | r* = r/D | Ranges from 0 to 0.5 |
| `z_star` | Non-dim axial coord | z* = z/D | Ranges from 0 to 50 |

### Non-Dimensional Flow Variables
| Variable | Name | Formula | Description |
|----------|------|---------|-------------|
| `w_star` | Non-dim axial velocity | w* = w/W_in | Normalized axial velocity |
| `u_star` | Non-dim radial velocity | u* = u/W_in | Normalized radial velocity |
| `p_star` | Non-dim pressure | p* = p/(ρ·W_in²) | Normalized pressure |
| `theta` | Non-dim temperature | θ = (T - T_in)·ρ·W_in·c_v/q_w | Normalized temperature difference |

---

## Grid Parameters

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `n_r` | Radial grid points | - | Number of grid cells in radial direction |
| `n_z` | Axial grid points | - | Number of grid cells in axial direction |
| `dr` | Radial spacing | m or - | Grid spacing in radial direction (uniform) |
| `dz` | Axial spacing | m or - | Grid spacing in axial direction (uniform) |
| `r_grid` | Radial grid array | m or - | 1D array of radial positions |
| `z_grid` | Axial grid array | m or - | 1D array of axial positions |

### Staggered Grid Notation
- Cell centers: Used for scalar quantities (p, T, grid indices i, j)
- Cell faces: Used for velocity components
  - `w` at axial faces (between cells in z-direction)
  - `u` at radial faces (between cells in r-direction)

---

## Time Parameters

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `t` | Time | s | Physical time |
| `dt` | Time step | s | Time step for explicit time integration |
| `n_steps` | Number of steps | - | Total time steps to simulate |
| `t_end` | End time | s | Final simulation time |

---

## Numerical Parameters

| Variable | Name | Units | Description |
|----------|------|-------|-------------|
| `tol_continuity` | Continuity tolerance | - | Convergence tolerance for continuity equation |
| `tol_velocity` | Velocity tolerance | - | Convergence tolerance for velocity field |
| `tol_temperature` | Temperature tolerance | - | Convergence tolerance for temperature field |
| `max_iter` | Max iterations | - | Maximum iterations for iterative solvers |
| `alpha_p` | Pressure relaxation | - | Under-relaxation factor for pressure correction |

---

## Analytical Solution (Fully Developed Flow)

### Velocity Profile (Fully Developed)
```
w_analytical = 2 * W_in * (1 - (2*r/D)²)
```

### Temperature Profile (Fully Developed)
```
theta_analytical = 4*(z/D) + Re*Pr*[ (1/2)*(2*r/D)² - (1/8)*(2*r/D)⁴ - 7/48 ]
```
NOTE: the original assignment PDF printed this with the radial coefficients doubled
(`(2r/D)² - (1/4)(2r/D)⁴ - 3/4`), i.e. missing the factor 1/2 on Re*Pr. That form is a
transcription error (it gives Nu = 24/11 instead of the classical 48/11 and violates
theta_bulk = 4z*). The corrected form above is the one used everywhere in this repo;
`Assignemnt/Assignment_WS25.pdf` is the corrected transcription, `Assignment_WS25_error.pdf`
the original. See `Report/nusselt_number_analysis.md` for the full derivation.

---

## Code Structure

### Calculation Functions
All calculation functions return computed results and do NOT create plots.

**Examples:**
- `calculate_velocity_field(u, w, p, ...)` → returns updated u, w
- `solve_poisson_equation(rhs, ...)` → returns pressure correction
- `calculate_temperature_field(T, u, w, ...)` → returns updated T
- `check_convergence(residuals, tol)` → returns boolean and residual norms
- `calculate_analytical_solution(r, z, Re, Pr)` → returns w_analytical, theta_analytical

### Plot Functions
All plot functions take pre-computed data and return matplotlib axes. No calculations inside.

**Examples:**
- `plot_velocity_contour(w_star, r_grid, z_grid)` → returns ax
- `plot_temperature_contour(theta, r_grid, z_grid)` → returns ax
- `plot_velocity_profile(w_star, r_grid, z_positions)` → returns ax
- `plot_temperature_profile(theta, r_grid, z_positions)` → returns ax
- `plot_convergence_history(residuals)` → returns ax
- `plot_comparison_with_analytical(w_numerical, w_analytical, r_grid)` → returns ax

---

## File Organization

```
Report_Incropmressible_Flow/
├── claude.md                    # This file
├── Heun/                       # Projected Heun/RK2 implementation
│   ├── src_python/
│   ├── __init__.py
│   ├── constants.py            # Physical constants, Re, Pr
│   ├── grid.py                 # Grid generation
│   ├── initial_conditions.py  # Set initial u, w, p, T
│   ├── boundary_conditions.py # Apply BC at each time step
│   ├── solver_momentum.py     # Solve momentum equations
│   ├── solver_poisson.py      # Poisson solver for pressure
│   ├── solver_temperature.py  # Solve energy equation
│   ├── time_integration.py    # Explicit 2nd order time stepping
│   ├── convergence.py         # Check convergence criteria
│   ├── analytical.py          # Analytical solutions
│   ├── plotting.py            # All plot functions
│   └── main.py                # Main simulation loop
│   ├── src_matlab/
│   ├── Plots_python/
│   └── Plots_matlab/
├── AdamsBashfort2/            # AB2 incremental pressure-correction variant
│   ├── src_python/
│   ├── src_matlab/
│   ├── Plots_python/
│   └── Plots_matlab/
└── Report/                    # LaTeX report files
```

---

## Naming Conventions

1. **Use lowercase with underscores** for Python variables: `w_star`, `n_r`, `Re`
2. **Greek letters**: Write out or use common abbreviations
   - ρ → `rho`
   - µ → `mu`
   - λ → `lambda_thermal` (avoid `lambda` keyword)
   - θ → `theta`
3. **Dimensional vs Non-dimensional**: Append `_star` for non-dimensional: `w` vs `w_star`
4. **Arrays**: Use `_grid` suffix for coordinate arrays: `r_grid`, `z_grid`
5. **Always use consistent terminology**:
   - "axial" (not "longitudinal")
   - "radial" (not "normal")
   - "wall heat flux" (not "heat flow")

---

## Notes

- Viscous dissipation is **neglected** in energy equation
- Axisymmetric flow: no azimuthal (φ) dependence
- Explicit time integration: 2nd order (e.g., RK2 or Adams-Bashforth)
- Spatial discretization: Finite Volume Method, 2nd order
- Solution method: Segregated approach with pressure correction (projection method)
