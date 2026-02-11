"""
Physical constants and simulation parameters
All variable names follow the conventions defined in claude.md
"""

# ============================================================================
# NON-DIMENSIONAL NUMBERS
# ============================================================================
Re = 100.0  # Reynolds number: ρ·W_in·D/µ
Pr = 5.0    # Prandtl number: µ·c_v/λ

# ============================================================================
# GEOMETRY (Non-dimensional)
# ============================================================================
L_over_D = 50.0  # Pipe length in diameters (L/D = 50)

# ============================================================================
# GRID PARAMETERS
# ============================================================================
n_r = 50   # Number of grid points in radial direction
n_z = 250  # Number of grid points in axial direction

# ============================================================================
# TIME PARAMETERS
# ============================================================================
dt = 0.001      # Time step (non-dimensional)
n_steps = 10000 # Number of time steps
output_interval = 100  # Save output every N steps

# ============================================================================
# NUMERICAL PARAMETERS
# ============================================================================
tol_continuity = 1e-6    # Convergence tolerance for continuity
tol_velocity = 1e-6      # Convergence tolerance for velocity
tol_temperature = 1e-6   # Convergence tolerance for temperature
max_iter = 1000          # Maximum iterations for Poisson solver
alpha_p = 0.7            # Pressure under-relaxation factor

# ============================================================================
# PHYSICAL PROPERTIES (Dimensional - for reference only)
# ============================================================================
# These are example values; the actual simulation uses non-dimensional equations
rho = 1000.0           # Density [kg/m³] (e.g., water at 20°C)
mu = 0.001             # Dynamic viscosity [Pa·s]
lambda_thermal = 0.6   # Thermal conductivity [W/(m·K)]
c_v = 4180.0           # Specific heat capacity [J/(kg·K)]
D = 0.01               # Pipe diameter [m] (1 cm)
W_in = 0.1             # Inlet velocity [m/s]
T_in = 293.15          # Inlet temperature [K] (20°C)
q_w = 1000.0           # Wall heat flux [W/m²]

# Computed dimensional parameters
L = L_over_D * D       # Pipe length [m]
