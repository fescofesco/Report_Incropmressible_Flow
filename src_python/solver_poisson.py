"""
Poisson solver for the pressure-correction equation and velocity/pressure update.

Projection method steps (non-dimensional, ρ*=1):
  Step 2: Solve  ∇²p' = (1/dt) · ∇·u*
  Step 3: Correct velocities
              u^{n+1}_{i+½,j} = u*_{i+½,j} − dt·(p'_{i+1,j} − p'_{i,j})/dr
              w^{n+1}_{i,j+½} = w*_{i,j+½} − dt·(p'_{i,j} − p'_{i,j-1})/dz
              (note: p' is at cell centres i,j; w face spans cells j-1..j)
  Step 4: Update pressure
              p^{n+1}_{i,j} = p^n_{i,j} + ζ·p'_{i,j}

Poisson equation in cylindrical coordinates (FVM, cell (i,j)):
  [r_{i+½}(p'_{i+1,j}−p'_{i,j}) − r_{i-½}(p'_{i,j}−p'_{i-1,j})] / (r_c[i]·dr²)
  + [p'_{i,j+1} − 2p'_{i,j} + p'_{i,j-1}] / dz²
  = b_{i,j}

where  b_{i,j} = (1/dt)·[(r_{i+½}·u*_{i+1,j} − r_{i-½}·u*_{i,j}) / (r_c[i]·dr)
                         + (w*_{i,j+1} − w*_{i,j}) / dz]

Boundary conditions for p':
  Centreline (i=0):      ∂p'/∂r = 0  → r_{-½}=0 term vanishes naturally
  Wall       (i=n_r-1):  ∂p'/∂r = 0  → outer flux = 0  (Neumann)
  Inlet      (j=0):      ∂p'/∂z = 0  → south flux = 0  (Neumann)
  Outlet     (j=n_z-1):  p' = 0      → Dirichlet (reference pressure)

The matrix is assembled ONCE and reused every time step.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ---------------------------------------------------------------------------
# Matrix assembly  (call once at start-up)
# ---------------------------------------------------------------------------

def build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z):
    """
    Assemble the sparse Poisson matrix A (size n_r*n_z × n_r*n_z).

    Flat index: k = i*n_z + j

    Parameters
    ----------
    r_c : (n_r,)   cell-centre radial coords
    r_f : (n_r+1,) r-face radial coords
    dr, dz : float
    n_r, n_z : int

    Returns
    -------
    A : scipy sparse CSR matrix  (n_r*n_z × n_r*n_z)
    """
    N = n_r * n_z
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(r)
        cols.append(c)
        vals.append(v)

    for i in range(n_r):
        for j in range(n_z):
            k = i * n_z + j

            # ---- Dirichlet at outlet (j = n_z-1): p' = 0 -------------------
            if j == n_z - 1:
                add(k, k, 1.0)
                continue

            # ---- coefficients from cylindrical Poisson ----------------------
            rf_o = r_f[i + 1]   # r_{i+½}
            rf_i = r_f[i]       # r_{i-½}
            rc   = r_c[i]

            coeff_E = rf_o / (rc * dr**2)       # p'[i+1, j]
            coeff_W = rf_i / (rc * dr**2)       # p'[i-1, j]
            coeff_N = 1.0 / dz**2               # p'[i, j+1]
            coeff_S = 1.0 / dz**2               # p'[i, j-1]

            # Neumann wall (i = n_r-1): outer flux = 0 → drop E coeff
            if i == n_r - 1:
                coeff_E = 0.0

            # Neumann inlet (j = 0): south flux = 0 → drop S coeff
            if j == 0:
                coeff_S = 0.0

            # Centre coefficient (sum of neighbour coeffs with negative sign)
            coeff_C = -(coeff_E + coeff_W + coeff_N + coeff_S)

            add(k, k, coeff_C)

            # East (i+1, j)
            if i < n_r - 1:
                add(k, (i + 1) * n_z + j, coeff_E)

            # West (i-1, j)  — for i=0: r_f[0]=0 → coeff_W=0, but add anyway
            if i > 0:
                add(k, (i - 1) * n_z + j, coeff_W)

            # North (i, j+1)
            if j < n_z - 1:
                add(k, i * n_z + (j + 1), coeff_N)

            # South (i, j-1)
            if j > 0:
                add(k, i * n_z + (j - 1), coeff_S)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    return A


# ---------------------------------------------------------------------------
# Solve the Poisson equation
# ---------------------------------------------------------------------------

def solve_poisson(A, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z):
    """
    Solve  A · p'_flat = b  for pressure correction p'.

    Parameters
    ----------
    A       : scipy sparse CSR  (n_r*n_z, n_r*n_z)
    u_star  : ndarray (n_r+1, n_z)  predicted radial velocity
    w_star  : ndarray (n_r, n_z+1)  predicted axial velocity
    r_c, r_f: arrays
    dr, dz, dt : float
    n_r, n_z   : int

    Returns
    -------
    p_prime : ndarray (n_r, n_z)
    """
    # RHS: b[i,j] = (1/dt) * div(u*)
    #   div(u*) = [r_{i+½}·u*[i+1,j] − r_{i-½}·u*[i,j]] / (r_c[i]·dr)
    #           + [w*[i,j+1] − w*[i,j]] / dz

    # Radial divergence
    rf_o = r_f[1:].reshape(-1, 1)    # (n_r, 1)
    rf_i = r_f[:-1].reshape(-1, 1)   # (n_r, 1)
    rc   = r_c.reshape(-1, 1)        # (n_r, 1)

    div_u = (rf_o * u_star[1:, :] - rf_i * u_star[:-1, :]) / (rc * dr)  # (n_r, n_z)
    div_w = (w_star[:, 1:] - w_star[:, :-1]) / dz                        # (n_r, n_z)

    b = (div_u + div_w) / dt   # (n_r, n_z)

    # Outlet row (j=n_z-1) has Dirichlet p'=0 → RHS = 0
    b[:, -1] = 0.0

    b_flat = b.ravel()

    # Direct solve
    p_prime_flat = spla.spsolve(A, b_flat)
    p_prime = p_prime_flat.reshape(n_r, n_z)

    return p_prime


# ---------------------------------------------------------------------------
# Velocity correction  (Step 3)
# ---------------------------------------------------------------------------

def correct_velocity(u_star, w_star, p_prime, r_f, dr, dz, dt, n_r, n_z):
    """
    Apply pressure-gradient correction to obtain divergence-free velocities.

    u^{n+1}[i, j] = u*[i, j] − dt · (p'[i, j] − p'[i-1, j]) / dr
                                        for i = 1..n_r-1
    w^{n+1}[i, j] = w*[i, j] − dt · (p'[i, j] − p'[i, j-1]) / dz
                                        for j = 1..n_z-1

    Parameters
    ----------
    u_star  : ndarray (n_r+1, n_z)
    w_star  : ndarray (n_r, n_z+1)
    p_prime : ndarray (n_r, n_z)
    dr, dz, dt : float
    n_r, n_z   : int

    Returns
    -------
    u_new : ndarray (n_r+1, n_z)
    w_new : ndarray (n_r, n_z+1)
    """
    u_new = u_star.copy()
    w_new = w_star.copy()

    # Radial velocity at interior r-faces (i=1..n_r-1):
    #   u[i,j] at r_f[i]; pressure gradient (p'[i,j] - p'[i-1,j]) / dr
    u_new[1:-1, :] = u_star[1:-1, :] - dt * (p_prime[1:, :] - p_prime[:-1, :]) / dr

    # Axial velocity at interior z-faces (j=1..n_z-1):
    #   w[i,j] at z_f[j]; pressure gradient (p'[i,j] - p'[i,j-1]) / dz
    #   p'[i,j] is the cell NORTH of face j  (cell index j, since z_f[j] borders cells j-1 and j)
    w_new[:, 1:-1] = w_star[:, 1:-1] - dt * (p_prime[:, 1:] - p_prime[:, :-1]) / dz

    return u_new, w_new


# ---------------------------------------------------------------------------
# Pressure update  (Step 4)
# ---------------------------------------------------------------------------

def update_pressure(p, p_prime, zeta):
    """
    p^{n+1} = p^n + ζ · p'

    Parameters
    ----------
    p       : ndarray (n_r, n_z)
    p_prime : ndarray (n_r, n_z)
    zeta    : float   under-relaxation factor (e.g. 0.5)

    Returns
    -------
    p_new : ndarray (n_r, n_z)
    """
    return p + zeta * p_prime
