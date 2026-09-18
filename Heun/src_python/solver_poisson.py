"""
Poisson solver for the pressure-correction equation and velocity/pressure update.

Projection method steps (non-dimensional, ρ*=1):
  Step 2: Solve  ∇²p' = (1/dt) · ∇·u*
  Step 3: Correct velocities
              u^{n+1}_{i+½,j} = u*_{i+½,j} − dt·(p'_{i+1,j} − p'_{i,j})/dr
              w^{n+1}_{i,j+½} = w*_{i,j+½} − dt·(p'_{i,j} − p'_{i,j-1})/dz
              (note: p' is at cell centres i,j; w face spans cells j-1..j)
  Step 4: Store a relaxed estimate of the non-incremental projection pressure
              p^{n+1}_{i,j} = (1-ζ)p^n_{i,j} + ζ·p'_{i,j}

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

The linear system is solved with a hand-written DIRECT method, as the
assignment requires: the operator is block tridiagonal (n_z blocks of size
n_r x n_r), and it is eliminated by the block Thomas algorithm on top of a
Doolittle LU kernel written out below. Because the operator is constant in
time it is factorised ONCE (factorize_poisson) before the time loop; each
time step then costs only the two sweeps in solve_poisson. No library linear
solve (scipy spsolve, numpy.linalg.solve, MATLAB backslash) is used.

build_poisson_matrix assembles the same operator as a sparse matrix. It is no
longer used by the time loop and is kept only so that verify_poisson_solver.py
can check the hand-written solver's residual against an independent reference.
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

    Not used by the solver itself -- see the module docstring.
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
# Dense LU factorisation with partial pivoting  (hand-written)
# ---------------------------------------------------------------------------
#
# The assignment requires a DIRECT solution method that is programmed, not a
# library call (no numpy.linalg.solve / scipy spsolve / MATLAB backslash).
# These two routines are the elimination kernel: classical Doolittle LU,
# P·A = L·U, written out explicitly. They operate on the small (n_r x n_r)
# blocks of the block-tridiagonal Poisson operator and are called ONLY during
# set-up (see factorize_poisson).
# ---------------------------------------------------------------------------

def lu_factor(A):
    """
    Doolittle LU factorisation with partial pivoting:  P·A = L·U.

    Parameters
    ----------
    A : ndarray (n, n)   dense matrix (not modified)

    Returns
    -------
    LU  : ndarray (n, n)
        Packed factors. The strict lower triangle holds L (its unit diagonal
        is implied), the upper triangle including the diagonal holds U.
    piv : ndarray (n,) int
        piv[k] is the row that was interchanged with row k at step k.
    """
    n = A.shape[0]
    LU = np.array(A, dtype=float)
    piv = np.zeros(n, dtype=int)

    for k in range(n):
        # --- partial pivoting: largest magnitude in column k, rows k..n-1
        m = k + int(np.argmax(np.abs(LU[k:, k])))
        piv[k] = m
        if LU[m, k] == 0.0:
            raise np.linalg.LinAlgError("zero pivot in column %d" % k)
        if m != k:
            LU[[k, m], :] = LU[[m, k], :]

        # --- elimination; the multipliers are stored in the created zeros
        LU[k + 1:, k] /= LU[k, k]
        LU[k + 1:, k + 1:] -= np.outer(LU[k + 1:, k], LU[k, k + 1:])

    return LU, piv


def lu_solve(LU, piv, B):
    """
    Solve A·X = B from the factors returned by lu_factor: forward substitution
    with L, then back substitution with U.

    Parameters
    ----------
    LU, piv : output of lu_factor
    B       : ndarray (n,) or (n, m)   one or several right-hand sides

    Returns
    -------
    X : ndarray, same shape as B
    """
    n = LU.shape[0]
    shape_in = B.shape
    X = np.array(B, dtype=float).reshape(n, -1)

    # --- apply the recorded row interchanges:  b <- P·b
    for k in range(n):
        m = piv[k]
        if m != k:
            X[[k, m], :] = X[[m, k], :]

    # --- forward substitution  L·Y = P·b   (L has unit diagonal)
    for k in range(n - 1):
        X[k + 1:, :] -= np.outer(LU[k + 1:, k], X[k, :])

    # --- back substitution  U·X = Y
    for k in range(n - 1, -1, -1):
        X[k, :] /= LU[k, k]
        X[:k, :] -= np.outer(LU[:k, k], X[k, :])

    return X.reshape(shape_in)


# ---------------------------------------------------------------------------
# Block-tridiagonal (block-Thomas) factorisation of the Poisson operator
# ---------------------------------------------------------------------------

def factorize_poisson(r_c, r_f, dr, dz, n_r, n_z):
    """
    Factorise the pressure-Poisson operator ONCE with a hand-written direct
    method: the block Thomas algorithm built on the dense LU kernel above.

    Structure
    ---------
    Group the unknowns p'[:, j] by axial station j. The 5-point stencil
    couples a cell only to its radial neighbours (same j) and to j ± 1, so the
    operator is BLOCK TRIDIAGONAL with n_z blocks of size n_r x n_r:

        S_j·p_{j-1} + D_j·p_j + N_j·p_{j+1} = b_j ,     j = 0 .. n_z-1

        D_j = tridiag( w_i , -(e_i + w_i + 2·a_z) , e_i )     (interior j)
        S_j = N_j = a_z·I ,     a_z = 1/dz²
        e_i = r_{i+½}/(r_i·dr²) ,    w_i = r_{i-½}/(r_i·dr²)

    The boundary conditions enter as: e_{n_r-1} = 0 (Neumann at the wall);
    w_0 = 0 automatically because r_{-½} = 0 (Neumann at the centreline); no
    south coupling at j = 0, so that diagonal carries a_z once instead of
    twice (Neumann at the inlet); and D = I, S = N = 0 at j = n_z-1
    (Dirichlet p' = 0 at the outlet).

    Block Thomas forward elimination
    --------------------------------
        M_0 = D_0
        M_j = D_j − S_j·M_{j-1}^{-1}·N_{j-1} = D_j − a_z²·M_{j-1}^{-1}

    which reduces the solve to the two sweeps in solve_poisson,
    y_j = M_j^{-1}(b_j − a_z·y_{j-1})  and  p_j = y_j − a_z·M_j^{-1}·p_{j+1}.

    The operator is constant in time, so the n_z Schur complements are
    inverted here, once. Every subsequent Poisson solve then costs only two
    (n_r x n_r) matrix-vector products per axial station, i.e. O(n_r²·n_z).

    Stability
    ---------
    −A is an irreducibly diagonally dominant M-matrix: |diagonal| equals the
    sum of the off-diagonal magnitudes in every interior row and is strictly
    greater in the Dirichlet outlet rows. Every Schur complement M_j inherits
    that property, so no pivoting is needed between blocks; partial pivoting
    is used inside lu_factor regardless.

    Returns
    -------
    fac : dict with keys
        'Minv' : ndarray (n_z, n_r, n_r)   inverted Schur complements M_j^{-1}
        'az'   : float                     a_z = 1/dz²
        'n_r', 'n_z' : int
    """
    az = 1.0 / dz**2

    # --- radial stencil coefficients (independent of j)
    e = (r_f[1:] / (r_c * dr**2)).copy()    # coupling to i+1
    w = (r_f[:-1] / (r_c * dr**2)).copy()   # coupling to i-1; w[0]=0 as r_f[0]=0
    e[-1] = 0.0                             # Neumann at the wall

    def tridiag(main_diag):
        M = np.zeros((n_r, n_r))
        idx = np.arange(n_r)
        M[idx, idx] = main_diag
        M[idx[:-1], idx[1:]] = e[:-1]       # super-diagonal
        M[idx[1:], idx[:-1]] = w[1:]        # sub-diagonal
        return M

    D_inlet = tridiag(-(e + w + az))        # j = 0, no south flux
    D_int = tridiag(-(e + w + 2.0 * az))    # 1 <= j <= n_z-2

    I_n = np.eye(n_r)
    Minv = np.zeros((n_z, n_r, n_r))

    # j = 0: no south coupling, so M_0 = D_0
    LU, piv = lu_factor(D_inlet)
    Minv[0] = lu_solve(LU, piv, I_n)

    # j = 1 .. n_z-2: Schur complement update, then invert
    for j in range(1, n_z - 1):
        LU, piv = lu_factor(D_int - az**2 * Minv[j - 1])
        Minv[j] = lu_solve(LU, piv, I_n)

    # j = n_z-1: Dirichlet row, D = I and no south coupling
    Minv[n_z - 1] = I_n

    return {'Minv': Minv, 'az': az, 'n_r': n_r, 'n_z': n_z}


# ---------------------------------------------------------------------------
# Solve the Poisson equation
# ---------------------------------------------------------------------------

def solve_poisson(fac, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z):
    """
    Solve  ∇²p' = (1/dt)·∇·u*  with the block-Thomas sweeps prepared by
    factorize_poisson. No library linear solve is involved.

    Parameters
    ----------
    fac : dict
        Output of factorize_poisson (built once, before the time loop).
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

    # Outlet column (j=n_z-1) has Dirichlet p'=0 → RHS = 0
    b[:, -1] = 0.0

    # Column j of b IS block j of the block-tridiagonal system, so no
    # flattening or flat-index bookkeeping is needed anywhere.
    Minv = fac['Minv']
    az = fac['az']

    # --- forward sweep:  y_j = M_j^{-1}·(b_j − a_z·y_{j-1})
    y = np.empty((n_r, n_z))
    y[:, 0] = Minv[0] @ b[:, 0]
    for j in range(1, n_z - 1):
        y[:, j] = Minv[j] @ (b[:, j] - az * y[:, j - 1])
    y[:, n_z - 1] = b[:, n_z - 1]     # Dirichlet row: M = I, no south coupling

    # --- back substitution:  p_j = y_j − a_z·M_j^{-1}·p_{j+1}
    p_prime = np.empty((n_r, n_z))
    p_prime[:, n_z - 1] = y[:, n_z - 1]
    for j in range(n_z - 2, -1, -1):
        p_prime[:, j] = y[:, j] - az * (Minv[j] @ p_prime[:, j + 1])

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
    p^{n+1} = p^n + ζ · (p' − p^n) = (1−ζ)·p^n + ζ·p'

    NOTE: this is a *relaxed blend towards p'*, not an accumulation of p'
    onto p^n. Because the momentum predictor here excludes the pressure
    gradient entirely (a non-incremental / Chorin-type projection, not the
    incremental-pressure-correction scheme), p' as computed by solve_poisson
    already IS the full physical pressure for this step (dimensionally
    consistent with du/dt = -grad(p) + ...), not a small correction that
    should accumulate. Naively doing p^{n+1} = p^n + ζ·p' (additive) makes
    p grow roughly linearly in time even after the flow reaches steady
    state, since p' does not vanish at steady state under this splitting
    (unlike the incremental scheme, where the predictor already carries the
    previous pressure gradient and p' -> 0 at convergence).

    Parameters
    ----------
    p       : ndarray (n_r, n_z)
    p_prime : ndarray (n_r, n_z)
    zeta    : float   relaxation factor towards the freshly solved p' (e.g. 0.7)

    Returns
    -------
    p_new : ndarray (n_r, n_z)
    """
    return p + zeta * (p_prime - p)
