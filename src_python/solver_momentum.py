"""
Predictor step for axial (w) and radial (u) momentum equations.

Uses the CONSERVATION form of the Navier-Stokes equations in
cylindrical coordinates (r, z) — axisymmetric, non-dimensional:

Axial momentum (w at z-faces):
  ∂w/∂t + (1/r)∂(r·u·w)/∂r + ∂(w²)/∂z
      = (1/Re)·[ (1/r)∂/∂r(r·∂w/∂r) + ∂²w/∂z² ]

Radial momentum (u at r-faces):
  ∂u/∂t + (1/r)∂(r·u²)/∂r + ∂(u·w)/∂z
      = (1/Re)·[ (1/r)∂/∂r(r·∂u/∂r) − u/r² + ∂²u/∂z² ]

The pressure gradient is EXCLUDED from the predictor; it is applied
in the corrector step of solver_poisson.py.

Convective terms use FIRST-ORDER UPWIND differencing to ensure
stability at high grid Peclet numbers (Pe_z = |w|·dz·Re ≈ 20).
Diffusive terms use standard second-order central differencing.

Staggered grid layout
---------------------
  w[i, j]  at (r_c[i], z_f[j])   shape (n_r,   n_z+1)
  u[i, j]  at (r_f[i], z_c[j])   shape (n_r+1, n_z)
  r_c[i] = (i+0.5)·dr,   r_f[i] = i·dr
  z_c[j] = (j+0.5)·dz,   z_f[j] = j·dz

FVM control volumes
-------------------
  w-CV: [r_f[i], r_f[i+1]] × [z_c[j-1], z_c[j]]   (for j = 1..n_z-1)
  u-CV: [r_c[i-1], r_c[i]] × [z_f[j], z_f[j+1]]   (for i = 1..n_r-1)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Axial velocity predictor  w*
# ---------------------------------------------------------------------------

def predictor_w(w, u, r_c, r_f, dr, dz, dt, Re):
    """
    Advance axial velocity by one explicit time step (no pressure).

    Convective fluxes use upwind differencing:
      At each face, the transported value of w is taken from the
      upwind side based on the sign of the transporting velocity.

    Parameters
    ----------
    w  : ndarray (n_r, n_z+1)   current axial velocity (BCs already applied)
    u  : ndarray (n_r+1, n_z)   current radial velocity (BCs already applied)
    r_c: ndarray (n_r,)
    r_f: ndarray (n_r+1,)
    dr, dz, dt, Re : float

    Returns
    -------
    w_star : ndarray (n_r, n_z+1)  predicted axial velocity (BCs NOT yet reapplied)
    """
    n_r, n_z1 = w.shape
    n_z = n_z1 - 1

    # ---- ghost extension in r for w ----------------------------------------
    # Wall ghost: no-slip w=0 at r=R → w_ghost = -w[n_r-1, :]
    # Axis ghost: symmetry ∂w/∂r=0  → w_ghost = w[0, :]
    w_ghost_wall = -w[-1:, :]        # shape (1, n_z+1)
    w_ghost_axis =  w[0:1, :]        # shape (1, n_z+1)
    # Extended array:  axis_ghost | w[0..n_r-1] | wall_ghost
    w_ext = np.concatenate([w_ghost_axis, w, w_ghost_wall], axis=0)  # (n_r+2, n_z+1)
    # w_ext[i+1, j] = w[i, j]   (offset of 1)

    # ---- interior slice j=1..n_z-1 ------------------------------------------
    w_int = w[:, 1:-1]   # shape (n_r, n_z-1)

    # ---- u interpolated to z_f[j] (j=1..n_z-1) -----------------------------
    u_at_zf = 0.5 * (u[:, :-1] + u[:, 1:])   # (n_r+1, n_z-1) at (r_f, z_f[1..n_z-1])

    # ---- face velocities and upwind values ----------------------------------

    # --- Radial convective flux  (1/r)∂(r·u·w)/∂r ---
    # Outer r-face (r_f[i+1]):
    u_E = u_at_zf[1:, :]                         # transport vel at outer face, (n_r, n_z-1)
    w_E_upwind = np.where(u_E > 0, w_int, w_ext[2:, 1:-1])  # upwind w

    # Inner r-face (r_f[i]):
    u_I = u_at_zf[:-1, :]                        # transport vel at inner face, (n_r, n_z-1)
    w_I_upwind = np.where(u_I > 0, w_ext[:-2, 1:-1], w_int)  # upwind w

    # --- Axial convective flux  ∂(w²)/∂z ---
    # North face at z_c[j]:  transport velocity = avg of w[i,j] and w[i,j+1]
    w_N = 0.5 * (w_int + w[:, 2:])
    w_N_upwind = np.where(w_N > 0, w_int, w[:, 2:])

    # South face at z_c[j-1]:  transport velocity = avg of w[i,j-1] and w[i,j]
    w_S = 0.5 * (w[:, :-2] + w_int)
    w_S_upwind = np.where(w_S > 0, w[:, :-2], w_int)

    # ---- reshaping for broadcasting (n_r,1) ---------------------------------
    rc  = r_c.reshape(-1, 1)
    rf_o = r_f[1:].reshape(-1, 1)    # r_f[i+1]
    rf_i = r_f[:-1].reshape(-1, 1)   # r_f[i]

    # ---- convective fluxes (upwind) -----------------------------------------
    conv_r = (rf_o * u_E * w_E_upwind - rf_i * u_I * w_I_upwind) / (rc * dr)
    conv_z = (w_N * w_N_upwind - w_S * w_S_upwind) / dz

    # ---- diffusive fluxes (central, 2nd order) ------------------------------
    # Radial:  [r_f[i+1]*(w[i+1]-w[i]) - r_f[i]*(w[i]-w[i-1])] / (r_c*dr^2)
    visc_r = (rf_o * (w_ext[2:, 1:-1] - w_int)
             - rf_i * (w_int - w_ext[:-2, 1:-1])) / (rc * dr**2) / Re

    # Axial:  (w[i,j+1] - 2*w[i,j] + w[i,j-1]) / dz^2
    visc_z = (w[:, 2:] - 2.0 * w_int + w[:, :-2]) / dz**2 / Re

    # ---- assemble predictor -------------------------------------------------
    w_star = w.copy()
    w_star[:, 1:-1] = w_int + dt * (-conv_r - conv_z + visc_r + visc_z)

    # Boundary faces are left as-is (inlet/outlet set by apply_bc_w)
    return w_star


# ---------------------------------------------------------------------------
# Radial velocity predictor  u*
# ---------------------------------------------------------------------------

def predictor_u(u, w, r_c, r_f, dr, dz, dt, Re):
    """
    Advance radial velocity by one explicit time step (no pressure).

    Convective fluxes use upwind differencing.

    Parameters
    ----------
    u  : ndarray (n_r+1, n_z)
    w  : ndarray (n_r,   n_z+1)
    r_c: ndarray (n_r,)
    r_f: ndarray (n_r+1,)
    dr, dz, dt, Re : float

    Returns
    -------
    u_star : ndarray (n_r+1, n_z)
    """
    n_rp1, n_z = u.shape

    # ---- ghost extension in z for u ----------------------------------------
    # Inlet (j=-1 ghost):  u=0 at inlet → u_ghost = -u[:,0] so avg = 0
    # Outlet (j=n_z ghost): zero-gradient → u_ghost = u[:,n_z-1]
    u_ghost_in  = -u[:, 0:1]           # shape (n_r+1, 1)
    u_ghost_out =  u[:, -1:]           # shape (n_r+1, 1)
    u_ext = np.concatenate([u_ghost_in, u, u_ghost_out], axis=1)  # (n_r+1, n_z+2)
    # u_ext[:, j+1] = u[:, j]

    # ---- interior slice i=1..n_r-1 ------------------------------------------
    u_int = u[1:-1, :]    # shape (n_r-1, n_z)

    # ---- face velocities for radial convection ------------------------------
    # Transport velocity at r_c[i] (outer face of u-CV):
    u_out = 0.5 * (u_int + u[2:, :])    # (n_r-1, n_z)
    # Transport velocity at r_c[i-1] (inner face):
    u_in  = 0.5 * (u[:-2, :] + u_int)   # (n_r-1, n_z)

    # Upwind u values for radial convection:
    u_out_upwind = np.where(u_out > 0, u_int, u[2:, :])
    u_in_upwind  = np.where(u_in > 0, u[:-2, :], u_int)

    # ---- face velocities for axial convection --------------------------------
    # w interpolated to r_f[i] (i=1..n_r-1):
    w_at_rf = 0.5 * (w[:-1, :] + w[1:, :])   # (n_r-1, n_z+1) at (r_f[1..n_r-1], z_f)

    # Transport velocity at top face z_f[j+1]:
    w_top = w_at_rf[:, 1:]     # (n_r-1, n_z)
    # Transport velocity at bottom face z_f[j]:
    w_bot = w_at_rf[:, :-1]    # (n_r-1, n_z)

    # u values for upwind in z-direction (using u_ext):
    # South of top face = u at z_c[j] = u_ext[1:-1, 1:-1] = u_int
    # North of top face = u at z_c[j+1] = u_ext[1:-1, 2:]
    u_top_upwind = np.where(w_top > 0, u_ext[1:-1, 1:-1], u_ext[1:-1, 2:])
    # South of bot face = u at z_c[j-1] = u_ext[1:-1, :-2]
    # North of bot face = u at z_c[j] = u_ext[1:-1, 1:-1] = u_int
    u_bot_upwind = np.where(w_bot > 0, u_ext[1:-1, :-2], u_ext[1:-1, 1:-1])

    # ---- reshaping for broadcasting (n_r-1, 1) ------------------------------
    rc_o  = r_c[1:  ].reshape(-1, 1)    # r_c[i],   r-outer face of u-CV
    rc_i  = r_c[:-1 ].reshape(-1, 1)    # r_c[i-1], r-inner face
    rf_int = r_f[1:-1].reshape(-1, 1)   # r_f[i],   where u is stored

    # ---- convective fluxes (upwind) -----------------------------------------
    conv_r = (rc_o * u_out * u_out_upwind - rc_i * u_in * u_in_upwind) / (rf_int * dr)
    conv_z = (w_top * u_top_upwind - w_bot * u_bot_upwind) / dz

    # ---- diffusive fluxes (central, 2nd order) ------------------------------
    # Radial (including −u/r² curvature term):
    visc_r = ((rc_o * (u[2:, :] - u_int) - rc_i * (u_int - u[:-2, :])) / (rf_int * dr**2)
              - u_int / rf_int**2) / Re

    # Axial:
    visc_z = (u_ext[1:-1, 2:] - 2.0 * u_int + u_ext[1:-1, :-2]) / dz**2 / Re

    # ---- assemble predictor -------------------------------------------------
    u_star = u.copy()
    u_star[1:-1, :] = u_int + dt * (-conv_r - conv_z + visc_r + visc_z)

    # Axis and wall faces stay 0 (enforced by apply_bc_u after this call)
    return u_star
