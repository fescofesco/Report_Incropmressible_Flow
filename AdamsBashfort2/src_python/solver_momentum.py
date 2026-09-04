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

Convective terms use a 2ND-ORDER UPWIND (linear-extrapolation / LUD) scheme
— see numerics.upwind2_face — which retains the upwind stability bias at
high grid Peclet numbers (Pe_z = |w|·dz·Re ≈ 20) while being formally
2nd-order accurate, matching the diffusive terms (standard 2nd-order
central differencing).

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

from numerics import upwind2_face


# ---------------------------------------------------------------------------
# Axial velocity predictor  w*
# ---------------------------------------------------------------------------

def compute_rhs_w(w, u, r_c, r_f, dr, dz, Re):
    """
    Evaluate the explicit RHS of the axial-momentum equation (no pressure,
    no time discretization):  RHS = -conv_r - conv_z + visc_r + visc_z

    Kept separate from time-stepping so the same evaluation can be reused
    as the current and stored previous RHS in the AB2 scheme.

    Convective fluxes use upwind differencing:
      At each face, the transported value of w is taken from the
      upwind side based on the sign of the transporting velocity.

    Parameters
    ----------
    w  : ndarray (n_r, n_z+1)   current axial velocity (BCs already applied)
    u  : ndarray (n_r+1, n_z)   current radial velocity (BCs already applied)
    r_c: ndarray (n_r,)
    r_f: ndarray (n_r+1,)
    dr, dz, Re : float

    Returns
    -------
    rhs : ndarray (n_r, n_z-1)  interior-face RHS (excludes j=0 and j=n_z faces)
    """
    n_r, n_z1 = w.shape
    n_z = n_z1 - 1

    # ---- ghost extension in r for w (depth 1, used by diffusion) -----------
    # Wall ghost: no-slip w=0 at r=R → w_ghost = -w[n_r-1, :]
    # Axis ghost: symmetry ∂w/∂r=0  → w_ghost = w[0, :]
    w_ghost_wall = -w[-1:, :]        # shape (1, n_z+1)
    w_ghost_axis =  w[0:1, :]        # shape (1, n_z+1)
    # Extended array:  axis_ghost | w[0..n_r-1] | wall_ghost
    w_ext = np.concatenate([w_ghost_axis, w, w_ghost_wall], axis=0)  # (n_r+2, n_z+1)
    # w_ext[i+1, j] = w[i, j]   (offset of 1)

    # ---- second ghost layer in r, needed by the 2nd-order-upwind stencil ----
    # Axis:  mirror one cell further in.  Wall: continue the antisymmetric
    # (odd) reflection about the no-slip wall face one cell further out.
    w_ghost_axis2 = w[1:2, :]
    w_ghost_wall2 = -w[-2:-1, :]
    w_ext2 = np.concatenate(
        [w_ghost_axis2, w_ghost_axis, w, w_ghost_wall, w_ghost_wall2], axis=0)  # (n_r+4, n_z+1)
    # w_ext2[i+2, j] = w[i, j]

    # ---- interior slice j=1..n_z-1 ------------------------------------------
    w_int = w[:, 1:-1]   # shape (n_r, n_z-1)

    # ---- u interpolated to z_f[j] (j=1..n_z-1) -----------------------------
    u_at_zf = 0.5 * (u[:, :-1] + u[:, 1:])   # (n_r+1, n_z-1) at (r_f, z_f[1..n_z-1])

    # --- Radial convective flux  (1/r)∂(r·u·w)/∂r ---
    w_r_im2 = w_ext2[0:n_r,     1:-1]
    w_r_im1 = w_ext2[1:n_r + 1, 1:-1]
    w_r_i   = w_ext2[2:n_r + 2, 1:-1]     # = w_int
    w_r_ip1 = w_ext2[3:n_r + 3, 1:-1]
    w_r_ip2 = w_ext2[4:n_r + 4, 1:-1]

    # Outer r-face (r_f[i+1]):
    u_E = u_at_zf[1:, :]                         # transport vel at outer face, (n_r, n_z-1)
    w_E_upwind = upwind2_face(w_r_im1, w_r_i, w_r_ip1, w_r_ip2, u_E)

    # Inner r-face (r_f[i]):
    u_I = u_at_zf[:-1, :]                        # transport vel at inner face, (n_r, n_z-1)
    w_I_upwind = upwind2_face(w_r_im2, w_r_im1, w_r_i, w_r_ip1, u_I)

    # --- Axial convective flux  ∂(w²)/∂z ---
    # w already has real values exactly AT the inlet/outlet faces (Dirichlet /
    # zero-gradient), so only ONE extra ghost point beyond each end is needed
    # for the 2nd-order-upwind stencil (constant extension: repeat the
    # boundary value, consistent with uniform inflow / zero-gradient outflow).
    w_zghost_before = w[:, 0:1]
    w_zghost_after  = w[:, -1:]
    w_zext1 = np.concatenate([w_zghost_before, w, w_zghost_after], axis=1)  # (n_r, n_z+3)
    # w_zext1[:, j+1] = w[:, j]

    w_z_jm2 = w_zext1[:, 0:n_z - 1]
    w_z_jm1 = w_zext1[:, 1:n_z]
    w_z_j   = w_zext1[:, 2:n_z + 1]      # = w_int
    w_z_jp1 = w_zext1[:, 3:n_z + 2]
    w_z_jp2 = w_zext1[:, 4:n_z + 3]

    # North face at z_c[j]:  transport velocity = avg of w[i,j] and w[i,j+1]
    w_N = 0.5 * (w_z_j + w_z_jp1)
    w_N_upwind = upwind2_face(w_z_jm1, w_z_j, w_z_jp1, w_z_jp2, w_N)

    # South face at z_c[j-1]:  transport velocity = avg of w[i,j-1] and w[i,j]
    w_S = 0.5 * (w_z_jm1 + w_z_j)
    w_S_upwind = upwind2_face(w_z_jm2, w_z_jm1, w_z_j, w_z_jp1, w_S)

    # ---- reshaping for broadcasting (n_r,1) ---------------------------------
    rc  = r_c.reshape(-1, 1)
    rf_o = r_f[1:].reshape(-1, 1)    # r_f[i+1]
    rf_i = r_f[:-1].reshape(-1, 1)   # r_f[i]

    # ---- convective fluxes (2nd-order upwind) --------------------------------
    conv_r = (rf_o * u_E * w_E_upwind - rf_i * u_I * w_I_upwind) / (rc * dr)
    conv_z = (w_N * w_N_upwind - w_S * w_S_upwind) / dz

    # ---- diffusive fluxes (central, 2nd order) ------------------------------
    # Radial:  [r_f[i+1]*(w[i+1]-w[i]) - r_f[i]*(w[i]-w[i-1])] / (r_c*dr^2)
    visc_r = (rf_o * (w_ext[2:, 1:-1] - w_int)
             - rf_i * (w_int - w_ext[:-2, 1:-1])) / (rc * dr**2) / Re

    # Axial:  (w[i,j+1] - 2*w[i,j] + w[i,j-1]) / dz^2
    visc_z = (w[:, 2:] - 2.0 * w_int + w[:, :-2]) / dz**2 / Re

    return -conv_r - conv_z + visc_r + visc_z


def predictor_w(w, u, r_c, r_f, dr, dz, dt, Re):
    """
    Advance axial velocity by one explicit Euler step (no pressure).
    Thin wrapper around compute_rhs_w; the main time loop uses
    compute_rhs_w directly to build the second-order AB2 update.

    Returns
    -------
    w_star : ndarray (n_r, n_z+1)  predicted axial velocity (BCs NOT yet reapplied)
    """
    rhs = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re)
    w_star = w.copy()
    w_star[:, 1:-1] = w[:, 1:-1] + dt * rhs
    # Boundary faces are left as-is (inlet/outlet set by apply_bc_w)
    return w_star


# ---------------------------------------------------------------------------
# Radial velocity predictor  u*
# ---------------------------------------------------------------------------

def compute_rhs_u(u, w, r_c, r_f, dr, dz, Re):
    """
    Evaluate the explicit RHS of the radial-momentum equation (no pressure,
    no time discretization). Reused as the current and previous AB2 RHS.

    Convective fluxes use upwind differencing.

    Parameters
    ----------
    u  : ndarray (n_r+1, n_z)
    w  : ndarray (n_r,   n_z+1)
    r_c: ndarray (n_r,)
    r_f: ndarray (n_r+1,)
    dr, dz, Re : float

    Returns
    -------
    rhs : ndarray (n_r-1, n_z)  interior-face RHS (excludes i=0 and i=n_r faces)
    """
    n_rp1, n_z = u.shape
    n_r = n_rp1 - 1

    # ---- ghost extension in z for u (depth 2), needed by diffusion and by ---
    # the 2nd-order-upwind stencil.
    # Inlet:  u=0 at inlet face → odd reflection, continued one cell further
    # Outlet: zero-gradient, continued one cell further
    u_ghost_in1  = -u[:, 0:1]
    u_ghost_in2  = -u[:, 1:2]
    u_ghost_out1 =  u[:, -1:]
    u_ghost_out2 =  u[:, -1:]
    u_ext = np.concatenate([u_ghost_in1, u, u_ghost_out1], axis=1)  # (n_r+1, n_z+2), depth 1 (diffusion)
    u_zext2 = np.concatenate(
        [u_ghost_in2, u_ghost_in1, u, u_ghost_out1, u_ghost_out2], axis=1)  # (n_r+1, n_z+4)
    # u_zext2[:, j+2] = u[:, j]

    # ---- second ghost layer in r, needed by the 2nd-order-upwind stencil ----
    # u=0 exactly at the axis (i=0) and wall (i=n_r) grid points already, so
    # the ghost one point beyond each is a simple linear extrapolation
    # (equivalently, odd reflection about the zero boundary value).
    u_ghost_axis_ext  = -u[1:2, :]      # represents "u[-1]"
    u_ghost_wall_ext  = -u[-2:-1, :]    # represents "u[n_r+1]"
    u_ext_r = np.concatenate([u_ghost_axis_ext, u, u_ghost_wall_ext], axis=0)  # (n_r+3, n_z)
    # u_ext_r[i+1, :] = u[i, :]

    # ---- interior slice i=1..n_r-1 ------------------------------------------
    u_int = u[1:-1, :]    # shape (n_r-1, n_z)

    # --- Radial convective flux  (1/r)∂(r·u²)/∂r ---
    u_r_im2 = u_ext_r[0:n_r - 1, :]
    u_r_im1 = u_ext_r[1:n_r,     :]
    u_r_i   = u_ext_r[2:n_r + 1, :]     # = u_int
    u_r_ip1 = u_ext_r[3:n_r + 2, :]
    u_r_ip2 = u_ext_r[4:n_r + 3, :]

    # Transport velocity at r_c[i] (outer face of u-CV):
    u_out = 0.5 * (u_int + u[2:, :])    # (n_r-1, n_z)
    # Transport velocity at r_c[i-1] (inner face):
    u_in  = 0.5 * (u[:-2, :] + u_int)   # (n_r-1, n_z)

    u_out_upwind = upwind2_face(u_r_im1, u_r_i, u_r_ip1, u_r_ip2, u_out)
    u_in_upwind  = upwind2_face(u_r_im2, u_r_im1, u_r_i, u_r_ip1, u_in)

    # ---- face velocities for axial convection --------------------------------
    # w interpolated to r_f[i] (i=1..n_r-1):
    w_at_rf = 0.5 * (w[:-1, :] + w[1:, :])   # (n_r-1, n_z+1) at (r_f[1..n_r-1], z_f)

    # Transport velocity at top face z_f[j+1]:
    w_top = w_at_rf[:, 1:]     # (n_r-1, n_z)
    # Transport velocity at bottom face z_f[j]:
    w_bot = w_at_rf[:, :-1]    # (n_r-1, n_z)

    # --- Axial convective flux  ∂(u·w)/∂z --- (2nd-order upwind, interior rows)
    u_z_jm2 = u_zext2[1:-1, 0:n_z]
    u_z_jm1 = u_zext2[1:-1, 1:n_z + 1]
    u_z_j   = u_zext2[1:-1, 2:n_z + 2]      # = u_int
    u_z_jp1 = u_zext2[1:-1, 3:n_z + 3]
    u_z_jp2 = u_zext2[1:-1, 4:n_z + 4]

    u_top_upwind = upwind2_face(u_z_jm1, u_z_j, u_z_jp1, u_z_jp2, w_top)
    u_bot_upwind = upwind2_face(u_z_jm2, u_z_jm1, u_z_j, u_z_jp1, w_bot)

    # ---- reshaping for broadcasting (n_r-1, 1) ------------------------------
    rc_o  = r_c[1:  ].reshape(-1, 1)    # r_c[i],   r-outer face of u-CV
    rc_i  = r_c[:-1 ].reshape(-1, 1)    # r_c[i-1], r-inner face
    rf_int = r_f[1:-1].reshape(-1, 1)   # r_f[i],   where u is stored

    # ---- convective fluxes (2nd-order upwind) --------------------------------
    conv_r = (rc_o * u_out * u_out_upwind - rc_i * u_in * u_in_upwind) / (rf_int * dr)
    conv_z = (w_top * u_top_upwind - w_bot * u_bot_upwind) / dz

    # ---- diffusive fluxes (central, 2nd order) ------------------------------
    # Radial (including −u/r² curvature term):
    visc_r = ((rc_o * (u[2:, :] - u_int) - rc_i * (u_int - u[:-2, :])) / (rf_int * dr**2)
              - u_int / rf_int**2) / Re

    # Axial:
    visc_z = (u_ext[1:-1, 2:] - 2.0 * u_int + u_ext[1:-1, :-2]) / dz**2 / Re

    return -conv_r - conv_z + visc_r + visc_z


def predictor_u(u, w, r_c, r_f, dr, dz, dt, Re):
    """
    Advance radial velocity by one explicit Euler step (no pressure).
    Thin wrapper around compute_rhs_u; the main time loop uses
    compute_rhs_u directly to build the second-order AB2 update.

    Returns
    -------
    u_star : ndarray (n_r+1, n_z)
    """
    rhs = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re)
    u_star = u.copy()
    u_star[1:-1, :] = u[1:-1, :] + dt * rhs
    # Axis and wall faces stay 0 (enforced by apply_bc_u after this call)
    return u_star
