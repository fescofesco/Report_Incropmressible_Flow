"""
Energy equation solver (non-dimensional temperature θ).

Non-dimensional energy equation (no viscous dissipation):
  ∂θ/∂t + (1/r)∂(r·u·θ)/∂r + ∂(w·θ)/∂z
      = (1/(Re·Pr))·[ (1/r)∂/∂r(r·∂θ/∂r) + ∂²θ/∂z² ]

Written in conservation form; identical FVM structure to the
axial momentum equation, with 1/Re → 1/(Re·Pr).

Convective terms use FIRST-ORDER UPWIND differencing to ensure
stability at high grid Peclet numbers (Pe_thermal = Re·Pr·dz = 100).
Diffusive terms use standard second-order central differencing.

θ is stored at cell centres (n_r, n_z).  Velocities u and w at the
faces interpolate θ to the face locations.

Wall BC:  ∂θ/∂r*|_{r*=0.5} = −1
          Implemented as a ghost cell:  θ_ghost[n_r, j] = θ[n_r-1, j] − dr  (Neumann)
Axis BC:  ∂θ/∂r|_{r*=0}   = 0
          Ghost: θ_ghost[-1, j] = θ[0, j]
Inlet:    θ[:, 0] = 0   (set in apply_bc_T)
Outlet:   ∂θ/∂z = 0     (ghost θ[:, n_z] = θ[:, n_z-1])
"""

import numpy as np


def advance_temperature(T, u, w, r_c, r_f, dr, dz, dt, Re, Pr):
    """
    Advance temperature by one explicit Euler step.

    Convective fluxes use upwind differencing.

    Parameters
    ----------
    T   : ndarray (n_r, n_z)    current θ  (BCs already applied)
    u   : ndarray (n_r+1, n_z)  radial velocity at r-faces
    w   : ndarray (n_r, n_z+1)  axial  velocity at z-faces
    r_c : ndarray (n_r,)
    r_f : ndarray (n_r+1,)
    dr, dz, dt : float
    Re, Pr : float

    Returns
    -------
    T_new : ndarray (n_r, n_z)
    """
    n_r, n_z = T.shape
    alpha = 1.0 / (Re * Pr)   # thermal diffusivity (non-dimensional)

    # ---- ghost extension in r -----------------------------------------------
    # Axis  ghost: ∂θ/∂r=0 → θ_ghost = θ[0,:]
    # Wall  ghost: ∂θ/∂r=-1 → θ_ghost = θ[n_r-1,:] - dr  (centred-difference)
    T_ghost_axis = T[0:1, :]                   # (1, n_z)
    T_ghost_wall = T[-1:, :] - dr              # (1, n_z),  gives gradient = -1/dr*dr = -1

    T_ext = np.concatenate([T_ghost_axis, T, T_ghost_wall], axis=0)  # (n_r+2, n_z)
    # T_ext[i+1, j] = T[i, j]

    # ---- ghost extension in z -----------------------------------------------
    # Inlet  (j=-1): θ=0 → ghost = -T[:,0]  so avg = 0 at inlet face
    # Outlet (j=n_z): zero-gradient → ghost = T[:,n_z-1]
    T_ghost_in  = -T[:, 0:1]
    T_ghost_out =  T[:, -1:]
    # Combine with T for axial stencil:  shape (n_r, n_z+2)
    T_zext = np.concatenate([T_ghost_in, T, T_ghost_out], axis=1)

    # ---- radial convective flux  (1/r)∂(r·u·θ)/∂r --------------------------
    # u at r-faces (directly from u array):
    u_E = u[1:, :]    # u at r_f[i+1], (n_r, n_z)
    u_I = u[:-1, :]   # u at r_f[i],   (n_r, n_z)

    # Upwind θ at outer r-face (r_f[i+1]):
    #   u_E > 0 (outward): T from inner cell → T[i,j]
    #   u_E < 0 (inward):  T from outer cell → T_ext[i+2,j] = T[i+1,j] or wall ghost
    T_rf_o_upwind = np.where(u_E > 0, T, T_ext[2:, :])

    # Upwind θ at inner r-face (r_f[i]):
    #   u_I > 0 (outward): T from further inner → T_ext[i,j] = T[i-1,j] or axis ghost
    #   u_I < 0 (inward):  T from this cell → T[i,j]
    T_rf_i_upwind = np.where(u_I > 0, T_ext[:-2, :], T)

    rf_o = r_f[1:].reshape(-1, 1)    # (n_r, 1)
    rf_i = r_f[:-1].reshape(-1, 1)   # (n_r, 1)
    rc   = r_c.reshape(-1, 1)        # (n_r, 1)

    conv_r = (rf_o * u_E * T_rf_o_upwind - rf_i * u_I * T_rf_i_upwind) / (rc * dr)

    # ---- axial convective flux  ∂(w·θ)/∂z -----------------------------------
    # w at z-faces (directly):
    w_tp = w[:, 1:]    # w at z_f[j+1], (n_r, n_z)
    w_bt = w[:, :-1]   # w at z_f[j],   (n_r, n_z)

    # Upwind θ at top face z_f[j+1]:
    #   w_tp > 0 (flow to right): T from south → T[i,j] = T_zext[:, 1:-1]
    #   w_tp < 0 (flow to left):  T from north → T[i,j+1] = T_zext[:, 2:]
    T_zf_p_upwind = np.where(w_tp > 0, T_zext[:, 1:-1], T_zext[:, 2:])

    # Upwind θ at bottom face z_f[j]:
    #   w_bt > 0: T from south → T[i,j-1] = T_zext[:, :-2]
    #   w_bt < 0: T from north → T[i,j] = T_zext[:, 1:-1]
    T_zf_m_upwind = np.where(w_bt > 0, T_zext[:, :-2], T_zext[:, 1:-1])

    conv_z = (w_tp * T_zf_p_upwind - w_bt * T_zf_m_upwind) / dz

    # ---- radial diffusion  (1/r)∂/∂r(r·∂θ/∂r) ------------------------------
    visc_r = alpha * (rf_o * (T_ext[2:, :] - T)
                    - rf_i * (T - T_ext[:-2, :])) / (rc * dr**2)

    # ---- axial diffusion  ∂²θ/∂z² -------------------------------------------
    visc_z = alpha * (T_zext[:, 2:] - 2.0 * T + T_zext[:, :-2]) / dz**2

    # ---- time step -----------------------------------------------------------
    T_new = T + dt * (-conv_r - conv_z + visc_r + visc_z)

    # Re-apply inlet BC (may be polluted by convective term at j=0)
    T_new[:, 0] = 0.0

    return T_new
