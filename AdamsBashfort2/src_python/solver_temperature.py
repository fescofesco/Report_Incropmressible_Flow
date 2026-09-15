"""
Energy equation solver (non-dimensional temperature θ).

Non-dimensional energy equation (no viscous dissipation):
  ∂θ/∂t + (1/r)∂(r·u·θ)/∂r + ∂(w·θ)/∂z
      = (1/(Re·Pr))·[ (1/r)∂/∂r(r·∂θ/∂r) + ∂²θ/∂z² ]

Written in conservation form; identical FVM structure to the
axial momentum equation, with 1/Re → 1/(Re·Pr).

Convective terms use a 2ND-ORDER UPWIND (linear-extrapolation / LUD) scheme
— see numerics.upwind2_face — which retains the upwind stability bias at
high grid Peclet numbers (Pe_thermal = Re·Pr·dz = 100) while being formally
2nd-order accurate, matching the diffusive terms (standard 2nd-order
central differencing).

θ is stored at cell centres (n_r, n_z).  Velocities u and w at the
faces interpolate θ to the face locations.

Wall BC:  ∂θ/∂r*|_{r*=0.5} = +Re·Pr
          Derivation: dimensionally, Fourier's law at the wall gives
          ∂T/∂r|_wall = q_w/λ (heat flows into the fluid, wall hotter
          than centerline). Non-dimensionalizing with r*=r/D and
          θ=(T-T_in)ρc_vW_in/q_w:
              ∂θ/∂r* = D·(ρc_vW_in/q_w)·∂T/∂r = ρc_vW_inD/λ = Re·Pr
          (NOT ±1 — that earlier value dropped the Re·Pr factor from
          the non-dimensionalization; see report Appendix B)
          Implemented as a ghost cell:  θ_ghost[n_r, j] = θ[n_r-1, j] + dr·Re·Pr  (Neumann)
Axis BC:  ∂θ/∂r|_{r*=0}   = 0
          Ghost: θ_ghost[-1, j] = θ[0, j]
Inlet:    θ = 0 at z*=0 (diffusive ghost plus exact prescribed inflow flux)
Outlet:   ∂θ/∂z = 0     (ghost θ[:, n_z] = θ[:, n_z-1])
"""

import numpy as np

from numerics import upwind2_face


def compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr):
    """
    Evaluate the explicit RHS of the energy equation (no time discretization):

        RHS = -conv_r - conv_z + visc_r + visc_z

    so the AB2 integrator can combine the current and previous evaluations.

    Convective fluxes use upwind differencing.

    Parameters
    ----------
    T   : ndarray (n_r, n_z)    current θ  (BCs already applied)
    u   : ndarray (n_r+1, n_z)  radial velocity at r-faces
    w   : ndarray (n_r, n_z+1)  axial  velocity at z-faces
    r_c : ndarray (n_r,)
    r_f : ndarray (n_r+1,)
    dr, dz : float
    Re, Pr : float

    Returns
    -------
    rhs : ndarray (n_r, n_z)
    """
    n_r, n_z = T.shape
    alpha = 1.0 / (Re * Pr)   # thermal diffusivity (non-dimensional)
    RePr = Re * Pr

    # ---- ghost extension in r -----------------------------------------------
    # Axis  ghost: ∂θ/∂r=0    → θ_ghost = θ[0,:]
    # Wall  ghost: ∂θ/∂r=+RePr → θ_ghost = θ[n_r-1,:] + dr*RePr  (centred-difference)
    T_ghost_axis = T[0:1, :]                        # (1, n_z)
    T_ghost_wall = T[-1:, :] + dr * RePr             # (1, n_z), gradient = (dr*RePr)/dr = RePr

    T_ext = np.concatenate([T_ghost_axis, T, T_ghost_wall], axis=0)  # (n_r+2, n_z)
    # T_ext[i+1, j] = T[i, j]

    # ---- ghost extension in z -----------------------------------------------
    # Inlet  (j=-1): θ=0 → ghost = -T[:,0]  so avg = 0 at inlet face
    # Outlet (j=n_z): zero-gradient → ghost = T[:,n_z-1]
    T_ghost_in  = -T[:, 0:1]
    T_ghost_out =  T[:, -1:]
    # Combine with T for axial stencil:  shape (n_r, n_z+2)  (used by diffusion)
    T_zext = np.concatenate([T_ghost_in, T, T_ghost_out], axis=1)

    # ---- second ghost layer, needed by the 2nd-order-upwind convective stencil
    # Axis:  mirror one cell further in  → θ_ghost2 = θ[1,:]
    # Wall:  continue the same linear (constant-gradient) extrapolation
    T_ghost_axis2 = T[1:2, :]
    T_ghost_wall2 = T[-1:, :] + 2.0 * dr * RePr
    T_ext2 = np.concatenate(
        [T_ghost_axis2, T_ghost_axis, T, T_ghost_wall, T_ghost_wall2], axis=0)  # (n_r+4, n_z)
    # T_ext2[i+2, j] = T[i, j]

    # Inlet:  continue the odd (Dirichlet) reflection one cell further
    # Outlet: continue the zero-gradient extrapolation
    T_ghost_in2  = -T[:, 1:2]
    T_ghost_out2 =  T[:, -1:]
    T_zext2 = np.concatenate(
        [T_ghost_in2, T_ghost_in, T, T_ghost_out, T_ghost_out2], axis=1)  # (n_r, n_z+4)
    # T_zext2[:, j+2] = T[:, j]

    # ---- radial convective flux  (1/r)∂(r·u·θ)/∂r --------------------------
    # u at r-faces (directly from u array):
    u_E = u[1:, :]    # u at r_f[i+1], (n_r, n_z)
    u_I = u[:-1, :]   # u at r_f[i],   (n_r, n_z)

    T_r_im2 = T_ext2[0:n_r, :]
    T_r_im1 = T_ext2[1:n_r + 1, :]
    T_r_i   = T_ext2[2:n_r + 2, :]     # = T
    T_r_ip1 = T_ext2[3:n_r + 3, :]
    T_r_ip2 = T_ext2[4:n_r + 4, :]

    # Outer r-face (r_f[i+1], between cells i and i+1)
    T_rf_o_upwind = upwind2_face(T_r_im1, T_r_i, T_r_ip1, T_r_ip2, u_E)
    # Inner r-face (r_f[i], between cells i-1 and i)
    T_rf_i_upwind = upwind2_face(T_r_im2, T_r_im1, T_r_i, T_r_ip1, u_I)

    rf_o = r_f[1:].reshape(-1, 1)    # (n_r, 1)
    rf_i = r_f[:-1].reshape(-1, 1)   # (n_r, 1)
    rc   = r_c.reshape(-1, 1)        # (n_r, 1)

    conv_r = (rf_o * u_E * T_rf_o_upwind - rf_i * u_I * T_rf_i_upwind) / (rc * dr)

    # ---- axial convective flux  ∂(w·θ)/∂z -----------------------------------
    # w at z-faces (directly):
    w_tp = w[:, 1:]    # w at z_f[j+1], (n_r, n_z)
    w_bt = w[:, :-1]   # w at z_f[j],   (n_r, n_z)

    T_z_jm2 = T_zext2[:, 0:n_z]
    T_z_jm1 = T_zext2[:, 1:n_z + 1]
    T_z_j   = T_zext2[:, 2:n_z + 2]    # = T
    T_z_jp1 = T_zext2[:, 3:n_z + 3]
    T_z_jp2 = T_zext2[:, 4:n_z + 4]

    # Top face z_f[j+1] (between cells j and j+1)
    T_zf_p_upwind = upwind2_face(T_z_jm1, T_z_j, T_z_jp1, T_z_jp2, w_tp)
    # Bottom face z_f[j] (between cells j-1 and j)
    T_zf_m_upwind = upwind2_face(T_z_jm2, T_z_jm1, T_z_j, T_z_jp1, w_bt)

    # At the physical inlet face the incoming transported value is prescribed,
    # not reconstructed: theta(z*=0)=0.  The odd ghosts above are still needed
    # for the second-order axial diffusion stencil.  Letting LUD extrapolate at
    # this inflow face would generally produce a nonzero boundary heat flux.
    T_zf_m_upwind[:, 0] = np.where(w_bt[:, 0] > 0.0,
                                    0.0, T_zf_m_upwind[:, 0])

    conv_z = (w_tp * T_zf_p_upwind - w_bt * T_zf_m_upwind) / dz

    # ---- radial diffusion  (1/r)∂/∂r(r·∂θ/∂r) ------------------------------
    visc_r = alpha * (rf_o * (T_ext[2:, :] - T)
                    - rf_i * (T - T_ext[:-2, :])) / (rc * dr**2)

    # ---- axial diffusion  ∂²θ/∂z² -------------------------------------------
    visc_z = alpha * (T_zext[:, 2:] - 2.0 * T + T_zext[:, :-2]) / dz**2

    return -conv_r - conv_z + visc_r + visc_z


def advance_temperature(T, u, w, r_c, r_f, dr, dz, dt, Re, Pr):
    """
    Advance temperature by one explicit Euler step (thin wrapper around
    compute_rhs_T). Kept for simple/one-off use; the main time loop uses
    compute_rhs_T directly to build the second-order AB2 update.

    Returns
    -------
    T_new : ndarray (n_r, n_z)
    """
    rhs = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr)
    T_new = T + dt * rhs

    # Inlet/outlet are enforced via ghost cells inside compute_rhs_T; no
    # direct cell-value overwrite is needed.
    return T_new
