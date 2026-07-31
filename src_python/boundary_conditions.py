"""
Boundary conditions for axisymmetric cylindrical pipe flow.

Ghost-cell approach: boundary values are set directly in the
face-centred velocity arrays and enforced via array slices.

Physical boundaries
-------------------
  Centreline  r*=0   (i_face = 0)
  Wall        r*=0.5 (i_face = n_r)
  Inlet       z*=0   (j_face = 0)
  Outlet      z*=L   (j_face = n_z)

Non-dimensional BCs
-------------------
  u[0,  :]  = 0            (symmetry: no radial velocity at axis)
  u[n_r,:]  = 0            (no-slip:  no radial velocity at wall)
  w[:,  0]  = 1            (inlet:    uniform axial velocity w*=W_in/W_in=1)
  w[:, n_z] = w[:,n_z-1]   (outlet:   zero-gradient ∂w/∂z = 0)

  For temperature θ:
    dθ/dr|wall = +Re·Pr   (wall heat flux: ∂θ/∂r* = D·(ρc_vW_in/q_w)·∂T/∂r|_wall
                            = D·(ρc_vW_in/q_w)·(q_w/λ) = ρc_vW_inD/λ = Re·Pr)
    dθ/dr|axis = 0    (symmetry)
    θ[:,0]    = 0     (inlet)
    dθ/dz|out = 0     (outlet)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Velocity BCs (applied directly to face arrays)
# ---------------------------------------------------------------------------

def apply_bc_w(w):
    """
    Enforce axial-velocity boundary conditions.

    Parameters
    ----------
    w : ndarray (n_r, n_z+1)  modified in-place

    BC summary
    ----------
    j=0   inlet:  w[:,0] = 1  (uniform inlet, W_in* = 1)
    j=n_z outlet: zero-gradient ∂w/∂z=0  → w[:,-1] = w[:,-2]
    r=0   axis:   ∂w/∂r=0 satisfied by symmetry (no explicit BC needed
                  for w since it lives at r_c, not r_f; handled in diffusion)
    r=R   wall:   w is defined at r_c, so the wall no-slip for w enters
                  through the ghost cell used in the momentum solver.
                  We do NOT set w here for the wall row; that is handled
                  inside solver_momentum via the ghost extension.
    """
    w[:, 0] = 1.0               # inlet: uniform w* = 1
    w[:, -1] = w[:, -2]         # outlet: zero-gradient


def apply_bc_u(u):
    """
    Enforce radial-velocity boundary conditions.

    Parameters
    ----------
    u : ndarray (n_r+1, n_z)  modified in-place

    BC summary
    ----------
    i=0    centreline:  u=0  (symmetry / no penetration at r=0)
    i=n_r  wall:        u=0  (no penetration at r=R)

    NOTE: u[:,0] is at z_c[0] (first cell centre), NOT at the inlet face z_f[0].
    The inlet BC u=0 at z=0 is enforced via the ghost cell in the predictor
    (u_ghost_in = -u[:,0], so the average at z_f[0] is zero).
    We must NOT override u[:,0] here, because the Poisson correction sets
    it to enforce div(u)=0 at the inlet cells.
    """
    u[0,  :] = 0.0              # axis (r_f[0] = 0)
    u[-1, :] = 0.0              # wall (r_f[n_r] = R)


# ---------------------------------------------------------------------------
# Temperature BCs
# ---------------------------------------------------------------------------

def apply_bc_T(T, dr, n_r, Re, Pr):
    """
    Enforce non-dimensional temperature boundary conditions.

    Non-dimensional wall heat-flux BC:
        ∂θ/∂r*|_{r*=0.5} = +Re·Pr
    This is implemented as a ghost-cell value inside solver_temperature.py:
        T_ghost_wall = T[n_r-1, :] + dr * Re*Pr

    Parameters
    ----------
    T    : ndarray (n_r, n_z)  modified in-place
    dr   : float
    n_r  : int
    Re, Pr : float  (not used for T array itself, but kept for signature clarity)

    Note: The wall BC is applied via ghost values inside solver_temperature;
    here we only enforce inlet and outlet conditions on the T array.
    """
    T[:, 0]  = 0.0              # inlet: θ = 0
    T[:, -1] = T[:, -2]         # outlet: zero-gradient ∂θ/∂z = 0


# ---------------------------------------------------------------------------
# Composite: apply all BCs at once
# ---------------------------------------------------------------------------

def apply_all_bc(u, w, T, dr, n_r, Re, Pr):
    """Apply all boundary conditions in one call."""
    apply_bc_u(u)
    apply_bc_w(w)
    apply_bc_T(T, dr, n_r, Re, Pr)
