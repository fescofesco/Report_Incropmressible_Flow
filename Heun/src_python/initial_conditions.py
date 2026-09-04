"""
Initial conditions for the pipe-flow simulation.

All fields are non-dimensional:
  u*  radial velocity  - shape (n_r+1, n_z)  at r-faces
  w*  axial  velocity  - shape (n_r, n_z+1)  at z-faces
  p*  pressure         - shape (n_r, n_z)    at cell centres
  T   non-dim temp θ   - shape (n_r, n_z)    at cell centres

Start-up state: uniform axial flow (w=1, u=0, p=0, θ=0).
"""

import numpy as np


def initialise_fields(n_r, n_z):
    """
    Allocate and initialise all field arrays.

    Parameters
    ----------
    n_r, n_z : int

    Returns
    -------
    u  : ndarray (n_r+1, n_z)   radial velocity  (zero)
    w  : ndarray (n_r,   n_z+1) axial  velocity  (uniform 1)
    p  : ndarray (n_r,   n_z)   pressure         (zero)
    T  : ndarray (n_r,   n_z)   temperature θ    (zero)
    """
    u = np.zeros((n_r + 1, n_z))           # u at r-faces
    w = np.ones((n_r, n_z + 1))            # w at z-faces, uniform inlet profile
    p = np.zeros((n_r, n_z))               # p at cell centres
    T = np.zeros((n_r, n_z))               # θ at cell centres

    return u, w, p, T
