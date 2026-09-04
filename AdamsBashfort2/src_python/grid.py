"""
Grid generation for 2D axisymmetric cylindrical pipe flow.

Staggered grid layout (non-dimensional, r* = r/D, z* = z/D):
  - Cell centres:  r_c[i], z_c[j]      p, T stored here
  - Radial faces:  r_f[i]               u stored here  (shape: n_r+1)
  - Axial faces:   z_f[j]               w stored here  (shape: n_z+1)

Index ranges
  i = 0 .. n_r-1   cell-centre radial index
  j = 0 .. n_z-1   cell-centre axial  index
  i = 0 .. n_r     r-face index  (i=0 → axis r*=0, i=n_r → wall r*=0.5)
  j = 0 .. n_z     z-face index  (j=0 → inlet z*=0, j=n_z → outlet z*=50)
"""

import numpy as np
from constants import n_r, n_z, L_over_D


def make_grid(n_r=n_r, n_z=n_z, L_over_D=L_over_D):
    """
    Build the staggered grid.

    Parameters
    ----------
    n_r, n_z : int
        Number of cells in radial and axial directions.
    L_over_D : float
        Non-dimensional pipe length (default 50).

    Returns
    -------
    r_c : ndarray (n_r,)   cell-centre radial coords  [0 < r* < 0.5]
    r_f : ndarray (n_r+1,) r-face coords              [0 <= r* <= 0.5]
    z_c : ndarray (n_z,)   cell-centre axial coords   [0 < z* < L_over_D]
    z_f : ndarray (n_z+1,) z-face coords              [0 <= z* <= L_over_D]
    dr  : float            uniform radial spacing
    dz  : float            uniform axial  spacing
    """
    dr = 0.5 / n_r          # r* ∈ [0, 0.5]
    dz = L_over_D / n_z     # z* ∈ [0, 50]

    r_f = np.linspace(0.0, 0.5, n_r + 1)          # r-faces
    r_c = 0.5 * (r_f[:-1] + r_f[1:])              # r-cell-centres

    z_f = np.linspace(0.0, L_over_D, n_z + 1)     # z-faces
    z_c = 0.5 * (z_f[:-1] + z_f[1:])              # z-cell-centres

    return r_c, r_f, z_c, z_f, dr, dz
