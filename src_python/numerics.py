"""
Small shared numerical building blocks used by more than one solver module.
"""

import numpy as np


def upwind2_face(phi_im1, phi_i, phi_ip1, phi_ip2, vel):
    """
    2nd-order upwind (linear-extrapolation / LUD) value of a transported
    quantity at a face, given the two cell values on each side of the face
    and the transport velocity at that face.

        face value = 1.5*phi_i   - 0.5*phi_im1   if vel > 0  (upwind: low-index side)
                   = 1.5*phi_ip1 - 0.5*phi_ip2    if vel < 0  (upwind: high-index side)

    This reduces to the classic 3-point linear extrapolation from the two
    cells upstream of the face, giving 2nd-order spatial accuracy while
    keeping the upwind bias needed for stability at high cell Peclet number
    (unlike pure central differencing).

    Parameters
    ----------
    phi_im1, phi_i, phi_ip1, phi_ip2 : ndarray, all the same shape
        Cell values at relative positions i-1, i, i+1, i+2 with respect to
        the face sitting between cells i and i+1.
    vel : ndarray, same shape
        Transport velocity at the face.

    Returns
    -------
    ndarray, same shape
    """
    return np.where(vel > 0,
                     1.5 * phi_i - 0.5 * phi_im1,
                     1.5 * phi_ip1 - 0.5 * phi_ip2)
