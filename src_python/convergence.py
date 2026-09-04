"""
Convergence monitoring for the pipe-flow simulation.

Three residuals are tracked:
  R_cont  : max |div(u)| over all cells     (continuity)
  R_vel   : max change rate of both u and w (momentum change rate)
  R_temp  : max |θ^{n+1} - θ^n| / dt       (temperature change rate)
"""

import numpy as np


def compute_divergence(u, w, r_c, r_f, dr, dz):
    """
    Compute the pointwise discrete divergence of the velocity field.

    div[i,j] = [r_{i+½}·u[i+1,j] − r_{i-½}·u[i,j]] / (r_c[i]·dr)
             + [w[i,j+1] − w[i,j]] / dz

    Parameters
    ----------
    u  : ndarray (n_r+1, n_z)
    w  : ndarray (n_r, n_z+1)
    r_c: ndarray (n_r,)
    r_f: ndarray (n_r+1,)
    dr, dz : float

    Returns
    -------
    div : ndarray (n_r, n_z)   pointwise divergence
    """
    rf_o = r_f[1:].reshape(-1, 1)
    rf_i = r_f[:-1].reshape(-1, 1)
    rc   = r_c.reshape(-1, 1)

    div_u = (rf_o * u[1:, :] - rf_i * u[:-1, :]) / (rc * dr)
    div_w = (w[:, 1:] - w[:, :-1]) / dz

    return div_u + div_w


def compute_residuals(u_new, u_old, w_new, w_old, T_new, T_old,
                      r_c, r_f, dr, dz, dt):
    """
    Compute all residuals for convergence checking.

    Parameters
    ----------
    u_new, u_old : ndarray (n_r+1, n_z)
    w_new, w_old : ndarray (n_r, n_z+1)
    T_new, T_old : ndarray (n_r, n_z)
    r_c, r_f, dr, dz, dt : as usual

    Returns
    -------
    residuals : dict with keys 'continuity', 'momentum', 'temperature'
    """
    R_cont = float(np.max(np.abs(compute_divergence(
        u_new, w_new, r_c, r_f, dr, dz))))
    R_u = np.max(np.abs(u_new[1:-1, :] - u_old[1:-1, :]))
    R_w = np.max(np.abs(w_new[:, 1:-1] - w_old[:, 1:-1]))
    R_vel = float(max(R_u, R_w) / dt)
    R_temp = float(np.max(np.abs(T_new - T_old)) / dt)

    return {
        'continuity' : R_cont,
        'momentum'   : R_vel,
        'temperature': R_temp,
    }


def check_convergence(residuals, tol_cont, tol_vel, tol_temp):
    """
    Return True when all residuals are below their tolerances.

    Parameters
    ----------
    residuals : dict from compute_residuals()
    tol_cont, tol_vel, tol_temp : float

    Returns
    -------
    converged : bool
    """
    return (residuals['continuity']  < tol_cont and
            residuals['momentum']    < tol_vel  and
            residuals['temperature'] < tol_temp)
