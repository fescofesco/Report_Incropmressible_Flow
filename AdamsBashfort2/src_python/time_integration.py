"""Second-order Adams-Bashforth incremental pressure-correction step."""

from boundary_conditions import apply_bc_u, apply_bc_w
from solver_momentum import compute_rhs_u, compute_rhs_w
from solver_poisson import correct_velocity, solve_poisson
from solver_temperature import compute_rhs_T


def ab2_step(u, w, p, T, rhs_previous, r_c, r_f, dr, dz, dt, Re, Pr,
             n_r, n_z, A_poisson, alpha_p):
    """Advance one step; use Euler when the AB2 history is unavailable."""
    rhs_u = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re)
    rhs_w = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re)
    rhs_T = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr)

    if rhs_previous is None:
        explicit_u, explicit_w, explicit_T = rhs_u, rhs_w, rhs_T
    else:
        rhs_u_old, rhs_w_old, rhs_T_old = rhs_previous
        explicit_u = 1.5 * rhs_u - 0.5 * rhs_u_old
        explicit_w = 1.5 * rhs_w - 0.5 * rhs_w_old
        explicit_T = 1.5 * rhs_T - 0.5 * rhs_T_old

    # Including p^n here makes the Poisson unknown a correction p'.
    u_star = u.copy()
    grad_p_r = (p[1:, :] - p[:-1, :]) / dr
    u_star[1:-1, :] = u[1:-1, :] + dt * (explicit_u - grad_p_r)

    w_star = w.copy()
    grad_p_z = (p[:, 1:] - p[:, :-1]) / dz
    w_star[:, 1:-1] = w[:, 1:-1] + dt * (explicit_w - grad_p_z)
    apply_bc_u(u_star)
    apply_bc_w(w_star)

    p_prime = solve_poisson(
        A_poisson, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
    u_new, w_new = correct_velocity(
        u_star, w_star, p_prime, r_f, dr, dz, dt, n_r, n_z)
    apply_bc_u(u_new)
    apply_bc_w(w_new)

    p_new = p + alpha_p * p_prime
    T_new = T + dt * explicit_T
    rhs_current = (rhs_u, rhs_w, rhs_T)
    return u_new, w_new, p_new, T_new, rhs_current
