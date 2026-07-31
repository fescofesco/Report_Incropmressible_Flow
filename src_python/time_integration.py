"""
Explicit 2nd-order (Heun / explicit trapezoidal RK2) time integration for the
pressure-projection pipe-flow solver.

Design choice: pressure is a constraint (Lagrange multiplier enforcing
div(u)=0), not an explicit RHS term, so it cannot be blended stage-by-stage
the way the convective/diffusive RHS terms are. We therefore:

  1. Evaluate the momentum RHS (no pressure) at the current state   -> k1
  2. Take an Euler sub-step to an intermediate state                -> stage 1
  3. Evaluate the momentum RHS again at that intermediate state     -> k2
  4. Combine k1, k2 with the Heun weights to get a 2nd-order-accurate
     "starred" (not yet divergence-free) velocity predictor
  5. Project ONCE (single Poisson solve) to make the combined predictor
     divergence-free

This is the standard simplification used for explicit RK2 projection
methods: continuity only needs to hold for the final combined step, so a
single end-of-step projection is sufficient and avoids solving the Poisson
equation twice per time step. Temperature is advanced the same way (Heun),
using the divergence-free velocity at the new time level for the second
stage's convective transport, consistent with how the corrected velocity is
already used in the original 1st-order scheme.
"""

from solver_momentum import compute_rhs_w, compute_rhs_u
from solver_temperature import compute_rhs_T
from solver_poisson import solve_poisson, correct_velocity, update_pressure
from boundary_conditions import apply_bc_u, apply_bc_w


def rk2_step(u, w, p, T, r_c, r_f, dr, dz, dt, Re, Pr, n_r, n_z, A_poisson, alpha_p):
    """
    Advance (u, w, p, T) by one explicit-trapezoidal (Heun / RK2) step.

    Parameters
    ----------
    u, w, p, T : current fields (BCs already applied)
    r_c, r_f, dr, dz, dt, Re, Pr : as usual
    n_r, n_z    : grid sizes
    A_poisson   : pre-assembled sparse Poisson matrix
    alpha_p     : pressure under-relaxation factor

    Returns
    -------
    u_new, w_new, p_new, T_new
    """
    # ---- Stage 1: RHS at the current state (k1) ------------------------------
    k1_w = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re)
    k1_u = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re)

    w_tilde = w.copy()
    w_tilde[:, 1:-1] = w[:, 1:-1] + dt * k1_w
    u_tilde = u.copy()
    u_tilde[1:-1, :] = u[1:-1, :] + dt * k1_u
    apply_bc_w(w_tilde)
    apply_bc_u(u_tilde)

    # ---- Stage 2: RHS at the intermediate (Euler) state (k2) -----------------
    k2_w = compute_rhs_w(w_tilde, u_tilde, r_c, r_f, dr, dz, Re)
    k2_u = compute_rhs_u(u_tilde, w_tilde, r_c, r_f, dr, dz, Re)

    # ---- Combine (Heun): 2nd-order predictor, pressure not yet applied -------
    w_star = w.copy()
    w_star[:, 1:-1] = w[:, 1:-1] + 0.5 * dt * (k1_w + k2_w)
    u_star = u.copy()
    u_star[1:-1, :] = u[1:-1, :] + 0.5 * dt * (k1_u + k2_u)
    apply_bc_w(w_star)
    apply_bc_u(u_star)

    # ---- Single pressure projection on the combined predictor ----------------
    p_prime = solve_poisson(A_poisson, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
    u_new, w_new = correct_velocity(u_star, w_star, p_prime, r_f, dr, dz, dt, n_r, n_z)
    apply_bc_u(u_new)
    apply_bc_w(w_new)
    p_new = update_pressure(p, p_prime, alpha_p)

    # ---- Temperature: Heun's method, using the new divergence-free velocity
    #      for the second stage's convective transport ------------------------
    k1_T = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr)
    T_tilde = T + dt * k1_T
    T_tilde[:, 0] = 0.0
    T_tilde[:, -1] = T_tilde[:, -2]

    k2_T = compute_rhs_T(T_tilde, u_new, w_new, r_c, r_f, dr, dz, Re, Pr)
    T_new = T + 0.5 * dt * (k1_T + k2_T)
    T_new[:, 0] = 0.0
    T_new[:, -1] = T_new[:, -2]

    return u_new, w_new, p_new, T_new
