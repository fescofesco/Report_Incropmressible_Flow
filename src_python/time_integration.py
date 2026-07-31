"""
Explicit 2nd-order (Heun / explicit trapezoidal RK2) time integration for the
pressure-projection pipe-flow solver.

Design: pressure is a constraint (Lagrange multiplier enforcing div(u)=0),
not an explicit RHS term, so it cannot be blended stage-by-stage the way
the convective/diffusive RHS terms are. Each RK stage's velocity therefore
gets its OWN pressure projection before being used to evaluate the next
stage's RHS (two Poisson solves per step) -- evaluating k2 on a velocity
field that doesn't satisfy continuity would corrupt the convective terms
(which assume a solenoidal transport velocity in this conservative FVM
formulation) and silently degrades the scheme's temporal accuracy below
2nd order.

  1. Evaluate the momentum RHS (no pressure) at the current, already
     divergence-free state                                    -> k1
  2. Euler sub-step to an intermediate predictor, then PROJECT it
     (Poisson solve #1) to get a divergence-free intermediate state u1,w1
  3. Evaluate the momentum RHS at that divergence-free intermediate state -> k2
  4. Combine k1, k2 with the Heun weights to get a 2nd-order-accurate
     "starred" (not yet divergence-free) velocity predictor
  5. Project the combined predictor (Poisson solve #2) to get the final
     divergence-free u^{n+1}, w^{n+1}

Temperature is advanced the same way (Heun), with its second stage using
the divergence-free INTERMEDIATE velocity (u1, w1) -- the same state used
for the momentum RHS's k2 -- rather than the final u^{n+1}, w^{n+1}, so
that both k2 evaluations are consistent with the same point in the
algorithm's stage structure.
"""

from solver_momentum import compute_rhs_w, compute_rhs_u
from solver_temperature import compute_rhs_T
from solver_poisson import solve_poisson, correct_velocity, update_pressure
from boundary_conditions import apply_bc_u, apply_bc_w


def rk2_step(u, w, p, T, r_c, r_f, dr, dz, dt, Re, Pr, n_r, n_z, A_poisson, alpha_p):
    """
    Advance (u, w, p, T) by one explicit-trapezoidal (Heun / RK2) step, with
    each velocity stage individually projected to remain divergence-free.

    Parameters
    ----------
    u, w, p, T : current fields (BCs already applied; u, w divergence-free)
    r_c, r_f, dr, dz, dt, Re, Pr : as usual
    n_r, n_z    : grid sizes
    A_poisson   : pre-assembled sparse Poisson matrix
    alpha_p     : pressure relaxation factor (blend towards the new p', not
                  an accumulation -- see solver_poisson.update_pressure)

    Returns
    -------
    u_new, w_new, p_new, T_new
    """
    # ---- Stage 1: RHS at the current (divergence-free) state (k1) ------------
    k1_w = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re)
    k1_u = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re)

    w1_star = w.copy()
    w1_star[:, 1:-1] = w[:, 1:-1] + dt * k1_w
    u1_star = u.copy()
    u1_star[1:-1, :] = u[1:-1, :] + dt * k1_u
    apply_bc_w(w1_star)
    apply_bc_u(u1_star)

    # Project the stage-1 predictor so k2 is evaluated on a divergence-free
    # intermediate velocity field (u1, w1), not a raw Euler predictor.
    p1_prime = solve_poisson(A_poisson, u1_star, w1_star, r_c, r_f, dr, dz, dt, n_r, n_z)
    u1, w1 = correct_velocity(u1_star, w1_star, p1_prime, r_f, dr, dz, dt, n_r, n_z)
    apply_bc_u(u1)
    apply_bc_w(w1)

    # ---- Stage 2: RHS at the divergence-free intermediate state (k2) ---------
    k2_w = compute_rhs_w(w1, u1, r_c, r_f, dr, dz, Re)
    k2_u = compute_rhs_u(u1, w1, r_c, r_f, dr, dz, Re)

    # ---- Combine (Heun): 2nd-order predictor, pressure not yet applied -------
    w_star = w.copy()
    w_star[:, 1:-1] = w[:, 1:-1] + 0.5 * dt * (k1_w + k2_w)
    u_star = u.copy()
    u_star[1:-1, :] = u[1:-1, :] + 0.5 * dt * (k1_u + k2_u)
    apply_bc_w(w_star)
    apply_bc_u(u_star)

    # ---- Final pressure projection on the combined predictor -----------------
    p_prime = solve_poisson(A_poisson, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
    u_new, w_new = correct_velocity(u_star, w_star, p_prime, r_f, dr, dz, dt, n_r, n_z)
    apply_bc_u(u_new)
    apply_bc_w(w_new)
    p_new = update_pressure(p, p_prime, alpha_p)

    # ---- Temperature: Heun's method, using the divergence-free INTERMEDIATE
    #      velocity (u1, w1) for the second stage's convective transport,
    #      consistent with the momentum k2 evaluation above -------------------
    k1_T = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr)
    T_tilde = T + dt * k1_T
    T_tilde[:, 0] = 0.0
    T_tilde[:, -1] = T_tilde[:, -2]

    k2_T = compute_rhs_T(T_tilde, u1, w1, r_c, r_f, dr, dz, Re, Pr)
    T_new = T + 0.5 * dt * (k1_T + k2_T)
    T_new[:, 0] = 0.0
    T_new[:, -1] = T_new[:, -2]

    return u_new, w_new, p_new, T_new
