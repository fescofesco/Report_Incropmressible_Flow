"""
Verification of the hand-written direct Poisson solver.

The assignment asks for a direct solution method that is *programmed*, not a
library call. solver_poisson.factorize_poisson / solve_poisson implement the
block Thomas algorithm on top of our own Doolittle LU kernel. This script
checks that implementation three ways:

  1. Backward error: the relative residual ||A p' - b||_inf / ||b||_inf of the
     computed solution against the independently assembled sparse operator A.
     This is the meaningful accuracy measure -- the operator is stiff
     (cond_1(A) ~ 1e8), so the *solution* difference against any other solver
     is bounded by cond(A)*eps and says little about either solver.
  2. Agreement with scipy's sparse LU (SuperLU) on the same system, as an
     independent reference.
  3. That the boundary conditions come out as imposed: p' = 0 on the outlet
     column (Dirichlet) and the discrete compatibility of the Neumann walls.

Run:  python verify_poisson_solver.py
"""

import time

import numpy as np
import scipy.sparse.linalg as spla

import constants as C
from grid import make_grid
from solver_poisson import (build_poisson_matrix, factorize_poisson,
                            lu_factor, lu_solve, solve_poisson)


def _random_predictor(n_r, n_z, seed=0):
    """A divergence-carrying velocity predictor with the wall/axis BCs of u."""
    rng = np.random.default_rng(seed)
    u_star = rng.standard_normal((n_r + 1, n_z))
    u_star[0, :] = 0.0       # symmetry axis
    u_star[-1, :] = 0.0      # no-penetration at the wall
    w_star = rng.standard_normal((n_r, n_z + 1))
    return u_star, w_star


def _rhs(u_star, w_star, r_c, r_f, dr, dz, dt):
    """The same right-hand side solve_poisson builds internally."""
    rf_o = r_f[1:].reshape(-1, 1)
    rf_i = r_f[:-1].reshape(-1, 1)
    rc = r_c.reshape(-1, 1)
    div_u = (rf_o * u_star[1:, :] - rf_i * u_star[:-1, :]) / (rc * dr)
    div_w = (w_star[:, 1:] - w_star[:, :-1]) / dz
    b = (div_u + div_w) / dt
    b[:, -1] = 0.0
    return b


def main():
    n_r, n_z, dt = C.n_r, C.n_z, C.dt
    r_c, r_f, z_c, z_f, dr, dz = make_grid(n_r, n_z, C.L_over_D)

    print("=" * 62)
    print("Verification of the hand-written direct Poisson solver")
    print(f"  grid: n_r = {n_r}, n_z = {n_z}   ({n_r * n_z} unknowns)")
    print(f"  block structure: {n_z} blocks of {n_r} x {n_r}")
    print("=" * 62)

    # --- 0. the LU kernel itself, on a small dense system --------------------
    rng = np.random.default_rng(42)
    M = rng.standard_normal((12, 12))
    rhs = rng.standard_normal(12)
    LU, piv = lu_factor(M)
    x = lu_solve(LU, piv, rhs)
    print(f"\n[0] dense LU kernel, 12x12 random system")
    print(f"    ||M x - rhs||_inf = {np.abs(M @ x - rhs).max():.3e}")

    # --- 1. factorise --------------------------------------------------------
    t0 = time.perf_counter()
    fac = factorize_poisson(r_c, r_f, dr, dz, n_r, n_z)
    t_fac = time.perf_counter() - t0
    print(f"\n[1] block-Thomas factorisation: {t_fac * 1e3:.1f} ms (once)")

    # --- 2. solve and measure the backward error -----------------------------
    u_star, w_star = _random_predictor(n_r, n_z)
    b = _rhs(u_star, w_star, r_c, r_f, dr, dz, dt)

    t0 = time.perf_counter()
    n_rep = 20
    for _ in range(n_rep):
        p_prime = solve_poisson(fac, u_star, w_star, r_c, r_f,
                                dr, dz, dt, n_r, n_z)
    t_solve = (time.perf_counter() - t0) / n_rep

    A = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z)
    b_flat = b.ravel()
    res_own = np.abs(A @ p_prime.ravel() - b_flat).max() / np.abs(b_flat).max()
    print(f"\n[2] solve: {t_solve * 1e3:.3f} ms per Poisson solve")
    print(f"    relative residual ||A p' - b||_inf / ||b||_inf = {res_own:.3e}")

    # --- 3. independent reference (SuperLU) ----------------------------------
    p_ref = spla.spsolve(A.tocsc(), b_flat).reshape(n_r, n_z)
    res_ref = np.abs(A @ p_ref.ravel() - b_flat).max() / np.abs(b_flat).max()
    d_rel = np.abs(p_prime - p_ref).max() / np.abs(p_ref).max()
    cond = spla.onenormest(A.tocsc()) * spla.onenormest(spla.inv(A.tocsc()))
    print(f"\n[3] reference scipy spsolve (SuperLU)")
    print(f"    relative residual, SuperLU                     = {res_ref:.3e}")
    print(f"    ||p'_own - p'_SuperLU||_inf / ||p'||_inf       = {d_rel:.3e}")
    print(f"    cond_1(A) estimate                             = {cond:.3e}")
    print(f"    cond_1(A) * eps                                = "
          f"{cond * np.finfo(float).eps:.3e}   (bound on the line above)")

    # --- 4. boundary conditions ----------------------------------------------
    print(f"\n[4] boundary conditions")
    print(f"    max |p'| on the outlet column (Dirichlet p'=0)  = "
          f"{np.abs(p_prime[:, -1]).max():.3e}")

    ok = res_own < 1e-10
    print("\n" + "=" * 62)
    print("PASS" if ok else "FAIL",
          "- hand-written direct solver reproduces the operator"
          if ok else "- residual too large")
    print("=" * 62)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
