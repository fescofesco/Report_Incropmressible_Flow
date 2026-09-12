"""
Main simulation loop for 2D axisymmetric incompressible pipe flow.

Algorithm: second-order Adams-Bashforth integration with an incremental
pressure correction and one Poisson solve per time step. Forward Euler is
used for the first step to initialise the multistep history.

Physical problem
----------------
  2D axisymmetric laminar flow in a heated cylindrical pipe.
  Re = 100, Pr = 5, L/D = 50.
  Governing equations in conservation form (see 2_2D_cylindrical_coordinates_uw.tex).
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend for saving figures
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla

# Ensure src_python is on the path when running from repo root
sys.path.insert(0, os.path.dirname(__file__))

from constants import (Re, Pr, n_r, n_z, dt, n_steps, output_interval,
                       tol_continuity, tol_velocity, tol_temperature, alpha_p)
from grid import make_grid
from initial_conditions import initialise_fields
from boundary_conditions import apply_all_bc
from solver_poisson import build_poisson_matrix
from time_integration import ab2_step
from convergence import compute_residuals, check_convergence
from analytical import (calculate_analytical_velocity_nondim,
                         calculate_analytical_temperature_nondim)
from plotting import (plot_velocity_contour_axial, plot_velocity_contour_radial,
                      plot_temperature_contour, plot_velocity_profile,
                      plot_temperature_profile, plot_convergence_history,
                      plot_comparison_with_analytical)


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
PLOTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'Plots_python')
os.makedirs(PLOTS_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Helper: interpolate staggered fields to cell centres for plotting
# ---------------------------------------------------------------------------

def w_to_cell_centre(w):
    """w at z-faces → average to cell centres.  Returns shape (n_r, n_z)."""
    return 0.5 * (w[:, :-1] + w[:, 1:])


def u_to_cell_centre(u):
    """u at r-faces → average to cell centres.  Returns shape (n_r, n_z)."""
    return 0.5 * (u[:-1, :] + u[1:, :])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Incompressible Pipe Flow Solver")
    print(f"  Re={Re}, Pr={Pr}, n_r={n_r}, n_z={n_z}")
    print(f"  dt={dt}, max steps={n_steps}")
    print("=" * 60)

    # ---- grid ---------------------------------------------------------------
    r_c, r_f, z_c, z_f, dr, dz = make_grid(n_r, n_z)

    # ---- CFL check ----------------------------------------------------------
    CFL_conv = dt / dr          # worst case w=1, CFL = dt/dz
    CFL_visc = dt * (1.0/Re) * (1.0/dr**2 + 1.0/dz**2)
    print(f"  CFL_conv ~ {max(CFL_conv, dt/dz):.4f}  (should be < 1)")
    print(f"  CFL_visc ~ {CFL_visc:.4f}  (should be < 0.5)")
    if CFL_visc > 0.5:
        print("  WARNING: viscous stability limit exceeded. Reduce dt or coarsen grid.")

    # ---- initialise fields --------------------------------------------------
    u, w, p, T = initialise_fields(n_r, n_z)
    apply_all_bc(u, w, T, dr, n_r, Re, Pr)

    # ---- build Poisson matrix (once) ----------------------------------------
    print("\nAssembling Poisson matrix ...")
    A_poisson = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z)
    print(f"  Matrix size: {A_poisson.shape}, nnz={A_poisson.nnz}")
    poisson_solve = spla.factorized(A_poisson.tocsc())

    # ---- history arrays for convergence -------------------------------------
    hist = {'continuity': [], 'momentum': [], 'temperature': []}

    # ---- time loop ----------------------------------------------------------
    print("\nStarting time integration ...")
    rhs_previous = None
    for step in range(1, n_steps + 1):

        # ---- AB2 predictor + one incremental pressure correction ------------
        u_new, w_new, p_new, T_new, rhs_previous = ab2_step(
            u, w, p, T, rhs_previous, r_c, r_f, dr, dz, dt, Re, Pr,
            n_r, n_z, poisson_solve, alpha_p)

        # ---- Convergence check ------------------------------------------------
        res = compute_residuals(u_new, u, w_new, w, T_new, T,
                                r_c, r_f, dr, dz, dt)
        for key in hist:
            hist[key].append(res[key])

        # ---- advance --------------------------------------------------------
        u, w, p, T = u_new, w_new, p_new, T_new

        # ---- progress print -------------------------------------------------
        # flush=True so progress is visible immediately when stdout is
        # redirected to a file/log (otherwise Python buffers stdout and
        # nothing appears until the process exits or the buffer fills).
        if step % output_interval == 0:
            print(f"  step {step:6d}/{n_steps}  "
                  f"R_cont={res['continuity']:.2e}  "
                  f"R_vel={res['momentum']:.2e}  "
                  f"R_T={res['temperature']:.2e}", flush=True)

        if check_convergence(res, tol_continuity, tol_velocity, tol_temperature):
            print(f"\n  Converged at step {step}!")
            break

    else:
        print("\nWARNING: maximum step count reached before all convergence "
              "tolerances were met.")

    print("\nTime integration finished.")

    # ---- post-processing: map to cell centres for plotting ------------------
    w_cc = w_to_cell_centre(w)    # (n_r, n_z)
    u_cc = u_to_cell_centre(u)    # (n_r, n_z)

    # ---- Plot 1: axial velocity contour -------------------------------------
    ax = plot_velocity_contour_axial(w_cc, r_c, z_c, u_star=u_cc,
                                     title='Axial velocity w/W_in and velocity vectors')
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'w_contour.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved w_contour.png")

    # ---- Plot 2: radial velocity contour ------------------------------------
    ax = plot_velocity_contour_radial(u_cc, r_c, z_c)
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'u_contour.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved u_contour.png")

    # ---- Plot 3: temperature contour ----------------------------------------
    ax = plot_temperature_contour(T, r_c, z_c)
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'theta_contour.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved theta_contour.png")

    # ---- Plot 4: velocity profiles at selected z-positions ------------------
    z_plot = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
    ax = plot_velocity_profile(w_cc, r_c, z_plot, z_c,
                               w_analytical=calculate_analytical_velocity_nondim)
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'w_profiles.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved w_profiles.png")

    # ---- Plot 5: temperature profiles ---------------------------------------
    def theta_an(r, z):
        return calculate_analytical_temperature_nondim(r, z, Re, Pr)

    ax = plot_temperature_profile(T, r_c, z_plot, z_c, theta_analytical=theta_an)
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'T_profiles.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved T_profiles.png")

    # ---- Plot 6: convergence history ----------------------------------------
    ax = plot_convergence_history(hist)
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'convergence.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved convergence.png")

    # ---- Plot 7: comparison with analytical at z*=50 (fully developed) ------
    z_fd_idx = np.argmin(np.abs(z_c - 50.0))
    w_num = w_cc[:, z_fd_idx]
    w_ana = calculate_analytical_velocity_nondim(r_c)
    ax = plot_comparison_with_analytical(w_num, w_ana, r_c,
                                         ylabel='w/W_in',
                                         title='Velocity at z*=50: Numerical vs Analytical')
    ax.get_figure().savefig(os.path.join(PLOTS_DIR, 'w_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close('all')
    print("  Saved w_comparison.png")

    print("\nAll plots saved to", PLOTS_DIR)
    return u, w, p, T, r_c, z_c, hist


if __name__ == '__main__':
    main()
