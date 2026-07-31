# TODO

## Verify before submitting

- [ ] **Double-check the corrected fully-developed temperature formula against lecture
      notes / course materials.** The assignment PDF states
      `θ = 4z* + Re·Pr[(2r*)² - (1/4)(2r*)⁴ - 3/4]`, but the code now uses
      `θ = 4z* + Re·Pr[½(2r*)² - ⅛(2r*)⁴ - 7/48]` (half the radial-term coefficients,
      constant re-derived from energy balance). This rests on 4 independent analytical
      derivations plus a properly-redone empirical check against the CFD solution
      (`Nu ≈ 4.6` at `z*=29.9`, within ~5% of the classical 48/11≈4.36, vs. ≈2.18
      implied by the sheet's literal coefficients) — see
      `Report/temperature_formula_review.md`. Still worth a sanity check against the
      professor's derivation/lecture notes before submitting.
      Affects: `src_python/analytical.py`, `src_matlab/analytical_temperature.m`,
      `Report/main.tex` (~line 633) — all already updated to the corrected formula.

## Done

- [x] Fix wall temperature BC (`+Re·Pr` instead of `-1`) in Python and MATLAB.
- [x] Upgrade both solvers to explicit 2nd-order (Heun/RK2) time integration.
- [x] Upgrade convective discretization to 2nd-order upwind (LUD) in both solvers.
- [x] Fix analytical volumetric flow rate formula (`pi*R*^2`, was `0.5*pi*R*^2` —
      exactly half the value confirmed by numerical integration of the velocity
      profile). `src_python/analytical.py::calculate_analytical_flow_rate_nondim`.
- [x] **Fix pressure update accumulating without bound.** `update_pressure` did
      `p_new = p + zeta*p_prime` (additive), but since the momentum predictor excludes
      the pressure gradient entirely (non-incremental projection), `p_prime` as solved
      each step already IS the full physical pressure, not a vanishing correction —
      additive accumulation made `p` grow roughly linearly forever (confirmed: max|p|
      reached ~30,000 after 2000 steps). Fixed to a relaxed blend:
      `p_new = p + zeta*(p_prime - p)`. Confirmed bounded/converging after the fix
      (max|p| ~17-33 and settling). Note: this bug never actually affected `u`, `w`, or
      `T` dynamics (they only ever consumed the fresh `p_prime`, not the accumulated
      `p`) — it was a real bug in the *reported/output* pressure field, not a silent
      corruption of the flow/temperature solution.
- [x] **Fix RK2 projection to project each stage, not just the final combined step.**
      The 2nd RK stage's RHS (`k2`) was evaluated on an un-projected Euler predictor
      that didn't satisfy continuity, which is inconsistent with the conservative FVM
      formulation (assumes solenoidal transport velocity) and can silently reduce the
      scheme below 2nd-order accuracy. Redesigned `rk2_step` (Python + MATLAB) to
      project the stage-1 predictor (extra Poisson solve) before using it to evaluate
      k2, for both momentum and temperature's second stage. Verified: MATLAB still
      matches Python bit-for-bit after the redesign. Runtime roughly doubles (extra
      Poisson solve per step, ~35 min for 30,000 steps instead of ~19 min).
- [x] Investigated the claimed temperature undershoot (`theta_min` significantly below
      0, unphysical given zero inlet temp and only heating). Could not reproduce with
      the corrected pressure/RK2 code over 6000 steps (worst case ~1e-6, floating-point
      noise) — likely was a downstream symptom of the pressure/RK2 bugs rather than an
      inherent LUD monotonicity issue for this problem's parameters. Not adding a
      TVD/limiter given no empirical evidence of need after the real bugs were fixed;
      revisit if a coarser grid, larger dt, or different Re/Pr is used later.
- [x] Cross-check MATLAB RHS functions against Python to floating-point precision.
- [x] Sync `Report/main.tex`'s numerical-method description (was still describing
      1st-order upwind / explicit Euler); report recompiles cleanly.
- [x] Fixed `Report/temperature_formula_review.md` section E, which had mixed up
      wall-to-*center* and wall-to-*bulk* temperature differences (Nu is defined via
      the flow-weighted bulk/mixing-cup average, not centerline) AND checked a
      z-position (`z*=44.9`) the thermal front hadn't convectively reached yet given
      the run length used (`t*=20` vs. the ~45 time units needed) — numbers from that
      check should not have been trusted.
- [x] Re-ran the full 30,000-step simulation in Python and MATLAB with the fixed
      pressure/RK2-projection code; `Plots_python/`/`Plots_matlab/` regenerated and
      confirmed consistent with the pre-fix runs (velocity/temperature results were
      barely affected in practice for this problem — the intermediate-stage
      divergence violation was apparently small — but the reported pressure field and
      the theoretical correctness of the scheme are now both fixed).
- [x] **Properly redid the empirical Nu check** using the full run (`t*=60`), computing
      `theta_bulk` via the flow-weighted average at three z-positions with decreasing
      settling margin (`t*-z*`). Global energy balance (`theta_bulk=4z*`) holds to ~1%
      at `z*=29.9` (best margin), and `Nu` is within ~5% of 48/11 there, degrading
      toward the outlet as expected. `Report/temperature_formula_review.md` section E
      and `Report/main.tex` updated with the corrected numbers.

## Still open

- [ ] Consider whether `n_steps`/`dt` need adjusting: momentum converges to machine
      precision well before step 15,000, but the temperature residual (`R_T`) never
      reaches `tol_temperature` because θ keeps rising with z at steady state (expected
      for constant wall heat flux — not a bug). Worth deciding/documenting what
      "converged" means for this problem before reporting a final answer.
