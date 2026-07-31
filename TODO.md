# TODO

## Verify before submitting

- [ ] **Double-check the corrected fully-developed temperature formula against lecture
      notes / course materials.** The assignment PDF states
      `θ = 4z* + Re·Pr[(2r*)² - (1/4)(2r*)⁴ - 3/4]`, but the code now uses
      `θ = 4z* + Re·Pr[½(2r*)² - ⅛(2r*)⁴ - 7/48]` (half the radial-term coefficients,
      constant re-derived from energy balance). This is based on 4 independent
      analytical derivations plus 2 empirical checks against the actual CFD solution in
      both Python and MATLAB (all agreeing with each other, none matching the literal
      PDF text) — see `Report/temperature_formula_review.md` for the full writeup.
      Still worth a sanity check against the professor's derivation/lecture notes before
      submitting, in case the PDF is correct and something else is being missed.
      Affects: `src_python/analytical.py`, `src_matlab/analytical_temperature.m`,
      `Report/main.tex` (~line 633 and the comparison-profile discussion) — all already
      updated to the corrected formula.

## Done

- [x] Fix wall temperature BC (`+Re·Pr` instead of `-1`) in Python and MATLAB.
- [x] Upgrade both solvers to explicit 2nd-order (Heun/RK2) time integration.
- [x] Upgrade convective discretization to 2nd-order upwind (LUD) in both solvers.
- [x] Cross-check MATLAB RHS functions against Python to floating-point precision, and
      a 300-step MATLAB stability run matches the Python smoke test exactly.
- [x] Re-run the full 30,000-step simulation in both Python and MATLAB with the
      corrected analytical formula; `Plots_python/` and `Plots_matlab/` regenerated.
      Velocity matches the analytical Hagen–Poiseuille profile to ~1e-3 at z*=50 in
      both. Temperature profiles now converge toward the analytical asymptote with
      increasing z*, as the assignment says they should.
- [x] Sync `Report/main.tex`'s numerical-method description (was still describing
      1st-order upwind / explicit Euler); report recompiles cleanly.

## Still open

- [ ] Consider whether `n_steps`/`dt` need adjusting: momentum converges to machine
      precision well before step 15,000, but the temperature residual (`R_T`) never
      reaches `tol_temperature` because θ keeps rising with z at steady state (expected
      for constant wall heat flux — not a bug). Worth deciding/documenting what
      "converged" means for this problem before reporting a final answer.
