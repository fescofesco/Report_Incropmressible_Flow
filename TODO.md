# TODO

## Verify before submitting

- [ ] **Double-check the corrected fully-developed temperature formula against lecture
      notes / course materials.** The assignment PDF states
      `θ = 4z* + Re·Pr[(2r*)² - (1/4)(2r*)⁴ - 3/4]`, but the code now uses
      `θ = 4z* + Re·Pr[½(2r*)² - ⅛(2r*)⁴ - 7/48]` (half the radial-term coefficients,
      constant re-derived from energy balance). This is based on 4 independent
      analytical derivations plus 1 empirical check against the actual CFD solution
      (all agreeing with each other, none matching the literal PDF text) — see
      `Report/temperature_formula_review.md` for the full writeup. Still worth a
      sanity check against the professor's derivation/lecture notes before submitting,
      in case the PDF is correct and something else is being missed.
      Affects: `src_python/analytical.py`, `src_matlab/analytical_temperature.m`,
      `Report/main.tex` (~line 604 and the comparison profile discussion).

## Remaining work

- [ ] **Port the Python fixes to MATLAB** (`src_matlab/`): wall temperature BC
      (`+Re·Pr` instead of `-1`), RK2/Heun 2nd-order time integration (currently
      1st-order Euler), 2nd-order-upwind convective scheme (currently 1st-order
      upwind), and the corrected analytical formula once confirmed above.
- [ ] **Update `Report/main.tex`** numerical-method sections that still describe the
      old scheme (1st-order upwind, explicit Euler) — search for "first-order upwind"
      and "Explicit Euler time stepping". Also re-run the `\lstinputlisting` code
      listings once the MATLAB source is updated, since the report embeds the MATLAB
      files directly.
- [ ] **Re-run the full 30,000-step simulation** with the corrected analytical formula
      in place and regenerate `Plots_python/` (the last full run used the fixed wall
      BC + RK2 + 2nd-order upwind, but still compared against the old/uncorrected
      analytical formula in the T_profiles/comparison plots).
- [ ] Consider whether `n_steps`/`dt` need adjusting: momentum converges to machine
      precision well before step 15,000, but the temperature residual (`R_T`) never
      reaches `tol_temperature` because θ keeps rising with z at steady state (expected
      for constant wall heat flux — not a bug). Worth deciding/documenting what
      "converged" means for this problem before reporting a final answer.
