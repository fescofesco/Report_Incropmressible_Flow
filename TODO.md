# TODO

## Verify before submitting

- [x] **Corrected fully-developed temperature formula — CONFIRMED.** The original assignment
      PDF stated `θ = 4z* + Re·Pr[(2r*)² - (1/4)(2r*)⁴ - 3/4]`; the code/report use
      `θ = 4z* + Re·Pr[½(2r*)² - ⅛(2r*)⁴ - 7/48]` (factor `½` on the radial terms, constant
      re-derived from the global energy balance so `θ_bulk = 4z*` for all `z* ≥ 0`). The
      course instructor confirmed by email (2026-08-31) that the assignment formula is wrong
      on both counts and that the corrected form is to be used; the still-present caveat is
      only that the analytical *radial shape* is not valid inside the thermal entry length
      (`z* ≲ 25`) — inherent to any fully-developed reference, not fixable. Full derivation,
      literature sources (Shah & London 1978; Incropera & DeWitt; Kays & Crawford), and the
      arithmetic for both candidate profiles are in Appendix B of `Report/main.pdf`
      (Temperature Formula and Nusselt Number Audit). New report subsection
      `Validity of the Fully-Developed Temperature Profile` (`\label{sec:fd-validity}`)
      explains bulk-level-exact vs. shape-only-asymptotic. `.tex` files carry no instructor
      attribution — the correction is presented as a from-first-principles derivation.
      The original PDF is kept as `Assignemnt/Assignment_WS25_error.pdf`; the canonical
      `Assignemnt/Assignment_WS25.pdf` is now a corrected transcription (the version the
      report embeds via `\includepdf`).

- Two-variant ("assignment-form") report: **dropped.** The instructor said to use the
  corrected form, so there is no reason to also build a report around the erroneous
  `Nu = 24/11` profile. Plan file `~/.claude/plans/do-you-see-the-resilient-neumann.md` is
  obsolete.

## Known limitations (deliberately not fixed this pass — documented, not silent)

- [ ] **Pressure outlet Dirichlet reference is at the last cell centre, not the outlet
      face.** `build_poisson_matrix` (Python) / MATLAB equivalent replaces the entire
      last-column equation with `p'=0`, which is a half-cell offset from the true outlet
      face. Behaves well on the assignment grid (`max|div|` ~1e-13 to 1e-15), but this is
      grid-sensitive: a coarser test grid showed outlet divergence ~1.6e-5 instead — the
      formulation is not discretely consistent as the grid refines, it just happens to be
      small enough on the assignment grid to not matter. A proper fix would use a
      ghost-cell/face-centred Dirichlet treatment analogous to the temperature/velocity
      BCs. Not attempted this pass — pressure is a gauge/reference quantity that doesn't
      feed back into the velocity/temperature dynamics here (only `p_prime`, freshly
      solved each step, does), so the practical impact is low, but it's a real
      discretization inconsistency worth fixing if the grid or grading rubric changes.
- [ ] **2nd-order upwind (LUD) convection is unbounded and can undershoot below zero on
      coarser grids/larger `dt`.** On the assignment grid (`n_r=50, n_z=250, dt=0.002`)
      the worst observed undershoot is negligible (~-7e-6, floating-point-noise scale, see
      the "Done" entry below) — but this is *specific to that grid/timestep*, not a
      general property of the scheme. A coarser grid was reported to show a much larger,
      non-negligible undershoot (~-0.70). A bounded/TVD-limited version of `upwind2_face`
      (e.g. minmod or similar limiter) would fix this generally; not implemented this pass
      since it isn't needed for the submitted grid, but don't assume the current scheme is
      safe if the grid, `dt`, or Re/Pr change.

## Done

- [x] **Wrote the Nusselt number analysis (now Appendix B of `Report/main.pdf`)**: the
      from-first-principles derivation that constant-wall-heat-flux fully-developed laminar
      pipe flow has `Nu = 48/11` exactly (linear BVP, no free parameters), literature sources
      (Shah & London 1978; Incropera & DeWitt; Kays & Crawford), the full derivation that the
      assignment-literal coefficients instead give `Nu = 24/11` (= ½·48/11, and violates
      `θ_bulk = 4z*`), the wall-slope contrast (`Re·Pr` vs `2·Re·Pr`), and the empirical
      CFD table. Arithmetic verified numerically + symbolically. The former standalone
      Markdown notes were merged into the report appendix and deleted. Expanded the Nusselt subsection (with derivation + `\cite`s +
      a new `thebibliography`) in `Report/main.tex` and `Report/All_derivations.tex`; both
      recompile clean (main.pdf 49 pp, All_derivations.pdf 40 pp, no undefined refs/cites).
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
- [x] **Fix temperature inlet/outlet BC to use ghost cells only, no direct cell-value
      overwrite.** `rk2_step`/`advance_temperature` were forcibly setting
      `T[:,0]=0`/`T[:,-1]=T[:,-2]` after every stage, even though `compute_rhs_T` already
      handles the Dirichlet inlet / zero-gradient outlet correctly via ghost cells
      (`T_ghost_in=-T[:,0]`, `T_ghost_out=T[:,-1]`). Since `T[:,0]`/`T[:,-1]` are cell
      centres a half-cell away from the actual boundary faces, the direct overwrite was
      over-constraining the PDE and suppressing real physics — e.g. it was forcing the
      *entire* first z-column (including cells right next to the heated wall) to exactly
      zero, when physically the near-wall region right at the inlet legitimately heats up
      quickly (wall flux applies from z=0 onward; near-wall axial velocity is ~0 due to
      no-slip, so local radial conduction dominates there). Fixed in `rk2_step.py/.m` and
      `advance_temperature.py/.m`; also removed the now-redundant/inconsistent overwrite
      from `apply_bc_T` (Python + MATLAB), which is now a no-op kept only for interface
      uniformity with `apply_all_bc`. Verified stable over 3000+ steps; MATLAB still
      matches Python bit-for-bit.
- [x] Investigated the claimed temperature undershoot on the assignment grid (`theta_min`
      significantly below 0, unphysical given zero inlet temp and only heating), on the
      assignment grid. **Correction from an earlier pass of this TODO**: sampling only
      every 100-500 steps missed a brief early transient — a fine-grained scan (every
      step) reproduces a real undershoot of ~-0.0053 during steps ~1-50 (`t*`~0.002-0.1)
      right at the inlet/wall corner (`r*`~0.47, `z*`~0.1), decaying to negligible
      (~-1e-5) by step ~200 and staying negligible for the remaining ~29,800 steps. This
      is a genuine LUD-overshoot artifact from the sharp initial transient at that corner,
      but it is small, localized, and transient (<1% of the simulated time span, doesn't
      recur or grow) — not treated as requiring a limiter for this submission, but this is
      the accurate, properly-measured magnitude (not "floating-point noise" as an earlier
      version of this entry claimed). See the "Known limitations" section above for the
      separate, more serious caveat that this does **not** generalize to other grids —
      confirmed non-negligible (~-0.70) on a coarser grid, and possibly not transient
      there either (not separately re-checked).
- [x] Cross-check MATLAB RHS functions against Python to floating-point precision.
- [x] Synced numerical-method descriptions across `Report/main.tex`, `Report/All_derivations.tex`,
      `src_python/main.py`, `src_matlab/main.m` (was still describing 1st-order upwind /
      explicit Euler / a single Poisson solve per step — the code does 2nd-order upwind,
      Heun/RK2, and *two* Poisson solves per step). Also fixed stale formula/BC text in
      `README.md`, `Report/README.md`, `src_python/analytical.py` docstrings, and
      `Report/All_derivations.tex` (wall BC sign/magnitude, temperature formula
      coefficients, centreline `u` BC wording). Report recompiles cleanly.
- [x] Fixed the empirical Nu check (former temperature-formula note), which had mixed up
      wall-to-*center* and wall-to-*bulk* temperature differences (Nu is defined via
      the flow-weighted bulk/mixing-cup average, not centerline) AND checked a
      z-position (`z*=44.9`) the thermal front hadn't convectively reached yet given
      the run length used (`t*=20` vs. the ~45 time units needed) — numbers from that
      check should not have been trusted. Also resolved the stale caveat about the 7/48
      constant "not being independently derived" (it since was, from the exact global
      energy balance requiring `θ_bulk=4z*`).
- [x] **Properly redid the empirical Nu check** using the full run (`t*=60`), computing
      `theta_bulk` via the flow-weighted average at three z-positions with decreasing
      settling margin (`t*-z*`). Global energy balance (`theta_bulk=4z*`) holds to ~1%
      at `z*=29.9` (best margin), and `Nu` is within ~5% of 48/11 there, degrading
      toward the outlet as expected. `Report/main.tex` (Section 5.4) updated with the
      corrected numbers.
- [x] **Recomputed the empirical Nu check from the converged field** (`t*=139.728`,
      step 69,864, `R_T<1e-6`). The `t*=60` table had been left in the report after the run
      was extended. Re-running reproduced the old values at `t*=60` (they used θ at the last
      cell centre as θ_wall), so the new table uses the wall-face value
      `θ_nr + ½·dr·Re·Pr` instead. Converged: `θ_bulk = 4z*` to 0.01 % at `z*=29.9/39.9/48.9`;
      Nu = 4.455 / 4.397 / 4.378, decreasing toward 48/11 (+0.3 % at the outlet), as expected
      for thermal entrance-region behaviour.
- [x] **Corrected the temperature-convergence reasoning.** The earlier claim ("R_T can't
      reach tolerance because θ keeps rising with z at steady state") is simply wrong — a
      steady field can have `dθ/dz≠0` while `dθ/dt=0`; spatial gradient does not prevent
      temporal convergence. The real reason: extrapolating the logged `R_T` decay curve
      (still monotonically decreasing at `t*=60`, local decay timescale ~26 non-dim time
      units) shows reaching the `1e-6` tolerance would need ~3.7 more hours of compute
      (thermal diffusion across the pipe radius has a timescale ~`Re·Pr·R*²`≈125, and the
      run only covers `t*=60`). Decision: report the actual achieved residuals honestly
      (`R_cont`, `R_vel` at machine precision; `R_T` still decaying, not at tolerance)
      rather than chase the multi-hour run or claim false convergence.
- [x] Added a note near the T_profiles figure explaining that the analytical (dashed)
      curves are the fully-developed asymptote and are expected to diverge — including
      going negative — at small z*, since the formula is only valid once the flow has had
      time to become thermally developed. Confirmed this isn't a bug: present in both the
      original PDF-literal formula (worse, e.g. -355 at z*=4.9) and the corrected one
      (-53 at z*=4.9), and `theta_contour.png` (numerical-only) correctly shows no
      negative values since the actual PDE solution is bounded below by zero given the
      zero inlet condition and heating-only wall BC.
- [x] **Fixed `plot_temperature_profile` (Python) to only overlay the analytical curve
      once, at the largest z-position (z*=50), instead of at every selected z.** MATLAB's
      `main.m` already only did this (`analytical_temperature(r_c, 50, ...)`, single black
      dashed curve) — Python's version was plotting the analytical (fully-developed-only)
      formula at every z-position including z*=4.9 and 9.9, which is what produced the
      confusing negative dashed curves discussed above. Now both languages present the
      same, less misleading plot: numerical profiles at all positions, one analytical
      reference curve at z*=50 only.
- [x] **Removed the obsolete erroneous-PDF comparison from the final plotting path.**
      The approved assignment now contains the corrected profile, so both Python and MATLAB
      plot that single approved fully-developed curve at the last actual scalar location
      (`z*=49.9`) against the numerical solution. The historical discrepancy remains documented
      in Appendix B of `Report/main.pdf`.

## Still open

- [x] Re-ran the full 30,000-step simulation in Python and MATLAB after the exact inlet
      convective-temperature boundary fix. Both runs give the same residual history
      (`R_T=2.19` at `t*=60`) and regenerated `Plots_python/`/`Plots_matlab/` using the
      approved assignment profile at the actual last scalar location `z*=49.9`.
