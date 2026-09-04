# Review: fully-developed temperature formula — moved

> **This note has been superseded by [`nusselt_number_analysis.md`](nusselt_number_analysis.md).**

The analysis of the assignment's fully-developed temperature profile — why its radial
coefficients appear to be doubled, the from-first-principles derivation that the constant-wall-
heat-flux Nusselt number must be exactly `Nu = 48/11`, the literature sources, the Nusselt
number each candidate profile implies (`48/11` corrected vs. `24/11` for the assignment-literal
form), and the empirical cross-check against the converged CFD solution — now lives in
`Report/nusselt_number_analysis.md`, in cleaned-up form.

The filename is kept so that existing `\texttt{Report/temperature\_formula\_review.md}`
cross-references in `main.tex`, `All_derivations.tex`, and the source docstrings still resolve;
those will be repointed at the next report edit.

**Bottom line (unchanged):** the code, `main.tex`, and `All_derivations.tex` use the corrected
profile
`θ = 4z* + Re·Pr[½(2r*)² − ⅛(2r*)⁴ − 7/48]`,
which gives `Nu = 48/11 ≈ 4.364`, satisfies the global energy balance `θ_bulk = 4z*`, and is
consistent with the assignment's stated wall condition `λ∂T/∂r|_wall = q_w`. The
earlier form `[(2r*)² − ¼(2r*)⁴ − ¾]` implies `Nu = 24/11 ≈ 2.18` (exactly half) and is a
transcription error in the radial terms (missing factor `½` on `Re·Pr`, and `−¾` where the
energy balance gives `−7/48`).
