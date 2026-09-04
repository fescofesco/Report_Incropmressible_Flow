# Assignment files

| File | What it is |
|---|---|
| `Assignment_WS25.pdf` | **Working / canonical version.** Faithful transcription of the original assignment with **one** correction in the fully-developed comparison temperature profile (see below). This is the file the report embeds via `\includepdf` and the one to read in future work. Built from `Assignment_WS25.tex`. |
| `Assignment_WS25.tex` | LaTeX source for the corrected transcription. |
| `Assignment_WS25_error.pdf` | The **original** assignment PDF exactly as distributed, kept unchanged for reference. Its fully-developed temperature formula is wrong. |
| `Assignment_WS25_blacked.pdf` | Original with personal data redacted (unchanged). |
| `Guidel_lines_for_report.pdf` | Report formatting guidelines (unchanged). |
| `Derivations_assignment` | Scratch notes on the derivations to produce (unchanged). |

## The correction

The original assignment prints the fully-developed profile as

```
(T - T_in)·ρ·W_in·c_v / q_w  =  4·z/D  +  Re·Pr·[ (2r/D)²  -  ¼(2r/D)⁴  -  ¾ ]
```

The radial terms are missing a factor `½`, and the additive constant `−¾` should be `−7/48`.
The corrected form (used everywhere in this repo) is

```
(T - T_in)·ρ·W_in·c_v / q_w  =  4·z/D  +  Re·Pr·[ ½(2r/D)²  -  ⅛(2r/D)⁴  -  7/48 ]
                             =  4·z/D  +  ½·Re·Pr·[ (2r/D)²  -  ¼(2r/D)⁴  -  7/24 ]   (same thing)
```

Why: integrating the non-dimensional energy equation for thermally fully developed flow with
constant wall heat flux gives these coefficients directly; the `−7/48` constant is fixed by the
global energy balance so that the bulk temperature `θ_bulk = 4·z/D` holds for **all** `z ≥ 0`.
The original form gives Nusselt number `24/11 ≈ 2.18` (half the classical `48/11 ≈ 4.36` for
this configuration) and violates `θ_bulk = 4z/D`.

Full derivation, the residual caveat (the analytical *radial shape* is only valid past the
thermal entry length `z ≳ 25D`; the *level* is always right), and literature sources are in
[`../Report/nusselt_number_analysis.md`](../Report/nusselt_number_analysis.md).
