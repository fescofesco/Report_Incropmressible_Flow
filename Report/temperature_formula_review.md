# Review: fully-developed temperature formula (assignment comparison profile)

## A) The formula in question

The assignment (Task 2 / comparison profiles, page 2 of `Assignemnt/Assignment_WS25.pdf`) gives
a fully-developed analytical temperature profile to be plotted alongside the numerical result.
As currently transcribed into the codebase — `src_python/analytical.py::calculate_analytical_temperature`,
`src_matlab/analytical_temperature.m`, and `Report/main.tex:604` — it reads:

```
θ(r*, z*) = 4·z* + Re·Pr · [ (2r*)² − (1/4)(2r*)⁴ − 3/4 ]
```

where `r* = r/D`, `z* = z/D`, and `θ = (T − T_in)·ρ·c_v·W_in/q_w` (the same non-dimensional
temperature used throughout the rest of the code).

## B) Why I think it's wrong

I believe the two radial (`r*`-dependent) terms should each carry **half** the coefficient
currently shown, i.e. the physically/mathematically consistent formula is:

```
θ(r*, z*) = 4·z* + Re·Pr · [ ½(2r*)² − (1/8)(2r*)⁴ − const ]
```

(the constant offset also needs re-deriving — see caveat at the end).

The `4·z*` term is **not** in question — it matches the standard integral energy balance and I
reproduced it independently.

The discrepancy is specifically in the coefficients multiplying `(2r*)²` and `(2r*)⁴`. Section D
below shows three independent derivations that all agree with each other and disagree with the
formula as transcribed by a consistent factor of 2 — and a fourth cross-check against a Nusselt
number your own codebase already hardcodes.

## C) What I read from the assignment PDF, and where

I read `Assignemnt/Assignment_WS25.pdf` directly (both pages) using Claude's PDF reader (which
renders/extracts the document, including the equation images/typesetting on page 2). The relevant
passage, page 2, under "Fügen Sie zum Vergleich die Profile für die voll entwickelte zylindrische
Rohrströmung hinzu, hier gegeben als":

```
w = 2 W_in [ 1 − (2r/D)² ] ,

(T − T_in)ϱ W_in c_v / q_w = 4 z/D + Re Pr [ (2r/D)² − (1/4)(2r/D)⁴ − 3/4 ] ,

welche für zunehmende Lauflänge z schlussendlich erreicht werden müssen.
```

This is an unambiguous match to what's in `analytical.py` / `analytical_temperature.m` / `main.tex` —
the code correctly transcribes what the PDF shows. **I do not have independent access to your
lecture notes or a second copy of the assignment**, so I cannot rule out that this is a typo
*in the assignment itself* (problem sheets do occasionally have transcription errors from
lecture-note derivations). I only have the one PDF you gave me, read through an automated PDF
renderer — I can't 100% rule out an extraction artifact on a stacked fraction/exponent, but the
text above rendered cleanly and consistently, so I don't think that's the issue here.

## D) What I did to determine it's likely wrong

I derived the fully-developed θ(r*) profile four independent ways. All four agree with each other
and disagree with the transcribed PDF formula by the same factor of 2 in the radial terms.

**Setup common to all derivations:** velocity profile as given, `w*(r*) = 2[1 − (2r*)²]`
(not in question), non-dimensional energy equation (also derivable independently, not in question):

```
w* ∂θ/∂z* = (1/(Re·Pr)) · (1/r*) ∂/∂r*(r* ∂θ/∂r*)
```

**Derivation 1 — direct BC non-dimensionalization (no fully-developed assumption needed).**
Fourier's law at the wall gives (dimensionally) `∂T/∂r|_wall = q_w/λ`. Converting to non-dimensional
variables (`r* = r/D`, `θ` as defined above):

```
∂θ/∂r*|_wall = (ρ c_v W_in / q_w) · D · (q_w/λ) = ρ c_v W_in D / λ = Re·Pr
```

This alone already shows the wall gradient magnitude must be `Re·Pr` (= 500 for Re=100, Pr=5),
not `O(1)`. This is a completely separate, simpler check that doesn't even depend on the velocity
profile.

**Derivation 2 — solving the fully-developed ODE in `r*` directly.**
Substituting `w*(r*)` and `∂θ/∂z* = 4` into the energy equation and integrating with the
regularity condition at `r*=0`:

```
dθ/dr* = Re·Pr · (4r* − 8r*³)
θ(r*)  = Re·Pr · (2r*² − 2r*⁴) + f(z*)
```

i.e. coefficients `2` and `−2`, **half** of the `4` and `−4` implied by `(2r*)² − ¼(2r*)⁴ = 4r*² − 4r*⁴`
in the transcribed formula.

**Derivation 3 — same ODE, solved in terms of `s = 2r*` instead (independent algebra as a
cross-check against a transcription slip in Derivation 2).** Converting the governing ODE to `s`
and integrating gives:

```
dθ/ds = Re·Pr · (s − s³/2)   →   dθ/dr* = Re·Pr·(2s − s³)
```

At the wall (`s=1`): `dθ/dr*|_wall = Re·Pr`, confirming Derivation 1 exactly. Integrating fully
and converting back to `r*` reproduces the same `2r*² − 2r*⁴` radial dependence as Derivation 2.

**Derivation 4 — cross-check against the Nusselt number your own code already uses.**
`analytical.py::calculate_nusselt_number_fully_developed()` returns the classic textbook result
`Nu = 48/11 ≈ 4.364` for laminar pipe flow with constant wall heat flux (this is a well-established
number, e.g. Kays & Crawford / Shah & London). I computed the flow-weighted (bulk) average of the
transcribed formula's radial bracket and got `θ_wall − θ_bulk = (11/24)·Re·Pr`, which corresponds to
`Nu = 24/11 ≈ 2.18` — exactly half of `48/11`. Using the corrected (halved) radial coefficients from
Derivations 2/3 instead gives `θ_wall − θ_bulk = (11/48)·Re·Pr`, which reproduces `Nu = 48/11`
exactly, consistent with the value already hardcoded elsewhere in your codebase.

**Caveat on the constant term.** I have *not* independently re-derived the `−3/4` (or, under the
correction, `−const`) offset from first principles matching the entrance condition at `z*=0` — that
requires matching the fully-developed asymptote to the developing solution, which is a longer
derivation I didn't attempt. If you correct the `(2r*)²`/`(2r*)⁴` coefficients, the constant should
be re-derived or checked against your lecture notes too, rather than assumed to just halve along
with the other terms.

## E) Empirical confirmation from the fixed CFD solver — CORRECTED

**An earlier version of this section had two errors, caught in review, both now fixed:**

1. **Wrong quantity.** It compared `θ_wall − θ_center` (centerline) against the Nu
   correlation, but the standard definition of the bulk/mean temperature used in
   `Nu = h·D/λ`, `h = q_w/(T_wall − T_bulk)` is the flow-weighted (mixing-cup) average
   over the cross-section, `θ_bulk = ∫w·θ·r dr / ∫w·r dr`, which is *not* the same as
   the centerline value. The correct check is `Nu = Re·Pr/(θ_wall − θ_bulk)` (derived
   in section D) — there is no factor of 2 in this formula.

2. **Wrong run length for the z-position checked.** The mean non-dimensional velocity
   is exactly 1, so a fluid parcel takes roughly `z*` non-dimensional time units to
   convect from the inlet to axial position `z*`. The original check ran only 10,000
   steps (`t*=20`) but read off `θ` at `z*=44.9` (needs `t*≈45` just to be reached at
   all) — the field there was still dominated by the initial condition heating in
   place, not by flow that had actually come from the inlet.

**Corrected check**, using the full 30,000-step run (`t*=60`) and the properly-defined
`θ_bulk`, at three axial positions with decreasing settling margin (`t* − z*`, the time
elapsed since the thermal front passed that location):

| `z*` | margin `t*−z*` | `θ_bulk` (numerical) | `4z*` (expected) | `θ_wall−θ_bulk` | `Nu` | classic `Nu=48/11` |
|---|---|---|---|---|---|---|
| 29.9 | 30.1 | 118.46 | 119.6 | 109.04 | **4.585** | 4.364 |
| 39.9 | 20.1 | 154.98 | 159.6 | 105.90 | **4.722** | 4.364 |
| 48.9 | 11.1 | 181.14 | 195.6 | 99.25  | **5.038** | 4.364 |

The global energy balance (`θ_bulk = 4z*` exactly, independent of profile shape) holds
to ~1% at `z*=29.9` where the flow has had the most settling time, and the numerical
`Nu` is within ~5% of the classical 4.364 there — both degrading somewhat closer to the
outlet, exactly as expected from having progressively less settling margin. This is now
a methodologically sound, independent confirmation (built only from the governing PDEs
and the wall BC, with no reference to the comparison formula) that lands close to
`Nu=48/11`, not the `≈2.18` implied by the PDF-literal (unhalved) coefficients — consistent
with, and now properly supporting, the four analytical derivations in section D.

## Bottom line

Four independent analytical derivations (one purely local/algebraic, not even requiring the
fully-developed assumption; one cross-checked against a Nusselt number your own code already
trusts), plus a properly-redone empirical check against the actual CFD solution (section E,
`Nu` within ~5% of the classical 4.364 with adequate settling time, versus ~2.18 implied by the
PDF-literal coefficients), all point to the same factor-of-2 discrepancy in the radial terms of
the given analytical formula. I'd treat this as "very likely a typo," but since I only have the
one PDF and no lecture notes to cross-reference, please sanity-check against your course
materials before treating this as final.
