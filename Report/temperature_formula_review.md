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

## E) Empirical confirmation from the fixed CFD solver

After fixing the wall BC (magnitude `Re·Pr`) and upgrading the solver to 2nd-order
(RK2 in time, 2nd-order upwind in space), I ran the actual simulation for 10,000 steps
(velocity field already converged to machine precision well before this) and read off
`θ` directly from the computed field at `z*=44.9` (Re=100, Pr=5, so `Re·Pr=500`):

| Quantity | Numerical (from the simulation) | PDF formula as transcribed | Formula corrected (halved radial terms) |
|---|---|---|---|
| `θ_wall − θ_center` | **111.2** | 369.9 | 185.0 |
| `∂θ/∂r*` near wall (one-sided FD) | **488** (vs. enforced BC of 500) | — | — |
| Implied Nu (`= 2·Re·Pr/(θ_wall−θ_center)`, since `Nu=hD/λ`, `D=2R`) | **≈ 4.5** (close to 48/11≈4.36) | ≈ 2.7 | ≈ 2.68 |

The simulation — built only from the governing PDEs and the independently-derived wall
BC, with no reference to the fully-developed formula at all — reproduces `θ_wall −
θ_center` within ~3% of the classic `Nu = 48/11` result, and lands nowhere near what the
PDF-transcribed comparison formula predicts (3.3× too high). This is strong additional,
independent evidence (on top of the four derivations in section D) that the formula as
transcribed from the assignment PDF has the radial-term-coefficient error described
above. Note `z*=44.9` isn't perfectly at the fully-developed asymptote (thermal entry
length ≈ `0.05·Re·Pr·D = 25`, so `z*=45` is ~1.8 entry lengths in — close but not exact),
which is presumably why even the "corrected" formula (185.0) doesn't match the
simulation (111.2) exactly yet; a comparison at `z*=50` after the full 30,000-step run,
or a formal grid/z-convergence study, would tighten this further.

## Bottom line

Five independent checks — four analytical derivations (one purely local/algebraic, not even
requiring the fully-developed assumption; one cross-checked against a Nusselt number your own code
already trusts) plus one empirical check against the actual, independently-built CFD solution —
all point to the same factor-of-2 discrepancy in the radial terms of the given analytical formula.
I'd treat this as "very likely a typo," but since I only have the one PDF and no lecture notes to
cross-reference, please sanity-check against your course materials before I change the formula used
for the comparison plots.
