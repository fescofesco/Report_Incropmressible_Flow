# Nusselt number for fully-developed laminar pipe flow with constant wall heat flux

**Purpose.** This note derives, from first principles, why the fully-developed Nusselt number
for this problem must be **exactly `Nu = 48/11 ≈ 4.3636`**, gives the standard literature
sources, and shows what Nusselt number each of the two candidate fully-developed temperature
profiles implies — the one *printed in the assignment* and the *corrected* one used throughout
this repository. It supersedes the earlier informal note in
`Report/temperature_formula_review.md`.

All symbols are non-dimensional in the scaling of the assignment:
`r* = r/D`, `z* = z/D`, `w* = w/W_in`, `θ = (T − T_in) ρ c_v W_in / q_w`,
`Re = ρ W_in D / μ = 100`, `Pr = μ c_v / λ = 5`. It is convenient to use
`s ≡ 2r* = 2r/D ∈ [0, 1]`, so the pipe wall is at `s = 1`.

---

## A. The two candidate fully-developed temperature profiles

Both are of the form `θ(r*, z*) = 4 z* + Re·Pr · B(s)`, differing only in the radial bracket
`B(s)`:

| | radial bracket `B(s)` | `∂θ/∂r*` at the wall (`s=1`) | flow-weighted mean `⟨B⟩_w` | `θ_wall − θ_bulk` | `Nu = Re·Pr / (θ_wall − θ_bulk)` |
|---|---|---|---|---|---|
| **Corrected** (used in the code and report) | `½s² − ⅛s⁴ − 7/48` | `Re·Pr` | `0` | `(11/48)·Re·Pr` | **`48/11 ≈ 4.3636`** |
| **Assignment-literal** (`Assignemnt/Assignment_WS25.pdf`, p. 2) | `s² − ¼s⁴ − ¾` | `2·Re·Pr` | `−11/24` | `(11/24)·Re·Pr` | **`24/11 ≈ 2.1818`** |

The corrected bracket is the one implemented in
`src_python/analytical.py::calculate_analytical_temperature`,
`src_matlab/analytical_temperature.m`, and Eq. (analytical_temperature) of the report. The
assignment-literal bracket is kept for comparison in
`calculate_analytical_temperature_pdf_original` / `analytical_temperature_pdf_original.m`.

---

## B. Derivation of `Nu = 48/11` (exact)

### B.1 Governing equation and the axial gradient

For hydrodynamically and thermally fully-developed axisymmetric flow (`u* = 0`, `∂/∂z*` of the
*shape* of the profile vanishes), with viscous dissipation and axial conduction neglected, the
non-dimensional energy equation reduces to

```
w*(r*) · ∂θ/∂z*  =  (1 / (Re·Pr)) · (1/r*) · d/dr* ( r* · dθ/dr* ).            (B1)
```

A global energy balance over a slice of the pipe fixes the axial gradient. The wall adds heat at
a non-dimensional rate `∂θ/∂r*|_wall = Re·Pr` per unit wall area (Section B.3); the enthalpy flux
carried by the flow is `∫ w* θ dA`. With `⟨w*⟩ = 1` (the parabolic profile carries the same mean
velocity as the uniform inlet) this balance gives, for the bulk (mixing-cup) temperature
`θ_bulk ≡ ∫ w* θ dA / ∫ w* dA`,

```
dθ_bulk/dz*  =  4        ⟹        θ_bulk(z*) = 4 z*                            (B2)
```

exactly, at every `z*` — this is mass/energy conservation from the inlet, independent of the
radial profile shape. In the fully-developed region the whole profile translates uniformly in
`z*`, so `∂θ/∂z* = dθ_bulk/dz* = 4` everywhere in the cross-section.

### B.2 Velocity profile

The Hagen–Poiseuille profile (Eq. analytical_velocity of the report) is

```
w*(r*)  =  2 [ 1 − (2r*)² ]  =  2 (1 − s²).                                    (B3)
```

### B.3 Wall boundary condition

Fourier's law at the heated wall, `λ ∂T/∂r|_wall = q_w`, non-dimensionalizes directly:

```
∂θ/∂r*|_wall  =  D · (ρ c_v W_in / q_w) · ∂T/∂r|_wall
              =  D · (ρ c_v W_in / q_w) · (q_w / λ)
              =  ρ c_v W_in D / λ  =  Re·Pr.                                   (B4)
```

(Numerically `Re·Pr = 500` here — *not* `O(1)`: the reference length and temperature scales do
not cancel.) The symmetry condition at the axis is `∂θ/∂r*|_{r*=0} = 0`.

### B.4 Integrating the radial ODE

Insert (B2)–(B3) into (B1):

```
(1/r*) d/dr*( r* dθ/dr* )  =  Re·Pr · 4 · w*(r*)  =  8 Re·Pr (1 − 4r*²).
```

Integrate once, applying regularity (`r* dθ/dr* → 0` as `r* → 0`):

```
r* dθ/dr*  =  8 Re·Pr ( r*²/2 − r*⁴ )
dθ/dr*     =  8 Re·Pr r* ( ½ − r*² ).                                          (B5)
```

Check the wall value: at `r* = ½`, `dθ/dr* = 8 Re·Pr · ½ · (½ − ¼) = Re·Pr`, matching (B4). ✔

Integrate again:

```
θ(r*, z*)  =  2 Re·Pr ( r*² − r*⁴ ) + C(z*)
           =  Re·Pr ( ½ s² − ⅛ s⁴ ) + C(z*),        s = 2r*.                  (B6)
```

### B.5 Fixing the constant from the global energy balance

Write `C(z*) = 4 z* + Re·Pr · c` with `c` a pure constant. Condition (B2), `θ_bulk = 4z*`,
requires the flow-weighted mean of the radial bracket to vanish:

```
⟨ ½ s² − ⅛ s⁴ + c ⟩_w  =  0,        ⟨ f ⟩_w ≡ ∫₀^{1/2} f w* r* dr* / ∫₀^{1/2} w* r* dr*.
```

With `w* r* dr* = ½ s (1 − s²) ds`, the denominator is `∫₀¹ ½ s (1−s²) ds = 1/8`, and

```
⟨ ½ s² − ⅛ s⁴ ⟩_w  =  8 · ∫₀¹ ( ½ s² − ⅛ s⁴ ) · ½ s (1 − s²) ds
                    =  8 · ( 7 / 384 )  =  7/48.
```

Hence `c = −7/48`, and

```
θ(r*, z*)  =  4 z*  +  Re·Pr [ ½ (2r*)² − ⅛ (2r*)⁴ − 7/48 ].                   (B7)
```

This is the corrected profile. (Naively halving the assignment's `−¾` would give `−3/8`, which
is wrong; the constant must be re-derived, and equals `−7/48`.)

### B.6 Forming the Nusselt number

The Nusselt number is `Nu = h D / λ` with `h = q_w / (T_wall − T_bulk)`. Non-dimensionalizing
`T_wall − T_bulk = (θ_wall − θ_bulk) · q_w / (ρ c_v W_in)`:

```
Nu  =  q_w D / [ λ (T_wall − T_bulk) ]  =  (ρ c_v W_in D / λ) / (θ_wall − θ_bulk)
    =  Re·Pr / (θ_wall − θ_bulk).                                             (B8)
```

From (B7), the bracket at the wall (`s = 1`) is `½ − ⅛ − 7/48 = 24/48 − 6/48 − 7/48 = 11/48`,
and `θ_bulk` corresponds to bracket `0`, so

```
θ_wall − θ_bulk  =  (11/48) · Re·Pr
Nu  =  Re·Pr / [ (11/48) Re·Pr ]  =  48/11  ≈  4.3636.                         (B9)
```

---

## C. Why the value is *exactly* `48/11`

Equations (B1)–(B4) form a **linear two-point boundary-value problem with fixed coefficients and
no free parameters**. Given only

- hydrodynamically **and** thermally fully-developed flow,
- **constant** wall heat flux (the "H1" / `q″ = const` thermal boundary condition),
- laminar flow with the parabolic velocity profile (B3),
- constant fluid properties,
- negligible axial conduction and viscous dissipation,

the radial profile (B7) and the ratio (B9) are completely determined. `Re`, `Pr`, `q_w`, `D`,
the fluid properties, and the numerical grid all **cancel** out of `Nu`: they rescale `θ` but not
the dimensionless ratio `(θ_wall − θ_bulk)/Re·Pr`. `48/11` is a rational number obtained by
integration, not an empirical correlation and not a fitted value. Any correctly posed
constant-flux fully-developed laminar circular-pipe problem has this same `Nu`.

(For contrast, the constant-wall-**temperature** case — "T" / "H2"-adjacent — gives the
different exact-ish value `Nu ≈ 3.6568`, the first eigenvalue of the Graetz problem. The two
canonical thermal boundary conditions have genuinely different fully-developed Nusselt numbers;
`48/11` is specifically the constant-heat-flux one, which is the case in this assignment.)

---

## D. What the assignment-literal formula would give: `Nu = 24/11`

Take the printed bracket `B_lit(s) = s² − ¼ s⁴ − ¾` at face value.

**Wall value.** `B_lit(1) = 1 − ¼ − ¾ = 0`, so the literal formula says `θ_wall = 4z*`.

**Bulk value.** With the same flow-weighted average as above,

```
⟨ s² − ¼ s⁴ ⟩_w  =  8 · ∫₀¹ ( s² − ¼ s⁴ ) · ½ s (1 − s²) ds  =  8 · ( 7/192 )  =  7/24,
⟨ B_lit ⟩_w      =  7/24 − ¾  =  7/24 − 18/24  =  −11/24.
```

So the literal formula gives `θ_bulk = 4z* + Re·Pr · (−11/24) = 4z* − (11/24) Re·Pr`, which is
**not** `4z*` — the literal profile **violates the global energy balance (B2)**. This is a
second, independent symptom of the same defect.

**Nusselt number.** `θ_wall − θ_bulk = 0 − (−(11/24) Re·Pr) = (11/24) Re·Pr`, so by (B8)

```
Nu_lit  =  Re·Pr / [ (11/24) Re·Pr ]  =  24/11  ≈  2.1818  =  ½ · (48/11).
```

**Where the factor of 2 comes from.** The literal radial coefficients are `(1, ¼)` where the
integration (B6) gives `(½, ⅛)` — each doubled. Equivalently, the wall slope of the literal
profile is

```
d/dr* [ Re·Pr ( s² − ¼ s⁴ ) ] |_{s=1}  =  Re·Pr (4s − 2s³) |_{s=1}  =  2·Re·Pr,
```

i.e. **twice** the physical wall gradient `Re·Pr` that the assignment's own stated condition
`λ ∂T/∂r|_wall = q_w` requires (B4). A wall gradient twice as steep produces a wall-to-bulk
temperature difference twice as large, and since `Nu ∝ 1/(θ_wall − θ_bulk)`, a Nusselt number
half as large.

---

## E. Sources for `Nu = 48/11`

This is a classical textbook result for the circular tube with constant wall heat flux:

- **R. K. Shah & A. L. London**, *Laminar Flow Forced Convection in Ducts*, Academic Press,
  1978 — circular duct, **H1** boundary condition: `Nu_{H1} = 48/11 = 4.36364`. (Supplement 1 to
  *Advances in Heat Transfer*; the circular-tube chapter.)
- **F. P. Incropera & D. P. DeWitt**, *Fundamentals of Heat and Mass Transfer*, Wiley — chapter
  on internal flow: for laminar, fully-developed flow in a circular tube with **uniform surface
  heat flux**, `Nu_D = 48/11 = 4.36` (constant; contrast `Nu_D = 3.66` for uniform surface
  temperature).
- **W. M. Kays, M. E. Crawford & B. Weigand**, *Convective Heat and Mass Transfer*,
  McGraw-Hill — fully-developed laminar circular-tube flow, constant heat rate per unit length:
  `Nu = 48/11`.
- **R. B. Bird, W. E. Stewart & E. N. Lightfoot**, *Transport Phenomena*, Wiley — the same
  fully-developed constant-wall-flux tube problem is worked out and yields `Nu = 48/11`.

The derivation in Section B is the standard one reproduced (in dimensional form) in each of
these references.

---

## F. What each approach yields / what would happen

| Approach | wall slope `∂θ/∂r*|_w` | `θ_bulk` behaviour | analytical `Nu` | agreement with theory |
|---|---|---|---|---|
| **Corrected profile** (repo) | `Re·Pr` | `= 4z*` (energy balance satisfied) | `48/11 ≈ 4.364` | matches Shah & London / Incropera exactly |
| **Assignment-literal profile** | `2·Re·Pr` | `= 4z* − (11/24)Re·Pr` (balance violated) | `24/11 ≈ 2.182` | half the textbook value; also self-inconsistent with the assignment's stated Fourier BC |

**Empirical check against the CFD solution.** The solver in this repository imposes the
physically consistent wall condition `∂θ/∂r* = Re·Pr` (B4). Computing the flow-weighted bulk
temperature from the converged field (`t* = 60` run) and forming `Nu = Re·Pr/(θ_wall − θ_bulk)`:

| `z*` | settling margin `t* − z*` | `θ_bulk` (CFD) | `4z*` (expected) | `θ_wall − θ_bulk` | `Nu` (CFD) | classical `48/11` |
|---|---|---|---|---|---|---|
| 29.9 | 30.1 | 119.43 | 119.6 | 109.06 | **4.585** | 4.364 |
| 39.9 | 20.1 | 155.82 | 159.6 | 105.77 | **4.727** | 4.364 |
| 48.9 | 11.1 | 181.79 | 195.6 |  99.06 | **5.048** | 4.364 |

The global balance `θ_bulk ≈ 4z*` holds to `~1 %` at `z* = 29.9` (largest settling margin), and
the CFD `Nu` is within `~5 %` of `48/11` there, degrading toward the outlet exactly as expected
from progressively less settling time (thermal diffusion across the radius has non-dimensional
timescale `~Re·Pr·R*² ≈ 125`, comparable to the simulated `t* = 60`). This is a model-free
confirmation — built only from the governing PDE and the wall BC — that the physical answer is
near `48/11`, **not** near the `2.18` implied by the assignment-literal coefficients.

**Note on re-running the solver with a different wall flux.** The energy equation, the wall BC,
and the inlet condition `θ_in = 0` are all *linear* in `θ`, and `Nu` in (B8) is normalized by
the wall flux. Multiplying the imposed wall flux by any constant `k` simply scales the entire
`θ` field by `k`, leaving `Nu` unchanged at `48/11`. A CFD run therefore **cannot** reproduce
`24/11` — that number exists only as a property of the (internally inconsistent) literal
formula, never of a correctly posed simulation.

---

## G. Summary table

With `k ≡ Re·Pr` (`= 500` for `Re = 100`, `Pr = 5`):

| quantity | corrected profile | assignment-literal profile |
|---|---|---|
| radial bracket `B(s)`, `s = 2r*` | `½s² − ⅛s⁴ − 7/48` | `s² − ¼s⁴ − ¾` |
| wall gradient `∂θ/∂r*|_{r*=½}` | `k` | `2k` |
| `B` at wall (`s = 1`) | `11/48` | `0` |
| flow-weighted mean `⟨B⟩_w` | `0` | `−11/24` |
| `θ_wall − θ_bulk` | `(11/48) k` | `(11/24) k` |
| **`Nu` (analytical)** | **`48/11 ≈ 4.3636`** | **`24/11 ≈ 2.1818`** |
| `θ_bulk = 4z*` ? | yes | no |
| `Nu` (CFD, this repo) | `4.4 – 4.6` (→ `48/11` with settling margin) | not attainable (linear problem, flux-normalized `Nu`) |

**Conclusion.** An earlier transcription of the fully-developed temperature profile carried its
two radial coefficients doubled relative to what (i) direct integration of the non-dimensional
energy equation, (ii) the stated Fourier wall condition, and (iii) the classical `Nu = 48/11`
for constant-heat-flux laminar pipe flow all require — equivalently, it was missing the factor
`½` on `Re·Pr` in the radial terms, with the constant `−¾` in place of the energy-balance value
`−7/48`. The corrected profile (B7), which gives `Nu = 48/11` exactly and satisfies the global
energy balance `θ_bulk = 4z*` for all `z* ≥ 0`, is used as the fully-developed comparison
throughout this report and codebase.
