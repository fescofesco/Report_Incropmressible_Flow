# Current Heun/RK2 Projection Algorithm

## Scope

The current Python and MATLAB solvers use the same explicit, second-order
Heun (explicit trapezoidal RK2) method. The velocity associated with every RK
stage is projected onto the divergence-free space. Consequently, one complete
time step requires two pressure-Poisson solves.

The velocity components are stored on a staggered cylindrical finite-volume
grid:

- radial velocity `u` is stored on radial cell faces;
- axial velocity `w` is stored on axial cell faces;
- pressure `p` and temperature `T` are stored at cell centres.

Let

\[
\mathcal{F}_u(u,w),\qquad \mathcal{F}_w(u,w)
\]

denote the spatially discretized convective and viscous momentum terms. They
do not contain a pressure gradient. Let \(\mathcal{F}_T(T,u,w)\) denote the
temperature convection-diffusion operator.

## One Heun/RK2 time step

### 1. First momentum stage

The first right-hand-side evaluations use the divergence-free fields at time
level \(n\):

\[
k_{1,u}=\mathcal{F}_u(u^n,w^n),\qquad
k_{1,w}=\mathcal{F}_w(u^n,w^n).
\]

An explicit Euler predictor is formed:

\[
u_1^*=u^n+\Delta t\,k_{1,u},\qquad
w_1^*=w^n+\Delta t\,k_{1,w}.
\]

Velocity boundary conditions are applied to this predictor.

### 2. First pressure projection

The intermediate predictor is generally not divergence-free. A projection
pressure \(\phi_1\) is obtained from

\[
\nabla^2\phi_1=\frac{1}{\Delta t}\nabla\cdot\boldsymbol{v}_1^*.
\]

The intermediate velocity is corrected:

\[
\boldsymbol{v}_1=\boldsymbol{v}_1^*-\Delta t\nabla\phi_1.
\]

Thus \(\boldsymbol{v}_1\) satisfies the discrete continuity equation before
it is used in the second RK stage.

### 3. Second momentum stage

The second right-hand-side evaluation uses the projected intermediate state:

\[
k_{2,u}=\mathcal{F}_u(u_1,w_1),\qquad
k_{2,w}=\mathcal{F}_w(u_1,w_1).
\]

The Heun predictor is then

\[
u^*=u^n+\frac{\Delta t}{2}(k_{1,u}+k_{2,u}),
\]

\[
w^*=w^n+\frac{\Delta t}{2}(k_{1,w}+k_{2,w}).
\]

### 4. Final pressure projection

A second Poisson equation is solved:

\[
\nabla^2\phi^{n+1}
=\frac{1}{\Delta t}\nabla\cdot\boldsymbol{v}^*.
\]

The final velocity is

\[
\boldsymbol{v}^{n+1}
=\boldsymbol{v}^*-\Delta t\nabla\phi^{n+1}.
\]

This second projection enforces discrete continuity at the new time level.

## Meaning of the stored pressure

The current momentum predictor omits the old pressure gradient. Therefore,
the Poisson solution \(\phi^{n+1}\) is the complete non-incremental projection
pressure for the current step, not a small correction to \(p^n\).

The stored pressure is relaxed toward this new projection pressure:

\[
p^{n+1}=(1-\alpha_p)p^n+\alpha_p\phi^{n+1}.
\]

It must not be updated with

\[
p^{n+1}=p^n+\alpha_p\phi^{n+1},
\]

because \(\phi^{n+1}\) does not tend to zero at steady state in this
non-incremental formulation. Adding it repeatedly would cause the stored
pressure to grow with time. Since the stored pressure is not included in the
next predictor, `alpha_p` affects the reported pressure history but not the
projected velocity solution.

## Temperature update

Temperature uses the same Heun weights. First,

\[
k_{1,T}=\mathcal{F}_T(T^n,u^n,w^n),\qquad
\widetilde{T}=T^n+\Delta t\,k_{1,T}.
\]

The second stage uses the same projected intermediate velocity as the
momentum equation:

\[
k_{2,T}=\mathcal{F}_T(\widetilde{T},u_1,w_1),
\]

\[
T^{n+1}=T^n+\frac{\Delta t}{2}(k_{1,T}+k_{2,T}).
\]

Thermal inlet and outlet conditions are imposed through ghost-cell values in
the temperature flux calculation; boundary-adjacent cell-centre values are
not overwritten with face boundary values.

## Properties and cost

- Explicit and second-order accurate in time when the stage projections and
  boundary conditions are applied consistently.
- Self-starting; no previous-step RHS is required.
- Requires two momentum RHS evaluations and two Poisson solves per step.
- Keeps the transport velocity divergence-free at both RK stages.
- Uses a non-incremental pressure projection, unlike the incremental
  Adams-Bashforth pressure-correction algorithm in `mitschrift.pdf`.

