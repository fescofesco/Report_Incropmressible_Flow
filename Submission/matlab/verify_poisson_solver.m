%% VERIFY_POISSON_SOLVER  Verification of the hand-written direct Poisson solver.
%
%  The assignment asks for a direct solution method that is PROGRAMMED, not a
%  library call. FACTORIZE_POISSON / SOLVE_POISSON implement the block Thomas
%  algorithm on top of our own Doolittle LU kernel (LU_FACTOR / LU_SOLVE).
%  This script checks that implementation:
%
%    [0] the dense LU kernel on a small random system;
%    [1] the backward error of the block-Thomas solution, i.e. the relative
%        residual ||A*p' - b||_inf / ||b||_inf against the independently
%        assembled sparse operator A from BUILD_POISSON_MATRIX. This is the
%        meaningful accuracy measure: the operator is stiff (cond_1(A) ~ 1e8),
%        so the *solution* difference against any other solver is bounded by
%        cond(A)*eps and says little about either solver;
%    [2] agreement with MATLAB's sparse direct solver (backslash) on the same
%        system, as an independent reference. This is the ONLY place in the
%        whole code where backslash appears -- it is a check, not the solver;
%    [3] that the Dirichlet outlet condition p' = 0 comes out as imposed.
%
%  Run:  verify_poisson_solver

clear; clc;

c = constants();
n_r = c.n_r;  n_z = c.n_z;  dt = c.dt;
[r_c, r_f, ~, ~, dr, dz] = make_grid(n_r, n_z, c.L_over_D);

fprintf('==============================================================\n');
fprintf('Verification of the hand-written direct Poisson solver\n');
fprintf('  grid: n_r = %d, n_z = %d   (%d unknowns)\n', n_r, n_z, n_r*n_z);
fprintf('  block structure: %d blocks of %d x %d\n', n_z, n_r, n_r);
fprintf('==============================================================\n');

%% [0] the dense LU kernel itself
rng(42);
M = randn(12);
rhs = randn(12, 1);
[LU, piv] = lu_factor(M);
x = lu_solve(LU, piv, rhs);
fprintf('\n[0] dense LU kernel, 12x12 random system\n');
fprintf('    ||M*x - rhs||_inf = %.3e\n', max(abs(M*x - rhs)));

%% [1] factorise and solve
t0 = tic;
P = factorize_poisson(r_c, r_f, dr, dz, n_r, n_z);
t_fac = toc(t0);
fprintf('\n[1] block-Thomas factorisation: %.1f ms (once)\n', t_fac*1e3);

% A divergence-carrying velocity predictor with the wall/axis BCs of u
rng(0);
u_star = randn(n_r+1, n_z);
u_star(1, :)   = 0.0;    % symmetry axis
u_star(end, :) = 0.0;    % no-penetration at the wall
w_star = randn(n_r, n_z+1);

n_rep = 20;
t0 = tic;
for k = 1:n_rep
    p_prime = solve_poisson(P, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z);
end
t_solve = toc(t0) / n_rep;

% Rebuild the same RHS that solve_poisson forms internally
div_u = (r_f(2:end) .* u_star(2:end, :) - r_f(1:end-1) .* u_star(1:end-1, :)) ./ (r_c * dr);
div_w = (w_star(:, 2:end) - w_star(:, 1:end-1)) / dz;
b = (div_u + div_w) / dt;
b(:, end) = 0.0;

% Flat ordering used by build_poisson_matrix: k = (i-1)*n_z + j
b_flat = reshape(b', [], 1);
p_flat = reshape(p_prime', [], 1);

A = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z);
res_own = max(abs(A*p_flat - b_flat)) / max(abs(b_flat));
fprintf('    solve: %.3f ms per Poisson solve\n', t_solve*1e3);
fprintf('    relative residual ||A p'' - b||_inf / ||b||_inf = %.3e\n', res_own);

%% [2] independent reference: MATLAB sparse backslash (check only)
p_ref_flat = A \ b_flat;
res_ref = max(abs(A*p_ref_flat - b_flat)) / max(abs(b_flat));
d_rel = max(abs(p_flat - p_ref_flat)) / max(abs(p_ref_flat));
cond_est = condest(A);
fprintf('\n[2] reference MATLAB backslash (used here only, as a check)\n');
fprintf('    relative residual, backslash                   = %.3e\n', res_ref);
fprintf('    ||p''_own - p''_backslash||_inf / ||p''||_inf      = %.3e\n', d_rel);
fprintf('    cond_1(A) estimate                             = %.3e\n', cond_est);
fprintf('    cond_1(A) * eps                                = %.3e   (bound on the line above)\n', ...
        cond_est * eps);

%% [3] boundary conditions
fprintf('\n[3] boundary conditions\n');
fprintf('    max |p''| on the outlet column (Dirichlet p''=0)  = %.3e\n', ...
        max(abs(p_prime(:, end))));

%% Verdict
fprintf('\n==============================================================\n');
if res_own < 1e-10
    fprintf('PASS - hand-written direct solver reproduces the operator\n');
else
    fprintf('FAIL - residual too large\n');
end
fprintf('==============================================================\n');
