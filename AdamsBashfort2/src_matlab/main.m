%% MAIN  Simulation loop for 2D axisymmetric incompressible pipe flow.
%
%  Second-order Adams-Bashforth integration with one incremental pressure
%  correction per step. Forward Euler initialises the multistep history.

clear; close all; clc;

%% Load parameters
c = constants();
Re = c.Re;  Pr = c.Pr;
n_r = c.n_r;  n_z = c.n_z;
dt = c.dt;  n_steps = c.n_steps;
alpha_p = c.alpha_p;

fprintf('============================================================\n');
fprintf('Incompressible Pipe Flow Solver (MATLAB)\n');
fprintf('  Re=%g, Pr=%g, n_r=%d, n_z=%d\n', Re, Pr, n_r, n_z);
fprintf('  dt=%g, max steps=%d\n', dt, n_steps);
fprintf('============================================================\n');

%% Grid
[r_c, r_f, z_c, z_f, dr, dz] = make_grid(n_r, n_z, c.L_over_D);

% CFL check
CFL_conv = max(dt/dr, dt/dz);
CFL_visc = dt * (1/Re) * (1/dr^2 + 1/dz^2);
fprintf('  CFL_conv = %.4f  (should be < 1)\n', CFL_conv);
fprintf('  CFL_visc = %.4f  (should be < 0.5)\n', CFL_visc);
if CFL_visc > 0.5
    warning('Viscous stability limit exceeded. Reduce dt or coarsen grid.');
end

%% Initialise fields
[u, w, p, T] = initialise_fields(n_r, n_z);
[u, w, T] = apply_all_bc(u, w, T);

%% Factorise the Poisson operator (once)
% Hand-written DIRECT solver: block-Thomas elimination of the
% block-tridiagonal pressure operator, built on our own dense LU kernel
% (lu_factor.m / lu_solve.m). No MATLAB backslash, decomposition() or inv()
% appears anywhere in the solve path -- backslash is used only inside
% verify_poisson_solver.m, as an independent check. The operator is constant
% in time, so it is factorised once here; each time step then costs only the
% two sweeps in solve_poisson.m.
fprintf('\nFactorising Poisson operator (block Thomas, direct) ...\n');
t_fac = tic;
poisson_fac = factorize_poisson(r_c, r_f, dr, dz, n_r, n_z);
fprintf('  %d blocks of %d x %d  ->  %d unknowns  (%.2f s)\n', ...
        n_z, n_r, n_r, n_r*n_z, toc(t_fac));

%% History arrays
hist_cont = zeros(n_steps, 1);
hist_vel  = zeros(n_steps, 1);
hist_temp = zeros(n_steps, 1);

%% Time loop
fprintf('\nStarting time integration ...\n');
converged_step = 0;
rhs_u_old = [];
rhs_w_old = [];
rhs_T_old = [];

for step = 1:n_steps

    % AB2 predictor + one incremental pressure correction
    [u_new, w_new, p_new, T_new, rhs_u, rhs_w, rhs_T] = ...
        ab2_step(u, w, p, T, rhs_u_old, rhs_w_old, rhs_T_old, ...
                 r_c, r_f, dr, dz, dt, Re, Pr, n_r, n_z, ...
                 poisson_fac, alpha_p);

    % Convergence check
    res = compute_residuals(u_new, u, w_new, w, T_new, T, r_c, r_f, dr, dz, dt);
    hist_cont(step) = res.continuity;
    hist_vel(step)  = res.momentum;
    hist_temp(step) = res.temperature;

    % Advance
    u = u_new;  w = w_new;  p = p_new;  T = T_new;
    rhs_u_old = rhs_u;
    rhs_w_old = rhs_w;
    rhs_T_old = rhs_T;

    % Progress print
    if mod(step, c.output_interval) == 0
        fprintf('  step %6d/%d  R_cont=%.2e  R_vel=%.2e  R_T=%.2e\n', ...
                step, n_steps, res.continuity, res.momentum, res.temperature);
    end

    if check_convergence(res, c.tol_continuity, c.tol_velocity, c.tol_temperature)
        fprintf('\n  Converged at step %d!\n', step);
        converged_step = step;
        break;
    end
end

if converged_step == 0
    converged_step = n_steps;
    warning('Maximum step count reached before all convergence tolerances were met.');
end
fprintf('\nTime integration finished.\n');

% Trim history
hist_cont = hist_cont(1:converged_step);
hist_vel  = hist_vel(1:converged_step);
hist_temp = hist_temp(1:converged_step);

%% Post-processing: map to cell centres for plotting
w_cc = 0.5 * (w(:, 1:end-1) + w(:, 2:end));   % (n_r, n_z)
u_cc = 0.5 * (u(1:end-1, :) + u(2:end, :));    % (n_r, n_z)

%% Output directory
plots_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Plots_matlab');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

%% Plot 1: axial velocity contour
light_figure();
contourf(z_c, r_c, w_cc, 20, 'LineStyle', 'none');
hold on;
% Decimated quiver layer so individual arrows stay readable; the underlying
% contour still contains every cell.
%
% Two corrections are needed because the axes are stretched ~100:1
% (z/D spans 50, r/D spans 0.5):
%
%  1. MATLAB scales quiver components in DATA units. A physically small
%     radial velocity (u/w ~ 0.06 at most) then covers a large fraction of
%     the short r-axis and the arrows render almost vertically, which is
%     visually wrong for a flow that is almost entirely axial. Multiplying
%     the radial component by the axis aspect ratio makes the arrow's
%     *on-screen* angle equal to the true atan(u/w).
%
%  2. Arrowhead size is likewise measured in data units, so the default
%     head puts wings tens of percent of the r-range wide -- they appear as
%     long diagonal streaks across the plot. The head is therefore scaled
%     down by the same aspect ratio.
z_span = 50.0;
r_span = 0.5;
aspect = r_span / z_span;

z_arrow = 1:8:length(z_c);
r_arrow = 1:3:length(r_c);
[Z_a, R_a] = meshgrid(z_c(z_arrow), r_c(r_arrow));

arrow_len = 1.6;                              % z/D length of the fastest arrow
a_scale = arrow_len / max(abs(w_cc(:)));

q = quiver(Z_a, R_a, ...
           w_cc(r_arrow, z_arrow) * a_scale, ...
           u_cc(r_arrow, z_arrow) * a_scale * aspect, ...
           0, 'k', 'LineWidth', 0.5);          % 0 -> autoscaling OFF
q.Clipping = 'on';
q.MaxHeadSize = 0.003;
xlim([0, 50]);
ylim([0, 0.5]);
cb = colorbar; cb.Label.String = 'w/W_{in}';
xlabel('z/D'); ylabel('r/D');
title('Axial velocity w/W_{in} and velocity vectors');
exportgraphics(gcf, fullfile(plots_dir, 'w_contour.png'), 'Resolution', 200);
fprintf('  Saved w_contour.png\n');

%% Plot 2: radial velocity contour
% u is ~0 over most of the domain (flow is fully developed for z* > ~5,
% i.e. 90% of the pipe length) except for a strong, spatially localised
% entrance-corner value near (r=0.5, z=0). A linear color scale spanning
% the raw min/max is dominated by that single feature and renders the
% entire developed region as one flat color. Instead use a diverging
% colormap centred at zero with a robust (98th-percentile) symmetric
% limit, so the near-zero bulk keeps visible contrast and the entrance
% feature simply saturates the color scale rather than washing it out.
light_figure();
u_abs_sorted = sort(abs(u_cc(:)));
u_lim = u_abs_sorted(max(1, round(0.98 * numel(u_abs_sorted))));
if u_lim <= 0
    u_lim = max(abs(u_cc(:)));
end
if u_lim <= 0
    u_lim = 1;
end
contourf(z_c, r_c, u_cc, linspace(-u_lim, u_lim, 21), 'LineStyle', 'none');
clim([-u_lim, u_lim]);
colorbar; colormap(gca, diverging_bwr(256)); xlabel('z/D'); ylabel('r/D');
title('Radial velocity u/W_{in}');
exportgraphics(gcf, fullfile(plots_dir, 'u_contour.png'), 'Resolution', 200);
fprintf('  Saved u_contour.png\n');

%% Plot 3: temperature contour
light_figure();
contourf(z_c, r_c, T, 20, 'LineStyle', 'none');
colorbar; colormap(gca, 'hot'); xlabel('z/D'); ylabel('r/D');
title('Non-dimensional temperature \theta');
exportgraphics(gcf, fullfile(plots_dir, 'theta_contour.png'), 'Resolution', 200);
fprintf('  Saved theta_contour.png\n');

%% Plot 4: velocity profiles at selected z-positions
z_plot = [5, 10, 20, 30, 40, 50];
light_figure(); hold on;
w_an = analytical_velocity(r_c);
for k = 1:length(z_plot)
    [~, jj] = min(abs(z_c - z_plot(k)));
    plot(r_c, w_cc(:, jj), '-o', 'MarkerSize', 3, 'DisplayName', sprintf('z*=%.0f', z_plot(k)));
end
plot(r_c, w_an, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('r/D'); ylabel('w/W_{in}');
title('Velocity profiles at selected z-positions');
legend('Location', 'best'); grid on;
exportgraphics(gcf, fullfile(plots_dir, 'w_profiles.png'), 'Resolution', 200);
fprintf('  Saved w_profiles.png\n');

%% Plot 5: temperature profiles
light_figure(); hold on;
for k = 1:length(z_plot)
    [~, jj] = min(abs(z_c - z_plot(k)));
    plot(r_c, T(:, jj), '-o', 'MarkerSize', 3, 'DisplayName', sprintf('z*=%.0f', z_plot(k)));
end
% Analytical fully-developed profile at z*=50 (approved form; see
% report, Appendix B)
[~, j_Tfd] = min(abs(z_c - z_plot(end)));
z_Tfd = z_c(j_Tfd);
theta_an = analytical_temperature(r_c, z_Tfd, Re, Pr);
plot(r_c, theta_an, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Analytical profile at z*=%.1f', z_Tfd));
xlabel('r/D'); ylabel('\theta');
title('Temperature profiles at selected z-positions');
legend('Location', 'best'); grid on;
exportgraphics(gcf, fullfile(plots_dir, 'T_profiles.png'), 'Resolution', 200);
fprintf('  Saved T_profiles.png\n');

%% Plot 6: convergence history
light_figure();
semilogy(1:converged_step, hist_cont, 'b-', 'DisplayName', 'Continuity');
hold on;
semilogy(1:converged_step, hist_vel, 'r-', 'DisplayName', 'Momentum');
semilogy(1:converged_step, hist_temp, 'g-', 'DisplayName', 'Temperature');
xlabel('Time step'); ylabel('Residual');
title('Convergence history');
legend('Location', 'best'); grid on;
exportgraphics(gcf, fullfile(plots_dir, 'convergence.png'), 'Resolution', 200);
fprintf('  Saved convergence.png\n');

%% Plot 7: comparison with analytical at z*=50
[~, j_fd] = min(abs(z_c - 50.0));
w_num = w_cc(:, j_fd);
w_ana = analytical_velocity(r_c);

light_figure();
subplot(2,1,1);
plot(r_c, w_num, 'bo-', 'MarkerSize', 4, 'DisplayName', 'Numerical');
hold on;
plot(r_c, w_ana, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('r/D'); ylabel('w/W_{in}');
title('Velocity at z*=50: Numerical vs Analytical');
legend('Location', 'best'); grid on;

subplot(2,1,2);
plot(r_c, abs(w_num - w_ana), 'k-', 'LineWidth', 1);
xlabel('r/D'); ylabel('|Error|');
title('Absolute error'); grid on;
exportgraphics(gcf, fullfile(plots_dir, 'w_comparison.png'), 'Resolution', 200);
fprintf('  Saved w_comparison.png\n');

fprintf('\nAll plots saved to %s\n', plots_dir);

%% Save final fields for offline inspection / diagnostics
% Save next to THIS script. The path is derived from mfilename rather than
% from plots_dir, because the containing folder is called src_matlab in the
% development tree but matlab in the submission bundle -- hard-coding
% '../src_matlab' made this line error out there, after the plots were
% already written.
save(fullfile(fileparts(mfilename('fullpath')), 'last_run.mat'), ...
     'u', 'w', 'p', 'T', 'u_cc', 'w_cc', 'r_c', 'z_c', 'hist_cont', 'hist_vel', 'hist_temp');

close all;
