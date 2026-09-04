%% MAIN  Simulation loop for 2D axisymmetric incompressible pipe flow.
%
%  Algorithm (per time step) -- see rk2_step.m:
%    Explicit 2nd-order (Heun / RK2) momentum predictor: TWO Poisson
%    solves/projections per step (one to make the stage-1 intermediate
%    velocity divergence-free before evaluating the 2nd RHS stage, one on
%    the final combined predictor to enforce div(u^{n+1})=0), a pressure
%    update, and a Heun step for the temperature equation.

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

%% Build Poisson matrix (once)
fprintf('\nAssembling Poisson matrix ...\n');
A_poisson = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z);
fprintf('  Matrix size: %d x %d, nnz = %d\n', size(A_poisson,1), size(A_poisson,2), nnz(A_poisson));
poisson_solver = decomposition(A_poisson, 'lu');

%% History arrays
hist_cont = zeros(n_steps, 1);
hist_vel  = zeros(n_steps, 1);
hist_temp = zeros(n_steps, 1);

%% Time loop
fprintf('\nStarting time integration ...\n');
converged_step = 0;

for step = 1:n_steps

    % 2nd-order explicit (Heun/RK2) step: momentum + Poisson + temperature
    [u_new, w_new, p_new, T_new] = rk2_step(u, w, p, T, r_c, r_f, dr, dz, dt, Re, Pr, n_r, n_z, poisson_solver, alpha_p);

    % Convergence check
    res = compute_residuals(u_new, u, w_new, w, T_new, T, r_c, r_f, dr, dz, dt);
    hist_cont(step) = res.continuity;
    hist_vel(step)  = res.momentum;
    hist_temp(step) = res.temperature;

    % Advance
    u = u_new;  w = w_new;  p = p_new;  T = T_new;

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
figure('Visible', 'off');
contourf(z_c, r_c, w_cc, 20, 'LineStyle', 'none');
colorbar; xlabel('z/D'); ylabel('r/D');
title('Axial velocity w/W_{in}');
saveas(gcf, fullfile(plots_dir, 'w_contour.png'));
fprintf('  Saved w_contour.png\n');

%% Plot 2: radial velocity contour
figure('Visible', 'off');
contourf(z_c, r_c, u_cc, 20, 'LineStyle', 'none');
colorbar; colormap(gca, 'cool'); xlabel('z/D'); ylabel('r/D');
title('Radial velocity u/W_{in}');
saveas(gcf, fullfile(plots_dir, 'u_contour.png'));
fprintf('  Saved u_contour.png\n');

%% Plot 3: temperature contour
figure('Visible', 'off');
contourf(z_c, r_c, T, 20, 'LineStyle', 'none');
colorbar; colormap(gca, 'hot'); xlabel('z/D'); ylabel('r/D');
title('Non-dimensional temperature \theta');
saveas(gcf, fullfile(plots_dir, 'theta_contour.png'));
fprintf('  Saved theta_contour.png\n');

%% Plot 4: velocity profiles at selected z-positions
z_plot = [5, 10, 20, 30, 40, 50];
figure('Visible', 'off'); hold on;
w_an = analytical_velocity(r_c);
for k = 1:length(z_plot)
    [~, jj] = min(abs(z_c - z_plot(k)));
    plot(r_c, w_cc(:, jj), '-o', 'MarkerSize', 3, 'DisplayName', sprintf('z*=%.0f', z_plot(k)));
end
plot(r_c, w_an, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('r/D'); ylabel('w/W_{in}');
title('Velocity profiles at selected z-positions');
legend('Location', 'best'); grid on;
saveas(gcf, fullfile(plots_dir, 'w_profiles.png'));
fprintf('  Saved w_profiles.png\n');

%% Plot 5: temperature profiles
figure('Visible', 'off'); hold on;
for k = 1:length(z_plot)
    [~, jj] = min(abs(z_c - z_plot(k)));
    plot(r_c, T(:, jj), '-o', 'MarkerSize', 3, 'DisplayName', sprintf('z*=%.0f', z_plot(k)));
end
% Analytical fully-developed profile at z*=50 (approved form; see
% Report/nusselt_number_analysis.md)
[~, j_Tfd] = min(abs(z_c - z_plot(end)));
z_Tfd = z_c(j_Tfd);
theta_an = analytical_temperature(r_c, z_Tfd, Re, Pr);
plot(r_c, theta_an, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('z*=%.1f', z_Tfd));
xlabel('r/D'); ylabel('\theta');
title('Temperature profiles at selected z-positions');
legend('Location', 'best'); grid on;
saveas(gcf, fullfile(plots_dir, 'T_profiles.png'));
fprintf('  Saved T_profiles.png\n');

%% Plot 6: convergence history
figure('Visible', 'off');
semilogy(1:converged_step, hist_cont, 'b-', 'DisplayName', 'Continuity');
hold on;
semilogy(1:converged_step, hist_vel, 'r-', 'DisplayName', 'Momentum');
semilogy(1:converged_step, hist_temp, 'g-', 'DisplayName', 'Temperature');
xlabel('Time step'); ylabel('Residual');
title('Convergence history');
legend('Location', 'best'); grid on;
saveas(gcf, fullfile(plots_dir, 'convergence.png'));
fprintf('  Saved convergence.png\n');

%% Plot 7: comparison with analytical at z*=50
[~, j_fd] = min(abs(z_c - 50.0));
w_num = w_cc(:, j_fd);
w_ana = analytical_velocity(r_c);

figure('Visible', 'off');
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
saveas(gcf, fullfile(plots_dir, 'w_comparison.png'));
fprintf('  Saved w_comparison.png\n');

fprintf('\nAll plots saved to %s\n', plots_dir);
close all;
