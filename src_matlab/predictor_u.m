function u_star = predictor_u(u, w, r_c, r_f, dr, dz, dt, Re)
% PREDICTOR_U  Advance radial velocity by one explicit step (no pressure).
%
%   Conservation-form FVM with upwind convection and central diffusion.
%
%   u      : (n_r+1, n_z)  current radial velocity
%   w      : (n_r, n_z+1)  current axial velocity
%   r_c    : (n_r, 1)
%   r_f    : (n_r+1, 1)
%   dr, dz, dt, Re : scalars
%
%   Returns u_star : (n_r+1, n_z)

    [~, n_z] = size(u);

    % Ghost extension in z for u
    % Inlet ghost:  u=0 at inlet -> u_ghost = -u(:,1)
    % Outlet ghost: zero-gradient -> u_ghost = u(:,end)
    u_ext = [-u(:,1), u, u(:,end)];   % (n_r+1, n_z+2)
    % u_ext(:, j+1) = u(:, j)  (1-based)

    % Interior slice i=2..n_r (rows 2:end-1 of u)
    u_int = u(2:end-1, :);   % (n_r-1, n_z)

    % --- Radial convection ---
    u_out = 0.5 * (u_int + u(3:end, :));     % transport vel at r_c(i)
    u_in  = 0.5 * (u(1:end-2, :) + u_int);   % transport vel at r_c(i-1)

    u_out_upwind = u_int .* (u_out > 0) + u(3:end, :) .* (u_out <= 0);
    u_in_upwind  = u(1:end-2, :) .* (u_in > 0) + u_int .* (u_in <= 0);

    % --- Axial convection ---
    % w interpolated to r_f(i) for i=2..n_r
    w_at_rf = 0.5 * (w(1:end-1, :) + w(2:end, :));   % (n_r-1, n_z+1)

    w_top = w_at_rf(:, 2:end);    % transport vel at z_f(j+1), (n_r-1, n_z)
    w_bot = w_at_rf(:, 1:end-1);  % transport vel at z_f(j),   (n_r-1, n_z)

    % u values for upwind in z (using u_ext, rows 2:end-1)
    u_south_top = u_ext(2:end-1, 2:end-1);   % u at z_c(j) = u_int
    u_north_top = u_ext(2:end-1, 3:end);      % u at z_c(j+1)
    u_top_upwind = u_south_top .* (w_top > 0) + u_north_top .* (w_top <= 0);

    u_south_bot = u_ext(2:end-1, 1:end-2);   % u at z_c(j-1)
    u_north_bot = u_ext(2:end-1, 2:end-1);   % u at z_c(j) = u_int
    u_bot_upwind = u_south_bot .* (w_bot > 0) + u_north_bot .* (w_bot <= 0);

    % Reshape for broadcasting
    rc_o   = r_c(2:end);       % r_c(i),   (n_r-1, 1)
    rc_i   = r_c(1:end-1);    % r_c(i-1), (n_r-1, 1)
    rf_int = r_f(2:end-1);    % r_f(i),   (n_r-1, 1)

    % Convective fluxes (upwind)
    conv_r = (rc_o .* u_out .* u_out_upwind - rc_i .* u_in .* u_in_upwind) ./ (rf_int .* dr);
    conv_z = (w_top .* u_top_upwind - w_bot .* u_bot_upwind) / dz;

    % Diffusive fluxes (central)
    visc_r = ((rc_o .* (u(3:end, :) - u_int) - rc_i .* (u_int - u(1:end-2, :))) ./ (rf_int * dr^2) ...
              - u_int ./ rf_int.^2) / Re;

    visc_z = (u_ext(2:end-1, 3:end) - 2*u_int + u_ext(2:end-1, 1:end-2)) / dz^2 / Re;

    % Assemble predictor
    u_star = u;
    u_star(2:end-1, :) = u_int + dt * (-conv_r - conv_z + visc_r + visc_z);
end
