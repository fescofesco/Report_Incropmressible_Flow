function rhs = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re)
% COMPUTE_RHS_U  Explicit RHS of the radial-momentum equation (no pressure,
%   no time discretization). Reused by the AB2 time integrator.
%
%   Convective fluxes use a 2nd-order upwind (linear-extrapolation / LUD)
%   scheme -- see upwind2_face.m.
%
%   u  : (n_r+1, n_z)
%   w  : (n_r,   n_z+1)
%   r_c, r_f, dr, dz, Re : as usual
%
%   Returns rhs : (n_r-1, n_z)  interior-face RHS (excludes i=1 and i=n_r+1)

    [n_rp1, n_z] = size(u);
    n_r = n_rp1 - 1;

    % ---- ghost extension in z for u (depth 2), needed by diffusion and by --
    % the 2nd-order-upwind stencil.
    u_ghost_in1  = -u(:, 1);
    u_ghost_in2  = -u(:, 2);
    u_ghost_out1 =  u(:, end);
    u_ghost_out2 =  u(:, end);
    u_ext = [u_ghost_in1, u, u_ghost_out1];                          % (n_r+1, n_z+2), depth 1 (diffusion)
    u_zext2 = [u_ghost_in2, u_ghost_in1, u, u_ghost_out1, u_ghost_out2];  % (n_r+1, n_z+4)

    % ---- second ghost layer in r, needed by the 2nd-order-upwind stencil ---
    u_ghost_axis_ext = -u(2, :);        % represents "u(0)" (one before axis)
    u_ghost_wall_ext = -u(end-1, :);    % represents "u(n_r+2)" (one beyond wall)
    u_ext_r = [u_ghost_axis_ext; u; u_ghost_wall_ext];   % (n_r+3, n_z)

    % ---- interior slice i=2..n_r ---------------------------------------------
    u_int = u(2:end-1, :);    % (n_r-1, n_z)

    % --- Radial convective flux (1/r) d(r*u^2)/dr ---
    u_r_im2 = u_ext_r(1:n_r-1, :);
    u_r_im1 = u_ext_r(2:n_r,   :);
    u_r_i   = u_ext_r(3:n_r+1, :);   % = u_int
    u_r_ip1 = u_ext_r(4:n_r+2, :);
    u_r_ip2 = u_ext_r(5:n_r+3, :);

    u_out = 0.5 * (u_int + u(3:end, :));
    u_in  = 0.5 * (u(1:end-2, :) + u_int);

    u_out_upwind = upwind2_face(u_r_im1, u_r_i, u_r_ip1, u_r_ip2, u_out);
    u_in_upwind  = upwind2_face(u_r_im2, u_r_im1, u_r_i, u_r_ip1, u_in);

    % ---- face velocities for axial convection ---------------------------------
    w_at_rf = 0.5 * (w(1:end-1, :) + w(2:end, :));   % (n_r-1, n_z+1)
    w_top = w_at_rf(:, 2:end);
    w_bot = w_at_rf(:, 1:end-1);

    % --- Axial convective flux d(u*w)/dz --- (2nd-order upwind, interior rows)
    u_z_jm2 = u_zext2(2:end-1, 1:n_z);
    u_z_jm1 = u_zext2(2:end-1, 2:n_z+1);
    u_z_j   = u_zext2(2:end-1, 3:n_z+2);   % = u_int
    u_z_jp1 = u_zext2(2:end-1, 4:n_z+3);
    u_z_jp2 = u_zext2(2:end-1, 5:n_z+4);

    u_top_upwind = upwind2_face(u_z_jm1, u_z_j, u_z_jp1, u_z_jp2, w_top);
    u_bot_upwind = upwind2_face(u_z_jm2, u_z_jm1, u_z_j, u_z_jp1, w_bot);

    % ---- reshaping for broadcasting -------------------------------------------
    rc_o   = r_c(2:end);
    rc_i   = r_c(1:end-1);
    rf_int = r_f(2:end-1);

    % ---- convective fluxes (2nd-order upwind) ----------------------------------
    conv_r = (rc_o .* u_out .* u_out_upwind - rc_i .* u_in .* u_in_upwind) ./ (rf_int .* dr);
    conv_z = (w_top .* u_top_upwind - w_bot .* u_bot_upwind) / dz;

    % ---- diffusive fluxes (central, 2nd order) ---------------------------------
    visc_r = ((rc_o .* (u(3:end, :) - u_int) - rc_i .* (u_int - u(1:end-2, :))) ./ (rf_int * dr^2) ...
              - u_int ./ rf_int.^2) / Re;

    visc_z = (u_ext(2:end-1, 3:end) - 2*u_int + u_ext(2:end-1, 1:end-2)) / dz^2 / Re;

    rhs = -conv_r - conv_z + visc_r + visc_z;
end
