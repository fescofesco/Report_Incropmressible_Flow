function rhs = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re)
% COMPUTE_RHS_W  Explicit RHS of the axial-momentum equation (no pressure,
%   no time discretization): rhs = -conv_r - conv_z + visc_r + visc_z
%
%   Kept separate from time-stepping for the AB2 history update.
%
%   Convective fluxes use a 2nd-order upwind (linear-extrapolation / LUD)
%   scheme -- see upwind2_face.m -- retaining the upwind stability bias at
%   high cell Peclet number while being formally 2nd-order accurate,
%   matching the diffusive terms (2nd-order central differencing).
%
%   w  : (n_r, n_z+1)   current axial velocity (BCs already applied)
%   u  : (n_r+1, n_z)   current radial velocity (BCs already applied)
%   r_c, r_f, dr, dz, Re : as usual
%
%   Returns rhs : (n_r, n_z-1)  interior-face RHS (excludes j=1 and j=n_z+1)

    [n_r, n_z1] = size(w);
    n_z = n_z1 - 1;

    % ---- ghost extension in r for w (depth 1, used by diffusion) ----------
    w_ghost_wall = -w(end, :);
    w_ghost_axis =  w(1, :);
    w_ext = [w_ghost_axis; w; w_ghost_wall];   % (n_r+2, n_z+1)

    % ---- second ghost layer in r, needed by the 2nd-order-upwind stencil --
    w_ghost_axis2 = w(2, :);
    w_ghost_wall2 = -w(end-1, :);
    w_ext2 = [w_ghost_axis2; w_ghost_axis; w; w_ghost_wall; w_ghost_wall2];  % (n_r+4, n_z+1)
    % w_ext2(i+2, :) = w(i, :)

    % ---- interior slice j=2..n_z --------------------------------------------
    w_int = w(:, 2:n_z);   % (n_r, n_z-1)

    % ---- u interpolated to z_f(j) for j=2..n_z ------------------------------
    u_at_zf = 0.5 * (u(:, 1:end-1) + u(:, 2:end));   % (n_r+1, n_z-1)

    % --- Radial convective flux (1/r) d(r*u*w)/dr ---
    w_r_im2 = w_ext2(1:n_r,     2:n_z);
    w_r_im1 = w_ext2(2:n_r+1,   2:n_z);
    w_r_i   = w_ext2(3:n_r+2,   2:n_z);   % = w_int
    w_r_ip1 = w_ext2(4:n_r+3,   2:n_z);
    w_r_ip2 = w_ext2(5:n_r+4,   2:n_z);

    u_E = u_at_zf(2:end, :);
    w_E_upwind = upwind2_face(w_r_im1, w_r_i, w_r_ip1, w_r_ip2, u_E);

    u_I = u_at_zf(1:end-1, :);
    w_I_upwind = upwind2_face(w_r_im2, w_r_im1, w_r_i, w_r_ip1, u_I);

    % --- Axial convective flux d(w^2)/dz ---
    % w already has real values exactly at the inlet/outlet faces, so only
    % one extra ghost point beyond each end is needed (constant extension).
    w_zghost_before = w(:, 1);
    w_zghost_after  = w(:, end);
    w_zext1 = [w_zghost_before, w, w_zghost_after];   % (n_r, n_z+3)

    w_z_jm2 = w_zext1(:, 1:n_z-1);
    w_z_jm1 = w_zext1(:, 2:n_z);
    w_z_j   = w_zext1(:, 3:n_z+1);        % = w_int
    w_z_jp1 = w_zext1(:, 4:n_z+2);
    w_z_jp2 = w_zext1(:, 5:n_z+3);

    w_N = 0.5 * (w_z_j + w_z_jp1);
    w_N_upwind = upwind2_face(w_z_jm1, w_z_j, w_z_jp1, w_z_jp2, w_N);

    w_S = 0.5 * (w_z_jm1 + w_z_j);
    w_S_upwind = upwind2_face(w_z_jm2, w_z_jm1, w_z_j, w_z_jp1, w_S);

    % ---- reshaping for broadcasting -----------------------------------------
    rc  = r_c;
    rf_o = r_f(2:end);
    rf_i = r_f(1:end-1);

    % ---- convective fluxes (2nd-order upwind) --------------------------------
    conv_r = (rf_o .* u_E .* w_E_upwind - rf_i .* u_I .* w_I_upwind) ./ (rc .* dr);
    conv_z = (w_N .* w_N_upwind - w_S .* w_S_upwind) / dz;

    % ---- diffusive fluxes (central, 2nd order) ------------------------------
    visc_r = (rf_o .* (w_ext(3:end, 2:n_z) - w_int) ...
            - rf_i .* (w_int - w_ext(1:end-2, 2:n_z))) ./ (rc * dr^2) / Re;

    visc_z = (w(:, 3:n_z+1) - 2*w_int + w(:, 1:n_z-1)) / dz^2 / Re;

    rhs = -conv_r - conv_z + visc_r + visc_z;
end
