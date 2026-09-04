function rhs = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr)
% COMPUTE_RHS_T  Explicit RHS of the energy equation (no time
%   discretization): rhs = -conv_r - conv_z + visc_r + visc_z
%
%   Reused by both RK2 stages in rk2_step.m.
%
%   Convective fluxes use a 2nd-order upwind (linear-extrapolation / LUD)
%   scheme -- see upwind2_face.m. Diffusive terms use standard 2nd-order
%   central differencing.
%
%   Wall BC: dtheta/dr* = +Re*Pr (see Report/nusselt_number_analysis.md
%   for the derivation -- NOT -1, which drops the Re*Pr factor from the
%   non-dimensionalization).
%   Implemented as a ghost cell: theta_ghost(n_r+1,:) = theta(n_r,:) + dr*Re*Pr
%
%   T   : (n_r, n_z)    current theta (BCs already applied)
%   u   : (n_r+1, n_z)  radial velocity at r-faces
%   w   : (n_r, n_z+1)  axial velocity at z-faces
%   r_c, r_f, dr, dz, Re, Pr : as usual
%
%   Returns rhs : (n_r, n_z)

    [n_r, n_z] = size(T);
    alpha = 1.0 / (Re * Pr);
    RePr = Re * Pr;

    % ---- ghost extension in r (depth 1, used by diffusion) -----------------
    T_ghost_axis = T(1, :);
    T_ghost_wall = T(end, :) + dr * RePr;
    T_ext = [T_ghost_axis; T; T_ghost_wall];   % (n_r+2, n_z)

    % ---- ghost extension in z (depth 1, used by diffusion) -----------------
    T_ghost_in  = -T(:, 1);
    T_ghost_out =  T(:, end);
    T_zext = [T_ghost_in, T, T_ghost_out];   % (n_r, n_z+2)

    % ---- second ghost layer, needed by the 2nd-order-upwind stencil --------
    T_ghost_axis2 = T(2, :);
    T_ghost_wall2 = T(end, :) + 2.0 * dr * RePr;
    T_ext2 = [T_ghost_axis2; T_ghost_axis; T; T_ghost_wall; T_ghost_wall2];  % (n_r+4, n_z)

    T_ghost_in2  = -T(:, 2);
    T_ghost_out2 =  T(:, end);
    T_zext2 = [T_ghost_in2, T_ghost_in, T, T_ghost_out, T_ghost_out2];  % (n_r, n_z+4)

    % ---- radial convective flux (1/r) d(r*u*theta)/dr -----------------------
    u_E = u(2:end, :);
    u_I = u(1:end-1, :);

    T_r_im2 = T_ext2(1:n_r,     :);
    T_r_im1 = T_ext2(2:n_r+1,   :);
    T_r_i   = T_ext2(3:n_r+2,   :);   % = T
    T_r_ip1 = T_ext2(4:n_r+3,   :);
    T_r_ip2 = T_ext2(5:n_r+4,   :);

    T_rf_o_upwind = upwind2_face(T_r_im1, T_r_i, T_r_ip1, T_r_ip2, u_E);
    T_rf_i_upwind = upwind2_face(T_r_im2, T_r_im1, T_r_i, T_r_ip1, u_I);

    rf_o = r_f(2:end);
    rf_i = r_f(1:end-1);
    rc   = r_c;

    conv_r = (rf_o .* u_E .* T_rf_o_upwind - rf_i .* u_I .* T_rf_i_upwind) ./ (rc .* dr);

    % ---- axial convective flux d(w*theta)/dz ---------------------------------
    w_tp = w(:, 2:end);
    w_bt = w(:, 1:end-1);

    T_z_jm2 = T_zext2(:, 1:n_z);
    T_z_jm1 = T_zext2(:, 2:n_z+1);
    T_z_j   = T_zext2(:, 3:n_z+2);    % = T
    T_z_jp1 = T_zext2(:, 4:n_z+3);
    T_z_jp2 = T_zext2(:, 5:n_z+4);

    T_zf_p_upwind = upwind2_face(T_z_jm1, T_z_j, T_z_jp1, T_z_jp2, w_tp);
    T_zf_m_upwind = upwind2_face(T_z_jm2, T_z_jm1, T_z_j, T_z_jp1, w_bt);

    % At the physical inlet face the incoming transported value is the
    % prescribed Dirichlet value theta=0, not a LUD reconstruction. The odd
    % ghosts remain necessary for the second-order axial diffusion stencil.
    inflow = w_bt(:, 1) > 0.0;
    T_zf_m_upwind(inflow, 1) = 0.0;

    conv_z = (w_tp .* T_zf_p_upwind - w_bt .* T_zf_m_upwind) / dz;

    % ---- radial diffusion (1/r) d/dr(r*dtheta/dr) ---------------------------
    visc_r = alpha * (rf_o .* (T_ext(3:end, :) - T) ...
                    - rf_i .* (T - T_ext(1:end-2, :))) ./ (rc * dr^2);

    % ---- axial diffusion d^2 theta/dz^2 --------------------------------------
    visc_z = alpha * (T_zext(:, 3:end) - 2.0 * T + T_zext(:, 1:end-2)) / dz^2;

    rhs = -conv_r - conv_z + visc_r + visc_z;
end
