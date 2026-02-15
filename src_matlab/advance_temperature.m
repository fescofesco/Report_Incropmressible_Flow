function T_new = advance_temperature(T, u, w, r_c, r_f, dr, dz, dt, Re, Pr)
% ADVANCE_TEMPERATURE  Advance temperature by one explicit Euler step.
%
%   Upwind convection, central diffusion.
%   Wall BC: dtheta/dr = -1 via ghost cell.
%
%   T   : (n_r, n_z)     current theta
%   u   : (n_r+1, n_z)   radial velocity
%   w   : (n_r, n_z+1)   axial velocity
%   Returns T_new : (n_r, n_z)

    [n_r, n_z] = size(T);
    alpha = 1.0 / (Re * Pr);

    % Ghost extension in r
    % Axis ghost:  dT/dr = 0 -> T_ghost = T(1,:)
    % Wall ghost:  dT/dr = -1 -> T_ghost = T(end,:) - dr
    T_ext = [T(1,:); T; T(end,:) - dr];   % (n_r+2, n_z)

    % Ghost extension in z
    % Inlet ghost:  T=0 -> ghost = -T(:,1) so avg = 0 at face
    % Outlet ghost: zero-gradient -> ghost = T(:,end)
    T_zext = [-T(:,1), T, T(:,end)];   % (n_r, n_z+2)

    % --- Radial convective flux (1/r) d(r*u*T)/dr ---
    u_E = u(2:end, :);     % u at r_f(i+1), (n_r, n_z)
    u_I = u(1:end-1, :);   % u at r_f(i),   (n_r, n_z)

    % Upwind T at outer r-face
    T_rf_o_up = T .* (u_E > 0) + T_ext(3:end, :) .* (u_E <= 0);
    % Upwind T at inner r-face
    T_rf_i_up = T_ext(1:end-2, :) .* (u_I > 0) + T .* (u_I <= 0);

    rf_o = r_f(2:end);      % (n_r, 1)
    rf_i = r_f(1:end-1);    % (n_r, 1)
    rc   = r_c;             % (n_r, 1)

    conv_r = (rf_o .* u_E .* T_rf_o_up - rf_i .* u_I .* T_rf_i_up) ./ (rc * dr);

    % --- Axial convective flux d(w*T)/dz ---
    w_tp = w(:, 2:end);    % w at z_f(j+1), (n_r, n_z)
    w_bt = w(:, 1:end-1);  % w at z_f(j),   (n_r, n_z)

    % T_zext(:, j+1) = T(:, j)
    T_south_top = T_zext(:, 2:end-1);   % T at z_c(j)
    T_north_top = T_zext(:, 3:end);     % T at z_c(j+1)
    T_zf_p_up = T_south_top .* (w_tp > 0) + T_north_top .* (w_tp <= 0);

    T_south_bot = T_zext(:, 1:end-2);   % T at z_c(j-1)
    T_north_bot = T_zext(:, 2:end-1);   % T at z_c(j)
    T_zf_m_up = T_south_bot .* (w_bt > 0) + T_north_bot .* (w_bt <= 0);

    conv_z = (w_tp .* T_zf_p_up - w_bt .* T_zf_m_up) / dz;

    % --- Radial diffusion (1/r) d/dr(r dT/dr) ---
    visc_r = alpha * (rf_o .* (T_ext(3:end,:) - T) ...
                    - rf_i .* (T - T_ext(1:end-2,:))) ./ (rc * dr^2);

    % --- Axial diffusion d^2 T/dz^2 ---
    visc_z = alpha * (T_zext(:,3:end) - 2*T + T_zext(:,1:end-2)) / dz^2;

    % Time step
    T_new = T + dt * (-conv_r - conv_z + visc_r + visc_z);

    % Re-apply inlet BC
    T_new(:, 1) = 0.0;
end
