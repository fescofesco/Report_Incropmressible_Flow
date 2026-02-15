function w_star = predictor_w(w, u, r_c, r_f, dr, dz, dt, Re)
% PREDICTOR_W  Advance axial velocity by one explicit step (no pressure).
%
%   Conservation-form FVM with upwind convection and central diffusion.
%
%   w      : (n_r, n_z+1)  current axial velocity
%   u      : (n_r+1, n_z)  current radial velocity
%   r_c    : (n_r, 1)      cell-centre radial coords
%   r_f    : (n_r+1, 1)    face radial coords
%   dr, dz, dt, Re : scalars
%
%   Returns w_star : (n_r, n_z+1)

    [n_r, n_z1] = size(w);
    n_z = n_z1 - 1;

    % Ghost extension in r for w
    % Wall ghost: no-slip w=0 -> w_ghost = -w(end,:)
    % Axis ghost: symmetry dw/dr=0 -> w_ghost = w(1,:)
    w_ext = [w(1,:); w; -w(end,:)];   % (n_r+2, n_z+1)
    % w_ext(i+1, j) = w(i, j)  (1-based offset)

    % Interior slice j=2..n_z (MATLAB 1-based -> columns 2:n_z)
    w_int = w(:, 2:n_z);   % (n_r, n_z-1)

    % u interpolated to z_f(j) for j=2..n_z
    u_at_zf = 0.5 * (u(:, 1:end-1) + u(:, 2:end));   % (n_r+1, n_z-1)

    % --- Radial convective flux (1/r) d(r*u*w)/dr ---
    % Outer r-face (r_f(i+1)):
    u_E = u_at_zf(2:end, :);                                   % (n_r, n_z-1)
    w_E_upwind = w_int .* (u_E > 0) + w_ext(3:end, 2:n_z) .* (u_E <= 0);

    % Inner r-face (r_f(i)):
    u_I = u_at_zf(1:end-1, :);                                 % (n_r, n_z-1)
    w_I_upwind = w_ext(1:end-2, 2:n_z) .* (u_I > 0) + w_int .* (u_I <= 0);

    % --- Axial convective flux d(w^2)/dz ---
    % North face at z_c(j): transport velocity
    w_N = 0.5 * (w_int + w(:, 3:n_z+1));
    w_N_upwind = w_int .* (w_N > 0) + w(:, 3:n_z+1) .* (w_N <= 0);

    % South face at z_c(j-1):
    w_S = 0.5 * (w(:, 1:n_z-1) + w_int);
    w_S_upwind = w(:, 1:n_z-1) .* (w_S > 0) + w_int .* (w_S <= 0);

    % Reshape for broadcasting (column vectors)
    rc  = r_c;                    % (n_r, 1)
    rf_o = r_f(2:end);           % (n_r, 1)  r_f(i+1)
    rf_i = r_f(1:end-1);        % (n_r, 1)  r_f(i)

    % Convective fluxes (upwind)
    conv_r = (rf_o .* u_E .* w_E_upwind - rf_i .* u_I .* w_I_upwind) ./ (rc .* dr);
    conv_z = (w_N .* w_N_upwind - w_S .* w_S_upwind) / dz;

    % Diffusive fluxes (central, 2nd order)
    visc_r = (rf_o .* (w_ext(3:end, 2:n_z) - w_int) ...
            - rf_i .* (w_int - w_ext(1:end-2, 2:n_z))) ./ (rc * dr^2) / Re;

    visc_z = (w(:, 3:n_z+1) - 2*w_int + w(:, 1:n_z-1)) / dz^2 / Re;

    % Assemble predictor
    w_star = w;
    w_star(:, 2:n_z) = w_int + dt * (-conv_r - conv_z + visc_r + visc_z);
end
