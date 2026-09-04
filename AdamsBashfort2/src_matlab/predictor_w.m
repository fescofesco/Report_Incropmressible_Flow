function w_star = predictor_w(w, u, r_c, r_f, dr, dz, dt, Re)
% PREDICTOR_W  Advance axial velocity by one explicit Euler step (no
%   pressure). Thin wrapper around compute_rhs_w; the main time loop uses
%   compute_rhs_w directly (via ab2_step) to build the AB2 update.
%
%   Returns w_star : (n_r, n_z+1)  predicted axial velocity (BCs NOT yet reapplied)

    rhs = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re);
    n_z1 = size(w, 2);
    n_z = n_z1 - 1;
    w_star = w;
    w_star(:, 2:n_z) = w(:, 2:n_z) + dt * rhs;
end
