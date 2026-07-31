function u_star = predictor_u(u, w, r_c, r_f, dr, dz, dt, Re)
% PREDICTOR_U  Advance radial velocity by one explicit Euler step (no
%   pressure). Thin wrapper around compute_rhs_u; the main time loop uses
%   compute_rhs_u directly (via rk2_step) to build the 2nd-order Heun update.
%
%   Returns u_star : (n_r+1, n_z)

    rhs = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re);
    u_star = u;
    u_star(2:end-1, :) = u(2:end-1, :) + dt * rhs;
end
