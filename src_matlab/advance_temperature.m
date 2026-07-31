function T_new = advance_temperature(T, u, w, r_c, r_f, dr, dz, dt, Re, Pr)
% ADVANCE_TEMPERATURE  Advance temperature by one explicit Euler step.
%   Thin wrapper around compute_rhs_T; the main time loop uses
%   compute_rhs_T directly (via rk2_step) to build the 2nd-order Heun update.
%
%   Returns T_new : (n_r, n_z)

    rhs = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr);
    T_new = T + dt * rhs;

    % Re-apply inlet BC (may be polluted by convective term at j=1)
    T_new(:, 1) = 0.0;
end
