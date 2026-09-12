function T_new = advance_temperature(T, u, w, r_c, r_f, dr, dz, dt, Re, Pr)
% ADVANCE_TEMPERATURE  Advance temperature by one explicit Euler step.
%   Thin wrapper around compute_rhs_T; the main time loop uses
%   compute_rhs_T directly (via ab2_step) to build the AB2 update.
%
%   Returns T_new : (n_r, n_z)

    rhs = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr);
    T_new = T + dt * rhs;

    % Inlet/outlet are enforced via ghost cells inside compute_rhs_T; no
    % direct cell-value overwrite is needed.
end
