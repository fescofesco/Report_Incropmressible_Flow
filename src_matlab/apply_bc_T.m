function T = apply_bc_T(T)
% APPLY_BC_T  Enforce temperature boundary conditions.
%
%   T(:, 1)   = 0             inlet: theta = 0
%   T(:, end) = T(:, end-1)  outlet: zero-gradient
%
%   Wall heat flux BC (dtheta/dr = +Re*Pr) is handled via ghost cells
%   inside the temperature solver (see compute_rhs_T.m and
%   Report/temperature_formula_review.md).

    T(:, 1)   = 0.0;          % inlet
    T(:, end) = T(:, end-1);  % outlet zero-gradient
end
