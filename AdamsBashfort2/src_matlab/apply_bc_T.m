function T = apply_bc_T(T)
% APPLY_BC_T  No-op for T -- kept only so apply_all_bc has a uniform
%   interface across u, w, T.
%
%   Inlet (Dirichlet theta=0) and outlet (zero-gradient) are enforced
%   purely through ghost cells inside compute_rhs_T (T_ghost_in=-T(:,1),
%   T_ghost_out=T(:,end)) -- NOT by overwriting the T(:,1)/T(:,end) cell
%   values here. Those are cell centres a half-cell away from the actual
%   inlet/outlet faces, so clamping them directly would over-constrain
%   the PDE and is inconsistent with the ghost-cell treatment used
%   everywhere else (see compute_rhs_T.m and
%   Appendix B of the report for the wall BC derivation).
%
%   Wall heat flux BC (dtheta/dr = +Re*Pr) is likewise handled entirely
%   via ghost cells inside compute_rhs_T.m.
end
