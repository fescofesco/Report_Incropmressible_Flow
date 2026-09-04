function [u, w, T] = apply_all_bc(u, w, T)
% APPLY_ALL_BC  Apply all boundary conditions in one call.

    u = apply_bc_u(u);
    w = apply_bc_w(w);
    T = apply_bc_T(T);
end
