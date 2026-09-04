function [u, w, p, T] = initialise_fields(n_r, n_z)
% INITIALISE_FIELDS  Allocate and set initial conditions.
%
%   [u, w, p, T] = initialise_fields(n_r, n_z)
%
%   u : (n_r+1, n_z) radial velocity at r-faces  -> zero
%   w : (n_r,   n_z+1) axial velocity at z-faces -> uniform 1
%   p : (n_r,   n_z) pressure at cell centres     -> zero
%   T : (n_r,   n_z) temperature at cell centres   -> zero

    u = zeros(n_r + 1, n_z);
    w = ones(n_r, n_z + 1);
    p = zeros(n_r, n_z);
    T = zeros(n_r, n_z);
end
