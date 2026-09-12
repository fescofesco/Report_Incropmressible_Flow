function [r_c, r_f, z_c, z_f, dr, dz] = make_grid(n_r, n_z, L_over_D)
% MAKE_GRID  Generate staggered grid for axisymmetric pipe flow.
%
%   [r_c, r_f, z_c, z_f, dr, dz] = make_grid(n_r, n_z, L_over_D)
%
%   Outputs (all non-dimensional, scaled by D):
%     r_c : (n_r x 1) cell-centre radial coords   [0 < r < 0.5]
%     r_f : (n_r+1 x 1) r-face coords             [0 <= r <= 0.5]
%     z_c : (n_z x 1) cell-centre axial coords
%     z_f : (n_z+1 x 1) z-face coords
%     dr  : radial spacing
%     dz  : axial spacing

    dr = 0.5 / n_r;
    dz = L_over_D / n_z;

    r_f = linspace(0, 0.5, n_r + 1)';       % column vector
    r_c = 0.5 * (r_f(1:end-1) + r_f(2:end));

    z_f = linspace(0, L_over_D, n_z + 1)';
    z_c = 0.5 * (z_f(1:end-1) + z_f(2:end));
end
