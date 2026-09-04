function w_an = analytical_velocity(r_star)
% ANALYTICAL_VELOCITY  Hagen-Poiseuille fully developed velocity profile.
%   w* = 2 * (1 - (2*r*)^2)

    w_an = 2.0 * (1.0 - (2.0 * r_star).^2);
end
