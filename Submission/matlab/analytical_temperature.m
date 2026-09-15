function theta_an = analytical_temperature(r_star, z_star, Re, Pr)
% ANALYTICAL_TEMPERATURE  Fully developed temperature profile (approved form).
%   theta = 4*z* + Re*Pr * [1/2*(2r*)^2 - 1/8*(2r*)^4 - 7/48]
%
%   Obtained by direct integration of the non-dimensional energy equation with
%   the wall condition dtheta/dr* = Re*Pr; the constant -7/48 follows from the
%   exact global energy balance (theta_bulk = 4*z* for all z* >= 0). Gives
%   Nu = 48/11. Full derivation and sources: report, Appendix B.
%   (An earlier transcription of the assignment omitted the 1/2 factor on the
%   radial terms and used -3/4; that form is not used anywhere in this codebase.)

    eta = 2.0 * r_star;
    theta_an = 4.0 * z_star + Re * Pr * (0.5 * eta.^2 - 0.125 * eta.^4 - 7.0/48.0);
end
