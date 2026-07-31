function theta_an = analytical_temperature(r_star, z_star, Re, Pr)
% ANALYTICAL_TEMPERATURE  Fully developed temperature profile.
%   theta = 4*z* + Re*Pr * [1/2*(2r*)^2 - 1/8*(2r*)^4 - 7/48]
%
%   NOTE: the assignment PDF states coefficients (2r*)^2 - 1/4*(2r*)^4 - 3/4
%   (double the radial-term coefficients used here). See
%   Report/temperature_formula_review.md for the re-derivation (4 independent
%   analytical methods + 1 empirical check against the CFD solution) showing
%   these halved coefficients are correct; TODO.md flags this for a final
%   sanity check against lecture notes before submission.

    eta = 2.0 * r_star;
    theta_an = 4.0 * z_star + Re * Pr * (0.5 * eta.^2 - 0.125 * eta.^4 - 7.0/48.0);
end
