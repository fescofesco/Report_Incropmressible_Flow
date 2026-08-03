function theta_an = analytical_temperature_pdf_original(r_star, z_star, Re, Pr)
% ANALYTICAL_TEMPERATURE_PDF_ORIGINAL  Fully developed temperature profile
%   EXACTLY as literally given in the assignment PDF (page 2), for direct
%   side-by-side comparison with analytical_temperature.m (the re-derived/
%   corrected version used elsewhere in this codebase).
%
%   theta = 4*z* + Re*Pr * [(2r*)^2 - 1/4*(2r*)^4 - 3/4]
%
%   This is the formula as transcribed from the assignment; it has NOT been
%   independently re-derived (unlike analytical_temperature.m) and four
%   analytical derivations plus an empirical check against the CFD solution
%   suggest its radial-term coefficients are double what they should be --
%   see Report/temperature_formula_review.md for the full analysis. Kept
%   here so both versions can be plotted/compared directly rather than
%   silently replacing one with the other.

    eta = 2.0 * r_star;
    theta_an = 4.0 * z_star + Re * Pr * (eta.^2 - 0.25 * eta.^4 - 0.75);
end
