function theta_an = analytical_temperature(r_star, z_star, Re, Pr)
% ANALYTICAL_TEMPERATURE  Fully developed temperature profile.
%   theta = 4*z* + Re*Pr * [(2r*)^2 - 1/4*(2r*)^4 - 3/4]

    eta = 2.0 * r_star;
    theta_an = 4.0 * z_star + Re * Pr * (eta.^2 - 0.25 * eta.^4 - 0.75);
end
