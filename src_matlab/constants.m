function c = constants()
% CONSTANTS  Physical constants and simulation parameters.
%   c = constants() returns a struct with all parameters.

    % Non-dimensional numbers
    c.Re = 100.0;       % Reynolds number
    c.Pr = 5.0;         % Prandtl number

    % Geometry (non-dimensional)
    c.L_over_D = 50.0;  % Pipe length in diameters

    % Grid parameters
    c.n_r = 50;         % Radial grid cells
    c.n_z = 250;        % Axial grid cells

    % Time parameters
    c.dt = 0.002;       % Time step
    c.n_steps = 30000;  % Max time steps
    c.output_interval = 1000;

    % Numerical parameters
    c.tol_continuity = 1e-6;
    c.tol_velocity   = 1e-6;
    c.tol_temperature = 1e-6;
    c.alpha_p = 0.7;    % Pressure under-relaxation
end
