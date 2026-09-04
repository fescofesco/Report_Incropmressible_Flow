function [u_new, w_new] = correct_velocity(u_star, w_star, p_prime, dr, dz, dt)
% CORRECT_VELOCITY  Apply pressure correction to obtain div-free velocities.
%
%   u^{n+1}(i,j) = u*(i,j) - dt * (p'(i,j) - p'(i-1,j)) / dr
%   w^{n+1}(i,j) = w*(i,j) - dt * (p'(i,j) - p'(i,j-1)) / dz

    u_new = u_star;
    w_new = w_star;

    % Radial velocity at interior r-faces (rows 2:end-1)
    u_new(2:end-1, :) = u_star(2:end-1, :) ...
        - dt * (p_prime(2:end, :) - p_prime(1:end-1, :)) / dr;

    % Axial velocity at interior z-faces (columns 2:end-1)
    w_new(:, 2:end-1) = w_star(:, 2:end-1) ...
        - dt * (p_prime(:, 2:end) - p_prime(:, 1:end-1)) / dz;
end
