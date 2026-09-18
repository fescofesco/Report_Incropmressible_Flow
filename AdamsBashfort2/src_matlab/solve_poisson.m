function p_prime = solve_poisson(P, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
% SOLVE_POISSON  Solve  lap(p') = (1/dt)*div(u*)  for the pressure correction,
%   by the block-Thomas sweeps prepared in FACTORIZE_POISSON.
%
%   p_prime = solve_poisson(P, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
%
%   This is a DIRECT method, programmed here rather than delegated: no
%   backslash, no decomposition(), no inv(). With the unknowns grouped by
%   axial station, column j of the (n_r x n_z) field IS block j of the
%   block-tridiagonal system, so no flattening or flat-index bookkeeping is
%   needed anywhere.
%
%       forward :  y_j = inv(M_j)*(b_j - a_z*y_{j-1})
%       back    :  p_j = y_j - a_z*inv(M_j)*p_{j+1}
%
%   Inputs:
%     P       : struct from FACTORIZE_POISSON (built once, before the time loop)
%     u_star  : (n_r+1 x n_z) predicted radial velocity
%     w_star  : (n_r x n_z+1) predicted axial velocity
%
%   Output:
%     p_prime : (n_r x n_z)
%
%   See also FACTORIZE_POISSON, LU_FACTOR, LU_SOLVE.

    % ---- right-hand side:  b(i,j) = (1/dt)*div(u*) -------------------------
    rf_o = r_f(2:end);      % (n_r, 1)
    rf_i = r_f(1:end-1);    % (n_r, 1)
    rc   = r_c;             % (n_r, 1)

    % Radial divergence
    div_u = (rf_o .* u_star(2:end, :) - rf_i .* u_star(1:end-1, :)) ./ (rc * dr);

    % Axial divergence
    div_w = (w_star(:, 2:end) - w_star(:, 1:end-1)) / dz;

    b = (div_u + div_w) / dt;   % (n_r, n_z)

    % Outlet column (j = n_z) has Dirichlet p' = 0 -> RHS = 0
    b(:, end) = 0.0;

    % ---- forward sweep ------------------------------------------------------
    Minv = P.Minv;
    az   = P.az;

    y = zeros(n_r, n_z);
    y(:, 1) = Minv(:, :, 1) * b(:, 1);
    for j = 2:(n_z - 1)
        y(:, j) = Minv(:, :, j) * (b(:, j) - az * y(:, j-1));
    end
    y(:, n_z) = b(:, n_z);      % Dirichlet row: M = I, no south coupling

    % ---- back substitution --------------------------------------------------
    p_prime = zeros(n_r, n_z);
    p_prime(:, n_z) = y(:, n_z);
    for j = (n_z - 1):-1:1
        p_prime(:, j) = y(:, j) - az * (Minv(:, :, j) * p_prime(:, j+1));
    end
end
