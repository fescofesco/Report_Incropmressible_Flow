function p_prime = solve_poisson(A, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
% SOLVE_POISSON  Solve A * p' = b for pressure correction.
%
%   p_prime = solve_poisson(A, u_star, w_star, r_c, r_f, dr, dz, dt, n_r, n_z)
%
%   RHS: b(i,j) = (1/dt) * div(u*)
%   Returns p_prime : (n_r, n_z)

    rf_o = r_f(2:end);      % (n_r, 1)
    rf_i = r_f(1:end-1);    % (n_r, 1)
    rc   = r_c;             % (n_r, 1)

    % Radial divergence
    div_u = (rf_o .* u_star(2:end, :) - rf_i .* u_star(1:end-1, :)) ./ (rc * dr);

    % Axial divergence
    div_w = (w_star(:, 2:end) - w_star(:, 1:end-1)) / dz;

    b = (div_u + div_w) / dt;   % (n_r, n_z)

    % Outlet row (j=n_z) has Dirichlet p'=0 -> RHS = 0
    b(:, end) = 0.0;

    % Flatten and solve
    b_flat = b(:);              % column-major flatten (same as reshape(b, [], 1))
    % MATLAB sparse direct solve (backslash)
    % Note: MATLAB stores column-major, Python row-major. The flat index
    % must match the matrix assembly order: k = (i-1)*n_z + j.
    % Reshaping b column-major: b(:) goes column by column (j varies fastest
    % for fixed i when b is (n_r, n_z)). But our matrix uses k=(i-1)*n_z+j,
    % where j varies fastest for fixed i. This matches b(:) since MATLAB
    % flattens column-major on (n_r, n_z) -> i varies fastest.
    % We need to transpose to match: b' is (n_z, n_r), then b'(:) has
    % j varying slowest... no. Let me be precise.
    %
    % Matrix assembly: k = (i-1)*n_z + j, so for i=1, j=1..n_z -> k=1..n_z.
    % MATLAB b(:) for b of size (n_r, n_z) gives b(1,1), b(2,1), ..., b(n_r,1),
    %   b(1,2), ..., i.e. i varies fastest. But we need j to vary fastest.
    % Solution: transpose b before flattening.
    b_flat = reshape(b', [], 1);   % b' is (n_z, n_r), flatten -> j varies fastest

    p_prime_flat = A \ b_flat;

    % Reshape back: p_prime_flat has k=(i-1)*n_z+j ordering
    p_prime = reshape(p_prime_flat, n_z, n_r)';   % (n_r, n_z)
end
