function P = factorize_poisson(r_c, r_f, dr, dz, n_r, n_z)
% FACTORIZE_POISSON  Block-Thomas factorisation of the pressure Poisson
%   operator, using the hand-written LU kernel in LU_FACTOR / LU_SOLVE.
%   Called ONCE, before the time loop.
%
%   P = factorize_poisson(r_c, r_f, dr, dz, n_r, n_z)
%
%   STRUCTURE
%   ---------
%   Group the unknowns p'(:, j) by axial station j. The 5-point stencil
%   couples a cell only to its radial neighbours (same j) and to j +/- 1, so
%   the operator is BLOCK TRIDIAGONAL with n_z blocks of size n_r x n_r:
%
%       S_j*p_{j-1} + D_j*p_j + N_j*p_{j+1} = b_j ,      j = 1 .. n_z
%
%       D_j = tridiag( w_i , -(e_i + w_i + 2*a_z) , e_i )     (interior j)
%       S_j = N_j = a_z*I ,     a_z = 1/dz^2
%       e_i = r_{i+1/2}/(r_i*dr^2) ,     w_i = r_{i-1/2}/(r_i*dr^2)
%
%   The boundary conditions enter as:
%     Wall       (i = n_r):  e_{n_r} = 0                          (Neumann)
%     Centreline (i = 1):    w_1 = 0 automatically, r_{1/2} = 0   (Neumann)
%     Inlet      (j = 1):    no south coupling, so that diagonal carries
%                            a_z once instead of twice            (Neumann)
%     Outlet     (j = n_z):  D = I, S = N = 0                     (Dirichlet)
%
%   BLOCK THOMAS FORWARD ELIMINATION
%   --------------------------------
%       M_1 = D_1
%       M_j = D_j - S_j*inv(M_{j-1})*N_{j-1} = D_j - a_z^2*inv(M_{j-1})
%
%   which reduces the solve to the two sweeps in SOLVE_POISSON,
%   y_j = inv(M_j)*(b_j - a_z*y_{j-1})  and  p_j = y_j - a_z*inv(M_j)*p_{j+1}.
%
%   The operator is constant in time, so the n_z Schur complements are
%   inverted here, once, via our own LU. Every subsequent Poisson solve then
%   costs only two (n_r x n_r) matrix-vector products per axial station, i.e.
%   O(n_r^2*n_z) -- cheaper than re-factorising anything.
%
%   STABILITY
%   ---------
%   -A is an irreducibly diagonally dominant M-matrix: |diagonal| equals the
%   sum of the off-diagonal magnitudes in every interior row and is strictly
%   greater in the Dirichlet outlet rows. Every Schur complement M_j inherits
%   that property, so no pivoting is needed between blocks. Partial pivoting
%   is used inside LU_FACTOR regardless.
%
%   Output:
%     P.Minv : (n_r x n_r x n_z) the inverted Schur complements inv(M_j)
%     P.az   : 1/dz^2
%     P.n_r, P.n_z
%
%   See also SOLVE_POISSON, LU_FACTOR, LU_SOLVE, VERIFY_POISSON_SOLVER.

    az = 1.0 / dz^2;

    % --- radial stencil coefficients (independent of j) ----------------------
    e = r_f(2:end)   ./ (r_c * dr^2);   % coupling to i+1
    w = r_f(1:end-1) ./ (r_c * dr^2);   % coupling to i-1; w(1)=0 as r_f(1)=0
    e(n_r) = 0.0;                       % Neumann at the wall

    off_super = diag(e(1:n_r-1),  1);
    off_sub   = diag(w(2:n_r),   -1);

    D_inlet = diag(-(e + w +       az)) + off_super + off_sub;  % j = 1
    D_int   = diag(-(e + w + 2.0 * az)) + off_super + off_sub;  % 2 <= j <= n_z-1

    I_n  = eye(n_r);
    Minv = zeros(n_r, n_r, n_z);

    % --- j = 1: no south coupling, so M_1 = D_1 -----------------------------
    [LU, piv] = lu_factor(D_inlet);
    Minv(:, :, 1) = lu_solve(LU, piv, I_n);

    % --- j = 2 .. n_z-1: Schur complement update, then invert ---------------
    for j = 2:(n_z - 1)
        [LU, piv] = lu_factor(D_int - az^2 * Minv(:, :, j-1));
        Minv(:, :, j) = lu_solve(LU, piv, I_n);
    end

    % --- j = n_z: Dirichlet row, D = I and no south coupling ----------------
    Minv(:, :, n_z) = I_n;

    P = struct('Minv', Minv, 'az', az, 'n_r', n_r, 'n_z', n_z);
end
