function [LU, piv] = lu_factor(A)
% LU_FACTOR  Doolittle LU factorisation with partial pivoting:  P*A = L*U.
%
%   [LU, piv] = lu_factor(A)
%
%   The assignment requires a DIRECT solution method for the pressure Poisson
%   equation that is programmed, not a library call. This routine and
%   LU_SOLVE are that elimination kernel, written out explicitly -- no
%   backslash, no lu(), no inv(), no decomposition(). They act on the small
%   (n_r x n_r) blocks of the block-tridiagonal Poisson operator and are
%   called only during set-up, from FACTORIZE_POISSON.
%
%   Inputs:
%     A   : (n x n) dense matrix (not modified)
%
%   Outputs:
%     LU  : (n x n) packed factors. The strict lower triangle holds L (its
%           unit diagonal is implied); the upper triangle including the
%           diagonal holds U.
%     piv : (n x 1) piv(k) is the row that was interchanged with row k at
%           elimination step k.
%
%   See also LU_SOLVE, FACTORIZE_POISSON, SOLVE_POISSON.

    n = size(A, 1);
    LU = A;
    piv = zeros(n, 1);

    for k = 1:n
        % --- partial pivoting: largest magnitude in column k, rows k..n
        [pivot_val, m] = max(abs(LU(k:n, k)));
        m = m + k - 1;
        piv(k) = m;
        if pivot_val == 0
            error('lu_factor:singular', 'Zero pivot in column %d.', k);
        end
        if m ~= k
            LU([k m], :) = LU([m k], :);
        end

        % --- elimination; the multipliers are stored in the created zeros
        rows = (k+1):n;
        LU(rows, k) = LU(rows, k) / LU(k, k);
        LU(rows, (k+1):n) = LU(rows, (k+1):n) - LU(rows, k) * LU(k, (k+1):n);
    end
end
