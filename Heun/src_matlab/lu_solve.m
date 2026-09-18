function X = lu_solve(LU, piv, B)
% LU_SOLVE  Solve A*X = B from the factors produced by LU_FACTOR.
%
%   X = lu_solve(LU, piv, B)
%
%   Forward substitution with the unit lower triangle L, then back
%   substitution with U. Written out explicitly -- no backslash.
%
%   Inputs:
%     LU, piv : output of LU_FACTOR
%     B       : (n x 1) or (n x m) one or several right-hand sides
%
%   Output:
%     X       : same size as B
%
%   See also LU_FACTOR, FACTORIZE_POISSON.

    n = size(LU, 1);
    X = B;

    % --- apply the recorded row interchanges:  b <- P*b
    for k = 1:n
        m = piv(k);
        if m ~= k
            X([k m], :) = X([m k], :);
        end
    end

    % --- forward substitution  L*Y = P*b   (L has unit diagonal)
    for k = 1:(n-1)
        X((k+1):n, :) = X((k+1):n, :) - LU((k+1):n, k) * X(k, :);
    end

    % --- back substitution  U*X = Y
    for k = n:-1:1
        X(k, :) = X(k, :) / LU(k, k);
        X(1:(k-1), :) = X(1:(k-1), :) - LU(1:(k-1), k) * X(k, :);
    end
end
