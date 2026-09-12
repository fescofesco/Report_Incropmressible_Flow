function A = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z)
% BUILD_POISSON_MATRIX  Assemble sparse Poisson matrix (called once).
%
%   A = build_poisson_matrix(r_c, r_f, dr, dz, n_r, n_z)
%
%   Returns A : sparse matrix (n_r*n_z x n_r*n_z)
%
%   Flat index: k = (i-1)*n_z + j   (1-based)
%
%   BCs for p':
%     Centreline (i=1):    Neumann (r_f(1)=0, flux vanishes naturally)
%     Wall (i=n_r):        Neumann (outer flux = 0)
%     Inlet (j=1):         Neumann (south flux = 0)
%     Outlet (j=n_z):      Dirichlet p' = 0

    N = n_r * n_z;

    % Pre-allocate triplets (at most 5 entries per row)
    rows = zeros(5*N, 1);
    cols = zeros(5*N, 1);
    vals = zeros(5*N, 1);
    idx = 0;

    for i = 1:n_r
        for j = 1:n_z
            k = (i-1)*n_z + j;

            % Dirichlet at outlet (j = n_z): p' = 0
            if j == n_z
                idx = idx + 1;
                rows(idx) = k; cols(idx) = k; vals(idx) = 1.0;
                continue;
            end

            % Coefficients from cylindrical Poisson
            rf_o = r_f(i + 1);   % r_{i+1/2}
            rf_i_val = r_f(i);   % r_{i-1/2}
            rc = r_c(i);

            coeff_E = rf_o / (rc * dr^2);       % p'(i+1, j)
            coeff_W = rf_i_val / (rc * dr^2);   % p'(i-1, j)
            coeff_N = 1.0 / dz^2;               % p'(i, j+1)
            coeff_S = 1.0 / dz^2;               % p'(i, j-1)

            % Neumann wall (i = n_r): outer flux = 0
            if i == n_r
                coeff_E = 0.0;
            end

            % Neumann inlet (j = 1): south flux = 0
            if j == 1
                coeff_S = 0.0;
            end

            % Centre coefficient
            coeff_C = -(coeff_E + coeff_W + coeff_N + coeff_S);

            % Diagonal
            idx = idx + 1;
            rows(idx) = k; cols(idx) = k; vals(idx) = coeff_C;

            % East (i+1, j)
            if i < n_r
                idx = idx + 1;
                rows(idx) = k; cols(idx) = i*n_z + j; vals(idx) = coeff_E;
            end

            % West (i-1, j)
            if i > 1
                idx = idx + 1;
                rows(idx) = k; cols(idx) = (i-2)*n_z + j; vals(idx) = coeff_W;
            end

            % North (i, j+1)
            if j < n_z
                idx = idx + 1;
                rows(idx) = k; cols(idx) = (i-1)*n_z + (j+1); vals(idx) = coeff_N;
            end

            % South (i, j-1)
            if j > 1
                idx = idx + 1;
                rows(idx) = k; cols(idx) = (i-1)*n_z + (j-1); vals(idx) = coeff_S;
            end
        end
    end

    % Trim and build sparse matrix
    rows = rows(1:idx);
    cols = cols(1:idx);
    vals = vals(1:idx);
    A = sparse(rows, cols, vals, N, N);
end
