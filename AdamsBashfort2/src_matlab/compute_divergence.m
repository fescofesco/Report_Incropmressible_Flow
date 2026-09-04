function div = compute_divergence(u, w, r_c, r_f, dr, dz)
% COMPUTE_DIVERGENCE  Pointwise discrete divergence of velocity field.
%
%   div(i,j) = [r_{i+1/2}*u(i+1,j) - r_{i-1/2}*u(i,j)] / (r_c(i)*dr)
%            + [w(i,j+1) - w(i,j)] / dz
%
%   Returns div : (n_r, n_z)

    rf_o = r_f(2:end);      % (n_r, 1)
    rf_i = r_f(1:end-1);    % (n_r, 1)
    rc   = r_c;             % (n_r, 1)

    div_u = (rf_o .* u(2:end, :) - rf_i .* u(1:end-1, :)) ./ (rc * dr);
    div_w = (w(:, 2:end) - w(:, 1:end-1)) / dz;

    div = div_u + div_w;
end
