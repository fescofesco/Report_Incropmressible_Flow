function res = compute_residuals(u_new, u_old, w_new, w_old, T_new, T_old, r_c, r_f, dr, dz, dt)
% COMPUTE_RESIDUALS  Compute convergence residuals.
%
%   Returns res struct with fields:
%     .continuity  : max |div(u)|
%     .momentum    : max change rate of both u and w (interior)
%     .temperature : max |T_new - T_old| / dt

    div = compute_divergence(u_new, w_new, r_c, r_f, dr, dz);
    res.continuity  = max(abs(div(:)));
    R_u = max(abs(u_new(2:end-1, :) - u_old(2:end-1, :)), [], 'all');
    R_w = max(abs(w_new(:, 2:end-1) - w_old(:, 2:end-1)), [], 'all');
    res.momentum    = max(R_u, R_w) / dt;
    res.temperature = max(abs(T_new(:) - T_old(:))) / dt;
end
