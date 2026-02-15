function res = compute_residuals(w_new, w_old, T_new, T_old, u, w, r_c, r_f, dr, dz, dt)
% COMPUTE_RESIDUALS  Compute convergence residuals.
%
%   Returns res struct with fields:
%     .continuity  : max |div(u)|
%     .momentum    : max |w_new - w_old| / dt  (interior)
%     .temperature : max |T_new - T_old| / dt

    div = compute_divergence(u, w, r_c, r_f, dr, dz);
    res.continuity  = max(abs(div(:)));
    res.momentum    = max(abs(w_new(:, 2:end-1) - w_old(:, 2:end-1)), [], 'all') / dt;
    res.temperature = max(abs(T_new(:) - T_old(:))) / dt;
end
