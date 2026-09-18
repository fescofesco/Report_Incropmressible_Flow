function [u_new, w_new, p_new, T_new, rhs_u, rhs_w, rhs_T] = ...
    ab2_step(u, w, p, T, rhs_u_old, rhs_w_old, rhs_T_old, ...
             r_c, r_f, dr, dz, dt, Re, Pr, n_r, n_z, poisson_fac, alpha_p)
% AB2_STEP  Adams-Bashforth 2 with one incremental pressure correction.
% The first call uses forward Euler when the previous RHS arrays are empty.

    rhs_u = compute_rhs_u(u, w, r_c, r_f, dr, dz, Re);
    rhs_w = compute_rhs_w(w, u, r_c, r_f, dr, dz, Re);
    rhs_T = compute_rhs_T(T, u, w, r_c, r_f, dr, dz, Re, Pr);

    if isempty(rhs_u_old)
        explicit_u = rhs_u;
        explicit_w = rhs_w;
        explicit_T = rhs_T;
    else
        explicit_u = 1.5 * rhs_u - 0.5 * rhs_u_old;
        explicit_w = 1.5 * rhs_w - 0.5 * rhs_w_old;
        explicit_T = 1.5 * rhs_T - 0.5 * rhs_T_old;
    end

    % The old pressure gradient makes the Poisson unknown a correction p'.
    u_star = u;
    grad_p_r = (p(2:end, :) - p(1:end-1, :)) / dr;
    u_star(2:end-1, :) = u(2:end-1, :) ...
        + dt * (explicit_u - grad_p_r);

    w_star = w;
    grad_p_z = (p(:, 2:end) - p(:, 1:end-1)) / dz;
    w_star(:, 2:end-1) = w(:, 2:end-1) ...
        + dt * (explicit_w - grad_p_z);
    u_star = apply_bc_u(u_star);
    w_star = apply_bc_w(w_star);

    p_prime = solve_poisson(poisson_fac, u_star, w_star, ...
                            r_c, r_f, dr, dz, dt, n_r, n_z);
    [u_new, w_new] = correct_velocity(u_star, w_star, ...
                                      p_prime, dr, dz, dt);
    u_new = apply_bc_u(u_new);
    w_new = apply_bc_w(w_new);

    p_new = update_pressure(p, p_prime, alpha_p);
    T_new = T + dt * explicit_T;
end
