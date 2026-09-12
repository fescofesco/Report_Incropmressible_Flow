function p_new = update_pressure(p, p_prime, zeta)
% UPDATE_PRESSURE  Accumulate relaxed correction: p_new = p + zeta*p_prime.

    p_new = p + zeta * p_prime;
end
