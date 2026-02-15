function p_new = update_pressure(p, p_prime, zeta)
% UPDATE_PRESSURE  p^{n+1} = p^n + zeta * p'

    p_new = p + zeta * p_prime;
end
