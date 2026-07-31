function p_new = update_pressure(p, p_prime, zeta)
% UPDATE_PRESSURE  p^{n+1} = p^n + zeta*(p' - p^n) = (1-zeta)*p^n + zeta*p'
%
%   NOTE: this is a relaxed blend towards p', not an accumulation of p'
%   onto p^n. The momentum predictor excludes the pressure gradient
%   entirely (non-incremental / Chorin-type projection), so p' as computed
%   by solve_poisson already IS the full physical pressure for this step,
%   not a vanishing correction -- accumulating it additively made p grow
%   roughly linearly in time even at steady state. See the corresponding
%   note in solver_poisson.py::update_pressure (Python) for the derivation.

    p_new = p + zeta * (p_prime - p);
end
