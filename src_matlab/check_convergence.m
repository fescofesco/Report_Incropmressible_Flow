function converged = check_convergence(res, tol_cont, tol_vel, tol_temp)
% CHECK_CONVERGENCE  Return true when all residuals are below tolerances.

    converged = (res.continuity < tol_cont) && ...
                (res.momentum < tol_vel) && ...
                (res.temperature < tol_temp);
end
