function u = apply_bc_u(u)
% APPLY_BC_U  Enforce radial velocity boundary conditions.
%
%   u(1, :)   = 0   centreline (r = 0, symmetry)
%   u(end, :) = 0   wall (r = R, no penetration)
%
%   NOTE: u(:,1) is at z_c(1), NOT the inlet face.
%   The inlet BC u=0 at z=0 is handled by the ghost cell in the predictor.

    u(1, :)   = 0;   % axis
    u(end, :) = 0;   % wall
end
