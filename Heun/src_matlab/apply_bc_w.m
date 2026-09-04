function w = apply_bc_w(w)
% APPLY_BC_W  Enforce axial velocity boundary conditions.
%
%   w(:, 1)   = 1   inlet: uniform w* = 1
%   w(:, end) = w(:, end-1)   outlet: zero-gradient

    w(:, 1)   = 1.0;          % inlet
    w(:, end) = w(:, end-1);  % outlet zero-gradient
end
