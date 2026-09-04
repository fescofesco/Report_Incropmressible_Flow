function face_val = upwind2_face(phi_im1, phi_i, phi_ip1, phi_ip2, vel)
% UPWIND2_FACE  2nd-order upwind (linear-extrapolation / LUD) face value.
%
%   face value = 1.5*phi_i   - 0.5*phi_im1   if vel > 0  (upwind: low-index side)
%              = 1.5*phi_ip1 - 0.5*phi_ip2   if vel < 0  (upwind: high-index side)
%
%   phi_im1, phi_i, phi_ip1, phi_ip2 : values at relative positions i-1, i,
%       i+1, i+2 with respect to the face sitting between cells i and i+1.
%   vel : transport velocity at the face.
%
%   Retains the upwind stability bias needed at high cell Peclet number
%   while being formally 2nd-order accurate (unlike pure central
%   differencing or 1st-order upwind).

    face_val = (1.5 * phi_i - 0.5 * phi_im1) .* (vel > 0) ...
              + (1.5 * phi_ip1 - 0.5 * phi_ip2) .* (vel <= 0);
end
