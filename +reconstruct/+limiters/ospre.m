function limited_slope = ospre(delta_minus, delta_plus, cfg) %#ok<INUSD> cfg unused
% OSPRE Smooth blended limiter for TVD schemes
%
% OSPRE limiter provides a smooth blend between minmod and superbee,
% giving good results for many practical problems. With
% r = delta_plus/delta_minus,
%     phi(r) = 1.5 (r^2 + r)/(r^2 + r + 1),
% and the limited slope phi(r) delta_minus is written below without the
% quotient r, symmetric in the two differences:
%     slope = 1.5 delta_minus delta_plus (delta_minus + delta_plus)
%             / (delta_minus^2 + delta_minus delta_plus + delta_plus^2)
% when delta_minus delta_plus > 0, and 0 otherwise.
%
% Syntax:
%   limited_slope = ospre(delta_minus, delta_plus)
%   limited_slope = ospre(delta_minus, delta_plus, cfg)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%   cfg         - (optional, unused) configuration, as the MUSCL
%                 reconstruction passes it to every limiter
%
% Output:
%   limited_slope - Limited slope
%
% References:
%   - The limiter of Waterson and Deconinck, in the form given by F. Kemm,
%     A comparative study of TVD-limiters: well-known limiters and an
%     introduction of new ones, Int. J. Numer. Meth. Fluids 67(4), 404-440
%     (2011).
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign: phi(r) delta_minus with r = delta_plus/delta_minus, in a form
    % that needs no division by delta_minus
    limited_slope = 1.5 * delta_minus * delta_plus * (delta_minus + delta_plus) ...
                    / (delta_minus^2 + delta_minus * delta_plus + delta_plus^2);
end
end
