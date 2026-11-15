function limited_slope = koren(delta_minus, delta_plus, cfg) %#ok<INUSD> cfg unused
% KOREN Third-order accurate TVD limiter
%
% Koren's limiter provides third-order accuracy in smooth regions
% while maintaining TVD properties. With r = delta_plus/delta_minus,
%     phi(r) = max(0, min(2 r, (1 + 2 r)/3, 2)),
% and the limited slope is phi(r) delta_minus, written below without the
% quotient r as
%     slope = s min(2 |delta_plus|, (|delta_minus| + 2 |delta_plus|)/3, 2 |delta_minus|)
% when delta_minus delta_plus > 0 (s the common sign), and 0 otherwise.
% Where (1 + 2 r)/3 is selected, the value q_i + slope/2 at the right edge
% of the cell is that of the third-order upwind-biased reconstruction
% (kappa = 1/3).
%
% Syntax:
%   limited_slope = koren(delta_minus, delta_plus)
%   limited_slope = koren(delta_minus, delta_plus, cfg)
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
%   - B. Koren (1993), the limiter in the form given by G. W. Griffiths and
%     W. E. Schiesser, Traveling Wave Analysis of Partial Differential
%     Equations, Academic Press 2011, p. 40.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign: phi(r) delta_minus with r = delta_plus/delta_minus, in a form
    % that needs no division
    a = abs(delta_minus);
    b = abs(delta_plus);
    limited_slope = sign(delta_minus) * min([2*b, (a + 2*b)/3, 2*a]);
end
end
