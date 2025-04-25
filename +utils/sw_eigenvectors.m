function [L, R, lambdas] = sw_eigenvectors(H, U, C, cfg)
% SW_EIGENVECTORS Calculates eigenvalues and eigenvectors for the shallow water equations.
%
% The Jacobian matrix A = [0 1; c^2-u^2 2u], where c = sqrt(gH).
%
% Syntax:
%   [L, R, lambdas] = sw_eigenvectors(H, U, C, cfg)
%
% Inputs:
%   H   - [scalar] Water depth [m]
%   U   - [scalar] Velocity [m/s]
%   C   - [scalar] Wave speed sqrt(gH) [m/s]
%   cfg - [struct] Configuration structure (not used currently, but for consistency)
%
% Outputs:
%   L       - [2x2] Matrix of left eigenvectors (rows)
%   R       - [2x2] Matrix of right eigenvectors (columns)
%   lambdas - [1x2] Eigenvalues [u-c, u+c]
%
% Reference:
%   LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%   Cambridge University Press. (Chapter 13)

    % Eigenvalues
    lambda1 = U - C;
    lambda2 = U + C;
    lambdas = [lambda1, lambda2];

    % Use global epsilon from config for division-by-zero protection
    epsilon = cfg.numerics.epsilon;

    % Right Eigenvectors (columns)
    R = [ 1,    1; ...
          U-C,  U+C ];

    % Left Eigenvectors (rows)
    % L = inv(R)
    if abs(C) < epsilon % Handle case C=0 (dry or critical flow)
        % Degenerate case, might need special handling depending on context.
        % For now, return identity matrices or raise an error.
        % Let's return a pseudo-inverse based on limits, but care is needed.
        % If C -> 0, R -> [1 1; U U]. This is singular.
        % A -> [0 1; -U^2 2U]. Eigenvalues are U, U. Eigenvector is [1; U].
        % For reconstruction purposes, maybe use identity? Or signal an issue.
        % Let's set L and R to identity for now if C is near zero. 
        % This will effectively revert to component-wise for that interface.
        L = eye(2);
        R = eye(2);
        % warning('sw_eigenvectors:ZeroWaveSpeed', 'Wave speed C is near zero (H=%g). Eigenvectors are ill-defined.', H);
    else
        inv_2C = 1.0 / (2.0 * C);
        L = inv_2C * [ U+C,  -1; ...
                      -U+C,   1 ];
    end

end
