%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/StegerWarming.m
%
% Purpose:
%   Computes the Steger-Warming numerical flux for the 1D Non-Linear Shallow
%   Water Equations (NSW). This is a classical Flux Vector Splitting (FVS)
%   scheme where the flux is split into positive and negative parts using the
%   spectral decomposition of the Jacobian matrix. The numerical flux at the
%   interface is then F_num = A+(wL)*wL + A-(wR)*wR.
%
% Syntax:
%   F_num = StegerWarming(wL, wR, cfg)
%
% Inputs:
%   wL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g (gravity), cfg.phys.dry_tolerance.
%
% Outputs:
%   F_num - [1 x 2, double] Steger-Warming numerical flux vector [Fh, Fq].
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g and cfg.phys.dry_tolerance).
%
% References:
%   - Steger, J. L., & Warming, R. F. (1981). Flux vector splitting of the
%     inviscid gasdynamic equations with application to finite-difference methods.
%     Journal of Computational Physics, 40(2), 263-293.
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 8)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F_num = StegerWarming(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Persistent Variable for Gravity                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    persistent g;
    if isempty(g)
        if isfield(cfg, 'phys') && isfield(cfg.phys, 'g')
            g = cfg.phys.g;     % [m/s^2] Gravity
        else
            g = 9.81; % Default gravity if not found in cfg
            warning('StegerWarming:UsingDefaultG', 'Using default g = 9.81 m/s^2');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested Helper Function: Calculate Split Jacobians           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A_plus, A_minus] = calculate_split_jacobian(w)
        % Ensure w is a column vector
        w_col = w(:);
        h = w_col(1);   % [m] Water depth
        q = w_col(2);   % [m^2/s] Discharge

        % --- Handle Dry States ---
        if h <= cfg.phys.dry_tolerance % Threshold for dry state
            A_plus = zeros(2, 2);
            A_minus = zeros(2, 2);
            return;
        end

        % --- Calculate Primitive Variables and Wave Speed ---
        u = q / h;          % [m/s] Velocity
        c = sqrt(g * h);    % [m/s] Wave speed (celerity)

        % --- Spectral Decomposition of Jacobian A(w) = R * Lambda * L ---
        % Eigenvalues (lambda_1, lambda_2)
        lam1 = u - c; % [m/s]
        lam2 = u + c; % [m/s]

        % Right eigenvectors (columns of R)
        R = [1,   1;...
             u-c, u+c];

        % Left eigenvectors (rows of L) - Note: Scaled by 1/(2c)
        L_scaled = [ u+c, -1;...
                    -u+c,  1];
        if abs(c) < 1e-10 % Avoid division by zero
            inv_2c = 0; % Or handle differently? Eigenvectors ill-defined if c=0
        else
            inv_2c = 1.0 / (2.0 * c);
        end
        L = inv_2c * L_scaled;

        % --- Split Eigenvalues ---
        % Lambda+ contains max(lambda_i, 0)
        % Lambda- contains min(lambda_i, 0)
        lam1p = max(lam1, 0);
        lam2p = max(lam2, 0);
        lam1m = min(lam1, 0);
        lam2m = min(lam2, 0);

        Lambda_p = diag([lam1p, lam2p]);
        Lambda_m = diag([lam1m, lam2m]);

        % --- Calculate Split Jacobian Matrices A+ and A- ---
        % A+ = R * Lambda+ * L
        % A- = R * Lambda- * L
        A_plus = R * Lambda_p * L;
        A_minus = R * Lambda_m * L;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of Nested Helper Function                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Steger-Warming Flux: F_num = A+(wL)*wL + A-(wR)*wR %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Ap_L, ~] = calculate_split_jacobian(wL); % Only need A+(wL)
    [~, Am_R] = calculate_split_jacobian(wR); % Only need A-(wR)

    F_num = Ap_L * wL(:) + Am_R * wR(:); % [m^2/s; m^3/s^2], ensure column vectors

    % Ensure output is ROW vector for compatibility with rhs function
    F_num = F_num(:)';

end