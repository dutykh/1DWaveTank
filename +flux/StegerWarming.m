function F_num = StegerWarming(wL, wR, cfg)

    % StegerWarming Calculate numerical flux using Steger-Warming Flux Vector Splitting.
    %
    %   F_num = StegerWarming(wL, wR, cfg)
    %
    %   Inputs:
    %       wL   : State vector [hL; qL] at the left side of the interface.
    %       wR   : State vector [hR; qR] at the right side of the interface.
    %       cfg  : Configuration structure (must contain cfg.param.g).
    %
    %   Outputs:
    %       F_num: Numerical flux vector [Fh; Fq] across the interface.
    %
    %   Method:
    %       Uses the Steger-Warming flux vector splitting scheme:
    %       F_num = A+(wL)*wL + A-(wR)*wR
    %       where A+ and A- are the positive and negative parts of the Jacobian A(w),
    %       calculated via spectral decomposition A = R * Lambda * L.
    %       A+ = R * max(Lambda, 0) * L
    %       A- = R * min(Lambda, 0) * L
    %
    %   Reference:
    %       Steger, J. L., & Warming, R. F. (1981). Flux vector splitting of the 
    %       inviscid gasdynamic equations with application to finite-difference methods. 
    %       Journal of Computational Physics, 40(2), 263-293.

    persistent g
    if isempty(g)
        if isfield(cfg, 'param') && isfield(cfg.param, 'g')
            g = cfg.param.g;
        else
            g = 9.81; % Default gravity
            warning('Using default g = 9.81 m/s^2');
        end
    end

    % --- Helper function to calculate A+, A- for a given state w --- 
    function [A_plus, A_minus] = calculate_split_jacobian(w, g_val)
        % Ensure w is a column vector
        w_col = w(:);
        h = w_col(1);
        q = w_col(2);

        % Handle dry states
        if h <= 1e-6 % Threshold for dry state
            A_plus = zeros(2, 2);
            A_minus = zeros(2, 2);
            return;
        end

        u = q / h;
        c = sqrt(g_val * h);

        % Eigenvalues
        lam1 = u - c;
        lam2 = u + c;

        % Right eigenvectors (columns)
        R = [1, 1; u - c, u + c];

        % Left eigenvectors (rows) - scaled by 2c for simplicity here
        L_scaled = [u + c, -1; 
                   -u + c,  1]; 
        inv_2c = 1.0 / (2.0 * c);
        L = inv_2c * L_scaled;

        % Positive and negative parts of Lambda
        lam1p = max(lam1, 0);
        lam2p = max(lam2, 0);
        lam1m = min(lam1, 0);
        lam2m = min(lam2, 0);
        
        Lambda_p = diag([lam1p, lam2p]);
        Lambda_m = diag([lam1m, lam2m]);

        % Positive and negative parts of Jacobian A = R*Lambda*L
        A_plus = R * Lambda_p * L;
        A_minus = R * Lambda_m * L;
    end

    % --- Calculate split Jacobians for Left and Right states --- 
    [Ap_L, ~] = calculate_split_jacobian(wL, g); % Only need A+(wL)
    [~, Am_R] = calculate_split_jacobian(wR, g); % Only need A-(wR)
    
    % --- Combine to get the numerical flux --- 
    % F_num = A+(wL)*wL + A-(wR)*wR
    F_num = Ap_L * wL(:) + Am_R * wR(:); % Ensure wL, wR are column vectors

    % Ensure output is ROW vector for compatibility with rhs_nsw_1st_order
    F_num = F_num(:).'; 
 
end