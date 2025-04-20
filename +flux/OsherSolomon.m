function F_num = OsherSolomon(wL, wR, cfg)

    % OsherSolomon Calculate numerical flux using Osher's scheme (Flux Vector Splitting variant).
    %
    %   F_num = OsherSolomon(wL, wR, cfg)
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
    %       Uses the flux vector splitting variant of Osher's scheme:
    %       F_num = F+(wL) + F-(wR)
    %       where F+(w) = A+(w)*w and F-(w) = A-(w)*w.
    %       A+ and A- are the positive and negative parts of the Jacobian A(w),
    %       calculated via spectral decomposition A = R * Lambda * L.
    %
    %   Reference:
    %       Osher, S., & Solomon, F. (1982). Upwind difference schemes for 
    %       hyperbolic systems of conservation laws. 
    %       Mathematics of Computation, 38(158), 339-374.

    persistent g

    if isempty(g)
        if isfield(cfg, 'param') && isfield(cfg.param, 'g')
            g = cfg.param.g;
        else
            g = 9.81; % Default gravity
            warning('Using default g = 9.81 m/s^2');
        end
    end

    % --- Helper function to calculate F+, F- for a given state w --- 
    function [Fp, Fm] = calculate_split_flux(w, g_val)
        % Ensure w is a column vector for calculations
        w_col = w(:); 
        h = w_col(1);
        q = w_col(2);

        % Handle dry states
        if h <= 1e-6 % Threshold for dry state
            Fp = [0; 0];
            Fm = [0; 0];
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

        % Split fluxes (using column vector w_col)
        Fp = A_plus * w_col;
        Fm = A_minus * w_col;
    end

    % --- Calculate split fluxes for Left and Right states --- 
    [Fp_L, ~] = calculate_split_flux(wL, g); % Only need F+(wL)
    [~, Fm_R] = calculate_split_flux(wR, g); % Only need F-(wR)
    
    % --- Combine to get the numerical flux --- 
    F_num = Fp_L + Fm_R;

    % Ensure output is ROW vector for compatibility with rhs_nsw_1st_order
    F_num = F_num(:).'; 
    
end