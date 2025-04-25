%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/colebrook_white.m
%
% Purpose:
%   Implements the Colebrook-White formula to calculate the Darcy friction 
%   factor f for the Darcy-Weisbach equation.
%   Formula: 1/√f = -2log₁₀(ks/(14.8R) + 2.51/(Re√f))
%
% Syntax:
%   f = colebrook_white(H, U, cfg)
%
% Inputs:
%   H    - [N x 1, double] Water depth at each cell [m].
%   U    - [N x 1, double] Velocity at each cell [m/s].
%   cfg  - [struct] Configuration structure containing:
%          cfg.phys.ks: [double] Equivalent sand roughness height [m].
%          cfg.phys.kinematic_viscosity: [double] Kinematic viscosity of water [m²/s].
%          cfg.phys.cw_iterations: [integer] Max iterations for solving implicit equation.
%          cfg.phys.cw_tolerance: [double] Convergence tolerance.
%
% Outputs:
%   f   - [N x 1, double] Darcy friction factor for each cell.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = colebrook_white(H, U, cfg)
    % COLEBROOK_WHITE Calculate Darcy friction factor 'f' using the Colebrook-White equation.
    %
    % Syntax:
    %   f = colebrook_white(H, U, cfg)
    %
    % Inputs:
    %   H   - Water depth [m] (scalar or vector)
    %   U   - Flow velocity [m/s] (scalar or vector, same size as H)
    %   cfg - Configuration structure containing:
    %         cfg.phys.ks (optional) - Equivalent sand roughness height [m].
    %                                  Default: 0.001 m
    %         cfg.phys.kinematic_viscosity (optional) - Kinematic viscosity of fluid [m^2/s].
    %                                                  Default: 1e-6 m^2/s (water at 20°C)
    %         cfg.solver.cw_iterations (optional) - Max iterations for Colebrook-White.
    %                                                Default: 100
    %         cfg.solver.cw_tolerance (optional) - Convergence tolerance for f.
    %                                              Default: 1e-6
    %
    % Output:
    %   f   - Darcy friction factor [-] (scalar or vector, same size as H and U)
    %
    % Description:
    %   This function calculates the Darcy friction factor 'f' based on the
    %   Colebrook-White equation, which is implicit and solved iteratively.
    %   It handles both laminar and turbulent flow regimes.
    %   For laminar flow (Re < 2300), f = 64/Re.
    %   For turbulent flow, the Colebrook-White equation is solved:
    %   1/sqrt(f) = -2 * log10( (ks/(3.7*Dh)) + (2.51/(Re*sqrt(f))) )
    %   where Dh is the hydraulic diameter (approximated as 4*H for wide channels).
    %
    %   The implementation uses an iterative fixed-point method to find f.
    %   The function is vectorized to handle array inputs for H and U efficiently.
    %
    % Author: Denys Dutykh
    % Date: 2024-04-23

    % Input validation
    if ~isequal(size(H), size(U))
        error('Inputs H and U must have the same size.');
    end

    % Get parameters from config or use defaults
    if isfield(cfg.phys, 'ks') && ~isempty(cfg.phys.ks)
        ks = cfg.phys.ks;
    else
        warning('Equivalent sand roughness not specified. Using default ks = 0.001 m');
        ks = 0.001; % Default value (1 mm roughness)
    end
    
    if isfield(cfg.phys, 'kinematic_viscosity') && ~isempty(cfg.phys.kinematic_viscosity)
        nu = cfg.phys.kinematic_viscosity;
    else
        nu = 1e-6; % Default value for water at 20°C [m²/s]
    end
    
    if isfield(cfg.phys, 'cw_iterations')
        max_iter = cfg.phys.cw_iterations;
    else
        max_iter = 100; % Default iteration limit
    end
    
    if isfield(cfg.phys, 'cw_tolerance')
        tol = cfg.phys.cw_tolerance;
    else
        tol = cfg.numerics.epsilon; % Use config epsilon as default convergence tolerance
    end

    % Initialize output array
    f = zeros(size(H));
    valid_indices = H > cfg.phys.dry_tolerance & abs(U) > cfg.numerics.epsilon; % Avoid division by zero or near-zero H/U

    if ~any(valid_indices)
        return; % No valid points to calculate friction for
    end

    % --- Calculations for valid indices --- 
    H_valid = H(valid_indices);
    U_valid = U(valid_indices);

    % Hydraulic Diameter (assuming wide channel, Dh approx 4*H)
    Dh = 4 * H_valid;

    % Reynolds number
    Re = abs(U_valid) .* Dh / nu;

    % Relative roughness
    rel_rough = ks ./ Dh;

    % Initialize friction factor array for valid points
    f_valid = zeros(size(H_valid));

    % --- Laminar Flow (Re < 2300) --- 
    laminar_mask = Re < 2300;
    if any(laminar_mask)
        f_valid(laminar_mask) = 64 ./ Re(laminar_mask);
    end

    % --- Turbulent Flow (Re >= 2300) --- 
    turbulent_mask = ~laminar_mask;
    if any(turbulent_mask)
        Re_turb = Re(turbulent_mask);
        rel_rough_turb = rel_rough(turbulent_mask);

        % Initial guess for f (e.g., using Haaland approximation or simpler guess)
        % Using a simple explicit approximation based on log laws for fully rough flow
        % f_guess = (1 ./ (-1.8 * log10(6.9 ./ Re_turb + (rel_rough_turb / 3.7).^1.11))).^2;
        % Or simpler guess based on typical values
        f_guess = 0.02 * ones(size(Re_turb)); % Initial guess

        f_iter = f_guess;
        converged = false(size(f_iter));
        
        for iter = 1:max_iter
            if all(converged)
                break;
            end
            
            f_old = f_iter;
            
            % Colebrook-White iteration (vectorized)
            % Avoid log10(0) or sqrt(negative) - although f should remain positive
            f_iter_safe = max(f_old, cfg.numerics.epsilon); 
            rhs = -2 * log10(rel_rough_turb / 3.7 + 2.51 ./ (Re_turb .* sqrt(f_iter_safe)));
            f_new = (1 ./ rhs).^2;
            
            % Update only non-converged elements
            f_iter(~converged) = f_new(~converged);
            
            % Check convergence for non-converged elements
            delta_f = abs(f_iter(~converged) - f_old(~converged));
            newly_converged = delta_f < tol;
            converged(~converged) = newly_converged;
        end

        if iter == max_iter && ~all(converged)
            warning('Colebrook-White iteration did not converge for %d elements after %d iterations.', sum(~converged), max_iter);
        end
        
        f_valid(turbulent_mask) = f_iter;
    end

    % Assign calculated f values back to the original full array
    f(valid_indices) = f_valid;

end