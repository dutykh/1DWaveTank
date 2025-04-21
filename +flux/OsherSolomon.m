%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/OsherSolomon.m
%
% Purpose:
%   Computes the genuine Osher–Solomon numerical flux for the 1-D Non-linear
%   Shallow Water Equations (SWE). This flux is based on the Flux Vector
%   Splitting concept, where the numerical flux across an interface is
%   calculated as F_num = F+(wL) + F-(wR). For the SWE, analytical formulas
%   exist for the split fluxes F+ and F- derived from path integration in
%   the state space, ensuring desirable properties like entropy consistency.
%
% Theory:
%   F_num = F^+(w_L) + F^-(w_R)
%
%   where, for a state w = [h, q] (water depth, discharge), the analytical
%   split fluxes are:
%
%   c  = sqrt(g*h)           (celerity)
%   u  = q / h               (velocity)
%
%   F^+(w) = 1/2 * [ q + h*c ;
%                     q*u + 0.5*g*h^2 + c*q ]
%
%   F^-(w) = 1/2 * [ q - h*c ;
%                     q*u + 0.5*g*h^2 - c*q ]
%
%   These expressions represent the exact path–integrated Osher–Solomon split
%   fluxes for the SWE.
%
% Syntax:
%   F_num = OsherSolomon(wL, wR, cfg)
%
% Input Arguments:
%   wL  - [Nx2 double] State vectors [h, q] on the left side of N interfaces.
%           h: Water depth [m]
%           q: Discharge (hu) [m^2/s]
%   wR  - [Nx2 double] State vectors [h, q] on the right side of N interfaces.
%   cfg - [struct] Configuration structure containing physical parameters.
%           Required field: cfg.phys.g (gravitational acceleration [m/s^2]).
%
% Output Arguments:
%   F_num - [Nx2 double] Osher-Solomon numerical flux vector [Fh, Fq]
%           across the N interfaces.
%           Fh: Flux component for the continuity equation (h) [m^2/s]
%           Fq: Flux component for the momentum equation (q) [m^3/s^2]
%
% Dependencies:
%   None (uses basic MATLAB operations and cfg structure).
%
% References:
%   - Osher, S., & Solomon, F. (1982). Upwind difference schemes for
%     hyperbolic systems of conservation laws.
%     Mathematics of Computation, 38(158), 339-374.
%
% Author:
%   Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
%
% Date:
%   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F_num = OsherSolomon(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Persistent Gravity & Parameter Extraction                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store gravity locally using a persistent variable for efficiency,
    % avoiding repeated structure lookups within loops or vectorized calls.
    persistent g
    if isempty(g)
        % If g is not already stored, retrieve it from the cfg structure.
        if isfield(cfg, 'phys') && isfield(cfg.phys, 'g')
            g = cfg.phys.g; % [m/s^2]
        else
            % Fallback to default SI value if not found in cfg.
            g = 9.81; % [m/s^2]
            warning('OsherSolomon:DefaultG', 'Using default g = 9.81 m/s^2');
        end
    end

    N = size(wL, 1);  % Determine the number of interfaces (fluxes to compute).
    eps_flux = 1e-10; % Small tolerance to handle dry states and avoid division by zero.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested Helper Function: Calculate Analytical Split Fluxes   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This nested function calculates F+ and F- for a given set of states 'w'.
    % It's vectorized to handle multiple states (interfaces) efficiently.
    function [Fp, Fm] = split_flux(w) % Input 'w' is an Nx2 matrix [h, q]
        
        % --- Extract State Variables --- 
        h = w(:,1); % [m] Water depth (Nx1 column vector)
        q = w(:,2); % [m^2/s] Discharge (Nx1 column vector)
        
        % --- Initialize Output Fluxes --- 
        % Pre-allocate with zeros. Dry states will retain these zero flux values.
        Fp = zeros(N, 2); % F+ flux [Fp_h, Fp_q] (Nx2 matrix)
        Fm = zeros(N, 2); % F- flux [Fm_h, Fm_q] (Nx2 matrix)

        % --- Identify Wet States --- 
        % Create a logical index for states where depth h is above the tolerance.
        % Calculations involving velocity (q/h) or celerity (sqrt(g*h)) are only
        % performed for these 'wet' states to avoid numerical issues (0/0, sqrt(0)).
        is_wet = h > eps_flux; % (Nx1 logical vector)

        % --- Calculate Fluxes for Wet States Only --- 
        if any(is_wet)
            % Extract h and q only for the wet states.
            h_wet = h(is_wet); % [m]
            q_wet = q(is_wet); % [m^2/s]

            % Calculate primitive variables and celerity for wet states.
            u_wet = q_wet ./ h_wet;       % [m/s] Velocity
            c_wet = sqrt(g * h_wet);    % [m/s] Wave speed (celerity)

            % Calculate the common physical flux term F2 = q*u + 0.5*g*h^2
            % This term appears in both F+ and F- expressions.
            common_wet = q_wet .* u_wet + 0.5 * g * h_wet.^2; % [m^3/s^2]

            % Calculate F+ components using the analytical formula
            % Fp_h = 0.5 * (q + h*c)
            % Fp_q = 0.5 * (common + c*q)
            Fp_wet1 = 0.5 * (q_wet + h_wet .* c_wet); % [m^2/s]
            Fp_wet2 = 0.5 * (common_wet + c_wet .* q_wet); % [m^3/s^2]

            % Calculate F- components using the analytical formula
            % Fm_h = 0.5 * (q - h*c)
            % Fm_q = 0.5 * (common - c*q)
            Fm_wet1 = 0.5 * (q_wet - h_wet .* c_wet); % [m^2/s]
            Fm_wet2 = 0.5 * (common_wet - c_wet .* q_wet); % [m^3/s^2]

            % --- Assign Results to Output Matrices --- 
            % Place the calculated wet-state fluxes back into the pre-allocated
            % Fp and Fm matrices using the logical index 'is_wet'.
            Fp(is_wet, 1) = Fp_wet1;
            Fp(is_wet, 2) = Fp_wet2;
            Fm(is_wet, 1) = Fm_wet1;
            Fm(is_wet, 2) = Fm_wet2;
        end
        % Note: Dry states (where is_wet is false) retain their initial zero flux.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of Nested Helper Function                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate Split Fluxes and Compute Final Numerical Flux      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % --- Evaluate Split Fluxes --- 
    % Call the helper function to get the split fluxes based on left (wL)
    % and right (wR) states at each interface.
    [Fp_L, ~] = split_flux(wL);    % Calculate F+(wL), ignore F-(wL) (result is Nx2)
    [~, Fm_R] = split_flux(wR);    % Calculate F-(wR), ignore F+(wR) (result is Nx2)

    % --- Compute Osher-Solomon Numerical Flux --- 
    % The final flux is the sum of the positive-going flux from the left state
    % and the negative-going flux from the right state.
    F_num = Fp_L + Fm_R; % Element-wise addition (Nx2 matrix)
    
    % --- Final Robustness Check --- 
    % Check for any NaN values that might have occurred despite safeguards
    % (e.g., potentially from sqrt of a tiny negative number if h wasn't checked).
    if any(isnan(F_num(:)))
        warning('OsherSolomon:NaNFlux', 'NaN detected in Osher-Solomon flux calculation. Replacing with 0.');
        F_num(isnan(F_num)) = 0; % Replace NaN with 0 as a basic fallback.
    end
    
end % Function OsherSolomon end