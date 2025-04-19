function F = physical_flux(w, cfg)

    %PHYSICAL_FLUX Calculates the physical flux vector for the Shallow Water Equations.
    %   F = PHYSICAL_FLUX(w, cfg) computes the flux F = [HU; HU^2 + 0.5*g*H^2]
    %   given the state vector w = [H; HU].
    %
    %   Inputs:
    %       w   - State vector(s) [H, HU] (can be N x 2 array).
    %       cfg - Configuration structure containing parameters like cfg.param.g.
    %
    %   Outputs:
    %       F   - Physical flux vector(s) [F_H, F_HU] (N x 2 array).
    
    % Extract parameters
    g = cfg.param.g;
    
    % Extract conserved variables H and HU
    H = w(:, 1);
    HU = w(:, 2);
        
    % Calculate primitive variable U (velocity)
    % Avoid division by zero in dry areas (H=0)
    U = zeros(size(H));
    wet_indices = H > 1e-10; % Define a tolerance for wet cells
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
        
    % Calculate flux components
    F_H = HU;                     % Flux for H equation
    F_HU = HU .* U + 0.5 * g * H.^2; % Flux for HU equation
        
    % Combine into flux vector:
    F = [F_H, F_HU];

end