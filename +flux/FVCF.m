function Phi = FVCF(vL, vR, cfg)

    %FVCF Roe-type numerical flux for the Nonlinear Shallow Water Equations.
    %   Phi = FVCF(vL, vR, cfg) calculates the numerical flux across an
    %   interface using a Roe-type scheme based on Roe averages.
    %
    %   NOTE: This function implements the logic from the user's legacy
    %   'NumFlux.m' file, which was commented as 'FVCF scheme' but uses
    %   Roe averages.
    %
    %   Inputs:
    %       vL  - State vector [H; HU] at the left of the interface.
    %       vR  - State vector [H; HU] at the right of the interface.
    %       cfg - Configuration structure (must contain param.g).
    %
    %   Outputs:
    %       Phi - Numerical flux vector [Phi_H; Phi_HU].
    
    % Extract gravity from config
    g = cfg.param.g;
    
    % Define epsilon for numerical stability (avoid division by zero)
    eps_flux = 1e-10;
        
    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2);
    Hr = vR(:,1); HuR = vR(:,2);
        
    % Calculate primitive variable U (velocity) on left and right
    uL = zeros(size(Hl)); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL);
    uR = zeros(size(Hr)); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR);
        
    % Calculate Roe averages (weights Hls, Hrs; averaged state mu1, mu2)
    Hls = sqrt(Hl); % weights for Roe averaging
    Hrs = sqrt(Hr); % weights for Roe averaging
        
    % Roe averaged depth (mu1 in legacy code) - Note: Simple average used here
    H_roe = 0.5 * (Hl + Hr);
        
    % Roe averaged velocity (mu2 in legacy code)
    u_roe = (Hls .* uL + Hrs .* uR) ./ (Hls + Hrs + eps_flux);
        
    % Roe averaged wave speed (cm in legacy code)
    c_roe = sqrt(g * H_roe);
        
    % Calculate signs of eigenvalues (lambda = u_roe +/- c_roe)
    lambda1 = u_roe + c_roe;
    lambda2 = u_roe - c_roe;
    s1 = sign(lambda1);
    s2 = sign(lambda2);
        
    % Calculate coefficients based on eigenvalue signs (sd, sp in legacy code)
    sd = 0.5 * (s2 - s1) ./ (c_roe + eps_flux); % Avoid division by zero if c_roe is near zero
    sp = 0.5 * (s1 + s2);
        
    % Calculate physical fluxes on the left and right
    fl = core.utils.physical_flux(vL, cfg);
    fr = core.utils.physical_flux(vR, cfg);
        
    % Calculate the numerical flux Phi using the Roe dissipation terms
    Phi = 0.5 * (fl + fr); % Start with the average (centered part)
    fd = 0.5 * (fl - fr);  % Difference part (used in dissipation)
        
    % Add Roe dissipation terms (matches legacy code logic)
    Phi(:,1) = Phi(:,1) + (sd .* u_roe + sp) .* fd(:,1) - sd .* fd(:,2);
    Phi(:,2) = Phi(:,2) + sd .* (u_roe.^2 - c_roe.^2) .* fd(:,1) + (sp - sd .* u_roe) .* fd(:,2);

end