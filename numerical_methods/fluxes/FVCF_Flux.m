function Phi = FVCF_Flux(vL, vR, g)
%FVCF_FLUX   Computes the FVCF numerical flux for the shallow water equations.
%
%   Phi = FVCF_Flux(vL, vR, g) returns the numerical flux at cell interfaces
%   using the FVCF scheme. vL and vR are the left and right states, and g is
the gravitational acceleration.
%
%   Inputs:
%       vL - left state vector [h, hu]
%       vR - right state vector [h, hu]
%       g  - gravitational acceleration
%
%   Output:
%       Phi - numerical flux vector
%
%   Reference: Adapted from legacy NumFlux.m

    % Roe-averaged weights
    Hls = sqrt(vL(:,1));
    Hrs = sqrt(vR(:,1));
    
    uL = vL(:,2) ./ (vL(:,1) + eps);
    uR = vR(:,2) ./ (vR(:,1) + eps);
    
    mu1 = 0.5 * (vL(:,1) + vR(:,1));
    mu2 = (Hls .* uL + Hrs .* uR) ./ (Hls + Hrs + eps);
    
    cm = sqrt(g * mu1);
    s1 = sign(mu2 + cm);
    s2 = sign(mu2 - cm);
    
    sd = 0.5 * (s2 - s1) ./ (cm + eps);
    sp = 0.5 * (s1 + s2);
    
    fl = PhysFlux(vL, g);
    fr = PhysFlux(vR, g);
    
    Phi = 0.5 * (fl + fr);
    fd = 0.5 * (fl - fr);
    
    Phi(:,1) = Phi(:,1) + (sd .* mu2 + sp) .* fd(:,1) - sd .* fd(:,2);
    Phi(:,2) = Phi(:,2) + sd .* (mu2.^2 - cm.^2) .* fd(:,1) + (sp - sd .* mu2) .* fd(:,2);
end