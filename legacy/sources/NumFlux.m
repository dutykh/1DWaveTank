%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, USMB           %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

%%%%% FVCF scheme
function Phi = NumFlux (vL, vR)

    global g
    
    Hls = sqrt(vL(:,1));    % weights for Roe averaging
    Hrs = sqrt(vR(:,1));    % weights for Roe averaging
    
    uL = vL(:,2)./(vL(:,1)+eps);
    uR = vR(:,2)./(vR(:,1)+eps);
    
    mu1 = 0.5*(vL(:,1) + vR(:,1));  % averaged states on faces
    mu2 = (Hls.*uL + Hrs.*uR)./(Hls + Hrs + eps);
    
    cm = sqrt(g*mu1);
    s1 = sign(mu2 + cm);
    s2 = sign(mu2 - cm);
    
    sd = 0.5*(s2 - s1)./(cm + eps);
    sp = 0.5*(s1 + s2);
    
    fl = PhysFlux(vL);
    fr = PhysFlux(vR);
    
    Phi = 0.5*(fl + fr);
    fd = 0.5*(fl - fr);
    
    Phi(:,1) = Phi(:,1) + (sd.*mu2 + sp).*fd(:,1) - sd.*fd(:,2);
    Phi(:,2) = Phi(:,2) + sd.*(mu2.^2 - cm.^2).*fd(:,1) + (sp - sd.*mu2).*fd(:,2);

end % NumFlux ()