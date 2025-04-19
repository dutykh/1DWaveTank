%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

% Right-hand side
function zrhs = RHS_NSWE(t, z)

    global cf g g2 dx h N

    w    = zeros(N,2);
    w(:,1) = z(1:N);
    w(:,2) = z(N+1:end);

    p = zeros(N, 2);
    
    p(:,1) = w(:,1);  % we will reconstruct conservative variables
    p(:,2) = w(:,2)./(w(:,1)+eps);
    
    vL = zeros(N-1,2);    % variables reconstructed from the left
    vR = zeros(N-1,2);    % variables reconstructed from the right
    
    dp = p(2:end,:) - p(1:end-1,:);
    Dp = p(3:end,:) - 2*p(2:end-1,:) + p(1:end-2,:);
    
    Dpl = mmod(Dp(1:end-1,:), Dp(2:end,:));
    
    s = mmod(dp(2:end-2,:)+0.5*Dpl(1:end-1,:), dp(3:end-1,:)-0.5*Dpl(2:end,:));
    
    vL(3:N-2,1) = p(3:N-2,1) + 0.5*s(:,1);
    vR(2:N-3,1) = p(3:N-2,1) - 0.5*s(:,1);
    
    vL(3:N-2,2) = p(3:N-2,2) + 0.5*s(:,2).*vR(2:N-3,1)./(w(3:N-2,1)+eps);
    vR(2:N-3,2) = p(3:N-2,2) - 0.5*s(:,2).*vL(3:N-2,1)./(w(3:N-2,1)+eps);
    
    % treatment of variables near from boundaries (cells: 1,2,N-2,N-1)
    vL(1,:) = p(1,:);     % we cannot reconstruct in the 1st cell from the left
    vR(N-1,:) = p(N,:);   % ... in the last cell from the right !
    
    vR(1,:) = 0.375*p(1,:) + 0.75*p(2,:) - 0.125*p(3,:);
    vL(2,:) = -0.125*p(1,:) + 0.75*p(2,:) + 0.375*p(3,:);
    vL(N-1,:) = -0.125*p(N-2,:) + 0.75*p(N-1,:) + 0.375*p(N,:);
    vR(N-2,:) = 0.375*p(N-2,:) + 0.75*p(N-1,:) - 0.125*p(N,:);
    
    %%% Hydrostatic reconstruction
    hL = h(1:end-1); hR = h(2:end);
    
    hi = min(hL, hR);
    HLs = max(w(1:end-1,1) - hL + hi, zeros(N-1,1));
    HRs = max(w(2:end,1) - hR + hi, zeros(N-1,1));
    
    vL(:,1) = HLs;
    vL(:,2) = HLs.*vL(:,2);
    
    vR(:,1) = HRs;
    vR(:,2) = HRs.*vR(:,2);
    
    %%% Computing the fluxes
    Fnum = NumFlux(vL(:,1:2), vR(:,1:2));
    
    % declaration of the result
    rhs = zeros(N,2);
    
    rhs(1:end-1,:) = rhs(1:end-1,:) - Fnum;
    rhs(2:end,:) = rhs(2:end,:) + Fnum;
    
    rhs(1:end-1,2) = rhs(1:end-1,2) - g2*(w(1:end-1,1).^2 - HLs.^2);
    rhs(2:end,2) = rhs(2:end,2) + g2*(w(2:end,1).^2 - HRs.^2);

    ind = p(:,1) > 1e-2;
    rhs(ind,2) = rhs(ind,2) - g*cf*(p(ind,2).*abs(p(ind,2)))./(p(ind,1).^(1.0/3.0));
    
    % left boundary
    bL = BoundaryValue(t);
    cs = sqrt (g*vL(1,1));
    u = vL(1,2)/(vL(1,1)+eps) + (1 - p(1,1)/bL)*cs;
    bFlux = [bL*u, bL*u^2 + g2*bL^2];
    rhs(1,:) = rhs(1,:) + bFlux;
    
    % right boundary
    bFlux = PhysFlux(w(end,:));
    rhs(N,:) = rhs(N,:) - bFlux;
    
    rhs = rhs/dx; % the value we return
    
    zrhs = zeros(2*N, 1);
	zrhs(1:N) = rhs(:,1);
    zrhs(N+1:end) = rhs(:,2);

end % RHS_NSWE ()