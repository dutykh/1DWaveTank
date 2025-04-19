
clear
close all
format longE
warning('on', 'verbose');

%%% Libraries we use:
addpath('data/');
addpath('sources/');
addpath('tools/');
addpath('tools/MagInset/')
addpath('tools/odetpbar/');
addpath('tools/export_fig/');
addpath('tools/boundedline/');

%%% Some abbreviations:
FS = 'FontSize';
IN = 'Interpreter';
LS = 'LineStyle';
LW = 'LineWidth';
MS = 'MarkerSize';

%%% Global variables:
global binterp cf d
global N g g2 dx d2x h nu

%%% Physical parameters:
g  = 9.81;  % gravity acceleration
g2 = 0.5*g; % g/2
cf = 0.0;   % friction coefficient
nu = 0.0;   % eddy viscosity
%
d  = 0.28;   % undisturbed water depth
%
a  = 0.0;   % the left boundary (incident wave)
st =0.41; b = 6.61;b1 = 01.61;% the point where the slope starts
td = 1/10; td1 = 1/15;%% Numerical parameters:
N  = 50;    x  = linspace(a, b, N+1)';          % cell interfaces
dx = x(2) - x(1);                   % spatial grid step
xc = 0.5*(x(1:end-1) + x(2:end));   % centers of cells
%% Bathymetry function:
height_1plane = (td*(b1 - st));

total_height = height_1plane + td1*(b - b1);

height_2plane = total_height - height_1plane;

h  = d*(xc >= a).*(xc <= st) + (d - td*(xc - st)).*(xc >= st).*(xc <= b1) + (d - total_height + abs((height_2plane).*(-(xc - b)/(b - b1)).^(4/3))).*(xc >= b1).*(xc <= b);

%%% Choice of the initial condition:
w  = zeros(2*N,1);
w(1:N) = max(h, eps+0*h);   % zero initial condition without velocities

%%% Read experimental data:
load dat.mat
fprintf(' Done\n');
t0  = 0.0;                              % initial time
Tef = 30;

%%% We get the first Gauge as boundary data:
irw = dat(:,1)-mean(dat(:,1));
load tim;    % Frequency is 10Hz

binterp = interp1(tim, irw, 'spline', 'pp');

% Time projection:
Tf    = Tef;  % final simulation time (the same as experiment)
M     = 200; % number of time instances where we project solution
tlist2022 = linspace(t0, Tf, M);
dt = (Tf-t0)/M
CFL = dt/dx*sqrt(g*(d + a))
% ODE solver options:
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar);
% [~, sol] = ode23(@RHS_NSWE, tlist2022, w, options);

% % [~, slb] = ode23(@RHS_mPer, tlist2022, w, options);
% % fprintf(' Done\n');

% Rup2022 = zeros(M,1);
% Rpb2022 = zeros(M,1);
% for t=1:M % loop in time
%   ind = find(sol(t,1:N) > 1e-2, 1, 'last');
%   Rup2022(t) = -h(ind);
% %   ind = find(slb(t,1:N) > 1e-2, 1, 'last');
%   Rpb2022(t) = -h(ind);
% end % for t

% % ind   = zeros(4,1);
% % Gauge = [0.6 1.38 2.14];
% % Data  = dat;
% % %%% We treat the wave gauges data (mPer)
% % fprintf('Plotting data...');
% % figure;
% % % set(gcf, 'pos', [1 1 1918 973]);
% % set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% % ha = tight_subplot(3,1, [0.045 0.019]);
% % % find positions of gauges in the mesh and plot it:
% % for i=2:4 
% %   [~, ind(i)] = min(abs(xc - Gauge(i)));
% %   % Data(:,i) = Data(:,i) - mean(Data(:,i));
% %   Data(:,i) = smooth(Data(:,i) - mean(Data(:,i)), 25, 'moving');
% % 
% %   axes(ha(i));
% %   plot(tim, Data(:,i), '.-'), grid off, hold on
% %   plot(tlist2022, slb(:,ind(i)) - h(ind(i))); hold on
% %   plot(tlist2022, sol(:,ind(i)) - h(ind(i))), hold off
% %   axis tight;
% %   xlim([t0 Tf]); % ylim([-0.23 0.26]);
% %   
% %   drawnow
% % end % for i
% % set(gcf, 'color', 'w');
