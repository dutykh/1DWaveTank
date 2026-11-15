% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the second figure of the README: a hump of the free surface leaving a flat tank through
% its two open ends. The tank 0 ⩽ x ⩽ 20 m holds water 0.5 m deep over the flat bottom z = 0,
% cut into 500 cells; the hump eta - eta0 = 0.1 exp(-2 (x - 10)^2) m starts at rest and splits
% into two waves that leave through the ends. It is computed with the open end of bc.open, which
% takes the outgoing Riemann invariant from the tank and the incoming one from still water at the
% level eta0, and with the extrapolation of the end cell that bc.open applies when no still-water
% level is given (the open end of versions 1.3.0 and earlier), by the second-order scheme (MUSCL
% with the van Leer limiter) and the first-order scheme, HLLE fluxes and the second-order
% strong-stability-preserving Runge-Kutta scheme at CFL = 0.5, up to t = 12 s.
%
% Run it from this folder in MATLAB or GNU Octave:
%   matlab -batch open_end_reflection_data        or        octave --no-gui open_end_reflection_data.m
% It writes open_end_reflection.mat beside itself (the free surface at t = 8 s and the largest
% elevation in the tank every 0.02 s); plot_open_end_reflection.py then draws
% open-end-reflection.png from it.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
clock0 = tic;

models = {'rhs_nsw_2nd_order', 'rhs_nsw_1st_order'};
model_names = {'second', 'first'};
ends = {'characteristic', 'extrapolation'};
eta0 = 0.5;                                              % [m] still-water level
T_snap = 8;                                              % [s] time of the free surface shown
data = struct();

for im = 1:numel(models)
  for ie = 1:numel(ends)
    config = cfg.default_config();
    config.mesh.domain.xmin = 0; config.mesh.domain.xmax = 20; config.mesh.N = 500;
    config.domain = config.mesh.domain;
    config.model = str2func(['core.' models{im}]);
    config.numFlux = @flux.HLLE;
    config.timeStepper = @time.integrate_ssp2_adaptive;
    config.time.cfl = 0.5;
    config.time.num_progress_reports = 1;
    config.reconstruction = [];
    config.reconstruct = struct('handle', @reconstruct.muscl, 'limiter', 'vanleer');
    config.bc.left.handle = @bc.open;
    config.bc.right.handle = @bc.open;
    config.bc.num_ghost_cells = 2;
    config.bathyHandle = @bathy.flat;
    config.bathy_params.flat_elevation = 0;
    config.ic_handle = @ic.gaussian_bump;
    config.ic_param = struct('a', 0.1, 'lambda', 2, 'x0', 10, 'H0', eta0);
    config.t0 = 0; config.tEnd = 12; config.vis.dt_plot = 0.02;
    config.tspan = config.t0:config.vis.dt_plot:config.tEnd;   % [s] output times
    config = cfg.validate_config(config);
    if strcmp(ends{ie}, 'extrapolation')
      % Without a still-water level, bc.open extrapolates the end cell
      config.bc.left.param = rmfield(config.bc.left.param, 'eta0');
      config.bc.right.param = rmfield(config.bc.right.param, 'eta0');
    end
    results = core.solver(config);
    t = results.t(:);
    eta = results.H - eta0;                              % elevation above the still-water level
    k = find(abs(t - T_snap) < 1e-9, 1);
    name = sprintf('%s_%s', model_names{im}, ends{ie});
    data.(['etamax_' name]) = max(abs(eta), [], 2);     % largest |eta - eta0| in the tank
    data.(['eta_snap_' name]) = eta(k, :).';            % eta - eta0 at t = T_snap
    data.(['volume_' name]) = sum(results.H, 2) * config.mesh.dx;
  end
end

data.t = t;
data.xc = config.mesh.xc(:);
data.T_snap = T_snap;
data.eta0 = eta0;
data.runtime = toc(clock0);
if exist('OCTAVE_VERSION', 'builtin') ~= 0
  data.engine = ['GNU Octave ' OCTAVE_VERSION];
else
  data.engine = ['MATLAB ' version];
end
save(fullfile(here, 'open_end_reflection.mat'), '-struct', 'data', '-v7');
fprintf('open_end_reflection.mat written: four runs up to t = %g s in %.0f s (%s)\n', ...
        data.t(end), data.runtime, data.engine);
