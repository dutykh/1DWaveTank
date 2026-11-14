% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%
% Data of the figure of the README: the default experiment of +cfg/simulation_config.m
% ('sloping_beach'), set up and solved exactly as run_simulation.m does it
% (cfg.simulation_config, then core.solver), without the animation. Sine waves of
% amplitude 0.2 m and period 4 s enter at the left end of a 40 m tank of still depth
% 1 m, whose bottom rises with slope 0.1 over its last third, and run up the beach
% until t = 25 s.
%
% Run it from this folder in MATLAB:
%   matlab -batch wavetank_readme_data        or, in the MATLAB window,   run wavetank_readme_data
% It writes wavetank_frames.mat beside itself (every output frame, 0.2 s apart, about
% 3 MB); plot_readme_figure.py then draws wavetank-sloping-beach.png from it. GNU Octave
% has no ode113: there the script integrates with ode45 at the same tolerances.
%   octave --no-gui wavetank_readme_data.m

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
clock0 = tic;

%%% Configuration of the default experiment, as run_simulation.m loads it:
config = cfg.simulation_config();
if exist('ode113') == 0
  config.time.matlab_solver = 'ode45';  % GNU Octave
end

%%% Solution, as run_simulation.m computes it:
results = core.solver(config);

%%% Output frames, bottom and parameters:
t  = results.t(:);                                  % output times [s]
xc = config.mesh.xc(:);                             % cell centres [m]
zb = config.bathyHandle(config, config.mesh.xc);    % bottom elevation z_b [m]
zb = zb(:);
H  = results.H;                                     % depth H, one row per output time [m]
HU = results.HU;                                    % discharge HU [m^2/s]
U  = results.U;                                     % velocity U = HU/H in wet cells [m/s]
a  = config.bc.left.param.a;                        % amplitude of the generated wave [m]
T  = config.bc.left.param.T;                        % period of the generated wave [s]
slope = config.phys.slope;                          % slope of the beach
dry_tolerance = config.phys.dry_tolerance;          % depth below which a cell is dry [m]
solver = config.time.matlab_solver;
runtime = toc(clock0);
if exist('OCTAVE_VERSION', 'builtin') ~= 0
  engine = ['GNU Octave ' OCTAVE_VERSION];
else
  engine = ['MATLAB ' version];
end

save(fullfile(here, 'wavetank_frames.mat'), 't', 'xc', 'zb', 'H', 'HU', 'U', 'a', 'T', ...
     'slope', 'dry_tolerance', 'solver', 'runtime', 'engine', '-v7');
fprintf('%d frames up to t = %.15g s with %s in %.0f s (%s)\n', numel(t), t(end), solver, ...
        runtime, engine);
