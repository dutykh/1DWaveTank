function config = default_config()
% config.default_config: Provides default values for the simulation configuration.
%   config = config.default_config() returns a structure 'config' containing default
%   parameters and function handles for the 1DWaveTank simulation.

config = struct(); % Initialize empty struct

% --- Default Physical Parameters ---
config.param.g = 9.81;      % Gravity acceleration [m/s^2]
config.param.H0 = 0.5;      % Default undisturbed water depth [m]
config.param.Cf = 0.0;      % Default friction coefficient [-]
config.param.nu = 0.0;      % Default eddy viscosity [m^2/s]

% --- Default Domain and Mesh ---
config.domain.xmin = 0.0;
config.domain.xmax = 10.0;   % Default domain length [m]
config.mesh.N = 100;         % Default number of cells

% --- Default Numerical Methods (Placeholders) ---
config.model = []; % To be set in specific config/run script
config.numFlux = @flux.placeholder_flux; % Default placeholder
config.reconstruction = []; % Default: No reconstruction (1st order)
config.well_balancing = []; % Default: No specific WB scheme
config.timeStepper = @time.placeholder_stepper; % Default placeholder

% --- Default Boundary Conditions (Placeholders) ---
config.bcL = @bc.placeholder_bc; % Default placeholder
config.bcR = @bc.placeholder_bc; % Default placeholder

% --- Default Initial Condition ---
config.ic = @ic.lake_at_rest; % Default to lake at rest

% --- Default Bathymetry ---
config.bathyHandle = @bathy.flat; % Default to flat bottom

% --- Default Run Control ---
config.t0 = 0.0;             % Default start time [s]
config.tEnd = 5.0;           % Default end time [s]
config.outputEvery = 10;     % Default output frequency (time steps)
config.useGpu = false;       % Default: Do not use GPU

% --- Default House-keeping ---
config.caseName = 'default_case';
config.outputPath = './results/'; % Default output path

end