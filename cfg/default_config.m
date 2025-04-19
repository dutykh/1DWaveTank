function cfg = default_config()

    %DEFAULT_CONFIG Provides default values for the simulation configuration.
    %   cfg = DEFAULT_CONFIG() returns a structure 'cfg' containing default
    %   parameters and function handles for the 1DWaveTank simulation.
    
    cfg = struct(); % Initialize empty struct
    
    % --- Default Physical Parameters ---
    cfg.param.g = 9.81;      % Gravity acceleration [m/s^2]
    cfg.param.H0 = 0.5;      % Default undisturbed water depth [m]
    cfg.param.Cf = 0.0;      % Default friction coefficient [-]
    cfg.param.nu = 0.0;      % Default eddy viscosity [m^2/s]
    
    % --- Default Domain and Mesh ---
    cfg.domain.xmin = 0.0;
    cfg.domain.xmax = 10.0;   % Default domain length [m]
    cfg.mesh.N = 100;         % Default number of cells
    
    % --- Default Numerical Methods (Placeholders) ---
    cfg.model = []; % To be set in specific config/run script
    cfg.numFlux = @flux.placeholder_flux; % Default placeholder
    cfg.reconstruction = []; % Default: No reconstruction (1st order)
    cfg.well_balancing = []; % Default: No specific WB scheme
    cfg.timeStepper = @time.placeholder_stepper; % Default placeholder
    
    % --- Default Boundary Conditions (Placeholders) ---
    cfg.bcL = @bc.placeholder_bc; % Default placeholder
    cfg.bcR = @bc.placeholder_bc; % Default placeholder
    
    % --- Default Initial Condition ---
    cfg.ic = @ic.lake_at_rest; % Default to lake at rest
    
    % --- Default Bathymetry ---
    cfg.bathyHandle = @cfg.bathy.flat; % Default to flat bottom
    
    % --- Default Run Control ---
    cfg.t0 = 0.0;             % Default start time [s]
    cfg.tEnd = 5.0;           % Default end time [s]
    cfg.outputEvery = 10;     % Default output frequency (time steps)
    cfg.useGpu = false;       % Default: Do not use GPU
    
    % --- Default House-keeping ---
    cfg.caseName = 'default_case';
    cfg.outputPath = './results/'; % Default output path
    
end