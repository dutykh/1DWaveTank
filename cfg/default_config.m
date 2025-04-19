function cfg = default_config()
    
    cfg = struct();

    %% ─────────────────────────────────────────────────────────────────────
    %% Physical domain & mesh
    cfg.domain.xmin = 0.0;         % left boundary [m]
    cfg.domain.xmax = 1.0;         % right boundary [m]

    cfg.mesh.N        = 100;       % number of control volumes
    cfg.mesh.generator = @grid.uniform;  % handle returning (xc,dx,x_edge)

    %% ─────────────────────────────────────────────────────────────────────
    %% Physical parameters (gravity, density, viscosity, …)
    cfg.param.g   = 9.81;          % gravitational acceleration [m s⁻²]
    cfg.param.H0  = 1.0;           % undisturbed water depth [m]
    cfg.param.nu  = 0.0;           % kinematic viscosity [m² s⁻¹] (0 → inviscid)

    %% ─────────────────────────────────────────────────────────────────────
    %% Governing PDE model & numerical choices
    cfg.model          = @rhs.nswe;     % nonlinear shallow‑water equations
    cfg.numFlux        = @flux.roe;     % Roe approximate Riemann solver
    cfg.reconstruction = @recon.uno2;   % 2‑nd‑order UNO (robust default)
    cfg.timeStepper    = @time.rk4;     % classical RK4 integrator

    cfg.wellBalancing  = false;         % lake‑at‑rest correction off by default

    %% ─────────────────────────────────────────────────────────────────────
    %% Bathymetry, ICs and BCs
    cfg.bathyHandle = @(x) bathy.flat(x, cfg);          % constant depth
    cfg.ic          = @(x) ic.lake_at_rest(x, cfg);     % η = 0, u = 0

    cfg.bcL = @bc.periodic;         % default: periodic domain
    cfg.bcR = @bc.periodic;

    %% ─────────────────────────────────────────────────────────────────────
    %% Time‑stepping & run control
    cfg.cfl          = 0.45;        % Courant number (dt chosen by solver)
    cfg.tEnd         = 1.0;         % final simulation time [s]
    cfg.outputEvery  = 1;           % diagnostics frequency (steps)

    %% The solver may populate these run‑time fields:
    cfg.Nt           = [];          % number of time steps actually taken
    cfg.stats        = struct();    % execution statistics (CPU, CFL, …)

    %% ─────────────────────────────────────────────────────────────────────
    %% I/O & book‑keeping
    cfg.caseName          = 'default_case';
    cfg.io.resultsDir     = fullfile('results', cfg.caseName);
    cfg.io.saveSolution   = false;  % turn on to dump snapshots to MAT files

    %% ─────────────────────────────────────────────────────────────────────
    %% Derived quantities: build mesh immediately so they are always present
    [cfg.xc, cfg.dx, cfg.x_edge] = cfg.mesh.generator(cfg.domain, cfg.mesh.N);

end