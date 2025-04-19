function cfg = flat_bottom_config()

    %FLAT_BOTTOM_CONFIG  Five‑metre, constant‑depth tank, N = 500 cells.
    %
    %   cfg = FLAT_BOTTOM_CONFIG() returns a fully‑initialised configuration
    %   structure that can be passed verbatim to CORE.SOLVER:
    %
    %       results = core.solver( flat_bottom_config() );
    
        cfg = default_config();        % start from safe defaults
                                   % ───────────────────────────
    
        %% Domain and spatial mesh
        cfg.domain.xmin = 0.0;             % left wall [m]
        cfg.domain.xmax = 5.0;             % right boundary [m]
        cfg.mesh.N      = 500;             % number of control volumes
        [cfg.xc, cfg.dx, cfg.x_edge] = uniform(cfg.domain, cfg.mesh.N);
    
        %% Bathymetry and initial condition
        cfg.param.H0    = 0.50;            % undisturbed water depth [m]
        cfg.bathyHandle = @(x) flat(x, cfg);
        cfg.ic          = @(x) lake_at_rest(x, cfg);   % η = 0, u = 0
    
        %% Governing model and numerical options
        % For direct use with ode45, always use a wrapper:
        % cfg.model = @(t, z) RHS_NSWE(t, z, params); % Not used directly in run_simulation.m
        cfg.model = []; % Placeholder: see run_simulation.m for actual wrapper
        cfg.numFlux        = @FVCF_Flux;   % FVCF numerical flux (new flux)
        cfg.reconstruction = [];           % no reconstruction - 1st order method
        cfg.timeStepper    = @integrate_ode45;  % Matlab's standard ode45 wrapper
    
        %% Boundary conditions
        cfg.bcL = @bc.wall;                % rigid wall at x = 0
        cfg.bcR = @bc.outgoing;            % weakly‑reflecting outflow at x = L
    
        %% Run‑control parameters
        cfg.tEnd        = 10.0;            % simulated time span [s]
        cfg.outputEvery = 100;             % snapshot frequency (time steps)
    
        %% House‑keeping
        cfg.caseName    = 'flat_bottom_L5m';      % useful for I/O folders

end