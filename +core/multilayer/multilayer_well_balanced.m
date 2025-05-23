function multilayer_well_balanced()

    % Start with a clean environment to avoid conflicts from previous runs.
    clc; clear; close all; format longE;
    rehash toolboxcache; % Force MATLAB to update its cache for new functions/packages

    % Add the project root directory and all its subdirectories to the MATLAB path.
    % This ensures all package functions (e.g., `+core`, `+cfg`) are accessible.
    addpath(genpath(pwd));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set random number generator (for the initial condition):
    if not(exist('OCTAVE_VERSION','builtin'))
        s = RandStream('mt19937ar','Seed',13);
        RandStream.setGlobalStream(s);
    else
        disp('Running in OCTAVE')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Environment Setup                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Optional: Suppress OpenGL warnings if they occur on specific systems.
    % warning('off', 'MATLAB:opengl:SwitchToSoftwareOpenGL');

    % Start with a clean environment to avoid conflicts from previous runs.
% % %     rehash toolboxcache; % Force MATLAB to update its cache for new functions/packages

    % Add the project root directory and all its subdirectories to the MATLAB path.
    % This ensures all package functions (e.g., `+core`, `+cfg`) are accessible.
    addpath(genpath(pwd));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status for user clarity (optional, can be commented)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('--- Setting up Simulation Configuration ---\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Load Default Configuration ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start with the baseline parameters defined in default_config.m
    config = cfg.default_config();
    fprintf('Default config loaded. Overriding for specific experiment...\n');

    % ======================================================================
    % --- Common Settings (can be overridden by specific setups below) ---
    % ======================================================================
    % These are the default values for domain, mesh, numerics, etc.
    % They can be overridden in the switch block for each experiment.

    config.caseName = 'gaussian_bump_rest_L20m_H0.5m_N500';
    
    % --- Domain ---
    config.domain.xmin   = km_to_m(0.00);    % [m] Left endpoint of domain
    config.domain.xmax   = km_to_m(35.0);    % [m] Right endpoint of domain
    
    % --- Mesh ---
    config.mesh.domain.xmin = config.domain.xmin; % [m] Left endpoint of domain
    config.mesh.domain.xmax = config.domain.xmax; % [m] Right endpoint of domain
    config.mesh.N           = 2^11;               % [integer] Number of spatial cells
    config.mesh.Nlayers     = 4;                  % [integer] Number of vertical layers
    config.mesh.xc          = linspace(config.mesh.domain.xmin + 0.5*(config.mesh.domain.xmax-config.mesh.domain.xmin)/config.mesh.N, config.mesh.domain.xmax - 0.5*(config.mesh.domain.xmax-config.mesh.domain.xmin)/config.mesh.N, config.mesh.N);     % Guarantee mesh.xc (cell centers) for ICs and bathymetry
    config.mesh.dx          = (config.mesh.domain.xmax - config.mesh.domain.xmin) / config.mesh.N;     % Ensure mesh.dx is set for CFL and solver compatibility
    
    % --- Common parameters ---
    config.L             = config.mesh.domain.xmax - config.mesh.domain.xmin; % Use mesh.domain
    config.param.H0      = 1.0;
    config.h0            = config.param.H0;
    
    % --- Vertical domain Mesh ---
    config.atmosphere.Pstart            = hectoPascal_to_pascal(900);
    config.atmosphere.Pfinal            = hectoPascal_to_pascal(500);
    config.atmosphere.Zj                = linspace(pressurealt(config.atmosphere.Pstart),pressurealt(config.atmosphere.Pfinal),config.mesh.Nlayers);  % [m] Altitude of the j-th level
    config.atmosphere.Zj                = [200 400 600 800];  % [m] Altitude of the j-th level
    config.atmosphere.Altitudes         = [0 config.atmosphere.Zj];    
    config.atmosphere.dH                = diff(config.atmosphere.Altitudes);
    [config.atmosphere.rho,~,~,~,~,~,~] = atmos(config.atmosphere.Altitudes(2:config.mesh.Nlayers + 1) ,'altType','geopotential','units','SI');

    % --- Time Integration ---
    config.timeStepper       = @integrate_rk4_adaptive_local; % Use RK4 for 4th order
    config.cfl_target        = 0.95;                          % CFL for RK4 (can be higher than SSP)
    config.time.cfl          = 0.75;                          % [unitless] CFL number for adaptive time stepping
    config.t0                = 0.0;
    config.tEnd              = 1500.0;
    config.vis.dt_plot       = 5;
    config.vis.plot_velocity = true;                          % [logical] Plot velocity in a subpanel
    config.vis.show_legend   = false;                         % [logical] Show legend in wave tank plot 
    
    % Ensure tspan always contains at least two points (start and end)
    if config.tEnd == config.t0
        config.tspan = [config.t0 config.tEnd];
    else
        config.tspan = config.t0:config.vis.dt_plot:config.tEnd;
    end
    
    % --- Output ---
    config.save_results = false;                   % [logical] Save results by default
    config.output_path = './results/';             % [char] Output directory   
    
    % --- Boundary condition ---
    config.bc.left.handle  = @bc.periodic;
    config.bc.right.handle = @bc.periodic;
    
    % Bathymetry: Gaussian bump
    config.bathyHandle = @bathy.gaussian_bump;
    config.bathy_bump_center = (config.mesh.domain.xmin + config.mesh.domain.xmax)/2;   % Use mesh.domain
    config.bathy_bump_height = 90.0;                                                     % [m] Bump height
    config.bathy_bump_width  = (config.mesh.domain.xmax - config.mesh.domain.xmin)/15;  % Use mesh.domain
    
    % --- Update RHS to high-order version ---
    config.model = @rhs_nsw_high_order_local;
    
    config.reconstruction.method = 'weno5';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- Initial Conditions ---
%     config.ic_handle = @(cfg) ic.gaussian_bump(cfg.mesh.xc, struct('a',2.5*cfg.bathy_bump_height,'lambda',1/(cfg.bathy_bump_width^2),'H0',cfg.h0,'x0',cfg.bathy_bump_center));
    
    config.ic_at_rest = 0;
    config.ic_handle  = @(cfg) traveling_wave(cfg.mesh.xc, config.atmosphere.dH, config.atmosphere.rho, config.phys.g, 7, 6.0, cfg, config.ic_at_rest);
    
    config.numFlux = @flux.HLLC;                   % Riemann solver (HLLC)
%     config.numFlux = @flux.HLL;                   % Riemann solver (HLLC)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Reconstruction Configuration ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch config.reconstruction.method
        
        case 'weno5'
            config.reconstruction.order          = 5;            
            config.reconstruction.handle         = @reconstruct.weno5;
            config.reconstruction.limiter_handle = reconstruct.limiters.limiter_selector(config.reconstruction.limiter);
            
            config.bc.num_ghost_cells            = 3;
            
        case 'muscl'

            config.reconstruction.order          = 2;
            config.reconstruction.handle         = @reconstruct.muscl;
            config.reconstruction.limiter_handle = reconstruct.limiters.limiter_selector(config.reconstruction.limiter);
            
            config.bc.num_ghost_cells            = 2;
            
        case default
            
            error('Reconstruction method!!!!!');
            
    end
    
    config.reconstruction.limiter        = 'vanleer';
    config.reconstruction.limiter        = 'minmod';
    
    config.reconstruction.characteristic = true; % Default: characteristic-wise for all methods

    config.reconstruction.mp5_mode  = 'characteristic';
    config.reconstruction.ppm_mode  = 'characteristic';
    config.reconstruction.weno_mode = 'characteristic';
    config.reconstruction.uno2_mode = 'characteristic';    
    
    config.reconstruct       = config.reconstruction;        % Ensure compatibility with core solver
    config.reconstruct.theta = 1/3;                          % Third-order accuracy in smooth regions
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output Directory and File Setup (Optional)                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Create output directory and define save path if saving results is enabled.
    if isfield(config, 'save_results') && config.save_results
        if ~isfield(config, 'outputPath') || isempty(config.outputPath)
            config.outputPath = './results/'; % Default output directory if not specified
            warning('Output path not specified in config, using default: %s', config.outputPath);
        end
        if ~isfolder(config.outputPath)
            mkdir(config.outputPath);
            fprintf('Created results directory: %s\n', config.outputPath);
        end
        % Generate a filename incorporating a timestamp for uniqueness.
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('results_%s.h5', timestamp);
        savePath = fullfile(config.outputPath, filename);
        fprintf('Results will be saved to: %s\n', savePath);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print Key Configuration Details                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display selected parameters to the command window for user verification.

    fprintf('  Experiment name: %s\n', config.caseName);
    
    % --- Reconstruction ---
    fprintf('  Reconstruction method: %s\n', config.reconstruction.method);    
    fprintf('  Reconstruction limiter: %s\n', config.reconstruction.limiter);
    fprintf('  Reconstruct characteristics: %s\n', is_true_input(config.reconstruction.characteristic));    
    
    % --- Boundary Conditions --- 
    % Print left BC handle and parameters (if any)
    left_bc_handle_str = func2str(config.bc.left.handle);
    if isfield(config.bc.left, 'param') && ~isempty(fieldnames(config.bc.left.param))
        params = config.bc.left.param;
        param_names = fieldnames(params);
        % Use cellfun to format each parameter as 'name=value'
        param_strs = cellfun(@(name) sprintf('%s=%.3g', name, params.(name)), param_names, 'UniformOutput', false);
        param_str = strjoin(param_strs, ', '); % Join into a single string
        fprintf('  BC Left: %s (Params: %s)\n', left_bc_handle_str, param_str);
    else
        fprintf('  BC Left: %s\n', left_bc_handle_str);
    end

    % Print right BC handle (assuming simple BCs like wall/open often don't need params printed)
    right_bc_handle_str = func2str(config.bc.right.handle);
    fprintf('  BC Right: %s\n', right_bc_handle_str);

    % --- Numerics and Time --- 
    fprintf('  Numerical Flux: %s\n', func2str(config.numFlux));
    fprintf('  Time Stepper: %s\n', func2str(config.timeStepper));
    if isfield(config, 'time') && isfield(config.time, 'cfl')
        fprintf('  Time Span: [%.2f, %.2f] s, CFL: %.2f\n', config.t0, config.tEnd, config.time.cfl);
    else
        fprintf('  Time Span: [%.2f, %.2f] s (CFL not applicable)\n', config.t0, config.tEnd);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the Core Solver                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The `core.solver` function encapsulates the main simulation logic.
    % It takes the configuration structure `config` and executes the time stepping
    % loop according to the specified numerical methods and parameters.
    % It returns a `results` structure containing the simulation output (time vector,
    % state variables H, HU, U, etc.) and statistics.

    fprintf('--- Running Core Solver ---\n');
    tic; % Start timer to measure solver execution time.

    results = solver_multilayer(config); % Execute the main simulation function.

    cpu_time = toc; % Stop timer and get elapsed time.
    fprintf('--- Core Solver Finished (CPU Time: %.3f s) ---\n', cpu_time);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print Execution Statistics                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('--- Simulation Statistics ---\n');
    fprintf('  Mesh Cells (N) : %d\n', config.mesh.N);
    if isfield(results, 'total_steps') && ~isnan(results.total_steps)
        fprintf('  Total Steps    : %d\n', results.total_steps);
    else
        fprintf('  Total Steps    : N/A (Using MATLAB ODE Solver or info missing)\n');
    end
    fprintf('  CPU Time       : %.3f s\n', cpu_time);
    fprintf('-----------------------------\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Visualization / Animation                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if the solver returned valid results before attempting to plot.
    if isfield(results, 't') && ~isempty(results.t) && isfield(results, 'H') && ~isempty(results.H)

        fprintf('--- Starting Visualization ---\n');

        %%% Get Bathymetry
        % Calculate bathymetry `h(x)` at cell centers `xc` using the function handle from config.
        % Ensure bathymetry is a row vector for consistent subtraction later.
        h_bathy = config.bathyHandle(config, config.mesh.xc); 
        if iscolumn(h_bathy); h_bathy = h_bathy'; end % Ensure row vector

        num_time_steps = length(results.t); % Number of output frames
        fig = []; % Initialize figure handle (plot_state will create if needed)

        %%% Compute Global Axis Limits
        % Calculate axis limits *before* the loop to ensure consistent axes across all animation frames.
        % This prevents the axes from rescaling dynamically, which can be distracting.

        % --- Y-Limits for Surface/Bathy Plot --- 
        eta_all = results.H + h_bathy; % Free surface = Total Water Depth H + Bottom Elevation z_b(x) (relative to z=0 datum)

        y_min_data = min(min(h_bathy(:)), min(eta_all(:))); % Min of bottom and surface
        y_max_data = max(max(h_bathy(:)), max(eta_all(:))); % Max of bottom and surface

        range = y_max_data - y_min_data;
        if range < 1e-6; range = 1; end % Avoid zero range for flat cases
        padding = 0.1 * range; % 10% padding

        y_limits = [y_min_data - padding, y_max_data + padding];

        % --- X-Limits --- 
        x_limits = [min(config.mesh.xc), max(config.mesh.xc)];

        % --- Y-Limits for Velocity Plot --- 
        if isfield(results, 'U') && ~isempty(results.U) % Check if velocity was calculated
            u_min = min(results.U(:));
            u_max = max(results.U(:));
            if u_min == u_max % Handle case of zero or constant velocity
                delta = max(abs(u_min), 1e-2) * 0.1; % 10% margin or at least 0.001
                u_limits = [u_min - delta, u_max + delta];
            else
                u_margin = 0.1 * (u_max - u_min); % 10% margin
                u_limits = [u_min - u_margin, u_max + u_margin];
            end
        else
            u_limits = [-1, 1]; % Default limits if U is not available
        end

        %%% Animation Loop
        % Iterate through each output time step stored in the results.
        for idx = 1:num_time_steps
            current_t = results.t(idx);         % [s] Time for the current frame
            current_H = results.H(idx, :)';     % [N x 1, m] Water depth (needs column vector for plot_state)
            if isfield(results, 'U') && ~isempty(results.U)
                 current_U = results.U(idx, :)'; % [N x 1, m/s] Velocity (needs column vector)
            else
                 current_U = nan(size(current_H)); % Use NaN if velocity is not plotted/available
            end
            current_h_bathy = h_bathy';         % [N x 1, m] Bathymetry (needs column vector)

            % Call the plotting function to draw/update the figure.
            % Pass the existing figure handle `fig` to update the same window.
            fig = vis.plot_state(config.mesh.xc, current_H, current_h_bathy, current_U, current_t, config, fig, x_limits, y_limits, u_limits);

            drawnow; % Force MATLAB to render the plot immediately.
            pause(0.05); % Pause briefly to control animation speed.

            % --- Optional: Save Frames for Movie --- 
            % Uncomment and configure this section to save each frame as an image file.
            % These frames can later be compiled into a movie using tools like ffmpeg.
            % movie_dir = fullfile(config.outputPath, 'frames');
            % if ~isfolder(movie_dir), mkdir(movie_dir); end
            % frame_filename = fullfile(movie_dir, sprintf('frame_%04d.png', idx));
            % saveas(fig, frame_filename); % Or use print() for more control
        end
        fprintf('--- Visualization Finished ---\n');
    else
        fprintf('--- Skipping Visualization (No valid results data) ---\n');
    end        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cleanup (Optional)                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove added paths if desired (might be useful in some contexts)
    % rmpath(genpath(pwd));

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    % End of Main  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    function Pa = hectoPascal_to_pascal(hPa)
        Pa = 100.0 .* hPa;
    end

    function m = km_to_m(km)
        m = 1000.0 .* km;
    end

    function Z_j = interaction(N, j, rho, h)
      % Input:
      %   N : Numero totale di livelli.
      %   j : Livello per cui calcolare la somma.
      %   rho: Array (1xN) contenente i valori di rho per ogni livello.
      %   h  : Matrice (NxM) contenente i valori di h (altezza) per ogni livello e punto spaziale.
      %
      % Output:
      %   Z_j : Array (1xM) contenente la somma calcolata per il livello j.

      M = size(h, 2);      % Dimensione spaziale
      Z = zeros(1, M);      % Z = zeros(1, M)

      somma1 = zeros(1, M);
      for k = 1:j-1
        if k >= 1
          somma1 = somma1 + h(k, :);
        end
      end

      somma2 = zeros(1, M);
      for k = j+1:N
        if k <= N
          somma2 = somma2 + (rho(k) / rho(j)) * h(k, :);
        end
      end

      Z_j = Z + somma1 + somma2;

    end

    function result = is_true_input(input)
    %IS_TRUE_INPUT Restituisce 'true' o 'false' come stringa
    %   Accetta input logici, numerici o stringhe come 'true', 'yes', 'on', ecc.

        if islogical(input)
            value = input;

        elseif isnumeric(input)
            value = input ~= 0;

        elseif ischar(input) || isstring(input)
            str = lower(string(input));
            true_values = ["true", "yes", "on", "1"];
            value = any(str == true_values);

        else
            value = false;
        end

        % Restituisci come stringa
        if value
            result = "true";
        else
            result = "false";
        end
    end

    function [w0] = traveling_wave(x, H0, rho, g, nwaves, max_u, cfg, rest)
        % Input:
        %   x        : Vettore coordinate spaziali (1xM).
        %   H0       : Vettore (1xN) altezze indisturbate.
        %   rho      : Vettore (1xN) densità strati.
        %   g        : Accelerazione gravità.
        %   nwaves   : Numero di onde nel dominio
        %   cfg      : Opzioni solver
        %   rest     : true per lake at rest condition
        %
        % Output:
        %   w0       : Matrice altezze e velocità iniziali (NlayersxNx2).
        %
        % Nota: Si assume che rho(1) sia la densità massima e che la densità diminuisca con l'aumentare di j.
        
        w0 = zeros(cfg.mesh.Nlayers, cfg.mesh.N, 2);         
        
        h0 = zeros(cfg.mesh.Nlayers, cfg.mesh.N);
        u0 = zeros(cfg.mesh.Nlayers, cfg.mesh.N);
        
        z_b = cfg.bathyHandle(cfg, x);
        
        if ~rest                        
            
            for j = 1:cfg.mesh.Nlayers
                
                fprintf('Generating layer %d of %d...', j, cfg.mesh.Nlayers);
                
                A = 0.00;           % Inizializza l'ampiezza
                Vref_layer = max_u; % Massima velocità da raggiungere
                
                lunghezza_onda = max(x) / nwaves; % Alcune onde nel dominio
                
                etaj = zeros(1, cfg.mesh.N);
                
                uj_max = 0;
                
                while uj_max <= Vref_layer
                    
                    % Definisci Ampiezza della perturbazione eta
                    A = A + 1e-3; % Incrementa l'ampiezza
                    
                    etaj = A * sin(2 * pi * x / lunghezza_onda);
                       
                    if j == 1
                        eta_zeta = etaj - z_b;
                    else
                        eta_zeta = etaj - 0*(h0(j-1, :) - H0(j-1)) - z_b;
                    end
                    h0(j, :) = H0(j) + eta_zeta + cfg.phys.dry_tolerance;
                    
                    w0(j,:,1) = h0(j, :);
                    
                    if j == 1
                        H_sotto = 0;
                        rho_sotto = rho(1); %atmos(0.0,'altType','geopotential','units','SI');
                    else
                        H_sotto = H0(j-1);
                        rho_sotto = rho(j-1);
                    end
                    
                    H_tot = sum(H0);
                    %                 H_tot = H_sotto + H0(j);
                    g_ridotto = g * (rho_sotto - rho(j)) / rho_sotto; % Differenza di densità con lo strato sottostante
                    radicando = 1 - (4 * g_ridotto * H_sotto * H0(j) / (g * H_tot^2));
                    c_quadro = g * H_tot * (0.5 + 0.5 * sqrt(radicando)); % Modo barotropico (+)
                    
                    %                 c_quadro = g.*H0(j);
                    c_j = sqrt(c_quadro);
                    u0(j, :) = c_j .* etaj ./ h0(j, :); % Modo barotropico (+)
                    
                    uj_max = max(abs(u0(j, :)));
                    
                end
                
                fprintf(' amplitude: %f [m] --- max |u|: %f [m/s] --- dH: %f [m]--- quota H: %f [m]\n', max(abs(etaj)), max(abs(u0(j, :))), H0(j), sum(H0(1:j)));
                
                w0(j, :, 2) = u0(j, :).*h0(j, :);
                
            end
            
        else
            
            for j = 1:cfg.mesh.Nlayers
                if j == 1
                    w0(j, :, 1) = H0(j) - z_b;
                else
                    w0(j, :, 1) = H0(j) - (w0(j-1, :, 1) - H0(j-1)) - z_b;
                end
                w0(j, :, 2) = 0.0;
            end
            
        end
        
        w0 = reshape(w0,cfg.mesh.Nlayers*cfg.mesh.N*2,1);
        
    end

    function results = solver_multilayer(cfg)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % +core/solver.m
        %
        % Purpose:
        %   Main simulation driver for the 1DWaveTank code. Orchestrates the
        %   simulation process: sets up initial conditions, selects the appropriate
        %   RHS function and time integrator, runs the simulation, and processes output.
        %
        % Syntax:
        %   results = solver(cfg)
        %
        % Inputs:
        %   cfg - [struct] Configuration structure containing all simulation parameters:
        %         cfg.mesh: Mesh details (N, x, xc, dx)
        %         cfg.time: Time integration parameters (t_span, integrator handle, dt_plot, cfl)
        %         cfg.phys: Physical parameters (g)
        %         cfg.prob: Problem-specific parameters (initial conditions handle, etc.)
        %         cfg.numerics: Numerical scheme details (rhs handle, flux handle, etc.)
        %         cfg.bc: Boundary condition handles
        %
        % Outputs:
        %   results - [struct] Structure containing the simulation output:
        %             results.t:    [M_out x 1] Column vector of output time points
        %             results.H:    [M_out x N] Matrix of water depth H at cell centers
        %             results.HU:   [M_out x N] Matrix of discharge HU at cell centers
        %             results.U:    [M_out x N] Matrix of velocity U at cell centers
        %             results.xc:   [N x 1] Vector of cell center coordinates
        %             results.cfg:  [struct] The configuration structure used for the run
        %             results.total_steps: [integer] Total number of time steps taken
        %
        % Dependencies:
        %   Expects properly configured cfg structure and function handles for IC, RHS, integrator, etc.
        %
        % References:
        %   - See 1DWaveTank UserGuide.md for structure and usage.
        %
        % Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
        % Date:   21 April 2025
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Input Validation (Basic)                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('--- Starting Core Solver ---\n');
        required_fields = {'tspan', 'model', 'timeStepper', 'mesh', 'ic_handle', 'ic_param'};
        for i = 1:length(required_fields)
            if ~isfield(cfg, required_fields{i})
                error('Configuration structure `cfg` is missing required field: %s', required_fields{i});
            end
        end
        if ~isfield(cfg.mesh, 'N')
            error('Configuration structure `cfg.mesh` is missing required field: N');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initial Condition Setup                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ensure mesh center coordinates are a row vector (1xN)
        if isfield(cfg.mesh, 'xc')
            if iscolumn(cfg.mesh.xc)
                cfg.mesh.xc = cfg.mesh.xc.'; % Transpose to row vector
                fprintf('Transposed cfg.mesh.xc to row vector for compatibility.\n');
            end
        end
        fprintf('Setting up initial condition...\n');
        % The initial condition function should return the state vector [H; HU]
        w_init = cfg.ic_handle(cfg); % Pass full config to IC handle (refactored for lake_at_rest)

        % --- Reshape Initial Condition --- 
        N = cfg.mesh.N;
        Nlayers = cfg.mesh.Nlayers;
        
        if isvector(w_init) && length(w_init) == 2*N*Nlayers
            % Already in [H; HU] format (flattened)
            w0 = w_init(:);
        elseif isvector(w_init) && length(w_init) == N*Nlayers
            % Only H is provided, assume HU = 0
            w0 = [w_init(:); zeros(N,1)];
        elseif size(w_init,2) >= 2 && size(w_init,1) == N
            % Use only the first two columns (H and HU provided, ignore extras)
            w0 = [w_init(:,1); w_init(:,2)];
        else
            error('Initial condition function must return N x 1, N x 2 (or more), or 2*N x 1 array.');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time Span Preparation                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tspan = cfg.tspan;
        if isempty(tspan) || length(tspan) < 2
            error('cfg.tspan must contain at least start and end times.');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare Function Handles                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the function handle for the time stepper
        time_stepper = cfg.timeStepper;

        % Handle for the RHS function (differs based on solver type)
        if isequal(time_stepper, @time.integrate_matlab_ode)
            % MATLAB solvers need f(t,w), so wrap cfg.model to include cfg
            rhs_handle = @(t, w) cfg.model(t, w, cfg);
        else
            % Our custom solvers will receive cfg separately and call f(t,w,cfg)
            rhs_handle = cfg.model; 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time Integration                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic; % Start timer for integration
        % Call the selected time stepper with the appropriate RHS handle
        [sol_out, t_out, stats] = time_stepper(rhs_handle, tspan, w0, cfg);
        integration_time = toc; % Stop timer
        total_steps = stats.nsteps;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Post-processing & Output                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Time integration completed in %.2f seconds.\n', integration_time);
        fprintf('Total time steps taken: %d\n', total_steps);

        % Check output dimensions (sol_out should have time points as rows)
        if length(t_out) ~= size(sol_out, 1)
            error('Dimension mismatch: length(t_out)=%d vs size(sol_out,1)=%d.', length(t_out), size(sol_out, 1));
        end

        % Prepare results structure
        fprintf('Processing results...\n');
        N = cfg.mesh.N; % Number of spatial cells
        Nlayers = cfg.mesh.Nlayers; % Number of spatial cells
        M_out = length(t_out); % Number of output time steps

        % Reshape the flat solution vector returned by the time stepper
        % sol_out is expected to be M_out x (2*N)
        if size(sol_out, 2) ~= 2*N*Nlayers
            error('Time stepper returned solution array with unexpected dimensions.');
        end

        results = struct();
        results.t = t_out(:); % Ensure time is a column vector

        % Extract H and HU (M_out x N)
        results.H = sol_out(:, 1:N);
        results.HU = sol_out(:, N+1:2*N);

        % Calculate velocity U (handle division by zero)
        results.U = zeros(M_out, N);
        wet_mask = results.H > 1e-10; % Avoid division by zero in dry cells
        results.U(wet_mask) = results.HU(wet_mask) ./ results.H(wet_mask);

        % Store the configuration used for this run
        results.cfg = cfg;

        % Store total steps and dt history (if available)
        results.total_steps = total_steps;
        % results.dt_history = stats.dt_history; % dt_history is no longer reliably available

        fprintf('--- Core Solver Finished ---\n');

    end

    function dwdt_flat = rhs_nsw_high_order_local(t, w_flat, cfg)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % +core/rhs_nsw_high_order_local.m
        %
        % Purpose:
        %   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
        %   equations using a high-order finite volume scheme. This version implements
        %   a well-balanced scheme based on hydrostatic reconstruction of linearly
        %   reconstructed variables and a centered source term, as described in
        %   Audusse et al. (2004) and related literature (e.g., "SchemaHydroWB.pdf").
        %
        % Syntax:
        %   dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)
        %
        % Inputs:
        %   t       - [scalar, double] Current simulation time [s].
        %   w_flat  - [2N x 1, double] Flattened state vector [H1;...;HN; HU1;...;HUN].
        %   cfg     - [struct] Configuration structure.
        %
        % Outputs:
        %   dwdt_flat - [2N x 1, double] Flattened time derivative vector [dH/dt; dHU/dt].
        %
        % Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
        % Date:   April 24, 2025
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Extract parameters
        ng = cfg.bc.num_ghost_cells;
        N = cfg.mesh.N;
        g = cfg.phys.g;
        dry_tol = cfg.phys.dry_tolerance;
        dx = cfg.mesh.dx;

        Nlayers = cfg.mesh.Nlayers;
        
        rho = cfg.atmosphere.rho;
        
        % Reshape state vector from flat to structured
        w_flat = reshape(w_flat,Nlayers,N,2);
        
        dwdt = zeros(Nlayers,N,2);
        
        W_padded_all = zeros(Nlayers,N+2*ng,2);
        
        for nl = 1:Nlayers
            w = squeeze(w_flat(nl,:,:));
            
            % Apply boundary conditions with ghost cells
            w_padded = zeros(N + 2*ng, 2);
            w_padded(ng+1 : N+ng, :) = w;
            w_padded = apply_boundary_conditions(w_padded, t, cfg);
            
            W_padded_all(nl,:,:) = w_padded(:,:);
            
        end
        
        for nl = 1:Nlayers
            w = squeeze(w_flat(nl,:,:));

            w_padded = squeeze(W_padded_all(nl,:,:));
                       
            % Interaction between layers
            z_interaction = interaction(Nlayers, nl, rho, squeeze(W_padded_all(:,:,1)));
            
            % Get bathymetry at all cell centers
            x_cell_centers = get_all_cell_centers(cfg, ng);
            z_b = z_interaction + cfg.bathyHandle(cfg, x_cell_centers);
            
            % Perform reconstruction to get interface values
            % Inline reconstruction logic (was perform_reconstruction)
            if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'handle') && ~isempty(cfg.reconstruct.handle)
                reconstruct_handle = cfg.reconstruct.handle;
            else
                reconstruct_handle = @reconstruct.none;
            end
            [wL_int, wR_int] = reconstruct_handle(w_padded, cfg);
            
            % Get additional topography at interfaces for source term
            z_interfaces = compute_interface_topo(z_b, cfg, ng, N);
            
            % Apply hydrostatic reconstruction
            [wL_hydro, wR_hydro] = hydrostatic_reconstruction(wL_int, wR_int, z_interfaces, dry_tol, g);
            
            % Calculate numerical fluxes
            F_num = cfg.numFlux(wL_hydro, wR_hydro, cfg);
            
            % Compute flux divergence term: -(F_{i+1/2} - F_{i-1/2})/dx
            dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx;
            
            % Compute well-balanced source term using exact discrete balancing
            dwdt_source = zeros(N, 2);
            for i = 1:N
                idx_nl = i + ng;
                h_i = w_padded(idx_nl, 1);

                % Skip dry cells
                if h_i <= dry_tol
                    continue;
                end

                % Source term that exactly balances the numerical flux for lake-at-rest
                z_diff = z_interfaces(i+1) - z_interfaces(i);
                dwdt_source(i, 2) = -g * h_i * z_diff / dx;  
            end
            
            % Add friction source term if specified
            if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
                H = w(:, 1);
                HU = w(:, 2);
                wet_indices = H > dry_tol;
                
                if any(wet_indices)
                    friction_term = cfg.phys.friction_model(H(wet_indices), HU(wet_indices), g, cfg);
                    dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term(:);
                end
            end
            
            % Combine flux and source terms
            dwdt(nl,:,:) = dwdt_flux + dwdt_source;
        end

        % Flatten output for ODE solver
        dwdt_flat = reshape(dwdt,Nlayers*N*2,1);

    end

    function [wL_hydro, wR_hydro] = hydrostatic_reconstruction(wL, wR, z_interfaces, dry_tol, g)
        % Initialize output with input values
        wL_hydro = wL;
        wR_hydro = wR;

        for i = 1:size(wL, 1)
            % Get elevations at interface
            z_left = z_interfaces(i);
            z_right = z_interfaces(i);

            % Maximum elevation at interface (exactly the same value for left and right)
            z_max = z_left;  % Since z_left = z_right for exact well-balancing

            % Extract values
            h_left = wL(i, 1);
            h_right = wR(i, 1);

            % Calculate velocities carefully
            u_left = 0;
            if h_left > dry_tol
                u_left = wL(i, 2) / h_left;
            end

            u_right = 0;
            if h_right > dry_tol
                u_right = wR(i, 2) / h_right;
            end

            % Hydrostatic reconstruction - formula from the paper
            h_left_recon = max(0, h_left + z_left - z_max);
            h_right_recon = max(0, h_right + z_right - z_max);

            % Update conserved variables
            wL_hydro(i, 1) = h_left_recon;
            wL_hydro(i, 2) = h_left_recon * u_left;

            wR_hydro(i, 1) = h_right_recon;
            wR_hydro(i, 2) = h_right_recon * u_right;
        end
    end

    function z_interfaces = compute_interface_topo(z_cell, cfg, ng, N)
        % Compute topography at interfaces using averaged cell values
        % This ensures z_i-1/2 is exactly the same when viewed from cells i-1 and i
        Nl = (cfg.mesh.Nlayers * 0) + 1/2;
        z_interfaces = zeros(1, N+1);
        z_interfaces(1:N+1) = Nl * (z_cell(ng+(1:N+1)-1) + z_cell(ng+(1:N+1)));
    end

    function x_centers = get_all_cell_centers(cfg, ng)
        % Get coordinates of cell centers including ghost cells
        dx = cfg.mesh.dx;
        x_interior = cfg.mesh.xc;

        % Ghost cells to the left
        x_left_ghost = zeros(1, ng);
        x_left_ghost(1:ng) = x_interior(1) - (ng-(1:ng)+1) * dx;

        % Ghost cells to the right
        x_right_ghost = zeros(1, ng);
        x_right_ghost(1:ng) = x_interior(end) + (1:ng) * dx;

        % Combine all cell centers
        x_centers = [x_left_ghost, x_interior, x_right_ghost];
    end

    function w_padded = apply_boundary_conditions(w_padded, t, cfg)
        % Apply boundary conditions to fill ghost cells
        is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);
        ng = cfg.bc.num_ghost_cells;

        if is_periodic
            w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, ng);
        else
            if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
                w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, ng);
            end
            if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
                w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, ng);
            end
        end
    end

    function [sol_out, t_out, stats] = integrate_rk4_adaptive_local(rhs_func, tspan, w0, cfg)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % +time/integrate_rk4_adaptive.m
        %
        % Purpose:
        %   Solves a system of ODEs dw/dt = rhs_func(t, w, cfg) using the classic
        %   explicit 4th-order Runge-Kutta (RK4) method with adaptive time stepping.
        %   While RK4 itself is fixed-order, the time step `dt` is adapted at each step
        %   based on the CFL condition to ensure numerical stability for hyperbolic problems.
        %   Solution is stored at user-specified output times, and the step size is
        %   adjusted to hit these times exactly.
        %
        % Syntax:
        %   [sol_out, t_out, stats] = integrate_rk4_adaptive(rhs_func, tspan, w0, cfg)
        %
        % Inputs:
        %   rhs_func - [function handle] RHS of the ODE system. Signature:
        %                f = rhs_func(t, w, cfg)
        %   tspan    - [vector, double] Time points [t0, t1, ..., tf] at which the
        %                solution output is requested. Must be monotonically increasing.
        %   w0       - [vector, double] Initial state vector (column vector) at time t0.
        %   cfg      - [struct] Configuration structure. Must contain:
        %                cfg.phys.g, cfg.time.cfl, cfg.mesh.N, cfg.mesh.dx,
        %                cfg.time.num_progress_reports (for progress bar).
        %
        % Outputs:
        %   sol_out  - [M x length(w0), double] Solution matrix. Each row `sol_out(i,:)`
        %                is the state vector corresponding to the time point `t_out(i)`.
        %   t_out    - [1 x M, double] Row vector of output times.
        %   stats    - [struct] Statistics:
        %                stats.nsteps:   Total number of internal RK4 time steps taken.
        %                stats.nfevals:  Total number of RHS evaluations (4 * nsteps for RK4).
        %
        % Dependencies:
        %   - utils.calculate_dt_cfl.m (for adaptive time step)
        %   - Progress bar utility (optional)
        %
        % References:
        %   - Butcher, J. C. (2008). Numerical Methods for Ordinary Differential Equations (2nd ed.). Wiley.
        %   - LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations.
        %     SIAM.
        %
        % Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
        % Date:   21 April 2025
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Input Validation and Setup                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin < 4
            error('integrate_rk4_adaptive:NotEnoughInputs', 'Not enough input arguments.');
        end
        if ~isfield(cfg, 'time') || ~isfield(cfg.time, 'cfl') || isempty(cfg.time.cfl) || cfg.time.cfl <= 0
            error('integrate_rk4_adaptive:MissingCFL', 'CFL number must be specified and positive in cfg.time.cfl');
        end
        if ~isvector(tspan) || ~issorted(tspan) || tspan(1) < 0 || length(tspan) < 2
            error('integrate_rk4_adaptive:InvalidTSPAN', 'TSPAN must be a monotonically increasing vector with at least two elements, starting from t0 >= 0.');
        end

        t0 = tspan(1);         % [s] Initial time
        tf = tspan(end);       % [s] Final time
        t_out_req = tspan(:)'; % Ensure requested output times is a row vector

        w = w0(:); % Ensure w0 is a column vector for internal calculations
        t = t0;    % [s] Current simulation time
        k = 0;     % Step counter
        nfevals = 0; % RHS evaluation counter

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preallocate Output Arrays                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_outputs = length(t_out_req);
        % Store solution as columns initially for easier concatenation if resize needed
        sol_out_internal = zeros(length(w0), num_outputs);
        t_out = zeros(1, num_outputs);
        % Estimate max steps for dt_history (can be resized if needed)
        max_diff_val = max(max(diff(t_out_req), 1e-6)); % Get the single maximum value
        estimated_steps = ceil(10 * (tf - t0) / max_diff_val) + 100; 
        dt_history = zeros(1, estimated_steps);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store Initial Condition                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        output_idx = 1;
        if abs(t - t_out_req(output_idx)) < 1e-12 % Check if t0 is the first output time
            sol_out_internal(:, output_idx) = w;
            t_out(output_idx) = t;
            output_idx = output_idx + 1;
        end
        if num_outputs >= output_idx
             t_next_plot = t_out_req(output_idx);
        else
             t_next_plot = tf + 1; % No more plotting needed
        end

        fprintf('Starting adaptive RK4 integration from t=%.3f to t=%.3f\n', t0, tf);
        fprintf('Output requested at %d time points.\n', num_outputs);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Progress Reporting Setup                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        last_report_time = t0;
        num_reports = 10; % Default number of reports
        if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
            num_reports = cfg.time.num_progress_reports;
        end
        report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times

        % --- aux vars
        t_plot = 0;
        ufind  = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main Time Stepping Loop                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        max_internal_steps = 1e7; % Safety break
        while t < tf
            if k >= max_internal_steps
                 warning('integrate_rk4_adaptive:MaxStepsExceeded', 'Maximum internal steps (%d) exceeded. Aborting.', max_internal_steps);
                 break;
            end

            % --- Calculate Adaptive Time Step Based on CFL ---
            dt = calculate_dt_cfl(w, cfg);

            % --- Adjust dt to Hit Output Times Exactly ---
            dt_to_tf = tf - t;
            dt_to_plot = t_next_plot - t;
            % Choose smallest of CFL dt, time to tf, time to next plot
            dt = min([dt, dt_to_tf, dt_to_plot]);

            % --- Safety Checks for dt ---
            if dt <= 1e-12 % Prevent excessively small steps
                if abs(t-tf) < 1e-9
                     fprintf('Reached final time tf=%.4f\n', tf);
                     break; % Exit loop if effectively at the end time
                else
                    warning('integrate_rk4_adaptive:SmallDt', 'Time step dt=%.3e is too small at t=%.3f. Aborting integration.', dt, t);
                    break;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Classic 4th-Order Runge-Kutta (RK4) Stages              %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Formula: w_{n+1} = w_n + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
            % where k_i are estimates of the slope at different points in the interval.
            k1 = rhs_func(t,             w,            cfg); % Slope at the beginning
            k2 = rhs_func(t + 0.5 * dt, w + 0.5 * dt*k1, cfg); % Slope at midpoint using k1
            k3 = rhs_func(t + 0.5 * dt, w + 0.5 * dt*k2, cfg); % Slope at midpoint using k2
            k4 = rhs_func(t + dt,       w + dt*k3,       cfg); % Slope at the end using k3
            nfevals = nfevals + 4; % Increment RHS evaluation count

            % --- Update Solution and Time ---
            w_new = w + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            t_new = t + dt;
            k = k + 1;
            
            % --- update plot & save time ---
            t_plot = t_plot + dt;

            % --- Store dt History (Resize if needed) ---
            if k > length(dt_history)
                 dt_chunk_size = length(dt_history); % Double the current size
                 warning('Time:Integrate:GrowStats', 'Growing dt_history size at step %d (t=%.3f)', k, t);
                 dt_history(end+1 : end+dt_chunk_size) = 0; % Grow using direct indexing
            end
            dt_history(k) = dt;

            % --- Progress Reporting ---
            if report_interval > 0 && t_new >= last_report_time + report_interval
                fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t_new, 100*(t_new-t0)/(tf-t0), dt);
                last_report_time = t_new;
            end

            % --- Output Handling: Store Solution at Requested Times ---
            % Check if the new time step landed exactly on a requested output time.
            if abs(t_new - t_next_plot) < 1e-12
                sol_out_internal(:, output_idx) = w_new;
                t_out(output_idx) = t_new;
                output_idx = output_idx + 1;
                if output_idx <= num_outputs
                     t_next_plot = t_out_req(output_idx);
                else
                     t_next_plot = tf + 1; % No more outputs needed
                end
            end

            % --- Prepare for Next Step ---
            w = w_new;
            t = t_new;

            if t_plot >= cfg.vis.dt_plot
                t_plot = 0;
                
                loc_plot = 0;
                
                Nplots_tot = 2;
                
                w_plt = reshape(w_new,cfg.mesh.Nlayers,cfg.mesh.N,2);                                
                
                z_b = cfg.bathyHandle(cfg, cfg.mesh.xc);
                
                loc_plot = loc_plot + 1;
                subplot(Nplots_tot,1,loc_plot)
                    for nl = 1:cfg.mesh.Nlayers
                        h_plt = sum(w_plt(1:nl,:,1), 1) + z_b;
                        plot(cfg.mesh.xc,h_plt)
                        hold on;
                    end
                    plot(cfg.mesh.xc,z_b,'k--')
                    hold off;
                    title(['Time:', blanks(1), num2str(t,'%6.3f')]);
                    ylim([0 cfg.atmosphere.Altitudes(cfg.mesh.Nlayers) + 2*max(cfg.atmosphere.dH)]);
                    drawnow;
                    
                loc_plot = loc_plot + 1;                    
                subplot(Nplots_tot,1,loc_plot)
                    u_plt = squeeze(w_plt(:,:,2)./w_plt(:,:,1));
                    for nl = 1:cfg.mesh.Nlayers                        
                        plot(cfg.mesh.xc,u_plt)
                        hold on;
                    end
                    
                    if ufind == 0
                        umax = 1.2*max(max(abs(u_plt)));
                        ufind = 0;
                    end
                    
                    hold off;
                    ylim([-umax umax]);
                    drawnow;
                    
                    hold off;
                    drawnow;                    
                    
            end
            
        end % End while loop

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Final Output Formatting and Statistics                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Trim unused preallocated space
        sol_out_internal = sol_out_internal(:, 1:output_idx-1);
        t_out = t_out(1:output_idx-1);
        dt_history = dt_history(1:k);

        % Transpose solution to match expected output format [M x length(w0)]
        sol_out = sol_out_internal';

        % --- Statistics ---
        stats.nsteps = k;
        stats.nfevals = nfevals;
        % stats.dt_history = dt_history; % Optionally return dt history

        fprintf('Integration finished at t = %.3f s after %d steps.\n', t_out(end), k);

    end % Function end

    function dt = calculate_dt_cfl(w, cfg)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % +core/+utils/calculate_dt_cfl.m
        %
        % Purpose:
        %   Computes the maximum stable time step for the 1D Non-Linear Shallow Water
        %   equations using the Courant-Friedrichs-Lewy (CFL) condition. Ensures
        %   numerical stability of explicit time-stepping schemes by limiting the
        %   time step according to the fastest wave speed in the domain.
        %
        % Syntax:
        %   dt = calculate_dt_cfl(w, cfg)
        %
        % Inputs:
        %   w   - [2N x 1, double] State vector at the current time, containing water depth (H)
        %         and discharge (HU) for all N cells: w = [H; HU].
        %   cfg - [struct] Configuration structure. Required fields:
        %         cfg.mesh.N:   [integer] Number of cells.
        %         cfg.time.cfl: [double] Desired CFL number (typically <= 1).
        %         cfg.mesh.dx:  [double] Spatial cell width [m].
        %         cfg.phys.g:   [double] Acceleration due to gravity [m/s^2].
        %
        % Outputs:
        %   dt  - [double] Calculated maximum stable time step according to the CFL condition [s].
        %
        % Dependencies:
        %   None (utility function, but expects correct cfg structure).
        %
        % References:
        %   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
        %   - Courant, R., Friedrichs, K., & Lewy, H. (1928). "On the Partial Difference Equations of Mathematical Physics."
        %
        % Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
        % Date:   21 April 2025
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Input Checks and Parameter Extraction                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isfield(cfg.time, 'cfl') || cfg.time.cfl <= 0
            error('cfg.time.cfl must be a positive value.');
        end
        if ~isfield(cfg.mesh, 'dx') || cfg.mesh.dx <= 0
            error('cfg.mesh.dx must be a positive value.');
        end
        if ~isfield(cfg.phys, 'g') || cfg.phys.g <= 0
            error('cfg.phys.g must be a positive value.');
        end
        if ~isfield(cfg.mesh, 'N') || cfg.mesh.N <= 0
            error('cfg.mesh.N must be a positive integer.');
        end

        N = cfg.mesh.N;
        Nlayers = cfg.mesh.Nlayers;
        
        if length(w) ~= 2*N*Nlayers
            error('Input state vector w must have length 2*N. Got length(w) = %d, N = %d, size(w) = [%d %d]', length(w), N, size(w,1), size(w,2));
        end

        g = cfg.phys.g;     % [m/s^2] Acceleration due to gravity
        dx = cfg.mesh.dx;   % [m] Spatial step size
        cfl = cfg.time.cfl; % [unitless] CFL number

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract H and HU from state vector w                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        w = reshape(w,Nlayers,N,2);
        
        H  = squeeze(w(Nlayers,N,1));         % [m] Water depth (first N elements)
        HU = squeeze(w(Nlayers,N,2));         % [m^2/s] Discharge (next N elements)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Time Step by CFL Condition                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The CFL condition for explicit FV schemes is:
        %   dt <= cfl * dx / max(|U| + c)
        % where U = HU/H (velocity), c = sqrt(g*H) (wave speed)
        %
        % To avoid instability in dry or nearly dry cells, a minimum threshold
        % is enforced for H. This prevents division by zero and ensures that
        % the wave speed is always well-defined.
        H_min = 1e-6; % [m] Minimum water depth threshold for stable velocity calculation
        H_floor = max(H, H_min); % Ensure H is above the minimum threshold

        % Calculate velocity U = HU / H, handling potential division by zero
        U = zeros(size(H)); % [m/s] Velocity vector
        wet_indices = H > H_min; % Find indices of wet cells
        U(wet_indices) = HU(wet_indices) ./ H(wet_indices);

        % Calculate the wave speed c = sqrt(g*H)
        wave_speed = sqrt(g * H_floor); % [m/s]

        % The maximum signal speed in each cell is |U| + c
        max_speed = abs(U) + wave_speed; % [m/s]
        max_signal_speed = max(max_speed); % [m/s] Maximum over all cells

        % If all cells are dry (H ~ 0), set dt to a large value (simulation will stop)
        if max_signal_speed < 1e-8
            dt = 1e6; % [s] Effectively disables time stepping (simulation should halt)
            return;
        end

        % CFL time step restriction
        dt = cfl * dx / max_signal_speed; % [s]

    end

    function b = gaussian_bump(cfg, x)
        
        %GAUSSIAN_BUMP Defines a Gaussian bump bathymetry.
        %   b = gaussian_bump(cfg, x) returns the bathymetry profile b(x)
        %   representing a flat bottom with a Gaussian bump centered in the domain.
        %   The bathymetry 'b' represents the bottom elevation relative to z=0.
        %   Still water level is at z = cfg.h0. Water depth is eta - b.
        %
        %   Inputs:
        %       cfg - Configuration structure containing parameters. Expected fields:
        %             cfg.h0                - Still water depth (required).
        %             cfg.L                 - Domain length (required).
        %             cfg.bathy_bump_center - Center position of the bump (optional, default: L/2).
        %             cfg.bathy_bump_height - Height of the bump (optional, default: 0.2 * h0).
        %             cfg.bathy_bump_width  - Characteristic width (std dev) of the bump (optional, default: L/10).
        %             cfg.min_depth         - Minimum allowed water depth (optional, used for warning/capping).
        %       x   - Vector of spatial coordinates.
        %
        %   Output:
        %       b   - Vector of bathymetry elevation values at coordinates x.
        %
        %   Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
        %   Date:   2025-05-12
        
        % --- Input Validation & Default Parameters ---
        % Robustly extract h0 and L (accepts config.h0 or config.param.H0, config.L or config.domain.xmax-xmin)
        if isfield(cfg, 'h0') && ~isempty(cfg.h0)
            h0 = cfg.h0;
        elseif isfield(cfg, 'param') && isfield(cfg.param, 'H0') && ~isempty(cfg.param.H0)
            h0 = cfg.param.H0;
        else
            error('gaussian_bump:MissingParameter', 'Still water depth cfg.h0 (or cfg.param.H0) must be provided.');
        end
        if isfield(cfg, 'L') && ~isempty(cfg.L)
            L = cfg.L;
        elseif isfield(cfg, 'domain') && isfield(cfg.domain, 'xmax') && isfield(cfg.domain, 'xmin')
            L = cfg.domain.xmax - cfg.domain.xmin;
        else
            error('gaussian_bump:MissingParameter', 'Domain length cfg.L (or cfg.domain.xmax-xmin) must be provided.');
        end
        
        % Use these robust values throughout
        
        % Default parameters
        default_center = L / 2;
        default_height = 0.2 * h0;
        default_width  = L / 10;
        default_min_depth = 1e-6; % Default minimum depth to avoid division by zero etc.
        
        bump_center = default_center;
        if isfield(cfg, 'bathy_bump_center') && ~isempty(cfg.bathy_bump_center)
            bump_center = cfg.bathy_bump_center;
        end
        
        bump_height = default_height;
        if isfield(cfg, 'bathy_bump_height') && ~isempty(cfg.bathy_bump_height)
            bump_height = cfg.bathy_bump_height;
            if bump_height < 0
                warning('gaussian_bump:InvalidHeight', 'Specified bump height is negative. Using absolute value.');
                bump_height = abs(bump_height);
            end
        end
        
        bump_width = default_width;
        if isfield(cfg, 'bathy_bump_width') && ~isempty(cfg.bathy_bump_width)
            bump_width = cfg.bathy_bump_width;
            if bump_width <= 0
                warning('gaussian_bump:InvalidWidth', 'Specified bump width must be positive. Using default.');
                bump_width = default_width;
            end
        end
        
        min_depth = default_min_depth;
        if isfield(cfg, 'min_depth') && ~isempty(cfg.min_depth)
            min_depth = cfg.min_depth;
        end
        % --- End Validation ---
        
        % --- Calculate Bathymetry ---
        % Gaussian bump centered at 'bump_center' with height 'bump_height'
        % and standard deviation 'bump_width'.
        % b(x) = bump_height * exp(-(x - bump_center)^2 / (2 * bump_width^2))
        b = bump_height * exp(-((x - bump_center).^2) / (2 * bump_width^2));
        % --- End Calculation ---
        
        % --- Check for excessive bump height ---
        % Water depth at rest over the bump is h0 - b. Ensure it's > min_depth.
        max_b = max(b);
        if max_b >= cfg.h0 - min_depth
            warning('gaussian_bump:PotentialDryArea', ...
                ['Maximum bump height (%.2f) results in water depth <= min_depth (%.2e) ',...
                'at some points for still water level h0=%.2f. Simulation might become unstable.'], ...
                max_b, min_depth, cfg.h0);
            % Optional: Cap the bump height to ensure minimum depth
            % scaling_factor = (cfg.h0 - min_depth) / max_b;
            % b = b * scaling_factor;
            % fprintf('Bump height capped to %.2f to maintain minimum depth.\n', max(b));
        end
        % --- End Check ---
        
    end

end















