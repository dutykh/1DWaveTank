function cfg = validate_config(cfg)
%VALIDATE_CONFIG Centralized validation for the simulation configuration structure.
%
% Syntax:
%   cfg_validated = validate_config(cfg_raw)
%
% Description:
%   This function takes a raw configuration structure (cfg_raw), validates
%   its fields against a defined schema, checks types and ranges, sets
%   default values for optional fields, and returns the validated and
%   potentially completed configuration structure (cfg_validated).
%
%   It calls individual validation functions for each major subsection
%   of the configuration (mesh, time, physics, bc, numerics, ic, output).
%
% Inputs:
%   cfg - [struct] The raw configuration structure loaded from
%          simulation_config.m and potentially merged with defaults.
%
% Outputs:
%   cfg - [struct] The validated and completed configuration structure.
%          An error is thrown if validation fails.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Assisted by: Cascade AI Assistant
% Date:   23 April 2025

    fprintf('--- Validating Configuration Structure ---\n');

    if ~isstruct(cfg)
        error('ConfigValidation:InvalidInput', 'Input configuration must be a structure.');
    end

    % Define the order of validation
    validation_functions = {
        @validate_mesh_config, ...
        @validate_physics_config, ...
        @validate_time_config, ...
        @validate_numerics_config, ...
        @validate_bc_config, ...
        @validate_ic_config, ...
        @validate_output_config ...
        % Add other validation functions here if needed
    };

    % Apply each validation function sequentially
    for i = 1:length(validation_functions)
        func = validation_functions{i};
        fprintf('  Validating section: %s\n', func2str(func));
        cfg = func(cfg);
    end

    fprintf('--- Configuration Validation Complete ---\n');

end

% =========================================================================
%                      Subsection Validation Functions
% =========================================================================

function cfg = validate_mesh_config(cfg)
    % Validates the [cfg.mesh] subsection
    if ~isfield(cfg, 'mesh'), cfg.mesh = struct(); end

    % --- Required Fields ---
    check_required_field(cfg.mesh, 'domain', 'mesh');
    check_required_field(cfg.mesh.domain, 'xmin', 'mesh.domain');
    check_required_field(cfg.mesh.domain, 'xmax', 'mesh.domain');
    check_required_field(cfg.mesh, 'N', 'mesh');

    % --- Type and Attribute Validation ---
    validateattributes(cfg.mesh.domain.xmin, {'numeric'}, {'scalar', 'finite'}, mfilename, 'cfg.mesh.domain.xmin');
    validateattributes(cfg.mesh.domain.xmax, {'numeric'}, {'scalar', 'finite', '>', cfg.mesh.domain.xmin}, mfilename, 'cfg.mesh.domain.xmax');
    validateattributes(cfg.mesh.N, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'cfg.mesh.N');

    % --- Optional Fields & Defaults ---
    % Calculate dx if not provided
    if ~isfield(cfg.mesh, 'dx') || isempty(cfg.mesh.dx)
        cfg.mesh.dx = (cfg.mesh.domain.xmax - cfg.mesh.domain.xmin) / cfg.mesh.N;
        fprintf('    Setting default: cfg.mesh.dx = %.4e\n', cfg.mesh.dx);
    else
        validateattributes(cfg.mesh.dx, {'numeric'}, {'scalar', 'positive', 'finite'}, mfilename, 'cfg.mesh.dx');
        % Optional: check consistency
        calculated_dx = (cfg.mesh.domain.xmax - cfg.mesh.domain.xmin) / cfg.mesh.N;
        if abs(cfg.mesh.dx - calculated_dx) > 1e-12 * abs(calculated_dx)
             warning('ConfigValidation:InconsistentDx', ...
                'Provided cfg.mesh.dx (%.4e) does not match value calculated from domain and N (%.4e). Using provided value.', ...
                cfg.mesh.dx, calculated_dx);
        end
    end

    % Calculate cell centers xc if not provided
    if ~isfield(cfg.mesh, 'xc') || isempty(cfg.mesh.xc)
        cfg.mesh.xc = linspace(cfg.mesh.domain.xmin + cfg.mesh.dx/2, ...
                               cfg.mesh.domain.xmax - cfg.mesh.dx/2, ...
                               cfg.mesh.N); % Row vector by default
        fprintf('    Setting default: cfg.mesh.xc (vector 1 x N)\n');
    else
        validateattributes(cfg.mesh.xc, {'numeric'}, {'vector', 'numel', cfg.mesh.N, 'finite'}, mfilename, 'cfg.mesh.xc');
        cfg.mesh.xc = cfg.mesh.xc(:).'; % Ensure row vector
    end

    % Calculate cell interfaces xf if not provided
    if ~isfield(cfg.mesh, 'xf') || isempty(cfg.mesh.xf)
        cfg.mesh.xf = linspace(cfg.mesh.domain.xmin, cfg.mesh.domain.xmax, cfg.mesh.N + 1)'; % Column vector
        fprintf('    Setting default: cfg.mesh.xf (vector (N+1) x 1)\n');
    else
        validateattributes(cfg.mesh.xf, {'numeric'}, {'vector', 'numel', cfg.mesh.N + 1, 'finite'}, mfilename, 'cfg.mesh.xf');
        cfg.mesh.xf = cfg.mesh.xf(:); % Ensure column vector
    end
end

% -------------------------------------------------------------------------

function cfg = validate_physics_config(cfg)
    % Validates the [cfg.phys] subsection
    if ~isfield(cfg, 'phys'), cfg.phys = struct(); end

    % --- Required Fields ---
    check_required_field(cfg.phys, 'g', 'phys');

    % --- Type and Attribute Validation ---
    validateattributes(cfg.phys.g, {'numeric'}, {'scalar', 'positive', 'finite'}, mfilename, 'cfg.phys.g');

    % --- Optional Fields & Defaults ---
    if ~isfield(cfg.phys, 'Cf') || isempty(cfg.phys.Cf)
        cfg.phys.Cf = 0.0; % Default: no friction
        fprintf('    Setting default: cfg.phys.Cf = %.4e\n', cfg.phys.Cf);
    else
        validateattributes(cfg.phys.Cf, {'numeric'}, {'scalar', 'nonnegative', 'finite'}, mfilename, 'cfg.phys.Cf');
    end
end

% -------------------------------------------------------------------------

function cfg = validate_time_config(cfg)
    % Validates the [cfg.time] subsection
    if ~isfield(cfg, 'time'), cfg.time = struct(); end

    % --- Required Fields ---
    check_required_field(cfg, 't0', 'global'); % Often set globally
    check_required_field(cfg, 'tEnd', 'global'); % Often set globally

    % --- Type and Attribute Validation ---
    validateattributes(cfg.t0, {'numeric'}, {'scalar', 'finite'}, mfilename, 'cfg.t0');
    validateattributes(cfg.tEnd, {'numeric'}, {'scalar', 'finite', '>', cfg.t0}, mfilename, 'cfg.tEnd');

    % --- Optional Fields & Defaults ---
    if ~isfield(cfg.time, 'cfl') || isempty(cfg.time.cfl)
        cfg.time.cfl = 0.95; % Default CFL number
        fprintf('    Setting default: cfg.time.cfl = %.2f\n', cfg.time.cfl);
    else
        validateattributes(cfg.time.cfl, {'numeric'}, {'scalar', 'positive', '<=', 1.0}, mfilename, 'cfg.time.cfl');
    end

    if ~isfield(cfg.time, 'dt_stable') || isempty(cfg.time.dt_stable)
        cfg.time.dt_stable = []; % Allow calculation by adaptive stepper
        fprintf('    Setting default: cfg.time.dt_stable = [] (will be calculated)\n');
    else
        validateattributes(cfg.time.dt_stable, {'numeric'}, {'scalar', 'positive', 'finite'}, mfilename, 'cfg.time.dt_stable');
    end

    % --- Output time vector --- 
    if ~isfield(cfg, 'vis') || ~isfield(cfg.vis, 'dt_plot') || isempty(cfg.vis.dt_plot)
        cfg.vis.dt_plot = (cfg.tEnd - cfg.t0) / 10; % Default: 10 output steps + initial
        fprintf('    Setting default: cfg.vis.dt_plot = %.4f\n', cfg.vis.dt_plot);
    else 
        validateattributes(cfg.vis.dt_plot, {'numeric'}, {'scalar', 'positive', 'finite'}, mfilename, 'cfg.vis.dt_plot');
    end
    
    if ~isfield(cfg.time, 't_out') || isempty(cfg.time.t_out)
        t_out = cfg.t0:cfg.vis.dt_plot:cfg.tEnd;
        if t_out(end) < cfg.tEnd - 1e-9 % Ensure tEnd is included
            t_out = [t_out, cfg.tEnd];
        end
        cfg.time.t_out = unique(t_out); % Ensure unique, sorted times
        fprintf('    Setting default: cfg.time.t_out (vector %d x 1)\n', numel(cfg.time.t_out));
    else
        validateattributes(cfg.time.t_out, {'numeric'}, {'vector', 'nondecreasing', 'finite'}, mfilename, 'cfg.time.t_out');
        cfg.time.t_out = unique(cfg.time.t_out(:)); % Ensure unique, sorted column vector
        if cfg.time.t_out(1) < cfg.t0 || cfg.time.t_out(end) > cfg.tEnd
            warning('ConfigValidation:InvalidTimeOut', 'cfg.time.t_out extends beyond [t0, tEnd]. Clamping.');
            cfg.time.t_out = cfg.time.t_out(cfg.time.t_out >= cfg.t0 & cfg.time.t_out <= cfg.tEnd);
        end
    end
    
    if ~isfield(cfg.time, 'num_progress_reports') || isempty(cfg.time.num_progress_reports)
        cfg.time.num_progress_reports = 10; % Default: 10 progress reports
        fprintf('    Setting default: cfg.time.num_progress_reports = %d\n', cfg.time.num_progress_reports);
    else
        validateattributes(cfg.time.num_progress_reports, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'cfg.time.num_progress_reports');
    end
end

% -------------------------------------------------------------------------

function cfg = validate_numerics_config(cfg)
    % Validates the numerical scheme settings (flux, stepper, etc.)
    if ~isfield(cfg, 'numerics'), cfg.numerics = struct(); end % Typically fields are global

    % --- Required Fields ---
    check_required_field(cfg, 'numFlux', 'global');
    check_required_field(cfg, 'timeStepper', 'global');
    check_required_field(cfg, 'model', 'global'); % The RHS function

    % --- Type and Attribute Validation ---
    validateattributes(cfg.numFlux, {'function_handle'}, {}, mfilename, 'cfg.numFlux');
    validateattributes(cfg.timeStepper, {'function_handle'}, {}, mfilename, 'cfg.timeStepper');
    validateattributes(cfg.model, {'function_handle'}, {}, mfilename, 'cfg.model');

    % --- Optional Fields & Defaults ---
    if ~isfield(cfg, 'reconstruction') || isempty(cfg.reconstruction)
        cfg.reconstruction = []; % Default: No reconstruction (1st order)
        fprintf('    Setting default: cfg.reconstruction = [] (1st order scheme)\n');
    else
        % Add validation if reconstruction options are implemented
        validateattributes(cfg.reconstruction, {'struct', 'function_handle'}, {}, mfilename, 'cfg.reconstruction');
    end
end

% -------------------------------------------------------------------------

function cfg = validate_bc_config(cfg)
    % Validates the [cfg.bc] subsection
    if ~isfield(cfg, 'bc'), cfg.bc = struct(); end

    % --- Required Fields ---
    check_required_field(cfg.bc, 'left', 'bc');
    check_required_field(cfg.bc, 'right', 'bc');
    check_required_field(cfg.bc.left, 'handle', 'bc.left');
    check_required_field(cfg.bc.right, 'handle', 'bc.right');

    % --- Type and Attribute Validation ---
    validateattributes(cfg.bc.left.handle, {'function_handle'}, {}, mfilename, 'cfg.bc.left.handle');
    validateattributes(cfg.bc.right.handle, {'function_handle'}, {}, mfilename, 'cfg.bc.right.handle');

    % --- Optional Fields & Defaults ---
    % Ensure 'param' field exists for BCs, even if empty
    if ~isfield(cfg.bc.left, 'param') || isempty(cfg.bc.left.param)
        cfg.bc.left.param = struct();
        fprintf('    Setting default: cfg.bc.left.param = struct()\n');
    else
        validateattributes(cfg.bc.left.param, {'struct'}, {}, mfilename, 'cfg.bc.left.param');
    end
    if ~isfield(cfg.bc.right, 'param') || isempty(cfg.bc.right.param)
        cfg.bc.right.param = struct();
        fprintf('    Setting default: cfg.bc.right.param = struct()\n');
    else
        validateattributes(cfg.bc.right.param, {'struct'}, {}, mfilename, 'cfg.bc.right.param');
    end

    % Default number of ghost cells (specific to implementation, assume 2 for now)
    if ~isfield(cfg.bc, 'num_ghost_cells') || isempty(cfg.bc.num_ghost_cells)
        cfg.bc.num_ghost_cells = 2; % Common choice for 2nd order schemes
        fprintf('    Setting default: cfg.bc.num_ghost_cells = %d\n', cfg.bc.num_ghost_cells);
    else
        validateattributes(cfg.bc.num_ghost_cells, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'cfg.bc.num_ghost_cells');
    end
end

% -------------------------------------------------------------------------

function cfg = validate_ic_config(cfg)
    % Validates the initial condition settings

    % --- Required Fields ---
    check_required_field(cfg, 'ic_handle', 'global');

    % --- Type and Attribute Validation ---
    validateattributes(cfg.ic_handle, {'function_handle'}, {}, mfilename, 'cfg.ic_handle');

    % --- Optional Fields & Defaults ---
    % Ensure 'ic_param' field exists, even if empty
    if ~isfield(cfg, 'ic_param') || isempty(cfg.ic_param)
        cfg.ic_param = struct();
        fprintf('    Setting default: cfg.ic_param = struct()\n');
    else
        validateattributes(cfg.ic_param, {'struct'}, {}, mfilename, 'cfg.ic_param');
    end

    % Optional: Validate bathymetry handle (often tied to IC)
    check_required_field(cfg, 'bathyHandle', 'global');
    validateattributes(cfg.bathyHandle, {'function_handle'}, {}, mfilename, 'cfg.bathyHandle');
end

% -------------------------------------------------------------------------

function cfg = validate_output_config(cfg)
    % Validates the [cfg.vis] and output settings
    if ~isfield(cfg, 'vis'), cfg.vis = struct(); end

    % --- Required Fields ---
    % (dt_plot is handled in time config)

    % --- Optional Fields & Defaults ---
    if ~isfield(cfg.vis, 'plot_velocity') || isempty(cfg.vis.plot_velocity)
        cfg.vis.plot_velocity = true; % Default: show velocity subplot
        fprintf('    Setting default: cfg.vis.plot_velocity = %s\n', mat2str(cfg.vis.plot_velocity));
    else
        validateattributes(cfg.vis.plot_velocity, {'logical', 'numeric'}, {'scalar', 'binary'}, mfilename, 'cfg.vis.plot_velocity');
    end
    
    if ~isfield(cfg.vis, 'show_legend') || isempty(cfg.vis.show_legend)
        cfg.vis.show_legend = false; % Default: hide legend
        fprintf('    Setting default: cfg.vis.show_legend = %s\n', mat2str(cfg.vis.show_legend));
    else
        validateattributes(cfg.vis.show_legend, {'logical', 'numeric'}, {'scalar', 'binary'}, mfilename, 'cfg.vis.show_legend');
    end

    if ~isfield(cfg, 'save_results') || isempty(cfg.save_results)
        cfg.save_results = false; % Default: don't save
        fprintf('    Setting default: cfg.save_results = %s\n', mat2str(cfg.save_results));
    else
        validateattributes(cfg.save_results, {'logical', 'numeric'}, {'scalar', 'binary'}, mfilename, 'cfg.save_results');
    end

    if cfg.save_results
        if ~isfield(cfg, 'outputPath') || isempty(cfg.outputPath)
            cfg.outputPath = './output/'; % Default output path
            fprintf('    Setting default: cfg.outputPath = ''%s''\n', cfg.outputPath);
        else
            validateattributes(cfg.outputPath, {'char', 'string'}, {'scalartext'}, mfilename, 'cfg.outputPath');
        end
        if ~isfield(cfg, 'caseName') || isempty(cfg.caseName)
             warning('ConfigValidation:MissingCaseName', 'cfg.save_results is true, but cfg.caseName is missing. Using default name.');
             cfg.caseName = sprintf('simulation_results_%s', datestr(now,'yyyymmdd_HHMMSS'));
        else
             validateattributes(cfg.caseName, {'char', 'string'}, {'scalartext'}, mfilename, 'cfg.caseName');
        end
    end
end

% =========================================================================
%                        Helper Functions
% =========================================================================

function check_required_field(structure, fieldname, struct_name)
    % Checks if a required field exists in a structure.
    if ~isfield(structure, fieldname) || isempty(structure.(fieldname))
        error('ConfigValidation:MissingField', ...
              'Required configuration field ''%s.%s'' is missing or empty.', ...
              struct_name, fieldname);
    end
end

% You can add more specific helper functions here if needed, e.g.,
% function check_function_handle(handle, expected_package, field_name)
%     % ... check if handle is a function handle and optionally if it belongs
%     % to an expected package (+bc, +ic, etc.) ...
% end
