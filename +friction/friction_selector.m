%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/friction_selector.m
%
% Purpose:
%   Helper function that returns the appropriate friction model function
%   handle based on the configuration. Serves as a central registry for
%   all available friction models.
%
% Syntax:
%   friction_model = friction_selector(friction_name)
%
% Inputs:
%   friction_name - [char] String identifier for the friction model:
%                   'none'   -> No friction
%                   'chezy'  -> Chézy model
%                   [future models will be added here]
%
% Outputs:
%   friction_model - [function_handle] Function handle to the selected friction model.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function friction_model = friction_selector(friction_name)

    % Default: no friction
    if nargin < 1 || isempty(friction_name)
        friction_name = 'none';
    end
    
    % Convert to lowercase for case-insensitive comparison
    friction_name = lower(friction_name);
    
    % Select appropriate friction model
    switch friction_name
        case {'none', 'no_friction', 'nofriction', 'off', '0'}
            friction_model = @friction.no_friction;
            
        case {'chezy', 'chézy', 'chezi'}
            friction_model = @friction.chezy;
            
        case {'manning', 'mannings', 'manning_n'}
            friction_model = @friction.manning;
            
        case {'darcy', 'darcy_weisbach', 'darcyweisbach'}
            friction_model = @friction.darcy_weisbach;
            
        % Note: Colebrook-White is not a standalone friction model but used with Darcy-Weisbach
            
        otherwise
            warning('Unknown friction model "%s". Using default (no friction).', friction_name);
            friction_model = @friction.no_friction;
    end

end