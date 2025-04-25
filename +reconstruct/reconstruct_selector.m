function reconstruct_handle = reconstruct_selector(method_name)
% RECONSTRUCT_SELECTOR Function to select a reconstruction method
%
% Purpose:
%   Returns the function handle for the specified reconstruction method.
%   This provides a central point for selecting and configuring
%   reconstruction methods.
%
% Syntax:
%   reconstruct_handle = reconstruct_selector(method_name)
%
% Inputs:
%   method_name - [char] String identifier for the reconstruction method:
%                 'none'   -> No reconstruction (1st order)
%                 'muscl'  -> MUSCL scheme (2nd order)
%                 'eno2'   -> ENO (2nd order)
%                 'uno2'   -> UNO (2nd order) 
%                 'weno5'  -> WENO (5th order)
%
% Outputs:
%   reconstruct_handle - [function_handle] Function handle to the selected method
%
% Dependencies:
%   The reconstruction method functions in the +reconstruct package
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Default: no reconstruction (1st order)
if nargin < 1 || isempty(method_name)
    method_name = 'none';
end

% Convert to lowercase for case-insensitive comparison
method_name = lower(method_name);

% Select appropriate reconstruction method
switch method_name
    case {'none', 'first_order', '1st_order', 'off', '0'}
        reconstruct_handle = @reconstruct.none;
        
    case {'muscl', 'muscl2', 'second_order', '2nd_order'}
        reconstruct_handle = @reconstruct.muscl;
        
    case {'eno', 'eno2'}
        reconstruct_handle = @reconstruct.eno2;
        
    case {'uno', 'uno2'}
        reconstruct_handle = @reconstruct.uno2;
        
    case {'weno', 'weno5'}
        reconstruct_handle = @reconstruct.weno5;
        
    case {'ppm'}
        reconstruct_handle = @reconstruct.ppm;
        
    otherwise
        warning('Unknown reconstruction method "%s". Using default (none/1st order).', method_name);
        reconstruct_handle = @reconstruct.none;
end
end
