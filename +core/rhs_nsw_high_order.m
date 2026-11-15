%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_high_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations by the second-order well-balanced scheme of
%   core.rhs_nsw_2nd_order, under a second name that configurations may
%   use. Under this name the friction coefficient cfg.phys.Cf is not
%   applied, friction coming from cfg.phys.friction_model alone, and the
%   results are those of core.rhs_nsw_2nd_order(t, w_flat, cfg, false), to
%   the last bit.
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

function dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)
    dwdt_flat = core.rhs_nsw_2nd_order(t, w_flat, cfg, false);
end
