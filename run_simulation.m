% run_simulation.m

clear; close all; format longE;
addpath(genpath(pwd)); % Add all subfolders

fprintf('Loading configuration...\n');
cfg = cfg.flat_bottom_config(); % Load specific experiment setup

fprintf('Starting simulation: %s\n', cfg.caseName);
results = core.solver(cfg); % Call the main solver function
fprintf('Simulation finished.\n');