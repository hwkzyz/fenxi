%% Run_Main_GradientXRange030_20250527
% One-command main pipeline for the 20250527 low-speed rotating template
% route with noise-aware gradient x-domain ratio = 0.30.
% Keep the original dynamic x-center convention. The GradientXRange030
% template is used for waveform fitting, while Step02 keeps the base
% DynamicMap naming/center used by the validated 20250527 route.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
run(fullfile(route_dir, 'Step01_Main_Build_GradientXRange030_Template_20250527.m'));
setenv('STEP02_TEMPLATE_SUFFIX', '');
setenv('STEP02_DYNAMIC_SUFFIX', '');
run(fullfile(route_dir, 'Step02_Build_Dynamic_Map_20250527.m'));
run(fullfile(route_dir, 'Step03_Main_Run_GradientXRange030_Identification_20250527.m'));
run(fullfile(route_dir, 'Step04_Main_Visualize_GradientXRange030_20250527.m'));
