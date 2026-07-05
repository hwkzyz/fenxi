%% Run_Main_GradientXRange030_20251222
% One-command main pipeline for the 20251222 low-speed rotating template
% route with noise-aware gradient x-domain ratio = 0.30.
% Keep the original validated DynamicMap naming/center convention. The
% GradientXRange030 template is used in Step03, while Step02 writes the
% base-named DynamicMap selected by the 50.2 s start-time scan.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
run(fullfile(route_dir, 'Step01_Main_Build_GradientXRange030_Template_20251222.m'));
setenv('STEP02_TEMPLATE_SUFFIX', '');
setenv('STEP02_DYNAMIC_SUFFIX', '');
run(fullfile(route_dir, 'Step02_Build_Dynamic_Map_20251222.m'));
run(fullfile(route_dir, 'Step03_Main_Run_GradientXRange030_Identification_20251222.m'));
run(fullfile(route_dir, 'Step04_Main_Visualize_GradientXRange030_20251222.m'));
