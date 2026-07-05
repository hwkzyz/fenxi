%% Run_Main_OPRCenterStd_20250527
% One-command OPRCenterStd direct-template pipeline for the 20250527 case.
% Figures are intentionally kept open for visual inspection.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
run(fullfile(routeDir, 'Step01_Main_Build_OPRCenterStd_Template_20250527.m'));
run(fullfile(routeDir, 'Step02_Main_Build_OPRCenterStd_DynamicMap_20250527.m'));
run(fullfile(routeDir, 'Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20250527.m'));
run(fullfile(routeDir, 'Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20250527.m'));
