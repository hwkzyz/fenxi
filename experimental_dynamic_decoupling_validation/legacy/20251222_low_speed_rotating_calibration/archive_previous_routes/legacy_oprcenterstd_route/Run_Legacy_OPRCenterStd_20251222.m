%% Run_Legacy_OPRCenterStd_20251222
% Legacy four-step OPRCenterStd route kept for old-result comparison.
%
% The current main route is Run_NewFlow_20251222 / Run_Main_OPRCenterStd_20251222.
% This legacy runner now lives under archive_previous_routes/legacy_oprcenterstd_route.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(routeDir);
run(fullfile(projectDir, 'Step01_Main_Build_OPRCenterStd_Template_20251222.m'));
run(fullfile(projectDir, 'Step02_Main_Build_OPRCenterStd_DynamicMap_20251222.m'));
run(fullfile(projectDir, 'Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20251222.m'));
run(fullfile(projectDir, 'Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222.m'));
