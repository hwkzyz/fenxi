%% Run_Legacy_OPRCenterStd_20251222
% Legacy four-step OPRCenterStd route kept for old-result comparison.
%
% This runner is kept as a compatibility alias for the restored root route.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
run(fullfile(routeDir, 'Step01_Main_Build_OPRCenterStd_Template_20251222.m'));
run(fullfile(routeDir, 'Step02_Main_Build_OPRCenterStd_DynamicMap_20251222.m'));
run(fullfile(routeDir, 'Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20251222.m'));
run(fullfile(routeDir, 'Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222.m'));
