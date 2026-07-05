%% Run_Main_OPRCenterStd_20251222
% Main 20251222 runner for the legacy OPRCenterStd direct-template route.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
cd(routeDir);

Step01_Main_Build_OPRCenterStd_Template_20251222;
Step02_Main_Build_OPRCenterStd_DynamicMap_20251222;
Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20251222;
Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222;
