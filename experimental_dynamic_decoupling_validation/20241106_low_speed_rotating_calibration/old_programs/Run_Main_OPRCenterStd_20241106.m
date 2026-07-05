%% Run_Main_OPRCenterStd_20241106
% One-command OPRCenterStd direct-template pipeline for the 20241106 case.
% The default run identifies all six blades around the manual 75 s window.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
addpath(routeDir);

Step04_Main_Run_AllBlades_DirectTemplate_20241106;
