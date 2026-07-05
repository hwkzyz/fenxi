%% Step02_Main_Build_OPRCenterStd_DynamicMap_20250527
% Build the dynamic waveform map in the same OPRCenterStd coordinate as the
% fixed direct-template low-speed template.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
coreDir = fullfile(routeDir, 'latest_flow_core');

setenv('STEP02_TEMPLATE_SUFFIX', 'GradientXRange030_OPRCenterStd');
setenv('STEP02_DYNAMIC_SUFFIX', 'Main20L_W3S1_GradientXRange030_OPRCenterStd');
setenv('STEP02_TARGET_LAPS', '20');
setenv('STEP02_ANALYSIS_WIN_SIZE', '3');
setenv('STEP02_SLIDING_STEP', '1');

run_script_local(fullfile(coreDir, 'Step02_Build_Dynamic_Map_20250527.m'));

function run_script_local(scriptPath)
run(scriptPath);
end
