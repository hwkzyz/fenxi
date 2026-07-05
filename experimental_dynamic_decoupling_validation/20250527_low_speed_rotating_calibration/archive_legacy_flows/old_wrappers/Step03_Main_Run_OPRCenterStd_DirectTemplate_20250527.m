%% Step03_Main_Run_OPRCenterStd_DirectTemplate_20250527
% Fixed direct low-speed-template identification route.
% No static gap library is used here.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
coreDir = fullfile(routeDir, 'latest_flow_core');
templateFile = fullfile(routeDir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S136_GradientXRange030_OPRCenterStd_20250527.mat');
dynamicFile = fullfile(routeDir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S136_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20250527.mat');

if exist(templateFile, 'file') ~= 2
    error(['OPRCenterStd template not found:\n  %s\n' ...
        'Run Step01_Main_Build_OPRCenterStd_Template_20250527 first.'], templateFile);
end
if exist(dynamicFile, 'file') ~= 2
    error(['OPRCenterStd DynamicMap not found:\n  %s\n' ...
        'Run Step02_Main_Build_OPRCenterStd_DynamicMap_20250527 first.'], dynamicFile);
end

setenv('STEP03_TEMPLATE_FILE', templateFile);
setenv('STEP03_DYNAMIC_MAP_FILE', dynamicFile);
setenv('STEP03_RESULT_SUFFIX', 'Main_DirectTemplate_OPRCenterStd');
setenv('STEP03_ANALYSIS_SENSORS', '1 3 6');
setenv('STEP03_MAIN_TOP_K_EO', '3');
setenv('STEP03D_DYNAMIC_EFFECTIVE_MODE', 'gradient');
setenv('STEP03D_PULSE_MODE', 'single');
setenv('STEP03D_DOMAIN_SELECTION_MODE', 'hard');
setenv('STEP03D_DOMAIN_SOFT_MARGIN_MM', '0');
setenv('STEP03D_QUERY_GUARD_MODE', 'adaptive');
setenv('STEP03D_QUERY_GUARD_QUANTILE', '95');
setenv('STEP03D_QUERY_GUARD_SAFETY_MM', '0.05');
setenv('STEP03D_QUERY_GUARD_MIN_MM', '0.12');
setenv('STEP03D_SENSOR_ETA_LIMIT_MM', '0.03');
setenv('STEP03D_SENSOR_ETA_REG_WEIGHT_V_PER_MM', '0');
setenv('STEP03D_OVERSHOOT_PENALTY_WEIGHT', '100');

run_script_local(fullfile(coreDir, 'Step03_Main_VPTop3SynchronousWaveform_20250527.m'));

function run_script_local(scriptPath)
run(scriptPath);
end
