%% Step03_Run_BaseFrame_Identification_20251222
% Theory-clean diagnostic route:
% GradientXRange030 waveform/weights represented in the base DynamicMap
% x-zero convention.  This does not overwrite the validated main route.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
template_file = fullfile(route_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_BaseFrame_20251222.mat');
dynamic_file = fullfile(route_dir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_20251222.mat');

if exist(template_file, 'file') ~= 2
    run(fullfile(route_dir, 'Step01B_Build_GradientXRange030_BaseFrame_Template_20251222.m'));
end
if exist(dynamic_file, 'file') ~= 2
    error(['Base-center dynamic map not found:\n  %s\n' ...
        'Run Step02_Build_Dynamic_Map_20251222 with the original main settings first.'], dynamic_file);
end

setenv('STEP03_TEMPLATE_FILE', template_file);
setenv('STEP03_DYNAMIC_MAP_FILE', dynamic_file);
setenv('STEP03_RESULT_SUFFIX', 'Main_GradientXRange030_BaseFrame');
setenv('STEP03_ANALYSIS_SENSORS', '1 2 3');
setenv('STEP03_MAIN_TOP_K_EO', '3');
setenv('STEP03D_DYNAMIC_EFFECTIVE_MODE', 'gradient');
setenv('STEP03D_PULSE_MODE', 'all');
setenv('STEP03D_DOMAIN_SELECTION_MODE', 'soft');
setenv('STEP03D_DOMAIN_SOFT_MARGIN_MM', '0.20');
setenv('STEP03D_SENSOR_ETA_LIMIT_MM', '0.03');
setenv('STEP03D_SENSOR_ETA_REG_WEIGHT_V_PER_MM', '0');
setenv('STEP03D_OVERSHOOT_PENALTY_WEIGHT', '100');
run(fullfile(route_dir, 'Step03_Main_VPTop3SynchronousWaveform_20251222.m'));
