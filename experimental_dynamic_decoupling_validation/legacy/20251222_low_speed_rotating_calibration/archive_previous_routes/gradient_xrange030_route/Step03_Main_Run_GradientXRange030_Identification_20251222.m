%% Step03_Main_Run_GradientXRange030_Identification_20251222
% Main 20251222 low-speed-template waveform identification.
% The main implementation is the VP top-3 synchronous waveform route.
% Original validated route: GradientXRange030 template, 50.2 s DynamicMap,
% all pulses, soft domain.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
template_suffix = strtrim(getenv('STEP01_CURRENT_TEMPLATE_SUFFIX'));
if isempty(template_suffix)
    template_suffix = 'GradientXRange030';
end
result_suffix = strtrim(getenv('STEP03_CURRENT_RESULT_SUFFIX'));
if isempty(result_suffix)
    result_suffix = 'Main_GradientXRange030';
end
template_file = fullfile(route_dir, 'output', 'templates', ...
    sprintf('Template_LowSpeedRotating_B1_S123_%s_20251222.mat', template_suffix));
dynamic_file = fullfile(route_dir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_20251222.mat');
if exist(template_file, 'file') ~= 2
    error(['Improved template not found:\n  %s\n' ...
        'Run Step01_Main_Build_GradientXRange030_Template_20251222 first.'], template_file);
end
if exist(dynamic_file, 'file') ~= 2
    error(['Main dynamic map not found:\n  %s\n' ...
        'Run Step02_Build_Dynamic_Map_20251222 with the original main settings first.'], dynamic_file);
end

setenv('STEP03_TEMPLATE_FILE', template_file);
setenv('STEP03_DYNAMIC_MAP_FILE', dynamic_file);
setenv('STEP03_RESULT_SUFFIX', result_suffix);
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
