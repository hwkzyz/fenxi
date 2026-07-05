%% Step01_Main_Build_GradientXRange030_Template_20251222
% Build the main low-speed template:
% coverage-first stable windows + SG-fit center + non-parametric spline
% waveform + noise-aware gradient x-domain threshold ratio = 0.30.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
stable_plan_file = fullfile(route_dir, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S123_20251222.csv');
reference_center_template_file = fullfile(route_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_20251222.mat');

if exist(stable_plan_file, 'file') ~= 2
    error(['Stable-window plan not found:\n  %s\n' ...
        'Run the archived Step01_00_Build_CoverageFirstStableWindowPlan_20251222 first.'], ...
        stable_plan_file);
end
if exist(reference_center_template_file, 'file') ~= 2
    error(['Reference center template not found:\n  %s\n' ...
        'Build the plain low-speed template once before rebuilding GradientXRange030.'], ...
        reference_center_template_file);
end

template_suffix = strtrim(getenv('STEP01_CURRENT_TEMPLATE_SUFFIX'));
if isempty(template_suffix)
    template_suffix = 'GradientXRange030_OPRAnchored';
end

setenv('STEP01_STABLE_WINDOW_PLAN', stable_plan_file);
setenv('STEP01_TEMPLATE_SUFFIX', template_suffix);
setenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE', reference_center_template_file);
setenv('STEP01_CENTER_MODE', 'sgfit');
setenv('STEP01_XRANGE_MODE', 'threshold');
setenv('STEP01_XRANGE_GRADIENT_MIN_RATIO', '0.30');
setenv('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', '0');
setenv('STEP01_XRANGE_MIN_HALF_WIDTH_MM', '2.5');
setenv('STEP01_XRANGE_MAX_HALF_WIDTH_MM', '4.2');
run(fullfile(route_dir, 'Step01_Build_Rotating_Template_20251222.m'));
