%% Step01_Main_Build_OPRCenterStd_Template_20250527
% Build the fixed direct-template low-speed template:
% multi-threshold OPR center + synchronized standard-angle reference,
% coverage-first stable windows, SG-fit center, and gradient x-domain.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
coreDir = fullfile(routeDir, 'latest_flow_core');
run_script_local(fullfile(coreDir, 'Step01_Extract_OPR_Blade_Timing_20250527.m'));

stablePlanFile = fullfile(routeDir, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S136_20250527.csv');
if exist(stablePlanFile, 'file') ~= 2
    error(['Stable-window plan not found:\n  %s\n' ...
        'Run the archived Step01_00_Build_CoverageFirstStableWindowPlan_20250527 first.'], ...
        stablePlanFile);
end

referenceCenterTemplateFile = fullfile(routeDir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S136_OPRCenterStdRef_20250527.mat');
if exist(referenceCenterTemplateFile, 'file') ~= 2
    fprintf('Building OPRCenterStd reference-center template first...\n');
    setenv('STEP01_STABLE_WINDOW_PLAN', '');
    setenv('STEP01_TEMPLATE_SUFFIX', 'OPRCenterStdRef');
    setenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE', '');
    setenv('STEP01_CENTER_MODE', 'sgfit');
    setenv('STEP01_XRANGE_MODE', 'energy');
    run_script_local(fullfile(coreDir, 'Step01_Build_Rotating_Template_20250527.m'));
end

setenv('STEP01_STABLE_WINDOW_PLAN', stablePlanFile);
setenv('STEP01_TEMPLATE_SUFFIX', 'GradientXRange030_OPRCenterStd');
setenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE', referenceCenterTemplateFile);
setenv('STEP01_CENTER_MODE', 'sgfit');
setenv('STEP01_XRANGE_MODE', 'threshold');
setenv('STEP01_XRANGE_GRADIENT_MIN_RATIO', '0.30');
setenv('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', '0');
setenv('STEP01_XRANGE_MIN_HALF_WIDTH_MM', '2.5');
setenv('STEP01_XRANGE_MAX_HALF_WIDTH_MM', '4.2');

run_script_local(fullfile(coreDir, 'Step01_Build_Rotating_Template_20250527.m'));

function run_script_local(scriptPath)
run(scriptPath);
end
