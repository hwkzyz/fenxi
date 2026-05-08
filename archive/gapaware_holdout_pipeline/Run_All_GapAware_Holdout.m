%% Run the main gap-aware holdout pipeline
% Main thread:
%   Step 3 raw-scan initialization
%   -> Step 4 joint_wide_free local joint refinement
%   -> Step 5/6 final visualization and analysis
clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

run(fullfile(scriptDir, 'Step_0_Config_GapAware_Holdout.m'));
run(fullfile(scriptDir, 'Step_1_Build_Holdout_TemplateLibrary.m'));
run(fullfile(scriptDir, 'Step_2_Simulate_Holdout_Data.m'));
run(fullfile(scriptDir, 'Step_3_Map_And_Estimate_StaticGap.m'));
run(fullfile(scriptDir, 'Step_4_Run_VARPRO_Comparison.m'));
run(fullfile(scriptDir, 'Step_5_Make_Figures.m'));
run(fullfile(scriptDir, 'Step_6_Analyze_Library_Influence.m'));

fprintf('\nAll main-pipeline step scripts finished. Current main thread = joint_wide_free.\n');
