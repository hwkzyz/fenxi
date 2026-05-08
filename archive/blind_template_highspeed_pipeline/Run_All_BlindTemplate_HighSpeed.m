%% Run the blind-template high-speed pipeline
% Main thread:
%   Step 1 simulate hidden-gap high-speed data
%   -> Step 2B keep all high-speed passing samples as a point cloud
%   -> Step 3B identify dual-frequency vibration with a waveform-level
%      blind shared-template profile optimization
%   -> Step 5B visualize and summarize
clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

run(fullfile(scriptDir, 'Step_0_Config_BlindTemplate_HighSpeed.m'));
run(fullfile(scriptDir, 'Step_1_Simulate_BlindTemplate_Data.m'));
run(fullfile(scriptDir, 'Step_2B_Build_HighSpeed_PointCloud.m'));
run(fullfile(scriptDir, 'Step_3B_Waveform_Profile_Identify.m'));
run(fullfile(scriptDir, 'Step_5B_Make_Profile_Figures.m'));

fprintf('\nAll waveform-level blind-template high-speed steps finished.\n');
