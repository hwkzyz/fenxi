%% Run all static-gap waveform pattern analysis steps
clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

run(fullfile(scriptDir, 'Step_0_Config_Load_StaticGap.m'));
run(fullfile(scriptDir, 'Step_1_Pointwise_Law_Analysis.m'));
run(fullfile(scriptDir, 'Step_2_Normalization_Collapse_Analysis.m'));
run(fullfile(scriptDir, 'Step_3_Feature_Law_Analysis.m'));
run(fullfile(scriptDir, 'Step_4_Fit_StandardShape_Model.m'));
run(fullfile(scriptDir, 'Step_5_LeaveOneOut_Reconstruction.m'));
run(fullfile(scriptDir, 'Step_6_Compare_Reconstruction_Models.m'));
run(fullfile(scriptDir, 'Step_7_LevelSet_Reconstruction.m'));
run(fullfile(scriptDir, 'Step_4_Summary_StaticGapPattern.m'));

fprintf('\nAll static-gap waveform pattern analysis steps finished.\n');
