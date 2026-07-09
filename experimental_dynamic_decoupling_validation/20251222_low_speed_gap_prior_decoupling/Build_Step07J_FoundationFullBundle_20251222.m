%% Build Step07J foundation-compatible full waveform bundle cache
% Rebuilds the matching 20251222 btt_data_foundation Step05 route with
% full bundle previews, so Step07J can use the same low-speed template,
% static eta, windows, masks, weights, and high-speed waveform points.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
foundationDir = fullfile(rootDir, '20251222_btt_data_foundation');
if exist(foundationDir, 'dir') ~= 7
    error('Foundation folder not found: %s', foundationDir);
end

envNames = {'STEP05_OUTPUT_TAG', 'STEP05_STORE_BUNDLE_PREVIEW_POINTS', ...
    'STEP05_ANALYSIS_START_TIME', 'STEP05_TARGET_LAPS', ...
    'STEP05_WINDOW_LAPS', 'STEP05_SLIDING_STEP_LAPS', ...
    'STEP05_SHOW_PLOTS', 'STEP05_SAVE_FIGURES', 'STEP05_FORCE_REBUILD'};
oldEnv = cell(size(envNames));
for i = 1:numel(envNames)
    oldEnv{i} = getenv(envNames{i});
end

cleanupObj = onCleanup(@() restore_env_local(envNames, oldEnv)); %#ok<NASGU>

C0 = CaseConfig();
Pflow = C0.flowConfig;

setenv('STEP05_OUTPUT_TAG', 'fullb_T50p2_L20');
setenv('STEP05_STORE_BUNDLE_PREVIEW_POINTS', 'inf');
setenv('STEP05_ANALYSIS_START_TIME', sprintf('%.15g', Pflow.identification.analysisStartTimeSec));
setenv('STEP05_TARGET_LAPS', sprintf('%d', Pflow.identification.targetBladePasses));
setenv('STEP05_WINDOW_LAPS', sprintf('%d', Pflow.identification.windowBladePasses));
setenv('STEP05_SLIDING_STEP_LAPS', sprintf('%d', Pflow.identification.slidingStepBladePasses));
setenv('STEP05_SHOW_PLOTS', '0');
setenv('STEP05_SAVE_FIGURES', '0');
setenv('STEP05_FORCE_REBUILD', '1');

oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir)); %#ok<NASGU>
cd(foundationDir);
Step05_SingleSync_DirectTemplate_Identification_20251222;

function restore_env_local(names, values)
for i = 1:numel(names)
    setenv(names{i}, values{i});
end
end
