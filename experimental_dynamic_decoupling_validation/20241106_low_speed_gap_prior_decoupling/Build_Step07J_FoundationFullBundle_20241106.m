%% Build Step07J foundation-compatible full waveform bundle cache
% Rebuilds the matching 20241106 btt_data_foundation Step05 route with
% full bundle previews, so Step07J can use the same low-speed template,
% static eta, windows, masks, weights, and high-speed waveform points.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
foundationDir = fullfile(rootDir, '20241106_btt_data_foundation');
if exist(foundationDir, 'dir') ~= 7
    error('Foundation folder not found: %s', foundationDir);
end

envNames = {'STEP05_OUTPUT_TAG', 'STEP05_STORE_BUNDLE_PREVIEW_POINTS', ...
    'STEP05_SHOW_PLOTS', 'STEP05_SAVE_FIGURES', 'STEP05_FORCE_REBUILD'};
oldEnv = cell(size(envNames));
for i = 1:numel(envNames)
    oldEnv{i} = getenv(envNames{i});
end

cleanupObj = onCleanup(@() restore_env_local(envNames, oldEnv));

setenv('STEP05_OUTPUT_TAG', 'fullbundle_for_gap_step07j');
setenv('STEP05_STORE_BUNDLE_PREVIEW_POINTS', 'inf');
setenv('STEP05_SHOW_PLOTS', '0');
setenv('STEP05_SAVE_FIGURES', '0');
setenv('STEP05_FORCE_REBUILD', '1');

oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(foundationDir);
Step05_SingleSync_DirectTemplate_Identification_20241106;

function restore_env_local(names, values)
for i = 1:numel(names)
    setenv(names{i}, values{i});
end
end
