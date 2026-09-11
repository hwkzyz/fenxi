function Run_NestedTransferHighSpeed_20250527()
% Run identical formal high-speed identification for each nested transfer model.

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
libraryDir = fullfile(thisDir, 'outputs', 'nested_libraries');
foundationFile = fullfile(caseDir, 'results', 'fixed_gap', ...
    '20250526_2500-3500_t400', 'FixedGap_B1_S136_T001p5s.mat');
assert(isfile(foundationFile), 'Missing fixed-gap foundation file: %s', foundationFile);

modelNames = {'P0_g0','P1_g0_gain','P2_registration','P3_full_transfer'};
oldEnv = capture_env_local();
cleanup = onCleanup(@() restore_env_local(oldEnv)); %#ok<NASGU>
for im = 1:numel(modelNames)
    modelName = modelNames{im};
    libraryFile = fullfile(libraryDir, sprintf('GapLibrary_20250527_%s.mat', modelName));
    assert(isfile(libraryFile), 'Missing nested transfer library: %s', libraryFile);
    setenv('STEP07J_CORRECTED_LIB_FILE', libraryFile);
    setenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE', foundationFile);
    setenv('STEP07J_RESULT_SUFFIX', sprintf('_transfer_%s_20260830', lower(modelName)));
    setenv('STEP07J_MAX_WINDOWS', 'inf');
    setenv('STEP07J_RUN_MODE', 'main');
    setenv('STEP07J_PHASE_SAFE_EXPANSION', '0');
    setenv('STEP07J_PREV_WINDOW_CANDIDATE_MODE', 'off');
    fprintf('\n=== High-speed propagation: %s ===\n', modelName);
    mainFile = fullfile(caseDir, 'Main10_GapAware_FullWave_Identification_20250527.m');
    % Main10 is a script beginning with clear. Execute it in the base
    % workspace so it cannot destroy this driver's loop or onCleanup state.
    evalin('base', sprintf('run(''%s'')', strrep(mainFile, '''', '''''')));
end
end

function state = capture_env_local()
names = {'STEP07J_CORRECTED_LIB_FILE','STEP07J_FOUNDATION_STEP05_RESULT_FILE', ...
    'STEP07J_RESULT_SUFFIX','STEP07J_MAX_WINDOWS','STEP07J_RUN_MODE', ...
    'STEP07J_PHASE_SAFE_EXPANSION','STEP07J_PREV_WINDOW_CANDIDATE_MODE'};
state = struct();
for i = 1:numel(names)
    state.(names{i}) = getenv(names{i});
end
end

function restore_env_local(state)
names = fieldnames(state);
for i = 1:numel(names)
    setenv(names{i}, state.(names{i}));
end
end
