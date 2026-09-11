function resultFile = Run_03_GapAware_20241106()
%RUN_03_GAPAWARE_20241106 Run the formal clearance-aware identification.
cfg = Config_20241106();
add_package_paths_local(cfg);
preparedFile = expected_prepared_file_local(cfg);
if ~isfile(preparedFile)
    Run_01_Prepare_Data_20241106();
end
foundationFile = expected_foundation_file_local(cfg);
if ~isfile(foundationFile)
    foundationFile = Run_02_Foundation_NoGap_20241106();
end
Check_Program_20241106('all', false);

envState = Apply_Config_Env_20241106(cfg);
cleanupEnv = onCleanup(@() Restore_Config_Env_20241106(envState));
setenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE', foundationFile);
setenv('STEP07J_LOCAL_BUNDLE_SOURCE', 'foundation_step05_bundle');

dynamicFile = expected_dynamic_map_file_local(cfg);
if ~isfile(dynamicFile)
    run_core_script_local(fullfile(cfg.paths.preparation, ...
        'Step06_BuildGapAwareDynamicMap_20241106.m'));
else
    fprintf('Reuse exact folder-local DynamicMap:\n  %s\n', dynamicFile);
end
assert(isfile(dynamicFile), 'Step06 did not create the expected DynamicMap: %s', dynamicFile);
setenv('STEP07J_DYNAMIC_MAP_FILE', dynamicFile);
run_core_script_local(fullfile(cfg.paths.gapAware, ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106.m'));

resultFile = expected_gap_result_file_local(cfg);
assert(isfile(resultFile), 'GapAware run did not create the expected result: %s', resultFile);
PackageConfig = cfg; %#ok<NASGU>
PackageProvenance = struct('packageVersion', cfg.packageVersion, ...
    'createdAt', datetime('now'), 'preparedCaseFile', preparedFile, ...
    'foundationResultFile', foundationFile, 'dynamicMapFile', dynamicFile, ...
    'resultFile', resultFile, 'analysisTimeSec', cfg.case.analysisStartTimeSec, ...
    'finalObjective', cfg.objective.finalType); %#ok<NASGU>
save(resultFile, 'PackageConfig', 'PackageProvenance', '-append');
fprintf('GapAware result complete:\n  %s\n', resultFile);
clear cleanupEnv;
end

function add_package_paths_local(cfg)
addpath(cfg.paths.root, cfg.paths.foundation, cfg.paths.preparation, ...
    cfg.paths.gapAware, cfg.paths.utilities);
end

function run_core_script_local(scriptFile)
run(scriptFile);
end

function file = expected_prepared_file_local(cfg)
file = fullfile(cfg.paths.prepared, ['PreparedCase_', case_tag_local(cfg), '.mat']);
end

function file = expected_foundation_file_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
file = fullfile(cfg.paths.foundationResults, ...
    run_tag_local(cfg), cfg.case.dynamicCase, ...
    sprintf('Result_Step05_FoundationMainPulseAdaptiveZeroEta_B%d_%s_20241106.mat', ...
    cfg.case.targetBlade, sensorTag));
end

function file = expected_dynamic_map_file_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
file = fullfile(cfg.paths.gapRuntime, sprintf( ...
    'Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s_W%dS%d.mat', ...
    cfg.case.targetBlade, sensorTag, timeTag, cfg.window.windowBladePasses, ...
    cfg.window.slidingStepBladePasses));
end

function file = expected_gap_result_file_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
suffix = sprintf('_main_gapfixedtilt_W%dS%d_F%dto%dHz_DFpm%gHz', ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)), ...
    cfg.frequency.refineHalfWidthHz);
file = fullfile(cfg.paths.gapResults, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B%d_%s%s.mat', ...
    cfg.case.targetBlade, sensorTag, suffix));
end

function tag = case_tag_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('B%d_%s_%s_W%dS%d_F%dto%d', cfg.case.targetBlade, sensorTag, ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end

function tag = run_tag_local(cfg)
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('%s_W%dS%d_F%d_%d', timeTag, ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end
