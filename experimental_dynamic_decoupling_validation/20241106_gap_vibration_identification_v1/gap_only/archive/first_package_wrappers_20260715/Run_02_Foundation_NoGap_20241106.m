function resultFile = Run_02_Foundation_NoGap_20241106()
%RUN_02_FOUNDATION_NOGAP_20241106 Run the independent fixed-gap baseline.
cfg = Config_20241106();
add_package_paths_local(cfg);
preparedFile = expected_prepared_file_local(cfg);
if ~isfile(preparedFile)
    Run_01_Prepare_Data_20241106();
end
Check_Program_20241106('all', false);

envState = Apply_Config_Env_20241106(cfg);
cleanupEnv = onCleanup(@() Restore_Config_Env_20241106(envState));
runTag = run_tag_local(cfg);
setenv('STEP05_OUTPUT_TAG', runTag);
run_core_script_local(fullfile(cfg.paths.foundation, ...
    'Step05_SingleSync_DirectTemplate_Identification_20241106.m'));

sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
resultFile = fullfile(cfg.paths.foundationResults, ...
    runTag, cfg.case.dynamicCase, ...
    sprintf('Result_Step05_FoundationMainPulseAdaptiveZeroEta_B%d_%s_20241106.mat', ...
    cfg.case.targetBlade, sensorTag));
assert(isfile(resultFile), 'Foundation run did not create the expected result: %s', resultFile);
PackageConfig = cfg; %#ok<NASGU>
PackageProvenance = provenance_local(cfg, preparedFile, resultFile); %#ok<NASGU>
save(resultFile, 'PackageConfig', 'PackageProvenance', '-append');
fprintf('Foundation no-gap result complete:\n  %s\n', resultFile);
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
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('B%d_%s_%s_W%dS%d_F%dto%d', cfg.case.targetBlade, sensorTag, ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
file = fullfile(cfg.paths.prepared, ['PreparedCase_', tag, '.mat']);
end

function tag = run_tag_local(cfg)
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('%s_W%dS%d_F%d_%d', timeTag, ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end

function P = provenance_local(cfg, preparedFile, resultFile)
P = struct('packageVersion', cfg.packageVersion, 'createdAt', datetime('now'), ...
    'preparedCaseFile', preparedFile, 'resultFile', resultFile, ...
    'analysisTimeSec', cfg.case.analysisStartTimeSec, ...
    'finalObjectiveRole', 'Foundation no-gap comparison method');
end
