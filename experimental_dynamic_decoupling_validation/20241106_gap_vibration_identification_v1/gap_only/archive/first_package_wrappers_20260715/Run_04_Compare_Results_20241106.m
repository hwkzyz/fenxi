function summaryFile = Run_04_Compare_Results_20241106()
%RUN_04_COMPARE_RESULTS_20241106 Compare exact formal result files only.
cfg = Config_20241106();
addpath(cfg.paths.root, cfg.paths.foundation, cfg.paths.preparation, ...
    cfg.paths.gapAware, cfg.paths.utilities);
gapResult = expected_gap_result_file_local(cfg);
if ~isfile(gapResult)
    Run_03_GapAware_20241106();
end

envState = Apply_Config_Env_20241106(cfg);
cleanupEnv = onCleanup(@() Restore_Config_Env_20241106(envState));
foundationTrend = expected_foundation_trend_local(cfg);
gapTrend = expected_gap_trend_local(cfg);
assert(isfile(foundationTrend), 'Missing exact Foundation trend: %s', foundationTrend);
assert(isfile(gapTrend), 'Missing exact GapAware trend: %s', gapTrend);
setenv('STEP07J_FVF_FOUNDATION_TREND_FILE', foundationTrend);
setenv('STEP07J_FVF_STEP07J_TREND_FILE', gapTrend);
setenv('STEP07J_FVF_OUTPUT_DIR', cfg.paths.comparison);
setenv('STEP07J_FVF_ALLOW_MAIN_FALLBACK', '1');
setenv('STEP07J_FVF_SUFFIX', 'foundation_vs_gapaware_plainrmse');
run_core_script_local(fullfile(cfg.paths.gapAware, ...
    'Compare_Step07J_FixedVsFoundation_20241106.m'));

caseTag = sprintf('B%d_S%s', cfg.case.targetBlade, ...
    sprintf('%d', cfg.case.analysisSensors));
summaryFile = fullfile(cfg.paths.comparison, sprintf( ...
    'Compare_Step07J_FixedVsFoundation_Summary_20241106_%s_foundation_vs_gapaware_plainrmse.csv', ...
    caseTag));
if ~isfile(summaryFile)
    hits = dir(fullfile(cfg.paths.comparison, ...
        'Compare_Step07J_FixedVsFoundation_Summary_20241106_*.csv'));
    assert(isscalar(hits), ...
        'Comparison did not create one unambiguous summary file.');
    summaryFile = fullfile(hits(1).folder, hits(1).name);
end
fprintf('Comparison complete:\n  %s\n', summaryFile);
clear cleanupEnv;
end

function file = expected_foundation_trend_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
file = fullfile(cfg.paths.foundationResults, ...
    run_tag_local(cfg), cfg.case.dynamicCase, ...
    sprintf('Trend_Step05_FoundationMainPulseAdaptiveZeroEta_B%d_%s_20241106.csv', ...
    cfg.case.targetBlade, sensorTag));
end

function file = expected_gap_result_file_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
file = fullfile(cfg.paths.gapResults, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B%d_%s%s.mat', ...
    cfg.case.targetBlade, sensorTag, result_suffix_local(cfg)));
end

function file = expected_gap_trend_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
file = fullfile(cfg.paths.gapResults, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_Trend_20241106_B%d_%s%s.csv', ...
    cfg.case.targetBlade, sensorTag, result_suffix_local(cfg)));
end

function suffix = result_suffix_local(cfg)
suffix = sprintf('_main_gapfixedtilt_W%dS%d_F%dto%dHz_DFpm%gHz', ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)), ...
    cfg.frequency.refineHalfWidthHz);
end

function tag = run_tag_local(cfg)
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('%s_W%dS%d_F%d_%d', timeTag, ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end

function run_core_script_local(scriptFile)
run(scriptFile);
end
