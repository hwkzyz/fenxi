function state = Apply_Config_Env_20241106(cfg)
%APPLY_CONFIG_ENV_20241106 Apply one configuration to both method programs.

names = {
    'STEP05_ANALYSIS_START_TIME_SEC','STEP05_TARGET_LAPS', ...
    'STEP05_WINDOW_LAPS','STEP05_SLIDING_STEP_LAPS', ...
    'STEP05_FREQ_SEARCH_HZ','STEP05_TARGET_BLADES','STEP05_ANALYSIS_SENSORS', ...
    'STEP05_SHOW_PLOTS','STEP05_SAVE_FIGURES','STEP05_FORCE_REBUILD', ...
    'STEP05_OUTPUT_TAG','STEP05_STATIC_ETA_FILE', ...
    'BLADE_CASE_BLADE_ID','BLADE_CASE_SENSOR_IDS', ...
    'STEP06G_TARGET_BLADE','STEP06G_ANALYSIS_SENSORS','STEP06G_GAP_SENSORS', ...
    'STEP06G_START_TIME_SEC','STEP06G_TARGET_LAPS','STEP06G_WINDOW_LAPS', ...
    'STEP06G_SLIDING_STEP_LAPS','STEP07J_ANALYSIS_SENSORS', ...
    'STEP07J_GAP_SENSORS','STEP07J_TOP_K_EO', ...
    'STEP07J_REFINE_EACH_CANDIDATE_FREQUENCY', ...
    'STEP07J_FREQ_REFINE_HALF_WIDTH_HZ','STEP07J_FINAL_OBJECTIVE', ...
    'STEP07J_FOUNDATION_STEP05_RESULT_FILE','STEP07J_LOCAL_BUNDLE_SOURCE', ...
    'STEP07J_RESULT_OUTPUT_DIR','STEP07J_RUN_MODE', ...
    'STEP07J_LOCK_GAP_EO_TO_FIXED','STEP07J_PREV_WINDOW_CANDIDATE_MODE', ...
    'STEP07J_FORCED_EO','STEP07J_STATIC_ETA_FILE','STEP07J_MAX_WINDOWS', ...
    'STEP07J_RESULT_SUFFIX','STEP07J_DYNAMIC_MAP_FILE', ...
    'STEP07J_TEMPLATE_FILE','STEP07J_CORRECTED_LIB_FILE', ...
    'STEP06G_OUTPUT_SUFFIX','STEP07J_FVF_FOUNDATION_TREND_FILE', ...
    'STEP07J_FVF_STEP07J_TREND_FILE','STEP07J_FVF_OUTPUT_DIR', ...
    'STEP07J_FVF_ALLOW_MAIN_FALLBACK','STEP07J_FVF_SUFFIX'};
state = struct('names', {names}, 'values', {cellfun(@getenv, names, 'UniformOutput', false)});

setenv('STEP05_ANALYSIS_START_TIME_SEC', num2str(cfg.case.analysisStartTimeSec, 15));
setenv('STEP05_TARGET_LAPS', sprintf('%d', cfg.window.targetBladePasses));
setenv('STEP05_WINDOW_LAPS', sprintf('%d', cfg.window.windowBladePasses));
setenv('STEP05_SLIDING_STEP_LAPS', sprintf('%d', cfg.window.slidingStepBladePasses));
setenv('STEP05_FREQ_SEARCH_HZ', sprintf('%.15g %.15g', cfg.frequency.searchHz));
setenv('STEP05_TARGET_BLADES', sprintf('%d ', cfg.case.targetBlade));
setenv('STEP05_ANALYSIS_SENSORS', sprintf('%d ', cfg.case.analysisSensors));
setenv('STEP05_SHOW_PLOTS', logical_text_local(cfg.run.showPlots));
setenv('STEP05_SAVE_FIGURES', logical_text_local(cfg.run.saveFigures));
setenv('STEP05_FORCE_REBUILD', logical_text_local(cfg.run.forceFoundationRebuild));
setenv('STEP05_OUTPUT_TAG', foundation_run_tag_local(cfg));
setenv('STEP05_STATIC_ETA_FILE', '');

setenv('BLADE_CASE_BLADE_ID', sprintf('%d', cfg.case.targetBlade));
setenv('BLADE_CASE_SENSOR_IDS', sprintf('%d ', cfg.case.analysisSensors));
setenv('STEP06G_TARGET_BLADE', sprintf('%d', cfg.case.targetBlade));
setenv('STEP06G_ANALYSIS_SENSORS', sprintf('%d ', cfg.case.analysisSensors));
setenv('STEP06G_GAP_SENSORS', sprintf('%d ', cfg.case.gapSensors));
setenv('STEP06G_START_TIME_SEC', num2str(cfg.case.analysisStartTimeSec, 15));
setenv('STEP06G_TARGET_LAPS', sprintf('%d', cfg.window.targetBladePasses));
setenv('STEP06G_WINDOW_LAPS', sprintf('%d', cfg.window.windowBladePasses));
setenv('STEP06G_SLIDING_STEP_LAPS', sprintf('%d', cfg.window.slidingStepBladePasses));
setenv('STEP07J_ANALYSIS_SENSORS', sprintf('%d ', cfg.case.analysisSensors));
setenv('STEP07J_GAP_SENSORS', sprintf('%d ', cfg.case.gapSensors));
setenv('STEP07J_TOP_K_EO', sprintf('%d', cfg.frequency.vpTopK));
setenv('STEP07J_REFINE_EACH_CANDIDATE_FREQUENCY', ...
    logical_text_local(cfg.frequency.refineEachCandidate));
setenv('STEP07J_FREQ_REFINE_HALF_WIDTH_HZ', ...
    num2str(cfg.frequency.refineHalfWidthHz, 15));
if ~strcmpi(cfg.objective.finalType, 'plain_rmse')
    error('The formal package requires cfg.objective.finalType = ''plain_rmse''.')
end
setenv('STEP07J_FINAL_OBJECTIVE', 'plain');
setenv('STEP07J_RESULT_OUTPUT_DIR', cfg.paths.gapResults);
setenv('STEP07J_RUN_MODE', 'main');
setenv('STEP07J_LOCK_GAP_EO_TO_FIXED', '0');
setenv('STEP07J_PREV_WINDOW_CANDIDATE_MODE', 'off');
setenv('STEP07J_FORCED_EO', '');
setenv('STEP07J_STATIC_ETA_FILE', '');
setenv('STEP07J_MAX_WINDOWS', '');
setenv('STEP07J_RESULT_SUFFIX', gap_result_suffix_local(cfg));
setenv('STEP07J_DYNAMIC_MAP_FILE', '');
setenv('STEP07J_TEMPLATE_FILE', '');
setenv('STEP07J_CORRECTED_LIB_FILE', '');
setenv('STEP06G_OUTPUT_SUFFIX', '');
setenv('STEP07J_FVF_FOUNDATION_TREND_FILE', '');
setenv('STEP07J_FVF_STEP07J_TREND_FILE', '');
setenv('STEP07J_FVF_OUTPUT_DIR', cfg.paths.comparison);
setenv('STEP07J_FVF_ALLOW_MAIN_FALLBACK', '1');
comparisonTimeTag = strrep(sprintf('T%07.3f', ...
    cfg.case.analysisStartTimeSec), '.', 'p');
setenv('STEP07J_FVF_SUFFIX', ...
    [comparisonTimeTag, '_foundation_vs_gapaware_plainrmse']);
end

function text = logical_text_local(value)
if value
    text = '1';
else
    text = '0';
end
end

function tag = foundation_run_tag_local(cfg)
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('%s_W%dS%d_F%d_%d', timeTag, ...
    cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end

function suffix = gap_result_suffix_local(cfg)
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
suffix = sprintf('_%s_main_gapfixedtilt_W%dS%d_F%dto%dHz_DFpm%gHz', ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)), ...
    cfg.frequency.refineHalfWidthHz);
end
