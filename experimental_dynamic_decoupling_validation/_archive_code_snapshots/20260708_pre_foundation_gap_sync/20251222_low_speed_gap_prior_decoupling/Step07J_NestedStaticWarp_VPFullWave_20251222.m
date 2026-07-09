%% Step07J: Nested high-speed static warp models
% The main gap_tilt route jointly identifies EO, vibration, static gap
% increment and a small high-speed tilt increment inside the direct-template
% + VP Top-K EO candidate set with the full waveform objective. VP narrows
% the EO candidates but does not fix the final EO. Comparison mode also
% computes nested static correction models on top of the low-speed template:
%
%   fixed             : no clearance increment
%   gap_only          : dg_s
%   gap_tilt          : dg_s + dmu_s*(x-tau_s)
%   gap_tilt_shift    : dg_s + dmu_s*(x-tau_s), plus sensor dTau_s
%
% Paper-facing main output uses gap_tilt. The other models are retained only
% as comparison/ablation methods.
%
% This is the current paper-facing identification route.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
C0 = CaseConfig();
Pflow = C0.flowConfig;
figDir = fullfile(outDir, 'figures_step07j_nested_static_warp_vp_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = C0.bladeId;
cfg.analysisSensors = C0.sensorIds;
cfg.freqSearchHz = [100, 1000];
cfg.eoPad = 2;
cfg.vpTopK = parse_positive_integer_env_local('STEP07J_TOP_K_EO', 3);
cfg.maxWindows = parse_positive_integer_env_local('STEP07J_MAX_WINDOWS', inf);
cfg.maxPointsPerWindow = parse_positive_integer_env_local('STEP07J_MAX_POINTS_PER_WINDOW', inf);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP07J_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DX_LIMIT_MM', 0.35);
cfg.deltaDxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_DX_LIMIT_MM', 0.05);
cfg.dxReferenceMode = lower(strtrim(getenv('STEP07J_DX_REFERENCE_MODE')));
if isempty(cfg.dxReferenceMode)
    cfg.dxReferenceMode = 'causal_vp_state';
end
if ~ismember(cfg.dxReferenceMode, {'seed_best', 'seed_median', 'causal_vp_state'})
    error('STEP07J_DX_REFERENCE_MODE must be "seed_best", "seed_median", or "causal_vp_state".');
end
cfg.causalDxRef = default_causal_dx_reference_cfg_local(cfg);
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_LIMIT_MM', 0.25);
cfg.deltaMuLimit = parse_nonnegative_numeric_env_local('STEP07J_DELTA_MU_LIMIT', 0.030);
cfg.deltaTauLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_TAU_LIMIT_MM', 0.16);
cfg.staticRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_STATIC_REG_WEIGHT_MV', 0.75);
cfg.mainModel = 'gap_tilt';
cfg.runMode = lower(strtrim(getenv('STEP07J_RUN_MODE')));
if isempty(cfg.runMode)
    cfg.runMode = 'main';
end
if ~ismember(cfg.runMode, {'main', 'comparison'})
    error('STEP07J_RUN_MODE must be "main" or "comparison".');
end
if strcmpi(cfg.runMode, 'comparison')
    cfg.modelNames = {'fixed','gap_only','gap_tilt','gap_tilt_shift'};
else
    cfg.modelNames = {'gap_tilt'};
end
cfg.mainEoScreenModel = 'none';
cfg.eoCandidateMode = lower(strtrim(getenv('STEP07J_EO_CANDIDATE_MODE')));
if isempty(cfg.eoCandidateMode)
    cfg.eoCandidateMode = 'vp';
end
if ~ismember(cfg.eoCandidateMode, {'vp', 'direct', 'direct_plus_vp'})
    error('STEP07J_EO_CANDIDATE_MODE must be "vp", "direct", or "direct_plus_vp".');
end
cfg.derivativeStepMm = 1e-3;
cfg.weightFloor = 0.05;
cfg.defaultSensorThreshold = 0.5;
cfg.domainMarginMm = parse_nonnegative_numeric_env_local('STEP07J_DOMAIN_MARGIN_MM', 0.02);
cfg.pulseSelectionMode = lower(strtrim(getenv('STEP07J_PULSE_MODE')));
if isempty(cfg.pulseSelectionMode)
    cfg.pulseSelectionMode = 'single';
end
if ~ismember(cfg.pulseSelectionMode, {'single', 'all'})
    error('STEP07J_PULSE_MODE must be "single" or "all".');
end
cfg.dynamicGradientMinRatio = parse_nonnegative_numeric_env_local('STEP07J_DYNAMIC_GRADIENT_MIN_RATIO', 0.08);
cfg.dynamicTimeGradientMinRatio = parse_nonnegative_numeric_env_local('STEP07J_DYNAMIC_TIME_GRADIENT_MIN_RATIO', 0.15);
cfg.dynamicPeakQuantile = parse_nonnegative_numeric_env_local('STEP07J_DYNAMIC_PEAK_QUANTILE', 85);
cfg.domainSelectionMode = lower(strtrim(getenv('STEP07J_DOMAIN_SELECTION_MODE')));
if isempty(cfg.domainSelectionMode)
    cfg.domainSelectionMode = 'hard';
end
if ~ismember(cfg.domainSelectionMode, {'hard', 'soft'})
    error('STEP07J_DOMAIN_SELECTION_MODE must be "hard" or "soft".');
end
cfg.domainSoftMarginMm = parse_nonnegative_numeric_env_local('STEP07J_DOMAIN_SOFT_MARGIN_MM', 0);
cfg.queryGuardFixedMm = parse_nonnegative_numeric_env_local( ...
    'STEP07J_QUERY_GUARD_MM', cfg.amplitudeLimitMm + cfg.dxLimitMm + 0.05);
cfg.queryGuardMode = lower(strtrim(getenv('STEP07J_QUERY_GUARD_MODE')));
if isempty(cfg.queryGuardMode)
    cfg.queryGuardMode = 'adaptive';
end
if ~ismember(cfg.queryGuardMode, {'fixed', 'adaptive'})
    error('STEP07J_QUERY_GUARD_MODE must be "fixed" or "adaptive".');
end
cfg.queryGuardQuantile = parse_nonnegative_numeric_env_local('STEP07J_QUERY_GUARD_QUANTILE', 95);
cfg.queryGuardSafetyMm = parse_nonnegative_numeric_env_local('STEP07J_QUERY_GUARD_SAFETY_MM', 0.05);
cfg.queryGuardMinMm = parse_nonnegative_numeric_env_local('STEP07J_QUERY_GUARD_MIN_MM', 0.12);
cfg.queryGuardMaxMm = parse_nonnegative_numeric_env_local('STEP07J_QUERY_GUARD_MAX_MM', cfg.queryGuardFixedMm);
cfg.phaseSafeExpansion = parse_logical_env_local('STEP07J_PHASE_SAFE_EXPANSION', true);
cfg.phaseSafeMarginMm = parse_nonnegative_numeric_env_local('STEP07J_PHASE_SAFE_MARGIN_MM', 0.03);
cfg.phaseSafeFallbackMode = lower(strtrim(getenv('STEP07J_PHASE_SAFE_FALLBACK_MODE')));
if isempty(cfg.phaseSafeFallbackMode)
    cfg.phaseSafeFallbackMode = 'linear_vp';
end
if ~ismember(cfg.phaseSafeFallbackMode, {'off', 'linear_vp'})
    error('STEP07J_PHASE_SAFE_FALLBACK_MODE must be "off" or "linear_vp".');
end
cfg.prevWindowRefreshEvery = floor(parse_nonnegative_numeric_env_local( ...
    'STEP07J_PREV_WINDOW_REFRESH_EVERY', 5));
cfg.overshootPenaltyWeight = parse_nonnegative_numeric_env_local('STEP07J_OVERSHOOT_PENALTY_WEIGHT', 100);
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 450, 'MaxFunEvals', 1800, ...
    'TolX', 1e-5, 'TolFun', 1e-5);

sensorOverride = strtrim(getenv('STEP07J_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07J_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
resultSuffix = sanitize_suffix_local(strtrim(getenv('STEP07J_RESULT_SUFFIX')));
if isempty(resultSuffix)
    if strcmpi(cfg.runMode, 'comparison')
        resultSuffix = '_method_compare';
    else
        resultSuffix = '_main_gaptilt';
    end
end

rotDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
rotTemplateDir = fullfile(rotDir, 'output', 'templates');
rotDynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');
gapDynamicDir = outDir;

correctedLibOverride = strtrim(getenv('STEP07J_CORRECTED_LIB_FILE'));
if ~isempty(correctedLibOverride)
    correctedLibFile = correctedLibOverride;
else
    correctedLibFile = find_corrected_gap_library_file_local(outDir, cfg.targetBlade, sensorTag);
end
templateOverride = strtrim(getenv('STEP07J_TEMPLATE_FILE'));
if ~isempty(templateOverride)
    templateFile = templateOverride;
else
    templateFile = find_low_speed_template_file_local(rotTemplateDir, cfg.targetBlade, cfg.analysisSensors);
end
dynamicOverride = strtrim(getenv('STEP07J_DYNAMIC_MAP_FILE'));
if ~isempty(dynamicOverride)
    if exist(dynamicOverride, 'file') ~= 2
        error('STEP07J_DYNAMIC_MAP_FILE does not exist: %s', dynamicOverride);
    end
    dynamicFile = dynamicOverride;
else
    dynamicFile = find_dynamic_map_file_local(gapDynamicDir, rotDynamicDir, ...
        cfg.targetBlade, cfg.analysisSensors, Pflow);
end
directMainFile = '';

if ~isfile(correctedLibFile)
    error('Run Step06I first. Missing file: %s', correctedLibFile);
end

St = load(templateFile, 'Template');
Sd = load(dynamicFile, 'DynamicMap');
Sc = load(correctedLibFile, 'CorrectedGapLibrary');
Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
DynamicMap = filter_dynamic_map_sensors_local(Sd.DynamicMap, cfg.analysisSensors, sensorTag);
coordinateCheck = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysisSensors, 1e-6);
CorrectedGapLibrary = Sc.CorrectedGapLibrary;
responseSurface = CorrectedGapLibrary.responseSurface;

numWindows = min(numel(DynamicMap.Window), cfg.maxWindows);
fprintf('\n=== Step07J: nested static warp VP/full-wave ===\n');
fprintf('Template: %s\n', templateFile);
fprintf('DynamicMap: %s\n', dynamicFile);
fprintf('Corrected gap library: %s\n', correctedLibFile);
fprintf('Active case: B%d, sensors %s, start %.6f s, laps %d, window %d, step %d.\n', ...
    cfg.targetBlade, mat2str(cfg.analysisSensors), ...
    Pflow.identification.analysisStartTimeSec, Pflow.identification.targetBladePasses, ...
    Pflow.identification.windowBladePasses, Pflow.identification.slidingStepBladePasses);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinateCheck.max_abs_delta_mm, coordinateCheck.status);
if ~coordinateCheck.is_consistent
    warning('Template/DynamicMap x-center mismatch. Do not use this run as the unified GradientXRange030 main result.');
end
fprintf('Windows: %d, sensors: %s, VP Top-K: %d, EO candidate mode: %s\n', ...
    numWindows, mat2str(cfg.analysisSensors), cfg.vpTopK, cfg.eoCandidateMode);
fprintf('Point selection: local Step03 direct-template bundle logic; no external Step03 result is loaded.\n');
fprintf('Phase-safe expansion: %d, margin %.3f mm, fallback %s.\n', ...
    cfg.phaseSafeExpansion, cfg.phaseSafeMarginMm, cfg.phaseSafeFallbackMode);
fprintf(['Pulse/domain: %s + %s, soft margin %.3f mm, query guard %s ' ...
    '(fixed %.3f, min %.3f, max %.3f, q%.1f + %.3f mm)\n'], ...
    cfg.pulseSelectionMode, cfg.domainSelectionMode, cfg.domainSoftMarginMm, ...
    cfg.queryGuardMode, cfg.queryGuardFixedMm, cfg.queryGuardMinMm, ...
    cfg.queryGuardMaxMm, cfg.queryGuardQuantile, cfg.queryGuardSafetyMm);
fprintf(['Models: %s | main %s, EO screen %s | dx limit %.3f mm, ' ...
    'dg %.3f mm, dmu %.3f, dtau %.3f mm\n'], ...
    strjoin(cfg.modelNames, ', '), cfg.mainModel, cfg.mainEoScreenModel, ...
    cfg.dxLimitMm, cfg.deltaGapLimitMm, cfg.deltaMuLimit, cfg.deltaTauLimitMm);

WindowResult = struct([]);
trendRows = cell(numWindows, 1);
bestIdx = 1;
bestRmse = inf;
prevLocalDirectPhaseInfo = empty_phase_safe_reference_local();
dxReferenceTracker = init_causal_dx_reference_tracker_local(cfg.causalDxRef);

for iw = 1:numWindows
    Wmap = DynamicMap.Window(iw);
    coreDirectBundle = build_local_direct_template_bundle_local(Wmap, Template, cfg, []);
    coreBundleFull = build_gap_bundle_from_direct_bundle_local( ...
        coreDirectBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
    eoCandidates = build_eo_candidates_local(coreBundleFull.rotFreqMeanHz, cfg.freqSearchHz, cfg.eoPad);
    localDirectCoreSeedTable = solve_local_direct_template_seed_scan_local(coreDirectBundle, Template, eoCandidates, cfg);
    phaseInfo = prevLocalDirectPhaseInfo;
    if cfg.prevWindowRefreshEvery > 0 && mod(iw - 1, cfg.prevWindowRefreshEvery) == 0
        phaseInfo = empty_phase_safe_reference_local();
        phaseInfo.reason = 'scheduled_refresh';
    end
    if cfg.phaseSafeExpansion && ~phaseInfo.usable && strcmpi(cfg.phaseSafeFallbackMode, 'linear_vp')
        phaseInfo = build_phase_safe_reference_local(localDirectCoreSeedTable, iw, cfg, phaseInfo.reason);
    end

    bundleFull = coreBundleFull;
    selectedDirectBundle = coreDirectBundle;
    if cfg.phaseSafeExpansion && phaseInfo.usable
        expandedDirectBundle = build_local_direct_template_bundle_local(Wmap, Template, cfg, phaseInfo);
        expandedBundleFull = build_gap_bundle_from_direct_bundle_local( ...
            expandedDirectBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        if expandedBundleFull.pointCount > coreBundleFull.pointCount
            bundleFull = expandedBundleFull;
            selectedDirectBundle = expandedDirectBundle;
        end
    end
    bundle = decimate_bundle_local(bundleFull, cfg.maxPointsPerWindow);
    seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg);
    [dxReferenceTracker, dxReferenceStep] = update_causal_dx_reference_tracker_local( ...
        dxReferenceTracker, seedTable, bundle, cfg.causalDxRef);
    selectedEO = select_eo_candidates_local(seedTable, cfg.vpTopK, [], iw, cfg.eoCandidateMode);
    mainCandidateEO = selectedEO;

    modelFits = struct();
    for im = 1:numel(cfg.modelNames)
        modeName = cfg.modelNames{im};
        cfgWindow = cfg;
        cfgWindow.CausalDxReference = dxReferenceStep;
        cfgWindow = configure_dx_reference_local(cfgWindow, seedTable, selectedEO, modeName);
        modelFits.(modeName) = refine_static_warp_fit_local(bundle, seedTable, selectedEO, cfgWindow, modeName);
    end
    wr = pack_window_result_local(iw, Wmap, bundleFull, seedTable, selectedEO, ...
        mainCandidateEO, modelFits, dxReferenceStep);
    if iw == 1
        WindowResult = repmat(wr, numWindows, 1);
    else
        WindowResult(iw) = wr;
    end
    trendRows{iw} = make_trend_row_local(wr, []);
    localDirectSelectedSeedTable = solve_local_direct_template_seed_scan_local(selectedDirectBundle, Template, eoCandidates, cfg);
    prevLocalDirectPhaseInfo = build_phase_safe_reference_local(localDirectSelectedSeedTable, iw, cfg, 'local_direct_selected');
    if modelFits.(cfg.mainModel).weightedRmseMv < bestRmse
        bestRmse = modelFits.(cfg.mainModel).weightedRmseMv;
        bestIdx = iw;
    end
    if strcmpi(cfg.runMode, 'comparison')
        fprintf('Window %02d/%02d laps %s: fixed %.2f | gap %.2f | gap+tilt %.2f | gap+tilt+shift %.2f mV\n', ...
            iw, numWindows, mat2str(Wmap.lap_range), modelFits.fixed.weightedRmseMv, ...
            modelFits.gap_only.weightedRmseMv, modelFits.gap_tilt.weightedRmseMv, ...
            modelFits.gap_tilt_shift.weightedRmseMv);
    else
        fprintf('Window %02d/%02d laps %s: gap_tilt %.2f mV\n', ...
            iw, numWindows, mat2str(Wmap.lap_range), modelFits.gap_tilt.weightedRmseMv);
    end
end

Trend = vertcat(trendRows{:});
Summary = build_summary_table_local(Trend, bestIdx, cfg);

Result = struct();
Result.dataset = '20251222';
if strcmpi(cfg.runMode, 'comparison')
    Result.method = 'gap_tilt_comparison_with_nested_static_warp_ablation_vp_full_wave';
else
    Result.method = 'gap_tilt_main_with_nested_static_warp_vp_full_wave';
end
Result.description = ['Main result uses gap_tilt: scalar gap increment plus a small high-speed x-dependent tilt increment. ' ...
    'Step07J builds the Step03-style direct-template point set locally from Template and DynamicMap, then converts those selected points to the gap-prior response model. No rotating_calibration Step03 result, EO trend, or stored bundle is loaded as prior information. ' ...
    'VP Top-K provides only the EO candidate set; gap_tilt jointly identifies EO, vibration, gap increment and tilt inside that set with the full waveform objective. ' ...
    'Fixed, gap_only and gap_tilt_shift are comparison/ablation models only.'];
Result.cfg = cfg;
Result.flowConfig = Pflow;
Result.templateFile = templateFile;
Result.dynamicFile = dynamicFile;
Result.correctedLibFile = correctedLibFile;
Result.directMainFile = directMainFile;
Result.CoordinateCheck = coordinateCheck;
Result.WindowResult = WindowResult;
Result.Trend = Trend;
Result.Summary = Summary;
Result.BestWindowIndex = bestIdx;
Result.MainModel = cfg.mainModel;
Result.BestWindow = WindowResult(bestIdx);

matFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20251222_B%d_%s.mat', ...
    cfg.targetBlade, [sensorTag resultSuffix]));
trendCsv = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_20251222_B%d_%s.csv', ...
    cfg.targetBlade, [sensorTag resultSuffix]));
summaryCsv = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Summary_20251222_B%d_%s.csv', ...
    cfg.targetBlade, [sensorTag resultSuffix]));
figTrend = fullfile(figDir, sprintf('Step07J_Trend_20251222_B%d_%s.png', cfg.targetBlade, [sensorTag resultSuffix]));
figBest = fullfile(figDir, sprintf('Step07J_BestWindow_Collapse_20251222_B%d_%s.png', cfg.targetBlade, [sensorTag resultSuffix]));
figGap = fullfile(figDir, sprintf('Step07J_StaticWarp_20251222_B%d_%s.png', cfg.targetBlade, [sensorTag resultSuffix]));

save(matFile, 'Result', 'Trend', 'Summary', '-v7.3');
writetable(Trend, trendCsv);
writetable(Summary, summaryCsv);
plot_trend_local(Trend, figTrend);
plot_best_window_collapse_local(Result.BestWindow, cfg, figBest);
plot_static_warp_local(Result, Trend, cfg, figGap);

fprintf('\nStep07J complete.\n');
disp(Summary);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, trendCsv, figTrend);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_BaseFrame_20251222.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
error('No low-speed template file found under %s for %s.', templateDir, sensorTag);
end

function dynamicFile = find_dynamic_map_file_local(gapDynamicDir, legacyDynamicDir, targetBlade, sensorIds, Pflow)
sensorTag = ['S', sprintf('%d', sensorIds)];
caseTag = sprintf('B%d_%s', targetBlade, sensorTag);
cacheKey = make_highmap_cache_key_step07j_local( ...
    Pflow.identification.resonanceRegionShortTag, ...
    Pflow.identification.analysisStartTimeSec, ...
    Pflow.identification.targetBladePasses, ...
    Pflow.identification.windowBladePasses, ...
    Pflow.identification.slidingStepBladePasses, ...
    Pflow.identification.pulseWindowSec);
exactStep06File = fullfile(gapDynamicDir, sprintf( ...
    'Step06_DynamicMap_20251222_%s_%s.mat', caseTag, cacheKey));
if isfile(exactStep06File)
    dynamicFile = exactStep06File;
    return;
end
patterns = {
    sprintf('Step06_DynamicMap_20251222_%s_%s*.mat', caseTag, cacheKey)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows*OPRCenterStd*_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    if startsWith(patterns{ip}, 'Step06_DynamicMap')
        files = dir(fullfile(gapDynamicDir, patterns{ip}));
    else
        files = dir(fullfile(legacyDynamicDir, patterns{ip}));
    end
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        dynamicFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
error(['No DynamicMap file found for %s. Run Step06_Build_HighMap_20251222 first ' ...
    'for start %.6f s, or set STEP07J_DYNAMIC_MAP_FILE explicitly.'], ...
    caseTag, Pflow.identification.analysisStartTimeSec);
end

function key = make_highmap_cache_key_step07j_local(regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec)
key = sprintf('%s_T%09.3f_L%02d_W%02d_S%02d_PW%07.5f', ...
    regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec);
key = strrep(key, '.', 'p');
key = strrep(key, '+', '');
key = strrep(key, '-', 'm');
end

function correctedLibFile = find_corrected_gap_library_file_local(outDir, targetBlade, sensorTag)
patterns = {
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_%s_OPRCenterStd.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_%s*.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        correctedLibFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
correctedLibFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_%s.mat', ...
    targetBlade, sensorTag));
end

function check = check_template_dynamic_xcenter_local(Template, DynamicMap, analysisSensors, toleranceMm)
rows = repmat(struct('sensor_id', NaN, 'template_xc_mm', NaN, ...
    'dynamic_xc_mm', NaN, 'delta_mm', NaN), numel(analysisSensors), 1);
for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    rows(is).sensor_id = sid;
    if ~isempty(Tpl) && isfield(Tpl, 'xc')
        rows(is).template_xc_mm = Tpl.xc;
    end
    rows(is).dynamic_xc_mm = infer_dynamic_xcenter_local(DynamicMap, sid);
    rows(is).delta_mm = rows(is).template_xc_mm - rows(is).dynamic_xc_mm;
end
delta = [rows.delta_mm];
check = struct();
check.table = struct2table(rows);
check.tolerance_mm = toleranceMm;
check.max_abs_delta_mm = max(abs(delta), [], 'omitnan');
check.is_consistent = isfinite(check.max_abs_delta_mm) && check.max_abs_delta_mm <= toleranceMm;
if check.is_consistent
    check.status = 'consistent';
else
    check.status = 'mismatch';
end
if isfield(DynamicMap, 'TemplateFileForXRel')
    check.dynamic_template_file_for_xrel = DynamicMap.TemplateFileForXRel;
end
end

function xc = infer_dynamic_xcenter_local(DynamicMap, sid)
xc = NaN;
for iw = 1:numel(DynamicMap.Window)
    S = DynamicMap.Window(iw).Sensor;
    idx = find([S.sensor_id] == sid, 1, 'first');
    if isempty(idx) || isempty(S(idx).x_abs) || isempty(S(idx).x_rel)
        continue;
    end
    v = S(idx).x_abs(:) - S(idx).x_rel(:);
    xc = median(v(isfinite(v)), 'omitnan');
    return;
end
end

function bundle = build_local_direct_template_bundle_local(Wmap, Template, cfg, phaseRef)
if nargin < 4
    phaseRef = [];
end
usePhaseSafe = isstruct(phaseRef) && isfield(phaseRef, 'usable') && ...
    phaseRef.usable && cfg.phaseSafeExpansion;
X = []; T = []; V = []; W = []; Theta = []; S = []; F0 = []; Fx = []; sensorIndex = [];
windowMeta = repmat(struct('sensor_id', NaN, 'valid_points', 0, ...
    'pulse_segment_count', 0, 'dynamic_effective_points', 0, ...
    'query_safe_points', 0, 'query_guard_mm', NaN, ...
    'phase_safe_expansion', false, 'phase_safe_points', 0, ...
    'template_domain', [NaN, NaN], 'selected_x_range', [NaN, NaN], ...
    'shape_rmse', NaN), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tpl = get_template_sensor_local(Template, sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    xRaw = D.x_rel(:);
    tRaw = D.t(:);
    thetaRaw = D.theta(:);
    vRaw = D.V(:);
    vRawMv = (vRaw - Tpl.baseline) * 1000;
    threshold = resolve_threshold_local(Tpl, cfg);
    thresholdMv = (threshold - Tpl.baseline) * 1000;
    f0Raw = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw, 'pchip', NaN);
    fxRaw = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw, 'pchip', NaN);
    xDomain = resolve_template_domain_local(Tpl);

    [maskPulse, pulseSegmentCount] = build_pulse_mask_with_count_local( ...
        vRaw, threshold, cfg.pulseSelectionMode);
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, cfg);
    maskDomain = xRaw >= xDomain(1) + cfg.domainMarginMm & ...
        xRaw <= xDomain(2) - cfg.domainMarginMm;
    maskFinite = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRaw) & isfinite(thetaRaw);
    maskBase = maskPulse & maskEffective & maskDomain & maskFinite;
    queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, maskBase, cfg);
    maskQuerySafe = xRaw >= xDomain(1) + queryGuardMm & ...
        xRaw <= xDomain(2) - queryGuardMm;
    maskPhaseSafe = false(size(maskQuerySafe));
    if usePhaseSafe
        maskPhaseSafe = build_phase_safe_query_mask_local(xRaw, thetaRaw, xDomain, is, phaseRef, cfg);
        if nnz(maskBase & maskPhaseSafe) >= 8
            maskQuerySafe = maskPhaseSafe;
        end
    end
    mask = maskBase;
    if strcmpi(cfg.domainSelectionMode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
    if nnz(mask) < 8
        mask = maskEffective & maskDomain & maskFinite;
        if strcmpi(cfg.domainSelectionMode, 'hard')
            safeMask = mask & maskQuerySafe;
            if nnz(safeMask) >= 8
                mask = safeMask;
            end
        end
    end
    if nnz(mask) < 8
        continue;
    end

    xSel = xRaw(mask);
    tSel = tRaw(mask);
    vSel = vRaw(mask);
    thetaSel = thetaRaw(mask);
    f0Sel = f0Raw(mask);
    fxSel = fxRaw(mask);
    wTotal = build_local_direct_weight_local(Tpl, xSel, tSel, vSel, xDomain, queryGuardMm, cfg);
    vStatic = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xSel, 'pchip', NaN);
    shapeRmse = sqrt(mean((vSel - vStatic).^2, 'omitnan'));

    X = [X; xSel]; %#ok<AGROW>
    T = [T; tSel]; %#ok<AGROW>
    V = [V; vSel]; %#ok<AGROW>
    W = [W; wTotal]; %#ok<AGROW>
    Theta = [Theta; thetaSel]; %#ok<AGROW>
    S = [S; repmat(sid, nnz(mask), 1)]; %#ok<AGROW>
    F0 = [F0; f0Sel]; %#ok<AGROW>
    Fx = [Fx; fxSel]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, nnz(mask), 1)]; %#ok<AGROW>

    windowMeta(is).sensor_id = sid;
    windowMeta(is).valid_points = nnz(mask);
    windowMeta(is).pulse_segment_count = pulseSegmentCount;
    windowMeta(is).dynamic_effective_points = nnz(maskEffective);
    windowMeta(is).query_safe_points = nnz(mask & maskQuerySafe);
    windowMeta(is).query_guard_mm = queryGuardMm;
    windowMeta(is).phase_safe_expansion = usePhaseSafe;
    windowMeta(is).phase_safe_points = nnz(maskBase & maskPhaseSafe);
    windowMeta(is).template_domain = xDomain;
    windowMeta(is).selected_x_range = [min(xSel), max(xSel)];
    windowMeta(is).shape_rmse = shapeRmse;
end
if isempty(T)
    error('No valid local direct-template points were constructed for window %d.', Wmap.window_id);
end
bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.sensor_index = sensorIndex;
bundle.sensor_ids = cfg.analysisSensors(:).';
bundle.window_meta = windowMeta([windowMeta.valid_points] > 0);
bundle.valid_segment_count = sum([bundle.window_meta.pulse_segment_count]);
bundle.point_count = numel(T);
if usePhaseSafe
    bundle.selection_pass = 'local_direct_phase_safe_expanded';
else
    bundle.selection_pass = 'local_direct_core_hard_adaptive';
end
bundle.phase_safe_margin_mm = cfg.phaseSafeMarginMm;
bundle.time_window = [min(T), max(T)];
bundle.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
bundle.rot_rpm_mean = Wmap.rot_rpm_mean;
bundle.window_id = Wmap.window_id;
bundle.lap_range = Wmap.lap_range;
bundle.query_guard_mm = max([bundle.window_meta.query_guard_mm], [], 'omitnan');
end

function seedTable = solve_local_direct_template_seed_scan_local(bundle, Template, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, ...
    'linearRmseMv', inf, 'fallbackRmseMv', inf), numel(eoCandidates), 1);
F0v = nan(size(bundle.X));
Fxv = nan(size(bundle.X));
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    mask = bundle.sensor_index == is;
    Tpl = get_template_sensor_local(Template, sid);
    F0v(mask) = interp1(Tpl.x_grid(:), Tpl.v_grid(:), bundle.X(mask), 'pchip', NaN);
    Fxv(mask) = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), bundle.X(mask), 'pchip', NaN);
end
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo .* bundle.Theta);
    c1 = cos(eo .* bundle.Theta);
    basis = [-Fxv(:), -Fxv(:) .* s1, -Fxv(:) .* c1];
    y = bundle.V(:) - F0v(:);
    valid = isfinite(y) & all(isfinite(basis), 2) & isfinite(bundle.W(:));
    sw = sqrt(max(bundle.W(valid), cfg.weightFloor));
    coeff = (basis(valid, :) .* sw) \ (y(valid) .* sw);
    dx = coeff(1);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y(valid) - basis(valid, :) * coeff;
    rmseMv = 1000 * sqrt(sum(bundle.W(valid) .* linRes.^2) / max(sum(bundle.W(valid)), eps));
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx = dx;
    rows(i).linearRmseMv = rmseMv;
    rows(i).fallbackRmseMv = rmseMv;
end
[~, order] = sortrows([[rows.fallbackRmseMv].', [rows.linearRmseMv].', [rows.EO].'], [1, 2, 3]);
seedTable = struct2table(rows(order));
end

function info = build_phase_safe_reference_local(seedTable, windowId, cfg, fallbackReason)
info = empty_phase_safe_reference_local();
if nargin < 4 || isempty(fallbackReason)
    fallbackReason = 'not_specified';
end
if isempty(seedTable) || height(seedTable) < 1
    info.reason = fallbackReason;
    return;
end
seed = seedTable(1, :);
if any(~isfinite([seed.EO, seed.A, seed.phi, seed.dx, seed.fallbackRmseMv]))
    info.reason = [fallbackReason, '_then_seed_nonfinite'];
    return;
end
info.usable = true;
info.EO = round(seed.EO);
info.A = min(abs(seed.A), cfg.amplitudeLimitMm);
info.phi = wrap_to_pi_local(seed.phi);
info.dx = max(min(seed.dx, cfg.dxLimitMm), -cfg.dxLimitMm);
info.deltaTauMm = zeros(1, numel(cfg.analysisSensors));
info.weightedRmseMv = seed.fallbackRmseMv;
info.source = 'local_direct_linear_vp_seed';
info.reason = [fallbackReason, '_then_local_direct_linear_vp'];
info.sourceWindowId = windowId;
info.usedPreviousWindow = false;
end

function info = empty_phase_safe_reference_local()
info = struct('usable', false, 'EO', NaN, 'A', NaN, 'phi', NaN, ...
    'dx', NaN, 'deltaTauMm', [], 'weightedRmseMv', NaN, ...
    'source', 'none', 'reason', 'not_evaluated', ...
    'sourceWindowId', NaN, 'usedPreviousWindow', false);
end

function xDomain = resolve_template_domain_local(Tpl)
if isfield(Tpl, 'x_domain') && numel(Tpl.x_domain) >= 2 && all(isfinite(Tpl.x_domain(1:2)))
    xDomain = Tpl.x_domain(1:2);
else
    xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
end
end

function [mask, segmentCount] = build_pulse_mask_with_count_local(v, threshold, modeName)
if strcmpi(modeName, 'single')
    [mask, segmentCount] = isolate_main_pulse_local(v, threshold);
else
    [mask, segmentCount] = isolate_all_pulses_local(v, threshold);
end
end

function maskSafe = build_phase_safe_query_mask_local(x, theta, xDomain, sensorLocalIndex, phaseRef, cfg)
if ~isstruct(phaseRef) || ~isfield(phaseRef, 'usable') || ~phaseRef.usable
    maskSafe = false(size(x));
    return;
end
dtau = 0;
if isfield(phaseRef, 'deltaTauMm') && numel(phaseRef.deltaTauMm) >= sensorLocalIndex
    dtau = phaseRef.deltaTauMm(sensorLocalIndex);
end
uEst = phaseRef.A .* sin(phaseRef.EO .* theta(:) + phaseRef.phi);
xQuery = x(:) - phaseRef.dx - dtau - uEst;
margin = cfg.phaseSafeMarginMm;
maskSafe = isfinite(xQuery) & ...
    xQuery >= xDomain(1) + margin & xQuery <= xDomain(2) - margin;
maskSafe = reshape(maskSafe, size(x));
end

function wTotal = build_local_direct_weight_local(Tpl, x, t, v, xDomain, queryGuardMm, cfg)
wEdge = build_edge_weight_local(t, v, cfg.weightFloor);
if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
    wTemplate = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x(:), 'linear', cfg.weightFloor);
else
    wTemplate = ones(size(x(:)));
end
wDomain = build_domain_soft_weight_local(x, xDomain, cfg.domainSoftMarginMm, cfg.weightFloor);
wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMm, cfg.domainSoftMarginMm, cfg.weightFloor);
wGradient = build_template_gradient_weight_local(Tpl, x, cfg.domainSelectionMode, cfg.weightFloor);
wTotal = max(cfg.weightFloor, wEdge(:) .* wTemplate(:) .* wDomain(:) .* wQuery(:) .* wGradient(:));
if max(wTotal) > 0
    wTotal = max(cfg.weightFloor, wTotal ./ max(wTotal));
end
end

function w = build_edge_weight_local(t, v, floorW)
t = t(:);
v = v(:);
w = ones(size(t));
if numel(t) >= 5
    dv = abs(gradient(v));
    scale = max(dv, [], 'omitnan');
    if isfinite(scale) && scale > 0
        w = max(floorW, dv ./ scale);
    end
end
end

function w = build_domain_soft_weight_local(x, xDomain, marginMm, floorW)
x = x(:);
w = ones(size(x));
if marginMm <= 0
    return;
end
dist = min(x - xDomain(1), xDomain(2) - x);
w = min(1, max(floorW, dist ./ marginMm));
end

function w = build_query_guard_soft_weight_local(x, xDomain, queryGuardMm, marginMm, floorW)
x = x(:);
w = ones(size(x));
softWidth = max(marginMm, 1e-6);
innerLo = xDomain(1) + queryGuardMm;
innerHi = xDomain(2) - queryGuardMm;
dist = min(x - innerLo, innerHi - x);
edge = dist < softWidth;
w(edge) = max(floorW, min(1, dist(edge) ./ softWidth));
w(~isfinite(w)) = floorW;
end

function w = build_template_gradient_weight_local(Tpl, x, modeName, floorW)
x = x(:);
if ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    w = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x, 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0 || strcmpi(modeName, 'soft')
    w = ones(size(x));
else
    w = max(floorW, g ./ gMax);
end
end

function bundle = build_gap_bundle_from_direct_bundle_local(Bsrc, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
if isfield(Bsrc, 'sensor_ids')
    sensorIds = Bsrc.sensor_ids(:).';
else
    sensorIds = cfg.analysisSensors(:).';
end
X = Bsrc.X(:);
T = Bsrc.T(:);
Theta = Bsrc.Theta(:);
W = max(Bsrc.W(:), cfg.weightFloor);
sensorIndex = Bsrc.sensor_index(:);
S = sensorIds(sensorIndex).';
V = nan(size(X));
F0 = nan(size(X));
Fx = nan(size(X));
g0BySensor = NaN(numel(sensorIds), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = sensorIndex == is;
    Tpl = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    V(mask) = (Bsrc.V(mask) - Tpl.baseline) * 1000;
    F0(mask) = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask));
    Fx(mask) = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask), cfg.derivativeStepMm);
    g0BySensor(is) = corr.g0Mm;
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    sensorInfo(is).pointCount = nnz(mask);
end
valid = isfinite(X) & isfinite(T) & isfinite(Theta) & isfinite(W) & ...
    isfinite(V) & isfinite(F0) & isfinite(Fx) & isfinite(sensorIndex);
bundle = struct();
bundle.X = X(valid);
bundle.T = T(valid);
bundle.TRel = T(valid) - min(T(valid));
bundle.V = V(valid);
bundle.W = W(valid);
bundle.Theta = Theta(valid);
bundle.S = S(valid);
bundle.sensorIndex = sensorIndex(valid);
bundle.sensorIds = sensorIds;
bundle.F0 = F0(valid);
bundle.Fx = Fx(valid);
bundle.Template = Template;
bundle.CorrectedGapLibrary = CorrectedGapLibrary;
bundle.responseSurface = responseSurface;
bundle.g0BySensor = g0BySensor;
bundle.sensorInfo = sensorInfo;
bundle.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
bundle.rotRpmMean = Wmap.rot_rpm_mean;
bundle.windowId = Wmap.window_id;
bundle.lapRange = Wmap.lap_range;
bundle.pointCount = numel(bundle.T);
end

function Template = filter_template_sensors_local(Template, sensorIds, sensorTag)
available = [Template.Sensor.sensor_id];
idx = zeros(size(sensorIds));
for i = 1:numel(sensorIds)
    hit = find(available == sensorIds(i), 1);
    if isempty(hit)
        error('Template does not contain CH%d.', sensorIds(i));
    end
    idx(i) = hit;
end
Template.Sensor = Template.Sensor(idx);
Template.SensorIDs = sensorIds;
Template.SensorTag = sensorTag;
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, sensorIds, sensorTag)
for iw = 1:numel(DynamicMap.Window)
    available = [DynamicMap.Window(iw).Sensor.sensor_id];
    idx = zeros(size(sensorIds));
    for i = 1:numel(sensorIds)
        hit = find(available == sensorIds(i), 1);
        if isempty(hit)
            error('DynamicMap window %d does not contain CH%d.', iw, sensorIds(i));
        end
        idx(i) = hit;
    end
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(idx);
end
DynamicMap.SensorIDs = sensorIds;
DynamicMap.SensorTag = sensorTag;
end

function bundle = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
X = []; T = []; V = []; W = []; Theta = []; S = []; sensorIndex = [];
F0 = []; Fx = []; g0BySensor = NaN(numel(cfg.analysisSensors), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tpl = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    xRaw = D.x_rel(:);
    tRaw = D.t(:);
    thetaRaw = D.theta(:);
    vRaw = D.V(:);
    vRawMv = (D.V(:) - Tpl.baseline) * 1000;
    wRaw = max(D.W(:), cfg.weightFloor);
    threshold = resolve_threshold_local(Tpl, cfg);
    thresholdMv = (threshold - Tpl.baseline) * 1000;

    f0Raw = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw);
    fxRaw = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw, cfg.derivativeStepMm);
    maskPulse = build_pulse_mask_local(vRawMv, thresholdMv, cfg.pulseSelectionMode);
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, cfg);
    maskDomain = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRawMv) & isfinite(thetaRaw);
    mask = maskPulse & maskEffective & maskDomain;
    queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, mask, cfg);
    xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    maskQuerySafe = xRaw >= xDomain(1) + queryGuardMm & xRaw <= xDomain(2) - queryGuardMm;
    if strcmpi(cfg.domainSelectionMode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
    if nnz(mask) < 8
        mask = maskEffective & maskDomain;
        if strcmpi(cfg.domainSelectionMode, 'hard')
            safeMask = mask & maskQuerySafe;
            if nnz(safeMask) >= 8
                mask = safeMask;
            end
        end
    end
    if nnz(mask) < 8
        mask = maskDomain;
        if strcmpi(cfg.domainSelectionMode, 'hard')
            safeMask = mask & maskQuerySafe;
            if nnz(safeMask) >= 8
                mask = safeMask;
            end
        end
    end
    if nnz(mask) < 8
        continue;
    end

    X = [X; xRaw(mask)]; %#ok<AGROW>
    T = [T; tRaw(mask)]; %#ok<AGROW>
    V = [V; vRawMv(mask)]; %#ok<AGROW>
    W = [W; normalize_weight_local(wRaw(mask), cfg.weightFloor)]; %#ok<AGROW>
    Theta = [Theta; thetaRaw(mask)]; %#ok<AGROW>
    S = [S; repmat(sid, nnz(mask), 1)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, nnz(mask), 1)]; %#ok<AGROW>
    F0 = [F0; f0Raw(mask)]; %#ok<AGROW>
    Fx = [Fx; fxRaw(mask)]; %#ok<AGROW>
    g0BySensor(is) = corr.g0Mm;
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = xDomain;
    sensorInfo(is).pointCount = nnz(mask);
    sensorInfo(is).queryGuardMm = queryGuardMm;
    sensorInfo(is).querySafeCount = nnz(maskQuerySafe);
end
if isempty(T)
    error('No valid points were constructed for window %d.', Wmap.window_id);
end
bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.TRel = T - min(T);
bundle.V = V;
bundle.W = W;
bundle.Theta = Theta;
bundle.S = S;
bundle.sensorIndex = sensorIndex;
bundle.sensorIds = cfg.analysisSensors(:).';
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.Template = Template;
bundle.CorrectedGapLibrary = CorrectedGapLibrary;
bundle.responseSurface = responseSurface;
bundle.g0BySensor = g0BySensor;
bundle.sensorInfo = sensorInfo;
bundle.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
bundle.rotRpmMean = Wmap.rot_rpm_mean;
bundle.windowId = Wmap.window_id;
bundle.lapRange = Wmap.lap_range;
bundle.pointCount = numel(T);
bundle.queryGuardMm = max([sensorInfo.queryGuardMm], [], 'omitnan');
end

function bundle = decimate_bundle_local(bundle, maxPoints)
if numel(bundle.T) <= maxPoints
    return;
end
idx = unique(round(linspace(1, numel(bundle.T), maxPoints)));
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    bundle.(f) = bundle.(f)(idx);
end
bundle.pointCount = numel(idx);
end

function seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, ...
    'linearRmseMv', inf, 'fallbackRmseMv', inf), numel(eoCandidates), 1);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo .* bundle.Theta);
    c1 = cos(eo .* bundle.Theta);
    basis = [-bundle.Fx(:), -bundle.Fx(:) .* s1, -bundle.Fx(:) .* c1];
    y = bundle.V(:) - bundle.F0(:);
    sw = sqrt(max(bundle.W(:), cfg.weightFloor));
    coeff = (basis .* sw) \ (y .* sw);
    dx = coeff(1);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y - basis * coeff;
    fixedEval = evaluate_static_warp_waveform_local([A, phi, dx], eo, bundle, cfg, 'fixed');
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx = dx;
    rows(i).linearRmseMv = sqrt(sum(bundle.W(:) .* linRes.^2) / max(numel(linRes), 1));
    rows(i).fallbackRmseMv = fixedEval.weightedRmseMv;
end
[~, order] = sortrows([[rows.fallbackRmseMv].', [rows.linearRmseMv].', [rows.EO].'], [1, 2, 3]);
seedTable = struct2table(rows(order));
end

function eoKeep = select_topk_eo_local(seedTable, topK)
topK = min(max(1, floor(topK)), height(seedTable));
eoKeep = unique(seedTable.EO(1:topK).', 'stable');
end

function eoKeep = select_eo_candidates_local(seedTable, topK, directTrend, windowId, mode)
vpKeep = select_topk_eo_local(seedTable, topK);
directEO = NaN;
if ~isempty(directTrend) && height(directTrend) >= windowId && ismember('EO_id', directTrend.Properties.VariableNames)
    directEO = directTrend.EO_id(windowId);
end
if strcmpi(mode, 'vp') || ~isfinite(directEO)
    eoKeep = vpKeep;
elseif strcmpi(mode, 'direct')
    eoKeep = unique(round(directEO));
else
    eoKeep = unique([round(directEO), vpKeep], 'stable');
end
end

function fit = refine_static_warp_fit_local(bundle, seedTable, eoCandidates, cfg, modeName)
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'dxReferenceMm', NaN, 'deltaDxMm', NaN, ...
    'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'deltaGapMm', [], 'deltaMuGapPerXMm', [], 'deltaTauMm', [], ...
    'VPred', [], 'uMm', [], 'score', inf), numel(eoCandidates), 1);
best = struct('score', inf);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    theta0 = initial_theta_local(seed.A, seed.phi, seed.dx, cfg, modeName, numel(bundle.sensorIds));
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfg, modeName);
    thetaOpt = fminsearch(fun, theta0, cfg.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfg, modeName, numel(bundle.sensorIds));
    detail = evaluate_static_warp_waveform_local(thetaOpt, eo, bundle, cfg, modeName);
    candidate(i).EO = eo;
    candidate(i).freqHz = eo * bundle.rotFreqMeanHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).dxReferenceMm = detail.dxReferenceMm;
    candidate(i).deltaDxMm = detail.deltaDxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    candidate(i).deltaGapMm = detail.deltaGapMm;
    candidate(i).deltaMuGapPerXMm = detail.deltaMuGapPerXMm;
    candidate(i).deltaTauMm = detail.deltaTauMm;
    candidate(i).VPred = detail.VPred;
    candidate(i).uMm = detail.uMm;
    candidate(i).score = detail.score;
    if detail.score < best.score
        best = candidate(i);
    end
end
fit = best;
fit.mode = modeName;
fit.CandidateTable = sortrows(struct2table(rmfield(candidate, {'VPred','uMm','deltaGapMm','deltaMuGapPerXMm','deltaTauMm'})), ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
end

function score = bounded_objective_local(thetaRaw, eo, bundle, cfg, modeName)
theta = bound_params_local(thetaRaw, cfg, modeName, numel(bundle.sensorIds));
detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName);
score = detail.score + 1e4 * sum((thetaRaw(:) - theta(:)).^2);
end

function theta0 = initial_theta_local(A, phi, dx, cfg, modeName, nSensor)
if use_dx_reference_local(cfg, modeName)
    dxParam = dx - cfg.dxReferenceMm;
else
    dxParam = dx;
end
theta0 = [A, phi, dxParam];
if mode_has_gap_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
if mode_has_tilt_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
if mode_has_shift_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
end

function theta = bound_params_local(thetaRaw, cfg, modeName, nSensor)
theta = thetaRaw(:).';
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
if use_dx_reference_local(cfg, modeName)
    theta(3) = max(min(theta(3), cfg.deltaDxLimitMm), -cfg.deltaDxLimitMm);
else
    theta(3) = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
end
need = 3 + mode_extra_count_local(modeName) * nSensor;
if numel(theta) < need
    theta(numel(theta)+1:need) = 0;
end
theta = theta(1:need);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
if ~isempty(idxDg)
    theta(idxDg) = max(min(theta(idxDg), cfg.deltaGapLimitMm), -cfg.deltaGapLimitMm);
end
if ~isempty(idxDmu)
    theta(idxDmu) = max(min(theta(idxDmu), cfg.deltaMuLimit), -cfg.deltaMuLimit);
end
if ~isempty(idxDtau)
    theta(idxDtau) = max(min(theta(idxDtau), cfg.deltaTauLimitMm), -cfg.deltaTauLimitMm);
end
end

function detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName)
theta = bound_params_local(theta, cfg, modeName, numel(bundle.sensorIds));
A = theta(1);
phi = theta(2);
if use_dx_reference_local(cfg, modeName)
    dxReference = cfg.dxReferenceMm;
    deltaDx = theta(3);
else
    dxReference = 0;
    deltaDx = theta(3);
end
dx = dxReference + deltaDx;
nSensor = numel(bundle.sensorIds);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
dg = zeros(nSensor, 1);
dmu = zeros(nSensor, 1);
dtau = zeros(nSensor, 1);
if ~isempty(idxDg), dg = theta(idxDg).'; end
if ~isempty(idxDmu), dmu = theta(idxDmu).'; end
if ~isempty(idxDtau), dtau = theta(idxDtau).'; end
u = A .* sin(eo .* bundle.Theta + phi);
vPred = nan(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensorIds)
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex == is;
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    xEval = bundle.X(mask) - dx - u(mask);
    [vPred(mask), overshoot(mask)] = eval_static_warp_template_local(Tpl, ...
        bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval);
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
pointCount = max(numel(bundle.V), 1);
if nnz(valid) < 8
    residualObj = inf;
    plainRmse = inf;
    overshootPenalty = inf;
else
    res = bundle.V(valid) - vPred(valid);
    ww = max(bundle.W(valid), cfg.weightFloor);
    residualObj = sum(ww .* res.^2);
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
    overshootPenalty = cfg.overshootPenaltyWeight * 1e6 * sum(ww .* overshoot(valid).^2);
end
regularizationPenalty = pointCount * (cfg.staticRegWeightMv * sqrt(( ...
    sum((dg ./ max(cfg.deltaGapLimitMm, eps)).^2) + ...
    sum((dmu ./ max(cfg.deltaMuLimit, eps)).^2) + ...
    sum((dtau ./ max(cfg.deltaTauLimitMm, eps)).^2)) / max(nSensor, 1))) ^ 2;
invalidPenalty = 1e12 * nnz(~valid);
objective = residualObj + overshootPenalty + regularizationPenalty + invalidPenalty;
weightedRmse = sqrt(objective / pointCount);
detail = struct();
detail.EO = eo;
detail.freqHz = eo * bundle.rotFreqMeanHz;
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx;
detail.dxReferenceMm = dxReference;
detail.deltaDxMm = deltaDx;
detail.deltaGapMm = dg(:);
detail.deltaMuGapPerXMm = dmu(:);
detail.deltaTauMm = dtau(:);
detail.weightedRmseMv = weightedRmse;
detail.plainRmseMv = plainRmse;
detail.residualObjectiveMv2 = residualObj;
detail.overshootPenalty = overshootPenalty;
detail.staticRegPenalty = regularizationPenalty;
detail.score = objective;
detail.VPred = vPred;
detail.uMm = u;
end

function cfgOut = configure_dx_reference_local(cfg, seedTable, eoCandidates, modeName)
cfgOut = cfg;
cfgOut.dxReferenceMm = 0;
if ~use_dx_reference_local(cfgOut, modeName)
    return;
end
if strcmpi(cfgOut.dxReferenceMode, 'causal_vp_state')
    if isfield(cfgOut, 'CausalDxReference') && isfield(cfgOut.CausalDxReference, 'dxRef') && ...
            isfinite(cfgOut.CausalDxReference.dxRef)
        cfgOut.dxReferenceMm = max(min(cfgOut.CausalDxReference.dxRef, cfgOut.dxLimitMm), -cfgOut.dxLimitMm);
    end
    return;
end
if strcmpi(cfgOut.dxReferenceMode, 'seed_best')
    if ~isempty(seedTable) && height(seedTable) >= 1 && ismember('dx', seedTable.Properties.VariableNames) && isfinite(seedTable.dx(1))
        cfgOut.dxReferenceMm = max(min(seedTable.dx(1), cfgOut.dxLimitMm), -cfgOut.dxLimitMm);
    end
    return;
end
eoCandidates = unique(round(eoCandidates(:).'));
dxSeed = nan(numel(eoCandidates), 1);
for ie = 1:numel(eoCandidates)
    idx = find(seedTable.EO == eoCandidates(ie), 1, 'first');
    if ~isempty(idx)
        dxSeed(ie) = seedTable.dx(idx);
    end
end
dxSeed = dxSeed(isfinite(dxSeed));
if ~isempty(dxSeed)
    cfgOut.dxReferenceMm = max(min(median(dxSeed), cfgOut.dxLimitMm), -cfgOut.dxLimitMm);
end
end

function cfgDx = default_causal_dx_reference_cfg_local(cfg)
cfgDx = struct();
cfgDx.K = cfg.vpTopK;
cfgDx.maxHypotheses = 6;
cfgDx.NConfirm = parse_positive_integer_env_local('STEP07J_DXREF_N_CONFIRM', 3);
cfgDx.RConfirm = parse_nonnegative_numeric_env_local('STEP07J_DXREF_R_CONFIRM', 5.0);
cfgDx.marginAbs = parse_nonnegative_numeric_env_local('STEP07J_DXREF_MARGIN_ABS', 0.20);
cfgDx.marginFrac = parse_nonnegative_numeric_env_local('STEP07J_DXREF_MARGIN_FRAC', 0.20);
cfgDx.sxMm = parse_nonnegative_numeric_env_local('STEP07J_DXREF_SX_MM', 0.04);
cfgDx.scMm = parse_nonnegative_numeric_env_local('STEP07J_DXREF_SC_MM', 0.08);
cfgDx.srMv = parse_nonnegative_numeric_env_local('STEP07J_DXREF_SR_MV', 12);
cfgDx.dxHardGateMm = parse_nonnegative_numeric_env_local('STEP07J_DXREF_HARD_GATE_MM', 0.08);
cfgDx.wRank = parse_nonnegative_numeric_env_local('STEP07J_DXREF_W_RANK', 0.30);
cfgDx.wRmse = parse_nonnegative_numeric_env_local('STEP07J_DXREF_W_RMSE', 0.70);
cfgDx.wDx = parse_nonnegative_numeric_env_local('STEP07J_DXREF_W_DX', 1.00);
cfgDx.wC = parse_nonnegative_numeric_env_local('STEP07J_DXREF_W_C', 0.25);
cfgDx.alphaDxCold = parse_nonnegative_numeric_env_local('STEP07J_DXREF_ALPHA_DX_COLD', 0.35);
cfgDx.alphaDxTrack = parse_nonnegative_numeric_env_local('STEP07J_DXREF_ALPHA_DX_TRACK', 0.25);
cfgDx.alphaCCold = parse_nonnegative_numeric_env_local('STEP07J_DXREF_ALPHA_C_COLD', 0.30);
cfgDx.alphaCTrack = parse_nonnegative_numeric_env_local('STEP07J_DXREF_ALPHA_C_TRACK', 0.20);
cfgDx.badCountLimit = parse_positive_integer_env_local('STEP07J_DXREF_BAD_COUNT_LIMIT', 2);
cfgDx.confRise = parse_nonnegative_numeric_env_local('STEP07J_DXREF_CONF_RISE', 0.15);
cfgDx.confDrop = parse_nonnegative_numeric_env_local('STEP07J_DXREF_CONF_DROP', 0.30);
end

function tracker = init_causal_dx_reference_tracker_local(cfgDx)
tracker = struct();
tracker.mode = 'cold_start';
tracker.hypotheses = empty_dx_reference_hypothesis_local();
tracker.dxRef = NaN;
tracker.cRef = [NaN, NaN];
tracker.confidence = 0;
tracker.badCount = 0;
tracker.cfg = cfgDx;
end

function h = empty_dx_reference_hypothesis_local()
h = struct('rootEO', {}, 'lastEO', {}, 'dxRef', {}, 'cRef', {}, ...
    'cost', {}, 'winCount', {}, 'Rcum', {}, 'maxDxInnovation', {});
end

function step = empty_causal_dx_reference_step_local()
step = struct('mode', 'cold_start', 'confirmedNow', false, 'resetNow', false, ...
    'dxRef', NaN, 'selectedRank', NaN, 'selectedEO', NaN, 'selectedDx', NaN, ...
    'selectedA', NaN, 'selectedFallbackRmse', NaN, 'dxInnovation', NaN, ...
    'bestCost', NaN, 'secondCost', NaN, 'costMargin', NaN, ...
    'winCount', 0, 'Rcum', 0, 'badCount', 0);
end

function [tracker, step] = update_causal_dx_reference_tracker_local(tracker, seedTable, bundle, cfgDx)
step = empty_causal_dx_reference_step_local();
seeds = extract_causal_dx_seed_rows_local(seedTable, cfgDx.K);
if isempty(seeds) || height(seeds) == 0
    step.mode = tracker.mode;
    step.dxRef = tracker.dxRef;
    return;
end
R_eff = estimate_effective_revolution_local(bundle);
if strcmp(tracker.mode, 'cold_start')
    [tracker, step] = update_causal_dx_cold_start_local(tracker, seeds, R_eff, cfgDx);
else
    [tracker, step] = update_causal_dx_confirmed_local(tracker, seeds, R_eff, cfgDx);
end
end

function seeds = extract_causal_dx_seed_rows_local(seedTable, K)
n = min(K, height(seedTable));
seeds = table();
if n < 1
    return;
end
seeds.rank = (1:n).';
seeds.EO = seedTable.EO(1:n);
seeds.A = seedTable.A(1:n);
seeds.phi = seedTable.phi(1:n);
seeds.dx = seedTable.dx(1:n);
if ismember('fallbackRmseMv', seedTable.Properties.VariableNames)
    seeds.fallbackRmse = seedTable.fallbackRmseMv(1:n);
else
    seeds.fallbackRmse = seedTable.linearRmseMv(1:n);
end
seeds.c1 = seeds.A .* cos(seeds.phi);
seeds.c2 = seeds.A .* sin(seeds.phi);
end

function R_eff = estimate_effective_revolution_local(bundle)
theta = bundle.Theta(:);
theta = theta(isfinite(theta));
if isempty(theta)
    R_eff = 0;
else
    R_eff = max(0, (max(theta) - min(theta)) / (2*pi));
end
end

function [tracker, step] = update_causal_dx_cold_start_local(tracker, seeds, R_eff, cfgDx)
step = empty_causal_dx_reference_step_local();
if isempty(tracker.hypotheses)
    newHyps = initialize_causal_dx_hypotheses_local(seeds, R_eff, cfgDx);
else
    newHyps = expand_causal_dx_hypotheses_local(tracker.hypotheses, seeds, R_eff, cfgDx, ...
        cfgDx.alphaDxCold, cfgDx.alphaCCold);
end
newHyps = sort_and_trim_causal_dx_hypotheses_local(newHyps, cfgDx.maxHypotheses);
tracker.hypotheses = newHyps;
[best, second] = best_two_causal_dx_hypotheses_local(newHyps);
margin = second.cost - best.cost;
isConfirmed = best.winCount >= cfgDx.NConfirm && best.Rcum >= cfgDx.RConfirm && ...
    margin >= max(cfgDx.marginAbs, cfgDx.marginFrac * max(abs(best.cost), eps)) && ...
    best.maxDxInnovation <= cfgDx.dxHardGateMm;
tracker.dxRef = best.dxRef;
tracker.cRef = best.cRef;
if isConfirmed
    tracker.mode = 'confirmed';
    tracker.confidence = 0.6;
end
step.mode = tracker.mode;
step.confirmedNow = isConfirmed;
step.dxRef = best.dxRef;
step.selectedEO = best.lastEO;
step.selectedA = hypot(best.cRef(1), best.cRef(2));
step.dxInnovation = best.maxDxInnovation;
step.bestCost = best.cost;
step.secondCost = second.cost;
step.costMargin = margin;
step.winCount = best.winCount;
step.Rcum = best.Rcum;
step.badCount = tracker.badCount;
end

function [tracker, step] = update_causal_dx_confirmed_local(tracker, seeds, R_eff, cfgDx)
step = empty_causal_dx_reference_step_local();
[idx, quality] = select_causal_dx_seed_measurement_local(seeds, tracker, cfgDx);
seed = seeds(idx, :);
dxInnovation = seed.dx - tracker.dxRef;
if quality.good
    tracker.dxRef = (1 - cfgDx.alphaDxTrack) * tracker.dxRef + cfgDx.alphaDxTrack * seed.dx;
    cSeed = [seed.c1, seed.c2];
    tracker.cRef = (1 - cfgDx.alphaCTrack) * tracker.cRef + cfgDx.alphaCTrack * cSeed;
    tracker.confidence = min(tracker.confidence + cfgDx.confRise, 1);
    tracker.badCount = max(tracker.badCount - 1, 0);
else
    tracker.confidence = max(tracker.confidence - cfgDx.confDrop, 0);
    tracker.badCount = tracker.badCount + 1;
end
resetNow = tracker.badCount >= cfgDx.badCountLimit;
if resetNow
    tracker.mode = 'cold_start';
    tracker.hypotheses = initialize_causal_dx_hypotheses_local(seeds, R_eff, cfgDx);
    tracker.badCount = 0;
end
step.mode = tracker.mode;
step.resetNow = resetNow;
step.dxRef = tracker.dxRef;
step.selectedRank = seed.rank;
step.selectedEO = seed.EO;
step.selectedDx = seed.dx;
step.selectedA = seed.A;
step.selectedFallbackRmse = seed.fallbackRmse;
step.dxInnovation = dxInnovation;
step.bestCost = quality.bestCost;
step.secondCost = quality.secondCost;
step.costMargin = quality.secondCost - quality.bestCost;
step.winCount = NaN;
step.Rcum = NaN;
step.badCount = tracker.badCount;
end

function hyps = initialize_causal_dx_hypotheses_local(seeds, R_eff, cfgDx)
hyps = empty_dx_reference_hypothesis_local();
bestRmse = min(seeds.fallbackRmse);
for i = 1:height(seeds)
    seed = seeds(i, :);
    h = struct();
    h.rootEO = seed.EO;
    h.lastEO = seed.EO;
    h.dxRef = seed.dx;
    h.cRef = [seed.c1, seed.c2];
    h.cost = causal_dx_seed_cost_local(seed, bestRmse, 0, 0, cfgDx);
    h.winCount = 1;
    h.Rcum = R_eff;
    h.maxDxInnovation = 0;
    hyps(end+1) = h; %#ok<AGROW>
end
end

function newHyps = expand_causal_dx_hypotheses_local(hyps, seeds, R_eff, cfgDx, alphaDx, alphaC)
newHyps = empty_dx_reference_hypothesis_local();
bestRmse = min(seeds.fallbackRmse);
for ih = 1:numel(hyps)
    h0 = hyps(ih);
    for is = 1:height(seeds)
        seed = seeds(is, :);
        dxInnovation = seed.dx - h0.dxRef;
        cInnovation = causal_dx_c_innovation_local(seed, h0);
        stepCost = causal_dx_seed_cost_local(seed, bestRmse, dxInnovation, cInnovation, cfgDx);
        h = h0;
        h.lastEO = seed.EO;
        h.dxRef = (1 - alphaDx) * h0.dxRef + alphaDx * seed.dx;
        h.cRef = (1 - alphaC) * h0.cRef + alphaC * [seed.c1, seed.c2];
        h.cost = h0.cost + stepCost;
        h.winCount = h0.winCount + 1;
        h.Rcum = h0.Rcum + R_eff;
        h.maxDxInnovation = max(h0.maxDxInnovation, abs(dxInnovation));
        newHyps(end+1) = h; %#ok<AGROW>
    end
end
end

function [idx, quality] = select_causal_dx_seed_measurement_local(seeds, tracker, cfgDx)
bestRmse = min(seeds.fallbackRmse);
cost = inf(height(seeds), 1);
for i = 1:height(seeds)
    seed = seeds(i, :);
    dxInnovation = seed.dx - tracker.dxRef;
    cInnovation = causal_dx_c_innovation_local(seed, struct('lastEO', seed.EO, 'cRef', tracker.cRef));
    cost(i) = causal_dx_seed_cost_local(seed, bestRmse, dxInnovation, cInnovation, cfgDx);
end
[costSorted, order] = sort(cost, 'ascend');
idx = order(1);
quality = struct();
quality.bestCost = costSorted(1);
if numel(costSorted) >= 2
    quality.secondCost = costSorted(2);
else
    quality.secondCost = inf;
end
quality.good = abs(seeds.dx(idx) - tracker.dxRef) <= cfgDx.dxHardGateMm;
end

function cost = causal_dx_seed_cost_local(seed, bestRmse, dxInnovation, cInnovation, cfgDx)
rankPenalty = seed.rank - 1;
rmsePenalty = max(0, seed.fallbackRmse - bestRmse) / max(cfgDx.srMv, eps);
dxPenalty = causal_dx_huber_local(dxInnovation / max(cfgDx.sxMm, eps));
cPenalty = causal_dx_huber_local(cInnovation / max(cfgDx.scMm, eps));
cost = cfgDx.wRank * rankPenalty + cfgDx.wRmse * rmsePenalty + ...
    cfgDx.wDx * dxPenalty + cfgDx.wC * cPenalty;
end

function cInnovation = causal_dx_c_innovation_local(seed, h)
if isfield(h, 'lastEO') && seed.EO == h.lastEO && all(isfinite(h.cRef))
    cInnovation = norm([seed.c1, seed.c2] - h.cRef);
else
    cInnovation = 0;
end
end

function v = causal_dx_huber_local(x)
a = abs(x);
if a <= 1
    v = 0.5 * a.^2;
else
    v = a - 0.5;
end
end

function hyps = sort_and_trim_causal_dx_hypotheses_local(hyps, maxN)
if isempty(hyps)
    return;
end
[~, order] = sort([hyps.cost], 'ascend');
hyps = hyps(order(1:min(maxN, numel(order))));
end

function [best, second] = best_two_causal_dx_hypotheses_local(hyps)
hyps = sort_and_trim_causal_dx_hypotheses_local(hyps, numel(hyps));
best = hyps(1);
if numel(hyps) >= 2
    second = hyps(2);
else
    second = best;
    second.cost = inf;
end
end

function tf = use_dx_reference_local(cfg, modeName)
tf = strcmpi(modeName, cfg.mainModel) && isfield(cfg, 'deltaDxLimitMm') && ...
    isfinite(cfg.deltaDxLimitMm) && cfg.deltaDxLimitMm > 0;
end

function n = mode_extra_count_local(modeName)
n = double(mode_has_gap_local(modeName)) + double(mode_has_tilt_local(modeName)) + ...
    double(mode_has_shift_local(modeName));
end

function tf = mode_has_gap_local(modeName)
tf = any(strcmpi(modeName, {'gap_only','gap_tilt','gap_tilt_shift'}));
end

function tf = mode_has_tilt_local(modeName)
tf = any(strcmpi(modeName, {'gap_tilt','gap_tilt_shift'}));
end

function tf = mode_has_shift_local(modeName)
tf = strcmpi(modeName, 'gap_tilt_shift');
end

function [idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor)
p = 3;
idxDg = [];
idxDmu = [];
idxDtau = [];
if mode_has_gap_local(modeName)
    idxDg = p + (1:nSensor);
    p = p + nSensor;
end
if mode_has_tilt_local(modeName)
    idxDmu = p + (1:nSensor);
    p = p + nSensor;
end
if mode_has_shift_local(modeName)
    idxDtau = p + (1:nSensor);
end
end

function [model, overshoot] = eval_static_warp_template_local(Tpl, responseSurface, corr, dg, dmu, dtau, xOpr)
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function [model, overshoot] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr)
xLowRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xLowRaw, xLo), xHi);
overshootLow = max(xLo - xLowRaw, 0) + max(xLowRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function dMdx = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, g, xOpr, h)
[mp, ~] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr + h);
[mm, ~] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr - h);
dMdx = (mp - mm) ./ (2 * h);
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr(:) - tau);
gEffRaw = g0 + mu .* (xOpr(:) - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
xLib = min(max(xLibRaw, xMin), xMax);
gEff = min(max(gEffRaw, gMin), gMax);
overshoot = hypot(max(xMin - xLibRaw, 0) + max(xLibRaw - xMax, 0), ...
    max(gMin - gEffRaw, 0) + max(gEffRaw - gMax, 0));
model = eval_response_surface_local(responseSurface, gEff, xLib);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1);
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function threshold = resolve_threshold_local(Tpl, cfg)
if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = cfg.defaultSensorThreshold;
end
end

function mask = build_pulse_mask_local(vMv, thresholdMv, modeName)
if strcmpi(modeName, 'single')
    mask = isolate_main_pulse_local(vMv, thresholdMv);
else
    mask = isolate_all_pulses_local(vMv, thresholdMv);
end
end

function [mask, segmentCount] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segmentCount = numel(starts);
end

function [mask, segmentCount] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segmentCount = 1;
    return;
end

segmentGap = idx(starts(2:end)) - idx(ends(1:end-1));
largeGapThreshold = max(10, round(0.02 * numel(v)));
pulseBreaks = find(segmentGap > largeGapThreshold);
groupStarts = [1; pulseBreaks + 1];
groupEnds = [pulseBreaks; numel(starts)];
for ig = 1:numel(groupStarts)
    segIds = groupStarts(ig):groupEnds(ig);
    bestSeg = segIds(1);
    bestPeak = -inf;
    for iseg = segIds
        seg = idx(starts(iseg):ends(iseg));
        peakVal = max(v(seg));
        if peakVal > bestPeak
            bestPeak = peakVal;
            bestSeg = iseg;
        end
    end
    mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end
segmentCount = numel(groupStarts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, cfg)
finite = isfinite(x) & isfinite(v) & isfinite(t);
gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gMax = max(gTpl(finite), [], 'omitnan');
if isfinite(gMax) && gMax > 0
    maskTpl = gTpl >= cfg.dynamicGradientMinRatio * gMax;
else
    maskTpl = finite;
end

gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if isfinite(gTimeMax) && gTimeMax > 0
    maskTime = gTime >= cfg.dynamicTimeGradientMinRatio * gTimeMax;
else
    maskTime = false(size(v));
end

peakQ = min(max(cfg.dynamicPeakQuantile, 0), 100);
peakLevel = prctile(v(finite), peakQ);
maskPeak = v >= max(threshold, peakLevel);

mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function queryGuardMm = resolve_query_guard_mm_local(Tpl, x, v, baseMask, cfg)
queryGuardMm = cfg.queryGuardFixedMm;
if ~strcmpi(cfg.queryGuardMode, 'adaptive') || nnz(baseMask) < 8
    return;
end
xSel = x(baseMask);
vSel = v(baseMask);
xStat = invert_template_voltage_local(Tpl, vSel, xSel);
uApp = abs(xSel(:) - xStat(:));
uApp = uApp(isfinite(uApp));
if isempty(uApp)
    return;
end
q = min(max(cfg.queryGuardQuantile, 0), 100);
adaptiveGuard = prctile(uApp, q) + cfg.queryGuardSafetyMm;
queryGuardMm = min(max(adaptiveGuard, cfg.queryGuardMinMm), cfg.queryGuardMaxMm);
end

function xStat = invert_template_voltage_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diffV = vGrid - vv;
    crossingX = [];
    exactIdx = find(abs(diffV) <= 1e-10);
    if ~isempty(exactIdx)
        crossingX = xGrid(exactIdx);
    end
    for k = 1:numel(diffV)-1
        if ~isfinite(diffV(k)) || ~isfinite(diffV(k+1))
            continue;
        end
        if diffV(k) == 0 || diffV(k) * diffV(k+1) > 0
            continue;
        end
        denom = vGrid(k+1) - vGrid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - vGrid(k)) / denom;
        crossingX(end+1, 1) = xGrid(k) + alpha * (xGrid(k+1) - xGrid(k)); %#ok<AGROW>
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function w = normalize_weight_local(wRaw, floorW)
w = max(wRaw(:), floorW);
if max(w) > 0
    w = max(floorW, w ./ max(w));
end
end

function eoCandidates = build_eo_candidates_local(rotFreqHz, freqSearchHz, eoPad)
rotFreqHz = max(rotFreqHz, eps);
eoMin = max(1, ceil(min(freqSearchHz) / rotFreqHz));
eoMax = max(eoMin, floor(max(freqSearchHz) / rotFreqHz));
eoCandidates = eoMin:eoMax;
if isempty(eoCandidates)
    eoCenter = max(1, round(mean(freqSearchHz) / rotFreqHz));
    eoCandidates = max(1, eoCenter - eoPad):max(1, eoCenter + eoPad);
end
freq = eoCandidates .* rotFreqHz;
eoCandidates = eoCandidates(freq >= min(freqSearchHz) & freq <= max(freqSearchHz));
end

function wr = pack_window_result_local(iw, Wmap, bundle, seedTable, selectedEO, mainCandidateEO, modelFits, dxReferenceStep)
wr = struct();
wr.windowId = iw;
wr.lapRange = Wmap.lap_range;
wr.timeWindow = Wmap.time_window;
wr.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
wr.rotRpmMean = Wmap.rot_rpm_mean;
wr.bundle = bundle;
wr.VPSeedTable = seedTable;
wr.VPSelectedEO = selectedEO;
wr.MainCandidateEO = mainCandidateEO;
wr.modelFits = modelFits;
wr.CausalDxReference = dxReferenceStep;
end

function row = make_trend_row_local(wr, directTrend)
directEO = NaN; directFreq = NaN; directAmp = NaN; directRmseMv = NaN;
if ~isempty(directTrend) && height(directTrend) >= wr.windowId
    directEO = directTrend.EO_id(wr.windowId);
    directFreq = directTrend.fn_id(wr.windowId);
    directAmp = directTrend.A_id(wr.windowId);
    directRmseMv = 1000 * directTrend.weighted_voltage_rmse(wr.windowId);
end
F = wr.modelFits;
fixedFit = get_model_fit_or_nan_local(F, 'fixed');
gapFit = get_model_fit_or_nan_local(F, 'gap_only');
tiltFit = get_model_fit_or_nan_local(F, 'gap_tilt');
shiftFit = get_model_fit_or_nan_local(F, 'gap_tilt_shift');
D = wr.CausalDxReference;
row = table(wr.windowId, wr.lapRange(1), wr.lapRange(end), wr.rotFreqMeanHz, ...
    D.dxRef, string(D.mode), D.confirmedNow, D.resetNow, D.selectedEO, ...
    D.selectedRank, D.selectedDx, D.dxInnovation, D.costMargin, ...
    directEO, directFreq, directAmp, directRmseMv, ...
    fixedFit.EO, fixedFit.freqHz, fixedFit.amplitudeMm, fixedFit.dxMm, fixedFit.weightedRmseMv, ...
    gapFit.EO, gapFit.freqHz, gapFit.amplitudeMm, gapFit.dxMm, gapFit.weightedRmseMv, ...
    mean(gapFit.deltaGapMm, 'omitnan'), ...
    tiltFit.EO, tiltFit.freqHz, tiltFit.amplitudeMm, tiltFit.dxMm, tiltFit.weightedRmseMv, ...
    mean(tiltFit.deltaGapMm, 'omitnan'), mean(tiltFit.deltaMuGapPerXMm, 'omitnan'), ...
    shiftFit.EO, shiftFit.freqHz, shiftFit.amplitudeMm, shiftFit.dxMm, shiftFit.weightedRmseMv, ...
    mean(shiftFit.deltaGapMm, 'omitnan'), mean(shiftFit.deltaMuGapPerXMm, 'omitnan'), ...
    mean(shiftFit.deltaTauMm, 'omitnan'), ...
    'VariableNames', {'window_id','lap_start','lap_end','rot_freq_hz', ...
    'dxref_mm','dxref_mode','dxref_confirmed_now','dxref_reset_now','dxref_selected_EO','dxref_selected_rank','dxref_selected_dx_mm','dxref_innovation_mm','dxref_cost_margin', ...
    'direct_EO','direct_frequency_hz','direct_amplitude_mm','direct_rmse_mV', ...
    'fixed_EO','fixed_frequency_hz','fixed_amplitude_mm','fixed_dx_mm','fixed_rmse_mV', ...
    'gap_EO','gap_frequency_hz','gap_amplitude_mm','gap_dx_mm','gap_rmse_mV','gap_mean_delta_gap_mm', ...
    'tilt_EO','tilt_frequency_hz','tilt_amplitude_mm','tilt_dx_mm','tilt_rmse_mV','tilt_mean_delta_gap_mm','tilt_mean_delta_mu', ...
    'shift_EO','shift_frequency_hz','shift_amplitude_mm','shift_dx_mm','shift_rmse_mV','shift_mean_delta_gap_mm','shift_mean_delta_mu','shift_mean_delta_tau_mm'});
end

function fit = get_model_fit_or_nan_local(F, name)
if isfield(F, name)
    fit = F.(name);
    return;
end
fit = struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'dxMm', NaN, 'weightedRmseMv', NaN, 'deltaGapMm', NaN, ...
    'deltaMuGapPerXMm', NaN, 'deltaTauMm', NaN);
end

function Summary = build_summary_table_local(Trend, bestIdx, cfg)
allRows = struct( ...
    'method', {'fixed','gap_only','gap_tilt','gap_tilt_shift'}, ...
    'eoCol', {'fixed_EO','gap_EO','tilt_EO','shift_EO'}, ...
    'freqCol', {'fixed_frequency_hz','gap_frequency_hz','tilt_frequency_hz','shift_frequency_hz'}, ...
    'ampCol', {'fixed_amplitude_mm','gap_amplitude_mm','tilt_amplitude_mm','shift_amplitude_mm'}, ...
    'rmseCol', {'fixed_rmse_mV','gap_rmse_mV','tilt_rmse_mV','shift_rmse_mV'});
keep = ismember({allRows.method}, cfg.modelNames);
rows = allRows(keep);
method = string({rows.method}).';
dominantEO = NaN(numel(rows), 1);
meanFreq = NaN(numel(rows), 1);
stdFreq = NaN(numel(rows), 1);
meanAmp = NaN(numel(rows), 1);
meanRmse = NaN(numel(rows), 1);
medianRmse = NaN(numel(rows), 1);
for i = 1:numel(rows)
    dominantEO(i) = mode_finite_local(Trend.(rows(i).eoCol));
    meanFreq(i) = mean(Trend.(rows(i).freqCol), 'omitnan');
    stdFreq(i) = std(Trend.(rows(i).freqCol), 'omitnan');
    meanAmp(i) = mean(Trend.(rows(i).ampCol), 'omitnan');
    meanRmse(i) = mean(Trend.(rows(i).rmseCol), 'omitnan');
    medianRmse(i) = median(Trend.(rows(i).rmseCol), 'omitnan');
end
bestWindowIndex = repmat(bestIdx, numel(rows), 1);
Summary = table(method, dominantEO, meanFreq, stdFreq, meanAmp, meanRmse, medianRmse, bestWindowIndex, ...
    'VariableNames', {'method','dominant_EO','mean_frequency_hz','std_frequency_hz', ...
    'mean_amplitude_mm','mean_rmse_mV','median_rmse_mV','best_free_gap_window'});
end

function plot_trend_local(T, figFile)
style = paper_style_local();
fig = figure('Name', 'Step07J nested static warp trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 13.0]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_EO, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_EO, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.gap_EO, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Gap');
plot(T.window_id, T.tilt_EO, '-', 'Color', style.green, 'LineWidth', 1.0, ...
    'Marker', '^', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt');
plot(T.window_id, T.shift_EO, '-', 'Color', style.orange, 'LineWidth', 1.0, ...
    'Marker', 'd', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt+shift');
ylabel('EO'); title('Order'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_amplitude_mm, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_amplitude_mm, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.gap_amplitude_mm, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Gap');
plot(T.window_id, T.tilt_amplitude_mm, '-', 'Color', style.green, 'LineWidth', 1.0, ...
    'Marker', '^', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt');
plot(T.window_id, T.shift_amplitude_mm, '-', 'Color', style.orange, 'LineWidth', 1.0, ...
    'Marker', 'd', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt+shift');
ylabel('A (mm)', 'Interpreter', 'tex'); title('Amplitude'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_rmse_mV, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_rmse_mV, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.gap_rmse_mV, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Gap');
plot(T.window_id, T.tilt_rmse_mV, '-', 'Color', style.green, 'LineWidth', 1.0, ...
    'Marker', '^', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt');
plot(T.window_id, T.shift_rmse_mV, '-', 'Color', style.orange, 'LineWidth', 1.0, ...
    'Marker', 'd', 'MarkerSize', 4.2, 'DisplayName', 'Gap+tilt+shift');
xlabel('Window'); ylabel('RMSE (mV)'); title('Voltage residual'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function plot_best_window_collapse_local(wr, cfg, figFile)
style = paper_style_local();
bundle = decimate_bundle_local(wr.bundle, cfg.maxPointsPerWindow);
fit = wr.modelFits.(cfg.mainModel);
fig = figure('Name', 'Step07J best-window dynamic collapse', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 4.8 * numel(cfg.analysisSensors)]);
tiledlayout(numel(cfg.analysisSensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = bundle.sensorIndex == is;
    xRaw = bundle.X(mask) - fit.dxMm;
    xComp = xRaw - fit.uMm(mask);
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    dg = fit.deltaGapMm(is);
    dmu = fit.deltaMuGapPerXMm(is);
    dtau = fit.deltaTauMm(is);
    [modelNo, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, dg, dmu, dtau, xRaw);
    [modelDyn, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, dg, dmu, dtau, xComp);
    [xs1, o1] = sort(xRaw);
    [xs2, o2] = sort(xComp);
    vObs = bundle.V(mask);

    nexttile; hold on; grid on; box on;
    plot(xRaw, vObs, '.', 'Color', style.gray, 'MarkerSize', 4.0, 'DisplayName', 'Samples');
    plot(xs1, modelNo(o1), '-', 'Color', style.red, 'LineWidth', 1.2, 'DisplayName', 'Static template');
    title(sprintf('CH%d before', sid));
    xlabel('x_{OPR} (mm)', 'Interpreter', 'tex'); ylabel('mV'); legend('Location', 'southeast', 'Box', 'off');
    format_axes_local(gca, style);

    nexttile; hold on; grid on; box on;
    plot(xComp, vObs, '.', 'Color', style.gray, 'MarkerSize', 4.0, 'DisplayName', 'Samples');
    plot(xs2, modelDyn(o2), '-', 'Color', style.red, 'LineWidth', 1.2, 'DisplayName', 'Dynamic template');
    title(sprintf('CH%d after, dg %.4f, dmu %.4f, dtau %.4f', sid, dg, dmu, dtau));
    xlabel('x_{OPR}-u(t) (mm)', 'Interpreter', 'tex'); ylabel('mV'); legend('Location', 'southeast', 'Box', 'off');
    format_axes_local(gca, style);
end
export_paper_figure_local(fig, figFile);
end

function plot_static_warp_local(Result, T, cfg, figFile)
style = paper_style_local();
fig = figure('Name', 'Step07J static warp parameters', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 13.0]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sensorIds = cfg.analysisSensors(:).';
[dgMat, dmuMat, dtauMat] = collect_static_warp_matrix_local(Result.WindowResult, sensorIds, cfg.mainModel);
colors = sensor_colors_local(numel(sensorIds), style);

nexttile; hold on; box on;
for is = 1:numel(sensorIds)
    plot(T.window_id, dgMat(:, is), '-', 'Color', colors(is, :), ...
        'LineWidth', 1.1, 'Marker', sensor_marker_local(is), ...
        'MarkerSize', 4.2, 'DisplayName', sprintf('CH%d', sensorIds(is)));
end
yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', 'dg = 0');
ylabel('\Delta g (mm)', 'Interpreter', 'tex');
title('Gap increment, gap-tilt main model');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
for is = 1:numel(sensorIds)
    plot(T.window_id, dmuMat(:, is), '-', 'Color', colors(is, :), ...
        'LineWidth', 1.1, 'Marker', sensor_marker_local(is), ...
        'MarkerSize', 4.2, 'DisplayName', sprintf('CH%d', sensorIds(is)));
end
yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', 'dmu = 0');
xlabel('Window');
ylabel('\Delta \mu', 'Interpreter', 'tex');
title('High-speed tilt increment, gap-tilt main model');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
for is = 1:numel(sensorIds)
    plot(T.window_id, dtauMat(:, is), '-', 'Color', colors(is, :), ...
        'LineWidth', 1.1, 'Marker', sensor_marker_local(is), ...
        'MarkerSize', 4.2, 'DisplayName', sprintf('CH%d', sensorIds(is)));
end
yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', 'dtau = 0');
xlabel('Window');
ylabel('\Delta \tau (mm)', 'Interpreter', 'tex');
title('Sensor x-shift increment (comparison model only)');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function [dgMat, dmuMat, dtauMat] = collect_static_warp_matrix_local(WindowResult, sensorIds, mainModel)
nWin = numel(WindowResult);
nSensor = numel(sensorIds);
dgMat = NaN(nWin, nSensor);
dmuMat = NaN(nWin, nSensor);
dtauMat = NaN(nWin, nSensor);
for iw = 1:nWin
    if ~isfield(WindowResult(iw), 'modelFits') || ...
            ~isfield(WindowResult(iw).modelFits, mainModel)
        continue;
    end
    fit = WindowResult(iw).modelFits.(mainModel);
    dg = fit.deltaGapMm(:).';
    dmu = fit.deltaMuGapPerXMm(:).';
    dtau = fit.deltaTauMm(:).';
    nUse = min(nSensor, numel(dg));
    dgMat(iw, 1:nUse) = dg(1:nUse);
    dmuMat(iw, 1:min(nSensor, numel(dmu))) = dmu(1:min(nSensor, numel(dmu)));
    dtauMat(iw, 1:min(nSensor, numel(dtau))) = dtau(1:min(nSensor, numel(dtau)));
end
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.tickFontSize = 8.2;
style.black = [0.08, 0.08, 0.08];
style.gray = [0.35, 0.35, 0.35];
style.blue = [0.00, 0.28, 0.70];
style.red = [0.82, 0.10, 0.10];
style.green = [0.00, 0.45, 0.28];
style.orange = [0.85, 0.37, 0.05];
end

function colors = sensor_colors_local(n, style)
base = [
    style.blue
    style.red
    style.green
    style.orange
    0.45, 0.20, 0.70
    0.15, 0.55, 0.70
    ];
if n <= size(base, 1)
    colors = base(1:n, :);
else
    colors = lines(n);
end
end

function marker = sensor_marker_local(idx)
markers = {'o', 's', '^', 'd', 'v', '>'};
marker = markers{mod(idx - 1, numel(markers)) + 1};
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on', ...
    'XMinorTick', 'off', 'YMinorTick', 'off');
grid(ax, 'off');
xlabel(ax, get(get(ax, 'XLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
ylabel(ax, get(get(ax, 'YLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
lgd = legend(ax);
if ~isempty(lgd) && isvalid(lgd)
    set(lgd, 'FontName', style.fontName, 'FontSize', style.tickFontSize);
end
end

function export_paper_figure_local(fig, pngFile)
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
pdfFile = fullfile(folder, [name, '.pdf']);
emfFile = fullfile(folder, [name, '.emf']);
try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
catch
    warning('Could not export PDF figure: %s', pdfFile);
end
try
    print(fig, emfFile, '-dmeta', '-r300');
catch
    warning('Could not export EMF figure: %s', emfFile);
end
end

function v = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    v = defaultValue;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = defaultValue;
else
    v = max(1, floor(tmp));
end
end

function v = parse_nonnegative_numeric_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    v = defaultValue;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp < 0
    error('%s must be a nonnegative numeric value.', name);
end
v = tmp;
end

function v = parse_logical_env_local(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    v = logical(defaultValue);
    return;
end
if ismember(raw, {'1', 'true', 'yes', 'on'})
    v = true;
elseif ismember(raw, {'0', 'false', 'no', 'off'})
    v = false;
else
    error('%s must be one of on/off, true/false, yes/no, or 1/0.', name);
end
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function suffix = sanitize_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_' suffix];
end
end

function m = mode_finite_local(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mode(x);
end
end

