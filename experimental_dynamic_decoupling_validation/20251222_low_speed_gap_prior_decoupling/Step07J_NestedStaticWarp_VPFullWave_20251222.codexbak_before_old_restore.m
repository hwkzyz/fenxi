%% Step07J: Nested high-speed static warp models
% The main gap_tilt route jointly identifies EO, vibration, static gap
% increment and a small high-speed tilt increment inside the direct-template
% + VP Top-K EO candidate set with the full waveform objective. VP narrows
% the EO candidates but does not fix the final EO. The gap-increment search
% half-width is set by Delta g_proj = epsilon_perp / S_g and
% Delta g_loc = alpha * Delta g_proj; STEP07J_DELTA_G_LIMIT_MM is only the
% projection cap/fallback unless STEP07J_USE_DELTA_G_PROJ=0 is used for an
% ablation run. The default main route computes only the paper-facing
% gap_tilt model. In comparison mode this script also computes nested static
% correction models on top of the low-speed template:
%
%   fixed             : no clearance increment
%   gap_only          : dg_s
%   gap_tilt          : dg_s + dmu_s*(x-tau_s)
%
% Paper-facing main output uses gap_tilt only. gap_only is retained only in
% the separate comparison run. gap_tilt_shift is excluded from the current
% comparison because the extra sensor-shift freedom is invalid for the
% paper-facing method.
%
% This is the current paper-facing identification route.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
C0 = CaseConfig();
[resonanceSelection, ~] = ResonanceRegionCatalog_20251222();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_nested_static_warp_vp_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = C0.bladeId;
cfg.analysisSensors = C0.sensorIds;
cfg.resonanceSelection = resonanceSelection;
cfg.resonanceRegionShortTag = resonanceSelection.shortTag;
cfg.freqSearchHz = [100, 1000];
cfg.eoPad = parse_positive_integer_env_local('STEP07J_EO_PAD', 2);
cfg.eoSearchMode = lower(strtrim(getenv('STEP07J_EO_SEARCH_MODE')));
if isempty(cfg.eoSearchMode)
    cfg.eoSearchMode = 'region_order';
end
if ~ismember(cfg.eoSearchMode, {'region_order', 'frequency_range'})
    error('STEP07J_EO_SEARCH_MODE must be "region_order" or "frequency_range".');
end
cfg.eoSearchCenterOrder = resonanceSelection.dominantOrder;
cfg.eoSearchHalfWidth = floor(parse_nonnegative_numeric_env_local( ...
    'STEP07J_EO_SEARCH_HALF_WIDTH', cfg.eoPad));
[cfg.freqSearchHz, hasFreqOverride] = parse_frequency_range_override_local(cfg.freqSearchHz);
if hasFreqOverride
    cfg.eoSearchMode = 'frequency_range';
end
cfg.vpTopK = parse_positive_integer_env_local('STEP07J_TOP_K_EO', 3);
cfg.maxWindows = parse_positive_integer_env_local('STEP07J_MAX_WINDOWS', inf);
cfg.maxPointsPerWindow = parse_positive_integer_env_local('STEP07J_MAX_POINTS_PER_WINDOW', inf);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP07J_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DX_LIMIT_MM', 0.35);
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_LIMIT_MM', 0.25);
cfg.useDeltaGapProjectionLimit = parse_logical_env_local('STEP07J_USE_DELTA_G_PROJ', true);
cfg.deltaGapProjAlpha = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_ALPHA', 2.5);
cfg.deltaGapProjMinLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_MIN_LIMIT_MM', 0.02);
cfg.deltaGapProjMaxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_MAX_LIMIT_MM', cfg.deltaGapLimitMm);
cfg.deltaGapProjFallbackLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_FALLBACK_LIMIT_MM', cfg.deltaGapLimitMm);
cfg.deltaGapDerivativeStepMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_DERIV_STEP_MM', 1e-4);
cfg.deltaGapProjSensitivityFloorMvPerMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_SG_FLOOR_MV_PER_MM', 1e-6);
cfg.deltaMuLimit = parse_nonnegative_numeric_env_local('STEP07J_DELTA_MU_LIMIT', 0.030);
cfg.deltaTauLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_TAU_LIMIT_MM', 0.16);
cfg.etaFallbackLimitMm = parse_nonnegative_numeric_env_local('STEP07J_ETA_FALLBACK_LIMIT_MM', 0.06);
cfg.etaRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_ETA_REG_WEIGHT_MV', 0.35);
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
    cfg.modelNames = {'fixed','gap_only','gap_tilt'};
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
cfg.prevWindowCandidateMode = lower(strtrim(getenv('STEP07J_PREV_WINDOW_CANDIDATE_MODE')));
if isempty(cfg.prevWindowCandidateMode)
    cfg.prevWindowCandidateMode = 'soft';
end
if ~ismember(cfg.prevWindowCandidateMode, {'off', 'soft'})
    error('STEP07J_PREV_WINDOW_CANDIDATE_MODE must be "off" or "soft".');
end
cfg.prevWindowRefreshEvery = floor(parse_nonnegative_numeric_env_local( ...
    'STEP07J_PREV_WINDOW_REFRESH_EVERY', 5));
cfg.phaseSafeExpansion = parse_logical_env_local('STEP07J_PHASE_SAFE_EXPANSION', true);
cfg.phaseSafeMarginMm = parse_nonnegative_numeric_env_local('STEP07J_PHASE_SAFE_MARGIN_MM', 0.03);
cfg.phaseSafeFallbackMode = lower(strtrim(getenv('STEP07J_PHASE_SAFE_FALLBACK_MODE')));
if isempty(cfg.phaseSafeFallbackMode)
    cfg.phaseSafeFallbackMode = 'linear_vp';
end
if ~ismember(cfg.phaseSafeFallbackMode, {'off', 'linear_vp'})
    error('STEP07J_PHASE_SAFE_FALLBACK_MODE must be "off" or "linear_vp".');
end
cfg.directBundleMode = lower(strtrim(getenv('STEP07J_DIRECT_BUNDLE_MODE')));
if isempty(cfg.directBundleMode)
    cfg.directBundleMode = 'off';
end
if ~ismember(cfg.directBundleMode, {'auto', 'off', 'off_when_phase_safe'})
    error('STEP07J_DIRECT_BUNDLE_MODE must be "auto", "off", or "off_when_phase_safe".');
end
cfg.derivativeStepMm = 1e-3;
cfg.weightFloor = 0.05;
cfg.defaultSensorThreshold = 0.5;
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
        resultSuffix = sprintf('_%s_method_compare', resonanceSelection.shortTag);
    else
        resultSuffix = sprintf('_%s_main_gaptilt', resonanceSelection.shortTag);
    end
end

rotDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
rotDynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');
rotResultDir = fullfile(rotDir, 'output', 'identification');
highMapCacheDir = outDir;

correctedLibOverride = strtrim(getenv('STEP07J_CORRECTED_LIB_FILE'));
if ~isempty(correctedLibOverride)
    correctedLibFile = correctedLibOverride;
    correctedLibMode = 'override';
else
    [CorrectedGapLibrary, correctedLibFile, correctedLibMode] = load_gap_calibration_library_local( ...
        outDir, C0, cfg);
    if isempty(CorrectedGapLibrary)
        correctedLibFile = find_corrected_gap_library_file_local(outDir, cfg.targetBlade, sensorTag);
        correctedLibMode = 'combined_fallback';
    end
end
templateOverride = strtrim(getenv('STEP07J_TEMPLATE_FILE'));
if ~isempty(templateOverride)
    templateFile = templateOverride;
    templateSourceMode = 'override';
else
    addpath(rotDir);
    templateCase = struct();
    templateCase.dataset = '20251222';
    templateCase.bladeId = cfg.targetBlade;
    templateCase.sensorIds = cfg.analysisSensors;
    templateCase.sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
    templateCase.caseTag = sprintf('B%d_%s', cfg.targetBlade, templateCase.sensorTag);
    [templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20251222(rotDir, templateCase, 'GradientXRange030_OPRCenterStd');
end
dynamicOverride = strtrim(getenv('STEP07J_DYNAMIC_MAP_FILE'));
highMapOverride = strtrim(getenv('STEP07J_HIGHMAP_FILE'));
dynamicSourceMode = 'step06_highmap_default';
if ~isempty(dynamicOverride)
    if exist(dynamicOverride, 'file') ~= 2
        error('STEP07J_DYNAMIC_MAP_FILE does not exist: %s', dynamicOverride);
    end
    dynamicFile = dynamicOverride;
    dynamicSourceMode = 'dynamic_map_override';
elseif ~isempty(highMapOverride)
    if exist(highMapOverride, 'file') ~= 2
        error('STEP07J_HIGHMAP_FILE does not exist: %s', highMapOverride);
    end
    dynamicFile = highMapOverride;
    dynamicSourceMode = 'highmap_override';
else
    try
        dynamicFile = ensure_dynamic_map_file_local(outDir, thisDir, C0, cfg);
        dynamicSourceMode = 'step06_dynamic_map_cache';
    catch highMapErr
        warning('Step06 DynamicMap unavailable (%s). Falling back to legacy DynamicMap.', highMapErr.message);
        dynamicFile = find_dynamic_map_file_local(rotDynamicDir, cfg.targetBlade, cfg.analysisSensors, resonanceSelection.shortTag);
        dynamicSourceMode = 'rotating_calibration_dynamic_map';
    end
end
directOverride = strtrim(getenv('STEP07J_DIRECT_RESULT_FILE'));
if ~isempty(directOverride)
    if exist(directOverride, 'file') ~= 2
        error('STEP07J_DIRECT_RESULT_FILE does not exist: %s', directOverride);
    end
    directMainFile = directOverride;
else
    directMainFile = '';
end

if ~isfile(correctedLibFile)
    error('Run Step06I first. Missing file: %s', correctedLibFile);
end

St = load(templateFile, 'Template');
Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
cfg.lowSpeedTemplate = Template;
[DynamicMap, highMapMeta] = load_dynamic_input_local(dynamicFile, dynamicSourceMode, cfg);
DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, cfg.analysisSensors, sensorTag);
coordinateCheck = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysisSensors, 1e-6);
if ~exist('CorrectedGapLibrary', 'var') || isempty(CorrectedGapLibrary)
    Sc = load(correctedLibFile, 'CorrectedGapLibrary');
    CorrectedGapLibrary = Sc.CorrectedGapLibrary;
    correctedLibMode = 'combined';
end
responseSurface = CorrectedGapLibrary.responseSurface;
directRef = load_direct_reference_local(directMainFile);
directTrend = directRef.Trend;

numWindows = min(numel(DynamicMap.Window), cfg.maxWindows);
fprintf('\n=== Step07J: nested static warp VP/full-wave ===\n');
fprintf('Template: %s [%s]\n', templateFile, templateSourceMode);
fprintf('Dynamic input: %s [%s]\n', dynamicFile, dynamicSourceMode);
fprintf('Corrected gap library: %s [%s]\n', correctedLibFile, correctedLibMode);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinateCheck.max_abs_delta_mm, coordinateCheck.status);
if ~coordinateCheck.is_consistent
    warning('Template/DynamicMap x-center mismatch. Do not use this run as the unified GradientXRange030 main result.');
end
if ~isempty(directMainFile)
    fprintf('Direct reference: %s\n', directMainFile);
else
    fprintf('Direct reference: none (local gap-prior route only)\n');
end
fprintf('Windows: %d, sensors: %s, VP Top-K: %d, EO candidate mode: %s\n', ...
    numWindows, mat2str(cfg.analysisSensors), cfg.vpTopK, cfg.eoCandidateMode);
fprintf('EO search: %s, center EO %g, half-width %d, frequency range [%.1f %.1f] Hz\n', ...
    cfg.eoSearchMode, cfg.eoSearchCenterOrder, cfg.eoSearchHalfWidth, ...
    cfg.freqSearchHz(1), cfg.freqSearchHz(2));
fprintf('Prev-window EO soft candidate: %s, refresh every %d window(s).\n', ...
    cfg.prevWindowCandidateMode, cfg.prevWindowRefreshEvery);
fprintf('Phase-safe expansion: %d, margin %.3f mm, fallback %s.\n', ...
    cfg.phaseSafeExpansion, cfg.phaseSafeMarginMm, cfg.phaseSafeFallbackMode);
fprintf('Direct bundle source mode: %s.\n', cfg.directBundleMode);
fprintf(['Pulse/domain: %s + %s, soft margin %.3f mm, query guard %s ' ...
    '(fixed %.3f, min %.3f, max %.3f, q%.1f + %.3f mm)\n'], ...
    cfg.pulseSelectionMode, cfg.domainSelectionMode, cfg.domainSoftMarginMm, ...
    cfg.queryGuardMode, cfg.queryGuardFixedMm, cfg.queryGuardMinMm, ...
    cfg.queryGuardMaxMm, cfg.queryGuardQuantile, cfg.queryGuardSafetyMm);
fprintf(['Models: %s | main %s, EO screen %s | dx limit %.3f mm, ' ...
    'dg %.3f mm, dmu %.3f, dtau %.3f mm\n'], ...
    strjoin(cfg.modelNames, ', '), cfg.mainModel, cfg.mainEoScreenModel, ...
    cfg.dxLimitMm, cfg.deltaGapLimitMm, cfg.deltaMuLimit, cfg.deltaTauLimitMm);
fprintf('Delta-g projection limit: use %d, alpha %.3g, min %.3f mm, max %.3f mm, fallback %.3f mm\n', ...
    cfg.useDeltaGapProjectionLimit, cfg.deltaGapProjAlpha, cfg.deltaGapProjMinLimitMm, ...
    cfg.deltaGapProjMaxLimitMm, cfg.deltaGapProjFallbackLimitMm);
fprintf('Eta policy: zero-start residual sensor shift, fallback limit %.3f mm, reg %.3f mV\n', ...
    cfg.etaFallbackLimitMm, cfg.etaRegWeightMv);

WindowResult = struct([]);
trendRows = cell(numWindows, 1);
bestIdx = 1;
bestRmse = inf;
prevWindowResult = [];

for iw = 1:numWindows
    Wmap = DynamicMap.Window(iw);
    prevPhaseInfo = build_phase_safe_reference_local(prevWindowResult, [], iw, cfg, 'prev_window');
    useDirectBundle = strcmpi(cfg.directBundleMode, 'auto') || ...
        (strcmpi(cfg.directBundleMode, 'off_when_phase_safe') && ~cfg.phaseSafeExpansion);
    if useDirectBundle && numel(directRef.WindowResult) >= iw && isfield(directRef.WindowResult(iw), 'bundle') && ...
            ~isempty(directRef.WindowResult(iw).bundle)
        coreBundleFull = build_gap_bundle_from_direct_bundle_local( ...
            directRef.WindowResult(iw).bundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
    else
        coreBundleFull = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
    end
    coreBundle = decimate_bundle_local(coreBundleFull, cfg.maxPointsPerWindow);
    eoCandidates = build_eo_candidates_local(coreBundle.rotFreqMeanHz, cfg);
    coreSeedTable = solve_vp_seed_scan_local(coreBundle, eoCandidates, cfg);
    phaseInfo = prevPhaseInfo;
    if cfg.phaseSafeExpansion && ~phaseInfo.usable && strcmpi(cfg.phaseSafeFallbackMode, 'linear_vp')
        phaseInfo = build_phase_safe_reference_local([], coreSeedTable, iw, cfg, prevPhaseInfo.reason);
    end

    bundleFull = coreBundleFull;
    usedPhaseSafeExpansion = false;
    expandedBundleFull = [];
    if cfg.phaseSafeExpansion && phaseInfo.usable
        expandedBundleFull = build_gap_observation_bundle_local( ...
            Wmap, Template, CorrectedGapLibrary, responseSurface, cfg, phaseInfo);
        if expandedBundleFull.pointCount > coreBundleFull.pointCount
            bundleFull = expandedBundleFull;
            usedPhaseSafeExpansion = true;
        end
    end
    bundle = decimate_bundle_local(bundleFull, cfg.maxPointsPerWindow);
    seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg);
    prevEoInfo = build_prev_window_eo_reference_local(prevWindowResult, iw, cfg);
    [selectedEO, selectionInfo] = select_eo_candidates_local( ...
        seedTable, cfg.vpTopK, directTrend, iw, cfg.eoCandidateMode, prevEoInfo);
    [coreSelectedEO, ~] = select_eo_candidates_local( ...
        coreSeedTable, cfg.vpTopK, directTrend, iw, cfg.eoCandidateMode, prevEoInfo);
    selectedEO = unique([selectedEO, coreSelectedEO], 'stable');
    if phaseInfo.usable
        selectedEO = unique([selectedEO, round(phaseInfo.EO)], 'stable');
    end
    selectionInfo.selectedEO = selectedEO;
    selectionInfo.coreSelectedEO = coreSelectedEO;
    selectionInfo.phaseSafeExpansionEnabled = logical(cfg.phaseSafeExpansion);
    selectionInfo.usedPhaseSafeExpansion = usedPhaseSafeExpansion;
    selectionInfo.phaseSafeReferenceEO = phaseInfo.EO;
    selectionInfo.phaseSafeReferenceSource = phaseInfo.source;
    selectionInfo.phaseSafeReferenceReason = phaseInfo.reason;
    selectionInfo.phaseSafeReferenceUsedPreviousWindow = phaseInfo.usedPreviousWindow;
    selectionInfo.phaseSafeCorePointCount = coreBundleFull.pointCount;
    selectionInfo.phaseSafeSelectedPointCount = bundleFull.pointCount;
    selectionInfo.phaseSafeCandidatePointCount = coreBundleFull.pointCount;
    if ~isempty(expandedBundleFull)
        selectionInfo.phaseSafeCandidatePointCount = expandedBundleFull.pointCount;
    end
    mainCandidateEO = selectedEO;

    cfgWindow = cfg;
    cfgWindow.sensorEtaLimitMm = collect_sensor_eta_limit_from_bundle_local(bundle);
    if cfgWindow.useDeltaGapProjectionLimit
        deltaGapProjection = estimate_delta_gap_projection_limits_local(bundle, seedTable, selectedEO, cfgWindow);
        cfgWindow.deltaGapLimitMm = deltaGapProjection.deltaGapLimitMm(:).';
        cfgWindow.DeltaGapProjection = deltaGapProjection;
        selectionInfo.DeltaGapProjection = deltaGapProjection;
    else
        selectionInfo.DeltaGapProjection = empty_delta_gap_projection_local(bundle, cfgWindow);
    end
    modelFits = struct();
    for im = 1:numel(cfg.modelNames)
        modeName = cfg.modelNames{im};
        modelFits.(modeName) = refine_static_warp_fit_local(bundle, seedTable, selectedEO, cfgWindow, modeName);
    end
    wr = pack_window_result_local(iw, Wmap, bundleFull, seedTable, ...
        mainCandidateEO, selectionInfo, modelFits);
    if iw == 1
        WindowResult = repmat(wr, numWindows, 1);
    else
        WindowResult(iw) = wr;
    end
    trendRows{iw} = make_trend_row_local(wr, directTrend, cfg);
    prevWindowResult = wr;
    if modelFits.(cfg.mainModel).weightedRmseMv < bestRmse
        bestRmse = modelFits.(cfg.mainModel).weightedRmseMv;
        bestIdx = iw;
    end
    fprintf('Window %02d/%02d laps %s: %s %.2f mV\n', ...
        iw, numWindows, mat2str(Wmap.lap_range), cfg.mainModel, ...
        modelFits.(cfg.mainModel).weightedRmseMv);
end

Trend = vertcat(trendRows{:});
Summary = build_summary_table_local(Trend, bestIdx, cfg);

Result = struct();
Result.dataset = '20251222';
if strcmpi(cfg.runMode, 'comparison')
    Result.method = 'gap_tilt_comparison_with_nested_static_warp_ablation_vp_full_wave';
else
    Result.method = 'gap_tilt_main_projection_limited_vp_full_wave';
end
Result.description = ['Current Step07J route uses the PrevWinPhaseSafe low-speed direct-template result as direct/fixed baseline and, by default, reuses its validated bundle. ' ...
    'The unified paper-facing gap-aware model is gap_tilt: VP Top-K plus direct/previous-window references provide only the EO candidate set, while gap_tilt jointly identifies EO, vibration, gap increment and a small x-dependent tilt increment with the full waveform objective. ' ...
    'The gap-increment search half-width is computed window-by-window from Delta g_proj = epsilon_perp / S_g and Delta g_loc = alpha Delta g_proj, with only explicit numerical floor/cap/fallback guards. ' ...
    'The default main run computes only gap_tilt. fixed and gap_only are computed only in STEP07J_RUN_MODE=comparison. gap_tilt_shift is excluded from the current comparison because the extra sensor-shift freedom is invalid for the paper-facing method.'];
Result.cfg = cfg;
Result.ResonanceSelection = resonanceSelection;
Result.templateFile = templateFile;
Result.dynamicFile = dynamicFile;
Result.dynamicSourceMode = dynamicSourceMode;
Result.highMapMeta = highMapMeta;
Result.correctedLibFile = correctedLibFile;
Result.correctedLibMode = correctedLibMode;
Result.gapCalibrationIndexFile = find_gap_calibration_index_file_local(outDir, C0);
Result.gapCalibrationSourceFile = correctedLibFile;
Result.directMainFile = directMainFile;
Result.CoordinateCheck = coordinateCheck;
Result.WindowResult = WindowResult;
Result.Trend = Trend;
Result.Summary = Summary;
Result.BestWindowIndex = bestIdx;
Result.MainModel = cfg.mainModel;
Result.BestWindow = WindowResult(bestIdx);

matFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s%s.mat', ...
    C0.dataset, C0.caseTag, resultSuffix));
trendCsv = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s%s.csv', ...
    C0.dataset, C0.caseTag, resultSuffix));
summaryCsv = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Summary_%s_%s%s.csv', ...
    C0.dataset, C0.caseTag, resultSuffix));
figTrend = fullfile(figDir, sprintf('Step07J_Trend_%s_%s%s.png', C0.dataset, C0.caseTag, resultSuffix));
figBest = fullfile(figDir, sprintf('Step07J_BestWindow_Collapse_%s_%s%s.png', C0.dataset, C0.caseTag, resultSuffix));
figGap = fullfile(figDir, sprintf('Step07J_StaticWarp_%s_%s%s.png', C0.dataset, C0.caseTag, resultSuffix));

save(matFile, 'Result', 'Trend', 'Summary', '-v7');
writetable(Trend, trendCsv);
writetable(Summary, summaryCsv);
plot_trend_local(Trend, figTrend, cfg);
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
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_OPRCenterStd_20251222.mat', ['B', num2str(targetBlade), '_', sensorTag])
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

function highMapFile = ensure_highmap_file_local(outDir, thisDir, C0, cfg)
expectedFile = expected_highmap_file_local(outDir, C0, cfg);
if isfile(expectedFile)
    highMapFile = expectedFile;
    return;
end
fprintf('Step07J highMap cache miss. Running Step06 to build:\n  %s\n', expectedFile);
STEP06_SKIP_CLEAR = true; %#ok<NASGU>
run(fullfile(thisDir, 'Step06_Build_HighMap_20251222.m'));
if ~isfile(expectedFile)
    error('Step06 finished but expected highMap cache was not created: %s', expectedFile);
end
highMapFile = expectedFile;
end

function dynamicMapFile = ensure_dynamic_map_file_local(outDir, thisDir, C0, cfg)
expectedFile = expected_dynamic_map_file_local(outDir, C0, cfg);
if isfile(expectedFile)
    dynamicMapFile = expectedFile;
    return;
end
highMapFile = ensure_highmap_file_local(outDir, thisDir, C0, cfg);
if ~isfile(expectedFile)
    fprintf('Step07J DynamicMap cache miss. Re-running Step06 to build:\n  %s\n', expectedFile);
    STEP06_SKIP_CLEAR = true; %#ok<NASGU>
    run(fullfile(thisDir, 'Step06_Build_HighMap_20251222.m'));
end
if ~isfile(expectedFile)
    error('Step06 finished but expected DynamicMap cache was not created: %s (highMap: %s)', ...
        expectedFile, highMapFile);
end
dynamicMapFile = expectedFile;
end

function highMapFile = expected_highmap_file_local(outDir, C0, cfg)
Pflow = C0.flowConfig;
cacheKey = make_highmap_cache_key_local(cfg.resonanceRegionShortTag, ...
    Pflow.identification.analysisStartTimeSec, ...
    Pflow.identification.targetBladePasses, ...
    Pflow.identification.windowBladePasses, ...
    Pflow.identification.slidingStepBladePasses, ...
    Pflow.identification.pulseWindowSec);
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_%s_%s_%s.mat', ...
    C0.dataset, C0.caseTag, cacheKey));
end

function dynamicMapFile = expected_dynamic_map_file_local(outDir, C0, cfg)
Pflow = C0.flowConfig;
cacheKey = make_highmap_cache_key_local(cfg.resonanceRegionShortTag, ...
    Pflow.identification.analysisStartTimeSec, ...
    Pflow.identification.targetBladePasses, ...
    Pflow.identification.windowBladePasses, ...
    Pflow.identification.slidingStepBladePasses, ...
    Pflow.identification.pulseWindowSec);
dynamicMapFile = fullfile(outDir, sprintf('Step06_DynamicMap_%s_%s_%s.mat', ...
    C0.dataset, C0.caseTag, cacheKey));
end

function key = make_highmap_cache_key_local(regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec)
key = sprintf('%s_T%09.3f_L%02d_W%02d_S%02d_PW%07.5f', ...
    regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec);
key = strrep(key, '.', 'p');
key = strrep(key, '+', '');
key = strrep(key, '-', 'm');
end

function [DynamicMap, meta] = load_dynamic_input_local(dynamicFile, sourceMode, cfg)
S = load(dynamicFile);
if isfield(S, 'DynamicMap')
    DynamicMap = S.DynamicMap;
    meta = struct('sourceType', 'DynamicMap', 'sourceFile', dynamicFile);
    return;
end
if ~isfield(S, 'highMap')
    error('Dynamic input file must contain DynamicMap or highMap: %s', dynamicFile);
end
highMap = S.highMap;
if isfield(cfg, 'lowSpeedTemplate') && ~isempty(cfg.lowSpeedTemplate)
    highMap.runtimeTemplate = cfg.lowSpeedTemplate;
end
DynamicMap = convert_highmap_to_dynamic_map_local(highMap, cfg, dynamicFile, sourceMode);
meta = struct();
meta.sourceType = 'highMap';
meta.sourceFile = dynamicFile;
meta.cacheKey = read_optional_field_local(highMap, 'cacheKey', '');
meta.analysisStartTime = read_optional_field_local(highMap, 'analysisStartTime', NaN);
meta.targetLaps = read_optional_field_local(highMap, 'targetLaps', NaN);
meta.analysisWinSize = read_optional_field_local(highMap, 'analysisWinSize', NaN);
meta.slidingStep = read_optional_field_local(highMap, 'slidingStep', NaN);
end

function DynamicMap = convert_highmap_to_dynamic_map_local(highMap, cfg, sourceFile, sourceMode)
if isfield(highMap, 'experiment_case') && ~isempty(highMap.experiment_case)
    DynamicMap = convert_experiment_case_to_dynamic_map_local(highMap, cfg, sourceFile, sourceMode);
    return;
end
requiredFields = {'t_v','x_v','x_fit_v','V_a','S_v','rev_v','W_v','theta_v'};
for i = 1:numel(requiredFields)
    if ~isfield(highMap, requiredFields{i})
        error('highMap is missing field %s.', requiredFields{i});
    end
end
windowLaps = read_optional_field_local(highMap, 'analysisWinSize', 3);
slidingStep = read_optional_field_local(highMap, 'slidingStep', 1);
targetLaps = read_optional_field_local(highMap, 'targetLaps', max(highMap.rev_v));
starts = 1:slidingStep:(targetLaps - windowLaps + 1);
Window = repmat(struct('window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, 'Sensor', []), 1, numel(starts));
for iw = 1:numel(starts)
    lapRange = [starts(iw), starts(iw) + windowLaps - 1];
    Sensor = repmat(struct('sensor_id', NaN, 't', [], 'x_abs', [], 'x_rel', [], ...
        'theta', [], 'V', [], 'V_unit', 'mV_baseline_corrected', 'W', []), 1, numel(cfg.analysisSensors));
    for is = 1:numel(cfg.analysisSensors)
        sid = cfg.analysisSensors(is);
        mask = highMap.S_v(:) == sid & highMap.rev_v(:) >= lapRange(1) & ...
            highMap.rev_v(:) <= lapRange(2);
        Sensor(is).sensor_id = sid;
        Sensor(is).t = highMap.t_v(mask);
        Sensor(is).x_rel = select_highmap_x_local(highMap, mask);
        Sensor(is).x_abs = select_highmap_x_abs_local(highMap, mask, Sensor(is).x_rel);
        Sensor(is).theta = highMap.theta_v(mask);
        Sensor(is).V = highMap.V_a(mask);
        Sensor(is).W = highMap.W_v(mask);
    end
    Window(iw).window_id = iw;
    Window(iw).lap_range = lapRange;
    Window(iw).time_window = infer_window_time_range_local(Sensor);
    Window(iw).rot_freq_mean_hz = read_optional_field_local(highMap, 'rotFreqHz', NaN);
    Window(iw).rot_rpm_mean = read_optional_field_local(highMap, 'rotRpm', 60 * Window(iw).rot_freq_mean_hz);
    Window(iw).Sensor = Sensor;
end
DynamicMap = struct();
DynamicMap.Dataset = read_optional_field_local(highMap, 'dataset', '20251222');
DynamicMap.SourceMode = sourceMode;
DynamicMap.SourceHighMapFile = sourceFile;
DynamicMap.TargetBlade = read_optional_field_local(highMap, 'targetBlade', cfg.targetBlade);
DynamicMap.SensorIDs = cfg.analysisSensors(:).';
DynamicMap.SensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
DynamicMap.RegionShortTag = read_optional_field_local(highMap, 'resonanceRegionShortTag', cfg.resonanceRegionShortTag);
DynamicMap.TemplateFileForXRel = read_optional_field_local(highMap, 'lowSpeedTemplateFile', '');
DynamicMap.Window = Window;
end

function DynamicMap = convert_experiment_case_to_dynamic_map_local(highMap, cfg, sourceFile, sourceMode)
E = highMap.experiment_case;
if ~isfield(E, 'Extracted_Data')
    error('highMap.experiment_case does not contain Extracted_Data.');
end
windowLaps = read_optional_field_local(highMap, 'analysisWinSize', 3);
slidingStep = read_optional_field_local(highMap, 'slidingStep', 1);
targetLaps = read_optional_field_local(highMap, 'targetLaps', inf);
if ~isfinite(targetLaps)
    targetLaps = infer_min_lap_count_local(E, cfg.analysisSensors);
end
starts = 1:slidingStep:(targetLaps - windowLaps + 1);
Window = repmat(struct('window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, 'Sensor', []), 1, numel(starts));
for iw = 1:numel(starts)
    lapRange = starts(iw):(starts(iw) + windowLaps - 1);
    Sensor = repmat(struct('sensor_id', NaN, 't', [], 'x_abs', [], ...
        'x_rel', [], 'V', [], 'W', [], 'theta', [], 'point_count', NaN), ...
        1, numel(cfg.analysisSensors));
    allT = [];
    for is = 1:numel(cfg.analysisSensors)
        sid = cfg.analysisSensors(is);
        dataIdx = find([E.Extracted_Data.sensor_id] == sid, 1, 'first');
        if isempty(dataIdx)
            error('experiment_case.Extracted_Data does not contain CH%d.', sid);
        end
        laps = E.Extracted_Data(dataIdx).Laps;
        Tpl = get_embedded_template_sensor_local(highMap, sid);
        t = [];
        xAbs = [];
        v = [];
        theta = [];
        for lapId = lapRange
            if lapId > numel(laps)
                continue;
            end
            D = laps(lapId);
            t = [t; D.t_points(:)]; %#ok<AGROW>
            xAbs = [xAbs; D.x_points(:)]; %#ok<AGROW>
            v = [v; D.v_points(:)]; %#ok<AGROW>
            theta = [theta; infer_lap_theta_local(E, D, lapId)]; %#ok<AGROW>
        end
        xRel = xAbs - Tpl.xc;
        W = build_simple_waveform_weight_local(v);
        Sensor(is).sensor_id = sid;
        Sensor(is).t = t(:);
        Sensor(is).x_abs = xAbs(:);
        Sensor(is).x_rel = xRel(:);
        Sensor(is).V = v(:);
        Sensor(is).W = W(:);
        Sensor(is).theta = theta(:);
        Sensor(is).point_count = numel(t);
        allT = [allT; t(:)]; %#ok<AGROW>
    end
    Window(iw).window_id = iw;
    Window(iw).lap_range = lapRange;
    finiteT = allT(isfinite(allT));
    if isempty(finiteT)
        Window(iw).time_window = [NaN, NaN];
    else
        Window(iw).time_window = [min(finiteT), max(finiteT)];
    end
    Window(iw).rot_freq_mean_hz = read_optional_field_local(highMap, 'rotFreqHz', NaN);
    Window(iw).rot_rpm_mean = read_optional_field_local(highMap, 'rotRpm', 60 * Window(iw).rot_freq_mean_hz);
    Window(iw).Sensor = Sensor;
end
DynamicMap = struct();
DynamicMap.Route = 'gap_prior_local_experiment_case_full_waveform_dynamic_map';
DynamicMap.SourceMode = sourceMode;
DynamicMap.SourceHighMapFile = sourceFile;
DynamicMap.TargetBlade = read_optional_field_local(highMap, 'targetBlade', cfg.targetBlade);
DynamicMap.SensorIDs = cfg.analysisSensors(:).';
DynamicMap.SensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
DynamicMap.RegionShortTag = read_optional_field_local(highMap, 'resonanceRegionShortTag', cfg.resonanceRegionShortTag);
DynamicMap.TemplateFileForXRel = read_optional_field_local(highMap, 'lowSpeedTemplateFile', '');
DynamicMap.ResonanceSelection = read_optional_field_local(highMap, 'region', table());
DynamicMap.GlobalTimeWindow = infer_dynamic_map_time_window_local(Window);
DynamicMap.Window = Window;
end

function nLap = infer_min_lap_count_local(E, sensorIds)
nLap = inf;
for sid = sensorIds(:).'
    dataIdx = find([E.Extracted_Data.sensor_id] == sid, 1, 'first');
    if isempty(dataIdx)
        continue;
    end
    nLap = min(nLap, numel(E.Extracted_Data(dataIdx).Laps));
end
if ~isfinite(nLap)
    error('Could not infer lap count from experiment_case.');
end
end

function Tpl = get_embedded_template_sensor_local(highMap, sid)
if isfield(highMap, 'runtimeTemplate') && isfield(highMap.runtimeTemplate, 'Sensor')
    idx = find([highMap.runtimeTemplate.Sensor.sensor_id] == sid, 1, 'first');
    if ~isempty(idx)
        Tpl = highMap.runtimeTemplate.Sensor(idx);
        return;
    end
end
if ~isfield(highMap, 'experiment_case') || ~isfield(highMap.experiment_case, 'LowSpeedTemplate')
    if isfield(highMap, 'lowSpeedTemplateFile') && exist(highMap.lowSpeedTemplateFile, 'file') == 2
        S = load(highMap.lowSpeedTemplateFile, 'Template');
        idx = find([S.Template.Sensor.sensor_id] == sid, 1, 'first');
        if ~isempty(idx)
            Tpl = S.Template.Sensor(idx);
            return;
        end
    end
    error('Cannot resolve low-speed template sensor CH%d from highMap.', sid);
end
Template = highMap.experiment_case.LowSpeedTemplate;
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Embedded low-speed template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function theta = infer_lap_theta_local(E, D, lapId)
t = D.t_points(:);
if isfield(D, 'theta_points') && numel(D.theta_points) == numel(t)
    theta = D.theta_points(:);
elseif isfield(E, 'Diagnostics') && isfield(E.Diagnostics, 'rot_freq_mean_hz') && isfinite(E.Diagnostics.rot_freq_mean_hz)
    theta = 2 * pi * E.Diagnostics.rot_freq_mean_hz .* (t - t(1));
else
    theta = 2 * pi * (lapId - 1) + linspace(0, 2*pi, numel(t)).';
end
end

function W = build_simple_waveform_weight_local(V)
V = V(:);
if isempty(V)
    W = V;
    return;
end
W = ones(size(V));
finite = isfinite(V);
if nnz(finite) >= 5
    g = abs(gradient(V));
    if max(g(finite), [], 'omitnan') > 0
        W = g ./ max(g(finite), [], 'omitnan');
    end
end
W = max(0.05, W);
end

function timeWindow = infer_dynamic_map_time_window_local(Window)
tMin = inf;
tMax = -inf;
for iw = 1:numel(Window)
    tw = Window(iw).time_window;
    if numel(tw) == 2 && all(isfinite(tw))
        tMin = min(tMin, tw(1));
        tMax = max(tMax, tw(2));
    end
end
if isfinite(tMin) && isfinite(tMax)
    timeWindow = [tMin, tMax];
else
    timeWindow = [NaN, NaN];
end
end

function x = select_highmap_x_local(highMap, mask)
if isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(mask) && any(isfinite(highMap.x_fit_v(mask)))
    x = highMap.x_fit_v(mask);
else
    x = highMap.x_v(mask);
end
end

function xAbs = select_highmap_x_abs_local(highMap, mask, fallbackX)
if isfield(highMap, 'x_mm_v') && numel(highMap.x_mm_v) == numel(mask) && any(isfinite(highMap.x_mm_v(mask)))
    xAbs = highMap.x_mm_v(mask);
else
    xAbs = fallbackX;
end
end

function timeWindow = infer_window_time_range_local(Sensor)
tAll = [];
for i = 1:numel(Sensor)
    tAll = [tAll; Sensor(i).t(:)]; %#ok<AGROW>
end
tAll = tAll(isfinite(tAll));
if isempty(tAll)
    timeWindow = [NaN, NaN];
else
    timeWindow = [min(tAll), max(tAll)];
end
end

function value = read_optional_field_local(S, fieldName, defaultValue)
if isstruct(S) && isfield(S, fieldName) && ~isempty(S.(fieldName))
    value = S.(fieldName);
else
    value = defaultValue;
end
end

function dynamicFile = find_dynamic_map_file_local(dynamicDir, targetBlade, sensorIds, regionShortTag)
sensorTag = ['S', sprintf('%d', sensorIds)];
caseTag = ['B', num2str(targetBlade), '_', sensorTag];
if nargin < 4
    regionShortTag = '';
end
patterns = {
    sprintf('DynamicMap_%s_SlidingWindows_%s_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', caseTag, regionShortTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_%s_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, sensorTag, regionShortTag)
    sprintf('DynamicMap_%s_SlidingWindows_%s*_20251222.mat', caseTag, regionShortTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_%s*_20251222.mat', targetBlade, sensorTag, regionShortTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', caseTag)
    sprintf('DynamicMap_%s_SlidingWindows*_20251222.mat', caseTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows*OPRCenterStd*_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20251222.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    files = dir(fullfile(dynamicDir, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        dynamicFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
error('No DynamicMap file found under %s for %s.', dynamicDir, sensorTag);
end

function correctedLibFile = find_corrected_gap_library_file_local(outDir, targetBlade, sensorTag)
patterns = {
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_%s.mat', ['B', num2str(targetBlade), '_', sensorTag])
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

function [CorrectedGapLibrary, sourceFile, sourceMode] = load_gap_calibration_library_local(outDir, C0, cfg)
CorrectedGapLibrary = [];
sourceFile = '';
sourceMode = 'missing';
if isfield(C0, 'flowConfig') && isfield(C0.flowConfig, 'calibration')
    libDirName = C0.flowConfig.calibration.libraryDirName;
    filePrefix = C0.flowConfig.calibration.filePrefix;
else
    libDirName = 'gap_calibration_library';
    filePrefix = 'GapCalib';
end
libDir = fullfile(outDir, libDirName);
if exist(libDir, 'dir') ~= 7
    return;
end
sensorIds = cfg.analysisSensors(:).';
sensorEntries = [];
baseLib = [];
sourceCombinedFile = '';
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    patterns = {
        sprintf('%s_%s_B%d_S%d.mat', filePrefix, C0.dataset, C0.bladeId, sid)
        sprintf('%s_%s_B%d_S%d*.mat', filePrefix, C0.dataset, C0.bladeId, sid)
        };
    sensorFile = '';
    for ip = 1:numel(patterns)
        files = dir(fullfile(libDir, patterns{ip}));
        if ~isempty(files)
            [~, idx] = max([files.datenum]);
            sensorFile = fullfile(files(idx).folder, files(idx).name);
            break;
        end
    end
    if isempty(sensorFile) || exist(sensorFile, 'file') ~= 2
        CorrectedGapLibrary = [];
        return;
    end
    S = load(sensorFile, 'GapCalib', 'SingleCorrectedGapLibrary');
    if isfield(S, 'GapCalib') && isfield(S.GapCalib, 'correctedGapLibrary')
        thisLib = S.GapCalib.correctedGapLibrary;
    elseif isfield(S, 'SingleCorrectedGapLibrary')
        thisLib = S.SingleCorrectedGapLibrary;
    else
        CorrectedGapLibrary = [];
        return;
    end
    if isempty(baseLib)
        baseLib = thisLib;
        if isfield(S, 'GapCalib') && isfield(S.GapCalib, 'sourceCombinedFile')
            sourceCombinedFile = S.GapCalib.sourceCombinedFile;
        elseif isfield(thisLib, 'sourceCombinedFile')
            sourceCombinedFile = thisLib.sourceCombinedFile;
        end
    end
    if ~isfield(thisLib, 'sensor') || isempty(thisLib.sensor)
        CorrectedGapLibrary = [];
        return;
    end
    sensorEntries = [sensorEntries; thisLib.sensor(1)]; %#ok<AGROW>
    if isempty(sourceFile)
        sourceFile = sensorFile;
    end
end
if ~isfield(baseLib, 'responseSurface') || isempty(baseLib.responseSurface)
    if isempty(sourceCombinedFile) || exist(sourceCombinedFile, 'file') ~= 2
        CorrectedGapLibrary = [];
        return;
    end
    C = load(sourceCombinedFile, 'CorrectedGapLibrary');
    if ~isfield(C, 'CorrectedGapLibrary') || ~isfield(C.CorrectedGapLibrary, 'responseSurface')
        CorrectedGapLibrary = [];
        return;
    end
    baseLib.responseSurface = C.CorrectedGapLibrary.responseSurface;
    baseLib.responseFile = C.CorrectedGapLibrary.responseFile;
    baseLib.templateFile = C.CorrectedGapLibrary.templateFile;
end
CorrectedGapLibrary = baseLib;
CorrectedGapLibrary.dataset = C0.dataset;
CorrectedGapLibrary.targetBlade = C0.bladeId;
CorrectedGapLibrary.analysisSensors = sensorIds;
CorrectedGapLibrary.sensor = sensorEntries;
CorrectedGapLibrary.libraryIndexFile = find_gap_calibration_index_file_local(outDir, C0);
sourceMode = 'calibration_library';
end

function indexFile = find_gap_calibration_index_file_local(outDir, C0)
if isfield(C0, 'flowConfig') && isfield(C0.flowConfig, 'calibration')
    libDirName = C0.flowConfig.calibration.libraryDirName;
    filePrefix = C0.flowConfig.calibration.filePrefix;
else
    libDirName = 'gap_calibration_library';
    filePrefix = 'GapCalib';
end
libDir = fullfile(outDir, libDirName);
patterns = {
    sprintf('%s_Index_%s_%s.csv', filePrefix, C0.dataset, C0.caseTag)
    sprintf('%s_Index_%s_%s*.csv', filePrefix, C0.dataset, C0.caseTag)
    };
indexFile = '';
for i = 1:numel(patterns)
    files = dir(fullfile(libDir, patterns{i}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        indexFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
end

function resultFile = find_direct_main_result_file_local(resultDir, targetBlade, sensorIds, regionShortTag)
sensorTag = ['S', sprintf('%d', sensorIds)];
caseTag = ['B', num2str(targetBlade), '_', sensorTag];
if nargin < 4
    regionShortTag = '';
end
preferred = {
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s_%s*DirectTemplate_OPRCenterStd*PrevWinPhaseSafe*_20251222.mat', caseTag, regionShortTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s*DirectTemplate_OPRCenterStd*PrevWinPhaseSafe*_20251222.mat', targetBlade, sensorTag, regionShortTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s*DirectTemplate_OPRCenterStd*PrevWinPhaseSafe*_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*DirectTemplate_OPRCenterStd*PrevWinPhaseSafe*_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s_Main_DirectTemplate_OPRCenterStd_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_DirectTemplate_OPRCenterStd_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s*DirectTemplate_OPRCenterStd*_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*DirectTemplate_OPRCenterStd*_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s*OPRCenterStd*_20251222.mat', caseTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*OPRCenterStd*_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_OPRAnchored_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_OPRAnchored*_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*GradientXRange030*_20251222.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(preferred)
    files = dir(fullfile(resultDir, preferred{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        resultFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
files = dir(fullfile(resultDir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*_20251222.mat', ...
    targetBlade, sensorTag)));
if isempty(files)
    files = dir(fullfile(resultDir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s*_20251222.mat', ...
        ['B', num2str(targetBlade), '_', sensorTag])));
end
if isempty(files)
    resultFile = '';
else
    [~, idx] = max([files.datenum]);
    resultFile = fullfile(files(idx).folder, files(idx).name);
end
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

function directRef = load_direct_reference_local(resultFile)
directRef = struct('Trend', table(), 'WindowResult', struct([]));
if isempty(resultFile) || ~isfile(resultFile)
    return;
end
S = load(resultFile, 'Result');
if isfield(S, 'Result') && isfield(S.Result, 'Trend')
    directRef.Trend = S.Result.Trend;
end
if isfield(S, 'Result') && isfield(S.Result, 'WindowResult')
    directRef.WindowResult = S.Result.WindowResult;
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
    sensorInfo(is).etaLimitMm = resolve_sensor_eta_limit_local(Tpl, corr, cfg);
    sensorInfo(is).etaMedianMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaMedianMm');
    sensorInfo(is).etaIqrMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaIqrMm');
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

function bundle = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg, phaseRef)
if nargin < 6
    phaseRef = [];
end
usePhaseSafe = isstruct(phaseRef) && isfield(phaseRef, 'usable') && ...
    phaseRef.usable && cfg.phaseSafeExpansion;
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
    if isfield(D, 'V_unit') && strcmpi(char(string(D.V_unit)), 'mV_baseline_corrected')
        vRawMv = D.V(:);
        vForInversion = Tpl.baseline + D.V(:) / 1000;
        thresholdMv = (resolve_threshold_local(Tpl, cfg) - Tpl.baseline) * 1000;
    else
        vRawMv = (D.V(:) - Tpl.baseline) * 1000;
        vForInversion = D.V(:);
        thresholdMv = (resolve_threshold_local(Tpl, cfg) - Tpl.baseline) * 1000;
    end
    wRaw = max(D.W(:), cfg.weightFloor);

    f0Raw = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw);
    fxRaw = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw, cfg.derivativeStepMm);
    maskPulse = build_pulse_mask_local(vRawMv, thresholdMv, cfg.pulseSelectionMode);
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRawMv, Tpl, thresholdMv, cfg);
    maskDomain = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRawMv) & isfinite(thetaRaw);
    mask = maskPulse & maskEffective & maskDomain;
    queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vForInversion, mask, cfg);
    xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    maskQuerySafe = xRaw >= xDomain(1) + queryGuardMm & xRaw <= xDomain(2) - queryGuardMm;
    maskPhaseSafe = false(size(maskQuerySafe));
    if usePhaseSafe
        maskPhaseSafe = build_phase_safe_query_mask_local(xRaw, thetaRaw, xDomain, is, phaseRef, cfg);
        if nnz(mask & maskPhaseSafe) >= 8
            maskQuerySafe = maskPhaseSafe;
        end
    end
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
    sensorInfo(is).phaseSafeExpansion = usePhaseSafe;
    sensorInfo(is).phaseSafeCount = nnz(mask & maskPhaseSafe);
    sensorInfo(is).etaLimitMm = resolve_sensor_eta_limit_local(Tpl, corr, cfg);
    sensorInfo(is).etaMedianMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaMedianMm');
    sensorInfo(is).etaIqrMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaIqrMm');
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
if usePhaseSafe
    bundle.selectionPass = 'phase_safe_expanded';
else
    bundle.selectionPass = 'core_hard_adaptive';
end
bundle.PhaseSafeReferenceInfo = phaseRef;
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

function info = build_phase_safe_reference_local(prevResult, seedTable, windowId, cfg, fallbackReason)
info = empty_phase_safe_reference_local();
if nargin < 5 || isempty(fallbackReason)
    fallbackReason = 'not_specified';
end
if ~isempty(prevResult) && isstruct(prevResult)
    if cfg.prevWindowRefreshEvery > 0 && mod(windowId - 1, cfg.prevWindowRefreshEvery) == 0
        info.reason = 'scheduled_refresh';
        return;
    end
    if ~isfield(prevResult, 'modelFits') || ~isfield(prevResult.modelFits, cfg.mainModel)
        info.reason = 'previous_missing_main_model';
        return;
    end
    prevFit = prevResult.modelFits.(cfg.mainModel);
    required = {'EO', 'amplitudeMm', 'phaseRad', 'dxMm', 'weightedRmseMv'};
    for i = 1:numel(required)
        if ~isfield(prevFit, required{i})
            info.reason = ['previous_missing_', required{i}];
            return;
        end
    end
    vals = [prevFit.EO, prevFit.amplitudeMm, prevFit.phaseRad, ...
        prevFit.dxMm, prevFit.weightedRmseMv];
    if any(~isfinite(vals))
        info.reason = 'previous_nonfinite_reference';
        return;
    end
    info.usable = true;
    info.EO = round(prevFit.EO);
    info.A = min(abs(prevFit.amplitudeMm), cfg.amplitudeLimitMm);
    info.phi = wrap_to_pi_local(prevFit.phaseRad);
    info.dx = max(min(prevFit.dxMm, cfg.dxLimitMm), -cfg.dxLimitMm);
    info.deltaTauMm = zeros(1, numel(cfg.analysisSensors));
    if isfield(prevFit, 'deltaTauMm') && ~isempty(prevFit.deltaTauMm)
        n = min(numel(info.deltaTauMm), numel(prevFit.deltaTauMm));
        info.deltaTauMm(1:n) = prevFit.deltaTauMm(1:n);
    end
    info.weightedRmseMv = prevFit.weightedRmseMv;
    info.source = 'prev_window_main_final';
    info.reason = 'prev_window_quality_pass';
    info.sourceWindowId = prevResult.windowId;
    info.usedPreviousWindow = true;
    return;
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
info.source = 'current_linear_vp_seed';
info.reason = [fallbackReason, '_then_linear_vp'];
info.sourceWindowId = windowId;
info.usedPreviousWindow = false;
end

function info = empty_phase_safe_reference_local()
info = struct('usable', false, 'EO', NaN, 'A', NaN, 'phi', NaN, ...
    'dx', NaN, 'deltaTauMm', [], 'weightedRmseMv', NaN, ...
    'source', 'none', 'reason', 'not_evaluated', ...
    'sourceWindowId', NaN, 'usedPreviousWindow', false);
end

function info = build_prev_window_eo_reference_local(prevResult, windowId, cfg)
info = struct('usable', false, 'EO', NaN, 'source', 'none', ...
    'reason', 'not_evaluated');
if strcmpi(cfg.prevWindowCandidateMode, 'off')
    info.reason = 'prev_window_candidate_off';
    return;
end
if isempty(prevResult) || ~isstruct(prevResult)
    info.reason = 'no_previous_window';
    return;
end
if cfg.prevWindowRefreshEvery > 0 && mod(windowId - 1, cfg.prevWindowRefreshEvery) == 0
    info.reason = 'scheduled_refresh';
    return;
end
if ~isfield(prevResult, 'modelFits') || ~isfield(prevResult.modelFits, cfg.mainModel)
    info.reason = 'previous_missing_main_model';
    return;
end
prevFit = prevResult.modelFits.(cfg.mainModel);
required = {'EO', 'weightedRmseMv'};
for i = 1:numel(required)
    if ~isfield(prevFit, required{i})
        info.reason = ['previous_missing_', required{i}];
        return;
    end
end
if ~isfinite(prevFit.EO) || ~isfinite(prevFit.weightedRmseMv)
    info.reason = 'previous_nonfinite_reference';
    return;
end
info.usable = true;
info.EO = round(prevFit.EO);
info.source = 'prev_window_main_final';
info.reason = 'prev_window_quality_pass';
end

function [eoKeep, info] = select_eo_candidates_local(seedTable, topK, directTrend, windowId, mode, prevEoInfo)
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
if nargin >= 6 && isstruct(prevEoInfo) && isfield(prevEoInfo, 'usable') && prevEoInfo.usable
    eoKeep = unique([eoKeep, round(prevEoInfo.EO)], 'stable');
end
info = struct();
info.mode = mode;
info.vpEO = vpKeep;
info.directEO = directEO;
info.selectedEO = eoKeep;
info.usedPrevWindow = nargin >= 6 && isstruct(prevEoInfo) && ...
    isfield(prevEoInfo, 'usable') && prevEoInfo.usable;
if nargin >= 6 && isstruct(prevEoInfo)
    info.prevEO = prevEoInfo.EO;
    info.prevSource = prevEoInfo.source;
    info.prevReason = prevEoInfo.reason;
else
    info.prevEO = NaN;
    info.prevSource = 'none';
    info.prevReason = 'not_provided';
end
end

function fit = refine_static_warp_fit_local(bundle, seedTable, eoCandidates, cfg, modeName)
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'residualObjectiveMv2', inf, 'overshootPenalty', inf, ...
    'staticRegPenalty', inf, 'etaRegPenalty', inf, 'invalidPenalty', inf, 'penaltyObjectiveMv2', inf, ...
    'penaltyShare', NaN, ...
    'sensorEtaMm', [], 'deltaGapMm', [], 'deltaMuGapPerXMm', [], 'deltaTauMm', [], ...
    'VPred', [], 'uMm', [], 'score', inf), numel(eoCandidates), 1);
best = struct('score', inf);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    theta0 = initial_theta_local(seed.A, seed.phi, seed.dx, modeName, numel(bundle.sensorIds));
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfg, modeName);
    thetaOpt = fminsearch(fun, theta0, cfg.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfg, modeName, numel(bundle.sensorIds));
    detail = evaluate_static_warp_waveform_local(thetaOpt, eo, bundle, cfg, modeName);
    candidate(i).EO = eo;
    candidate(i).freqHz = eo * bundle.rotFreqMeanHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    candidate(i).residualObjectiveMv2 = detail.residualObjectiveMv2;
    candidate(i).overshootPenalty = detail.overshootPenalty;
    candidate(i).staticRegPenalty = detail.staticRegPenalty;
    candidate(i).etaRegPenalty = detail.etaRegPenalty;
    candidate(i).invalidPenalty = detail.invalidPenalty;
    candidate(i).penaltyObjectiveMv2 = detail.penaltyObjectiveMv2;
    candidate(i).penaltyShare = detail.penaltyShare;
    candidate(i).sensorEtaMm = detail.sensorEtaMm;
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
fit.CandidateTable = sortrows(struct2table(rmfield(candidate, {'VPred','uMm','sensorEtaMm','deltaGapMm','deltaMuGapPerXMm','deltaTauMm'})), ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
end

function proj = estimate_delta_gap_projection_limits_local(bundle, seedTable, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
nSensor = numel(bundle.sensorIds);
rows = repmat(struct('sensorId', NaN, 'EO', NaN, 'epsilonPerpMv', NaN, ...
    'SgMvPerMm', NaN, 'deltaGapProjMm', NaN, 'deltaGapLsMm', NaN, ...
    'deltaGapLocMm', NaN, 'deltaGapLimitMm', NaN, 'usedFallback', true, ...
    'clippedMin', false, 'clippedMax', false), nSensor, 1);
for is = 1:nSensor
    bestRow = rows(is);
    bestLimit = -inf;
    for ie = 1:numel(eoCandidates)
        eo = eoCandidates(ie);
        seedIdx = find(seedTable.EO == eo, 1, 'first');
        if isempty(seedIdx)
            seedIdx = 1;
        end
        seed = seedTable(seedIdx, :);
        row = estimate_delta_gap_projection_one_sensor_local(bundle, seed.A, seed.phi, seed.dx, eo, is, cfg);
        if row.deltaGapLimitMm > bestLimit
            bestLimit = row.deltaGapLimitMm;
            bestRow = row;
        end
    end
    rows(is) = bestRow;
end
proj = struct();
proj.alpha = cfg.deltaGapProjAlpha;
proj.sensorId = [rows.sensorId].';
proj.EO = [rows.EO].';
proj.epsilonPerpMv = [rows.epsilonPerpMv].';
proj.SgMvPerMm = [rows.SgMvPerMm].';
proj.deltaGapProjMm = [rows.deltaGapProjMm].';
proj.deltaGapLsMm = [rows.deltaGapLsMm].';
proj.deltaGapLocMm = [rows.deltaGapLocMm].';
proj.deltaGapLimitMm = [rows.deltaGapLimitMm].';
proj.usedFallback = [rows.usedFallback].';
proj.clippedMin = [rows.clippedMin].';
proj.clippedMax = [rows.clippedMax].';
proj.Table = struct2table(rows);
end

function row = estimate_delta_gap_projection_one_sensor_local(bundle, A, phi, dx, eo, sensorLocalIndex, cfg)
sid = bundle.sensorIds(sensorLocalIndex);
mask = bundle.sensorIndex == sensorLocalIndex;
Tpl = get_template_sensor_local(bundle.Template, sid);
corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
theta = bundle.Theta(mask);
u = A .* sin(eo .* theta + phi);
xEval = bundle.X(mask) - dx - u;
[v0, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval);
hG = max(cfg.deltaGapDerivativeStepMm, eps);
[vGp, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, hG, 0, 0, xEval);
[vGm, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, -hG, 0, 0, xEval);
hX = max(cfg.derivativeStepMm, eps);
[vXp, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval + hX);
[vXm, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval - hX);
r = bundle.V(mask) - v0;
jg = (vGp - vGm) ./ (2 * hG);
dMdx = (vXp - vXm) ./ (2 * hX);
Jnu = [-dMdx, -dMdx .* sin(eo .* theta), -dMdx .* cos(eo .* theta)];
valid = isfinite(r) & isfinite(jg) & isfinite(bundle.W(mask)) & all(isfinite(Jnu), 2);
sw = sqrt(max(bundle.W(mask), cfg.weightFloor));
rw = sw(valid) .* r(valid);
jgw = sw(valid) .* jg(valid);
Jw = Jnu(valid, :) .* sw(valid);
rPerp = orthogonal_residual_local(rw, Jw);
jgPerp = orthogonal_residual_local(jgw, Jw);
epsilonPerp = norm(rPerp);
Sg = norm(jgPerp);
deltaLs = NaN;
if isfinite(Sg) && Sg > cfg.deltaGapProjSensitivityFloorMvPerMm
    deltaProj = epsilonPerp / Sg;
    deltaLs = abs((jgPerp(:).' * rPerp(:)) / max(Sg^2, eps));
    usedFallback = false;
else
    deltaProj = cfg.deltaGapProjFallbackLimitMm / max(cfg.deltaGapProjAlpha, eps);
    usedFallback = true;
end
deltaLoc = cfg.deltaGapProjAlpha * deltaProj;
limit = deltaLoc;
clippedMin = false;
clippedMax = false;
if ~isfinite(limit) || limit <= 0
    limit = cfg.deltaGapProjFallbackLimitMm;
    usedFallback = true;
end
if limit < cfg.deltaGapProjMinLimitMm
    limit = cfg.deltaGapProjMinLimitMm;
    clippedMin = true;
end
if limit > cfg.deltaGapProjMaxLimitMm
    limit = cfg.deltaGapProjMaxLimitMm;
    clippedMax = true;
end
row = struct('sensorId', sid, 'EO', eo, 'epsilonPerpMv', epsilonPerp, ...
    'SgMvPerMm', Sg, 'deltaGapProjMm', deltaProj, 'deltaGapLsMm', deltaLs, ...
    'deltaGapLocMm', deltaLoc, 'deltaGapLimitMm', limit, 'usedFallback', usedFallback, ...
    'clippedMin', clippedMin, 'clippedMax', clippedMax);
end

function proj = empty_delta_gap_projection_local(bundle, cfg)
nSensor = numel(bundle.sensorIds);
proj = struct();
proj.alpha = cfg.deltaGapProjAlpha;
proj.sensorId = bundle.sensorIds(:);
proj.EO = NaN(nSensor, 1);
proj.epsilonPerpMv = NaN(nSensor, 1);
proj.SgMvPerMm = NaN(nSensor, 1);
proj.deltaGapProjMm = NaN(nSensor, 1);
proj.deltaGapLsMm = NaN(nSensor, 1);
proj.deltaGapLocMm = NaN(nSensor, 1);
proj.deltaGapLimitMm = repmat(cfg.deltaGapLimitMm, nSensor, 1);
proj.usedFallback = true(nSensor, 1);
proj.clippedMin = false(nSensor, 1);
proj.clippedMax = false(nSensor, 1);
proj.Table = table(proj.sensorId, proj.EO, proj.epsilonPerpMv, proj.SgMvPerMm, ...
    proj.deltaGapProjMm, proj.deltaGapLsMm, proj.deltaGapLocMm, proj.deltaGapLimitMm, ...
    proj.usedFallback, proj.clippedMin, proj.clippedMax, ...
    'VariableNames', {'sensorId','EO','epsilonPerpMv','SgMvPerMm','deltaGapProjMm', ...
    'deltaGapLsMm','deltaGapLocMm','deltaGapLimitMm','usedFallback','clippedMin','clippedMax'});
end

function yPerp = orthogonal_residual_local(y, J)
y = y(:);
if isempty(J)
    yPerp = y;
    return;
end
keep = all(isfinite(J), 2) & isfinite(y);
Jk = J(keep, :);
yk = y(keep);
colKeep = vecnorm(Jk, 2, 1) > eps;
Jk = Jk(:, colKeep);
yPerp = y;
if isempty(Jk)
    yPerp(keep) = yk;
else
    yPerp(keep) = yk - Jk * (pinv(Jk) * yk);
end
end

function score = bounded_objective_local(thetaRaw, eo, bundle, cfg, modeName)
theta = bound_params_local(thetaRaw, cfg, modeName, numel(bundle.sensorIds));
detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName);
score = detail.score + 1e4 * sum((thetaRaw(:) - theta(:)).^2);
end

function theta0 = initial_theta_local(A, phi, dx, modeName, nSensor)
theta0 = [A, phi, dx];
theta0 = [theta0, zeros(1, nSensor)];
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
theta(3) = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
need = 3 + nSensor + mode_extra_count_local(modeName) * nSensor;
if numel(theta) < need
    theta(numel(theta)+1:need) = 0;
end
theta = theta(1:need);
[idxEta, idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
if ~isempty(idxEta)
    etaLimit = cfg.etaFallbackLimitMm * ones(1, nSensor);
    if isfield(cfg, 'sensorEtaLimitMm') && numel(cfg.sensorEtaLimitMm) >= nSensor
        etaLimit = cfg.sensorEtaLimitMm(:).';
    end
    etaLimit = max(etaLimit, eps);
    theta(idxEta) = max(min(theta(idxEta), etaLimit), -etaLimit);
end
if ~isempty(idxDg)
    deltaGapLimit = expand_limit_vector_local(cfg.deltaGapLimitMm, nSensor);
    theta(idxDg) = max(min(theta(idxDg), deltaGapLimit), -deltaGapLimit);
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
dx = theta(3);
nSensor = numel(bundle.sensorIds);
[idxEta, idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
eta = zeros(nSensor, 1);
dg = zeros(nSensor, 1);
dmu = zeros(nSensor, 1);
dtau = zeros(nSensor, 1);
if ~isempty(idxEta), eta = theta(idxEta).'; end
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
    xEval = bundle.X(mask) - dx - u(mask) - eta(is);
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
etaLimit = cfg.etaFallbackLimitMm * ones(nSensor, 1);
if isfield(cfg, 'sensorEtaLimitMm') && numel(cfg.sensorEtaLimitMm) >= nSensor
    etaLimit = cfg.sensorEtaLimitMm(:);
end
etaPenalty = pointCount * (cfg.etaRegWeightMv * sqrt(sum((eta ./ max(etaLimit, eps)).^2) / max(nSensor, 1))) ^ 2;
deltaGapLimit = expand_limit_vector_local(cfg.deltaGapLimitMm, nSensor).';
regularizationPenalty = pointCount * (cfg.staticRegWeightMv * sqrt(( ...
    sum((dg ./ max(deltaGapLimit, eps)).^2) + ...
    sum((dmu ./ max(cfg.deltaMuLimit, eps)).^2) + ...
    sum((dtau ./ max(cfg.deltaTauLimitMm, eps)).^2)) / max(nSensor, 1))) ^ 2;
invalidPenalty = 1e12 * nnz(~valid);
objective = residualObj + overshootPenalty + regularizationPenalty + etaPenalty + invalidPenalty;
weightedRmse = sqrt(objective / pointCount);
penaltyObjective = overshootPenalty + regularizationPenalty + etaPenalty + invalidPenalty;
if isfinite(objective) && objective > 0
    penaltyShare = penaltyObjective / objective;
else
    penaltyShare = NaN;
end
detail = struct();
detail.EO = eo;
detail.freqHz = eo * bundle.rotFreqMeanHz;
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx;
detail.sensorEtaMm = eta(:);
detail.deltaGapMm = dg(:);
detail.deltaMuGapPerXMm = dmu(:);
detail.deltaTauMm = dtau(:);
detail.weightedRmseMv = weightedRmse;
detail.plainRmseMv = plainRmse;
detail.residualObjectiveMv2 = residualObj;
detail.overshootPenalty = overshootPenalty;
detail.staticRegPenalty = regularizationPenalty;
detail.etaRegPenalty = etaPenalty;
detail.invalidPenalty = invalidPenalty;
detail.penaltyObjectiveMv2 = penaltyObjective;
detail.penaltyShare = penaltyShare;
detail.score = objective;
detail.VPred = vPred;
detail.uMm = u;
end

function n = mode_extra_count_local(modeName)
n = double(mode_has_gap_local(modeName)) + double(mode_has_tilt_local(modeName)) + ...
    double(mode_has_shift_local(modeName));
end

function tf = mode_has_gap_local(modeName)
tf = any(strcmpi(modeName, {'gap_only','gap_tilt'}));
end

function tf = mode_has_tilt_local(modeName)
tf = strcmpi(modeName, 'gap_tilt');
end

function tf = mode_has_shift_local(~)
tf = false;
end

function [idxEta, idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor)
p = 3;
idxEta = p + (1:nSensor);
p = p + nSensor;
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

function etaLimit = collect_sensor_eta_limit_from_bundle_local(bundle)
nSensor = numel(bundle.sensorIds);
etaLimit = nan(nSensor, 1);
for is = 1:nSensor
    if isfield(bundle, 'sensorInfo') && numel(bundle.sensorInfo) >= is && ...
            isfield(bundle.sensorInfo(is), 'etaLimitMm') && isfinite(bundle.sensorInfo(is).etaLimitMm)
        etaLimit(is) = bundle.sensorInfo(is).etaLimitMm;
    end
end
etaLimit(~isfinite(etaLimit)) = 0.06;
etaLimit = max(etaLimit, 1e-3);
end

function value = resolve_sensor_eta_limit_local(Tpl, corr, cfg)
value = resolve_sensor_eta_stat_local(Tpl, corr, 'etaLimitMm');
if ~isfinite(value)
    value = cfg.etaFallbackLimitMm;
end
value = max(value, 1e-3);
end

function value = resolve_sensor_eta_stat_local(Tpl, corr, fieldName)
value = NaN;
if isstruct(corr) && isfield(corr, fieldName)
    value = first_finite_numeric_local(corr.(fieldName));
end
if ~isfinite(value)
    tplField = lower_camel_to_template_field_local(fieldName);
    if isstruct(Tpl) && isfield(Tpl, tplField)
        value = first_finite_numeric_local(Tpl.(tplField));
    end
end
end

function fieldName = lower_camel_to_template_field_local(name)
switch name
    case 'etaLimitMm'
        fieldName = 'eta_limit_mm';
    case 'etaMedianMm'
        fieldName = 'eta_median_mm';
    case 'etaIqrMm'
        fieldName = 'eta_iqr_mm';
    otherwise
        fieldName = name;
end
end

function value = first_finite_numeric_local(raw)
value = NaN;
if isempty(raw) || ~isnumeric(raw)
    return;
end
raw = raw(:);
raw = raw(isfinite(raw));
if ~isempty(raw)
    value = raw(1);
end
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

function mask = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
end

function mask = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
bestSeg = 1;
bestPeak = -inf;
for i = 1:numel(starts)
    seg = idx(starts(i):ends(i));
    peakVal = max(v(seg));
    if peakVal > bestPeak
        bestPeak = peakVal;
        bestSeg = i;
    end
end
mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end

function mask = build_dynamic_effective_mask_local(x, vMv, Tpl, thresholdMv, cfg)
finite = isfinite(x) & isfinite(vMv);
gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:) * 1000, x(:), 'linear', 0));
end
gMax = max(gTpl(finite), [], 'omitnan');
if isfinite(gMax) && gMax > 0
    maskGradient = gTpl >= cfg.dynamicGradientMinRatio * gMax;
else
    maskGradient = false(size(x));
end
maskPeak = vMv >= thresholdMv;
mask = finite & (maskGradient | maskPeak);
if nnz(mask) < 8
    mask = finite;
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

function eoCandidates = build_eo_candidates_local(rotFreqHz, cfg)
if strcmpi(cfg.eoSearchMode, 'region_order')
    eoCenter = max(1, round(cfg.eoSearchCenterOrder));
    halfWidth = max(0, floor(cfg.eoSearchHalfWidth));
    eoCandidates = max(1, eoCenter - halfWidth):(eoCenter + halfWidth);
    return;
end
rotFreqHz = max(rotFreqHz, eps);
freqSearchHz = cfg.freqSearchHz;
eoMin = max(1, ceil(min(freqSearchHz) / rotFreqHz));
eoMax = max(eoMin, floor(max(freqSearchHz) / rotFreqHz));
eoCandidates = eoMin:eoMax;
if isempty(eoCandidates)
    eoCenter = max(1, round(mean(freqSearchHz) / rotFreqHz));
    eoCandidates = max(1, eoCenter - cfg.eoPad):max(1, eoCenter + cfg.eoPad);
end
freq = eoCandidates .* rotFreqHz;
eoCandidates = eoCandidates(freq >= min(freqSearchHz) & freq <= max(freqSearchHz));
end

function wr = pack_window_result_local(iw, Wmap, bundle, seedTable, mainCandidateEO, selectionInfo, modelFits)
wr = struct();
wr.windowId = iw;
wr.lapRange = Wmap.lap_range;
wr.timeWindow = Wmap.time_window;
wr.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
wr.rotRpmMean = Wmap.rot_rpm_mean;
wr.bundle = bundle;
wr.VPSeedTable = seedTable;
wr.VPSelectedEO = selectionInfo.vpEO;
wr.MainCandidateEO = mainCandidateEO;
wr.EOSelectionInfo = selectionInfo;
wr.modelFits = modelFits;
end

function limit = expand_limit_vector_local(value, n)
if isscalar(value)
    limit = repmat(value, 1, n);
else
    limit = value(:).';
    if numel(limit) < n
        limit(end+1:n) = limit(end);
    end
    limit = limit(1:n);
end
limit = max(limit, eps);
end

function row = make_trend_row_local(wr, directTrend, cfg)
directEO = NaN; directFreq = NaN; directAmp = NaN; directRmseMv = NaN;
if ~isempty(directTrend) && height(directTrend) >= wr.windowId
    directEO = directTrend.EO_id(wr.windowId);
    directFreq = directTrend.fn_id(wr.windowId);
    directAmp = directTrend.A_id(wr.windowId);
    directRmseMv = 1000 * directTrend.weighted_voltage_rmse(wr.windowId);
end
selectionInfo = wr.EOSelectionInfo;
selectionMode = string(selectionInfo.mode);
selectedEOText = string(mat2str(selectionInfo.selectedEO));
vpEOText = string(mat2str(selectionInfo.vpEO));
directCandidateEO = selectionInfo.directEO;
prevCandidateEO = selectionInfo.prevEO;
prevReferenceSource = string(selectionInfo.prevSource);
prevReferenceReason = string(selectionInfo.prevReason);
usedPrevWindowReference = logical(selectionInfo.usedPrevWindow);
phaseSafeReferenceEO = selectionInfo.phaseSafeReferenceEO;
phaseSafeReferenceSource = string(selectionInfo.phaseSafeReferenceSource);
phaseSafeReferenceReason = string(selectionInfo.phaseSafeReferenceReason);
usedPhaseSafeExpansion = logical(selectionInfo.usedPhaseSafeExpansion);
phaseSafeCorePointCount = selectionInfo.phaseSafeCorePointCount;
phaseSafeExpandedPointCount = selectionInfo.phaseSafeSelectedPointCount;
phaseSafeAddedPointCount = phaseSafeExpandedPointCount - phaseSafeCorePointCount;
gapProj = selectionInfo.DeltaGapProjection;
F = wr.modelFits;
row = table(wr.windowId, wr.lapRange(1), wr.lapRange(end), wr.rotFreqMeanHz, ...
    selectionMode, selectedEOText, vpEOText, directCandidateEO, ...
    prevCandidateEO, prevReferenceSource, prevReferenceReason, usedPrevWindowReference, ...
    phaseSafeReferenceEO, phaseSafeReferenceSource, phaseSafeReferenceReason, ...
    usedPhaseSafeExpansion, phaseSafeCorePointCount, phaseSafeExpandedPointCount, phaseSafeAddedPointCount, ...
    'VariableNames', {'window_id','lap_start','lap_end','rot_freq_hz', ...
    'eo_candidate_mode','selected_eo','vp_candidate_eo','direct_candidate_eo', ...
    'prev_candidate_eo','prev_reference_source','prev_reference_reason','used_prev_window_reference', ...
    'phase_safe_reference_eo','phase_safe_reference_source','phase_safe_reference_reason', ...
    'used_phase_safe_expansion','core_point_count','phase_safe_point_count','phase_safe_added_point_count'});
if strcmpi(cfg.runMode, 'comparison')
    row.direct_EO = directEO;
    row.direct_frequency_hz = directFreq;
    row.direct_amplitude_mm = directAmp;
    row.direct_rmse_mV = directRmseMv;
end
if isfield(F, 'fixed')
    row.fixed_EO = F.fixed.EO;
    row.fixed_frequency_hz = F.fixed.freqHz;
    row.fixed_amplitude_mm = F.fixed.amplitudeMm;
    row.fixed_dx_mm = F.fixed.dxMm;
    row.fixed_rmse_mV = F.fixed.weightedRmseMv;
end
if isfield(F, 'gap_only')
    row.gap_EO = F.gap_only.EO;
    row.gap_frequency_hz = F.gap_only.freqHz;
    row.gap_amplitude_mm = F.gap_only.amplitudeMm;
    row.gap_dx_mm = F.gap_only.dxMm;
    row.gap_rmse_mV = F.gap_only.weightedRmseMv;
    row.gap_mean_delta_gap_mm = mean(F.gap_only.deltaGapMm, 'omitnan');
end
if isfield(F, 'gap_tilt')
    row.tilt_EO = F.gap_tilt.EO;
    row.tilt_frequency_hz = F.gap_tilt.freqHz;
    row.tilt_amplitude_mm = F.gap_tilt.amplitudeMm;
    row.tilt_dx_mm = F.gap_tilt.dxMm;
    row.tilt_rmse_mV = F.gap_tilt.weightedRmseMv;
    row.tilt_mean_delta_gap_mm = mean(F.gap_tilt.deltaGapMm, 'omitnan');
    row.tilt_mean_delta_mu = mean(F.gap_tilt.deltaMuGapPerXMm, 'omitnan');
    row.tilt_mean_eta_mm = mean(F.gap_tilt.sensorEtaMm, 'omitnan');
    row.tilt_max_abs_eta_mm = max_abs_finite_local(F.gap_tilt.sensorEtaMm);
    row.delta_g_proj_alpha = gapProj.alpha;
    row.delta_g_proj_mean_mm = mean(gapProj.deltaGapProjMm, 'omitnan');
    row.delta_g_proj_max_mm = max(gapProj.deltaGapProjMm, [], 'omitnan');
    row.delta_g_loc_mean_mm = mean(gapProj.deltaGapLocMm, 'omitnan');
    row.delta_g_limit_mean_mm = mean(gapProj.deltaGapLimitMm, 'omitnan');
    row.delta_g_limit_max_mm = max(gapProj.deltaGapLimitMm, [], 'omitnan');
    row.epsilon_perp_mean_mV = mean(gapProj.epsilonPerpMv, 'omitnan');
    row.S_g_min_mV_per_mm = min(gapProj.SgMvPerMm, [], 'omitnan');
    row.delta_g_proj_fallback_count = sum(gapProj.usedFallback);
    row.delta_g_proj_cap_count = sum(gapProj.clippedMax);
end
end

function Summary = build_summary_table_local(Trend, bestIdx, cfg)
if strcmpi(cfg.runMode, 'comparison')
    spec = {
    "direct_low_template_main", "direct_EO", "direct_frequency_hz", "direct_amplitude_mm", "direct_rmse_mV"
    "fixed", "fixed_EO", "fixed_frequency_hz", "fixed_amplitude_mm", "fixed_rmse_mV"
    "gap_only", "gap_EO", "gap_frequency_hz", "gap_amplitude_mm", "gap_rmse_mV"
    "gap_tilt", "tilt_EO", "tilt_frequency_hz", "tilt_amplitude_mm", "tilt_rmse_mV"
    };
else
    spec = {
    "gap_tilt", "tilt_EO", "tilt_frequency_hz", "tilt_amplitude_mm", "tilt_rmse_mV"
    };
end
keep = false(size(spec, 1), 1);
for i = 1:size(spec, 1)
    keep(i) = all(ismember(cellstr(spec(i, 2:5)), Trend.Properties.VariableNames));
end
spec = spec(keep, :);
n = size(spec, 1);
method = strings(n, 1);
dominantEO = NaN(n, 1);
meanFreq = NaN(n, 1);
stdFreq = NaN(n, 1);
meanAmp = NaN(n, 1);
meanRmse = NaN(n, 1);
medianRmse = NaN(n, 1);
for i = 1:n
    method(i) = spec{i, 1};
    dominantEO(i) = mode_finite_local(Trend.(spec{i, 2}));
    meanFreq(i) = mean(Trend.(spec{i, 3}), 'omitnan');
    stdFreq(i) = std(Trend.(spec{i, 3}), 'omitnan');
    meanAmp(i) = mean(Trend.(spec{i, 4}), 'omitnan');
    meanRmse(i) = mean(Trend.(spec{i, 5}), 'omitnan');
    medianRmse(i) = median(Trend.(spec{i, 5}), 'omitnan');
end
bestWindowIndex = repmat(bestIdx, n, 1);
Summary = table(method, dominantEO, meanFreq, stdFreq, meanAmp, meanRmse, medianRmse, bestWindowIndex, ...
    'VariableNames', {'method','dominant_EO','mean_frequency_hz','std_frequency_hz', ...
    'mean_amplitude_mm','mean_rmse_mV','median_rmse_mV','best_free_gap_window'});
end

function plot_trend_local(T, figFile, cfg)
style = paper_style_local();
fig = figure('Name', 'Step07J nested static warp trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 13.0]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on; box on;
if strcmpi(cfg.runMode, 'comparison')
    plot_if_column_local(T, 'direct_EO', style.black, '.', 'Direct');
end
plot_if_column_local(T, 'fixed_EO', style.blue, 'o', 'Fixed gap');
plot_if_column_local(T, 'gap_EO', style.red, 's', 'Gap');
plot_if_column_local(T, 'tilt_EO', style.green, '^', 'Gap+tilt');
ylabel('EO'); title('Order'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
if strcmpi(cfg.runMode, 'comparison')
    plot_if_column_local(T, 'direct_amplitude_mm', style.black, '.', 'Direct');
end
plot_if_column_local(T, 'fixed_amplitude_mm', style.blue, 'o', 'Fixed gap');
plot_if_column_local(T, 'gap_amplitude_mm', style.red, 's', 'Gap');
plot_if_column_local(T, 'tilt_amplitude_mm', style.green, '^', 'Gap+tilt');
ylabel('A (mm)', 'Interpreter', 'tex'); title('Amplitude'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
if strcmpi(cfg.runMode, 'comparison')
    plot_if_column_local(T, 'direct_rmse_mV', style.black, '.', 'Direct');
end
plot_if_column_local(T, 'fixed_rmse_mV', style.blue, 'o', 'Fixed gap');
plot_if_column_local(T, 'gap_rmse_mV', style.red, 's', 'Gap');
plot_if_column_local(T, 'tilt_rmse_mV', style.green, '^', 'Gap+tilt');
xlabel('Window'); ylabel('RMSE (mV)'); title('Voltage residual'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function plot_if_column_local(T, columnName, colorValue, markerValue, displayName)
if ~ismember(columnName, T.Properties.VariableNames)
    return;
end
plot(T.window_id, T.(columnName), '-', 'Color', colorValue, 'LineWidth', 1.0, ...
    'Marker', markerValue, 'MarkerSize', 4.2, 'DisplayName', displayName);
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
    xRaw = bundle.X(mask) - fit.dxMm - fit.sensorEtaMm(is);
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
title(sprintf('Gap increment, %s main model', cfg.mainModel), 'Interpreter', 'none');
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
title(sprintf('High-speed tilt increment, %s main model', cfg.mainModel), 'Interpreter', 'none');
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
title(sprintf('Sensor x-shift increment, %s main model', cfg.mainModel), 'Interpreter', 'none');
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

function [freqRange, hasOverride] = parse_frequency_range_override_local(defaultRange)
rawMin = strtrim(getenv('STEP07J_FREQ_MIN_HZ'));
rawMax = strtrim(getenv('STEP07J_FREQ_MAX_HZ'));
hasOverride = ~isempty(rawMin) || ~isempty(rawMax);
freqRange = defaultRange;
if ~hasOverride
    return;
end
if isempty(rawMin) || isempty(rawMax)
    error('Set both STEP07J_FREQ_MIN_HZ and STEP07J_FREQ_MAX_HZ, or neither.');
end
fMin = str2double(rawMin);
fMax = str2double(rawMax);
if ~isfinite(fMin) || ~isfinite(fMax) || fMin <= 0 || fMax <= fMin
    error('STEP07J_FREQ_MIN_HZ/MAX_HZ must define a positive increasing frequency range.');
end
freqRange = [fMin, fMax];
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

function value = max_abs_finite_local(x)
x = abs(x(:));
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end

