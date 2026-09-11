%% Step07J: Nested high-speed static warp models
% The main gap_fixed_tilt route jointly identifies EO, vibration and static
% gap increment inside the direct-template + VP Top-K EO candidate set with
% the full waveform objective. The low-speed calibrated sensor tilt slope is
% kept fixed in the corrected gap library; no residual high-speed dmu is
% identified in the paper-facing main route. VP narrows the EO candidates
% but does not fix the final EO. The default main route computes only the
% paper-facing gap_fixed_tilt model. In comparison mode this script also
% computes nested static correction models on top of the low-speed template:
%
%   fixed             : direct-EO baseline with no clearance increment
%   gap_only          : dg_s with low-speed calibrated mu fixed (gap_fixed_tilt)
%   gap_tilt          : dg_s + residual dmu_s*(x-tau_s)
%
% Paper-facing main output uses gap_only/gap_fixed_tilt. gap_tilt is retained
% only in the separate comparison run to quantify sensitivity to a residual
% high-speed tilt increment. gap_tilt_shift is excluded from the current
% comparison because the extra sensor-shift freedom is invalid for the
% paper-facing method.
%
% This is the current paper-facing identification route. Each three-lap
% sliding window independently selects an EO from the VP Top-K set and then
% refines that candidate frequency locally. No frequency, phase or EO is
% shared across a multi-window region.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
packageCfg = Config_20241106();
rootDir = packageCfg.paths.root;
addpath(thisDir);
flowCfg = ProjectionFlow_Config_20241106();
flowCfg = apply_case_overrides_to_flow_config_local(flowCfg);
outDir = packageCfg.paths.gapRuntime;
resultOutDir = strtrim(getenv('STEP07J_RESULT_OUTPUT_DIR'));
if isempty(resultOutDir)
    resultOutDir = packageCfg.paths.gapResults;
end
if exist(resultOutDir, 'dir') ~= 7
    mkdir(resultOutDir);
end
figDir = fullfile(resultOutDir, 'figures_step07j_nested_static_warp_vp_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = flowCfg.identification.targetBlade;
cfg.analysisSensors = flowCfg.identification.analysisSensors;
cfg.gapSensors = flowCfg.identification.gapSensors;
cfg.analysisStartTimeSec = flowCfg.identification.analysisStartTimeSec;
cfg.windowBladePasses = flowCfg.identification.windowBladePasses;
cfg.slidingStepBladePasses = flowCfg.identification.slidingStepBladePasses;
cfg.directOnlySensors = setdiff(cfg.analysisSensors(:).', cfg.gapSensors(:).', 'stable');
cfg.freqSearchHz = flowCfg.identification.freqSearchHz;
cfg.eoPad = 2;
cfg.vpTopK = parse_positive_integer_env_local('STEP07J_TOP_K_EO', 3);
cfg.maxWindowsConfig = inf;
if isfield(flowCfg, 'identification') && isfield(flowCfg.identification, 'maxRunWindows')
    cfg.maxWindowsConfig = flowCfg.identification.maxRunWindows;
end
cfg.maxWindowsEnvIgnored = strtrim(getenv('STEP07J_MAX_WINDOWS'));
cfg.maxWindows = cfg.maxWindowsConfig;
if ~(isscalar(cfg.maxWindows) && isnumeric(cfg.maxWindows)) || ...
        (~isfinite(cfg.maxWindows) && ~isinf(cfg.maxWindows)) || cfg.maxWindows <= 0
    error('ProjectionFlow_Config_20241106 identification.maxRunWindows must be positive or inf.');
end
if isfinite(cfg.maxWindows)
    cfg.maxWindows = max(1, floor(cfg.maxWindows));
end
cfg.maxPointsPerWindow = parse_positive_integer_env_local('STEP07J_MAX_POINTS_PER_WINDOW', inf);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP07J_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DX_LIMIT_MM', 0.35);
% Diagnostic-only controls.  The normal paper-facing route leaves both
% switches inactive: EO remains selected from VP/full-wave screening and
% the matching Foundation Step05 fixed result remains the fixed baseline.
if ~isempty(strtrim(getenv('STEP07J_FORCED_EO')))
    error('Final 20241106 Step07J does not allow forced EO. Clear STEP07J_FORCED_EO.');
end
cfg.forcedEO = 0;
cfg.skipFoundationFixedResult = false;
cfg.keepConfiguredStaticEta = true;
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_LIMIT_MM', flowCfg.model.deltaGapLimitMm);
% Low-speed calibration already carries the fixed per-sensor tilt slope.
% Residual dmu is tightly bounded and used only in comparison mode.
cfg.deltaMuLimit = parse_nonnegative_numeric_env_local('STEP07J_DELTA_MU_LIMIT', 0.005);
cfg.deltaTauLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_TAU_LIMIT_MM', flowCfg.model.deltaTauLimitMm);
cfg.etaFallbackLimitMm = parse_nonnegative_numeric_env_local('STEP07J_ETA_FALLBACK_LIMIT_MM', 0.06);
cfg.etaRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_ETA_REG_WEIGHT_MV', 0.35);
cfg.staticRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_STATIC_REG_WEIGHT_MV', 0.75);
cfg.vpSeedObservationMode = 'gradient_displacement';
cfg.vpGradientMinRatio = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_MIN_RATIO', 0.10);
cfg.vpGradientReferenceQuantile = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_REFERENCE_QUANTILE', 95);
cfg.vpGradientWeightPower = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_WEIGHT_POWER', 2);
cfg.finalObjective = lower(strtrim(getenv('STEP07J_FINAL_OBJECTIVE')));
if isempty(cfg.finalObjective)
    cfg.finalObjective = 'plain';
end
if ~ismember(cfg.finalObjective, {'plain','weighted'})
    error('STEP07J_FINAL_OBJECTIVE must be "plain" or "weighted".');
end
% Every retained EO receives the same local continuous-frequency refinement
% before final selection. Set the environment switch to 0 only for the
% integer-frequency ablation.
cfg.continuousFrequencyRefine = parse_logical_env_local( ...
    'STEP07J_REFINE_EACH_CANDIDATE_FREQUENCY', true);
cfg.frequencyRefineHalfWidthHz = parse_nonnegative_numeric_env_local( ...
    'STEP07J_FREQ_REFINE_HALF_WIDTH_HZ', 2.0);
cfg.lowConfidenceEoMarginPercent = parse_nonnegative_numeric_env_local( ...
    'STEP07J_LOW_CONFIDENCE_EO_MARGIN_PERCENT', 0.5);
cfg.activeEtaModel = 'theoretical_zero_eta';
cfg.sensorEtaMode = 'fixed_zero';
cfg.fitSensorEta = false;
cfg.noGapKernelMode = 'foundation_core';
cfg.gapUpdateMode = 'vp_topk_delayed_gap_forward';
cfg.lockGapModelsToFixedEO = parse_logical_env_local( ...
    'STEP07J_LOCK_GAP_EO_TO_FIXED', false);
% A weighted fixed-model floor cannot constrain an unweighted final
% objective. In plain-RMSE mode every Top-K candidate is optimized and
% ranked directly by the same final voltage objective.
cfg.initializeGapModelsFromFixed = ~strcmpi(cfg.finalObjective, 'plain');
cfg.mainModel = 'gap_only';
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
    cfg.modelNames = {'gap_only'};
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
    cfg.prevWindowCandidateMode = 'off';
end
if ~ismember(cfg.prevWindowCandidateMode, {'off', 'soft'})
    error('STEP07J_PREV_WINDOW_CANDIDATE_MODE must be "off" or "soft".');
end
cfg.prevWindowRefreshEvery = floor(parse_nonnegative_numeric_env_local( ...
    'STEP07J_PREV_WINDOW_REFRESH_EVERY', 5));
cfg.phaseSafeExpansion = parse_logical_env_local('STEP07J_PHASE_SAFE_EXPANSION', true);
cfg.phaseSafeMarginMm = parse_nonnegative_numeric_env_local('STEP07J_PHASE_SAFE_MARGIN_MM', 0.02);
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
cfg.fixedEoMode = lower(strtrim(getenv('STEP07J_FIXED_EO_MODE')));
if isempty(cfg.fixedEoMode)
    cfg.fixedEoMode = 'selected';
end
if ~ismember(cfg.fixedEoMode, {'direct', 'selected'})
    error('STEP07J_FIXED_EO_MODE must be "direct" or "selected".');
end
cfg.useLegacyDirectFixed = parse_logical_env_local('STEP07J_USE_LEGACY_DIRECT_FIXED', false);
cfg.alignDirectOnlyXCenter = parse_logical_env_local('STEP07J_ALIGN_DIRECT_ONLY_XCENTER', true);
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
cfg.domainSoftMarginMm = parse_nonnegative_numeric_env_local('STEP07J_DOMAIN_SOFT_MARGIN_MM', 0.10);
cfg.domainMarginMm = parse_nonnegative_numeric_env_local('STEP07J_DOMAIN_MARGIN_MM', 0.02);
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
cfg.sensorQueryGuardScaleTable = [5 0.80; 7 0.75];
cfg.localBundleSource = lower(strtrim(getenv('STEP07J_LOCAL_BUNDLE_SOURCE')));
if isempty(cfg.localBundleSource)
    cfg.localBundleSource = 'foundation_step05_bundle';
end
if ~ismember(cfg.localBundleSource, {'compact_case', 'dynamic_map', 'foundation_step05_bundle'})
    error('STEP07J_LOCAL_BUNDLE_SOURCE must be "compact_case", "dynamic_map", or "foundation_step05_bundle".');
end
cfg.allowDynamicMapFallback = parse_logical_env_local('STEP07J_ALLOW_DYNAMIC_MAP_FALLBACK', false);
cfg.allowFoundationTimeMismatch = parse_logical_env_local( ...
    'STEP07J_ALLOW_FOUNDATION_TIME_MISMATCH', false);
if strcmpi(cfg.localBundleSource, 'dynamic_map') && ~cfg.allowDynamicMapFallback
    error(['STEP07J_LOCAL_BUNDLE_SOURCE=dynamic_map is the current 20241106 foundation-core route. ' ...
        'Keep STEP07J_ALLOW_DYNAMIC_MAP_FALLBACK=1, or set STEP07J_LOCAL_BUNDLE_SOURCE=compact_case with a valid compact waveform case.']);
end
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
gapSensorOverride = strtrim(getenv('STEP07J_GAP_SENSORS'));
if ~isempty(gapSensorOverride)
    parsedGapSensors = sscanf(gapSensorOverride, '%d').';
    if isempty(parsedGapSensors)
        error('STEP07J_GAP_SENSORS must contain integer sensor IDs, for example "5 7".');
    end
    cfg.gapSensors = parsedGapSensors;
end
cfg.gapSensors = intersect(cfg.gapSensors(:).', cfg.analysisSensors(:).', 'stable');
assert_capacitive_gap_sensors_local(cfg.gapSensors, 'STEP07J_ALLOW_NONCAP_GAP_SENSORS');
cfg.directOnlySensors = setdiff(cfg.analysisSensors(:).', cfg.gapSensors(:).', 'stable');
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
gapSensorTag = ['S', sprintf('%d', cfg.gapSensors)];
caseTag = sprintf('B%d_%s', cfg.targetBlade, sensorTag);
staticEtaOverride = strtrim(getenv('STEP07J_STATIC_ETA_FILE'));
if ~isempty(staticEtaOverride)
    error('Final 20241106 Step07J enforces eta=0; STEP07J_STATIC_ETA_FILE is no longer accepted.');
end
cfg.fixedSensorEtaMm = zeros(1, numel(cfg.analysisSensors));
cfg.fixedSensorEtaInfo = struct('file', 'built_in_theoretical_zero', ...
    'source', 'theoretical_zero');
resultSuffix = sanitize_suffix_local(strtrim(getenv('STEP07J_RESULT_SUFFIX')));
if isempty(resultSuffix)
    if strcmpi(cfg.runMode, 'comparison')
        resultSuffix = '_method_compare';
    else
        if cfg.continuousFrequencyRefine
            frequencyTag = sprintf('_DFpm%gHz', cfg.frequencyRefineHalfWidthHz);
        else
            frequencyTag = '_integerFrequency';
        end
        resultSuffix = sprintf('_main_gapfixedtilt_W%dS%d_F%dto%dHz%s', ...
            flowCfg.identification.windowBladePasses, ...
            flowCfg.identification.slidingStepBladePasses, ...
            round(cfg.freqSearchHz(1)), round(cfg.freqSearchHz(2)), frequencyTag);
    end
end

rotDir = fullfile(rootDir, 'archive', 'unused_rotating_calibration');
rotTemplateDir = fullfile(rotDir, 'output', 'templates');
rotDynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');
rotResultDir = fullfile(rotDir, 'output', 'new_flow', '06_identification');
gapTemplateDir = outDir;
gapDynamicDir = outDir;
templateBankFile = fullfile(packageCfg.paths.calibrationGap, ...
    flowCfg.calibration.templateBankFile);
gapBankFile = fullfile(packageCfg.paths.calibrationGap, ...
    flowCfg.calibration.gapBankFile);

correctedLibOverride = strtrim(getenv('STEP07J_CORRECTED_LIB_FILE'));
if ~isempty(correctedLibOverride)
    correctedLibFile = correctedLibOverride;
else
    correctedLibFile = find_corrected_gap_library_file_local(outDir, cfg.targetBlade, gapSensorTag, gapBankFile);
end
templateOverride = strtrim(getenv('STEP07J_TEMPLATE_FILE'));
if ~isempty(templateOverride)
    templateFile = templateOverride;
else
    templateFile = '';
end
dynamicOverride = strtrim(getenv('STEP07J_DYNAMIC_MAP_FILE'));
if ~isempty(dynamicOverride)
    if exist(dynamicOverride, 'file') ~= 2
        error('STEP07J_DYNAMIC_MAP_FILE does not exist: %s', dynamicOverride);
    end
    dynamicFile = dynamicOverride;
else
    dynamicFile = find_dynamic_map_file_local({gapDynamicDir, rotDynamicDir}, ...
        cfg.targetBlade, cfg.analysisSensors, flowCfg.identification.analysisStartTimeSec, ...
        flowCfg.identification.windowBladePasses, ...
        flowCfg.identification.slidingStepBladePasses);
end
directOverride = strtrim(getenv('STEP07J_DIRECT_RESULT_FILE'));
if ~isempty(directOverride)
    if exist(directOverride, 'file') ~= 2
        error('STEP07J_DIRECT_RESULT_FILE does not exist: %s', directOverride);
    end
    directMainFile = directOverride;
else
    directMainFile = find_direct_main_result_file_local(rotResultDir, cfg.targetBlade, ...
        cfg.analysisSensors, flowCfg.identification.analysisStartTimeSec);
end
compactOverride = strtrim(getenv('STEP07J_COMPACT_WAVEFORM_CASE_FILE'));
if ~isempty(compactOverride)
    if exist(compactOverride, 'file') ~= 2
        error('STEP07J_COMPACT_WAVEFORM_CASE_FILE does not exist: %s', compactOverride);
    end
    compactWaveformCaseFile = compactOverride;
else
    compactWaveformCaseFile = find_compact_waveform_case_file_local( ...
        rotDir, cfg.targetBlade, cfg.analysisSensors, flowCfg.identification.analysisStartTimeSec);
end
foundationStep05File = '';
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    foundationStep05Override = strtrim(getenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE'));
    if ~isempty(foundationStep05Override)
        if exist(foundationStep05Override, 'file') ~= 2
            error('STEP07J_FOUNDATION_STEP05_RESULT_FILE does not exist: %s', foundationStep05Override);
        end
        foundationStep05File = foundationStep05Override;
    else
        foundationStep05File = find_foundation_step05_result_file_local( ...
            rootDir, cfg.targetBlade, cfg.analysisSensors, ...
            flowCfg.identification.analysisStartTimeSec, ...
            flowCfg.identification.windowBladePasses, ...
            flowCfg.identification.slidingStepBladePasses, cfg.freqSearchHz);
        if isempty(foundationStep05File)
            error(['No exact Foundation Step05 result exists for T=%.6f s, W%dS%d, F%d-%d Hz. ' ...
                'Run Step00_Run_ProjectionFlow_20241106.m to build the matched source.'], ...
                flowCfg.identification.analysisStartTimeSec, ...
                flowCfg.identification.windowBladePasses, ...
                flowCfg.identification.slidingStepBladePasses, ...
                round(cfg.freqSearchHz(1)), round(cfg.freqSearchHz(2)));
        end
    end
end

if ~isfile(correctedLibFile)
    error('Run Step06I first. Missing file: %s', correctedLibFile);
end

Sd = load(dynamicFile, 'DynamicMap');
% The selected DynamicMap is authoritative for the analyzed time segment.
% This also keeps metadata correct when an explicit tagged map is supplied.
if isfield(Sd.DynamicMap, 'SourceSettings') && isstruct(Sd.DynamicMap.SourceSettings)
    sourceSettings = Sd.DynamicMap.SourceSettings;
    if isfield(sourceSettings, 'analysisStartTimeSec')
        flowCfg.identification.analysisStartTimeSec = sourceSettings.analysisStartTimeSec;
    end
    if isfield(sourceSettings, 'targetLaps')
        flowCfg.identification.targetBladePasses = sourceSettings.targetLaps;
    end
    if isfield(sourceSettings, 'windowLaps') && ...
            sourceSettings.windowLaps ~= flowCfg.identification.windowBladePasses
        error(['DynamicMap uses %d-lap windows, but the active main configuration requires %d laps. ' ...
            'Rerun Step06_BuildGapAwareDynamicMap_20241106.m.'], ...
            sourceSettings.windowLaps, flowCfg.identification.windowBladePasses);
    end
    if isfield(sourceSettings, 'slidingStepLaps') && ...
            sourceSettings.slidingStepLaps ~= flowCfg.identification.slidingStepBladePasses
        error(['DynamicMap uses a %d-lap step, but the active main configuration requires %d. ' ...
            'Rerun Step06_BuildGapAwareDynamicMap_20241106.m.'], ...
            sourceSettings.slidingStepLaps, flowCfg.identification.slidingStepBladePasses);
    end
end
if isempty(Sd.DynamicMap.Window) || ~isfield(Sd.DynamicMap.Window, 'lap_range') || ...
        any(arrayfun(@(w) numel(w.lap_range), Sd.DynamicMap.Window) ~= ...
        flowCfg.identification.windowBladePasses)
    error(['DynamicMap window contents do not match the requested %d-lap main window. ' ...
        'Rerun Step06_BuildGapAwareDynamicMap_20241106.m.'], ...
        flowCfg.identification.windowBladePasses);
end
FoundationStep05Source = [];
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    FoundationStep05Source = load_foundation_step05_source_local(foundationStep05File, cfg);
end
[foundationTimeMatch, foundationTimeCheck] = foundation_time_matches_dynamic_map_local( ...
    FoundationStep05Source, Sd.DynamicMap);
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle') && ~foundationTimeMatch
    if cfg.allowFoundationTimeMismatch
        warning(['Foundation Step05 bundle time span %s does not match DynamicMap span %s. ' ...
            'Keeping the Foundation bundle because STEP07J_ALLOW_FOUNDATION_TIME_MISMATCH=1.'], ...
            mat2str(foundationTimeCheck.foundation_span_s, 9), ...
            mat2str(foundationTimeCheck.dynamic_span_s, 9));
        foundationTimeCheck.action = 'override_keep_foundation_bundle';
    elseif cfg.allowDynamicMapFallback
        warning(['Foundation Step05 bundle time span %s does not match DynamicMap span %s. ' ...
            'Using the current-time DynamicMap waveform bundle instead.'], ...
            mat2str(foundationTimeCheck.foundation_span_s, 9), ...
            mat2str(foundationTimeCheck.dynamic_span_s, 9));
        cfg.localBundleSource = 'dynamic_map';
        foundationTimeCheck.action = 'fallback_to_dynamic_map';
    else
        error(['Foundation Step05 bundle and DynamicMap use different time spans. ' ...
            'Rebuild the Step05 source for the requested time or enable DynamicMap fallback.']);
    end
end
Template = [];
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle') && ...
        ~isempty(FoundationStep05Source) && isfield(FoundationStep05Source, 'Template') && ...
        isstruct(FoundationStep05Source.Template)
    Template = filter_template_sensors_local(FoundationStep05Source.Template, cfg.analysisSensors, sensorTag);
    templateFile = '[foundation Step05 Result.Template]';
elseif strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
        rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag);
elseif isempty(templateOverride) && isfile(templateBankFile)
    Template = load_template_from_bank_local(templateBankFile, cfg.targetBlade, cfg.analysisSensors, sensorTag);
    templateFile = templateBankFile;
end
CorrectedGapLibrary = load_corrected_gap_library_for_run_local(correctedLibFile, cfg, Template);
if isempty(Template) && isempty(templateOverride) && ...
        isfield(CorrectedGapLibrary, 'lowSpeedTemplate') && isstruct(CorrectedGapLibrary.lowSpeedTemplate)
    if isempty(cfg.directOnlySensors)
        Template = filter_template_sensors_local(CorrectedGapLibrary.lowSpeedTemplate, cfg.analysisSensors, sensorTag);
    else
        baseTemplateFile = find_low_speed_template_file_local({rotTemplateDir, gapTemplateDir}, cfg.targetBlade, cfg.analysisSensors);
        Template = load_low_speed_template_local(baseTemplateFile, cfg.analysisSensors, sensorTag);
        Template = merge_gap_template_sensors_local(Template, CorrectedGapLibrary.lowSpeedTemplate, cfg.gapSensors, sensorTag);
    end
    templateFile = '[CorrectedGapLibrary.lowSpeedTemplate]';
else
    if isempty(Template) && isempty(templateFile)
        templateFile = find_low_speed_template_file_local({rotTemplateDir, gapTemplateDir}, cfg.targetBlade, cfg.analysisSensors);
    end
    if isempty(Template)
        St = load(templateFile, 'Template');
        Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
    end
end
DynamicMap = filter_dynamic_map_sensors_local(Sd.DynamicMap, cfg.analysisSensors, sensorTag);
if ~strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    Template = align_direct_only_template_xcenter_local(Template, DynamicMap, cfg);
end
Template = sanitize_template_bundle_local(Template);
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    coordinateCheck = struct('table', table(), 'tolerance_mm', 1e-6, ...
        'max_abs_delta_mm', 0, 'is_consistent', true, ...
        'status', 'foundation_step05_template_bundle_coordinates');
else
    coordinateCheck = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysisSensors, 1e-6);
end
responseSurface = CorrectedGapLibrary.responseSurface;
directRef = load_direct_reference_local(directMainFile);
directTrend = directRef.Trend;
CompactWaveformCase = load_compact_waveform_case_local(compactWaveformCaseFile, cfg);

availableWindows = numel(DynamicMap.Window);
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    if isempty(FoundationStep05Source)
        error(['STEP07J_LOCAL_BUNDLE_SOURCE=foundation_step05_bundle requires a Step05 foundation result. ' ...
            'Set STEP07J_FOUNDATION_STEP05_RESULT_FILE or run the matching btt_data_foundation Step05 first.']);
    end
    availableWindows = min(availableWindows, numel(FoundationStep05Source.Window));
end
numWindows = min(availableWindows, cfg.maxWindows);
fprintf('\n=== Step07J: nested static warp VP/full-wave ===\n');
fprintf('Template: %s\n', templateFile);
fprintf('DynamicMap: %s\n', dynamicFile);
fprintf('Corrected gap library: %s\n', correctedLibFile);
fprintf('Active case: B%d, sensors %s, gap sensors %s, start %.6f s, laps %d, window %d, step %d.\n', ...
    cfg.targetBlade, mat2str(cfg.analysisSensors), mat2str(cfg.gapSensors), ...
    flowCfg.identification.analysisStartTimeSec, flowCfg.identification.targetBladePasses, ...
    flowCfg.identification.windowBladePasses, flowCfg.identification.slidingStepBladePasses);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinateCheck.max_abs_delta_mm, coordinateCheck.status);
if ~coordinateCheck.is_consistent
    warning('Template/DynamicMap x-center mismatch. Do not use this run as the unified GradientXRange030 main result.');
end
if ~isempty(directMainFile)
    fprintf('Direct reference: %s\n', directMainFile);
end
if ~isempty(compactWaveformCaseFile)
    fprintf('Compact waveform case: %s\n', compactWaveformCaseFile);
end
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    fprintf('Foundation Step05 bundle source: %s\n', foundationStep05File);
end
fprintf('Windows: run %d/%d available, sensors: %s, VP Top-K: %d, EO candidate mode: %s\n', ...
    numWindows, availableWindows, mat2str(cfg.analysisSensors), cfg.vpTopK, cfg.eoCandidateMode);
if ~isempty(cfg.maxWindowsEnvIgnored)
    fprintf(['Ignored STEP07J_MAX_WINDOWS=%s; run-window limit is controlled by ' ...
        'ProjectionFlow_Config_20241106 identification.maxRunWindows.\n'], ...
        cfg.maxWindowsEnvIgnored);
end
if numWindows < availableWindows
    fprintf('Window run limit active: ProjectionFlow_Config_20241106 identification.maxRunWindows = %d.\n', ...
        cfg.maxWindows);
end
fprintf('Gap-aware sensors: %s; direct-only vibration sensors: %s.\n', ...
    mat2str(cfg.gapSensors), mat2str(cfg.directOnlySensors));
fprintf('Prev-window EO soft candidate: %s, refresh every %d window(s).\n', ...
    cfg.prevWindowCandidateMode, cfg.prevWindowRefreshEvery);
fprintf('Phase-safe expansion: %d, margin %.3f mm, fallback %s.\n', ...
    cfg.phaseSafeExpansion, cfg.phaseSafeMarginMm, cfg.phaseSafeFallbackMode);
fprintf('Direct bundle source mode: %s.\n', cfg.directBundleMode);
fprintf('Local foundation-core bundle source: %s (DynamicMap route enabled: %d).\n', ...
    cfg.localBundleSource, cfg.allowDynamicMapFallback);
if cfg.allowFoundationTimeMismatch
    fprintf('Diagnostic Foundation-time override is active: retain matching Foundation bundles.\n');
end
fprintf('Fixed baseline EO mode: %s (legacy direct fixed override: %d).\n', ...
    cfg.fixedEoMode, cfg.useLegacyDirectFixed);
if cfg.skipFoundationFixedResult
    fprintf('Fixed baseline mode: refit current bundle; Foundation Step05 Result is not copied.\n');
end
fprintf('Direct-only x-center alignment: %d.\n', cfg.alignDirectOnlyXCenter);
fprintf(['Pulse/domain: %s + %s, soft margin %.3f mm, query guard %s ' ...
    '(fixed %.3f, min %.3f, max %.3f, q%.1f + %.3f mm)\n'], ...
    cfg.pulseSelectionMode, cfg.domainSelectionMode, cfg.domainSoftMarginMm, ...
    cfg.queryGuardMode, cfg.queryGuardFixedMm, cfg.queryGuardMinMm, ...
    cfg.queryGuardMaxMm, cfg.queryGuardQuantile, cfg.queryGuardSafetyMm);
fprintf(['Models: %s | main %s, EO screen %s | dx limit %.3f mm, ' ...
    'dg %.3f mm, dmu %.3f, dtau %.3f mm\n'], ...
    strjoin(cfg.modelNames, ', '), cfg.mainModel, cfg.mainEoScreenModel, ...
    cfg.dxLimitMm, cfg.deltaGapLimitMm, cfg.deltaMuLimit, cfg.deltaTauLimitMm);
fprintf('VP seed: %s, gradient gate %.3f*q%.1f, weight power %.3f.\n', ...
    cfg.vpSeedObservationMode, cfg.vpGradientMinRatio, ...
    cfg.vpGradientReferenceQuantile, cfg.vpGradientWeightPower);
fprintf('Final forward-fit objective: %s voltage RMSE.\n', cfg.finalObjective);
if cfg.lockGapModelsToFixedEO
    fprintf('Gap-stage EO decision: locked to fixed-model EO.\n');
else
    fprintf('Gap-stage EO decision: all retained VP Top-K candidates enter final gap-aware fitting.\n');
end
fprintf('Per-candidate continuous frequency refinement: %d (half-width %.3f Hz).\n', ...
    cfg.continuousFrequencyRefine, cfg.frequencyRefineHalfWidthHz);
fprintf('Low-confidence EO margin threshold: %.3f%%.\n', ...
    cfg.lowConfidenceEoMarginPercent);
fprintf('Eta policy: theoretical zero, eta=%s mm.\n', ...
    mat2str(cfg.fixedSensorEtaMm, 6));

WindowResult = struct([]);
trendRows = cell(numWindows, 1);
bestIdx = 1;
bestRmse = inf;
prevWindowResult = [];

for iw = 1:numWindows
    Wmap = DynamicMap.Window(iw);
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        Wmap = foundation_wmap_for_window_local(FoundationStep05Source.Window(iw), Wmap);
    end
    % Keep the formal main route window-local. Expanded points may use this
    % window's linear VP fallback, but never a preceding-window phase model.
    prevPhaseInfo = build_phase_safe_reference_local([], [], iw, cfg, 'window_local_only');
    useDirectBundle = strcmpi(cfg.directBundleMode, 'auto') || ...
        (strcmpi(cfg.directBundleMode, 'off_when_phase_safe') && ~cfg.phaseSafeExpansion);
    directBundle = direct_window_bundle_local(directRef, iw);
    useCompactFallback = false;
    if useDirectBundle && ~isempty(directBundle)
        coreBundleFull = build_gap_bundle_from_direct_bundle_local( ...
            directBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        bundleSource = 'direct_bundle';
    elseif strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        foundationWindow = FoundationStep05Source.Window(iw);
        coreBundleFull = build_gap_bundle_from_foundation_bundle_local( ...
            foundationWindow.CoreBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        coreBundleFull.localSource = 'foundation_step05_core_bundle';
        bundleSource = 'foundation_step05_bundle';
    elseif strcmpi(cfg.localBundleSource, 'compact_case') && ~isempty(CompactWaveformCase)
        compactCore = compact_case_window_bundle_local(CompactWaveformCase, Wmap, cfg, []);
        coreBundleFull = build_gap_bundle_from_direct_bundle_local( ...
            compactCore, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        coreBundleFull.localSource = 'compact_case_core';
        useCompactFallback = true;
        bundleSource = 'compact_case';
    elseif strcmpi(cfg.localBundleSource, 'compact_case')
        error(['Compact waveform case is required for compact_case reconstruction but is missing. ' ...
            'Generate the matching CompactWaveformCase, use direct bundle auto mode with a valid direct result, ' ...
            'or use the current STEP07J_LOCAL_BUNDLE_SOURCE=dynamic_map foundation-core route.']);
    elseif strcmpi(cfg.localBundleSource, 'dynamic_map') && cfg.allowDynamicMapFallback
        coreBundleFull = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        coreBundleFull.localSource = 'dynamic_map_foundation_core';
        bundleSource = 'dynamic_map_foundation_core';
    else
        error('No valid bundle source is available for window %d.', iw);
    end
    coreBundle = decimate_bundle_local(coreBundleFull, cfg.maxPointsPerWindow);
    eoCandidates = build_eo_candidates_local(coreBundle.rotFreqMeanHz, cfg.freqSearchHz, cfg.eoPad);
    coreSeedTable = solve_vp_seed_scan_local(coreBundle, eoCandidates, cfg);
    phaseInfo = prevPhaseInfo;
    if cfg.phaseSafeExpansion && ~phaseInfo.usable && strcmpi(cfg.phaseSafeFallbackMode, 'linear_vp')
        phaseInfo = build_phase_safe_reference_local([], coreSeedTable, iw, cfg, prevPhaseInfo.reason);
    end

    bundleFull = coreBundleFull;
    usedPhaseSafeExpansion = false;
    expandedBundleFull = [];
    if cfg.phaseSafeExpansion && phaseInfo.usable
        if useCompactFallback
            compactExpanded = compact_case_window_bundle_local(CompactWaveformCase, Wmap, cfg, phaseInfo);
            expandedBundleFull = build_gap_bundle_from_direct_bundle_local( ...
                compactExpanded, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
            expandedBundleFull.localSource = 'compact_case_phase_safe';
        elseif strcmpi(bundleSource, 'direct_bundle') && ~isempty(CompactWaveformCase)
            compactExpanded = compact_case_window_bundle_local(CompactWaveformCase, Wmap, cfg, phaseInfo);
            expandedBundleFull = build_gap_bundle_from_direct_bundle_local( ...
                compactExpanded, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
            expandedBundleFull.localSource = 'compact_case_phase_safe_from_direct_core';
        elseif strcmpi(bundleSource, 'dynamic_map_foundation_core') && cfg.allowDynamicMapFallback
            expandedBundleFull = build_gap_observation_bundle_local( ...
                Wmap, Template, CorrectedGapLibrary, responseSurface, cfg, phaseInfo);
            expandedBundleFull.localSource = 'dynamic_map_foundation_core_phase_safe';
        elseif strcmpi(bundleSource, 'foundation_step05_bundle')
            foundationWindow = FoundationStep05Source.Window(iw);
            expandedBundleFull = build_gap_bundle_from_foundation_bundle_local( ...
                foundationWindow.FinalBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
            expandedBundleFull.localSource = 'foundation_step05_final_bundle';
        end
        if ~isempty(expandedBundleFull) && expandedBundleFull.pointCount > coreBundleFull.pointCount
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
    selectionInfo.forcedEO = NaN;
    selectionInfo.forceReason = '';
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
    modelFits = struct();
    fixedCandidateEO = model_eo_candidates_local('fixed', selectedEO, selectionInfo, cfgWindow);
    modelFits.fixed = [];
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle') && ...
            ~cfg.skipFoundationFixedResult
        modelFits.fixed = fixed_fit_from_foundation_step05_local( ...
            FoundationStep05Source.Window(iw), bundle, cfgWindow);
    end
    if isempty(modelFits.fixed)
        modelFits.fixed = refine_static_warp_fit_local(bundle, seedTable, fixedCandidateEO, cfgWindow, 'fixed');
    end
    fitModelNames = setdiff(cfg.modelNames, {'fixed'}, 'stable');
    for im = 1:numel(fitModelNames)
        modeName = fitModelNames{im};
        modelCandidateEO = model_eo_candidates_local(modeName, selectedEO, selectionInfo, cfgWindow);
        if cfg.lockGapModelsToFixedEO && mode_has_gap_local(modeName) && isfinite(modelFits.fixed.EO)
            modelCandidateEO = modelFits.fixed.EO;
        end
        modelFits.(modeName) = refine_static_warp_fit_local( ...
            bundle, seedTable, modelCandidateEO, cfgWindow, modeName, modelFits.fixed);
    end
    modelFits.(cfg.mainModel).qualityStatus = validate_formal_fit_quality_local( ...
        modelFits.(cfg.mainModel), cfgWindow, iw, cfg.mainModel);
    if cfg.useLegacyDirectFixed && strcmpi(cfg.runMode, 'comparison') && isfield(modelFits, 'fixed')
        directFixedFit = fixed_fit_from_direct_reference_local(directRef, iw, bundle, cfgWindow);
        if ~isempty(directFixedFit)
            modelFits.fixed = directFixedFit;
        end
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
    fprintf('%s\n', format_window_identification_progress_local( ...
        iw, numWindows, Wmap, bundleFull, modelFits.(cfg.mainModel), cfgWindow));
end

Trend = vertcat(trendRows{:});
Summary = build_summary_table_local(Trend, bestIdx, cfg);

Result = struct();
Result.dataset = '20241106';
if strcmpi(cfg.runMode, 'comparison')
    Result.method = 'gap_tilt_comparison_with_nested_static_warp_ablation_vp_full_wave';
else
    Result.method = 'gap_fixed_tilt_main_vp_full_wave';
end
Result.description = ['Current Step07J route enforces theoretical zero eta, gradient-displacement VP seeds, free EO selection, and full-waveform RMSE refinement/ranking. ' ...
    'Each three-lap sliding window is identified separately within 300-1000 Hz; no multi-window shared-frequency or phase-continuity layer is used. ' ...
    'The formal route requires an exact time/window/frequency-matched Foundation Step05 source and refuses DynamicMap fallback. VP scans all EO cheaply, then only its adaptive Top-K/neighbor set enters nonlinear full-wave fitting. ' ...
    'The fixed-gap model no longer locks the gap stage to one EO: every retained VP Top-K candidate enters the final gap-aware forward fit, and the configured final voltage objective selects EO. The package default is ordinary unweighted voltage RMSE; weighting is retained only for VP initialization. No EO value is hard-coded. Boundary-hitting solutions are rejected. ' ...
    'Every retained candidate is independently refined with f=EO*fr+delta_f inside the configured local bound before final EO selection; set STEP07J_REFINE_EACH_CANDIDATE_FREQUENCY=0 only for the integer-frequency ablation. ' ...
    'The paper-facing gap_only/gap_fixed_tilt model identifies only the per-sensor gap increment dg; the low-speed calibrated sensor slope mu remains fixed and no residual high-speed dmu is identified. gap_only with zero dg reduces to the fixed no-gap kernel. ' ...
    'The default main run computes only gap_only. fixed and gap_tilt are computed only in STEP07J_RUN_MODE=comparison. gap_tilt is a sensitivity ablation because dmu can compete with vibration amplitude; gap_tilt_shift remains excluded because the extra sensor-shift freedom is invalid for the paper-facing method.'];
Result.cfg = cfg;
Result.flowConfig = flowCfg;
Result.templateFile = templateFile;
Result.dynamicFile = dynamicFile;
Result.correctedLibFile = correctedLibFile;
Result.directMainFile = directMainFile;
Result.compactWaveformCaseFile = compactWaveformCaseFile;
Result.foundationStep05File = foundationStep05File;
Result.FoundationStep05TimeCheck = foundationTimeCheck;
Result.CoordinateCheck = coordinateCheck;
Result.WindowResult = WindowResult;
Result.Trend = Trend;
Result.Summary = Summary;
Result.BestWindowIndex = bestIdx;
Result.MainModel = cfg.mainModel;
Result.BestWindow = WindowResult(bestIdx);

matFile = fullfile(resultOutDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s%s.mat', ...
    flowCfg.dataset, caseTag, resultSuffix));
trendCsv = fullfile(resultOutDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s%s.csv', ...
    flowCfg.dataset, caseTag, resultSuffix));
summaryCsv = fullfile(resultOutDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_Summary_%s_%s%s.csv', ...
    flowCfg.dataset, caseTag, resultSuffix));
figTrend = fullfile(figDir, sprintf('Step07J_Trend_%s_%s%s.png', flowCfg.dataset, caseTag, resultSuffix));
figBest = fullfile(figDir, sprintf('Step07J_BestWindow_Collapse_%s_%s%s.png', flowCfg.dataset, caseTag, resultSuffix));
figGap = fullfile(figDir, sprintf('Step07J_StaticWarp_%s_%s%s.png', flowCfg.dataset, caseTag, resultSuffix));

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
function P = apply_case_overrides_to_flow_config_local(P)
P.identification.targetBlade = parse_positive_integer_override_env_local( ...
    'BLADE_CASE_BLADE_ID', P.identification.targetBlade);
P.identification.analysisSensors = parse_int_vector_env_local( ...
    'BLADE_CASE_SENSOR_IDS', P.identification.analysisSensors);
P.identification.analysisSensors = parse_int_vector_env_local( ...
    'STEP07J_ANALYSIS_SENSORS', P.identification.analysisSensors);
P.identification.gapSensors = parse_int_vector_env_local( ...
    'STEP06G_GAP_SENSORS', P.identification.gapSensors);
P.identification.gapSensors = parse_int_vector_env_local( ...
    'STEP07J_GAP_SENSORS', P.identification.gapSensors);
P.identification.analysisStartTimeSec = parse_positive_numeric_env_local( ...
    'STEP06G_START_TIME_SEC', P.identification.analysisStartTimeSec);
P.identification.targetBladePasses = parse_positive_integer_override_env_local( ...
    'STEP06G_TARGET_LAPS', P.identification.targetBladePasses);
P.identification.windowBladePasses = parse_positive_integer_override_env_local( ...
    'STEP06G_WINDOW_LAPS', P.identification.windowBladePasses);
P.identification.slidingStepBladePasses = parse_positive_integer_override_env_local( ...
    'STEP06G_SLIDING_STEP_LAPS', P.identification.slidingStepBladePasses);
P.identification.analysisSensors = reshape(P.identification.analysisSensors, 1, []);
P.identification.gapSensors = intersect(P.identification.gapSensors(:).', ...
    P.identification.analysisSensors, 'stable');
P.identification.directOnlySensors = setdiff(P.identification.analysisSensors, ...
    P.identification.gapSensors, 'stable');
end

function value = parse_positive_integer_override_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 1 || abs(value - round(value)) > eps(value)
    error('%s must be a positive integer.', name);
end
value = round(value);
end

function value = parse_positive_numeric_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value <= 0
    error('%s must be a positive number.', name);
end
end

function values = parse_int_vector_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValue;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values) || any(values < 1)
    error('%s must contain positive integer IDs separated by spaces.', name);
end
end

function templateFile = find_low_speed_template_file_local(templateDirs, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
if ischar(templateDirs) || isstring(templateDirs)
    templateDirs = cellstr(templateDirs);
end
patterns = {
    sprintf('TemplateBundle_LowSpeedRotating_B%d_S2357_GradientXRange030_OPRCenterStd_20241106.mat', targetBlade)
    sprintf('TemplateBundle_LowSpeedRotating_B%d_S2357_OPRCenterStdRef_20241106.mat', targetBlade)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_OPRCenterStdRef_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_OPRCenterStd_20241106.mat', ['B', num2str(targetBlade), '_', sensorTag])
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_BaseFrame_20241106.mat', targetBlade, sensorTag)
    sprintf('Step06I_Generated_LowSpeed_Template_20241106_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06I_Generated_LowSpeed_Template_20241106_%s.mat', ['B', num2str(targetBlade), '_', sensorTag])
    sprintf('Step06I_Generated_LowSpeed_Template_20241106_B%d_%s*.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    for id = 1:numel(templateDirs)
        candidate = fullfile(templateDirs{id}, patterns{i});
        if isfile(candidate)
            templateFile = candidate;
            return;
        end
        files = dir(fullfile(templateDirs{id}, '**', patterns{i}));
        if ~isempty(files)
            [~, idx] = max([files.datenum]);
            templateFile = fullfile(files(idx).folder, files(idx).name);
            return;
        end
    end
end
error('No low-speed template file found under %s for %s.', strjoin(templateDirs, '; '), sensorTag);
end

function Template = load_low_speed_template_local(templateFile, sensorIds, sensorTag)
S = load(templateFile);
if isfield(S, 'Template')
    Template = S.Template;
elseif isfield(S, 'sensor_template')
    Template = S.sensor_template;
elseif isfield(S, 'TemplateBundle')
    Template = load_template_bundle_index_local(S.TemplateBundle, sensorIds);
else
    error('Unsupported template file %s.', templateFile);
end
if isfield(Template, 'Sensor')
    Template = filter_template_sensors_local(Template, sensorIds, sensorTag);
else
    error('Template file %s does not provide Sensor entries.', templateFile);
end
end

function [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
    rootDir, targetBlade, sensorIds, sensorTag)
packageCfg = Config_20241106();
templateDir = fullfile(packageCfg.paths.calibrationFoundation, ...
    'step04_low_speed_template');
preferred = fullfile(templateDir, ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat');
if isfile(preferred)
    templateFile = preferred;
else
    files = dir(fullfile(templateDir, ...
        'Template_OPRCenterStd_LowSpeed_AllBlades_*_20241106.mat'));
    if isempty(files)
        error('No foundation Step04 template found under %s.', templateDir);
    end
    [~, idx] = max([files.datenum]);
    templateFile = fullfile(files(idx).folder, files(idx).name);
end
S = load(templateFile, 'Template');
if ~isfield(S, 'Template') || ~isfield(S.Template, 'SensorBlade')
    error('Foundation Step04 template must contain Template.SensorBlade: %s', templateFile);
end
Sensor = repmat(struct('sensor_id', NaN, 'x_grid', [], 'v_grid', [], ...
    'dv_dx', [], 'bin_weight', [], 'x_domain', [NaN NaN], ...
    'threshold', NaN, 'xc', 0, 'baseline', 0), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    hit = [];
    for i = 1:numel(S.Template.SensorBlade)
        if S.Template.SensorBlade(i).sensor_id == sid && ...
                S.Template.SensorBlade(i).blade_id == targetBlade
            hit = S.Template.SensorBlade(i);
            break;
        end
    end
    if isempty(hit)
        error('Foundation Step04 template %s does not contain B%d CH%d.', ...
            templateFile, targetBlade, sid);
    end
    Sensor(is).sensor_id = sid;
    Sensor(is).x_grid = hit.x_grid(:);
    Sensor(is).v_grid = hit.v_grid(:);
    Sensor(is).dv_dx = hit.dv_dx(:);
    if isfield(hit, 'weight_grid') && ~isempty(hit.weight_grid)
        Sensor(is).bin_weight = hit.weight_grid(:);
    else
        Sensor(is).bin_weight = ones(size(Sensor(is).x_grid));
    end
    Sensor(is).x_domain = hit.x_domain(:).';
    Sensor(is).threshold = hit.threshold;
    Sensor(is).xc = hit.xc;
    if isfield(hit, 'baseline') && isfinite(hit.baseline)
        Sensor(is).baseline = hit.baseline;
    end
end
Template = struct();
Template.Dataset = '20241106';
Template.SourceFile = templateFile;
Template.SourceMode = 'foundation_step04_sensorblade';
Template.Sensor = Sensor;
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = sensorTag;
Template.TargetBlade = targetBlade;
end

function Template = load_template_bundle_index_local(bundleIndex, sensorIds)
if ~isfield(bundleIndex, 'SensorFiles') || ~istable(bundleIndex.SensorFiles)
    error('TemplateBundle index is missing SensorFiles.');
end
Template = rmfield_if_present_local(bundleIndex, {'SensorFiles', 'BundleFile'});
Template.Sensor = struct([]);
for sid = sensorIds(:).'
    idx = find(bundleIndex.SensorFiles.sensor_id == sid, 1);
    if isempty(idx)
        error('TemplateBundle index does not contain CH%d.', sid);
    end
    sensorFile = char(bundleIndex.SensorFiles.template_file(idx));
    S = load(sensorFile);
    if isfield(S, 'sensor_template')
        sensorTemplate = S.sensor_template;
    elseif isfield(S, 'Template')
        sensorTemplate = S.Template;
    else
        error('Unsupported sensor template file %s.', sensorFile);
    end
    if ~isfield(sensorTemplate, 'Sensor') || isempty(sensorTemplate.Sensor)
        error('Sensor template file %s does not contain Sensor.', sensorFile);
    end
    sensor = sensorTemplate.Sensor;
    if numel(sensor) > 1
        sensor = sensor([sensor.sensor_id] == sid);
    end
    if isempty(sensor)
        error('Sensor template file %s does not contain CH%d.', sensorFile, sid);
    end
    if isempty(Template.Sensor)
        Template.Sensor = sensor;
    else
        Template.Sensor(end + 1) = sensor; %#ok<AGROW>
    end
end
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = ['S', sprintf('%d', sensorIds)];
end

function Template = merge_gap_template_sensors_local(Template, gapTemplate, gapSensors, sensorTag)
for sid = gapSensors(:).'
    idxDst = find([Template.Sensor.sensor_id] == sid, 1);
    idxSrc = find([gapTemplate.Sensor.sensor_id] == sid, 1);
    if isempty(idxDst) || isempty(idxSrc)
        error('Cannot merge gap template for CH%d.', sid);
    end
    [Template.Sensor, gapSensor] = align_struct_array_fields_local(Template.Sensor, gapTemplate.Sensor(idxSrc));
    Template.Sensor(idxDst) = gapSensor;
end
Template.SensorTag = sensorTag;
end

function Template = load_template_from_bank_local(templateBankFile, targetBlade, sensorIds, sensorTag)
S = load(templateBankFile, 'LowSpeedTemplateBank');
bank = S.LowSpeedTemplateBank;
Template = struct();
Template.dataset = bank.dataset;
Template.targetBlade = targetBlade;
Template.analysisSensors = sensorIds(:).';
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = sensorTag;
Template.Sensor = struct([]);
for sid = sensorIds(:).'
    idx = find([bank.entry.bladeId] == targetBlade & [bank.entry.sensorId] == sid, 1, 'first');
    if isempty(idx)
        error('LowSpeedTemplateBank does not contain B%d CH%d.', targetBlade, sid);
    end
    Template.Sensor = append_sensor_struct_local(Template.Sensor, bank.entry(idx).sensorTemplate);
end
end

function CorrectedGapLibrary = load_corrected_gap_library_for_run_local(correctedLibFile, cfg, Template)
S = load(correctedLibFile);
if isfield(S, 'CorrectedGapLibrary')
    CorrectedGapLibrary = S.CorrectedGapLibrary;
    return;
end
if ~isfield(S, 'GapCalibrationBank')
    error('Unsupported corrected gap source: %s', correctedLibFile);
end
bank = S.GapCalibrationBank;
sensorCorr = struct([]);
for is = 1:numel(cfg.gapSensors)
    sid = cfg.gapSensors(is);
    idx = find([bank.entry.bladeId] == cfg.targetBlade & [bank.entry.sensorId] == sid, 1, 'first');
    if isempty(idx)
        error('GapCalibrationBank does not contain B%d CH%d.', cfg.targetBlade, sid);
    end
    sensorCorr = append_sensor_struct_local(sensorCorr, bank.entry(idx).sensor);
end
if isempty(Template)
    Template = load_template_from_bank_local(bank.templateBankFile, cfg.targetBlade, cfg.gapSensors, ...
        ['S', sprintf('%d', cfg.gapSensors)]);
else
    Template = filter_template_sensors_local(Template, cfg.gapSensors, ['S', sprintf('%d', cfg.gapSensors)]);
end
CorrectedGapLibrary = struct();
CorrectedGapLibrary.dataset = bank.dataset;
CorrectedGapLibrary.method = 'gap_bank_view_for_step07j';
CorrectedGapLibrary.description = 'Per-run view assembled from all-blade GapCalibrationBank.';
CorrectedGapLibrary.responseFile = bank.responseFile;
CorrectedGapLibrary.templateBankFile = bank.templateBankFile;
CorrectedGapLibrary.templateFile = string(bank.templateBankFile);
CorrectedGapLibrary.templateSourceMode = 'low_speed_template_bank';
CorrectedGapLibrary.templateSourceFile = string(bank.templateBankFile);
CorrectedGapLibrary.targetBlade = cfg.targetBlade;
CorrectedGapLibrary.analysisSensors = cfg.gapSensors(:).';
CorrectedGapLibrary.gapSensors = cfg.gapSensors(:).';
CorrectedGapLibrary.directOnlySensors = cfg.directOnlySensors(:).';
CorrectedGapLibrary.responseSurface = bank.responseSurface;
CorrectedGapLibrary.lowSpeedTemplate = Template;
CorrectedGapLibrary.sensor = sensorCorr;
CorrectedGapLibrary.formula = 'T_low(x)+a*(F_raw(g+mu*(x-tau),k*(x-tau))-F_raw(g0+mu*(x-tau),k*(x-tau)))';
CorrectedGapLibrary.cfg = bank.cfg;
end

function sensorArray = append_sensor_struct_local(sensorArray, sensor)
if isempty(sensorArray)
    sensorArray = sensor;
    return;
end
[sensorArray, sensor] = align_struct_array_fields_local(sensorArray, sensor);
sensorArray(end + 1) = sensor; %#ok<AGROW>
end

function [arrayOut, itemOut] = align_struct_array_fields_local(arrayIn, itemIn)
arrayOut = arrayIn;
itemOut = itemIn;
fieldsAll = union(fieldnames(arrayOut), fieldnames(itemOut), 'stable');
for i = 1:numel(fieldsAll)
    name = fieldsAll{i};
    if ~isfield(arrayOut, name)
        [arrayOut.(name)] = deal([]);
    end
    if ~isfield(itemOut, name)
        itemOut.(name) = [];
    end
end
arrayOut = orderfields(arrayOut, fieldsAll);
itemOut = orderfields(itemOut, fieldsAll);
end

function dynamicFile = find_dynamic_map_file_local(dynamicDirs, targetBlade, sensorIds, ...
        analysisStartTimeSec, windowLaps, slidingStepLaps)
sensorTag = ['S', sprintf('%d', sensorIds)];
if ischar(dynamicDirs) || isstring(dynamicDirs)
    dynamicDirs = cellstr(dynamicDirs);
end
patterns = {
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s_W%dS%d.mat', ...
        targetBlade, sensorTag, strrep(sprintf('T%07.3f', analysisStartTimeSec), '.', 'p'), ...
        windowLaps, slidingStepLaps)
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s_W%dS%d*.mat', ...
        targetBlade, sensorTag, strrep(sprintf('T%07.3f', analysisStartTimeSec), '.', 'p'), ...
        windowLaps, slidingStepLaps)
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s.mat', ...
        targetBlade, sensorTag, strrep(sprintf('T%07.3f', analysisStartTimeSec), '.', 'p'))
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_%s.mat', ['B', num2str(targetBlade), '_', sensorTag])
    sprintf('Step06_BuildGapAwareDynamicMap_20241106_B%d_%s*.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_T075s_GradientXRange030_OPRCenterStd_20241106.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20241106.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_%s_SlidingWindows*_20241106.mat', ['B', num2str(targetBlade), '_', sensorTag])
    sprintf('DynamicMap_B%d_%s_SlidingWindows*OPRCenterStd*_20241106.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20241106.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20241106.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    for id = 1:numel(dynamicDirs)
        files = dir(fullfile(dynamicDirs{id}, patterns{ip}));
        if ~isempty(files)
            [~, idx] = max([files.datenum]);
            dynamicFile = fullfile(files(idx).folder, files(idx).name);
            return;
        end
    end
end
error('No DynamicMap file found under %s for %s.', strjoin(dynamicDirs, '; '), sensorTag);
end

function correctedLibFile = find_corrected_gap_library_file_local(outDir, targetBlade, sensorTag, gapBankFile)
if nargin >= 4 && ~isempty(gapBankFile) && isfile(gapBankFile)
    correctedLibFile = gapBankFile;
    return;
end
patterns = {
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_%s.mat', ['B', num2str(targetBlade), '_', sensorTag])
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s*.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        correctedLibFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
correctedLibFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
    targetBlade, sensorTag));
end

function resultFile = find_direct_main_result_file_local(resultDir, targetBlade, sensorIds, startTimeSec)
sensorTag = ['S', sprintf('%d', sensorIds)];
timeTag = time_label_step07j_local(startTimeSec);
legacyResultDir = strrep(resultDir, fullfile('output', 'new_flow', '06_identification'), fullfile('output', 'identification'));
searchDirs = {
    fullfile(resultDir, sensorTag)
    resultDir
    legacyResultDir
    };
preferred = {
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta200.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta000.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_*.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_only.mat', timeTag, targetBlade)
    };
for id = 1:numel(searchDirs)
    if exist(searchDirs{id}, 'dir') ~= 7
        continue;
    end
    for ip = 1:numel(preferred)
        files = dir(fullfile(searchDirs{id}, preferred{ip}));
        if ~isempty(files)
            [~, idx] = max([files.datenum]);
            resultFile = fullfile(files(idx).folder, files(idx).name);
            return;
        end
    end
end
resultFile = '';
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
S = load(resultFile);
if isfield(S, 'Result') && isfield(S.Result, 'Trend')
    directRef.Trend = S.Result.Trend;
end
if isfield(S, 'Result') && isfield(S.Result, 'WindowResult')
    directRef.WindowResult = S.Result.WindowResult;
end
if isfield(S, 'IdentificationResult') && isfield(S.IdentificationResult, 'Trend')
    directRef.Trend = S.IdentificationResult.Trend;
end
if isfield(S, 'IdentificationResult') && isfield(S.IdentificationResult, 'WindowResult')
    directRef.WindowResult = S.IdentificationResult.WindowResult;
end
end

function resultFile = find_foundation_step05_result_file_local(rootDir, targetBlade, sensorIds, ...
        analysisStartTimeSec, windowLaps, slidingStepLaps, freqSearchHz)
resultFile = '';
sensorTag = ['S', sprintf('%d', sensorIds)];
packageCfg = Config_20241106();
foundationOutputRoot = packageCfg.paths.foundationResults;
timeTag = strrep(sprintf('T%07.3f', analysisStartTimeSec), '.', 'p');
runTag = sprintf('step05_single_sync_direct_template_%s_W%dS%d_F%d_%d', ...
    timeTag, windowLaps, slidingStepLaps, ...
    round(freqSearchHz(1)), round(freqSearchHz(2)));
patterns = {
    fullfile(foundationOutputRoot, runTag, '*', sprintf( ...
        'Result_Step05_FoundationMainPulseAdaptiveZeroEta_B%d_%s_20241106.mat', ...
        targetBlade, sensorTag))
    };
for ip = 1:numel(patterns)
    files = dir(patterns{ip});
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        resultFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
end

function Source = load_foundation_step05_source_local(resultFile, cfg)
Source = [];
if isempty(resultFile) || exist(resultFile, 'file') ~= 2
    return;
end
S = load(resultFile, 'Result');
if ~isfield(S, 'Result') || ~isfield(S.Result, 'WindowResult')
    warning('Unsupported foundation Step05 result: %s', resultFile);
    return;
end
if ~isfield(S.Result, 'AnalysisSettings')
    error('Foundation Step05 result lacks AnalysisSettings: %s', resultFile);
end
A = S.Result.AnalysisSettings;
requiredSettings = {'analysis_start_time_s','window_laps','sliding_step_laps', ...
    'freq_search_hz','eo_candidate_policy'};
if ~all(isfield(A, requiredSettings)) || ...
        abs(A.analysis_start_time_s - cfg.analysisStartTimeSec) > 1e-6 || ...
        A.window_laps ~= cfg.windowBladePasses || ...
        A.sliding_step_laps ~= cfg.slidingStepBladePasses || ...
        ~isequal(double(A.freq_search_hz(:).'), double(cfg.freqSearchHz(:).')) || ...
        ~strcmpi(A.eo_candidate_policy, 'adaptive_full_fit')
    error(['Foundation Step05 settings do not exactly match the active case. ' ...
        'Expected T=%.6f, W%dS%d, F%s. File: %s'], ...
        cfg.analysisStartTimeSec, cfg.windowBladePasses, cfg.slidingStepBladePasses, ...
        mat2str(cfg.freqSearchHz), resultFile);
end
if ~isfield(S.Result, 'TargetBlade') || S.Result.TargetBlade ~= cfg.targetBlade || ...
        ~isfield(S.Result, 'SensorIDs') || ...
        ~isequal(double(S.Result.SensorIDs(:).'), double(cfg.analysisSensors(:).'))
    error('Foundation Step05 blade/sensor selection does not match the active case: %s', resultFile);
end
Wsrc = S.Result.WindowResult;
Window = repmat(struct('FinalBundle', [], 'CoreBundle', [], 'FoundationResult', [], ...
    'lap_range', [NaN NaN], 'time_window_s', [NaN NaN]), numel(Wsrc), 1);
for iw = 1:numel(Wsrc)
    wr = Wsrc(iw);
    finalBundle = foundation_pick_bundle_local(wr, ...
        {'BundlePreview', 'ExpandedBundlePreview', 'CoreBundlePreview', 'CandidateBundlePreview'});
    coreBundle = foundation_pick_bundle_local(wr, ...
        {'CoreBundlePreview', 'BundlePreview', 'ExpandedBundlePreview', 'CandidateBundlePreview'});
    if isempty(finalBundle)
        continue;
    end
    if isfield(wr, 'Result')
        Window(iw).FoundationResult = wr.Result;
    end
    Window(iw).FinalBundle = foundation_filter_bundle_sensors_local(finalBundle, cfg.analysisSensors);
    if isempty(coreBundle)
        coreBundle = finalBundle;
    end
    Window(iw).CoreBundle = foundation_filter_bundle_sensors_local(coreBundle, cfg.analysisSensors);
    if isfield(wr, 'lap_range')
        Window(iw).lap_range = wr.lap_range;
    end
    if isfield(wr, 'time_window_s')
        Window(iw).time_window_s = wr.time_window_s;
    end
end
Source = struct();
Source.file = resultFile;
Source.Window = Window;
if isfield(S.Result, 'SensorIDs')
    Source.sensorIds = S.Result.SensorIDs(:).';
end
if isfield(S.Result, 'Template')
    Source.Template = S.Result.Template;
end
if isfield(S.Result, 'Trend')
    Source.Trend = S.Result.Trend;
else
    Source.Trend = table();
end
end

function [isMatch, info] = foundation_time_matches_dynamic_map_local(Source, DynamicMap)
info = struct('foundation_span_s', [NaN NaN], ...
    'dynamic_span_s', [NaN NaN], 'tolerance_s', 0.1, 'action', 'none');
isMatch = true;
if isempty(Source) || ~isfield(Source, 'Window') || isempty(Source.Window) || ...
        ~isfield(DynamicMap, 'Window') || isempty(DynamicMap.Window)
    return;
end

foundation_windows = reshape([Source.Window.time_window_s], 2, []).';
dynamic_windows = reshape([DynamicMap.Window.time_window], 2, []).';
info.foundation_span_s = [min(foundation_windows(:, 1)), max(foundation_windows(:, 2))];
info.dynamic_span_s = [min(dynamic_windows(:, 1)), max(dynamic_windows(:, 2))];
isMatch = all(isfinite([info.foundation_span_s, info.dynamic_span_s])) && ...
    max(abs(info.foundation_span_s - info.dynamic_span_s)) <= info.tolerance_s;
if isMatch
    info.action = 'foundation_step05_bundle';
end
end

function B = foundation_pick_bundle_local(wr, names)
B = [];
for i = 1:numel(names)
    name = names{i};
    if isfield(wr, name) && isstruct(wr.(name)) && ...
            isfield(wr.(name), 'x') && isfield(wr.(name), 'v')
        B = wr.(name);
        return;
    end
end
end

function B = foundation_filter_bundle_sensors_local(B, sensorIds)
if isempty(B) || ~isfield(B, 'sensor_id')
    return;
end
keep = ismember(B.sensor_id(:), sensorIds(:));
fields = fieldnames(B);
for i = 1:numel(fields)
    f = fields{i};
    v = B.(f);
    if isnumeric(v) || islogical(v)
        if numel(v) == numel(keep)
            B.(f) = v(keep);
        end
    end
end
B.point_count = nnz(keep);
if isfield(B, 'core_mask') && numel(B.core_mask) == B.point_count
    B.core_point_count = nnz(B.core_mask);
else
    B.core_point_count = B.point_count;
end
end

function [etaMm, info] = fixed_eta_from_foundation_step05_source_local(Source, sensorIds)
etaMm = [];
info = struct('file', Source.file, 'source', 'foundation_step05_result_sensor_eta_id');
if ~isstruct(Source) || ~isfield(Source, 'Window')
    return;
end
for iw = 1:numel(Source.Window)
    if ~isfield(Source.Window(iw), 'FoundationResult') || ...
            ~isstruct(Source.Window(iw).FoundationResult)
        continue;
    end
    R = Source.Window(iw).FoundationResult;
    if isfield(R, 'sensor_eta_id') && numel(R.sensor_eta_id) >= numel(sensorIds)
        eta = R.sensor_eta_id(:).';
        sourceSensorIds = [];
        if isfield(R, 'sensor_ids')
            sourceSensorIds = R.sensor_ids(:).';
        elseif isfield(Source, 'sensorIds')
            sourceSensorIds = Source.sensorIds(:).';
        end
        if numel(sourceSensorIds) == numel(eta) && all(ismember(sensorIds, sourceSensorIds))
            [~, loc] = ismember(sensorIds, sourceSensorIds);
            eta = eta(loc);
        else
            eta = eta(1:numel(sensorIds));
        end
        if all(isfinite(eta))
            etaMm = eta;
            info.window_id = iw;
            info.sensorIds = sensorIds(:).';
            return;
        end
    end
end
end

function Wmap = foundation_wmap_for_window_local(Fw, Wfallback)
Wmap = Wfallback;
if isfield(Fw, 'FinalBundle') && isstruct(Fw.FinalBundle)
    B = Fw.FinalBundle;
    if isfield(B, 'rot_freq_mean_hz') && isfinite(B.rot_freq_mean_hz)
        Wmap.rot_freq_mean_hz = B.rot_freq_mean_hz;
    end
    if isfield(B, 'rot_rpm_mean') && isfinite(B.rot_rpm_mean)
        Wmap.rot_rpm_mean = B.rot_rpm_mean;
    elseif isfield(Wmap, 'rot_freq_mean_hz')
        Wmap.rot_rpm_mean = 60 .* Wmap.rot_freq_mean_hz;
    end
end
if isfield(Fw, 'lap_range') && numel(Fw.lap_range) == 2 && all(isfinite(Fw.lap_range))
    Wmap.lap_range = Fw.lap_range;
end
if isfield(Fw, 'time_window_s') && numel(Fw.time_window_s) == 2
    Wmap.time_window_s = Fw.time_window_s;
    Wmap.time_window = Fw.time_window_s;
end
end

function bundle = direct_window_bundle_local(directRef, windowId)
bundle = [];
if ~isstruct(directRef) || ~isfield(directRef, 'WindowResult') || ...
        numel(directRef.WindowResult) < windowId
    return;
end
wr = directRef.WindowResult(windowId);
candidateNames = {'bundle', 'Bundle'};
for ic = 1:numel(candidateNames)
    name = candidateNames{ic};
    if isfield(wr, name) && ~isempty(wr.(name))
        bundle = wr.(name);
        return;
    end
end
end

function caseFile = find_compact_waveform_case_file_local(rotDir, targetBlade, sensorIds, startTimeSec)
caseFile = '';
sensorTag = ['S', sprintf('%d', sensorIds)];
timeTag = time_label_step07j_local(startTimeSec);
searchDir = fullfile(rotDir, 'output', 'new_flow', '05_waveform_library', sensorTag);
patterns = {
    sprintf('Step06_CompactWaveformCase_B%d_%s_20241106.mat', targetBlade, timeTag)
    };
for ip = 1:numel(patterns)
    files = dir(fullfile(searchDir, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        caseFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
end

function label = time_label_step07j_local(tSec)
if abs(tSec - round(tSec)) < 1e-9
    label = sprintf('T%03ds', round(tSec));
else
    label = sprintf('T%07.3fs', tSec);
    label = strrep(label, '.', 'p');
end
end

function CompactWaveformCase = load_compact_waveform_case_local(caseFile, cfg)
CompactWaveformCase = [];
if isempty(caseFile) || exist(caseFile, 'file') ~= 2
    if strcmpi(cfg.localBundleSource, 'compact_case')
        warning(['Compact waveform case is not available; Step07J will require a matching compact case ', ...
            'unless a valid direct bundle is used.']);
    end
    return;
end
S = load(caseFile, 'CompactWaveformCase');
if ~isfield(S, 'CompactWaveformCase') || ~isfield(S.CompactWaveformCase, 'WaveformLibrary')
    warning('Unsupported compact waveform case: %s', caseFile);
    return;
end
CompactWaveformCase = S.CompactWaveformCase;
end

function Bsrc = compact_case_window_bundle_local(CompactWaveformCase, Wmap, cfg, phaseRef)
if nargin < 4
    phaseRef = [];
end
WaveformLibrary = CompactWaveformCase.WaveformLibrary;
bladeIdx = find([WaveformLibrary.Blade.blade_id] == cfg.targetBlade, 1, 'first');
if isempty(bladeIdx)
    error('Compact waveform case does not contain B%d.', cfg.targetBlade);
end
B = WaveformLibrary.Blade(bladeIdx);
lapRange = Wmap.lap_range;
X = []; T = []; V = []; W = []; Theta = []; F0 = []; Fx = []; sensorIndex = [];
sensorEtaLimit = nan(1, numel(cfg.analysisSensors));
sourceMeta = repmat(struct('sensor_id', NaN, 'raw_points', 0, 'pulse_points', 0, ...
    'dynamic_effective_points', 0, 'base_points', 0, 'query_guard_mm', NaN, ...
    'query_safe_points', 0, 'selected_points', 0), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    sidx = find([B.Sensor.sensor_id] == sid, 1, 'first');
    if isempty(sidx)
        continue;
    end
    Sdata = B.Sensor(sidx);
    Tpl = compact_sensor_template_local(Sdata);
    [tRaw, xRaw, vRaw, thetaRaw, f0Raw, fxRaw] = concatenate_compact_laps_local(Sdata.Lap, lapRange, Tpl);
    [mask, queryGuardMm, meta] = select_compact_step06_points_local( ...
        tRaw, xRaw, vRaw, thetaRaw, f0Raw, fxRaw, Tpl, cfg, sid, is, phaseRef); %#ok<ASGLU>
    sourceMeta(is) = meta;
    if nnz(mask) < 8
        continue;
    end
    wSel = build_step06_filtered_weight_local(tRaw(mask), xRaw(mask), vRaw(mask), Tpl, queryGuardMm, cfg);
    X = [X; xRaw(mask)]; %#ok<AGROW>
    T = [T; tRaw(mask)]; %#ok<AGROW>
    V = [V; vRaw(mask)]; %#ok<AGROW>
    W = [W; wSel(:)]; %#ok<AGROW>
    Theta = [Theta; thetaRaw(mask)]; %#ok<AGROW>
    F0 = [F0; f0Raw(mask)]; %#ok<AGROW>
    Fx = [Fx; fxRaw(mask)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, nnz(mask), 1)]; %#ok<AGROW>
    if isfield(Tpl, 'eta_limit_mm') && isfinite(Tpl.eta_limit_mm) && Tpl.eta_limit_mm > 0
        sensorEtaLimit(is) = Tpl.eta_limit_mm;
    end
end
if isempty(T)
    error('Compact waveform case produced no valid points for window %d.', Wmap.window_id);
end
Bsrc = struct();
Bsrc.X = X;
Bsrc.T = T;
Bsrc.T_rel = T - min(T);
Bsrc.V = V;
Bsrc.W = W;
Bsrc.Theta = Theta;
Bsrc.F0 = F0;
Bsrc.Fx = Fx;
Bsrc.sensor_index = sensorIndex;
Bsrc.sensor_ids = cfg.analysisSensors(:).';
Bsrc.point_count = numel(T);
Bsrc.time_window = [min(T), max(T)];
Bsrc.selection_pass = 'step06_compact_case_local_filtering';
Bsrc.sensor_eta_limit_mm = sensorEtaLimit;
Bsrc.forceDirectBaseCurve = true;
Bsrc.window_meta = sourceMeta;
Bsrc.query_guard_mm = NaN;
Bsrc.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
Bsrc.rot_rpm_mean = Wmap.rot_rpm_mean;
end

function Tpl = compact_sensor_template_local(Sdata)
Tpl = struct();
Tpl.sensor_id = Sdata.sensor_id;
Tpl.x_grid = Sdata.template_x(:);
Tpl.v_grid = Sdata.template_v(:);
Tpl.dv_dx = Sdata.template_dv_dx(:);
Tpl.xc = Sdata.template_xc;
Tpl.x_domain = Sdata.template_domain(:).';
Tpl.baseline = 0;
Tpl.threshold = NaN;
if isfield(Sdata, 'template_eta_limit_mm')
    Tpl.eta_limit_mm = Sdata.template_eta_limit_mm;
else
    Tpl.eta_limit_mm = NaN;
end
end

function [tRaw, xRaw, vRaw, thetaRaw, f0Raw, fxRaw] = concatenate_compact_laps_local(Lap, lapRange, Tpl)
tRaw = []; xRaw = []; vRaw = []; thetaRaw = []; f0Raw = []; fxRaw = [];
for k = lapRange(:).'
    if k < 1 || k > numel(Lap)
        continue;
    end
    L = Lap(k);
    t = L.t(:);
    x = L.x_rel(:);
    v = L.V(:);
    theta = L.theta(:);
    f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x, 'pchip', NaN);
    fx = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x, 'pchip', NaN);
    keep = isfinite(t) & isfinite(x) & isfinite(v) & isfinite(theta);
    tRaw = [tRaw; t(keep)]; %#ok<AGROW>
    xRaw = [xRaw; x(keep)]; %#ok<AGROW>
    vRaw = [vRaw; v(keep)]; %#ok<AGROW>
    thetaRaw = [thetaRaw; theta(keep)]; %#ok<AGROW>
    f0Raw = [f0Raw; f0(keep)]; %#ok<AGROW>
    fxRaw = [fxRaw; fx(keep)]; %#ok<AGROW>
end
end

function [mask, queryGuardMm, meta] = select_compact_step06_points_local( ...
        tRaw, xRaw, vRaw, thetaRaw, f0Raw, fxRaw, Tpl, cfg, sensorId, sensorLocalIndex, phaseRef)
meta = struct('sensor_id', sensorId, 'raw_points', numel(tRaw), 'pulse_points', 0, ...
    'dynamic_effective_points', 0, 'base_points', 0, 'query_guard_mm', NaN, ...
    'query_safe_points', 0, 'selected_points', 0);
threshold = cfg.defaultSensorThreshold;
if isfield(Tpl, 'threshold') && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
end
maskPulse = build_pulse_mask_local(vRaw, threshold, cfg.pulseSelectionMode);
maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, cfg);
xDomainUse = expand_sensor_domain_step06_local(Tpl.x_domain, cfg, sensorId);
maskDomain = xRaw >= xDomainUse(1) + cfg.domainMarginMm & ...
    xRaw <= xDomainUse(2) - cfg.domainMarginMm;
maskFinite = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRaw) & isfinite(thetaRaw);
maskBase = maskPulse & maskEffective & maskDomain & maskFinite;
queryGuardMm = resolve_query_guard_mm_sensor_scaled_local(Tpl, xRaw, vRaw, maskBase, cfg, sensorId);
xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
maskQuerySafe = xRaw >= xDomain(1) + queryGuardMm & xRaw <= xDomain(2) - queryGuardMm;
if isstruct(phaseRef) && isfield(phaseRef, 'usable') && phaseRef.usable
    maskPhaseSafe = build_phase_safe_query_mask_local(xRaw, thetaRaw, xDomain, sensorLocalIndex, phaseRef, cfg);
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
meta.pulse_points = nnz(maskPulse);
meta.dynamic_effective_points = nnz(maskEffective);
meta.base_points = nnz(maskBase);
meta.query_guard_mm = queryGuardMm;
meta.query_safe_points = nnz(maskQuerySafe & maskBase);
meta.selected_points = nnz(mask);
end

function queryGuardMm = resolve_query_guard_mm_sensor_scaled_local(Tpl, x, v, baseMask, cfg, sensorId)
queryGuardMm = cfg.queryGuardFixedMm;
if strcmpi(cfg.queryGuardMode, 'adaptive') && nnz(baseMask) >= 8
    xSel = x(baseMask);
    vSel = v(baseMask);
    xStat = invert_template_voltage_step06_local(Tpl, vSel, xSel);
    uApp = abs(xSel(:) - xStat(:));
    uApp = uApp(isfinite(uApp));
    if ~isempty(uApp)
        adaptiveGuard = prctile(uApp, min(max(cfg.queryGuardQuantile, 0), 100)) + ...
            cfg.queryGuardSafetyMm;
        queryGuardMm = min(max(adaptiveGuard, cfg.queryGuardMinMm), cfg.queryGuardMaxMm);
    end
end
if isempty(cfg.sensorQueryGuardScaleTable) || size(cfg.sensorQueryGuardScaleTable, 2) < 2
    return;
end
idx = find(cfg.sensorQueryGuardScaleTable(:, 1) == sensorId, 1, 'first');
if isempty(idx)
    return;
end
scale = cfg.sensorQueryGuardScaleTable(idx, 2);
if isfinite(scale) && scale > 0
    queryGuardMm = queryGuardMm * scale;
end
end

function xStat = invert_template_voltage_step06_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = nan(size(v));
for i = 1:numel(v)
    vv = v(i);
    if ~isfinite(vv)
        continue;
    end
    diffV = vGrid - vv;
    crossingX = [];
    for k = 1:(numel(vGrid) - 1)
        if diffV(k) == 0
            crossingX(end + 1, 1) = xGrid(k); %#ok<AGROW>
        elseif diffV(k) * diffV(k + 1) <= 0
            denom = vGrid(k + 1) - vGrid(k);
            if abs(denom) < eps
                continue;
            end
            alpha = (vv - vGrid(k)) / denom;
            crossingX(end + 1, 1) = xGrid(k) + alpha * (xGrid(k + 1) - xGrid(k)); %#ok<AGROW>
        end
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

function xDomainUse = expand_sensor_domain_step06_local(xDomain, cfg, ~)
xDomainUse = xDomain;
if isfield(cfg, 'sensorDomainExpandMmTable') && ~isempty(cfg.sensorDomainExpandMmTable)
    % Reserved for parity with old Step06; current 20241106 setting is empty.
end
end

function wTotal = build_step06_filtered_weight_local(t, xComp, v, Tpl, queryGuardMm, cfg)
wEdge = build_edge_weight_step06_local(t, v, cfg.weightFloor);
wDomain = build_domain_soft_weight_step06_local(xComp, Tpl.x_domain, cfg.domainSoftMarginMm, cfg.weightFloor);
wQuery = build_query_guard_soft_weight_step06_local(xComp, [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))], ...
    queryGuardMm, cfg.domainSoftMarginMm, cfg.weightFloor);
wGradient = build_template_gradient_weight_step06_local(Tpl, xComp, cfg.domainSelectionMode, cfg.weightFloor);
wTotal = max(cfg.weightFloor, wEdge(:) .* wDomain(:) .* wQuery(:) .* wGradient(:));
if max(wTotal) > 0
    wTotal = max(cfg.weightFloor, wTotal ./ max(wTotal));
end
end

function wEdge = build_edge_weight_step06_local(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_step06_local(x, xDomain, marginMm, floorW)
if marginMm <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMm, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_step06_local(x, xDomain, queryGuardMm, marginMm, floorW)
if queryGuardMm <= 0 || marginMm <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMm - marginMm, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMm, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_step06_local(Tpl, x, domainSelectionMode, floorW)
if ~strcmpi(domainSelectionMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end

function vAligned = align_voltage_to_template_units_local(vRaw, Tpl)
vAligned = vRaw(:);
finiteV = isfinite(vAligned);
finiteTpl = isfinite(Tpl.v_grid(:));
if ~any(finiteV) || ~any(finiteTpl)
    return;
end
scaleV = median(abs(vAligned(finiteV)), 'omitnan');
scaleTpl = median(abs(Tpl.v_grid(finiteTpl) - Tpl.baseline), 'omitnan');
if isfinite(scaleV) && isfinite(scaleTpl) && scaleV > 0 && scaleTpl / scaleV > 100
    vAligned = vAligned * 1000;
end
end

function vRel = direct_bundle_relative_voltage_local(Bsrc, mask, Tpl)
vRaw = Bsrc.V(mask);
if isfield(Bsrc, 'forceDirectBaseCurve') && Bsrc.forceDirectBaseCurve && ...
        isfield(Bsrc, 'F0') && numel(Bsrc.F0) == numel(Bsrc.V)
    vRel = align_relative_voltage_to_template_units_local(vRaw, Tpl);
    return;
end
if isfield(Bsrc, 'F0') && numel(Bsrc.F0) == numel(Bsrc.V)
    f0Raw = Bsrc.F0(mask);
    finiteVF = isfinite(vRaw(:)) & isfinite(f0Raw(:));
    if nnz(finiteVF) >= 8
        vScale = median(abs(vRaw(finiteVF)), 'omitnan');
        f0Scale = median(abs(f0Raw(finiteVF)), 'omitnan');
        diffScale = sqrt(mean((vRaw(finiteVF) - f0Raw(finiteVF)).^2, 'omitnan'));
        scaleRef = max([vScale, f0Scale, eps]);
        if isfinite(vScale) && isfinite(f0Scale) && isfinite(diffScale) && ...
                vScale > 0 && f0Scale > 0 && ...
                max(vScale, f0Scale) / max(min(vScale, f0Scale), eps) < 10 && ...
                diffScale / scaleRef < 0.25
            vRel = align_relative_voltage_to_template_units_local(vRaw, Tpl);
            return;
        end
    end
end
vAligned = align_voltage_to_template_units_local(vRaw, Tpl);
vRel = vAligned - Tpl.baseline;
end

function vRel = align_relative_voltage_to_template_units_local(vRaw, Tpl)
vRel = vRaw(:);
finiteV = isfinite(vRel);
finiteTpl = isfinite(Tpl.v_grid(:));
if ~any(finiteV) || ~any(finiteTpl)
    return;
end
scaleV = median(abs(vRel(finiteV)), 'omitnan');
scaleTpl = median(abs(Tpl.v_grid(finiteTpl) - Tpl.baseline), 'omitnan');
if isfinite(scaleV) && isfinite(scaleTpl) && scaleV > 0 && scaleTpl / scaleV > 100
    vRel = vRel * 1000;
end
end

function tf = direct_bundle_uses_relative_voltage_local(Bsrc, mask)
tf = false;
if ~isfield(Bsrc, 'F0') || numel(Bsrc.F0) ~= numel(Bsrc.V)
    return;
end
if isfield(Bsrc, 'forceDirectBaseCurve') && Bsrc.forceDirectBaseCurve
    tf = true;
    return;
end
vRaw = Bsrc.V(mask);
f0Raw = Bsrc.F0(mask);
finiteVF = isfinite(vRaw(:)) & isfinite(f0Raw(:));
if nnz(finiteVF) < 8
    return;
end
vScale = median(abs(vRaw(finiteVF)), 'omitnan');
f0Scale = median(abs(f0Raw(finiteVF)), 'omitnan');
diffScale = sqrt(mean((vRaw(finiteVF) - f0Raw(finiteVF)).^2, 'omitnan'));
scaleRef = max([vScale, f0Scale, eps]);
tf = isfinite(vScale) && isfinite(f0Scale) && isfinite(diffScale) && ...
    vScale > 0 && f0Scale > 0 && ...
    max(vScale, f0Scale) / max(min(vScale, f0Scale), eps) < 10 && ...
    diffScale / scaleRef < 0.25;
end

function yRel = direct_bundle_relative_series_local(yRaw, Tpl)
yRel = align_relative_voltage_to_template_units_local(yRaw, Tpl);
end

function [baseX, baseF0] = make_direct_base_curve_local(xRaw, f0Raw)
xRaw = xRaw(:);
f0Raw = f0Raw(:);
valid = isfinite(xRaw) & isfinite(f0Raw);
if nnz(valid) < 8
    baseX = [];
    baseF0 = [];
    return;
end
xRaw = xRaw(valid);
f0Raw = f0Raw(valid);
[baseX, order] = sort(xRaw);
baseF0 = f0Raw(order);
[baseX, ~, groupIdx] = unique(baseX, 'stable');
if numel(baseX) < numel(baseF0)
    baseF0 = accumarray(groupIdx, baseF0, [], @median);
end
end

function etaLimit = direct_bundle_sensor_eta_limit_local(Bsrc, sensorPos, cfg)
etaLimit = cfg.etaFallbackLimitMm;
if ~isfield(Bsrc, 'sensor_eta_limit_mm') || isempty(Bsrc.sensor_eta_limit_mm)
    return;
end
raw = Bsrc.sensor_eta_limit_mm;
if isscalar(raw)
    etaLimit = raw;
elseif numel(raw) >= sensorPos
    etaLimit = raw(sensorPos);
end
if ~isfinite(etaLimit) || etaLimit <= 0
    etaLimit = cfg.etaFallbackLimitMm;
end
end

function bundle = build_gap_bundle_from_foundation_bundle_local(Bsrc, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
if isempty(Bsrc) || ~isstruct(Bsrc) || ~isfield(Bsrc, 'x') || ~isfield(Bsrc, 'v')
    error('Foundation Step05 bundle is missing x/v observations.');
end
sensorIds = cfg.analysisSensors(:).';
sensorId = Bsrc.sensor_id(:);
[present, sensorIndex] = ismember(sensorId, sensorIds);
X = Bsrc.x(:);
T = Bsrc.t(:);
Theta = Bsrc.theta(:);
V = Bsrc.v(:);
if isfield(Bsrc, 'fit_weight') && numel(Bsrc.fit_weight) == numel(V)
    W = Bsrc.fit_weight(:);
elseif isfield(Bsrc, 'legacy_fit_weight') && numel(Bsrc.legacy_fit_weight) == numel(V)
    W = Bsrc.legacy_fit_weight(:);
elseif isfield(Bsrc, 'template_weight') && numel(Bsrc.template_weight) == numel(V)
    W = Bsrc.template_weight(:);
else
    W = ones(size(V));
end
if isfield(Bsrc, 'template_v') && numel(Bsrc.template_v) == numel(V)
    F0 = Bsrc.template_v(:);
else
    F0 = nan(size(V));
end
if isfield(Bsrc, 'template_dv_dx') && numel(Bsrc.template_dv_dx) == numel(V)
    Fx = Bsrc.template_dv_dx(:);
else
    Fx = nan(size(V));
end
voltageScale = foundation_voltage_scale_to_template_local(Bsrc, Template, sensorIds);
V = V .* voltageScale;
F0 = F0 .* voltageScale;
Fx = Fx .* voltageScale;

valid = present(:) & isfinite(X) & isfinite(T) & isfinite(Theta) & ...
    isfinite(V) & isfinite(W) & isfinite(F0) & isfinite(Fx);
X = X(valid);
T = T(valid);
Theta = Theta(valid);
V = V(valid);
W = max(W(valid), cfg.weightFloor);
F0 = F0(valid);
Fx = Fx(valid);
sensorIndex = sensorIndex(valid);
S = sensorIds(sensorIndex).';

sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], ...
    'pointCount', 0, 'hasGapCorrection', false, 'directBaseX', [], ...
    'directBaseF0', [], 'etaLimitMm', cfg.etaFallbackLimitMm, ...
    'etaMedianMm', NaN, 'etaIqrMm', NaN), numel(sensorIds), 1);
g0BySensor = NaN(numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = sensorIndex == is;
    Tpl = get_template_sensor_local(Template, sid);
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    sensorInfo(is).pointCount = nnz(mask);
    sensorInfo(is).hasGapCorrection = is_gap_sensor_local(sid, cfg);
    [baseX, baseF0] = make_direct_base_curve_local(X(mask), F0(mask));
    sensorInfo(is).directBaseX = baseX;
    sensorInfo(is).directBaseF0 = baseF0;
    if sensorInfo(is).hasGapCorrection
        corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
        g0BySensor(is) = corr.g0Mm;
        sensorInfo(is).etaLimitMm = resolve_sensor_eta_limit_local(Tpl, corr, cfg);
        sensorInfo(is).etaMedianMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaMedianMm');
        sensorInfo(is).etaIqrMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaIqrMm');
    end
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
bundle.sensorIds = sensorIds;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.useDirectBaseCurve = false;
bundle.sourceVoltageScale = voltageScale;
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

function scale = foundation_voltage_scale_to_template_local(Bsrc, Template, sensorIds)
scale = 1;
if ~isfield(Bsrc, 'v') || isempty(Bsrc.v)
    return;
end
obs = abs(Bsrc.v(:));
obs = obs(isfinite(obs));
tplVals = [];
for sid = sensorIds(:).'
    Tpl = get_template_sensor_local(Template, sid);
    tplVals = [tplVals; abs(Tpl.v_grid(:))]; %#ok<AGROW>
end
tplVals = tplVals(isfinite(tplVals));
if isempty(obs) || isempty(tplVals)
    return;
end
obsRef = prctile(obs, 95);
tplRef = prctile(tplVals, 95);
if isfinite(obsRef) && isfinite(tplRef) && obsRef < 10 && tplRef > 100
    scale = 1000;
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
hasDirectBaseCurve = false;
g0BySensor = NaN(numel(sensorIds), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0, ...
    'hasGapCorrection', false, 'directBaseX', [], 'directBaseF0', []), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = sensorIndex == is;
    Tpl = get_template_sensor_local(Template, sid);
    hasGapCorrection = is_gap_sensor_local(sid, cfg);
    V(mask) = direct_bundle_relative_voltage_local(Bsrc, mask, Tpl);
    usesDirectBase = direct_bundle_uses_relative_voltage_local(Bsrc, mask);
    if hasGapCorrection
        corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
        F0(mask) = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask));
        Fx(mask) = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask), cfg.derivativeStepMm);
        g0BySensor(is) = corr.g0Mm;
        sensorInfo(is).etaLimitMm = resolve_sensor_eta_limit_local(Tpl, corr, cfg);
        sensorInfo(is).etaMedianMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaMedianMm');
        sensorInfo(is).etaIqrMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaIqrMm');
    else
        F0(mask) = eval_base_template_local(Tpl, X(mask));
        Fx(mask) = eval_base_template_derivative_local(Tpl, X(mask), cfg.derivativeStepMm);
        g0BySensor(is) = NaN;
        sensorInfo(is).etaLimitMm = cfg.etaFallbackLimitMm;
        sensorInfo(is).etaMedianMm = NaN;
        sensorInfo(is).etaIqrMm = NaN;
    end
    if usesDirectBase
        F0(mask) = direct_bundle_relative_series_local(Bsrc.F0(mask), Tpl);
        if isfield(Bsrc, 'Fx') && numel(Bsrc.Fx) == numel(Bsrc.V)
            Fx(mask) = direct_bundle_relative_series_local(Bsrc.Fx(mask), Tpl);
        end
        [baseX, baseF0] = make_direct_base_curve_local(X(mask), F0(mask));
        sensorInfo(is).directBaseX = baseX;
        sensorInfo(is).directBaseF0 = baseF0;
        sensorInfo(is).etaLimitMm = direct_bundle_sensor_eta_limit_local(Bsrc, is, cfg);
        hasDirectBaseCurve = true;
    end
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    sensorInfo(is).pointCount = nnz(mask);
    sensorInfo(is).hasGapCorrection = hasGapCorrection;
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
bundle.useDirectBaseCurve = hasDirectBaseCurve;
bundle.Template = Template;
bundle.CorrectedGapLibrary = CorrectedGapLibrary;
bundle.responseSurface = responseSurface;
if isfield(Bsrc, 'window_meta')
    bundle.sourceWindowMeta = Bsrc.window_meta;
end
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

function Template = align_direct_only_template_xcenter_local(Template, DynamicMap, cfg)
if ~cfg.alignDirectOnlyXCenter || isempty(cfg.directOnlySensors)
    return;
end
for sid = cfg.directOnlySensors(:).'
    idx = find([Template.Sensor.sensor_id] == sid, 1);
    if isempty(idx)
        continue;
    end
    dynamicXc = infer_dynamic_xcenter_local(DynamicMap, sid);
    if ~isfinite(dynamicXc)
        continue;
    end
    if ~isfield(Template.Sensor(idx), 'xc') || ~isfinite(Template.Sensor(idx).xc)
        Template.Sensor(idx).xc = dynamicXc;
        continue;
    end
    delta = Template.Sensor(idx).xc - dynamicXc;
    if abs(delta) <= 1e-9
        continue;
    end
    Template.Sensor(idx).x_grid = Template.Sensor(idx).x_grid(:) - delta;
    if isfield(Template.Sensor(idx), 'x_domain') && numel(Template.Sensor(idx).x_domain) == 2
        Template.Sensor(idx).x_domain = Template.Sensor(idx).x_domain - delta;
    end
    Template.Sensor(idx).xc = dynamicXc;
    Template.Sensor(idx).xc_reference_shift_mm = delta;
end
end

function bundle = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg, phaseRef)
if nargin < 6
    phaseRef = [];
end
usePhaseSafe = isstruct(phaseRef) && isfield(phaseRef, 'usable') && ...
    phaseRef.usable && cfg.phaseSafeExpansion;
X = []; T = []; V = []; W = []; Theta = []; S = []; sensorIndex = [];
F0 = []; Fx = []; g0BySensor = NaN(numel(cfg.analysisSensors), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0, ...
    'hasGapCorrection', false), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tpl = get_template_sensor_local(Template, sid);
    hasGapCorrection = is_gap_sensor_local(sid, cfg);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    xRaw = D.x_rel(:);
    tRaw = D.t(:);
    thetaRaw = D.theta(:);
    vRaw = D.V(:);
    vRawMv = D.V(:);
    wRaw = max(D.W(:), cfg.weightFloor);
    thresholdMv = resolve_threshold_local(Tpl, cfg);

    if hasGapCorrection
        corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
        f0Raw = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw);
        fxRaw = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw, cfg.derivativeStepMm);
    else
        corr = [];
        f0Raw = eval_base_template_local(Tpl, xRaw);
        fxRaw = eval_base_template_derivative_local(Tpl, xRaw, cfg.derivativeStepMm);
    end
    maskPulse = build_pulse_mask_local(vRawMv, thresholdMv, cfg.pulseSelectionMode);
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRawMv, tRaw, Tpl, thresholdMv, cfg);
    maskDomain = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRawMv) & isfinite(thetaRaw);
    mask = maskPulse & maskEffective & maskDomain;
    queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, mask, cfg);
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
    if hasGapCorrection
        g0BySensor(is) = corr.g0Mm;
    end
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = xDomain;
    sensorInfo(is).pointCount = nnz(mask);
    sensorInfo(is).queryGuardMm = queryGuardMm;
    sensorInfo(is).querySafeCount = nnz(maskQuerySafe);
    sensorInfo(is).phaseSafeExpansion = usePhaseSafe;
    sensorInfo(is).phaseSafeCount = nnz(mask & maskPhaseSafe);
    sensorInfo(is).hasGapCorrection = hasGapCorrection;
    if hasGapCorrection
        sensorInfo(is).etaLimitMm = resolve_sensor_eta_limit_local(Tpl, corr, cfg);
        sensorInfo(is).etaMedianMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaMedianMm');
        sensorInfo(is).etaIqrMm = resolve_sensor_eta_stat_local(Tpl, corr, 'etaIqrMm');
    else
        sensorInfo(is).etaLimitMm = cfg.etaFallbackLimitMm;
        sensorInfo(is).etaMedianMm = NaN;
        sensorInfo(is).etaIqrMm = NaN;
    end
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
fullEval = @(theta, eo, b, c) evaluate_static_warp_waveform_local(theta, eo, b, c, 'fixed');
seedTable = step07jcore.solve_gradient_displacement_vp_seed(bundle, eoCandidates, cfg, fullEval);
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
info.source = 'current_gradient_displacement_vp_seed';
info.reason = [fallbackReason, '_then_gradient_displacement_vp'];
info.sourceWindowId = windowId;
info.usedPreviousWindow = false;
end

function info = empty_phase_safe_reference_local()
info = struct('usable', false, 'EO', NaN, 'A', NaN, 'phi', NaN, ...
    'dx', NaN, 'deltaTauMm', [], 'weightedRmseMv', NaN, ...
    'source', 'none', 'reason', 'not_evaluated', ...
    'sourceWindowId', NaN, 'usedPreviousWindow', false);
end

function line = format_window_identification_progress_local(iw, numWindows, Wmap, bundle, fit, cfg)
pointCountText = 'NA';
if isstruct(bundle) && isfield(bundle, 'pointCount')
    pointCountText = sprintf('%d', bundle.pointCount);
end
line = sprintf(['Window %02d/%02d laps %s | %s: EO=%g, f=%.3f Hz, ' ...
    'A=%.4f mm, phi=%.3f rad, dx=%.4f mm, RMSE=%.2f mV, points=%s'], ...
    iw, numWindows, mat2str(Wmap.lap_range), cfg.mainModel, ...
    fit.EO, fit.freqHz, fit.amplitudeMm, wrap_to_pi_local(fit.phaseRad), ...
    fit.dxMm, fit.weightedRmseMv, pointCountText);

if isfield(fit, 'deltaGapMm') && ~isempty(fit.deltaGapMm)
    line = sprintf('%s, |dg|max=%.4f mm', line, max_abs_finite_local(fit.deltaGapMm));
end
if isfield(fit, 'deltaMuGapPerXMm') && ~isempty(fit.deltaMuGapPerXMm)
    line = sprintf('%s, |dmu|max=%.4f', line, max_abs_finite_local(fit.deltaMuGapPerXMm));
end
if isfield(fit, 'sensorEtaMm') && ~isempty(fit.sensorEtaMm)
    line = sprintf('%s, |eta|max=%.4f mm', line, max_abs_finite_local(fit.sensorEtaMm));
end

nearText = near_limit_text_local(fit, bundle, cfg);
if ~isempty(nearText)
    line = sprintf('%s  near-limit: %s', line, nearText);
end
end

function text = near_limit_text_local(fit, bundle, cfg)
items = {};
tol = 0.95;
if isfield(fit, 'amplitudeMm') && isfinite(fit.amplitudeMm) && ...
        cfg.amplitudeLimitMm > 0 && abs(fit.amplitudeMm) >= tol * cfg.amplitudeLimitMm
    items{end+1} = sprintf('A %.3g/%.3g', abs(fit.amplitudeMm), cfg.amplitudeLimitMm); %#ok<AGROW>
end
if isfield(fit, 'dxMm') && isfinite(fit.dxMm) && ...
        cfg.dxLimitMm > 0 && abs(fit.dxMm) >= tol * cfg.dxLimitMm
    items{end+1} = sprintf('dx %.3g/%.3g', abs(fit.dxMm), cfg.dxLimitMm); %#ok<AGROW>
end
if isfield(fit, 'deltaGapMm') && cfg.deltaGapLimitMm > 0 && ...
        max_abs_finite_local(fit.deltaGapMm) >= tol * cfg.deltaGapLimitMm
    items{end+1} = sprintf('|dg|max %.3g/%.3g', ...
        max_abs_finite_local(fit.deltaGapMm), cfg.deltaGapLimitMm); %#ok<AGROW>
end
if isfield(fit, 'deltaMuGapPerXMm') && cfg.deltaMuLimit > 0 && ...
        max_abs_finite_local(fit.deltaMuGapPerXMm) >= tol * cfg.deltaMuLimit
    items{end+1} = sprintf('|dmu|max %.3g/%.3g', ...
        max_abs_finite_local(fit.deltaMuGapPerXMm), cfg.deltaMuLimit); %#ok<AGROW>
end
if isfield(fit, 'sensorEtaMm') && ~isempty(fit.sensorEtaMm) && ...
        isfield(cfg, 'fitSensorEta') && cfg.fitSensorEta
    etaLimit = collect_sensor_eta_limit_from_bundle_local(bundle);
    eta = abs(fit.sensorEtaMm(:));
    etaLimit = etaLimit(:);
    n = min(numel(eta), numel(etaLimit));
    if n > 0
        ratio = eta(1:n) ./ max(etaLimit(1:n), eps);
        ratio = ratio(isfinite(ratio));
        if ~isempty(ratio) && max(ratio) >= tol
            items{end+1} = sprintf('eta %.1f%%', 100 * max(ratio)); %#ok<AGROW>
        end
    end
end
text = strjoin(items, ', ');
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
if ~isempty(directTrend) && height(directTrend) >= windowId
    directEO = table_value_by_names_local(directTrend, windowId, {'EO_id', 'EO'});
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

function eoCandidates = model_eo_candidates_local(modeName, selectedEO, selectionInfo, cfg)
eoCandidates = selectedEO;
if strcmpi(modeName, 'fixed') && strcmpi(cfg.fixedEoMode, 'direct') && ...
        isstruct(selectionInfo) && isfield(selectionInfo, 'directEO') && isfinite(selectionInfo.directEO)
    eoCandidates = round(selectionInfo.directEO);
end
eoCandidates = unique(round(eoCandidates(:).'), 'stable');
end

function fit = refine_static_warp_fit_local(bundle, seedTable, eoCandidates, cfg, modeName, fixedBaseFit)
if nargin < 6
    fixedBaseFit = [];
end
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'integerFrequencyHz', NaN, 'frequencyOffsetHz', 0, ...
    'integerWeightedRmseMv', NaN, 'rmseDeltaMv', 0, ...
    'frequencyRefinementExitFlag', NaN, 'frequencyRefinementIterations', NaN, ...
    'residualObjectiveMv2', inf, 'overshootPenalty', inf, ...
    'staticRegPenalty', inf, 'etaRegPenalty', inf, 'invalidPenalty', inf, 'penaltyObjectiveMv2', inf, ...
    'penaltyShare', NaN, ...
    'sensorEtaMm', [], 'deltaGapMm', [], 'deltaMuGapPerXMm', [], 'deltaTauMm', [], ...
    'VPred', [], 'uMm', [], 'score', inf), numel(eoCandidates), 1);
best = struct('score', inf);
cfgEval = attach_nested_fixed_base_local(cfg, modeName, fixedBaseFit);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    theta0 = initial_theta_local(seed.A, seed.phi, seed.dx, modeName, numel(bundle.sensorIds));
    theta0 = apply_fixed_kernel_initial_theta_local(theta0, fixedBaseFit, cfg, modeName, eo);
    theta0 = apply_gap_ls_initial_theta_local(theta0, eo, bundle, cfgEval, modeName);
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfgEval, modeName);
    thetaOpt = fminsearch(fun, theta0, cfg.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfg, modeName, numel(bundle.sensorIds), bundle.sensorIds);
    detail = evaluate_static_warp_waveform_local(thetaOpt, eo, bundle, cfgEval, modeName);
    if cfg.continuousFrequencyRefine && mode_has_gap_local(modeName)
        detail = refine_continuous_frequency_fit_local(bundle, detail, cfgEval, modeName);
    end
    candidate(i).EO = eo;
    candidate(i).freqHz = detail.freqHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    if isfield(detail, 'integerFrequencyHz')
        candidate(i).integerFrequencyHz = detail.integerFrequencyHz;
        candidate(i).frequencyOffsetHz = detail.frequencyOffsetHz;
        candidate(i).integerWeightedRmseMv = detail.integerWeightedRmseMv;
        candidate(i).rmseDeltaMv = detail.rmseDeltaMv;
        candidate(i).frequencyRefinementExitFlag = detail.frequencyRefinementExitFlag;
        candidate(i).frequencyRefinementIterations = detail.frequencyRefinementIterations;
    else
        candidate(i).integerFrequencyHz = eo * bundle.rotFreqMeanHz;
        candidate(i).integerWeightedRmseMv = detail.weightedRmseMv;
    end
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
if use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName) && ...
        all(isfinite([fixedBaseFit.EO, fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, fixedBaseFit.dxMm]))
    eo = round(fixedBaseFit.EO);
    thetaBase = initial_theta_local(fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, ...
        fixedBaseFit.dxMm, modeName, numel(bundle.sensorIds));
    thetaBase = bound_params_local(thetaBase, cfg, modeName, numel(bundle.sensorIds), bundle.sensorIds);
    detail = evaluate_static_warp_waveform_local(thetaBase, eo, bundle, cfgEval, modeName);
    baseCandidate = candidate(1);
    baseCandidate.EO = eo;
    baseCandidate.freqHz = eo * bundle.rotFreqMeanHz;
    baseCandidate.integerFrequencyHz = baseCandidate.freqHz;
    baseCandidate.frequencyOffsetHz = 0;
    baseCandidate.integerWeightedRmseMv = detail.weightedRmseMv;
    baseCandidate.rmseDeltaMv = 0;
    baseCandidate.frequencyRefinementExitFlag = NaN;
    baseCandidate.frequencyRefinementIterations = NaN;
    baseCandidate.amplitudeMm = detail.amplitudeMm;
    baseCandidate.phaseRad = detail.phaseRad;
    baseCandidate.dxMm = detail.dxMm;
    baseCandidate.weightedRmseMv = detail.weightedRmseMv;
    baseCandidate.plainRmseMv = detail.plainRmseMv;
    baseCandidate.residualObjectiveMv2 = detail.residualObjectiveMv2;
    baseCandidate.overshootPenalty = detail.overshootPenalty;
    baseCandidate.staticRegPenalty = detail.staticRegPenalty;
    baseCandidate.etaRegPenalty = detail.etaRegPenalty;
    baseCandidate.invalidPenalty = detail.invalidPenalty;
    baseCandidate.penaltyObjectiveMv2 = detail.penaltyObjectiveMv2;
    baseCandidate.penaltyShare = detail.penaltyShare;
    baseCandidate.sensorEtaMm = detail.sensorEtaMm;
    baseCandidate.deltaGapMm = detail.deltaGapMm;
    baseCandidate.deltaMuGapPerXMm = detail.deltaMuGapPerXMm;
    baseCandidate.deltaTauMm = detail.deltaTauMm;
    baseCandidate.VPred = detail.VPred;
    baseCandidate.uMm = detail.uMm;
    baseCandidate.score = detail.score;
    % The zero-gap nested candidate must be re-evaluated with the same
    % formal objective as every optimized gap candidate. Never transplant
    % the Foundation weighted score into a plain-RMSE candidate table.
    baseCandidate.deltaGapMm = zeros(numel(bundle.sensorIds), 1);
    baseCandidate.deltaMuGapPerXMm = zeros(numel(bundle.sensorIds), 1);
    baseCandidate.deltaTauMm = zeros(numel(bundle.sensorIds), 1);
    candidate(end+1) = baseCandidate; %#ok<AGROW>
    if baseCandidate.score <= best.score
        best = baseCandidate;
    end
end
fit = best;
fit.mode = modeName;
fit.nestedBaseMode = 'none';
if use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName)
    fit.nestedBaseMode = 'fixed_kernel_delta';
    fit.nestedBaseEO = fixedBaseFit.EO;
    fit.nestedBaseAmplitudeMm = fixedBaseFit.amplitudeMm;
    fit.nestedBaseDxMm = fixedBaseFit.dxMm;
end
fit.CandidateTable = sortrows(struct2table(rmfield(candidate, {'VPred','uMm','sensorEtaMm','deltaGapMm','deltaMuGapPerXMm','deltaTauMm'})), ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
fit = enforce_nested_fixed_floor_local(fit, fixedBaseFit, cfg, modeName, bundle);
fit = attach_eo_ambiguity_local(fit);
end

function fit = attach_eo_ambiguity_local(fit)
fit.distinctCandidateCount = NaN;
fit.secondBestEO = NaN;
fit.secondBestWeightedRmseMv = NaN;
fit.eoMarginMv = NaN;
fit.eoMarginPercent = NaN;
if ~isfield(fit, 'CandidateTable') || isempty(fit.CandidateTable)
    return;
end
T = fit.CandidateTable;
eoValues = unique(T.EO, 'stable');
bestScore = nan(numel(eoValues), 1);
bestRmse = nan(numel(eoValues), 1);
for i = 1:numel(eoValues)
    rows = find(T.EO == eoValues(i));
    [bestScore(i), k] = min(T.score(rows));
    bestRmse(i) = T.weightedRmseMv(rows(k));
end
[~, order] = sort(bestScore, 'ascend');
fit.distinctCandidateCount = numel(order);
if numel(order) < 2
    return;
end
fit.secondBestEO = eoValues(order(2));
fit.secondBestWeightedRmseMv = bestRmse(order(2));
fit.eoMarginMv = bestRmse(order(2)) - bestRmse(order(1));
fit.eoMarginPercent = 100 * fit.eoMarginMv / max(bestRmse(order(1)), eps);
end

function status = validate_formal_fit_quality_local(fit, cfg, windowId, modeName)
status = 'pass';
required = {'EO','freqHz','amplitudeMm','dxMm','weightedRmseMv','invalidPenalty'};
for i = 1:numel(required)
    if ~isfield(fit, required{i}) || ~isfinite(fit.(required{i}))
        error('Window %d %s fit lacks finite %s; refusing a formal EO result.', ...
            windowId, modeName, required{i});
    end
end
if fit.freqHz < cfg.freqSearchHz(1) || fit.freqHz > cfg.freqSearchHz(2)
    error('Window %d %s frequency %.6f Hz is outside %s.', ...
        windowId, modeName, fit.freqHz, mat2str(cfg.freqSearchHz));
end
tol = 0.995;
reasons = strings(0, 1);
if cfg.amplitudeLimitMm > 0 && fit.amplitudeMm >= tol * cfg.amplitudeLimitMm
    reasons(end+1) = sprintf('A=%.6g/%.6g mm', fit.amplitudeMm, cfg.amplitudeLimitMm); %#ok<AGROW>
end
if cfg.dxLimitMm > 0 && abs(fit.dxMm) >= tol * cfg.dxLimitMm
    reasons(end+1) = sprintf('|dx|=%.6g/%.6g mm', abs(fit.dxMm), cfg.dxLimitMm); %#ok<AGROW>
end
if mode_has_gap_local(modeName) && cfg.deltaGapLimitMm > 0 && ...
        isfield(fit, 'deltaGapMm') && ...
        max(abs(fit.deltaGapMm), [], 'omitnan') >= tol * cfg.deltaGapLimitMm
    reasons(end+1) = sprintf('|dg|max=%.6g/%.6g mm', ...
        max(abs(fit.deltaGapMm), [], 'omitnan'), cfg.deltaGapLimitMm); %#ok<AGROW>
end
if cfg.continuousFrequencyRefine && mode_has_gap_local(modeName) && ...
        cfg.frequencyRefineHalfWidthHz > 0 && isfield(fit, 'frequencyOffsetHz') && ...
        abs(fit.frequencyOffsetHz) >= tol * cfg.frequencyRefineHalfWidthHz
    reasons(end+1) = sprintf('|delta_f|=%.6g/%.6g Hz', ...
        abs(fit.frequencyOffsetHz), cfg.frequencyRefineHalfWidthHz); %#ok<AGROW>
end
if fit.invalidPenalty > 0
    reasons(end+1) = sprintf('invalidPenalty=%.6g', fit.invalidPenalty); %#ok<AGROW>
end
if ~isempty(reasons)
    error(['Window %d %s solution reached a hard boundary (%s). ' ...
        'The program refuses to report this EO as valid.'], ...
        windowId, modeName, strjoin(cellstr(reasons), ', '));
end
if isfield(fit, 'eoMarginPercent') && isfinite(fit.eoMarginPercent) && ...
        fit.eoMarginPercent < cfg.lowConfidenceEoMarginPercent
    status = 'low_eo_margin';
end
end

function cfgOut = attach_nested_fixed_base_local(cfg, modeName, fixedBaseFit)
cfgOut = cfg;
if use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName)
    cfgOut.NestedFixedBaseFit = fixedBaseFit;
end
end

function fit = enforce_nested_fixed_floor_local(fit, fixedBaseFit, cfg, modeName, bundle)
if isfield(cfg, 'finalObjective') && strcmpi(cfg.finalObjective, 'plain')
    % The formal package already inserts a zero-gap candidate evaluated on
    % the current bundle with the plain objective. A weighted Foundation
    % floor is not comparable and must not replace the selected result.
    return;
end
if ~use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName) || ...
        ~isfield(fixedBaseFit, 'weightedRmseMv') || ~isfinite(fixedBaseFit.weightedRmseMv)
    return;
end
if isfield(fit, 'weightedRmseMv') && isfinite(fit.weightedRmseMv) && ...
        fit.weightedRmseMv <= fixedBaseFit.weightedRmseMv
    return;
end
candidateTable = [];
if isfield(fit, 'CandidateTable')
    candidateTable = fit.CandidateTable;
end
fit = fixedBaseFit;
fit.mode = modeName;
fit.nestedBaseMode = 'fixed_kernel_delta_floor';
fit.nestedBaseEO = fixedBaseFit.EO;
fit.nestedBaseAmplitudeMm = fixedBaseFit.amplitudeMm;
fit.nestedBaseDxMm = fixedBaseFit.dxMm;
nSensor = numel(bundle.sensorIds);
fit.deltaGapMm = zeros(nSensor, 1);
fit.deltaMuGapPerXMm = zeros(nSensor, 1);
fit.deltaTauMm = zeros(nSensor, 1);
if ~isempty(candidateTable)
    fit.CandidateTable = candidateTable;
end
end

function theta0 = apply_fixed_kernel_initial_theta_local(theta0, fixedBaseFit, cfg, modeName, eo)
if ~use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName)
    return;
end
if ~all(isfinite([fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, fixedBaseFit.dxMm]))
    return;
end
if nargin >= 5 && isfinite(eo) && isfinite(fixedBaseFit.EO) && round(eo) ~= round(fixedBaseFit.EO)
    return;
end
theta0(1:3) = [fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, fixedBaseFit.dxMm];
end

function theta0 = apply_gap_ls_initial_theta_local(theta0, eo, bundle, cfg, modeName)
if ~mode_has_gap_local(modeName)
    return;
end
nSensor = numel(bundle.sensorIds);
[idxDg, ~, ~] = mode_param_indices_local(modeName, nSensor);
if isempty(idxDg)
    return;
end
A = theta0(1);
phi = theta0(2);
dx = theta0(3);
eta = step07jcore.fixed_sensor_eta_by_position(bundle.sensorIds, cfg);
nestedFixedBase = nested_fixed_base_from_cfg_local(cfg, modeName, bundle);
if ~isempty(nestedFixedBase)
    A = nestedFixedBase.amplitudeMm;
    phi = nestedFixedBase.phaseRad;
    dx = nestedFixedBase.dxMm;
    if isfield(nestedFixedBase, 'sensorEtaMm') && numel(nestedFixedBase.sensorEtaMm) >= nSensor
        eta = nestedFixedBase.sensorEtaMm(:);
        eta = eta(1:nSensor);
    end
end
dgSeed = zeros(1, nSensor);
hG = 1e-4;
for is = 1:nSensor
    if ~bundle_sensor_has_gap_correction_local(bundle, is)
        continue;
    end
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex == is;
    if nnz(mask) < 8
        continue;
    end
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    if ~isempty(nestedFixedBase) && isfield(nestedFixedBase, 'uMm') && numel(nestedFixedBase.uMm) == numel(bundle.V)
        u = nestedFixedBase.uMm(mask);
    else
        u = A .* sin(eo .* bundle.Theta(mask) + phi);
    end
    xEval = bundle.X(mask) - dx - u - eta(is);
    if ~isempty(nestedFixedBase)
        v0 = nestedFixedBase.VPred(mask);
        [vGp, ~] = eval_static_warp_delta_local(bundle.responseSurface, corr, hG, 0, 0, xEval);
        [vGm, ~] = eval_static_warp_delta_local(bundle.responseSurface, corr, -hG, 0, 0, xEval);
    else
        [v0, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval, bundle, is);
        [vGp, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, hG, 0, 0, xEval, bundle, is);
        [vGm, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, -hG, 0, 0, xEval, bundle, is);
    end
    r = bundle.V(mask) - v0;
    jg = (vGp - vGm) ./ (2 * hG);
    valid = isfinite(r) & isfinite(jg) & isfinite(bundle.W(mask));
    if nnz(valid) < 8
        continue;
    end
    ww = max(bundle.W(mask), cfg.weightFloor);
    denom = sum(ww(valid) .* jg(valid).^2);
    if isfinite(denom) && denom > eps
        dgSeed(is) = sum(ww(valid) .* jg(valid) .* r(valid)) / denom;
    end
end
dgSeed(~isfinite(dgSeed)) = 0;
theta0(idxDg) = max(min(dgSeed, cfg.deltaGapLimitMm), -cfg.deltaGapLimitMm);
end

function tf = use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName)
tf = mode_has_gap_local(modeName) && ...
    isfield(cfg, 'initializeGapModelsFromFixed') && cfg.initializeGapModelsFromFixed && ...
    isstruct(fixedBaseFit) && isfield(fixedBaseFit, 'EO') && isfinite(fixedBaseFit.EO);
end

function fixedBaseFit = nested_fixed_base_from_cfg_local(cfg, modeName, bundle)
fixedBaseFit = [];
if ~mode_has_gap_local(modeName) || ~isfield(cfg, 'NestedFixedBaseFit') || ...
        ~isstruct(cfg.NestedFixedBaseFit)
    return;
end
candidate = cfg.NestedFixedBaseFit;
if ~use_fixed_kernel_base_local(candidate, cfg, modeName) || ...
        ~isfield(candidate, 'VPred') || numel(candidate.VPred) ~= numel(bundle.V)
    return;
end
fixedBaseFit = candidate;
end

function fit = fixed_fit_from_direct_reference_local(directRef, windowId, bundle, ~)
fit = [];
if ~isstruct(directRef) || ~isfield(directRef, 'WindowResult') || numel(directRef.WindowResult) < windowId
    return;
end
wr = directRef.WindowResult(windowId);
if ~isfield(wr, 'Result') || isempty(wr.Result)
    return;
end
D = wr.Result;
required = {'EO_id','fn_id','A_id','phi_id_wrapped','dx_c_id','weighted_voltage_rmse'};
for i = 1:numel(required)
    if ~isfield(D, required{i})
        return;
    end
end
nSensor = numel(bundle.sensorIds);
eta = zeros(nSensor, 1);
if isfield(D, 'sensor_eta_id') && numel(D.sensor_eta_id) >= nSensor
    eta = D.sensor_eta_id(:);
    eta = eta(1:nSensor);
end
rmseMv = D.weighted_voltage_rmse;
if isfinite(rmseMv) && rmseMv < 1
    rmseMv = rmseMv * 1000;
end
fit = struct();
fit.EO = round(D.EO_id);
fit.freqHz = D.fn_id;
fit.amplitudeMm = D.A_id;
fit.phaseRad = D.phi_id_wrapped;
fit.dxMm = D.dx_c_id;
fit.weightedRmseMv = rmseMv;
fit.plainRmseMv = NaN;
if isfield(D, 'plain_voltage_rmse')
    fit.plainRmseMv = D.plain_voltage_rmse;
    if isfinite(fit.plainRmseMv) && fit.plainRmseMv < 1
        fit.plainRmseMv = fit.plainRmseMv * 1000;
    end
end
weightNormalizer = bundle_weight_normalizer_local(bundle, 0.05);
fit.weightNormalizer = weightNormalizer;
fit.residualObjectiveMv2 = (rmseMv.^2) * weightNormalizer;
fit.overshootPenalty = 0;
fit.staticRegPenalty = 0;
fit.etaRegPenalty = 0;
fit.invalidPenalty = 0;
fit.penaltyObjectiveMv2 = 0;
fit.penaltyShare = 0;
fit.sensorEtaMm = eta(:);
fit.deltaGapMm = zeros(nSensor, 1);
fit.deltaMuGapPerXMm = zeros(nSensor, 1);
fit.deltaTauMm = zeros(nSensor, 1);
if isfield(D, 'V_pred')
    fit.VPred = D.V_pred;
else
    fit.VPred = [];
end
if isfield(D, 'u_est')
    fit.uMm = D.u_est;
else
    fit.uMm = [];
end
fit.score = fit.residualObjectiveMv2;
fit.mode = 'fixed';
fit.CandidateTable = table(fit.EO, fit.freqHz, fit.amplitudeMm, fit.dxMm, rmseMv, ...
    'VariableNames', {'EO','freqHz','amplitudeMm','dxMm','weightedRmseMv'});
fit.source = 'direct_reference';
fit.metricNote = 'Fixed baseline copied from rotating_calibration direct result for apples-to-apples comparison.';
end

function fit = fixed_fit_from_foundation_step05_local(foundationWindow, bundle, cfg)
fit = [];
if ~isstruct(foundationWindow) || ~isfield(foundationWindow, 'FoundationResult') || ...
        ~isstruct(foundationWindow.FoundationResult)
    return;
end
D = foundationWindow.FoundationResult;
required = {'EO_id','fn_id','A_id','phi_id_wrapped','dx_c_id','weighted_voltage_rmse'};
for i = 1:numel(required)
    if ~isfield(D, required{i}) || ~isfinite(D.(required{i}))
        return;
    end
end
nSensor = numel(bundle.sensorIds);
eta = step07jcore.fixed_sensor_eta_by_position(bundle.sensorIds, cfg);
if isfield(D, 'sensor_eta_id') && numel(D.sensor_eta_id) >= nSensor
    eta = D.sensor_eta_id(:);
    sourceSensorIds = [];
    if isfield(D, 'sensor_ids')
        sourceSensorIds = D.sensor_ids(:).';
    end
    if numel(sourceSensorIds) == numel(eta) && all(ismember(bundle.sensorIds, sourceSensorIds))
        [~, loc] = ismember(bundle.sensorIds, sourceSensorIds);
        eta = eta(loc);
    else
        eta = eta(1:nSensor);
    end
end
rmseMv = D.weighted_voltage_rmse;
if isfinite(rmseMv) && rmseMv < 1
    rmseMv = 1000 * rmseMv;
end
plainRmseMv = NaN;
if isfield(D, 'plain_voltage_rmse') && isfinite(D.plain_voltage_rmse)
    plainRmseMv = D.plain_voltage_rmse;
    if plainRmseMv < 1
        plainRmseMv = 1000 * plainRmseMv;
    end
end
theta = [D.A_id, D.phi_id_wrapped, D.dx_c_id];
detail = evaluate_static_warp_waveform_local(theta, round(D.EO_id), bundle, cfg, 'fixed');
if isfield(detail, 'weightNormalizer') && isfinite(detail.weightNormalizer) && detail.weightNormalizer > 0
    weightNormalizer = detail.weightNormalizer;
else
    weightNormalizer = bundle_weight_normalizer_local(bundle, cfg.weightFloor);
end
fit = detail;
fit.EO = round(D.EO_id);
fit.freqHz = D.fn_id;
fit.amplitudeMm = D.A_id;
fit.phaseRad = D.phi_id_wrapped;
fit.dxMm = D.dx_c_id;
fit.sensorEtaMm = eta(:);
fit.deltaGapMm = zeros(nSensor, 1);
fit.deltaMuGapPerXMm = zeros(nSensor, 1);
fit.deltaTauMm = zeros(nSensor, 1);
fit.weightedRmseMv = rmseMv;
fit.plainRmseMv = plainRmseMv;
fit.weightNormalizer = weightNormalizer;
fit.residualObjectiveMv2 = (rmseMv.^2) * weightNormalizer;
fit.overshootPenalty = 0;
fit.staticRegPenalty = 0;
fit.etaRegPenalty = 0;
fit.invalidPenalty = 0;
fit.penaltyObjectiveMv2 = 0;
fit.penaltyShare = 0;
fit.score = fit.residualObjectiveMv2;
fit.mode = 'fixed';
fit.nestedBaseMode = 'foundation_step05_canonical';
fit.CandidateTable = table(fit.EO, fit.freqHz, fit.amplitudeMm, fit.dxMm, rmseMv, ...
    'VariableNames', {'EO','freqHz','amplitudeMm','dxMm','weightedRmseMv'});
fit.source = 'foundation_step05_result';
fit.metricNote = ['Canonical fixed no-gap kernel copied from the matching ' ...
    'btt_data_foundation Step05 result; gap models are nested updates from this fixed point.'];
if isfield(D, 'V_pred') && numel(D.V_pred) == numel(bundle.V)
    vPredDirect = D.V_pred(:);
    if isfinite(D.weighted_voltage_rmse) && D.weighted_voltage_rmse < 1
        vPredDirect = 1000 * vPredDirect;
    end
    fit.VPred = vPredDirect;
end
if isfield(D, 'u_est') && numel(D.u_est) == numel(bundle.V)
    fit.uMm = D.u_est(:);
end
end

function fit = refine_continuous_frequency_fit_local(bundle, integerFit, cfg, modeName)
halfWidthHz = cfg.frequencyRefineHalfWidthHz;
if ~(isfinite(halfWidthHz) && halfWidthHz > 0) || ...
        ~isstruct(integerFit) || ~isfinite(integerFit.EO) || ~isfinite(integerFit.freqHz)
    fit = integerFit;
    return;
end

nSensor = numel(bundle.sensorIds);
theta0 = initial_theta_local(integerFit.amplitudeMm, integerFit.phaseRad, ...
    integerFit.dxMm, modeName, nSensor);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
if ~isempty(idxDg) && numel(integerFit.deltaGapMm) >= nSensor
    theta0(idxDg) = integerFit.deltaGapMm(1:nSensor);
end
if ~isempty(idxDmu) && numel(integerFit.deltaMuGapPerXMm) >= nSensor
    theta0(idxDmu) = integerFit.deltaMuGapPerXMm(1:nSensor);
end
if ~isempty(idxDtau) && numel(integerFit.deltaTauMm) >= nSensor
    theta0(idxDtau) = integerFit.deltaTauMm(1:nSensor);
end

lbTheta = -inf(size(theta0));
ubTheta = inf(size(theta0));
lbTheta(1:3) = [0, -pi, -cfg.dxLimitMm];
ubTheta(1:3) = [cfg.amplitudeLimitMm, pi, cfg.dxLimitMm];
if ~isempty(idxDg)
    lbTheta(idxDg) = -cfg.deltaGapLimitMm;
    ubTheta(idxDg) = cfg.deltaGapLimitMm;
    gapMask = gap_sensor_mask_for_ids_local(cfg, bundle.sensorIds, nSensor);
    lbTheta(idxDg(~gapMask)) = 0;
    ubTheta(idxDg(~gapMask)) = 0;
end
if ~isempty(idxDmu)
    lbTheta(idxDmu) = -cfg.deltaMuLimit;
    ubTheta(idxDmu) = cfg.deltaMuLimit;
end
if ~isempty(idxDtau)
    lbTheta(idxDtau) = -cfg.deltaTauLimitMm;
    ubTheta(idxDtau) = cfg.deltaTauLimitMm;
end

f0 = integerFit.freqHz;
z0 = [theta0, f0];
lb = [lbTheta, max(cfg.freqSearchHz(1), f0 - halfWidthHz)];
ub = [ubTheta, min(cfg.freqSearchHz(2), f0 + halfWidthHz)];
fun = @(z) continuous_frequency_objective_local(z, integerFit.EO, bundle, cfg, modeName);
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
    'MaxIterations', 100, 'MaxFunctionEvaluations', 2000, ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-9, 'ScaleProblem', true);
[zOpt, ~, exitFlag, output] = fmincon(fun, z0, [], [], [], [], lb, ub, [], options);
thetaOpt = bound_params_local(zOpt(1:end-1), cfg, modeName, nSensor, bundle.sensorIds);
cfgEval = cfg;
cfgEval.activeContinuousFrequencyHz = zOpt(end);
detail = evaluate_static_warp_waveform_local(thetaOpt, integerFit.EO, bundle, cfgEval, modeName);

if ~isfinite(detail.score) || detail.score > integerFit.score * (1 + 1e-9)
    fit = integerFit;
    fit.integerFrequencyHz = integerFit.freqHz;
    fit.integerWeightedRmseMv = integerFit.weightedRmseMv;
    fit.frequencyOffsetHz = 0;
    fit.rmseDeltaMv = 0;
    fit.frequencyRefinementExitFlag = -10;
    fit.frequencyRefinementIterations = output.iterations;
    return;
end

fit = detail;
fit.EO = integerFit.EO;
fit.mode = modeName;
fit.integerFrequencyHz = integerFit.freqHz;
fit.integerWeightedRmseMv = integerFit.weightedRmseMv;
fit.frequencyOffsetHz = fit.freqHz - integerFit.freqHz;
fit.rmseDeltaMv = fit.weightedRmseMv - integerFit.weightedRmseMv;
fit.frequencyRefinementExitFlag = exitFlag;
fit.frequencyRefinementIterations = output.iterations;
end

function score = continuous_frequency_objective_local(z, eo, bundle, cfg, modeName)
cfg.activeContinuousFrequencyHz = z(end);
detail = evaluate_static_warp_waveform_local(z(1:end-1), eo, bundle, cfg, modeName);
score = detail.score;
end

function phaseArg = vibration_phase_argument_local(eo, bundle, cfg)
if isfield(cfg, 'activeContinuousFrequencyHz') && isfinite(cfg.activeContinuousFrequencyHz)
    if ~isfield(bundle, 'T') || numel(bundle.T) ~= numel(bundle.Theta)
        error('Continuous-frequency refinement requires bundle.T aligned with bundle.Theta.');
    end
    t = bundle.T(:);
    t0 = mean(t, 'omitnan');
    orderLineHz = eo * bundle.rotFreqMeanHz;
    detuningHz = cfg.activeContinuousFrequencyHz - orderLineHz;
    phaseArg = eo .* bundle.Theta + 2 * pi * detuningHz .* (t - t0);
else
    phaseArg = eo .* bundle.Theta;
end
end

function freqHz = vibration_frequency_hz_local(eo, bundle, cfg)
if isfield(cfg, 'activeContinuousFrequencyHz') && isfinite(cfg.activeContinuousFrequencyHz)
    freqHz = cfg.activeContinuousFrequencyHz;
else
    freqHz = eo * bundle.rotFreqMeanHz;
end
end

function score = bounded_objective_local(thetaRaw, eo, bundle, cfg, modeName)
theta = bound_params_local(thetaRaw, cfg, modeName, numel(bundle.sensorIds), bundle.sensorIds);
detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName);
score = detail.score + 1e4 * sum((thetaRaw(:) - theta(:)).^2);
end

function weightNormalizer = bundle_weight_normalizer_local(bundle, weightFloor)
valid = isfinite(bundle.V) & isfinite(bundle.W);
if nnz(valid) < 1
    weightNormalizer = max(numel(bundle.V), 1);
else
    weightNormalizer = max(numel(bundle.V), 1);
end
end

function theta0 = initial_theta_local(A, phi, dx, modeName, nSensor)
theta0 = [A, phi, dx];
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

function theta = bound_params_local(thetaRaw, cfg, modeName, nSensor, sensorIds)
if nargin < 5 || isempty(sensorIds)
    sensorIds = cfg.analysisSensors(1:nSensor);
end
theta = thetaRaw(:).';
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
theta(3) = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
need = 3 + mode_extra_count_local(modeName) * nSensor;
if numel(theta) < need
    theta(numel(theta)+1:need) = 0;
end
theta = theta(1:need);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
if ~isempty(idxDg)
    theta(idxDg) = max(min(theta(idxDg), cfg.deltaGapLimitMm), -cfg.deltaGapLimitMm);
    gapMask = gap_sensor_mask_for_ids_local(cfg, sensorIds, nSensor);
    theta(idxDg(~gapMask)) = 0;
end
if ~isempty(idxDmu)
    theta(idxDmu) = max(min(theta(idxDmu), cfg.deltaMuLimit), -cfg.deltaMuLimit);
    gapMask = gap_sensor_mask_for_ids_local(cfg, sensorIds, nSensor);
    theta(idxDmu(~gapMask)) = 0;
end
if ~isempty(idxDtau)
    theta(idxDtau) = max(min(theta(idxDtau), cfg.deltaTauLimitMm), -cfg.deltaTauLimitMm);
    gapMask = gap_sensor_mask_for_ids_local(cfg, sensorIds, nSensor);
    theta(idxDtau(~gapMask)) = 0;
end
end

function detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName)
theta = bound_params_local(theta, cfg, modeName, numel(bundle.sensorIds), bundle.sensorIds);
A = theta(1);
phi = theta(2);
dx = theta(3);
nSensor = numel(bundle.sensorIds);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
eta = step07jcore.fixed_sensor_eta_by_position(bundle.sensorIds, cfg);
dg = zeros(nSensor, 1);
dmu = zeros(nSensor, 1);
dtau = zeros(nSensor, 1);
if ~isempty(idxDg), dg = theta(idxDg).'; end
if ~isempty(idxDmu), dmu = theta(idxDmu).'; end
if ~isempty(idxDtau), dtau = theta(idxDtau).'; end
nestedFixedBase = nested_fixed_base_from_cfg_local(cfg, modeName, bundle);
if ~isempty(nestedFixedBase)
    if isfield(nestedFixedBase, 'sensorEtaMm') && numel(nestedFixedBase.sensorEtaMm) >= nSensor
        eta = nestedFixedBase.sensorEtaMm(:);
        eta = eta(1:nSensor);
    end
    u = A .* sin(vibration_phase_argument_local(eo, bundle, cfg) + phi);
    if isfield(nestedFixedBase, 'uMm') && numel(nestedFixedBase.uMm) == numel(bundle.V)
        uBase = nestedFixedBase.uMm(:);
    else
        uBase = nestedFixedBase.amplitudeMm .* sin( ...
            nestedFixedBase.EO .* bundle.Theta + nestedFixedBase.phaseRad);
    end
    vPred = nestedFixedBase.VPred(:);
    overshoot = zeros(size(bundle.V));
    for is = 1:numel(bundle.sensorIds)
        sid = bundle.sensorIds(is);
        mask = bundle.sensorIndex == is;
        xEval = bundle.X(mask) - dx - u(mask) - eta(is);
        xBase = bundle.X(mask) - nestedFixedBase.dxMm - uBase(mask) - eta(is);
        Tpl = get_template_sensor_local(bundle.Template, sid);
        if bundle_sensor_has_gap_correction_local(bundle, is)
            corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
            [vNoGap, overNoGap] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, 0, 0, 0, xEval, bundle, is);
            [vBaseNoGap, overBaseNoGap] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, 0, 0, 0, xBase, bundle, is);
            [deltaV, overshoot(mask)] = eval_static_warp_delta_local( ...
                bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval);
            vPred(mask) = vPred(mask) + (vNoGap - vBaseNoGap) + deltaV;
            overshoot(mask) = max(max(overshoot(mask), overNoGap), overBaseNoGap);
        else
            [vNoGap, overNoGap] = eval_direct_or_base_template_with_overshoot_local(bundle, is, Tpl, xEval);
            [vBaseNoGap, overBaseNoGap] = eval_direct_or_base_template_with_overshoot_local(bundle, is, Tpl, xBase);
            vPred(mask) = vPred(mask) + (vNoGap - vBaseNoGap);
            overshoot(mask) = max(overNoGap, overBaseNoGap);
            dg(is) = 0;
            dmu(is) = 0;
            dtau(is) = 0;
        end
    end
else
    u = A .* sin(vibration_phase_argument_local(eo, bundle, cfg) + phi);
    vPred = nan(size(bundle.V));
    overshoot = zeros(size(bundle.V));
    for is = 1:numel(bundle.sensorIds)
        sid = bundle.sensorIds(is);
        mask = bundle.sensorIndex == is;
        Tpl = get_template_sensor_local(bundle.Template, sid);
        xEval = bundle.X(mask) - dx - u(mask) - eta(is);
        if bundle_sensor_has_gap_correction_local(bundle, is)
            corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
            [vPred(mask), overshoot(mask)] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval, bundle, is);
        else
            [vPred(mask), overshoot(mask)] = eval_direct_or_base_template_with_overshoot_local(bundle, is, Tpl, xEval);
            dg(is) = 0;
            dmu(is) = 0;
            dtau(is) = 0;
        end
    end
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
pointCount = max(numel(bundle.V), 1);
if nnz(valid) < 8
    residualObj = inf;
    plainRmse = inf;
    overshootPenalty = inf;
    weightNormalizer = pointCount;
else
    res = bundle.V(valid) - vPred(valid);
    ww = max(bundle.W(valid), cfg.weightFloor);
    if strcmpi(cfg.finalObjective, 'plain')
        residualObj = sum(res.^2);
    else
        residualObj = sum(ww .* res.^2);
    end
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
    overshootPenalty = cfg.overshootPenaltyWeight * 1e6 * sum(ww .* overshoot(valid).^2);
    weightNormalizer = pointCount;
end
etaPenalty = 0;
regularizationPenalty = weightNormalizer * (cfg.staticRegWeightMv * sqrt(( ...
    sum((dg ./ max(cfg.deltaGapLimitMm, eps)).^2) + ...
    sum((dmu ./ max(cfg.deltaMuLimit, eps)).^2) + ...
    sum((dtau ./ max(cfg.deltaTauLimitMm, eps)).^2)) / max(count_gap_corrected_sensors_local(bundle), 1))) ^ 2;
invalidPenalty = 1e12 * nnz(~valid);
objective = residualObj + overshootPenalty + regularizationPenalty + etaPenalty + invalidPenalty;
weightedRmse = sqrt(objective / weightNormalizer);
penaltyObjective = overshootPenalty + regularizationPenalty + etaPenalty + invalidPenalty;
if isfinite(objective) && objective > 0
    penaltyShare = penaltyObjective / objective;
else
    penaltyShare = NaN;
end
detail = struct();
detail.EO = eo;
detail.freqHz = vibration_frequency_hz_local(eo, bundle, cfg);
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx;
detail.sensorEtaMm = eta(:);
detail.deltaGapMm = dg(:);
detail.deltaMuGapPerXMm = dmu(:);
detail.deltaTauMm = dtau(:);
detail.weightedRmseMv = weightedRmse;
detail.plainRmseMv = plainRmse;
detail.weightNormalizer = weightNormalizer;
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
if ~isempty(nestedFixedBase)
    detail.nestedDeltaMode = 'additive_fixed_vpred';
    detail.nestedBaseMaxAbsDiffMv = max(abs(vPred(:) - nestedFixedBase.VPred(:)), [], 'omitnan');
else
    detail.nestedDeltaMode = 'direct_template';
    detail.nestedBaseMaxAbsDiffMv = NaN;
end
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

function [model, overshoot] = eval_static_warp_template_local(Tpl, responseSurface, corr, dg, dmu, dtau, xOpr, bundle, sensorPos)
if nargin < 8
    bundle = [];
    sensorPos = NaN;
end
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
[vLow, overDirect] = eval_direct_base_curve_local(bundle, sensorPos, xWarpRaw);
if isempty(vLow)
    vLow = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xLow, 'pchip', NaN);
    overDirect = zeros(size(overshootLow));
end
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(max(overDyn, overBase), overshootLow), overDirect);
end

function [deltaV, overshoot] = eval_static_warp_delta_local(responseSurface, corr, dg, dmu, dtau, xOpr)
xDyn = xOpr(:) - dtau;
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xDyn);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr(:));
deltaV = corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(overDyn, overBase);
end

function [model, overshoot] = eval_direct_or_base_template_with_overshoot_local(bundle, sensorPos, Tpl, xOpr)
[model, overshoot] = eval_direct_base_curve_local(bundle, sensorPos, xOpr);
if isempty(model)
    [model, overshoot] = eval_base_template_with_overshoot_local(Tpl, xOpr);
end
end

function [model, overshoot] = eval_direct_base_curve_local(bundle, sensorPos, xOpr)
model = [];
overshoot = [];
if ~isstruct(bundle) || ~isfield(bundle, 'useDirectBaseCurve') || ~bundle.useDirectBaseCurve || ...
        ~isfield(bundle, 'sensorInfo') || sensorPos < 1 || sensorPos > numel(bundle.sensorInfo)
    return;
end
info = bundle.sensorInfo(sensorPos);
if ~isfield(info, 'directBaseX') || ~isfield(info, 'directBaseF0') || ...
        numel(info.directBaseX) < 8 || numel(info.directBaseF0) ~= numel(info.directBaseX)
    return;
end
xRaw = xOpr(:);
xLo = min(info.directBaseX(:));
xHi = max(info.directBaseX(:));
xClip = min(max(xRaw, xLo), xHi);
model = interp1(info.directBaseX(:), info.directBaseF0(:), xClip, 'pchip', NaN);
overshoot = max(xLo - xRaw, 0) + max(xRaw - xHi, 0);
end

function [model, overshoot] = eval_base_template_with_overshoot_local(Tpl, xOpr)
xRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xEval = min(max(xRaw, xLo), xHi);
overshoot = max(xLo - xRaw, 0) + max(xRaw - xHi, 0);
model = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xEval, 'pchip', NaN);
end

function model = eval_base_template_local(Tpl, xOpr)
[model, ~] = eval_base_template_with_overshoot_local(Tpl, xOpr);
end

function dMdx = eval_base_template_derivative_local(Tpl, xOpr, h)
mp = eval_base_template_local(Tpl, xOpr + h);
mm = eval_base_template_local(Tpl, xOpr - h);
dMdx = (mp - mm) ./ (2 * h);
end

function [model, overshoot] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr)
xLowRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xLowRaw, xLo), xHi);
overshootLow = max(xLo - xLowRaw, 0) + max(xLowRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xLow, 'pchip', NaN);
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
B0 = interp_response_coeff_finite_local(responseSurface.xGrid, responseSurface.coeff(:, 1), xq);
B1 = interp_response_coeff_finite_local(responseSurface.xGrid, responseSurface.coeff(:, 2), xq);
B2 = interp_response_coeff_finite_local(responseSurface.xGrid, responseSurface.coeff(:, 3), xq);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function yq = interp_response_coeff_finite_local(xGrid, yGrid, xq)
x = xGrid(:);
y = yGrid(:);
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 2
    yq = NaN(size(xq));
    return;
end
xOk = x(ok);
yOk = y(ok);
xqClip = min(max(xq, min(xOk)), max(xOk));
yq = interp1(xOk, yOk, xqClip, 'linear', NaN);
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
if isfield(Tpl, 'x_grid') && isfield(Tpl, 'v_grid')
    [Tpl.x_grid, Tpl.v_grid, Tpl.baseline] = clean_template_curve_local(Tpl.x_grid, Tpl.v_grid, Tpl.baseline);
end
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1);
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function Template = sanitize_template_bundle_local(Template)
for is = 1:numel(Template.Sensor)
    if isfield(Template.Sensor(is), 'x_grid') && isfield(Template.Sensor(is), 'v_grid')
        scale = infer_template_voltage_scale_local(Template.Sensor(is));
        Template.Sensor(is).v_grid = Template.Sensor(is).v_grid(:) * scale;
        if isfield(Template.Sensor(is), 'dv_dx') && ~isempty(Template.Sensor(is).dv_dx)
            Template.Sensor(is).dv_dx = Template.Sensor(is).dv_dx(:) * scale;
        end
        if isfield(Template.Sensor(is), 'threshold') && isfinite(Template.Sensor(is).threshold)
            Template.Sensor(is).threshold = Template.Sensor(is).threshold * scale;
        end
        [Template.Sensor(is).x_grid, Template.Sensor(is).v_grid, Template.Sensor(is).baseline] = ...
            clean_template_curve_local(Template.Sensor(is).x_grid, Template.Sensor(is).v_grid, Template.Sensor(is).baseline);
        Template.Sensor(is).voltageUnit = 'mV';
    end
end
Template.voltageUnit = 'mV';
end

function [xGrid, vGrid, baseline] = clean_template_curve_local(xGrid, vGrid, baseline)
xGrid = xGrid(:);
vGrid = vGrid(:);
finiteMask = isfinite(xGrid) & isfinite(vGrid);
if ~any(finiteMask)
    error('Template curve has no finite points.');
end
xGrid = xGrid(finiteMask);
vGrid = vGrid(finiteMask);
[xGrid, ia] = unique(xGrid, 'stable');
vGrid = vGrid(ia);
if nargin < 3 || ~isfinite(baseline)
    baseline = 0;
end
if numel(xGrid) < 2
    error('Template curve needs at least two finite support points.');
end
end

function scale = infer_template_voltage_scale_local(Tpl)
scale = 1;
if isfield(Tpl, 'voltageUnit') && ~isempty(Tpl.voltageUnit)
    unit = lower(strtrim(char(string(Tpl.voltageUnit))));
    if strcmp(unit, 'v')
        scale = 1000;
        return;
    elseif strcmp(unit, 'mv')
        scale = 1;
        return;
    end
end
v = abs(Tpl.v_grid(:));
v = v(isfinite(v));
if isempty(v)
    return;
end
if median(v, 'omitnan') < 20
    scale = 1000;
end
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

function tf = is_gap_sensor_local(sid, cfg)
tf = isfield(cfg, 'gapSensors') && any(cfg.gapSensors(:).' == sid);
end

function mask = gap_sensor_mask_for_ids_local(cfg, sensorIds, nSensor)
mask = true(1, nSensor);
if nargin >= 2 && numel(sensorIds) >= nSensor && isfield(cfg, 'gapSensors')
    mask = ismember(sensorIds(1:nSensor), cfg.gapSensors);
end
end

function tf = bundle_sensor_has_gap_correction_local(bundle, sensorIndex)
tf = true;
if isfield(bundle, 'sensorInfo') && numel(bundle.sensorInfo) >= sensorIndex && ...
        isfield(bundle.sensorInfo(sensorIndex), 'hasGapCorrection')
    tf = logical(bundle.sensorInfo(sensorIndex).hasGapCorrection);
end
end

function n = count_gap_corrected_sensors_local(bundle)
n = numel(bundle.sensorIds);
if isfield(bundle, 'sensorInfo') && isfield(bundle.sensorInfo, 'hasGapCorrection')
    n = sum([bundle.sensorInfo.hasGapCorrection]);
end
n = max(n, 1);
end

function S = rmfield_if_present_local(S, fieldNames)
if ischar(fieldNames) || isstring(fieldNames)
    fieldNames = cellstr(fieldNames);
end
for i = 1:numel(fieldNames)
    name = fieldNames{i};
    if isfield(S, name)
        S = rmfield(S, name);
    end
end
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
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
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
end

function mask = build_dynamic_effective_mask_local(x, vMv, t, Tpl, thresholdMv, cfg)
finite = isfinite(x) & isfinite(vMv) & isfinite(t);
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

gTime = zeros(size(vMv));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(vMv(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if isfinite(gTimeMax) && gTimeMax > 0
    maskTime = gTime >= cfg.dynamicTimeGradientMinRatio * gTimeMax;
else
    maskTime = false(size(vMv));
end

peakQ = min(max(cfg.dynamicPeakQuantile, 0), 100);
peakLevel = prctile(vMv(finite), peakQ);
maskPeak = vMv >= max(thresholdMv, peakLevel);

mask = finite & (maskTpl | maskTime | maskPeak) & vMv >= 0.5 * thresholdMv;
if nnz(mask) < 8
    mask = finite & vMv >= thresholdMv;
end
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
eta = 0;
if isfield(cfg, 'fixedSensorEtaMm') && numel(cfg.fixedSensorEtaMm) >= sensorLocalIndex
    eta = cfg.fixedSensorEtaMm(sensorLocalIndex);
end
uEst = phaseRef.A .* sin(phaseRef.EO .* theta(:) + phaseRef.phi);
xQuery = x(:) - phaseRef.dx - eta - dtau - uEst;
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

function row = make_trend_row_local(wr, directTrend, cfg)
directEO = NaN; directFreq = NaN; directAmp = NaN; directRmseMv = NaN;
if ~isempty(directTrend) && height(directTrend) >= wr.windowId
    directEO = table_value_by_names_local(directTrend, wr.windowId, {'EO_id', 'EO'});
    directFreq = table_value_by_names_local(directTrend, wr.windowId, {'fn_id', 'FrequencyHz'});
    directAmp = table_value_by_names_local(directTrend, wr.windowId, {'A_id', 'AmplitudeMM'});
    directRmseMv = direct_rmse_mv_from_trend_local(directTrend, wr.windowId);
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
    row.fixed_objective_rmse_mV = F.fixed.weightedRmseMv;
    row.fixed_plain_rmse_mV = F.fixed.plainRmseMv;
end
if isfield(F, 'gap_only')
    row.gap_EO = F.gap_only.EO;
    row.gap_frequency_hz = F.gap_only.freqHz;
    row.gap_amplitude_mm = F.gap_only.amplitudeMm;
    row.gap_dx_mm = F.gap_only.dxMm;
    row.gap_rmse_mV = F.gap_only.weightedRmseMv;
    row.gap_objective_rmse_mV = F.gap_only.weightedRmseMv;
    row.gap_plain_rmse_mV = F.gap_only.plainRmseMv;
    row.gap_mean_delta_gap_mm = mean(F.gap_only.deltaGapMm, 'omitnan');
    row.gap_distinct_candidate_count = F.gap_only.distinctCandidateCount;
    row.gap_second_best_EO = F.gap_only.secondBestEO;
    row.gap_second_best_rmse_mV = F.gap_only.secondBestWeightedRmseMv;
    row.gap_EO_margin_mV = F.gap_only.eoMarginMv;
    row.gap_EO_margin_percent = F.gap_only.eoMarginPercent;
    row.gap_quality_status = string(F.gap_only.qualityStatus);
    if isfield(F.gap_only, 'frequencyOffsetHz')
        row.gap_frequency_offset_hz = F.gap_only.frequencyOffsetHz;
    else
        row.gap_frequency_offset_hz = 0;
    end
    row.gap_frequency_offset_at_bound = cfg.continuousFrequencyRefine && ...
        cfg.frequencyRefineHalfWidthHz > 0 && ...
        abs(row.gap_frequency_offset_hz) >= 0.995 * cfg.frequencyRefineHalfWidthHz;
end
if isfield(F, 'gap_tilt')
    row.tilt_EO = F.gap_tilt.EO;
    row.tilt_frequency_hz = F.gap_tilt.freqHz;
    row.tilt_amplitude_mm = F.gap_tilt.amplitudeMm;
    row.tilt_dx_mm = F.gap_tilt.dxMm;
    row.tilt_rmse_mV = F.gap_tilt.weightedRmseMv;
    row.tilt_objective_rmse_mV = F.gap_tilt.weightedRmseMv;
    row.tilt_plain_rmse_mV = F.gap_tilt.plainRmseMv;
    row.tilt_mean_delta_gap_mm = mean(F.gap_tilt.deltaGapMm, 'omitnan');
    row.tilt_mean_delta_mu = mean(F.gap_tilt.deltaMuGapPerXMm, 'omitnan');
    row.tilt_mean_eta_mm = mean(F.gap_tilt.sensorEtaMm, 'omitnan');
    row.tilt_max_abs_eta_mm = max_abs_finite_local(F.gap_tilt.sensorEtaMm);
end
end

function value = direct_rmse_mv_from_trend_local(T, rowIdx)
value = table_value_by_names_local(T, rowIdx, {'weighted_voltage_rmse_mV', 'WeightedVoltageRMSEmV'});
if isfinite(value)
    return;
end
valueV = table_value_by_names_local(T, rowIdx, {'weighted_voltage_rmse', 'WeightedVoltageRMSE'});
if isfinite(valueV)
    value = 1000 * valueV;
end
end

function value = residual_rmse_mv_local(fit, bundle)
pointCount = max(numel(bundle.V), 1);
value = sqrt(fit.residualObjectiveMv2 / pointCount);
end

function value = table_value_by_names_local(T, rowIdx, names)
value = NaN;
if isempty(T) || rowIdx > height(T)
    return;
end
for i = 1:numel(names)
    name = names{i};
    if ~ismember(name, T.Properties.VariableNames)
        continue;
    end
    raw = T.(name)(rowIdx);
    if iscell(raw)
        raw = raw{1};
    end
    if isstring(raw) || ischar(raw)
        raw = str2double(raw);
    end
    if isnumeric(raw) || islogical(raw)
        value = double(raw(1));
    end
    return;
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
    "gap_only", "gap_EO", "gap_frequency_hz", "gap_amplitude_mm", "gap_rmse_mV"
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
windowCount = zeros(n, 1);
dominantWindowCount = zeros(n, 1);
dominantWindowFraction = NaN(n, 1);
excludedMixedEoWindows = zeros(n, 1);
rawAllEoMeanFreq = NaN(n, 1);
rawAllEoMeanAmp = NaN(n, 1);
for i = 1:n
    method(i) = spec{i, 1};
    eoValues = Trend.(spec{i, 2});
    freqValues = Trend.(spec{i, 3});
    ampValues = Trend.(spec{i, 4});
    rmseValues = Trend.(spec{i, 5});
    dominantEO(i) = mode_finite_local(eoValues);
    validEo = isfinite(eoValues);
    dominantMask = validEo & eoValues == dominantEO(i);
    windowCount(i) = nnz(validEo);
    dominantWindowCount(i) = nnz(dominantMask);
    excludedMixedEoWindows(i) = windowCount(i) - dominantWindowCount(i);
    dominantWindowFraction(i) = dominantWindowCount(i) / max(windowCount(i), 1);
    meanFreq(i) = mean(freqValues(dominantMask), 'omitnan');
    stdFreq(i) = std(freqValues(dominantMask), 'omitnan');
    meanAmp(i) = mean(ampValues(dominantMask), 'omitnan');
    meanRmse(i) = mean(rmseValues(dominantMask), 'omitnan');
    medianRmse(i) = median(rmseValues(dominantMask), 'omitnan');
    rawAllEoMeanFreq(i) = mean(freqValues(validEo), 'omitnan');
    rawAllEoMeanAmp(i) = mean(ampValues(validEo), 'omitnan');
end
bestWindowIndex = repmat(bestIdx, n, 1);
Summary = table(method, dominantEO, windowCount, dominantWindowCount, ...
    dominantWindowFraction, excludedMixedEoWindows, meanFreq, stdFreq, meanAmp, ...
    meanRmse, medianRmse, rawAllEoMeanFreq, rawAllEoMeanAmp, bestWindowIndex, ...
    'VariableNames', {'method','dominant_EO','window_count','dominant_EO_window_count', ...
    'dominant_EO_window_fraction','excluded_mixed_EO_windows', ...
    'mean_frequency_hz','std_frequency_hz','mean_amplitude_mm', ...
    'mean_rmse_mV','median_rmse_mV','raw_all_EO_mean_frequency_hz', ...
    'raw_all_EO_mean_amplitude_mm','best_free_gap_window'});
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
etaVec = step07jcore.fixed_sensor_eta_by_position(bundle.sensorIds, cfg);
if isfield(fit, 'sensorEtaMm') && numel(fit.sensorEtaMm) >= numel(cfg.analysisSensors)
    etaVec = fit.sensorEtaMm(:);
end
fig = figure('Name', 'Step07J best-window dynamic collapse', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 4.8 * numel(cfg.analysisSensors)]);
tiledlayout(numel(cfg.analysisSensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = bundle.sensorIndex == is;
    xRaw = bundle.X(mask) - fit.dxMm - etaVec(is);
    xComp = xRaw - fit.uMm(mask);
    Tpl = get_template_sensor_local(bundle.Template, sid);
    dg = fit.deltaGapMm(is);
    dmu = fit.deltaMuGapPerXMm(is);
    dtau = fit.deltaTauMm(is);
    if bundle_sensor_has_gap_correction_local(bundle, is)
        corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
        [modelNo, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, dg, dmu, dtau, xRaw, bundle, is);
        [modelDyn, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, dg, dmu, dtau, xComp, bundle, is);
        titleAfter = sprintf('CH%d after, dg %.4f, dmu %.4f, dtau %.4f', sid, dg, dmu, dtau);
    else
        [modelNo, ~] = eval_direct_or_base_template_with_overshoot_local(bundle, is, Tpl, xRaw);
        [modelDyn, ~] = eval_direct_or_base_template_with_overshoot_local(bundle, is, Tpl, xComp);
        titleAfter = sprintf('CH%d after, direct-only', sid);
    end
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
    title(titleAfter);
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

function assert_capacitive_gap_sensors_local(sensorIds, overrideEnv)
allowed = [5 7];
if all(ismember(sensorIds(:).', allowed))
    return;
end
if parse_logical_env_local(overrideEnv, false)
    warning('Using non-capacitive sensors %s in Step07J gap fitting because %s is enabled.', ...
        mat2str(sensorIds), overrideEnv);
    return;
end
error(['The 20241106 gap-prior branch uses a capacitive-probe gap library for sensors [5 7]. ' ...
    'Requested sensors %s include non-capacitive probes. Set %s=1 only for diagnostics.'], ...
    mat2str(sensorIds), overrideEnv);
end

