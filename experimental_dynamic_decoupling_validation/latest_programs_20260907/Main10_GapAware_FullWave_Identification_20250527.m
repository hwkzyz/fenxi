%% Step07J: Nested high-speed static warp models
% Current route, 2026-06-30:
%   1) Step07J builds the Step03-style direct-template point set locally
%      from the current Template and DynamicMap.
%   2) No rotating_calibration Step03 result, EO trend, or stored bundle is
%      used as prior information in the gap-prior main route.
%   3) The default main route computes only the paper-facing gap_fixed_tilt
%      model (internal field name: gap_only).
%      Other methods are computed only with STEP07J_RUN_MODE=comparison.
%
% The main gap_fixed_tilt route jointly identifies EO, vibration and static
% gap increment inside the direct-template + VP Top-K EO candidate set with
% the full waveform objective. The low-speed calibrated sensor tilt slope is
% kept fixed in the corrected gap library; no residual high-speed dmu is
% identified in the paper-facing main route. VP narrows
% the EO candidates but does not fix the final EO. The gap-increment search
% half-width is set by Delta g_proj = epsilon_perp / S_g and
% Delta g_loc = alpha * Delta g_proj; STEP07J_DELTA_G_LIMIT_MM is only the
% projection cap/fallback unless STEP07J_USE_DELTA_G_PROJ=0 is used for an
% ablation run. In comparison mode this script also computes nested static
% correction models on top of the low-speed template:
%
%   fixed             : no clearance increment
%   gap_only          : dg_s with low-speed calibrated mu fixed (gap_fixed_tilt)
%   gap_tilt          : dg_s + residual dmu_s*(x-tau_s)
%
% Paper-facing main output uses gap_only/gap_fixed_tilt. gap_tilt is retained
% only in the separate comparison run to quantify sensitivity to a residual
% high-speed tilt increment. gap_tilt_shift is not part of the current
% comparison.
%
% This is the current paper-facing identification route.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20250527();
rootDir = packageCfg.paths.root;
C0 = CaseConfig();
outDir = packageCfg.paths.gapResults;
figDir = fullfile(outDir, 'figures_step07j_nested_static_warp_vp_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = C0.bladeId;
cfg.analysisSensors = C0.sensorIds;
if isfield(C0, 'gapSensors')
    cfg.gapSensors = double(C0.gapSensors(:).');
else
    cfg.gapSensors = cfg.analysisSensors;
end
cfg.gapActiveSensorMask = ismember(double(cfg.analysisSensors(:).'), cfg.gapSensors);
cfg.directOnlySensors = cfg.analysisSensors(~cfg.gapActiveSensorMask);
cfg.freqSearchHz = packageCfg.frequency.searchHz;
cfg.structureMode = packageCfg.frequency.structureMode;
cfg.candidateMethod = packageCfg.frequency.candidateMethod;
cfg.eoPad = 2;
cfg.vpTopK = parse_positive_integer_env_local('STEP07J_TOP_K_EO', packageCfg.frequency.vpTopK);
cfg.maxWindows = parse_positive_integer_env_local('STEP07J_MAX_WINDOWS', inf);
cfg.maxPointsPerWindow = parse_positive_integer_env_local('STEP07J_MAX_POINTS_PER_WINDOW', inf);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP07J_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DX_LIMIT_MM', 0.35);
cfg.deltaDxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_DX_LIMIT_MM', 0.05);
cfg.dxReferenceMode = 'seed_median';
cfg.causalDxRef = default_causal_dx_reference_cfg_local(cfg);
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_LIMIT_MM', 0.25);
cfg.deltaGapPriorScaleMm = cfg.deltaGapLimitMm;
cfg.useDeltaGapProjectionLimit = parse_logical_env_local('STEP07J_USE_DELTA_G_PROJ', true);
cfg.deltaGapProjAlpha = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_ALPHA', 2.5);
cfg.deltaGapProjMinLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_MIN_LIMIT_MM', 0.02);
cfg.deltaGapProjMaxLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_MAX_LIMIT_MM', cfg.deltaGapLimitMm);
cfg.deltaGapProjFallbackLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_FALLBACK_LIMIT_MM', cfg.deltaGapLimitMm);
cfg.deltaGapDerivativeStepMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_DERIV_STEP_MM', 1e-4);
cfg.deltaGapProjSensitivityFloorMvPerMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_SG_FLOOR_MV_PER_MM', 1e-6);
cfg.deltaGapProjEtaFloor = parse_nonnegative_numeric_env_local('STEP07J_DELTA_G_PROJ_ETA_FLOOR', 0.02);
% The calibrated low-speed library already carries the per-sensor tilt
% slope corr.muGapPerXMm.  Step07J therefore only estimates a residual
% high-speed tilt increment dmu.  Keep this increment tightly constrained;
% otherwise dmu can explain EO-synchronous waveform changes and pull the
% identified vibration amplitude downward.
cfg.deltaMuLimit = parse_nonnegative_numeric_env_local('STEP07J_DELTA_MU_LIMIT', 0.005);
cfg.deltaTauLimitMm = parse_nonnegative_numeric_env_local('STEP07J_DELTA_TAU_LIMIT_MM', 0.16);
cfg.etaFallbackLimitMm = parse_nonnegative_numeric_env_local('STEP07J_ETA_FALLBACK_LIMIT_MM', 0.06);
cfg.etaRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_ETA_REG_WEIGHT_MV', 0.35);
cfg.staticRegWeightMv = parse_nonnegative_numeric_env_local('STEP07J_STATIC_REG_WEIGHT_MV', 0.75);
cfg.vpSeedObservationMode = 'gradient_displacement';
cfg.vpGradientMinRatio = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_MIN_RATIO', 0.10);
cfg.vpGradientReferenceQuantile = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_REFERENCE_QUANTILE', 95);
cfg.vpGradientWeightPower = parse_nonnegative_numeric_env_local('STEP07J_VP_GRADIENT_WEIGHT_POWER', 2);
cfg.finalObjective = 'plain';
cfg.continuousFrequencyRefine = packageCfg.frequency.refineEachCandidate;
cfg.frequencyRefineHalfWidthHz = packageCfg.frequency.refineHalfWidthHz;
cfg.activeEtaModel = 'theoretical_zero_eta';
cfg.sensorEtaMode = 'fixed_zero';
cfg.fitSensorEta = false;
cfg.noGapKernelMode = 'foundation_core';
cfg.gapUpdateMode = 'vp_topk_delayed_gap_forward';
cfg.lockGapModelsToFixedEO = parse_logical_env_local('STEP07J_LOCK_GAP_EO_TO_FIXED', false);
cfg.initializeGapModelsFromFixed = false;
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
assert(strcmpi(cfg.structureMode, 'single_sync') && ...
    strcmpi(cfg.candidateMethod, 'structured_voltage_vp') && ...
    strcmpi(cfg.eoCandidateMode, 'vp'), ...
    'Formal unified route requires single_sync structured_voltage_vp candidates.');
cfg.prevWindowCandidateMode = 'off';
cfg.prevWindowRefreshEvery = floor(parse_nonnegative_numeric_env_local( ...
    'STEP07J_PREV_WINDOW_REFRESH_EVERY', 5));
cfg.phaseSafeExpansion = false;
cfg.phaseSafeMarginMm = parse_nonnegative_numeric_env_local('STEP07J_PHASE_SAFE_MARGIN_MM', 0.03);
cfg.phaseSafeFallbackMode = lower(strtrim(getenv('STEP07J_PHASE_SAFE_FALLBACK_MODE')));
if isempty(cfg.phaseSafeFallbackMode)
    cfg.phaseSafeFallbackMode = 'linear_vp';
end
if ~ismember(cfg.phaseSafeFallbackMode, {'off', 'linear_vp'})
    error('STEP07J_PHASE_SAFE_FALLBACK_MODE must be "off" or "linear_vp".');
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
cfg.localBundleSource = lower(strtrim(getenv('STEP07J_LOCAL_BUNDLE_SOURCE')));
if isempty(cfg.localBundleSource)
    cfg.localBundleSource = 'foundation_step05_bundle';
end
if ~ismember(cfg.localBundleSource, {'foundation_step05_bundle', 'local_direct'})
    error('STEP07J_LOCAL_BUNDLE_SOURCE must be "foundation_step05_bundle" or "local_direct".');
end
cfg.overshootPenaltyWeight = parse_nonnegative_numeric_env_local('STEP07J_OVERSHOOT_PENALTY_WEIGHT', 100);
cfg.domainToleranceMm = 1e-8;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 450, 'MaxFunEvals', 1800, ...
    'TolX', 1e-5, 'TolFun', 1e-5);
cfg.syntheticBundleFile = strtrim(getenv('STEP07J_SYNTHETIC_BUNDLE_FILE'));
if ~isempty(cfg.syntheticBundleFile)
    if exist(cfg.syntheticBundleFile, 'file') ~= 2
        error('STEP07J_SYNTHETIC_BUNDLE_FILE does not exist: %s', cfg.syntheticBundleFile);
    end
    cfg.maxWindows = 1;
    cfg.phaseSafeExpansion = false;
    cfg.prevWindowCandidateMode = 'off';
end

sensorOverride = strtrim(getenv('STEP07J_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07J_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
staticEtaOverride = strtrim(getenv('STEP07J_STATIC_ETA_FILE'));
if ~isempty(staticEtaOverride)
    error('Formal V1 enforces eta=0; STEP07J_STATIC_ETA_FILE is not accepted.');
end
cfg.fixedSensorEtaMm = zeros(1, numel(cfg.analysisSensors));
cfg.fixedSensorEtaInfo = struct('file', 'built_in_theoretical_zero', ...
    'source', 'theoretical_zero');
resultSuffixOverride = strtrim(getenv('STEP07J_RESULT_SUFFIX'));
resultSuffix = sanitize_suffix_local(resultSuffixOverride);
if isempty(resultSuffix)
    timeTag = format_analysis_time_tag_local(C0.flowConfig.identification.analysisStartTimeSec);
    if strcmpi(cfg.runMode, 'comparison')
        resultSuffix = ['_', timeTag, '_method_compare'];
    else
        resultSuffix = sprintf('_%s_main_gaponly_W%dS%d_F%dto%dHz_DFpm%gHz', ...
            timeTag, C0.flowConfig.identification.windowBladePasses, ...
            C0.flowConfig.identification.slidingStepBladePasses, ...
            round(cfg.freqSearchHz(1)), round(cfg.freqSearchHz(2)), ...
            cfg.frequencyRefineHalfWidthHz);
    end
end

rotDir = packageCfg.paths.calibration;
rotDynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');
expectedAnalysisStartTime = C0.flowConfig.identification.analysisStartTimeSec;

correctedLibOverride = strtrim(getenv('STEP07J_CORRECTED_LIB_FILE'));
if ~isempty(correctedLibOverride)
    correctedLibFile = correctedLibOverride;
    correctedLibMode = 'override';
else
    correctedLibFile = fullfile(packageCfg.paths.calibrationInputs, ...
        sprintf('Step06I_OffsetTiltShared_GapLibrary_%s_%s.mat', C0.dataset, C0.caseTag));
    CorrectedGapLibrary = [];
    correctedLibMode = 'bundled_formal_input';
end
templateOverride = strtrim(getenv('STEP07J_TEMPLATE_FILE'));
if isempty(templateOverride)
    templateOverride = packageCfg.files.lowSpeedTemplate;
end
if strcmpi(templateOverride, packageCfg.files.lowSpeedTemplate) && ~isfile(templateOverride)
    Prepare_AdaptiveSG_Template_20250527(false);
end
if ~isempty(templateOverride)
    templateFile = templateOverride;
    templateSourceMode = 'override';
elseif strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    templateFile = '';
    templateSourceMode = 'foundation_step05_bundle';
else
    addpath(rotDir);
    templateCase = struct();
    templateCase.dataset = '20250527';
    templateCase.bladeId = cfg.targetBlade;
    templateCase.sensorIds = cfg.analysisSensors;
    templateCase.sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
    templateCase.caseTag = sprintf('B%d_%s', cfg.targetBlade, templateCase.sensorTag);
    [templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20250527(rotDir, templateCase, 'GradientXRange030_OPRCenterStd');
end
foundationStep05Override = strtrim(getenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE'));
if ~isempty(foundationStep05Override)
    if exist(foundationStep05Override, 'file') ~= 2
        error('STEP07J_FOUNDATION_STEP05_RESULT_FILE does not exist: %s', foundationStep05Override);
    end
    foundationStep05File = foundationStep05Override;
else
    foundationStep05File = packageCfg.files.formalFoundationResult;
end
dynamicOverride = strtrim(getenv('STEP07J_DYNAMIC_MAP_FILE'));
highMapOverride = strtrim(getenv('STEP07J_HIGHMAP_FILE'));
if ~isempty(highMapOverride)
    if exist(highMapOverride, 'file') ~= 2
        error('STEP07J_HIGHMAP_FILE does not exist: %s', highMapOverride);
    end
    dynamicFile = highMapOverride;
    dynamicSourceMode = 'step06_highmap_override';
elseif ~isempty(dynamicOverride)
    if exist(dynamicOverride, 'file') ~= 2
        error('STEP07J_DYNAMIC_MAP_FILE does not exist: %s', dynamicOverride);
    end
    dynamicFile = dynamicOverride;
    dynamicSourceMode = 'dynamic_map_override';
elseif strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    dynamicFile = '[foundation Step05 bundle windows]';
    dynamicSourceMode = 'foundation_step05_bundle';
else
    dynamicFile = find_step06_highmap_file_local(outDir, C0.dataset, C0.caseTag, expectedAnalysisStartTime);
    dynamicSourceMode = 'step06_highmap';
end
directMainFile = '';

if ~isfile(correctedLibFile)
    error('Run Step06I first. Missing file: %s', correctedLibFile);
end

if strcmpi(dynamicSourceMode, 'foundation_step05_bundle')
    FoundationStep05Source = load_foundation_step05_source_local(foundationStep05File, cfg);
    if isempty(FoundationStep05Source)
        error(['STEP07J_LOCAL_BUNDLE_SOURCE=foundation_step05_bundle requires a Step05 foundation result. ' ...
            'Run Build_Step07J_FoundationFullBundle_20250527 first or set STEP07J_FOUNDATION_STEP05_RESULT_FILE.']);
    end
    if ~isempty(templateOverride)
        [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
            rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag, templateOverride);
    elseif isfield(FoundationStep05Source, 'Template') && isstruct(FoundationStep05Source.Template)
        Template = filter_template_sensors_local(FoundationStep05Source.Template, cfg.analysisSensors, sensorTag);
        templateFile = '[foundation Step05 Result.Template]';
    else
        [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
            rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag);
    end
    DynamicMap = dynamic_map_from_foundation_source_local(FoundationStep05Source, cfg);
elseif strcmpi(dynamicSourceMode, 'step06_highmap') || strcmpi(dynamicSourceMode, 'step06_highmap_override')
    Sd = load(dynamicFile, 'highMap');
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        FoundationStep05Source = load_foundation_step05_source_local(foundationStep05File, cfg);
        if isempty(FoundationStep05Source)
            error(['STEP07J_LOCAL_BUNDLE_SOURCE=foundation_step05_bundle requires a Step05 foundation result. ' ...
                'Run Build_Step07J_FoundationFullBundle_20250527 first or set STEP07J_FOUNDATION_STEP05_RESULT_FILE.']);
        end
        if ~isempty(templateOverride)
            [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
                rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag, templateOverride);
        elseif isfield(FoundationStep05Source, 'Template') && isstruct(FoundationStep05Source.Template)
            Template = filter_template_sensors_local(FoundationStep05Source.Template, cfg.analysisSensors, sensorTag);
            templateFile = '[foundation Step05 Result.Template]';
        else
            [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
                rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag);
        end
        DynamicMap = build_dynamic_map_from_step06_highmap_local(Sd.highMap, Template, cfg, expectedAnalysisStartTime);
    else
        St = load(templateFile, 'Template');
        Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
        DynamicMap = build_dynamic_map_from_step06_highmap_local(Sd.highMap, Template, cfg, expectedAnalysisStartTime);
    end
else
    Sd = load(dynamicFile, 'DynamicMap');
    DynamicMap = assert_dynamic_map_time_matches_local(Sd.DynamicMap, expectedAnalysisStartTime, dynamicFile);
end
if ~exist('FoundationStep05Source', 'var')
    FoundationStep05Source = load_foundation_step05_source_local(foundationStep05File, cfg);
end
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    if isempty(FoundationStep05Source)
        error(['STEP07J_LOCAL_BUNDLE_SOURCE=foundation_step05_bundle requires a Step05 foundation result. ' ...
            'Run Build_Step07J_FoundationFullBundle_20250527 first or set STEP07J_FOUNDATION_STEP05_RESULT_FILE.']);
    end
    % The Foundation bundle supplies points and templates only. Formal V1
    % never imports its historical nonzero eta values.
    cfg.fixedSensorEtaMm = zeros(1, numel(cfg.analysisSensors));
    if ~isempty(templateOverride)
        [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
            rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag, templateOverride);
    elseif isfield(FoundationStep05Source, 'Template') && isstruct(FoundationStep05Source.Template)
        Template = filter_template_sensors_local(FoundationStep05Source.Template, cfg.analysisSensors, sensorTag);
        templateFile = '[foundation Step05 Result.Template]';
    else
        [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
            rootDir, cfg.targetBlade, cfg.analysisSensors, sensorTag);
    end
elseif ~exist('Template', 'var') || isempty(Template)
    St = load(templateFile, 'Template');
    Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
end
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    coordinateCheck = struct('table', table(), 'tolerance_mm', 1e-6, ...
        'max_abs_delta_mm', 0, 'is_consistent', true, ...
        'status', 'foundation_step05_template_bundle_coordinates');
else
    DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, cfg.analysisSensors, sensorTag);
    coordinateCheck = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysisSensors, 1e-6);
end
if ~exist('CorrectedGapLibrary', 'var') || isempty(CorrectedGapLibrary)
    Sc = load(correctedLibFile, 'CorrectedGapLibrary');
    CorrectedGapLibrary = Sc.CorrectedGapLibrary;
    correctedLibMode = 'combined';
end
responseSurface = CorrectedGapLibrary.responseSurface;

availableWindows = numel(DynamicMap.Window);
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    availableWindows = min(availableWindows, numel(FoundationStep05Source.Window));
end
numWindows = min(availableWindows, cfg.maxWindows);
SyntheticInput = struct();
if ~isempty(cfg.syntheticBundleFile)
    SyntheticInput = load(cfg.syntheticBundleFile, 'SyntheticBundle', 'SyntheticTruth');
    if ~isfield(SyntheticInput, 'SyntheticBundle') || ~isstruct(SyntheticInput.SyntheticBundle)
        error('Synthetic input must contain a SyntheticBundle struct: %s', cfg.syntheticBundleFile);
    end
    numWindows = 1;
end
fprintf('\n=== Step07J: nested static warp VP/full-wave ===\n');
fprintf('Template: %s [%s]\n', templateFile, templateSourceMode);
fprintf('Dynamic source: %s [%s]\n', dynamicFile, dynamicSourceMode);
fprintf('Corrected gap library: %s [%s]\n', correctedLibFile, correctedLibMode);
if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
    fprintf('Foundation Step05 bundle source: %s\n', foundationStep05File);
end
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinateCheck.max_abs_delta_mm, coordinateCheck.status);
legacyBaseXcGradientTpl = contains(templateFile, 'GradientXRange030') && ...
    ~contains(templateFile, 'OPRAnchored') && ~contains(dynamicFile, 'GradientXRange030');
if ~coordinateCheck.is_consistent && legacyBaseXcGradientTpl
    fprintf('Coordinate note: legacy 20250527 route uses base-xc DynamicMap with GradientXRange030 template.\n');
elseif ~coordinateCheck.is_consistent
    warning('Template/DynamicMap x-center mismatch. Do not use this run as the unified GradientXRange030 main result.');
end
fprintf('Windows: %d, sensors: %s, VP Top-K: %d, EO candidate mode: %s\n', ...
    numWindows, mat2str(cfg.analysisSensors), cfg.vpTopK, cfg.eoCandidateMode);
fprintf('Prev-window EO soft candidate: %s, refresh every %d window(s).\n', ...
    cfg.prevWindowCandidateMode, cfg.prevWindowRefreshEvery);
fprintf('Phase-safe expansion: %d, margin %.3f mm, fallback %s.\n', ...
    cfg.phaseSafeExpansion, cfg.phaseSafeMarginMm, cfg.phaseSafeFallbackMode);
fprintf('Local foundation-core bundle source: %s.\n', cfg.localBundleSource);
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
fprintf('VP seed: %s, gradient gate %.3f*q%.1f, weight power %.3f.\n', ...
    cfg.vpSeedObservationMode, cfg.vpGradientMinRatio, ...
    cfg.vpGradientReferenceQuantile, cfg.vpGradientWeightPower);
fprintf('Eta policy: theoretical zero from %s, eta=%s mm.\n', ...
    cfg.fixedSensorEtaInfo.file, mat2str(cfg.fixedSensorEtaMm, 6));

WindowResult = struct([]);
trendRows = cell(numWindows, 1);
bestIdx = 1;
bestRmse = inf;
prevWindowResult = [];
prevLocalDirectPhaseInfo = empty_phase_safe_reference_local();
dxReferenceTracker = init_causal_dx_reference_tracker_local(cfg.causalDxRef);

for iw = 1:numWindows
    Wmap = DynamicMap.Window(iw);
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        Wmap = foundation_wmap_for_window_local(FoundationStep05Source.Window(iw), Wmap);
        foundationWindow = FoundationStep05Source.Window(iw);
        coreBundleFull = build_gap_bundle_from_foundation_bundle_local( ...
            foundationWindow.CoreBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        coreBundleFull.localSource = 'foundation_step05_core_bundle';
        coreDirectBundle = [];
    else
        coreDirectBundle = build_local_direct_template_bundle_local(Wmap, Template, cfg, []);
        coreBundleFull = build_gap_bundle_from_direct_bundle_local( ...
            coreDirectBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
        coreBundleFull.localSource = 'local_direct_core';
    end
    if ~isempty(cfg.syntheticBundleFile)
        coreBundleFull = SyntheticInput.SyntheticBundle;
        coreBundleFull.localSource = 'synthetic_bundle_override';
        coreBundleFull.pointCount = numel(coreBundleFull.V);
        coreDirectBundle = [];
    end
    coreBundle = decimate_bundle_local(coreBundleFull, cfg.maxPointsPerWindow);
    eoCandidates = build_eo_candidates_local(coreBundle.rotFreqMeanHz, cfg.freqSearchHz, cfg.eoPad);
    coreSeedTable = solve_vp_seed_scan_local(coreBundle, eoCandidates, cfg);
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        localDirectCoreSeedTable = coreSeedTable;
    else
        localDirectCoreSeedTable = solve_local_direct_template_seed_scan_local(coreDirectBundle, Template, eoCandidates, cfg);
    end
    phaseInfo = prevLocalDirectPhaseInfo;
    if cfg.prevWindowRefreshEvery > 0 && mod(iw - 1, cfg.prevWindowRefreshEvery) == 0
        phaseInfo = empty_phase_safe_reference_local();
        phaseInfo.reason = 'scheduled_refresh';
    end
    if cfg.phaseSafeExpansion && ~phaseInfo.usable && strcmpi(cfg.phaseSafeFallbackMode, 'linear_vp')
        phaseInfo = build_phase_safe_reference_local([], localDirectCoreSeedTable, iw, cfg, phaseInfo.reason);
    end

    bundleFull = coreBundleFull;
    usedPhaseSafeExpansion = false;
    expandedBundleFull = [];
    selectedDirectBundle = coreDirectBundle;
    if cfg.phaseSafeExpansion && phaseInfo.usable
        if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
            expandedBundleFull = build_gap_bundle_from_foundation_bundle_local( ...
                foundationWindow.FinalBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
            expandedBundleFull.localSource = 'foundation_step05_final_bundle';
        else
            expandedDirectBundle = build_local_direct_template_bundle_local(Wmap, Template, cfg, phaseInfo);
            expandedBundleFull = build_gap_bundle_from_direct_bundle_local( ...
                expandedDirectBundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
            selectedDirectBundle = expandedDirectBundle;
        end
        if expandedBundleFull.pointCount > coreBundleFull.pointCount
            bundleFull = expandedBundleFull;
            usedPhaseSafeExpansion = true;
        end
    end
    bundle = decimate_bundle_local(bundleFull, cfg.maxPointsPerWindow);
    seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg);
    [dxReferenceTracker, dxReferenceStep] = update_causal_dx_reference_tracker_local( ...
        dxReferenceTracker, seedTable, bundle, cfg.causalDxRef);
    prevEoInfo = build_prev_window_eo_reference_local(prevWindowResult, iw, cfg);
    [selectedEO, selectionInfo] = select_eo_candidates_local( ...
        seedTable, cfg.vpTopK, [], iw, cfg.eoCandidateMode, prevEoInfo);
    [coreSelectedEO, ~] = select_eo_candidates_local( ...
        coreSeedTable, cfg.vpTopK, [], iw, cfg.eoCandidateMode, prevEoInfo);
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
    if ~isempty(selectedDirectBundle) && isfield(selectedDirectBundle, 'selection_pass')
        selectionInfo.localDirectSelectionPass = selectedDirectBundle.selection_pass;
    else
        selectionInfo.localDirectSelectionPass = 'foundation_step05_bundle';
    end
    selectionInfo.CausalDxReference = dxReferenceStep;
    mainCandidateEO = selectedEO;

    cfgWindow = cfg;
    cfgWindow.CausalDxReference = dxReferenceStep;
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
    cfgFixed = configure_dx_reference_local(cfgWindow, seedTable, selectedEO, 'fixed');
    modelFits.fixed = [];
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle') && isempty(cfg.syntheticBundleFile)
        modelFits.fixed = fixed_fit_from_foundation_step05_local( ...
            FoundationStep05Source.Window(iw), bundle, cfgFixed);
    end
    if isempty(modelFits.fixed)
        modelFits.fixed = refine_static_warp_fit_local(bundle, seedTable, selectedEO, cfgFixed, 'fixed');
    end
    fitModelNames = setdiff(cfg.modelNames, {'fixed'}, 'stable');
    for im = 1:numel(fitModelNames)
        modeName = fitModelNames{im};
        fitEO = selectedEO;
        if cfg.lockGapModelsToFixedEO && mode_has_gap_local(modeName) && isfinite(modelFits.fixed.EO)
            fitEO = modelFits.fixed.EO;
        end
        cfgFit = configure_dx_reference_local(cfgWindow, seedTable, fitEO, modeName);
        modelFits.(modeName) = refine_static_warp_fit_local( ...
            bundle, seedTable, fitEO, cfgFit, modeName, modelFits.fixed);
    end
    wr = pack_window_result_local(iw, Wmap, bundleFull, seedTable, ...
        mainCandidateEO, selectionInfo, modelFits);
    if iw == 1
        WindowResult = repmat(wr, numWindows, 1);
    else
        WindowResult(iw) = wr;
    end
    trendRows{iw} = make_trend_row_local(wr, [], cfg);
    prevWindowResult = wr;
    if strcmpi(cfg.localBundleSource, 'foundation_step05_bundle')
        prevLocalDirectPhaseInfo = build_phase_safe_reference_local( ...
            [], seedTable, iw, cfg, 'foundation_step05_selected');
    else
        localDirectSelectedSeedTable = solve_local_direct_template_seed_scan_local(selectedDirectBundle, Template, eoCandidates, cfg);
        prevLocalDirectPhaseInfo = build_phase_safe_reference_local([], localDirectSelectedSeedTable, iw, cfg, 'local_direct_selected');
    end
    if modelFits.(cfg.mainModel).plainRmseMv < bestRmse
        bestRmse = modelFits.(cfg.mainModel).plainRmseMv;
        bestIdx = iw;
    end
    fprintf('Window %02d/%02d laps %s: %s %.2f mV\n', ...
        iw, numWindows, mat2str(Wmap.lap_range), cfg.mainModel, ...
        modelFits.(cfg.mainModel).plainRmseMv);
end

Trend = vertcat(trendRows{:});
Summary = build_summary_table_local(Trend, bestIdx, cfg);

Result = struct();
Result.dataset = '20250527';
if strcmpi(cfg.runMode, 'comparison')
    Result.method = 'gap_tilt_comparison_with_nested_static_warp_ablation_vp_full_wave';
else
    Result.method = 'unified_structured_voltage_single_sync_full_wave';
end
Result.description = ['Formal V1 uses theoretical eta=0, gradient-displacement VP Top-K seeds, independent windows, per-candidate continuous-frequency refinement, and ordinary voltage RMSE ranking. ' ...
    'The gap-aware models are nested updates of this fixed kernel: gap_only and gap_tilt keep the fixed solution as a nested zero-update baseline but may reselect EO inside the same full waveform objective. ' ...
    'The paper-facing gap_only/gap_fixed_tilt model identifies only the per-sensor gap increment dg; the low-speed calibrated sensor slope mu remains fixed and no residual high-speed dmu is identified. gap_only with zero dg reduces to the fixed no-gap kernel. ' ...
    'The gap-increment search half-width is computed window-by-window from Delta g_proj = epsilon_perp / S_g and Delta g_loc = alpha Delta g_proj, with only explicit numerical floor/cap/fallback guards. ' ...
    'The default main run computes only gap_only. fixed and gap_tilt are computed only in STEP07J_RUN_MODE=comparison. gap_tilt is a sensitivity ablation because dmu can compete with vibration amplitude; gap_tilt_shift remains excluded because the extra sensor-shift freedom is invalid for the paper-facing method.'];
Result.cfg = cfg;
Result.MethodContract = struct( ...
    'schema_version', 'method_contract_v2', ...
    'stage', 'gap_aware', ...
    'structure_mode', cfg.structureMode, ...
    'structure_selection', packageCfg.frequency.structureSelection, ...
    'unknown_structure_criterion', packageCfg.frequency.unknownStructureCriterion, ...
    'candidate_method', cfg.candidateMethod, ...
    'eo_top_k', 3, ...
    'strict_top_k', struct('top_k_min', 3, 'top_k_max', 3, ...
        'top_rank_keep', 3, 'gap_ratio_keep', 0, 'include_neighbor_eo', false), ...
    'frequency_search_hz', [300 1000], ...
    'continuous_frequency_refinement', true, ...
    'frequency_refinement_half_width_hz', 2, ...
    'final_objective', 'plain_unweighted_voltage_rmse', ...
    'low_speed_template_method', packageCfg.lowSpeed.templateMethod, ...
    'template_interpolation', packageCfg.lowSpeed.interpolationMethod, ...
    'forward_model', packageCfg.model.forwardModel, ...
    'eta_policy', 'zero', ...
    'model', 'gap_only', ...
    'bundle_policy', 'core_only_bundle', ...
    'phase_safe_expansion', false, ...
    'previous_window_information_used', false);
Result.templateFile = templateFile;
Result.dynamicFile = dynamicFile;
Result.correctedLibFile = correctedLibFile;
Result.correctedLibMode = correctedLibMode;
Result.gapCalibrationIndexFile = find_gap_calibration_index_file_local(outDir, C0);
Result.gapCalibrationSourceFile = correctedLibFile;
Result.directMainFile = directMainFile;
Result.foundationStep05File = foundationStep05File;
Result.FoundationStep05SourceMode = cfg.localBundleSource;
Result.CoordinateCheck = coordinateCheck;
Result.WindowResult = WindowResult;
Result.Trend = Trend;
Result.Summary = Summary;
Result.BestWindowIndex = bestIdx;
Result.MainModel = cfg.mainModel;
Result.BestWindow = WindowResult(bestIdx);
if ~isempty(cfg.syntheticBundleFile)
    Result.syntheticBundleFile = cfg.syntheticBundleFile;
    if isfield(SyntheticInput, 'SyntheticTruth')
        Result.SyntheticTruth = SyntheticInput.SyntheticTruth;
    end
end

matFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s%s.mat', ...
    C0.dataset, C0.caseTag, resultSuffix));
if strcmpi(cfg.runMode, 'main') && isinf(cfg.maxWindows) && isempty(resultSuffixOverride)
    matFile = packageCfg.files.formalGapResult;
end
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
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_20250527.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20250527.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_BaseFrame_20250527.mat', targetBlade, sensorTag)
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

function dynamicFile = find_dynamic_map_file_local(dynamicDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows*OPRCenterStd*_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Diag80L_W3S1_BaseXc_GradientTpl_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Diag80L_W3S1_GradientXRange030_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Diag80L_W3S1_OPRAnchored_20250527.mat', targetBlade, sensorTag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20250527.mat', targetBlade, sensorTag)
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

function highMapFile = find_step06_highmap_file_local(outDir, dataset, caseTag, analysisStartTime)
timeTag = format_analysis_time_tag_local(analysisStartTime);
candidate = fullfile(outDir, sprintf('Step06_HighMap_%s_%s_%s.mat', dataset, caseTag, timeTag));
if isfile(candidate)
    highMapFile = candidate;
    return;
end
legacyCandidate = fullfile(outDir, sprintf('Step06_HighMap_%s_%s.mat', dataset, caseTag));
if isfile(legacyCandidate)
    loaded = load(legacyCandidate, 'highMap');
    if isfield(loaded, 'highMap') && highmap_time_matches_local(loaded.highMap, analysisStartTime)
        highMapFile = legacyCandidate;
        return;
    end
end
error(['No Step06 highMap matching analysisStartTime %.3f s was found. ' ...
    'Run Step06_Build_HighMap_20250527.m with BLADE_ANALYSIS_START_TIME_SEC set to the same value.'], ...
    analysisStartTime);
end

function DynamicMap = assert_dynamic_map_time_matches_local(DynamicMap, expectedStartTime, dynamicFile)
actualStartTime = NaN;
if isfield(DynamicMap, 'SourceSettings') && isfield(DynamicMap.SourceSettings, 'analysis_start_time')
    actualStartTime = DynamicMap.SourceSettings.analysis_start_time;
elseif isfield(DynamicMap, 'analysisStartTime')
    actualStartTime = DynamicMap.analysisStartTime;
end
if isfinite(expectedStartTime) && isfinite(actualStartTime) && abs(actualStartTime - expectedStartTime) > 1e-9
    error(['DynamicMap analysis_start_time %.3f s does not match current config %.3f s. ' ...
        'File: %s'], actualStartTime, expectedStartTime, dynamicFile);
end
end

function tf = highmap_time_matches_local(highMap, expectedStartTime)
tf = isfield(highMap, 'analysisStartTime') && ...
    isfinite(highMap.analysisStartTime) && ...
    abs(highMap.analysisStartTime - expectedStartTime) <= 1e-9;
end

function DynamicMap = build_dynamic_map_from_step06_highmap_local(highMap, Template, cfg, expectedStartTime)
if ~highmap_time_matches_local(highMap, expectedStartTime)
    error('Step06 highMap analysisStartTime does not match current config %.3f s.', expectedStartTime);
end
if ~isfield(highMap, 'targetLaps') || ~isfield(highMap, 'analysisWinSize') || ~isfield(highMap, 'slidingStep')
    error('Step06 highMap is missing sliding-window settings. Re-run Step06.');
end
if highMap.targetLaps < highMap.analysisWinSize
    error('Step06 highMap targetLaps (%d) is smaller than analysisWinSize (%d).', ...
        highMap.targetLaps, highMap.analysisWinSize);
end

numWindows = floor((highMap.targetLaps - highMap.analysisWinSize) / highMap.slidingStep) + 1;
Window = repmat(struct( ...
    'window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, ...
    'Sensor', []), numWindows, 1);

for iw = 1:numWindows
    lapStart = 1 + (iw - 1) * highMap.slidingStep;
    lapRange = lapStart:(lapStart + highMap.analysisWinSize - 1);
    Sensor = repmat(struct( ...
        'sensor_id', NaN, 't', [], 'x_abs', [], 'x_rel', [], ...
        'V', [], 'W', [], 'theta', [], 'point_count', NaN, ...
        'x_center_reference_mm', NaN), numel(cfg.analysisSensors), 1);
    allT = [];
    for is = 1:numel(cfg.analysisSensors)
        sid = cfg.analysisSensors(is);
        mask = highMap.S_v == sid & ismember(highMap.rev_v, lapRange);
        Tpl = get_template_sensor_local(Template, sid);
        Sensor(is).sensor_id = sid;
        Sensor(is).t = highMap.t_v(mask);
        Sensor(is).x_abs = highMap.x_mm_v(mask);
        Sensor(is).x_rel = highMap.x_v(mask);
        Sensor(is).V = highMap.V_raw_mV(mask) ./ 1000;
        Sensor(is).W = highMap.W_v(mask);
        Sensor(is).theta = highMap.theta_v(mask);
        Sensor(is).point_count = nnz(mask);
        Sensor(is).x_center_reference_mm = Tpl.xc;
        allT = [allT; Sensor(is).t(:)]; %#ok<AGROW>
    end
    Window(iw).window_id = iw;
    Window(iw).lap_range = lapRange;
    Window(iw).time_window = [min(allT), max(allT)];
    Window(iw).rot_freq_mean_hz = highMap.rotFreqHz;
    Window(iw).rot_rpm_mean = highMap.rotRpm;
    Window(iw).Sensor = Sensor;
end

DynamicMap = struct();
DynamicMap.Dataset = highMap.dataset;
DynamicMap.CaseName = highMap.caseName;
DynamicMap.SensorIDs = cfg.analysisSensors;
DynamicMap.SensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
DynamicMap.SourceSettings = struct( ...
    'analysis_start_time', highMap.analysisStartTime, ...
    'target_laps', highMap.targetLaps, ...
    'analysis_win_size', highMap.analysisWinSize, ...
    'sliding_step', highMap.slidingStep);
DynamicMap.TemplateFileForXRel = highMap.lowSpeedTemplateFile;
DynamicMap.Window = Window;
end

function correctedLibFile = find_corrected_gap_library_file_local(outDir, targetBlade, sensorTag)
patterns = {
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s_OPRCenterStd.mat', targetBlade, sensorTag)
    };
for ip = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        correctedLibFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
correctedLibFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s.mat', ...
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

function [Template, templateFile] = load_foundation_step04_template_for_blade_local( ...
    rootDir, targetBlade, sensorIds, sensorTag, explicitTemplateFile)
if nargin < 5
    explicitTemplateFile = '';
end
packageCfg = Config_20250527();
templateDir = fullfile(packageCfg.paths.preparedInputs, 'foundation', ...
    'step04_low_speed_template');
preferred = fullfile(templateDir, ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S136_20250527.mat');
if ~isempty(explicitTemplateFile)
    templateFile = explicitTemplateFile;
elseif isfile(preferred)
    templateFile = preferred;
else
    files = dir(fullfile(templateDir, ...
        'Template_OPRCenterStd_LowSpeed_AllBlades_*_20250527.mat'));
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
    if isfield(hit, 'threshold') && isfinite(hit.threshold)
        Sensor(is).threshold = hit.threshold;
    end
    Sensor(is).xc = hit.xc;
    if isfield(hit, 'baseline') && isfinite(hit.baseline)
        Sensor(is).baseline = hit.baseline;
    end
end
Template = struct();
Template.Dataset = '20250527';
Template.SourceFile = templateFile;
Template.SourceMode = 'foundation_step04_sensorblade';
Template.Sensor = Sensor;
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = sensorTag;
Template.TargetBlade = targetBlade;
end

function resultFile = find_foundation_step05_result_file_local(rootDir, targetBlade, sensorIds)
packageCfg = Config_20250527();
resultFile = fullfile(packageCfg.paths.preparedInputs, 'foundation', ...
    sprintf('Result_Step05_SingleSyncDirectTemplate_B%d_S%s_20250527.mat', ...
    targetBlade, sprintf('%d', sensorIds)));
if isfile(resultFile)
    return;
end
resultFile = '';
sensorTag = ['S', sprintf('%d', sensorIds)];
foundationOutputRoot = fullfile(packageCfg.paths.preparedInputs, 'foundation');
patterns = {
    fullfile(foundationOutputRoot, 'step05_single_sync_direct_template_fullbundle_gap_step07j', '*', sprintf( ...
        'Result_Step05_SingleSyncDirectTemplate_B%d_%s_20250527.mat', targetBlade, sensorTag))
    fullfile(foundationOutputRoot, 'step05_single_sync_direct_template_fullbundle_for_gap_step07j', '*', sprintf( ...
        'Result_Step05_SingleSyncDirectTemplate_B%d_%s_20250527.mat', targetBlade, sensorTag))
    fullfile(foundationOutputRoot, 'step05_single_sync_direct_template*', '*', sprintf( ...
        'Result_Step05_SingleSyncDirectTemplate_B%d_%s_20250527.mat', targetBlade, sensorTag))
    fullfile(foundationOutputRoot, 'step05_single_sync_direct_template*', '*', sprintf( ...
        'Result_Step05_*FixedJoint*B%d_%s_20250527.mat', targetBlade, sensorTag))
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
Wsrc = S.Result.WindowResult;
Window = repmat(struct('FinalBundle', [], 'CoreBundle', [], ...
    'FoundationResult', [], 'lap_range', [NaN NaN], ...
    'time_window_s', [NaN NaN]), numel(Wsrc), 1);
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
if isfield(S.Result, 'Template')
    Source.Template = S.Result.Template;
end
if isfield(S.Result, 'Trend')
    Source.Trend = S.Result.Trend;
else
    Source.Trend = table();
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
    if (isnumeric(v) || islogical(v)) && numel(v) == numel(keep)
        B.(f) = v(keep);
    end
end
B.point_count = nnz(keep);
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
        eta = eta(1:numel(sensorIds));
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
end
end

function DynamicMap = dynamic_map_from_foundation_source_local(Source, cfg)
Window = repmat(struct('window_id', NaN, 'lap_range', [], ...
    'time_window', [], 'time_window_s', [], 'rot_freq_mean_hz', NaN, ...
    'rot_rpm_mean', NaN, 'Sensor', []), numel(Source.Window), 1);
for iw = 1:numel(Source.Window)
    Fw = Source.Window(iw);
    B = Fw.FinalBundle;
    if isempty(B)
        B = Fw.CoreBundle;
    end
    Window(iw).window_id = iw;
    Window(iw).lap_range = Fw.lap_range;
    Window(iw).time_window_s = Fw.time_window_s;
    Window(iw).time_window = Fw.time_window_s;
    if isstruct(B) && isfield(B, 'rot_freq_mean_hz')
        Window(iw).rot_freq_mean_hz = B.rot_freq_mean_hz;
        Window(iw).rot_rpm_mean = 60 .* B.rot_freq_mean_hz;
    end
    if isstruct(B) && isfield(B, 'rot_rpm_mean') && isfinite(B.rot_rpm_mean)
        Window(iw).rot_rpm_mean = B.rot_rpm_mean;
    end
end
DynamicMap = struct();
DynamicMap.Dataset = '20250527';
DynamicMap.SensorIDs = cfg.analysisSensors;
DynamicMap.SensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
DynamicMap.SourceSettings = struct('analysis_start_time', NaN);
DynamicMap.Window = Window;
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
F0v = nan(size(bundle.X));
Fxv = nan(size(bundle.X));
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    mask = bundle.sensor_index == is;
    Tpl = get_template_sensor_local(Template, sid);
    F0v(mask) = interp1(Tpl.x_grid(:), Tpl.v_grid(:), bundle.X(mask), 'pchip', NaN);
    Fxv(mask) = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), bundle.X(mask), 'pchip', NaN);
end
seedBundle = bundle;
seedBundle.F0 = F0v;
seedBundle.Fx = Fxv;
seedBundle.sensorIds = bundle.sensor_ids;
seedBundle.sensorIndex = bundle.sensor_index;
seedBundle.pointCount = bundle.point_count;
fullEval = @(theta, eo, b, c) evaluate_local_direct_template_waveform_local(theta, eo, b, Template, c);
seedTable = step07jcore.solve_gradient_displacement_vp_seed(seedBundle, eoCandidates, cfg, fullEval);
end

function detail = evaluate_local_direct_template_waveform_local(theta, eo, bundle, Template, cfg)
theta = theta(:).';
A = min(abs(theta(1)), cfg.amplitudeLimitMm);
phi = wrap_to_pi_local(theta(2));
dx = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
eta = step07jcore.fixed_sensor_eta_by_position(bundle.sensor_ids, cfg);
u = A .* sin(vibration_phase_argument_local(eo, bundle, cfg) + phi);
vPred = nan(size(bundle.V));
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    mask = bundle.sensor_index == is;
    Tpl = get_template_sensor_local(Template, sid);
    xEval = bundle.X(mask) - dx - u(mask) - eta(is);
    xEval = min(max(xEval, min(Tpl.x_grid(:))), max(Tpl.x_grid(:)));
    vPred(mask) = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xEval, 'pchip', NaN);
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
if nnz(valid) < 8
    weightedRmse = inf;
    plainRmse = inf;
else
    res = bundle.V(valid) - vPred(valid);
    ww = max(bundle.W(valid), cfg.weightFloor);
    weightedRmse = 1000 * sqrt(sum(ww .* res.^2) / max(sum(ww), eps));
    plainRmse = 1000 * sqrt(mean(res.^2, 'omitnan'));
end
detail = struct('weightedRmseMv', weightedRmse, 'plainRmseMv', plainRmse);
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

function bundle = build_gap_bundle_from_foundation_bundle_local(Bsrc, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
if isempty(Bsrc) || ~isstruct(Bsrc) || ~isfield(Bsrc, 'x') || ~isfield(Bsrc, 'v')
    error('Foundation Step05 bundle is missing x/v observations.');
end
sensorIds = cfg.analysisSensors(:).';
[present, sensorIndex] = ismember(Bsrc.sensor_id(:), sensorIds);
X = Bsrc.x(:);
T = Bsrc.t(:);
Theta = Bsrc.theta(:);
Vraw = Bsrc.v(:);
if isfield(Bsrc, 'fit_weight') && numel(Bsrc.fit_weight) == numel(Vraw)
    W = Bsrc.fit_weight(:);
elseif isfield(Bsrc, 'legacy_fit_weight') && numel(Bsrc.legacy_fit_weight) == numel(Vraw)
    W = Bsrc.legacy_fit_weight(:);
elseif isfield(Bsrc, 'template_weight') && numel(Bsrc.template_weight) == numel(Vraw)
    W = Bsrc.template_weight(:);
else
    W = ones(size(Vraw));
end
V = nan(size(Vraw));
F0 = nan(size(Vraw));
Fx = nan(size(Vraw));
g0BySensor = NaN(numel(sensorIds), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = sensorIndex == is;
    Tpl = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    V(mask) = voltage_to_template_mv_relative_local(Vraw(mask), Tpl);
    % Always rebuild the candidate linearization from the active template
    % and the same corrected forward model used by the final fit. Cached
    % Foundation template_v/template_dv_dx may belong to an older template.
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
valid = present(:) & isfinite(X) & isfinite(T) & isfinite(Theta) & ...
    isfinite(V) & isfinite(W) & isfinite(F0) & isfinite(Fx) & isfinite(sensorIndex);
bundle = struct();
bundle.X = X(valid);
bundle.T = T(valid);
bundle.TRel = T(valid) - min(T(valid));
bundle.V = V(valid);
bundle.W = max(W(valid), cfg.weightFloor);
bundle.Theta = Theta(valid);
bundle.S = sensorIds(sensorIndex(valid)).';
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

function vMv = voltage_to_template_mv_relative_local(vRaw, Tpl)
vRaw = vRaw(:);
tplRef = max(abs(Tpl.v_grid(:)), [], 'omitnan');
rawRef = max(abs(vRaw(:)), [], 'omitnan');
if ~isfinite(tplRef), tplRef = 0; end
if ~isfinite(rawRef), rawRef = 0; end
if tplRef < 10 && rawRef < 10
    vMv = (vRaw - Tpl.baseline) .* 1000;
elseif tplRef < 10
    vMv = vRaw - Tpl.baseline .* 1000;
elseif rawRef < 10
    vMv = vRaw .* 1000 - Tpl.baseline;
else
    vMv = vRaw - Tpl.baseline;
end
end

function dvMv = voltage_gradient_to_mv_per_mm_local(dvRaw, Tpl)
dvRaw = dvRaw(:);
tplRef = max(abs(Tpl.v_grid(:)), [], 'omitnan');
if isfinite(tplRef) && tplRef < 10
    dvMv = dvRaw .* 1000;
else
    dvMv = dvRaw;
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

function fit = refine_static_warp_fit_local(bundle, seedTable, eoCandidates, cfg, modeName, fixedBaseFit)
if nargin < 6
    fixedBaseFit = [];
end
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'dxReferenceMm', NaN, 'deltaDxMm', NaN, ...
    'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'integerFrequencyHz', NaN, 'frequencyOffsetHz', 0, ...
    'integerPlainRmseMv', NaN, 'rmseDeltaMv', 0, ...
    'frequencyRefinementExitFlag', NaN, 'frequencyRefinementIterations', NaN, ...
    'residualObjectiveMv2', inf, 'overshootPenalty', inf, ...
    'staticRegPenalty', inf, 'etaRegPenalty', inf, 'invalidPenalty', inf, 'penaltyObjectiveMv2', inf, ...
    'penaltyShare', NaN, 'domainFeasible', false, 'domainMaxViolationMm', inf, ...
    'sensorEtaMm', [], 'deltaGapMm', [], 'deltaMuGapPerXMm', [], 'deltaTauMm', [], ...
    'VPred', [], 'uMm', [], 'score', inf), numel(eoCandidates), 1);
best = struct('score', inf);
cfgEval = attach_nested_fixed_base_local(cfg, modeName, fixedBaseFit);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    cfgCandidate = configure_gap_projection_candidate_local(cfgEval, eo, bundle.sensorIds);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    theta0 = initial_theta_local(seed.A, seed.phi, seed.dx, cfgCandidate, modeName, numel(bundle.sensorIds));
    theta0 = apply_fixed_kernel_initial_theta_local(theta0, fixedBaseFit, cfgCandidate, modeName, eo);
    theta0 = apply_delta_gap_projection_initial_theta_local(theta0, cfgCandidate, modeName, bundle.sensorIds, eo);
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfgCandidate, modeName);
    thetaOpt = fminsearch(fun, theta0, cfgCandidate.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfgCandidate, modeName, numel(bundle.sensorIds));
    detail = evaluate_static_warp_waveform_local(thetaOpt, eo, bundle, cfgCandidate, modeName);
    detail = refine_continuous_frequency_fit_local(bundle, detail, cfgCandidate, modeName);
    candidate(i).EO = eo;
    candidate(i).freqHz = detail.freqHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).dxReferenceMm = detail.dxReferenceMm;
    candidate(i).deltaDxMm = detail.deltaDxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    candidate(i).integerFrequencyHz = detail.integerFrequencyHz;
    candidate(i).frequencyOffsetHz = detail.frequencyOffsetHz;
    candidate(i).integerPlainRmseMv = detail.integerPlainRmseMv;
    candidate(i).rmseDeltaMv = detail.rmseDeltaMv;
    candidate(i).frequencyRefinementExitFlag = detail.frequencyRefinementExitFlag;
    candidate(i).frequencyRefinementIterations = detail.frequencyRefinementIterations;
    candidate(i).residualObjectiveMv2 = detail.residualObjectiveMv2;
    candidate(i).overshootPenalty = detail.overshootPenalty;
    candidate(i).staticRegPenalty = detail.staticRegPenalty;
    candidate(i).etaRegPenalty = detail.etaRegPenalty;
    candidate(i).invalidPenalty = detail.invalidPenalty;
    candidate(i).penaltyObjectiveMv2 = detail.penaltyObjectiveMv2;
    candidate(i).penaltyShare = detail.penaltyShare;
    candidate(i).domainFeasible = detail.domainFeasible;
    candidate(i).domainMaxViolationMm = detail.domainMaxViolationMm;
    candidate(i).sensorEtaMm = detail.sensorEtaMm;
    candidate(i).deltaGapMm = detail.deltaGapMm;
    candidate(i).deltaMuGapPerXMm = detail.deltaMuGapPerXMm;
    candidate(i).deltaTauMm = detail.deltaTauMm;
    candidate(i).VPred = detail.VPred;
    candidate(i).uMm = detail.uMm;
    candidate(i).score = detail.score;
    if detail.domainFeasible && detail.score < best.score
        best = candidate(i);
    end
end
if use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName) && ...
        all(isfinite([fixedBaseFit.EO, fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, fixedBaseFit.dxMm]))
    eo = round(fixedBaseFit.EO);
    thetaBase = initial_theta_local(fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, ...
        fixedBaseFit.dxMm, cfg, modeName, numel(bundle.sensorIds));
    thetaBase = bound_params_local(thetaBase, cfg, modeName, numel(bundle.sensorIds));
    detail = evaluate_static_warp_waveform_local(thetaBase, eo, bundle, cfgEval, modeName);
    baseCandidate = candidate(1);
    baseCandidate.EO = eo;
    baseCandidate.freqHz = eo * bundle.rotFreqMeanHz;
    baseCandidate.amplitudeMm = detail.amplitudeMm;
    baseCandidate.phaseRad = detail.phaseRad;
    baseCandidate.dxMm = detail.dxMm;
    baseCandidate.dxReferenceMm = detail.dxReferenceMm;
    baseCandidate.deltaDxMm = detail.deltaDxMm;
    baseCandidate.weightedRmseMv = detail.weightedRmseMv;
    baseCandidate.plainRmseMv = detail.plainRmseMv;
    baseCandidate.residualObjectiveMv2 = detail.residualObjectiveMv2;
    baseCandidate.overshootPenalty = detail.overshootPenalty;
    baseCandidate.staticRegPenalty = detail.staticRegPenalty;
    baseCandidate.etaRegPenalty = detail.etaRegPenalty;
    baseCandidate.invalidPenalty = detail.invalidPenalty;
    baseCandidate.penaltyObjectiveMv2 = detail.penaltyObjectiveMv2;
    baseCandidate.penaltyShare = detail.penaltyShare;
    baseCandidate.domainFeasible = detail.domainFeasible;
    baseCandidate.domainMaxViolationMm = detail.domainMaxViolationMm;
    baseCandidate.sensorEtaMm = detail.sensorEtaMm;
    baseCandidate.deltaGapMm = detail.deltaGapMm;
    baseCandidate.deltaMuGapPerXMm = detail.deltaMuGapPerXMm;
    baseCandidate.deltaTauMm = detail.deltaTauMm;
    baseCandidate.VPred = detail.VPred;
    baseCandidate.uMm = detail.uMm;
    baseCandidate.score = detail.score;
    baseCandidate.integerFrequencyHz = baseCandidate.freqHz;
    baseCandidate.frequencyOffsetHz = 0;
    baseCandidate.integerPlainRmseMv = detail.plainRmseMv;
    baseCandidate.rmseDeltaMv = 0;
    baseCandidate.deltaGapMm = zeros(numel(bundle.sensorIds), 1);
    baseCandidate.deltaMuGapPerXMm = zeros(numel(bundle.sensorIds), 1);
    baseCandidate.deltaTauMm = zeros(numel(bundle.sensorIds), 1);
    candidate(end+1) = baseCandidate; %#ok<AGROW>
    if baseCandidate.domainFeasible && baseCandidate.score <= best.score
        best = baseCandidate;
    end
end
if ~isfield(best, 'domainFeasible') || ~best.domainFeasible
    error('No physically feasible %s candidate was found for this window.', modeName);
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
candidateTable = struct2table(rmfield(candidate, ...
    {'VPred','uMm','sensorEtaMm','deltaGapMm','deltaMuGapPerXMm','deltaTauMm'}));
for is = 1:numel(bundle.sensorIds)
    values = NaN(numel(candidate), 1);
    for ic = 1:numel(candidate)
        if numel(candidate(ic).deltaGapMm) >= is
            values(ic) = candidate(ic).deltaGapMm(is);
        end
    end
    candidateTable.(sprintf('deltaGapS%dMm', bundle.sensorIds(is))) = values;
end
fit.CandidateTable = sortrows(candidateTable, ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
fit = enforce_nested_fixed_floor_local(fit, fixedBaseFit, cfg, modeName, bundle);
end

function fit = refine_continuous_frequency_fit_local(bundle, integerFit, cfg, modeName)
fit = integerFit;
fit.integerFrequencyHz = integerFit.freqHz;
fit.frequencyOffsetHz = 0;
fit.integerPlainRmseMv = integerFit.plainRmseMv;
fit.rmseDeltaMv = 0;
fit.frequencyRefinementExitFlag = 0;
fit.frequencyRefinementIterations = 0;
if ~cfg.continuousFrequencyRefine || cfg.frequencyRefineHalfWidthHz <= 0
    return;
end

nSensor = numel(bundle.sensorIds);
theta0 = initial_theta_local(integerFit.amplitudeMm, integerFit.phaseRad, ...
    integerFit.dxMm, cfg, modeName, nSensor);
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

z0 = [theta0, integerFit.freqHz];
fun = @(z) continuous_frequency_objective_local(z, integerFit.EO, bundle, cfg, modeName);
[zOpt, ~, exitFlag, output] = fminsearch(fun, z0, cfg.fminOptions);
thetaOpt = bound_params_local(zOpt(1:end-1), cfg, modeName, nSensor);
freqOpt = min(max(zOpt(end), max(cfg.freqSearchHz(1), ...
    integerFit.freqHz - cfg.frequencyRefineHalfWidthHz)), ...
    min(cfg.freqSearchHz(2), integerFit.freqHz + cfg.frequencyRefineHalfWidthHz));
cfgEval = cfg;
cfgEval.activeContinuousFrequencyHz = freqOpt;
detail = evaluate_static_warp_waveform_local(thetaOpt, integerFit.EO, bundle, cfgEval, modeName);
if ~detail.domainFeasible || ~isfinite(detail.score) || detail.score > integerFit.score * (1 + 1e-9)
    fit.frequencyRefinementExitFlag = -10;
    fit.frequencyRefinementIterations = output.iterations;
    return;
end

fit = detail;
fit.mode = modeName;
fit.integerFrequencyHz = integerFit.freqHz;
fit.frequencyOffsetHz = fit.freqHz - integerFit.freqHz;
fit.integerPlainRmseMv = integerFit.plainRmseMv;
fit.rmseDeltaMv = fit.plainRmseMv - integerFit.plainRmseMv;
fit.frequencyRefinementExitFlag = exitFlag;
fit.frequencyRefinementIterations = output.iterations;
end

function score = continuous_frequency_objective_local(z, eo, bundle, cfg, modeName)
nSensor = numel(bundle.sensorIds);
theta = bound_params_local(z(1:end-1), cfg, modeName, nSensor);
f0 = eo * bundle.rotFreqMeanHz;
freq = min(max(z(end), max(cfg.freqSearchHz(1), ...
    f0 - cfg.frequencyRefineHalfWidthHz)), ...
    min(cfg.freqSearchHz(2), f0 + cfg.frequencyRefineHalfWidthHz));
cfg.activeContinuousFrequencyHz = freq;
detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName);
boundPenalty = 1e6 * (sum((z(1:end-1) - theta).^2) + (z(end) - freq)^2);
score = detail.score + boundPenalty;
end

function phaseArg = vibration_phase_argument_local(eo, bundle, cfg)
phaseArg = eo .* bundle.Theta;
if ~isfield(cfg, 'activeContinuousFrequencyHz') || ...
        ~isfinite(cfg.activeContinuousFrequencyHz)
    return;
end
if ~isfield(bundle, 'T') || numel(bundle.T) ~= numel(bundle.Theta)
    error('Continuous-frequency refinement requires bundle.T aligned with bundle.Theta.');
end
t = bundle.T(:);
t0 = mean(t, 'omitnan');
detuningHz = cfg.activeContinuousFrequencyHz - eo * bundle.rotFreqMeanHz;
phaseArg = eo .* bundle.Theta + 2 * pi * detuningHz .* (t - t0);
end

function freqHz = vibration_frequency_hz_local(eo, bundle, cfg)
freqHz = eo * bundle.rotFreqMeanHz;
if isfield(cfg, 'activeContinuousFrequencyHz') && ...
        isfinite(cfg.activeContinuousFrequencyHz)
    freqHz = cfg.activeContinuousFrequencyHz;
end
end

function cfgOut = attach_nested_fixed_base_local(cfg, modeName, fixedBaseFit)
cfgOut = cfg;
if use_fixed_kernel_base_local(fixedBaseFit, cfg, modeName)
    cfgOut.NestedFixedBaseFit = fixedBaseFit;
end
end

function fit = enforce_nested_fixed_floor_local(fit, fixedBaseFit, cfg, modeName, bundle)
if strcmpi(cfg.finalObjective, 'plain')
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
if use_dx_reference_local(cfg, modeName)
    dxParam = fixedBaseFit.dxMm - cfg.dxReferenceMm;
else
    dxParam = fixedBaseFit.dxMm;
end
theta0(1:3) = [fixedBaseFit.amplitudeMm, fixedBaseFit.phaseRad, dxParam];
end

function cfgOut = configure_gap_projection_candidate_local(cfg, eo, sensorIds)
cfgOut = cfg;
nSensor = numel(sensorIds);
cfgOut.deltaGapCenterMm = zeros(1, nSensor);
limit = expand_limit_vector_local(cfg.deltaGapPriorScaleMm, nSensor);
cfgOut.deltaGapLowerMm = -limit;
cfgOut.deltaGapUpperMm = limit;
if ~isfield(cfg, 'DeltaGapProjection') || isempty(cfg.DeltaGapProjection) || ...
        ~isfield(cfg.DeltaGapProjection, 'AllTable') || isempty(cfg.DeltaGapProjection.AllTable)
    return;
end
T = cfg.DeltaGapProjection.AllTable;
required = {'sensorId','EO','deltaGapCenterMm','deltaGapLowerMm','deltaGapUpperMm'};
if ~all(ismember(required, T.Properties.VariableNames))
    return;
end
for is = 1:nSensor
    row = find(T.sensorId == sensorIds(is) & T.EO == round(eo), 1, 'first');
    if isempty(row)
        continue;
    end
    values = [T.deltaGapCenterMm(row), T.deltaGapLowerMm(row), T.deltaGapUpperMm(row)];
    if all(isfinite(values)) && values(2) <= values(1) && values(1) <= values(3)
        cfgOut.deltaGapCenterMm(is) = values(1);
        cfgOut.deltaGapLowerMm(is) = values(2);
        cfgOut.deltaGapUpperMm(is) = values(3);
    end
end
cfgOut.deltaGapLimitMm = max(abs([cfgOut.deltaGapLowerMm; cfgOut.deltaGapUpperMm]), [], 1);
end

function theta0 = apply_delta_gap_projection_initial_theta_local(theta0, cfg, modeName, sensorIds, ~)
nSensor = numel(sensorIds);
if ~mode_has_gap_local(modeName) || ~isfield(cfg, 'DeltaGapProjection') || ...
        isempty(cfg.DeltaGapProjection)
    return;
end
[idxDg, ~, ~] = mode_param_indices_local(modeName, nSensor);
if isempty(idxDg)
    return;
end
if isfield(cfg, 'deltaGapCenterMm')
    dgSeed = cfg.deltaGapCenterMm(:).';
elseif isfield(cfg.DeltaGapProjection, 'deltaGapLsSignedMm')
    dgSeed = cfg.DeltaGapProjection.deltaGapLsSignedMm(:).';
elseif isfield(cfg.DeltaGapProjection, 'deltaGapLsMm')
    dgSeed = cfg.DeltaGapProjection.deltaGapLsMm(:).';
else
    return;
end
if numel(dgSeed) < nSensor
    dgSeed(end+1:nSensor) = 0;
end
dgSeed = dgSeed(1:nSensor);
dgSeed(~isfinite(dgSeed)) = 0;
if isfield(cfg, 'deltaGapLowerMm') && isfield(cfg, 'deltaGapUpperMm')
    lower = expand_limit_vector_local(cfg.deltaGapLowerMm, nSensor);
    upper = expand_limit_vector_local(cfg.deltaGapUpperMm, nSensor);
    theta0(idxDg) = max(min(dgSeed, upper), lower);
else
    deltaGapLimit = expand_limit_vector_local(cfg.deltaGapLimitMm, nSensor);
    theta0(idxDg) = max(min(dgSeed, deltaGapLimit), -deltaGapLimit);
end
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
    eta = eta(1:nSensor);
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
theta = initial_theta_local(D.A_id, D.phi_id_wrapped, D.dx_c_id, ...
    cfg, 'fixed', nSensor);
detail = evaluate_static_warp_waveform_local(theta, round(D.EO_id), bundle, cfg, 'fixed');
if isfield(detail, 'weightNormalizer') && isfinite(detail.weightNormalizer) && detail.weightNormalizer > 0
    weightNormalizer = detail.weightNormalizer;
else
    weightNormalizer = max(numel(bundle.V), 1);
end
fit = detail;
fit.EO = round(D.EO_id);
fit.freqHz = D.fn_id;
fit.amplitudeMm = D.A_id;
fit.phaseRad = D.phi_id_wrapped;
fit.dxMm = D.dx_c_id;
fit.dxReferenceMm = 0;
fit.deltaDxMm = D.dx_c_id;
fit.sensorEtaMm = eta(:);
fit.deltaGapMm = zeros(nSensor, 1);
fit.deltaMuGapPerXMm = zeros(nSensor, 1);
fit.deltaTauMm = zeros(nSensor, 1);
fit.weightedRmseMv = rmseMv;
fit.plainRmseMv = plainRmseMv;
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

function proj = estimate_delta_gap_projection_limits_local(bundle, seedTable, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
nSensor = numel(bundle.sensorIds);
rowTemplate = struct('sensorId', NaN, 'EO', NaN, 'epsilonPerpMv', NaN, ...
    'SgMvPerMm', NaN, 'deltaGapProjMm', NaN, 'deltaGapLsMm', NaN, ...
    'deltaGapLsSignedMm', NaN, 'deltaGapLocMm', NaN, ...
    'deltaGapLimitMm', NaN, 'usedFallback', true, ...
    'clippedMin', false, 'clippedMax', false, ...
    'deltaGapCenterMm', NaN, 'deltaGapLowerMm', NaN, ...
    'deltaGapUpperMm', NaN, 'etaGap', NaN, ...
    'nuisanceRank', NaN, 'projectedGapRank', NaN, ...
    'nuisanceCondition', NaN, 'projectedGapCondition', NaN);
allRows = repmat(rowTemplate, 0, 1);
for ie = 1:numel(eoCandidates)
    eo = eoCandidates(ie);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    candidateRows = estimate_delta_gap_projection_candidate_local( ...
        bundle, seed.A, seed.phi, seed.dx, eo, cfg, rowTemplate);
    allRows = [allRows; candidateRows(:)]; %#ok<AGROW>
end
rows = repmat(rowTemplate, nSensor, 1);
for is = 1:nSensor
    sensorRows = allRows([allRows.sensorId] == bundle.sensorIds(is));
    if isempty(sensorRows)
        rows(is) = rowTemplate;
        rows(is).sensorId = bundle.sensorIds(is);
        continue;
    end
    [~, idx] = max([sensorRows.deltaGapLimitMm]);
    rows(is) = sensorRows(idx);
end
proj = struct();
proj.alpha = cfg.deltaGapProjAlpha;
proj.sensorId = [rows.sensorId].';
proj.EO = [rows.EO].';
proj.epsilonPerpMv = [rows.epsilonPerpMv].';
proj.SgMvPerMm = [rows.SgMvPerMm].';
proj.deltaGapProjMm = [rows.deltaGapProjMm].';
proj.deltaGapLsMm = [rows.deltaGapLsMm].';
proj.deltaGapLsSignedMm = [rows.deltaGapLsSignedMm].';
proj.deltaGapLocMm = [rows.deltaGapLocMm].';
proj.deltaGapLimitMm = [rows.deltaGapLimitMm].';
proj.deltaGapCenterMm = [rows.deltaGapCenterMm].';
proj.deltaGapLowerMm = [rows.deltaGapLowerMm].';
proj.deltaGapUpperMm = [rows.deltaGapUpperMm].';
proj.etaGap = [rows.etaGap].';
proj.nuisanceRank = [rows.nuisanceRank].';
proj.projectedGapRank = [rows.projectedGapRank].';
proj.nuisanceCondition = [rows.nuisanceCondition].';
proj.projectedGapCondition = [rows.projectedGapCondition].';
proj.usedFallback = [rows.usedFallback].';
proj.clippedMin = [rows.clippedMin].';
proj.clippedMax = [rows.clippedMax].';
proj.Table = struct2table(rows);
proj.AllTable = struct2table(allRows);
end

function rows = estimate_delta_gap_projection_candidate_local(bundle, A, phi, dx, eo, cfg, rowTemplate)
nSensor = numel(bundle.sensorIds);
nPoint = numel(bundle.V);
v0 = nan(nPoint, 1);
Jg = zeros(nPoint, nSensor);
Jnu = nan(nPoint, 4);
phase = eo .* bundle.Theta(:) + phi;
u = A .* sin(phase);
t = bundle.T(:);
tCenter = mean(t, 'omitnan');
eta = step07jcore.fixed_sensor_eta_by_position(bundle.sensorIds, cfg);
hG = max(cfg.deltaGapDerivativeStepMm, eps);
hX = max(cfg.derivativeStepMm, eps);
for is = 1:nSensor
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex == is;
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    xEval = bundle.X(mask) - dx - u(mask) - eta(is);
    [v0(mask), ~] = eval_static_warp_template_local( ...
        Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval);
    [vGp, ~] = eval_static_warp_template_local( ...
        Tpl, bundle.responseSurface, corr, hG, 0, 0, xEval);
    [vGm, ~] = eval_static_warp_template_local( ...
        Tpl, bundle.responseSurface, corr, -hG, 0, 0, xEval);
    [vXp, ~] = eval_static_warp_template_local( ...
        Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval + hX);
    [vXm, ~] = eval_static_warp_template_local( ...
        Tpl, bundle.responseSurface, corr, 0, 0, 0, xEval - hX);
    dMdx = (vXp - vXm) ./ (2 * hX);
    Jg(mask, is) = (vGp - vGm) ./ (2 * hG);
    Jnu(mask, 1) = -dMdx;
    Jnu(mask, 2) = -dMdx .* sin(eo .* bundle.Theta(mask));
    Jnu(mask, 3) = -dMdx .* cos(eo .* bundle.Theta(mask));
    Jnu(mask, 4) = -dMdx .* A .* cos(phase(mask)) .* 2 .* pi .* (t(mask) - tCenter);
end
r = bundle.V(:) - v0;
valid = isfinite(r) & isfinite(bundle.W(:)) & all(isfinite(Jg), 2) & all(isfinite(Jnu), 2);
rows = repmat(rowTemplate, nSensor, 1);
if nnz(valid) < max(8, size(Jnu, 2) + nSensor)
    for is = 1:nSensor
        rows(is) = fallback_projection_row_local(rowTemplate, bundle.sensorIds(is), eo, cfg);
    end
    return;
end
sw = sqrt(max(bundle.W(valid), cfg.weightFloor));
rw = sw .* r(valid);
Jgw = Jg(valid, :) .* sw;
Jnuw = Jnu(valid, :) .* sw;
[Unu, nuisanceRank, nuisanceCondition] = stable_column_basis_local(Jnuw);
rPerp = rw - Unu * (Unu.' * rw);
G = Jgw - Unu * (Unu.' * Jgw);
[~, projectedGapRank, projectedGapCondition] = stable_column_basis_local(G);
dgJoint = stable_least_squares_local(G, rPerp);
globalLimit = expand_limit_vector_local(cfg.deltaGapProjMaxLimitMm, nSensor);
for is = 1:nSensor
    others = setdiff(1:nSensor, is);
    [Uother, ~, ~] = stable_column_basis_local(G(:, others));
    q = G(:, is) - Uother * (Uother.' * G(:, is));
    rConditional = rPerp - Uother * (Uother.' * rPerp);
    qNorm = norm(q);
    rawNorm = norm(Jgw(:, is));
    etaGap = qNorm / max(rawNorm, eps);
    epsilonPerp = norm(rConditional);
    weak = ~isfinite(qNorm) || qNorm <= cfg.deltaGapProjSensitivityFloorMvPerMm || ...
        ~isfinite(etaGap) || etaGap < cfg.deltaGapProjEtaFloor || ...
        projectedGapRank < nSensor;
    if weak
        rows(is) = fallback_projection_row_local(rowTemplate, bundle.sensorIds(is), eo, cfg);
        rows(is).etaGap = etaGap;
    else
        signedLs = dgJoint(is);
        deltaProj = epsilonPerp / qNorm;
        halfWidthRaw = cfg.deltaGapProjAlpha * deltaProj;
        halfWidth = min(max(halfWidthRaw, cfg.deltaGapProjMinLimitMm), ...
            min(cfg.deltaGapProjMaxLimitMm, globalLimit(is)));
        center = min(max(signedLs, -globalLimit(is)), globalLimit(is));
        lower = max(-globalLimit(is), center - halfWidth);
        upper = min(globalLimit(is), center + halfWidth);
        effectiveHalfWidth = max(center - lower, upper - center);
        rows(is) = rowTemplate;
        rows(is).sensorId = bundle.sensorIds(is);
        rows(is).EO = eo;
        rows(is).epsilonPerpMv = epsilonPerp;
        rows(is).SgMvPerMm = qNorm;
        rows(is).deltaGapProjMm = deltaProj;
        rows(is).deltaGapLsMm = abs(signedLs);
        rows(is).deltaGapLsSignedMm = signedLs;
        rows(is).deltaGapLocMm = halfWidthRaw;
        rows(is).deltaGapLimitMm = effectiveHalfWidth;
        rows(is).usedFallback = false;
        rows(is).clippedMin = halfWidthRaw < cfg.deltaGapProjMinLimitMm;
        rows(is).clippedMax = halfWidthRaw > min(cfg.deltaGapProjMaxLimitMm, globalLimit(is));
        rows(is).deltaGapCenterMm = center;
        rows(is).deltaGapLowerMm = lower;
        rows(is).deltaGapUpperMm = upper;
        rows(is).etaGap = etaGap;
    end
    rows(is).nuisanceRank = nuisanceRank;
    rows(is).projectedGapRank = projectedGapRank;
    rows(is).nuisanceCondition = nuisanceCondition;
    rows(is).projectedGapCondition = projectedGapCondition;
end
end

function row = fallback_projection_row_local(rowTemplate, sensorId, eo, cfg)
row = rowTemplate;
limit = min(cfg.deltaGapProjFallbackLimitMm, cfg.deltaGapProjMaxLimitMm);
row.sensorId = sensorId;
row.EO = eo;
row.deltaGapProjMm = limit / max(cfg.deltaGapProjAlpha, eps);
row.deltaGapLocMm = limit;
row.deltaGapLimitMm = limit;
row.deltaGapLsMm = 0;
row.deltaGapLsSignedMm = 0;
row.deltaGapCenterMm = 0;
row.deltaGapLowerMm = -limit;
row.deltaGapUpperMm = limit;
row.usedFallback = true;
end

function [U, rankA, condA] = stable_column_basis_local(A)
if isempty(A)
    U = zeros(size(A, 1), 0);
    rankA = 0;
    condA = 1;
    return;
end
[Uall, S, ~] = svd(A, 'econ');
s = diag(S);
if isempty(s) || ~isfinite(s(1)) || s(1) <= 0
    U = zeros(size(A, 1), 0);
    rankA = 0;
    condA = inf;
    return;
end
tol = max(size(A)) * eps(s(1));
rankA = sum(s > tol);
U = Uall(:, 1:rankA);
if rankA > 0
    condA = s(1) / s(rankA);
else
    condA = inf;
end
end

function x = stable_least_squares_local(A, b)
if isempty(A)
    x = zeros(size(A, 2), 1);
    return;
end
[U, S, V] = svd(A, 'econ');
s = diag(S);
if isempty(s) || ~isfinite(s(1)) || s(1) <= 0
    x = zeros(size(A, 2), 1);
    return;
end
tol = max(size(A)) * eps(s(1));
keep = s > tol;
if ~any(keep)
    x = zeros(size(A, 2), 1);
    return;
end
x = V(:, keep) * ((U(:, keep).' * b) ./ s(keep));
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
proj.deltaGapLsSignedMm = NaN(nSensor, 1);
proj.deltaGapLocMm = NaN(nSensor, 1);
proj.deltaGapLimitMm = repmat(cfg.deltaGapLimitMm, nSensor, 1);
proj.deltaGapCenterMm = zeros(nSensor, 1);
proj.deltaGapLowerMm = -proj.deltaGapLimitMm;
proj.deltaGapUpperMm = proj.deltaGapLimitMm;
proj.etaGap = NaN(nSensor, 1);
proj.nuisanceRank = NaN(nSensor, 1);
proj.projectedGapRank = NaN(nSensor, 1);
proj.nuisanceCondition = NaN(nSensor, 1);
proj.projectedGapCondition = NaN(nSensor, 1);
proj.usedFallback = true(nSensor, 1);
proj.clippedMin = false(nSensor, 1);
proj.clippedMax = false(nSensor, 1);
proj.Table = table(proj.sensorId, proj.EO, proj.epsilonPerpMv, proj.SgMvPerMm, ...
    proj.deltaGapProjMm, proj.deltaGapLsMm, proj.deltaGapLsSignedMm, ...
    proj.deltaGapLocMm, proj.deltaGapLimitMm, proj.deltaGapCenterMm, ...
    proj.deltaGapLowerMm, proj.deltaGapUpperMm, proj.etaGap, ...
    proj.nuisanceRank, proj.projectedGapRank, ...
    proj.nuisanceCondition, proj.projectedGapCondition, ...
    proj.usedFallback, proj.clippedMin, proj.clippedMax, ...
    'VariableNames', {'sensorId','EO','epsilonPerpMv','SgMvPerMm','deltaGapProjMm', ...
    'deltaGapLsMm','deltaGapLsSignedMm','deltaGapLocMm','deltaGapLimitMm', ...
    'deltaGapCenterMm','deltaGapLowerMm','deltaGapUpperMm','etaGap', ...
    'nuisanceRank','projectedGapRank','nuisanceCondition','projectedGapCondition', ...
    'usedFallback','clippedMin','clippedMax'});
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
    if isfield(cfg, 'deltaGapLowerMm') && isfield(cfg, 'deltaGapUpperMm')
        lower = expand_limit_vector_local(cfg.deltaGapLowerMm, nSensor);
        upper = expand_limit_vector_local(cfg.deltaGapUpperMm, nSensor);
        theta(idxDg) = max(min(theta(idxDg), upper), lower);
    else
        deltaGapLimit = expand_limit_vector_local(cfg.deltaGapLimitMm, nSensor);
        theta(idxDg) = max(min(theta(idxDg), deltaGapLimit), -deltaGapLimit);
    end
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
        uBase = nestedFixedBase.amplitudeMm .* sin(nestedFixedBase.EO .* bundle.Theta + nestedFixedBase.phaseRad);
    end
    vPred = nestedFixedBase.VPred(:);
    overshoot = zeros(size(bundle.V));
    for is = 1:numel(bundle.sensorIds)
        sid = bundle.sensorIds(is);
        mask = bundle.sensorIndex == is;
        Tpl = get_template_sensor_local(bundle.Template, sid);
        xEval = bundle.X(mask) - dx - u(mask) - eta(is);
        xBase = bundle.X(mask) - nestedFixedBase.dxMm - uBase(mask) - eta(is);
        if is_gap_sensor_local(sid, cfg)
            corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
            [vNoGap, overNoGap] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, 0, 0, 0, xEval);
            [vBaseNoGap, overBaseNoGap] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, 0, 0, 0, xBase);
            [deltaV, overshoot(mask)] = eval_static_warp_delta_local( ...
                bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval);
            vPred(mask) = vPred(mask) + (vNoGap - vBaseNoGap) + deltaV;
            overshoot(mask) = max(max(overshoot(mask), overNoGap), overBaseNoGap);
        else
            [vNoGap, overNoGap] = eval_direct_or_base_template_with_overshoot_local(Tpl, xEval);
            [vBaseNoGap, overBaseNoGap] = eval_direct_or_base_template_with_overshoot_local(Tpl, xBase);
            vPred(mask) = vPred(mask) + (vNoGap - vBaseNoGap);
            overshoot(mask) = max(overNoGap, overBaseNoGap);
            dg(is) = 0; dmu(is) = 0; dtau(is) = 0;
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
        if is_gap_sensor_local(sid, cfg)
            corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
            [vPred(mask), overshoot(mask)] = eval_static_warp_template_local(Tpl, ...
                bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval);
        else
            [vPred(mask), overshoot(mask)] = eval_direct_or_base_template_with_overshoot_local(Tpl, xEval);
            dg(is) = 0; dmu(is) = 0; dtau(is) = 0;
        end
    end
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
domainMaxViolationMm = max(overshoot, [], 'omitnan');
if isempty(domainMaxViolationMm) || ~isfinite(domainMaxViolationMm)
    domainMaxViolationMm = inf;
end
domainFeasible = domainMaxViolationMm <= cfg.domainToleranceMm;
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
deltaGapPriorScale = expand_limit_vector_local(cfg.deltaGapPriorScaleMm, nSensor).';
regularizationPenalty = weightNormalizer * (cfg.staticRegWeightMv * sqrt(( ...
    sum((dg ./ max(deltaGapPriorScale, eps)).^2) + ...
    sum((dmu ./ max(cfg.deltaMuLimit, eps)).^2) + ...
    sum((dtau ./ max(cfg.deltaTauLimitMm, eps)).^2)) / max(nSensor, 1))) ^ 2;
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
detail.dxReferenceMm = dxReference;
detail.deltaDxMm = deltaDx;
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
detail.domainFeasible = domainFeasible;
detail.domainMaxViolationMm = domainMaxViolationMm;
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

function [model, overshoot] = eval_static_warp_template_local(Tpl, responseSurface, corr, dg, dmu, dtau, xOpr)
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
if abs(dg) <= eps && abs(dmu) <= eps
    model = vLow;
    overshoot = overshootLow;
    return;
end
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
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

function [model, overshoot] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr)
xLowRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xLowRaw, xLo), xHi);
overshootLow = max(xLo - xLowRaw, 0) + max(xLowRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
if abs(g - corr.g0Mm) <= eps
    model = vLow;
    overshoot = overshootLow;
    return;
end
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

function tf = is_gap_sensor_local(sid, cfg)
if isfield(cfg, 'gapSensors')
    tf = ismember(double(sid), double(cfg.gapSensors(:).'));
else
    tf = true;
end
end

function [model, overshoot] = eval_direct_or_base_template_with_overshoot_local(Tpl, xOpr)
xRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xq = min(max(xRaw, xLo), xHi);
overshoot = max(xLo - xRaw, 0) + max(xRaw - xHi, 0);
model = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xq, 'pchip', NaN);
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
dxRef = selectionInfo.CausalDxReference;
F = wr.modelFits;
row = table(wr.windowId, wr.lapRange(1), wr.lapRange(end), wr.rotFreqMeanHz, ...
    selectionMode, selectedEOText, vpEOText, directCandidateEO, ...
    prevCandidateEO, prevReferenceSource, prevReferenceReason, usedPrevWindowReference, ...
    phaseSafeReferenceEO, phaseSafeReferenceSource, phaseSafeReferenceReason, ...
    usedPhaseSafeExpansion, phaseSafeCorePointCount, phaseSafeExpandedPointCount, phaseSafeAddedPointCount, ...
    dxRef.dxRef, string(dxRef.mode), dxRef.confirmedNow, dxRef.resetNow, ...
    dxRef.selectedEO, dxRef.selectedRank, dxRef.selectedDx, dxRef.dxInnovation, dxRef.costMargin, ...
    'VariableNames', {'window_id','lap_start','lap_end','rot_freq_hz', ...
    'eo_candidate_mode','selected_eo','vp_candidate_eo','direct_candidate_eo', ...
    'prev_candidate_eo','prev_reference_source','prev_reference_reason','used_prev_window_reference', ...
    'phase_safe_reference_eo','phase_safe_reference_source','phase_safe_reference_reason', ...
    'used_phase_safe_expansion','core_point_count','phase_safe_point_count','phase_safe_added_point_count', ...
    'dxref_mm','dxref_mode','dxref_confirmed_now','dxref_reset_now','dxref_selected_EO','dxref_selected_rank','dxref_selected_dx_mm','dxref_innovation_mm','dxref_cost_margin'});
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
    row.gap_waveform_rmse_mV = F.gap_only.plainRmseMv;
    row.gap_objective = F.gap_only.score;
    row.gap_domain_feasible = F.gap_only.domainFeasible;
    row.gap_domain_max_violation_mm = F.gap_only.domainMaxViolationMm;
    row.gap_integer_frequency_hz = F.gap_only.integerFrequencyHz;
    row.gap_frequency_offset_hz = F.gap_only.frequencyOffsetHz;
    row.gap_mean_delta_gap_mm = mean(F.gap_only.deltaGapMm, 'omitnan');
    selectedGapProj = projection_table_for_eo_local(gapProj, F.gap_only.EO);
    row.delta_g_proj_alpha = gapProj.alpha;
    row.delta_g_proj_mean_mm = mean(selectedGapProj.deltaGapProjMm, 'omitnan');
    row.delta_g_proj_max_mm = max(selectedGapProj.deltaGapProjMm, [], 'omitnan');
    row.delta_g_loc_mean_mm = mean(selectedGapProj.deltaGapLocMm, 'omitnan');
    row.delta_g_limit_mean_mm = mean(selectedGapProj.deltaGapLimitMm, 'omitnan');
    row.delta_g_limit_max_mm = max(selectedGapProj.deltaGapLimitMm, [], 'omitnan');
    row.epsilon_perp_mean_mV = mean(selectedGapProj.epsilonPerpMv, 'omitnan');
    row.S_g_min_mV_per_mm = min(selectedGapProj.SgMvPerMm, [], 'omitnan');
    row.eta_g_min = min(selectedGapProj.etaGap, [], 'omitnan');
    row.delta_g_proj_fallback_count = sum(selectedGapProj.usedFallback);
    row.delta_g_proj_cap_count = sum(selectedGapProj.clippedMax);
    row.delta_g_center_mm = string(mat2str(selectedGapProj.deltaGapCenterMm.', 8));
    row.delta_g_lower_mm = string(mat2str(selectedGapProj.deltaGapLowerMm.', 8));
    row.delta_g_upper_mm = string(mat2str(selectedGapProj.deltaGapUpperMm.', 8));
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
end
end

function T = projection_table_for_eo_local(proj, eo)
T = proj.Table;
if ~isfield(proj, 'AllTable') || isempty(proj.AllTable) || ~isfinite(eo)
    return;
end
candidate = proj.AllTable(proj.AllTable.EO == round(eo), :);
if ~isempty(candidate)
    T = candidate;
end
end

function Summary = build_summary_table_local(Trend, bestIdx, cfg)
if strcmpi(cfg.runMode, 'comparison')
    spec = {
    "direct_low_template_main", "direct_EO", "direct_frequency_hz", "direct_amplitude_mm", "direct_rmse_mV"
    "fixed", "fixed_EO", "fixed_frequency_hz", "fixed_amplitude_mm", "fixed_rmse_mV"
    "gap_only", "gap_EO", "gap_frequency_hz", "gap_amplitude_mm", "gap_waveform_rmse_mV"
    "gap_tilt", "tilt_EO", "tilt_frequency_hz", "tilt_amplitude_mm", "tilt_rmse_mV"
    };
else
    spec = {
    "gap_only", "gap_EO", "gap_frequency_hz", "gap_amplitude_mm", "gap_waveform_rmse_mV"
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
    'mean_amplitude_mm','mean_waveform_rmse_mV','median_waveform_rmse_mV','best_free_gap_window'});
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
plot_if_column_local(T, 'gap_waveform_rmse_mV', style.red, 's', 'Gap');
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
    if mode_has_tilt_local(cfg.mainModel)
        title(sprintf('CH%d after, dg %.4f, dmu %.4f', sid, dg, dmu));
    else
        title(sprintf('CH%d after, dg %.4f', sid, dg));
    end
    xlabel('x_{OPR}-u(t) (mm)', 'Interpreter', 'tex'); ylabel('mV'); legend('Location', 'southeast', 'Box', 'off');
    format_axes_local(gca, style);
end
export_paper_figure_local(fig, figFile);
end

function plot_static_warp_local(Result, T, cfg, figFile)
style = paper_style_local();
fig = figure('Name', 'Step07J static warp parameters', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 5.0]);
plotTilt = mode_has_tilt_local(cfg.mainModel);
plotShift = mode_has_shift_local(cfg.mainModel);
nRows = 1 + double(plotTilt) + double(plotShift);
fig.Position(4) = 4.5 * nRows;
tiledlayout(nRows, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
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
xlabel('Window');
title(sprintf('Gap increment, %s main model', cfg.mainModel), 'Interpreter', 'none');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

if plotTilt
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
end

if plotShift
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
end
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

function suffix = sanitize_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function tag = format_analysis_time_tag_local(tSec)
if ~isfinite(tSec)
    tag = 'tNaNs';
    return;
end
tag = sprintf('t%07.3fs', tSec);
tag = strrep(tag, '.', 'p');
tag = regexprep(tag, '0+s$', 's');
tag = regexprep(tag, 'p+s$', 's');
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function out = ternary_local(cond, a, b)
if cond
    out = a;
else
    out = b;
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
