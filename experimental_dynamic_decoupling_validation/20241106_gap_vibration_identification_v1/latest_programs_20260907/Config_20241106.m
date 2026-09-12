function cfg = Config_20241106()
%CONFIG_20241106 Single user-facing configuration for this package.

rootDir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.dataset = '20241106';
cfg.packageVersion = '20241106_gap_vibration_identification_v1';
cfg.methodFreezeId = '20241106_UNIFIED_STRUCTURED_VOLTAGE_20260811';
cfg.methodStatus = 'unified_simulation_experiment_formal_route';
cfg.paths.root = rootDir;
cfg.paths.rawDatasetRoot = fullfile('D:\', '博士-国科', '试验台数据', '新试验', '20241106_2');
cfg.paths.functions = fullfile(rootDir, 'functions');
% Repeat the Unicode path after the legacy assignment above so MATLAB uses
% the intended public raw-data directory regardless of console encoding.
cfg.paths.rawDatasetRoot = fullfile('D:\', '博士-国科', ...
    '试验台数据', '新试验', '20241106_2');
cfg.paths.foundation = fullfile(cfg.paths.functions, 'foundation');
cfg.paths.preparation = fullfile(cfg.paths.functions, 'preparation');
cfg.paths.gapAware = fullfile(cfg.paths.functions, 'gap_aware');
cfg.paths.utilities = fullfile(cfg.paths.functions, 'utilities');
cfg.paths.diagnostics = fullfile(rootDir, 'diagnostics');
cfg.paths.docs = fullfile(rootDir, 'docs');
cfg.paths.tools = fullfile(rootDir, 'tools');
cfg.paths.inputs = fullfile(rootDir, 'inputs');
cfg.paths.calibrationFoundation = fullfile(cfg.paths.inputs, 'calibration', 'foundation');
cfg.paths.calibrationGap = fullfile(cfg.paths.inputs, 'calibration', 'gap');
cfg.paths.calibrationInputs = fullfile(cfg.paths.inputs, 'calibration');
cfg.paths.staticGapWaveformLibrary = fullfile(cfg.paths.inputs, 'calibration', ...
    'static_gap_waveform_library_20260524_2000Hz');
cfg.paths.bladeOffsetFile = fullfile(cfg.paths.inputs, 'calibration', ...
    'blade_offset', 'analysis_04_blade_offset_summary_matched.csv');
cfg.paths.bundledPreparedFoundation = fullfile(cfg.paths.inputs, 'prepared', 'foundation');
cfg.paths.bundledPreparedGap = fullfile(cfg.paths.inputs, 'prepared', 'gap');
cfg.paths.results = fullfile(rootDir, 'results');
cfg.paths.prepared = fullfile(cfg.paths.results, 'prepared');
cfg.paths.preparedFoundation = fullfile(cfg.paths.prepared, 'foundation');
cfg.paths.gapRuntime = fullfile(cfg.paths.prepared, 'gap_runtime');
cfg.paths.foundationResults = fullfile(cfg.paths.results, 'foundation_no_gap');
cfg.paths.gapResults = fullfile(cfg.paths.results, 'gap_aware');
cfg.paths.comparison = fullfile(cfg.paths.results, 'comparison');
cfg.paths.figures = fullfile(cfg.paths.results, 'figures');
cfg.paths.calibrationRuntime = fullfile(cfg.paths.results, 'calibration');
cfg.paths.strainEvidence = fullfile(cfg.paths.results, 'strain_evidence');
cfg.paths.strainValidation = fullfile(cfg.paths.results, 'strain_validation');
cfg.paths.analysis = fullfile(rootDir, 'analysis');
cfg.paths.common = fullfile(fileparts(fileparts(rootDir)), 'common');

cfg.case.targetBlade = 4;
cfg.case.analysisSensors = [2 5 7];
cfg.case.gapSensors = [5 7];
cfg.case.analysisStartTimeSec = 86.6;
cfg.case.dynamicCase = '3000_3150';

cfg.window.targetBladePasses = 20;
cfg.window.windowBladePasses = 3;
cfg.window.slidingStepBladePasses = 1;
cfg.window.count = floor((cfg.window.targetBladePasses - ...
    cfg.window.windowBladePasses) / cfg.window.slidingStepBladePasses) + 1;

cfg.frequency.searchHz = [300 1000];
cfg.frequency.vpTopK = 3;
cfg.frequency.candidateMethod = 'structured_voltage_vp';
cfg.frequency.structureMode = 'single_sync';
cfg.frequency.structureSelection = 'known_from_experimental_design';
cfg.frequency.unknownStructureCriterion = 'bic_after_complete_voltage_optimization';
cfg.frequency.candidateTopK = 3;
cfg.frequency.candidateMaxTopK = 3;
cfg.frequency.refineEachCandidate = true;
cfg.frequency.refineHalfWidthHz = 2.0;

% The formal low-speed templates use one common physical lap segment for
% every participating sensor.  The bundled template was generated from
% laps 1--30; keep the provenance in the public configuration so a rebuild
% cannot silently mix sensor-specific windows.
cfg.calibration.commonLowSpeedLapRange = [1 30];
cfg.calibration.requireCommonLowSpeedWindow = true;
cfg.calibration.templateSource = 'bundled_step04_template';

cfg.files.lowSpeedTemplateSource = fullfile(cfg.paths.calibrationFoundation, ...
    'step04_low_speed_template', 'LowSpeed_Template_SourceData_20241106.mat');
cfg.files.lowSpeedTemplateLegacy = fullfile(cfg.paths.calibrationFoundation, ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat');
cfg.files.lowSpeedTemplate = fullfile(cfg.paths.calibrationFoundation, ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S2357_20241106.mat');
cfg.files.responseSurface = fullfile(cfg.paths.inputs, 'calibration', ...
    'Step05_Response_Surface_20241106.mat');
cfg.files.gapLibrary = fullfile(cfg.paths.calibrationGap, ...
    'GapCalibrationBank_20241106_B1toB6_S57.mat');
cfg.lowSpeed.templateMethod = 'support_aware_grouped_lap_cv_sg_pchip';
cfg.lowSpeed.interpolationMethod = 'pchip_no_extrapolation';
cfg.lowSpeed.candidateWindowMm = [0.10 0.18 0.34 0.50 0.82 1.22 1.62 2.02 2.42];
cfg.lowSpeed.groupedFoldCount = 5;
cfg.lowSpeed.minBinCount = 5;

% Weighting is retained only for pseudo-displacement/VP initialization.
% The final forward voltage fit and EO ranking use ordinary voltage RMSE.
cfg.objective.finalType = 'plain_rmse';
cfg.objective.vpType = 'gradient_voltage_linearized_displacement';
cfg.objective.candidateType = 'structured_weighted_voltage_candidate';
cfg.objective.finalResidualWeighting = 'none';

cfg.model.etaPolicy = 'zero';
cfg.model.etaMm = 0;
cfg.model.mainGapModel = 'gap_only';
cfg.model.forwardModel = 'low_template_plus_static_gap_increment';
cfg.model.deltaGapLimitMm = 0.25;
cfg.model.deltaMuLimit = 0.060;
cfg.model.deltaTauLimitMm = 0.16;

% R5 keeps the legacy gap_only dynamic contract.  The anchor registration
% and sensor-conditioned latent are frozen before the window fits; only the
% vibration state and independent dg_s are window-varying.
cfg.r5 = struct();
cfg.r5.anchorBlade = 2;
cfg.r5.latentMode = 'sensor_conditioned';
% Keep the low-speed calibrated gap path in the frozen response library.
% No high-speed tilt parameter is fitted; the only clearance variable in a
% dynamic window is the independent increment dg_s.
% The low-speed gap path is part of the frozen observation operator.  It is
% not a dynamic tilt degree of freedom: only the independent dg_s increment
% is fitted in the high-speed window.
% Localization identified a scalar sensor-conditioned effective gap.
% Legacy mu/tau is not part of this response surface coordinate.
cfg.r5.gapBaselineModel = 'sensor_conditioned_scalar';
cfg.r5.requireStaticAudit = true;
cfg.r5.minSupportFraction = 0.98;
cfg.r5.amplitudeLimitMm = 0.50;
cfg.r5.dxLimitMm = 0.35;
cfg.r5.deltaGapLimitMm = cfg.model.deltaGapLimitMm;
cfg.r5.useLocalGapProjection = false;
cfg.r5.gapProjectionAlpha = 2.5;
cfg.r5.gapProjectionMinLimitMm = 0.02;
cfg.r5.gapProjectionSensitivityFloorMvPerMm = 0.05;
cfg.r5.gapProjectionStepMm = 1e-3;
cfg.r5.useFoundationVibrationSeedOnly = true;
% Diagnostic switch for the legacy nested Foundation anchor route.
cfg.r5.useNestedFoundationAnchor = false;
cfg.r5.useV1DynamicIncrement = false;
cfg.r5.auditFoundationReplay = true;
cfg.r5.foundationReplayToleranceMv = 1e-3;
% The target offset is a localization/audit quantity.  The bundled template
% and response surface already share the registered OPR coordinate.
cfg.r5.multiStartCount = 5;
cfg.audit.referenceEO = 12;
% EO remains a data-driven candidate/refinement variable in the production
% route. A fixed EO may be used only by separate diagnostic scripts.
cfg.r5.lockReferenceEO = false;
% Optional diagnostic: refit all EO candidates on selected windows to audit
% whether the production Top-K screen recalls the full-wave winner. It is
% off for formal runs and never constrains the reported EO.
cfg.r5.auditAllEO = false;
cfg.r5.auditAllEOWindowIds = [];

cfg.run.showPlots = false;
cfg.run.saveFigures = false;
cfg.run.forceFoundationRebuild = false;
cfg.preparation.rebuildFromRaw = false;

% The strain chain is a formal physical-validation stage, not an EO prior.
cfg.strain.sampleRateHz = 10000;
cfg.strain.primaryChannel = 'AI1-04';
cfg.strain.timeOffsetPriorSec = -102.6;
cfg.strain.frequencyBandHz = [300 1000];
cfg.strain.alignmentCasesSec = [75.0 86.6];

% Optional diagnostic-only override. Formal runs leave this environment
% variable empty and therefore retain the frozen configuration above.
overrideFile = strtrim(getenv('BTT_PARAMETER_STUDY_OVERRIDE_FILE'));
if ~isempty(overrideFile)
    assert(isfile(overrideFile), ...
        'Parameter-study override file does not exist: %s', overrideFile);
    O = load(overrideFile, 'StudyOverride', 'StudyMetadata');
    assert(isfield(O, 'StudyOverride') && isstruct(O.StudyOverride), ...
        'Override file must contain a StudyOverride struct.');
    cfg = merge_struct_local(cfg, O.StudyOverride);
    cfg.window.count = floor((cfg.window.targetBladePasses - ...
        cfg.window.windowBladePasses) / ...
        cfg.window.slidingStepBladePasses) + 1;
    assert(cfg.window.count >= 1, ...
        'Parameter-study window settings produce no valid windows.');
    assert(strcmpi(cfg.model.etaPolicy, 'zero') && cfg.model.etaMm == 0, ...
        'Parameter studies may not change the frozen eta=0 policy.');
    assert(strcmpi(cfg.objective.finalType, 'plain_rmse'), ...
        'Parameter studies may not change the formal plain-RMSE objective.');
    cfg.parameterStudy.overrideFile = overrideFile;
    if isfield(O, 'StudyMetadata')
        cfg.parameterStudy.metadata = O.StudyMetadata;
    end
end
end


function output = merge_struct_local(base, override)
output = base;
names = fieldnames(override);
for i = 1:numel(names)
    name = names{i};
    if isfield(output, name) && isstruct(output.(name)) && ...
            isstruct(override.(name)) && isscalar(output.(name)) && ...
            isscalar(override.(name))
        output.(name) = merge_struct_local(output.(name), override.(name));
    else
        output.(name) = override.(name);
    end
end
end
