function cfg = Config_20251222()
%CONFIG_20251222 Single user-facing configuration for the formal V1 package.

rootDir = fileparts(mfilename('fullpath'));
utilitiesDir = fullfile(rootDir, 'functions', 'utilities');
if exist(utilitiesDir, 'dir') == 7
    addpath(utilitiesDir);
end

cfg = struct();
cfg.dataset = '20251222';
cfg.packageVersion = '20251222_gap_vibration_identification_v1';
cfg.methodFreezeId = '20251222_UNIFIED_STRUCTURED_VOLTAGE_20260811';
cfg.methodStatus = 'unified_simulation_experiment_formal_route';

cfg.paths.root = rootDir;
cfg.paths.rawDatasetRoot = fullfile('E:\', '试验数据', '20251222', '传感器数据');
cfg.paths.strainRoot = fullfile('E:\', '试验数据', '20251222', '应变片数据');
% All executable dependencies are package-local. These legacy fields are
% retained only for schema compatibility and are intentionally empty.
cfg.paths.staticGapLibrary = '';
cfg.paths.referenceBladeGapAnalysis = '';
cfg.paths.functions = fullfile(rootDir, 'functions');
cfg.paths.preparation = fullfile(cfg.paths.functions, 'preparation');
cfg.paths.foundation = fullfile(cfg.paths.functions, 'foundation');
cfg.paths.calibration = fullfile(cfg.paths.functions, 'calibration');
cfg.paths.gapAware = fullfile(cfg.paths.functions, 'gap_aware');
cfg.paths.utilities = fullfile(cfg.paths.functions, 'utilities');
cfg.paths.inputs = fullfile(rootDir, 'inputs');
cfg.paths.calibrationInputs = fullfile(cfg.paths.inputs, 'calibration');
cfg.paths.preparedInputs = fullfile(cfg.paths.inputs, 'prepared');
cfg.paths.preparedFoundation = fullfile(cfg.paths.preparedInputs, 'foundation');
cfg.paths.results = fullfile(rootDir, 'results');
cfg.paths.prepared = fullfile(cfg.paths.results, 'prepared');
cfg.paths.preparedLegacy = fullfile(cfg.paths.prepared, 'legacy_btt');
cfg.paths.calibrationWork = fullfile(cfg.paths.prepared, 'calibration');
cfg.paths.foundationResults = fullfile(cfg.paths.results, 'fixed_gap');
cfg.paths.gapResults = fullfile(cfg.paths.results, 'gap_aware');
cfg.paths.comparison = fullfile(cfg.paths.results, 'comparison');
cfg.paths.strainValidation = fullfile(cfg.paths.results, 'strain_validation');
cfg.paths.figures = fullfile(cfg.paths.results, 'figures');
cfg.paths.analysis = fullfile(rootDir, 'analysis');

cfg.case.dynamicCase = '1000_2500_3500';
cfg.case.lowSpeedCase = '1000rpm无振动';
[selection, cfg.case.resonanceCatalog] = ResonanceRegionCatalog_20251222();
cfg.case.resonance = selection;
cfg.case.targetBlade = selection.representativeBladeId;
cfg.case.analysisSensors = [1 2 3];
cfg.case.gapSensors = [1 2 3];
cfg.case.analysisStartTimeSec = selection.analysisStartTimeSec;

cfg.machine.bladeCount = 6;
cfg.machine.oprPulsesPerRevolution = 6;
cfg.machine.sampleRateHz = 5e6;
cfg.machine.tipRadiusMm = 62.0;

cfg.window.targetBladePasses = selection.targetBladePasses;
cfg.window.windowBladePasses = selection.windowBladePasses;
cfg.window.slidingStepBladePasses = selection.slidingStepBladePasses;
cfg.window.count = floor((cfg.window.targetBladePasses - ...
    cfg.window.windowBladePasses) / cfg.window.slidingStepBladePasses) + 1;
cfg.window.pulseWindowSec = 6e-4;

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
cfg.frequency.foundationCandidatePolicy = 'adaptive_full_fit';
cfg.frequency.gapCandidatePolicy = 'gapaware_vp_top3_only';

cfg.objective.vpType = 'gradient_voltage_linearized_displacement';
cfg.objective.candidateType = 'structured_weighted_voltage_candidate';
cfg.objective.finalType = 'plain_rmse';
cfg.objective.finalResidualWeighting = 'none';
cfg.model.etaPolicy = 'zero';
cfg.model.etaMm = 0;
cfg.model.mainGapModel = 'gap_only';
cfg.model.forwardModel = 'low_template_plus_static_gap_increment';
cfg.model.deltaGapLimitMm = 0.25;
cfg.model.useDeltaGapProjectionLimit = true;
cfg.model.deltaGapProjectionAlpha = 2.5;
cfg.model.deltaGapProjectionMinMm = 0.02;
cfg.model.deltaMuLimit = 0.005;
cfg.model.deltaTauLimitMm = 0.16;
cfg.independence.usePreviousWindowCandidate = false;
cfg.independence.useCausalDxState = false;
cfg.independence.dxReferenceMode = 'seed_median';
cfg.bundle.foundationRoleInGapAware = 'common_raw_waveform_bundle_only';
cfg.bundle.inheritFoundationEOAmplitudePhaseDx = false;

cfg.opr.lowSpeedLapRange = [2 21];
cfg.opr.slotAngleDeg = [59.6191671794336, 60.4147103395678, ...
    59.9241580729946, 59.6366741842815, 60.3594476244228, 60.0458425992998];
cfg.opr.highToLowSlotOffset = -2;
cfg.files.oprSlotCalibration = fullfile(cfg.paths.calibrationInputs, ...
    'OPR6_SlotAngleCalibration_Laps2_21_20251222.mat');
cfg.files.lowSpeedTemplateSource = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', 'LowSpeed_Template_SourceData_20251222.mat');
cfg.files.lowSpeedTemplateLegacy = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');
cfg.files.lowSpeedTemplate = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', 'Template_AdaptiveSG.mat');
cfg.files.responseSurface = fullfile(cfg.paths.calibrationInputs, ...
    'Step05_Response_Surface_20251222.mat');
cfg.files.sharedResponseSurface = fullfile(cfg.paths.calibrationInputs, ...
    'Step05I_OffsetTilt_Shared_Response_Surface_20251222.mat');
cfg.files.gapLibrary = fullfile(cfg.paths.calibrationInputs, sprintf( ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_S123.mat', ...
    cfg.case.targetBlade));

cfg.strain.primaryChannel = 'AI1-03';
cfg.strain.frequencyBandHz = [520 640];
cfg.strain.transferRatio1PerM = 1.266;
cfg.lowSpeed.templateMethod = 'support_aware_grouped_lap_cv_sg_pchip';
cfg.lowSpeed.interpolationMethod = 'pchip_no_extrapolation';
cfg.lowSpeed.candidateWindowMm = [0.10 0.18 0.34 0.50 0.82 1.22 1.62 2.02 2.42];
cfg.lowSpeed.groupedFoldCount = 5;
cfg.lowSpeed.minBinCount = 5;
cfg.run.showPlots = false;
cfg.run.saveFigures = false;
cfg.run.forceRebuild = false;
% The low-speed calibration contains the explicit frozen spatial gap profile.
% It is used in both zero-dynamic and dynamic evaluations; high-speed fits
% may change dg only, never mu or tau.
cfg.r5.useStaticGapSlopeInDynamic = true;
cfg.r5.staticTiltCarrier = 'explicit_profile';
cfg.r5.staticTiltSource = 'low_speed_sensor_conditioned_calibration';
end
