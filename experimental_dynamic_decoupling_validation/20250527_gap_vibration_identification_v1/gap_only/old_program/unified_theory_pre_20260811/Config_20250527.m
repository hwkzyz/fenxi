function cfg = Config_20250527()
%CONFIG_20250527 Single user-facing configuration for the formal V1 package.

rootDir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.dataset = '20250527';
cfg.packageVersion = '20250527_gap_vibration_identification_v1';

cfg.paths.root = rootDir;
cfg.paths.rawDatasetRoot = fullfile('E:\', '试验数据', '20250527', '试验20250527');
cfg.paths.staticGapLibrary = fullfile(fileparts(rootDir), ...
    'static_gap_waveform_library_20260524_2000Hz');
cfg.paths.referenceBladeGapAnalysis = fullfile(fileparts(rootDir), ...
    'reference_blade_gap_analysis');
cfg.paths.functions = fullfile(rootDir, 'functions');
cfg.paths.preparation = fullfile(cfg.paths.functions, 'preparation');
cfg.paths.foundation = fullfile(cfg.paths.functions, 'foundation');
cfg.paths.calibration = fullfile(cfg.paths.functions, 'calibration');
cfg.paths.gapAware = fullfile(cfg.paths.functions, 'gap_aware');
cfg.paths.utilities = fullfile(cfg.paths.functions, 'utilities');
cfg.paths.inputs = fullfile(rootDir, 'inputs');
cfg.paths.calibrationInputs = fullfile(cfg.paths.inputs, 'calibration');
cfg.paths.preparedInputs = fullfile(cfg.paths.inputs, 'prepared');
cfg.paths.results = fullfile(rootDir, 'results');
cfg.paths.prepared = fullfile(cfg.paths.results, 'prepared');
cfg.paths.foundationResults = fullfile(cfg.paths.results, 'fixed_gap');
cfg.paths.gapResults = fullfile(cfg.paths.results, 'gap_aware');
cfg.paths.comparison = fullfile(cfg.paths.results, 'comparison');
cfg.paths.strainValidation = fullfile(cfg.paths.results, 'strain_validation');
cfg.paths.figures = fullfile(cfg.paths.results, 'figures');
cfg.paths.analysis = fullfile(rootDir, 'analysis');

cfg.files.formalGapResult = fullfile(cfg.paths.gapResults, ...
    'GapAware_B1_S136_T001p5s_W3S1_DF2Hz.mat');
cfg.files.formalFoundationResult = fullfile(cfg.paths.foundationResults, ...
    '20250526_2500-3500_t400', ...
    'FixedGap_B1_S136_T001p5s.mat');

cfg.case.dynamicCase = '20250526_2500-3500_t400';
cfg.case.lowSpeedCase = '20250526_910';
cfg.case.targetBlade = 1;
cfg.case.analysisSensors = [1 3 6];
cfg.case.gapSensors = [1 3 6];
cfg.case.analysisStartTimeSec = 1.5;
analysisStartOverride = strtrim(getenv('BTT_ANALYSIS_START_TIME_SEC'));
if ~isempty(analysisStartOverride)
    parsedStart = str2double(analysisStartOverride);
    if ~isfinite(parsedStart) || parsedStart < 0
        error('BTT_ANALYSIS_START_TIME_SEC must be a finite nonnegative number.');
    end
    cfg.case.analysisStartTimeSec = parsedStart;
end

cfg.machine.bladeCount = 6;
cfg.machine.oprPulsesPerRevolution = 6;
cfg.machine.sampleRateHz = 5e6;
cfg.machine.tipRadiusMm = 62.0;

cfg.window.targetBladePasses = 20;
cfg.window.windowBladePasses = 3;
cfg.window.slidingStepBladePasses = 1;
cfg.window.count = floor((cfg.window.targetBladePasses - ...
    cfg.window.windowBladePasses) / cfg.window.slidingStepBladePasses) + 1;
cfg.window.pulseWindowSec = 6e-4;

cfg.frequency.searchHz = [300 1000];
cfg.frequency.vpTopK = 3;
cfg.frequency.refineEachCandidate = true;
cfg.frequency.refineHalfWidthHz = 2.0;

cfg.objective.vpType = 'weighted_displacement';
cfg.objective.finalType = 'plain_rmse';

cfg.model.etaPolicy = 'zero';
cfg.model.etaMm = 0;
cfg.model.mainGapModel = 'gap_only';
cfg.model.deltaGapLimitMm = 0.25;
cfg.model.useDeltaGapProjectionLimit = true;
cfg.model.deltaGapProjectionAlpha = 2.5;
cfg.model.deltaGapProjectionMinMm = 0.02;
cfg.model.deltaMuLimit = 0.005;
cfg.model.deltaTauLimitMm = 0.16;

cfg.independence.usePreviousWindowCandidate = false;
cfg.independence.useCausalDxState = false;

cfg.strain.primaryChannel = 'AI1-04';
cfg.strain.frequencyBandHz = [520 640];
cfg.strain.transferRatio1PerM = 1.010;
cfg.strain.transferRatioStatus = 'provisional_580Hz_FE_gauge_center';

cfg.run.showPlots = false;
cfg.run.saveFigures = false;
cfg.run.forceRebuild = false;
end
