function cfg = Config_20241106()
%CONFIG_20241106 Single user-facing configuration for this package.

rootDir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.dataset = '20241106';
cfg.packageVersion = '20241106_gap_vibration_identification_v1';
cfg.paths.root = rootDir;
cfg.paths.rawDatasetRoot = fullfile('D:\', '博士-国科', '试验台数据', '新试验', '20241106_2');
cfg.paths.src = fullfile(rootDir, 'src');
% Repeat the Unicode path after the legacy assignment above so MATLAB uses
% the intended public raw-data directory regardless of console encoding.
cfg.paths.rawDatasetRoot = fullfile('D:\', '博士-国科', ...
    '试验台数据', '新试验', '20241106_2');
cfg.paths.foundation = fullfile(cfg.paths.src, 'foundation');
cfg.paths.preparation = fullfile(cfg.paths.src, 'preparation');
cfg.paths.gapAware = fullfile(cfg.paths.src, 'gap_aware');
cfg.paths.utilities = fullfile(cfg.paths.src, 'utilities');
cfg.paths.inputs = fullfile(rootDir, 'inputs');
cfg.paths.calibrationFoundation = fullfile(cfg.paths.inputs, 'calibration', 'foundation');
cfg.paths.calibrationGap = fullfile(cfg.paths.inputs, 'calibration', 'gap');
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
cfg.frequency.refineEachCandidate = true;
cfg.frequency.refineHalfWidthHz = 2.0;

% Weighting is retained only for pseudo-displacement/VP initialization.
% The final forward voltage fit and EO ranking use ordinary voltage RMSE.
cfg.objective.finalType = 'plain_rmse';

cfg.model.etaPolicy = 'zero';
cfg.model.etaMm = 0;
cfg.model.deltaGapLimitMm = 0.25;
cfg.model.deltaMuLimit = 0.060;
cfg.model.deltaTauLimitMm = 0.16;

cfg.run.showPlots = false;
cfg.run.saveFigures = false;
cfg.run.forceFoundationRebuild = false;
cfg.preparation.rebuildFromRaw = false;
end
