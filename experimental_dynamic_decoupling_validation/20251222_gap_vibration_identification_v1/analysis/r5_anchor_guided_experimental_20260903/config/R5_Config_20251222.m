function cfg = R5_Config_20251222()
%R5_CONFIG_20251222 Read-only configuration for the R5 sidecar analysis.
% The production package is loaded only to obtain paths and frozen contracts.

thisDir = fileparts(mfilename('fullpath'));
cfg.r5Root = fileparts(thisDir);
cfg.analysisRoot = fileparts(cfg.r5Root);
cfg.productionRoot = fileparts(cfg.analysisRoot);

addpath(cfg.productionRoot);
cleanupPath = onCleanup(@() rmpath(cfg.productionRoot));
prod = Config_20251222(); %#ok<NASGU>
clear cleanupPath

cfg.dataset = prod.dataset;
cfg.production = prod;
cfg.case = prod.case;
cfg.paths = prod.paths;
cfg.paths.r5Root = cfg.r5Root;
cfg.paths.r5Results = fullfile(cfg.r5Root, 'results');
cfg.paths.r5Tables = fullfile(cfg.paths.r5Results, 'tables');
cfg.paths.r5Figures = fullfile(cfg.paths.r5Results, 'figures');
cfg.paths.r5Logs = fullfile(cfg.paths.r5Results, 'logs');
cfg.outputs.root = cfg.paths.r5Results;

cfg.inputs.responseSurface = fullfile(cfg.paths.prepared, 'calibration', ...
    'Step05_Response_Surface_20251222_RefBladeAnchor.mat');
cfg.inputs.responseSurfaceFallback = prod.files.responseSurface;
cfg.inputs.lowSpeedTemplate = prod.files.lowSpeedTemplate;
cfg.inputs.dynamicMap = '';
cfg.inputs.correctedGapLibrary = prod.files.gapLibrary;

cfg.state.b2BladeId = 2;
cfg.state.nonB2AxisName = 'effective_state';
cfg.state.maxCandidateCount = 9;
cfg.state.allowContinuousRefine = false;
cfg.state.continuityLambda = 0;

cfg.donor.topK = 3;
cfg.donor.distanceFloor = 1e-9;
cfg.donor.weights = struct('voltage', 1.0, 'gradient', 0.25, ...
    'phase', 0.0, 'amplitude', 0.25);
cfg.donor.useFixedTrainingScale = true;
cfg.donor.manifoldMaxDistance = inf;

cfg.profile.frequencyHz = prod.frequency.searchHz;
cfg.profile.frequencyGridHz = linspace(prod.frequency.searchHz(1), ...
    prod.frequency.searchHz(2), 41); % coarse state scan; refine selected EO separately
cfg.profile.eoPad = 0;
cfg.profile.dxBoundsMm = [-0.25 0.25];
cfg.profile.amplitudeBoundsMm = [0 0.50];
cfg.profile.phaseBoundsRad = [-pi pi];
cfg.profile.maxIterations = 100;
cfg.profile.objective = 'unweighted_voltage_sse';

cfg.gates.maxAbsRhoSx = 0.97;
cfg.gates.maxJacobianCondition = 1e8;
cfg.gates.minProfileMarginMv2 = 0;
cfg.gates.maxStateJump = inf;
cfg.gates.minSupportFraction = 0.60;

cfg.opr = prod.opr;
cfg.window = prod.window;
cfg.machine = prod.machine;
end
