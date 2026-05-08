%% Step 0: configure the blind-template high-speed pipeline
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

cfg = struct();
cfg.rootDir = rootDir;
cfg.scriptDir = scriptDir;
cfg.outDir = outDir;
cfg.makeFigures = true;

txtCandidates = dir(fullfile(rootDir, '**', '*2mm*.txt'));
txtNames = lower(string({txtCandidates.name}));
isTarget = ~startsWith(txtNames, "wave");
targetIdx = find(isTarget, 1, 'first');
if isempty(targetIdx)
    error('Could not find the stacked straight-blade 2 mm gap data file under: %s', rootDir);
end
cfg.dataFile = fullfile(txtCandidates(targetIdx).folder, txtCandidates(targetIdx).name);

cfg.g_low = 0.8;
cfg.g_holdout = 1.0;
cfg.R_tip = 62;
cfg.fs = 2e5;
cfg.RPM_high = 3000;
cfg.NumRevs_high = 8;
cfg.alpha_k = deg2rad([55, 116.76, 178.52]);
cfg.noiseRatio = 0.002;

cfg.A_true = [0.25, 0.15];
cfg.f_true = [500, 1300];
cfg.phi_true = [pi/4, -pi/3];

cfg.domainPadRatio = 0.05;
cfg.segment = struct();
cfg.segment.xGridN = 801;
cfg.segment.activeLevel = 0.08;
cfg.segment.minValidFrac = 0.55;
cfg.segment.minPointCount = 40;
cfg.segment.weightFloor = 0.20;

cfg.blind = struct();
cfg.blind.useGainBias = false;
cfg.blind.maxIter = 12;
cfg.blind.tolTemplate = 1e-4;
cfg.blind.smoothWindow = 11;
cfg.blind.dGridCoarse = -0.45:0.01:0.45;
cfg.blind.dFineHalfWidth = 0.015;
cfg.blind.dFineStep = 0.001;
cfg.blind.seedCount = 6;

cfg.vibrationFit = struct();
cfg.vibrationFit.f1Grid = 400:5:650;
cfg.vibrationFit.f2Grid = 1100:5:1400;
cfg.vibrationFit.minFreqSeparation = 20;
cfg.vibrationFit.ampUpperBound = 0.8;
cfg.vibrationFit.maxIter = 400;
cfg.vibrationFit.shiftRoundTol = 1e-6;

cfg.profile = struct();
cfg.profile.method = 'low_template_guided_blind_refine';
cfg.profile.xGridN = 1201;
cfg.profile.smoothWindow = 13;
cfg.profile.maxOuterIter = 3;
cfg.profile.tolParam = 1e-5;
cfg.profile.tolTemplate = 1e-5;
cfg.profile.numInitCandidates = 120;
cfg.profile.numProfileStarts = 5;
cfg.profile.numLinearFullScaleStarts = 0;
cfg.profile.numLinearRobustStarts = 20;
cfg.profile.linearRobustAmpScales = [10, 15];
cfg.profile.numLinearDiverseStarts = 80;
cfg.profile.initAmpScales = [1, 3, 5, 7, 10, 15, 20];
cfg.profile.derivativeFloorFrac = 0.03;
cfg.profile.f1Grid = 400:5:650;
cfg.profile.f2Grid = 1100:5:1400;
cfg.profile.minFreqSeparation = 20;
cfg.profile.ampUpperBound = 0.8;
cfg.profile.maxLsqIter = 300;
cfg.profile.useLowTemplateInit = true;
cfg.profile.lowTemplateGap = cfg.g_low;
cfg.profile.lowInitKeep = 6;
cfg.profile.lowInitRefineCount = 14;
cfg.profile.lowInitUseGainBias = true;
cfg.profile.lowInitAmpScales = [1.0];
cfg.profile.lowInitFreqJitters = 0;
cfg.profile.useLowToHighTemplateTransfer = true;
cfg.profile.transferUseGainBias = true;
cfg.profile.transferResidualSmoothWindow = 21;
cfg.profile.useLocalRefineBounds = true;
cfg.profile.localFreqHalfWidth = 20;
cfg.profile.localAmpScaleLo = 0.60;
cfg.profile.localAmpScaleHi = 1.40;
cfg.profile.edgeGuardHz = 8;
cfg.profile.edgePenaltyWeight = 0.45;
cfg.profile.fitWeightFloor = 0.15;
cfg.profile.fitWeightPower = 1.0;
cfg.profile.transferResidualClipFrac = 0.45;
cfg.profile.selectionCvWeight = 0.65;
cfg.profile.cvCandidateCount = Inf;
cfg.profile.useGlobalProfileSearch = false;
cfg.profile.globalSwarmSize = 160;
cfg.profile.globalMaxIter = 120;
cfg.profile.globalTemplateGridN = 181;
cfg.profile.globalNumRuns = 3;
cfg.profile.globalKeepCount = 6;
cfg.profile.globalUsePatternSearch = true;
cfg.profile.globalPatternMaxIter = 80;
cfg.profile.globalSmoothPenalty = 1e-2;
cfg.profile.globalUseSurrogateSearch = true;
cfg.profile.globalSurrogateMaxEvals = 1800;

overridePath = fullfile(outDir, '__blind_template_config_override.mat');
if isfile(overridePath)
    S = load(overridePath, 'overrideCfg');
    cfg = merge_struct_recursive(cfg, S.overrideCfg);
end

save(fullfile(outDir, 'stage0_config.mat'), 'cfg');

if cfg.makeFigures
    figure('Name', 'Blind-Template Step 0 - Configuration', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [3, 3, 18, 9]);
    axis off;
    text(0.02, 0.88, 'Blind-Template High-Speed Decoupling', 'FontSize', 15, 'FontWeight', 'bold');
    text(0.02, 0.72, sprintf('Hidden high-speed gap: %.1f mm', cfg.g_holdout), 'FontSize', 11);
    text(0.02, 0.62, sprintf('High-speed RPM: %d', cfg.RPM_high), 'FontSize', 11);
    text(0.02, 0.52, sprintf('True frequencies: %.0f Hz, %.0f Hz', cfg.f_true(1), cfg.f_true(2)), 'FontSize', 11);
    text(0.02, 0.42, sprintf('True amplitudes: %.2f mm, %.2f mm', cfg.A_true(1), cfg.A_true(2)), 'FontSize', 11);
    text(0.02, 0.32, sprintf('Segment grid points: %d', cfg.segment.xGridN), 'FontSize', 11);
    text(0.02, 0.22, sprintf('Blind-template max iterations: %d', cfg.blind.maxIter), 'FontSize', 11);
    text(0.02, 0.10, 'This branch uses only high-speed local passing segments to reconstruct the shared template.', 'FontSize', 10);
end

fprintf('[Step 0] Blind-template config saved to: %s\n', fullfile(outDir, 'stage0_config.mat'));

function cfg = merge_struct_recursive(cfg, overrideCfg)
names = fieldnames(overrideCfg);
for i = 1:numel(names)
    name = names{i};
    if isstruct(overrideCfg.(name)) && isfield(cfg, name) && isstruct(cfg.(name))
        cfg.(name) = merge_struct_recursive(cfg.(name), overrideCfg.(name));
    else
        cfg.(name) = overrideCfg.(name);
    end
end
end
