%% Step 0: configure the holdout validation
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
if ~isfolder(outDir)
    mkdir(outDir);
end

cfg = struct();
cfg.rootDir = rootDir;
cfg.scriptDir = scriptDir;
cfg.outDir = outDir;
cfg.makeFigures = true;
cfg.computeStep4Baselines = true;
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
cfg.RPM_low = 900;
cfg.RPM_high = 3000;
cfg.NumRevs_low = 4;
cfg.NumRevs_high = 8;
cfg.alpha_k = deg2rad([55, 116.76, 178.52]);
cfg.noiseRatio = 0.002;

cfg.A_true = [0.25, 0.15];
cfg.f_true = [500, 1300];
cfg.phi_true = [pi/4, -pi/3];

cfg.xGridN = 2201;
cfg.gapQueryN = 241;
cfg.dxStaticGrid = linspace(-0.6, 0.6, 161);
cfg.fitActiveLevel = 0.08;
cfg.f1Grid = 400:5:650;
cfg.f2Grid = 1100:5:1400;
cfg.numVarproCandidates = 10;
cfg.rawGapSearchMargin = 0.25;
cfg.rawProfileGridN = 401;
cfg.rawProfilePercentile = 65;
cfg.rawProfileMinCount = 2;
cfg.rawProfileSmoothWindow = 9;
cfg.rawHuberDelta = 1.5;
cfg.rawWeightFloor = 0.08;
cfg.jointRefine.gHalfWidth = 0.05;
cfg.jointRefine.dxHalfWidth = 0.05;
cfg.jointRefine.fHalfWidth = 20;
cfg.jointRefine.ampScaleLo = 0.50;
cfg.jointRefine.ampScaleHi = 1.50;

if ~isfile(cfg.dataFile)
    error('Data file not found: %s', cfg.dataFile);
end

overridePath = fullfile(outDir, '__gapaware_config_override.mat');
if isfile(overridePath)
    S = load(overridePath, 'overrideCfg');
    cfg = merge_struct_recursive(cfg, S.overrideCfg);
end

save(fullfile(outDir, 'stage0_config.mat'), 'cfg');

if cfg.makeFigures
    figure('Name', 'Step 0 - Configuration', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [3, 3, 18, 9]);
    axis off;
    text(0.02, 0.88, 'Gap-Aware Holdout Validation', 'FontSize', 15, 'FontWeight', 'bold');
    text(0.02, 0.72, sprintf('Low-speed reference gap: %.1f mm', cfg.g_low), 'FontSize', 11);
    text(0.02, 0.62, sprintf('Hidden high-speed gap: %.1f mm', cfg.g_holdout), 'FontSize', 11);
    text(0.02, 0.52, sprintf('RPM low / high: %d / %d', cfg.RPM_low, cfg.RPM_high), 'FontSize', 11);
    text(0.02, 0.42, sprintf('True frequencies: %.0f Hz, %.0f Hz', cfg.f_true(1), cfg.f_true(2)), 'FontSize', 11);
    text(0.02, 0.32, sprintf('True amplitudes: %.2f mm, %.2f mm', cfg.A_true(1), cfg.A_true(2)), 'FontSize', 11);
    text(0.02, 0.22, sprintf('Equivalent sensor phases: %s deg', mat2str(rad2deg(cfg.alpha_k), 4)), 'FontSize', 11);
    text(0.02, 0.10, 'This step only defines the blind-test design and saves cfg for later steps.', 'FontSize', 10);
end

fprintf('[Step 0] Config saved to: %s\n', fullfile(outDir, 'stage0_config.mat'));

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

