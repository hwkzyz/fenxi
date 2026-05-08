%% Explore possible refinements for no-vibration reference-gap reconstruction
% This script does not replace the main method. It checks whether the
% current dense equal-level registration can be improved by changing:
%   1) reference gap;
%   2) equal-level anchor range and density;
%   3) interpolation scheme;
%   4) amplitude reconstruction rule.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = resolve_stacked_curve_file(scriptDir);
outDir = fullfile(scriptDir, 'no_vibration_reconstruction_refinement_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Figure defaults
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.15, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% Load and preprocess no-vibration waveforms
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);

curve = repmat(struct( ...
    'gap', 0, 'x', [], 'yRaw', [], 'yBc', [], ...
    'base', 0, 'amp', 0, 'xPeak', 0, 'xCenter', [], ...
    'yNorm', []), nGap, 1);

for ig = 1:nGap
    x = xCell{ig}(:);
    y = yCell{ig}(:);
    [base, amp, yNorm] = normalize_waveform(y);
    yBc = y - base;
    [~, iPeak] = max(yBc);
    xPeak = x(iPeak);

    curve(ig).gap = gapList(ig);
    curve(ig).x = x;
    curve(ig).yRaw = y;
    curve(ig).yBc = yBc;
    curve(ig).base = base;
    curve(ig).amp = amp;
    curve(ig).xPeak = xPeak;
    curve(ig).xCenter = x - xPeak;
    curve(ig).yNorm = yNorm;
end

%% Baseline configuration: current main method
nU = 1401;
baselineLevels = (0.02:0.02:0.98)';
[~, defaultRefIdx] = max(gapList);
defaultInterp = 'pchip';
defaultU = make_reference_grid(curve(defaultRefIdx), baselineLevels, nU);

baselineC = register_all_curves(curve, defaultRefIdx, defaultU, baselineLevels, defaultInterp);
baselineMetrics = evaluate_reconstruction(gapList, baselineC, defaultRefIdx, 'peak', ...
    "baseline_current", "Current dense equal-level");
writetable(baselineMetrics.globalTable, fullfile(outDir, 'baseline_current_global.csv'));
writetable(baselineMetrics.summaryTable, fullfile(outDir, 'baseline_current_by_gap.csv'));

%% Experiment 1: reference-gap sweep
refRows(nGap, 1) = struct( ...
    'reference_gap_mm', NaN, ...
    'reference_amplitude_pF', NaN, ...
    'u_span_mm', NaN, ...
    'mean_shape_rmse_to_ref', NaN, ...
    'max_shape_rmse_to_ref', NaN, ...
    'mean_reconstructed_rmse_to_ref', NaN, ...
    'max_reconstructed_rmse_to_ref', NaN, ...
    'mean_reconstructed_nrmse_to_ref', NaN, ...
    'max_reconstructed_nrmse_to_ref', NaN);

for ir = 1:nGap
    uGrid = make_reference_grid(curve(ir), baselineLevels, nU);
    Cmat = register_all_curves(curve, ir, uGrid, baselineLevels, defaultInterp);
    metrics = evaluate_reconstruction(gapList, Cmat, ir, 'peak', ...
        "ref_gap_sweep", sprintf('Reference gap %.1f mm', gapList(ir)));

    refRows(ir).reference_gap_mm = gapList(ir);
    refRows(ir).reference_amplitude_pF = max(Cmat(ir, :), [], 2);
    refRows(ir).u_span_mm = max(uGrid) - min(uGrid);
    refRows(ir).mean_shape_rmse_to_ref = metrics.globalTable.mean_shape_rmse_to_ref;
    refRows(ir).max_shape_rmse_to_ref = metrics.globalTable.max_shape_rmse_to_ref;
    refRows(ir).mean_reconstructed_rmse_to_ref = metrics.globalTable.mean_reconstructed_rmse_to_ref;
    refRows(ir).max_reconstructed_rmse_to_ref = metrics.globalTable.max_reconstructed_rmse_to_ref;
    refRows(ir).mean_reconstructed_nrmse_to_ref = metrics.globalTable.mean_reconstructed_nrmse_to_ref;
    refRows(ir).max_reconstructed_nrmse_to_ref = metrics.globalTable.max_reconstructed_nrmse_to_ref;
end

refSweepTable = struct2table(refRows);
writetable(refSweepTable, fullfile(outDir, 'reference_gap_sweep.csv'));

%% Experiment 2: equal-level anchor sweep on the same evaluation grid
levelSets = { ...
    make_level_set(0.005, 0.995, 199), ...
    make_level_set(0.010, 0.990, 99), ...
    make_level_set(0.020, 0.980, 49), ...
    make_level_set(0.030, 0.970, 63), ...
    make_level_set(0.050, 0.950, 45), ...
    make_level_set(0.100, 0.900, 41)};
levelNames = [ ...
    "0.005-0.995, 199 anchors"; ...
    "0.010-0.990, 99 anchors"; ...
    "0.020-0.980, 49 anchors"; ...
    "0.030-0.970, 63 anchors"; ...
    "0.050-0.950, 45 anchors"; ...
    "0.100-0.900, 41 anchors"];

levelTables = cell(numel(levelSets), 1);
for iset = 1:numel(levelSets)
    levels = levelSets{iset};
    Cmat = register_all_curves(curve, defaultRefIdx, defaultU, levels, defaultInterp);
    levelTables{iset} = evaluate_reconstruction(gapList, Cmat, defaultRefIdx, 'peak', ...
        sprintf("level_set_%02d", iset), levelNames(iset)).globalTable;
end
levelSweepTable = vertcat(levelTables{:});
writetable(levelSweepTable, fullfile(outDir, 'level_set_sweep.csv'));

%% Experiment 3: interpolation scheme sweep
interpMethods = ["linear"; "pchip"; "makima"; "spline"];
interpTables = {};
for im = 1:numel(interpMethods)
    methodName = interpMethods(im);
    try
        Cmat = register_all_curves(curve, defaultRefIdx, defaultU, baselineLevels, char(methodName));
        interpTables{end+1, 1} = evaluate_reconstruction(gapList, Cmat, defaultRefIdx, 'peak', ...
            sprintf("interp_%s", methodName), sprintf('Interpolation: %s', methodName)).globalTable; %#ok<SAGROW>
    catch ME
        warning('Interpolation method %s failed: %s', methodName, ME.message);
    end
end
interpSweepTable = vertcat(interpTables{:});
writetable(interpSweepTable, fullfile(outDir, 'interpolation_sweep.csv'));

%% Experiment 4: combined equal-level and interpolation sweep
comboTables = {};
for iset = 1:3
    for im = 1:numel(interpMethods)
        levels = levelSets{iset};
        methodName = interpMethods(im);
        try
            Cmat = register_all_curves(curve, defaultRefIdx, defaultU, levels, char(methodName));
            label = sprintf('%s + %s', levelNames(iset), methodName);
            comboTables{end+1, 1} = evaluate_reconstruction(gapList, Cmat, defaultRefIdx, 'peak', ...
                sprintf("combo_%02d_%s", iset, methodName), label).globalTable; %#ok<SAGROW>
        catch ME
            warning('Combined setting %d / %s failed: %s', iset, methodName, ME.message);
        end
    end
end
comboSweepTable = vertcat(comboTables{:});
writetable(comboSweepTable, fullfile(outDir, 'combined_level_interpolation_sweep.csv'));

[~, bestComboIdx] = min(comboSweepTable.mean_reconstructed_rmse_to_ref);
[bestLevelIdx, bestInterp] = parse_combo_key(comboSweepTable.method_key(bestComboIdx));
bestLevels = levelSets{bestLevelIdx};
bestComboC = register_all_curves(curve, defaultRefIdx, defaultU, bestLevels, char(bestInterp));
bestComboMetrics = evaluate_reconstruction(gapList, bestComboC, defaultRefIdx, 'peak', ...
    "refined_dense_equal_level", comboSweepTable.method_label(bestComboIdx));
writetable(bestComboMetrics.summaryTable, fullfile(outDir, 'best_refined_by_gap.csv'));
writetable(bestComboMetrics.globalTable, fullfile(outDir, 'best_refined_global.csv'));
export_surface_table(defaultU, gapList, bestComboC, defaultRefIdx, ...
    bestComboMetrics.CrefRecon, fullfile(outDir, 'best_refined_surfaces.csv'));

%% Experiment 5: amplitude reconstruction strategy
ampStrategies = ["peak"; "top1pct"; "top5pct"; "least_squares_to_ref"; "template_only"];
ampLabels = [ ...
    "Peak amplitude replacement"; ...
    "Top 1% mean amplitude"; ...
    "Top 5% mean amplitude"; ...
    "Least-squares scale to reference"; ...
    "Template-only lower bound"];
ampTables = cell(numel(ampStrategies), 1);
for ia = 1:numel(ampStrategies)
    ampTables{ia} = evaluate_reconstruction(gapList, baselineC, defaultRefIdx, char(ampStrategies(ia)), ...
        sprintf("amp_%s", ampStrategies(ia)), ampLabels(ia)).globalTable;
end
ampSweepTable = vertcat(ampTables{:});
writetable(ampSweepTable, fullfile(outDir, 'amplitude_strategy_sweep.csv'));

%% Experiment 6: small-noise robustness check
noiseLevels = [0, 1e-5, 5e-5, 1e-4, 2e-4]';
nMonteCarlo = 25;
noiseRows(numel(noiseLevels) * 2, 1) = struct( ...
    'method_key', "", ...
    'method_label', "", ...
    'noise_sigma_pF', NaN, ...
    'mean_rmse_pF', NaN, ...
    'std_rmse_pF', NaN, ...
    'mean_max_rmse_pF', NaN, ...
    'std_max_rmse_pF', NaN);
rowId = 0;
refinedLevels = bestLevels;
refinedInterp = char(bestInterp);
rng(7);
for in = 1:numel(noiseLevels)
    sigmaNoise = noiseLevels(in);
    baselineMean = zeros(nMonteCarlo, 1);
    baselineMax = zeros(nMonteCarlo, 1);
    refinedMean = zeros(nMonteCarlo, 1);
    refinedMax = zeros(nMonteCarlo, 1);

    for imc = 1:nMonteCarlo
        noisyCurve = add_noise_and_reprocess(curve, sigmaNoise);
        Cbase = register_all_curves(noisyCurve, defaultRefIdx, defaultU, baselineLevels, defaultInterp);
        Mbase = evaluate_reconstruction(gapList, Cbase, defaultRefIdx, 'peak', ...
            "noise_baseline", "Current dense equal-level").globalTable;
        Crefined = register_all_curves(noisyCurve, defaultRefIdx, defaultU, refinedLevels, refinedInterp);
        Mrefined = evaluate_reconstruction(gapList, Crefined, defaultRefIdx, 'peak', ...
            "noise_refined", "Refined dense equal-level").globalTable;

        baselineMean(imc) = Mbase.mean_reconstructed_rmse_to_ref;
        baselineMax(imc) = Mbase.max_reconstructed_rmse_to_ref;
        refinedMean(imc) = Mrefined.mean_reconstructed_rmse_to_ref;
        refinedMax(imc) = Mrefined.max_reconstructed_rmse_to_ref;
    end

    rowId = rowId + 1;
    noiseRows(rowId).method_key = "baseline_current";
    noiseRows(rowId).method_label = "Current dense equal-level";
    noiseRows(rowId).noise_sigma_pF = sigmaNoise;
    noiseRows(rowId).mean_rmse_pF = mean(baselineMean);
    noiseRows(rowId).std_rmse_pF = std(baselineMean);
    noiseRows(rowId).mean_max_rmse_pF = mean(baselineMax);
    noiseRows(rowId).std_max_rmse_pF = std(baselineMax);

    rowId = rowId + 1;
    noiseRows(rowId).method_key = "refined_dense_equal_level";
    noiseRows(rowId).method_label = sprintf('Refined, %s + %s', levelNames(1), refinedInterp);
    noiseRows(rowId).noise_sigma_pF = sigmaNoise;
    noiseRows(rowId).mean_rmse_pF = mean(refinedMean);
    noiseRows(rowId).std_rmse_pF = std(refinedMean);
    noiseRows(rowId).mean_max_rmse_pF = mean(refinedMax);
    noiseRows(rowId).std_max_rmse_pF = std(refinedMax);
end
noiseRobustnessTable = struct2table(noiseRows);
writetable(noiseRobustnessTable, fullfile(outDir, 'noise_robustness_check.csv'));

%% Aggregate compact best-candidate table
overviewTable = [ ...
    baselineMetrics.globalTable; ...
    best_ref_row(refSweepTable, "reference_gap_sweep", "Best reference gap by absolute RMSE"); ...
    best_row(levelSweepTable, "level_set_sweep", "Best level setting"); ...
    best_row(interpSweepTable, "interpolation_sweep", "Best interpolation"); ...
    best_row(comboSweepTable, "combined_sweep", "Best combined setting"); ...
    best_row(ampSweepTable, "amplitude_strategy_sweep", "Best amplitude strategy")];
writetable(overviewTable, fullfile(outDir, 'refinement_overview.csv'));

%% Figures
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 10]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(categorical(string(refSweepTable.reference_gap_mm)), refSweepTable.mean_reconstructed_rmse_to_ref, ...
    0.65, 'FaceColor', [0.58 0.75 0.88], 'EdgeColor', 'k', 'LineWidth', 0.5);
xlabel('Reference gap g_0 (mm)');
ylabel('Mean reconstructed RMSE (pF)');
title('(a) Reference-gap sweep');

nexttile;
bar(categorical(levelSweepTable.method_label, levelSweepTable.method_label, 'Ordinal', true), ...
    levelSweepTable.mean_reconstructed_rmse_to_ref, ...
    0.65, 'FaceColor', [0.74 0.86 0.66], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean reconstructed RMSE (pF)');
title('(b) Equal-level setting');
xtickangle(25);

nexttile;
bar(categorical(interpSweepTable.method_label, interpSweepTable.method_label, 'Ordinal', true), ...
    interpSweepTable.mean_reconstructed_rmse_to_ref, ...
    0.65, 'FaceColor', [0.95 0.78 0.50], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean reconstructed RMSE (pF)');
title('(c) Interpolation scheme');
xtickangle(25);

nexttile;
bar(categorical(ampSweepTable.method_label, ampSweepTable.method_label, 'Ordinal', true), ...
    ampSweepTable.mean_reconstructed_rmse_to_ref, ...
    0.65, 'FaceColor', [0.88 0.68 0.70], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean reconstructed RMSE (pF)');
title('(d) Amplitude strategy');
xtickangle(25);

title(tl1, sprintf('Refinement checks for no-vibration reconstruction, default g_0 = %.1f mm', ...
    gapList(defaultRefIdx)), 'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_refinement_sweeps.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_refinement_sweeps.pdf'), 'ContentType', 'vector');

fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

[~, hardIdx] = max(baselineMetrics.summaryTable.reconstructed_rmse_to_ref);
Cref = baselineC(defaultRefIdx, :).';
baselineRec = reconstruct_curves(baselineC, defaultRefIdx, 'peak');
lsRec = reconstruct_curves(baselineC, defaultRefIdx, 'least_squares_to_ref');

nexttile; hold on;
plot(defaultU, Cref, 'k--', 'DisplayName', 'True reference');
plot(defaultU, baselineRec(hardIdx, :).', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Current');
plot(defaultU, lsRec(hardIdx, :).', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Least-squares scale');
xlabel('Common coordinate u (mm)');
ylabel('Capacitance (pF)');
title(sprintf('(a) Hardest gap %.1f mm', gapList(hardIdx)));
legend('Location', 'best', 'FontSize', 7.0);

nexttile; hold on;
plot(defaultU, baselineRec(hardIdx, :).' - Cref, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Current residual');
plot(defaultU, lsRec(hardIdx, :).' - Cref, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Least-squares residual');
yline(0, 'k--');
xlabel('Common coordinate u (mm)');
ylabel('Residual (pF)');
title('(b) Residual comparison');
legend('Location', 'best', 'FontSize', 7.0);

title(tl2, 'Diagnostic amplitude-scaling lower bound', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_current_vs_lsq_residual.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_current_vs_lsq_residual.pdf'), 'ContentType', 'vector');

fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8]);
tiledlayout(fig3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
for methodKey = ["baseline_current", "refined_dense_equal_level"]
    mask = noiseRobustnessTable.method_key == methodKey;
    errorbar(noiseRobustnessTable.noise_sigma_pF(mask), ...
        noiseRobustnessTable.mean_rmse_pF(mask), ...
        noiseRobustnessTable.std_rmse_pF(mask), '-o', 'LineWidth', 1.1, ...
        'DisplayName', noiseRobustnessTable.method_label(find(mask, 1, 'first')));
end
xlabel('Added noise standard deviation (pF)');
ylabel('Mean reconstructed RMSE (pF)');
title('(a) Mean error robustness');
legend('Location', 'northwest', 'FontSize', 7.0);

nexttile; hold on;
for methodKey = ["baseline_current", "refined_dense_equal_level"]
    mask = noiseRobustnessTable.method_key == methodKey;
    errorbar(noiseRobustnessTable.noise_sigma_pF(mask), ...
        noiseRobustnessTable.mean_max_rmse_pF(mask), ...
        noiseRobustnessTable.std_max_rmse_pF(mask), '-o', 'LineWidth', 1.1, ...
        'DisplayName', noiseRobustnessTable.method_label(find(mask, 1, 'first')));
end
xlabel('Added noise standard deviation (pF)');
ylabel('Mean worst-gap RMSE (pF)');
title('(b) Worst-gap robustness');

exportgraphics(fig3, fullfile(outDir, 'fig3_noise_robustness.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_noise_robustness.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== No-vibration reconstruction refinement check ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Default reference gap g0 = %.1f mm\n\n', gapList(defaultRefIdx));
disp(overviewTable(:, {'method_key', 'method_label', ...
    'mean_reconstructed_rmse_to_ref', 'max_reconstructed_rmse_to_ref', ...
    'mean_shape_rmse_to_ref'}));
fprintf('\nBest combined setting:\n');
disp(best_row(comboSweepTable, "combined_sweep", "Best combined setting"));
fprintf('\nNoise robustness check:\n');
disp(noiseRobustnessTable(:, {'method_key', 'noise_sigma_pF', 'mean_rmse_pF', 'std_rmse_pF'}));
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function dataFile = resolve_stacked_curve_file(scriptDir)
    files = dir(fullfile(scriptDir, '*.txt'));
    bestFile = "";
    bestScore = -inf;
    for i = 1:numel(files)
        filePath = fullfile(files(i).folder, files(i).name);
        try
            data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
            if size(data, 2) < 2
                continue;
            end
            data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
            if size(data, 1) < 100
                continue;
            end
            x = data(:, 1);
            nCurve = numel([find(diff(x) < 0); numel(x)]);
            score = nCurve * 1e6 + size(data, 1);
            if nCurve >= 5 && score > bestScore
                bestScore = score;
                bestFile = string(filePath);
            end
        catch
        end
    end
    if strlength(bestFile) == 0
        error('No stacked no-vibration curve file was found in: %s', scriptDir);
    end
    dataFile = char(bestFile);
end

function [gapList, xCell, yCell] = load_stacked_curves(filePath)
    data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
    x = data(:, 1);
    y = data(:, 2);
    breakIdx = [find(diff(x) < 0); numel(x)];
    startIdx = [1; breakIdx(1:end-1) + 1];
    nCurve = numel(breakIdx);
    if nCurve == 7
        gapList = (0.2:0.2:1.4)';
    else
        gapList = (1:nCurve)';
    end
    xCell = cell(nCurve, 1);
    yCell = cell(nCurve, 1);
    for i = 1:nCurve
        idx = startIdx(i):breakIdx(i);
        xCell{i} = x(idx);
        yCell{i} = y(idx);
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function levels = make_level_set(low, high, nLevel)
    levels = linspace(low, high, nLevel)';
end

function uGrid = make_reference_grid(refCurve, levels, nU)
    [xLeft, xRight] = level_crossings_normalized(refCurve.xCenter, refCurve.yNorm, levels);
    valid = isfinite(xLeft) & isfinite(xRight);
    if nnz(valid) < 20
        error('Not enough valid anchors for the reference gap %.3f mm.', refCurve.gap);
    end
    uMin = xLeft(find(valid, 1, 'first'));
    uMax = xRight(find(valid, 1, 'first'));
    uGrid = linspace(uMin, uMax, nU)';
end

function Cmat = register_all_curves(curve, refIdx, uGrid, levels, interpMethod)
    nGap = numel(curve);
    Cmat = nan(nGap, numel(uGrid));
    refCurve = curve(refIdx);
    for ig = 1:nGap
        Cmat(ig, :) = register_by_equal_level(curve(ig), refCurve, uGrid, levels, interpMethod).';
    end
end

function yOnU = register_by_equal_level(cur, ref, uGrid, levels, interpMethod)
    [xLeftCur, xRightCur] = level_crossings_normalized(cur.xCenter, cur.yNorm, levels);
    [xLeftRef, xRightRef] = level_crossings_normalized(ref.xCenter, ref.yNorm, levels);
    valid = isfinite(xLeftCur) & isfinite(xRightCur) & ...
            isfinite(xLeftRef) & isfinite(xRightRef);
    if nnz(valid) < max(20, ceil(0.35 * numel(levels)))
        yOnU = interp1(cur.xCenter, cur.yBc, uGrid, interpMethod, 'extrap');
        yOnU = yOnU(:);
        return;
    end

    xLeftAnchor = [xLeftCur(valid); 0];
    uLeftAnchor = [xLeftRef(valid); 0];
    [xLeftAnchor, leftIdx] = unique(xLeftAnchor, 'stable');
    uLeftAnchor = uLeftAnchor(leftIdx);

    xRightAnchor = [0; flipud(xRightCur(valid))];
    uRightAnchor = [0; flipud(xRightRef(valid))];
    [xRightAnchor, rightIdx] = unique(xRightAnchor, 'stable');
    uRightAnchor = uRightAnchor(rightIdx);

    xCenter = cur.xCenter(:);
    yBc = cur.yBc(:);
    [~, iPeak] = max(yBc);

    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), interpMethod, 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), interpMethod, 'extrap');
    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];

    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    [uSorted, uniqIdx] = unique(uSorted, 'stable');
    ySorted = ySorted(uniqIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, interpMethod, 'extrap');
    yOnU = yOnU(:);
end

function [xLeft, xRight] = level_crossings_normalized(x, yNorm, levels)
    x = x(:);
    yNorm = yNorm(:);
    levels = levels(:);
    [~, iPeak] = max(yNorm);
    xRise = x(1:iPeak);
    yRise = yNorm(1:iPeak);
    xFall = x(iPeak:end);
    yFall = yNorm(iPeak:end);
    [yRiseU, idxRiseU] = unique(yRise, 'stable');
    xRiseU = xRise(idxRiseU);
    [yFallU, idxFallU] = unique(yFall, 'stable');
    xFallU = xFall(idxFallU);
    xLeft = nan(size(levels));
    xRight = nan(size(levels));
    for k = 1:numel(levels)
        lv = levels(k);
        if lv >= min(yRiseU) && lv <= max(yRiseU) && lv >= min(yFallU) && lv <= max(yFallU)
            xLeft(k) = interp1(yRiseU, xRiseU, lv, 'linear');
            xRight(k) = interp1(flipud(yFallU), flipud(xFallU), lv, 'linear');
        end
    end
end

function metrics = evaluate_reconstruction(gapList, Cmat, refIdx, ampMode, methodKey, methodLabel)
    Cref = Cmat(refIdx, :).';
    CrefRecon = reconstruct_curves(Cmat, refIdx, ampMode);
    Ameas = max(Cmat, [], 2);
    Hmat = Cmat ./ max(Ameas, eps);
    Href = Hmat(refIdx, :).';

    nGap = numel(gapList);
    rows(nGap, 1) = struct( ...
        'gap_mm', NaN, ...
        'registered_rmse_to_ref', NaN, ...
        'registered_corr_to_ref', NaN, ...
        'shape_rmse_to_ref', NaN, ...
        'shape_corr_to_ref', NaN, ...
        'reconstructed_rmse_to_ref', NaN, ...
        'reconstructed_corr_to_ref', NaN, ...
        'reconstructed_nrmse_to_ref', NaN);

    ampRef = max(abs(Cref), [], 'omitnan');
    for ig = 1:nGap
        regWave = Cmat(ig, :).';
        recWave = CrefRecon(ig, :).';
        shapeWave = Hmat(ig, :).';

        rows(ig).gap_mm = gapList(ig);
        rows(ig).registered_rmse_to_ref = sqrt(mean((regWave - Cref).^2, 'omitnan'));
        rows(ig).registered_corr_to_ref = corr(regWave, Cref, 'Rows', 'complete');
        rows(ig).shape_rmse_to_ref = sqrt(mean((shapeWave - Href).^2, 'omitnan'));
        rows(ig).shape_corr_to_ref = corr(shapeWave, Href, 'Rows', 'complete');
        rows(ig).reconstructed_rmse_to_ref = sqrt(mean((recWave - Cref).^2, 'omitnan'));
        rows(ig).reconstructed_corr_to_ref = corr(recWave, Cref, 'Rows', 'complete');
        rows(ig).reconstructed_nrmse_to_ref = rows(ig).reconstructed_rmse_to_ref / max(ampRef, eps);
    end

    summaryTable = struct2table(rows);
    globalTable = table( ...
        string(methodKey), string(methodLabel), gapList(refIdx), ...
        mean(summaryTable.registered_rmse_to_ref), ...
        max(summaryTable.registered_rmse_to_ref), ...
        mean(summaryTable.shape_rmse_to_ref), ...
        max(summaryTable.shape_rmse_to_ref), ...
        mean(summaryTable.reconstructed_rmse_to_ref), ...
        max(summaryTable.reconstructed_rmse_to_ref), ...
        mean(summaryTable.reconstructed_nrmse_to_ref), ...
        max(summaryTable.reconstructed_nrmse_to_ref), ...
        'VariableNames', {'method_key', 'method_label', 'reference_gap_mm', ...
        'mean_registered_rmse_to_ref', 'max_registered_rmse_to_ref', ...
        'mean_shape_rmse_to_ref', 'max_shape_rmse_to_ref', ...
        'mean_reconstructed_rmse_to_ref', 'max_reconstructed_rmse_to_ref', ...
        'mean_reconstructed_nrmse_to_ref', 'max_reconstructed_nrmse_to_ref'});

    metrics.summaryTable = summaryTable;
    metrics.globalTable = globalTable;
    metrics.CrefRecon = CrefRecon;
end

function CrefRecon = reconstruct_curves(Cmat, refIdx, ampMode)
    nGap = size(Cmat, 1);
    Cref = Cmat(refIdx, :).';
    CrefRecon = nan(size(Cmat));

    switch char(ampMode)
        case 'peak'
            A = max(Cmat, [], 2);
            Aref = A(refIdx);
            for ig = 1:nGap
                CrefRecon(ig, :) = (Aref / max(A(ig), eps)) * Cmat(ig, :);
            end

        case 'top1pct'
            A = top_fraction_mean(Cmat, 0.01);
            Aref = A(refIdx);
            for ig = 1:nGap
                CrefRecon(ig, :) = (Aref / max(A(ig), eps)) * Cmat(ig, :);
            end

        case 'top5pct'
            A = top_fraction_mean(Cmat, 0.05);
            Aref = A(refIdx);
            for ig = 1:nGap
                CrefRecon(ig, :) = (Aref / max(A(ig), eps)) * Cmat(ig, :);
            end

        case 'least_squares_to_ref'
            for ig = 1:nGap
                c = Cmat(ig, :).';
                scale = (c' * Cref) / max(c' * c, eps);
                CrefRecon(ig, :) = (scale * c).';
            end

        case 'template_only'
            for ig = 1:nGap
                CrefRecon(ig, :) = Cref.';
            end

        otherwise
            error('Unknown amplitude mode: %s', ampMode);
    end
end

function A = top_fraction_mean(Cmat, frac)
    nGap = size(Cmat, 1);
    A = zeros(nGap, 1);
    nTop = max(1, round(frac * size(Cmat, 2)));
    for ig = 1:nGap
        values = sort(Cmat(ig, :), 'descend');
        A(ig) = mean(values(1:nTop), 'omitnan');
    end
end

function row = best_row(tableIn, methodKey, methodLabel)
    [~, idx] = min(tableIn.mean_reconstructed_rmse_to_ref);
    row = tableIn(idx, :);
    row.method_key = string(methodKey);
    row.method_label = string(methodLabel);
end

function row = best_ref_row(tableIn, methodKey, methodLabel)
    [~, idx] = min(tableIn.mean_reconstructed_rmse_to_ref);
    row = table( ...
        string(methodKey), string(methodLabel), tableIn.reference_gap_mm(idx), ...
        NaN, NaN, ...
        tableIn.mean_shape_rmse_to_ref(idx), tableIn.max_shape_rmse_to_ref(idx), ...
        tableIn.mean_reconstructed_rmse_to_ref(idx), tableIn.max_reconstructed_rmse_to_ref(idx), ...
        tableIn.mean_reconstructed_nrmse_to_ref(idx), tableIn.max_reconstructed_nrmse_to_ref(idx), ...
        'VariableNames', {'method_key', 'method_label', 'reference_gap_mm', ...
        'mean_registered_rmse_to_ref', 'max_registered_rmse_to_ref', ...
        'mean_shape_rmse_to_ref', 'max_shape_rmse_to_ref', ...
        'mean_reconstructed_rmse_to_ref', 'max_reconstructed_rmse_to_ref', ...
        'mean_reconstructed_nrmse_to_ref', 'max_reconstructed_nrmse_to_ref'});
end

function [levelIdx, interpName] = parse_combo_key(methodKey)
    token = regexp(char(methodKey), '^combo_(\d+)_(.+)$', 'tokens', 'once');
    if isempty(token)
        error('Unexpected combo method key: %s', methodKey);
    end
    levelIdx = str2double(token{1});
    interpName = string(token{2});
end

function noisyCurve = add_noise_and_reprocess(curve, sigmaNoise)
    noisyCurve = curve;
    for ig = 1:numel(curve)
        yNoisy = curve(ig).yRaw + sigmaNoise * randn(size(curve(ig).yRaw));
        [base, amp, yNorm] = normalize_waveform(yNoisy);
        yBc = yNoisy - base;
        [~, iPeak] = max(yBc);
        xPeak = curve(ig).x(iPeak);
        noisyCurve(ig).yRaw = yNoisy;
        noisyCurve(ig).yBc = yBc;
        noisyCurve(ig).base = base;
        noisyCurve(ig).amp = amp;
        noisyCurve(ig).xPeak = xPeak;
        noisyCurve(ig).xCenter = curve(ig).x - xPeak;
        noisyCurve(ig).yNorm = yNorm;
    end
end

function export_surface_table(uGrid, gapList, Cmat, refIdx, CrefRecon, filePath)
    surfaceTable = table();
    surfaceTable.u_mm = uGrid(:);
    Ameas = max(Cmat, [], 2);
    Hmat = Cmat ./ max(Ameas, eps);
    surfaceTable.C_ref = Cmat(refIdx, :).';
    surfaceTable.H_ref = Hmat(refIdx, :).';
    for ig = 1:numel(gapList)
        suffix = sprintf('%03d', round(10 * gapList(ig)));
        surfaceTable.(sprintf('C_g_%s', suffix)) = Cmat(ig, :).';
        surfaceTable.(sprintf('H_g_%s', suffix)) = Hmat(ig, :).';
        surfaceTable.(sprintf('Cref_from_g_%s', suffix)) = CrefRecon(ig, :).';
    end
    writetable(surfaceTable, filePath);
end
