%% No-vibration reference-gap reconstruction
% Processing logic for no-vibration waveforms only:
%   1) subtract baseline;
%   2) register all waveforms to a common coordinate u using equal-level
%      left/right branch mapping to the largest-gap reference;
%   3) decompose C(u,g) = A(g) * H(u,g);
%   4) remove gap influence by rebuilding each waveform at the reference
%      gap amplitude while preserving its registered normalized shape:
%
%         C_ref_from_g(u) = A(g0) * H(u,g)
%
% The effect is evaluated by how close C_ref_from_g(u) is to the true
% reference-gap waveform C(u,g0).

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    dataFile = resolve_data_file(scriptDir);
end
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'no_vibration_reference_gap_reconstruction_results');
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

%% Settings
levelMap = linspace(0.005, 0.995, 199)';
nU = 1401;

%% Load no-vibration waveforms
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);
gRef = gapList(refIdx);

%% Preprocess and extract level anchors
curve = repmat(struct( ...
    'gap', 0, 'x', [], 'yRaw', [], 'yBc', [], ...
    'base', 0, 'amp', 0, 'xPeak', 0, 'xCenter', [], ...
    'yNorm', [], 'xLeftMap', [], 'xRightMap', []), nGap, 1);

for ig = 1:nGap
    x = xCell{ig}(:);
    y = yCell{ig}(:);
    [base, amp, yNorm] = normalize_waveform(y);
    yBc = y - base;
    [~, iPeak] = max(yBc);
    xPeak = x(iPeak);
    xCenter = x - xPeak;

    [xLeftMap, xRightMap] = level_crossings_normalized(xCenter, yNorm, levelMap);
    validMap = isfinite(xLeftMap) & isfinite(xRightMap);
    if nnz(validMap) < 20
        error('Not enough valid equal-level anchors for gap %.3f mm.', gapList(ig));
    end

    curve(ig).gap = gapList(ig);
    curve(ig).x = x;
    curve(ig).yRaw = y;
    curve(ig).yBc = yBc;
    curve(ig).base = base;
    curve(ig).amp = amp;
    curve(ig).xPeak = xPeak;
    curve(ig).xCenter = xCenter;
    curve(ig).yNorm = yNorm;
    curve(ig).xLeftMap = xLeftMap;
    curve(ig).xRightMap = xRightMap;
end

%% Build common coordinate u from the largest-gap reference
refCurve = curve(refIdx);
validRef = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
uMin = refCurve.xLeftMap(find(validRef, 1, 'first'));
uMax = refCurve.xRightMap(find(validRef, 1, 'first'));
uGrid = linspace(uMin, uMax, nU)';

%% Map all waveforms into common domain
Cmat = nan(nGap, nU);
for ig = 1:nGap
    [~, yOnU] = map_waveform_to_reference_domain(curve(ig), refCurve, uGrid);
    Cmat(ig, :) = yOnU(:)';
end

%% Amplitude-shape decomposition
Ameas = max(Cmat, [], 2);
Hmat = Cmat ./ max(Ameas, eps);

Cref = Cmat(refIdx, :).';
Aref = Ameas(refIdx);
Href = Hmat(refIdx, :).';

%% Reconstruct every waveform to the reference gap
CrefRecon = nan(size(Cmat));
for ig = 1:nGap
    CrefRecon(ig, :) = (Aref * Hmat(ig, :).').';
end

%% Metrics
rows(nGap, 1) = struct( ...
    'gap_mm', NaN, ...
    'registered_rmse_to_ref', NaN, ...
    'registered_corr_to_ref', NaN, ...
    'shape_rmse_to_ref', NaN, ...
    'shape_corr_to_ref', NaN, ...
    'reconstructed_rmse_to_ref', NaN, ...
    'reconstructed_corr_to_ref', NaN);

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
end

summaryTable = struct2table(rows);
writetable(summaryTable, fullfile(outDir, 'no_vibration_reference_gap_reconstruction_summary.csv'));

globalSummary = table( ...
    gRef, ...
    mean(summaryTable.registered_rmse_to_ref), ...
    max(summaryTable.registered_rmse_to_ref), ...
    mean(summaryTable.shape_rmse_to_ref), ...
    max(summaryTable.shape_rmse_to_ref), ...
    mean(summaryTable.reconstructed_rmse_to_ref), ...
    max(summaryTable.reconstructed_rmse_to_ref), ...
    'VariableNames', {'reference_gap_mm', ...
    'mean_registered_rmse_to_ref', 'max_registered_rmse_to_ref', ...
    'mean_shape_rmse_to_ref', 'max_shape_rmse_to_ref', ...
    'mean_reconstructed_rmse_to_ref', 'max_reconstructed_rmse_to_ref'});
writetable(globalSummary, fullfile(outDir, 'no_vibration_reference_gap_reconstruction_global.csv'));

%% Export calibrated data
surfaceTable = table();
surfaceTable.u_mm = uGrid;
surfaceTable.C_ref = Cref;
surfaceTable.H_ref = Href;
for ig = 1:nGap
    surfaceTable.(sprintf('C_g_%03d', round(10 * gapList(ig)))) = Cmat(ig, :).';
    surfaceTable.(sprintf('H_g_%03d', round(10 * gapList(ig)))) = Hmat(ig, :).';
    surfaceTable.(sprintf('Cref_from_g_%03d', round(10 * gapList(ig)))) = CrefRecon(ig, :).';
end
writetable(surfaceTable, fullfile(outDir, 'calibrated_registered_and_reconstructed_surfaces.csv'));

%% Figure 1: original, registered, normalized shapes, reconstructed
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(nGap);

nexttile; hold on;
for ig = 1:nGap
    plot(curve(ig).xCenter, curve(ig).yBc, '-', 'Color', colors(ig, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(ig)));
end
xlabel('Peak-centered coordinate x - x_p (mm)');
ylabel('Baseline-corrected capacitance (pF)');
title('(a) Original no-vibration waveforms');
legend('Location', 'eastoutside', 'FontSize', 7.0);

nexttile; hold on;
for ig = 1:nGap
    plot(uGrid, Cmat(ig, :).', '-', 'Color', colors(ig, :));
end
xlabel('Common coordinate u (mm)');
ylabel('Registered waveform C(u,g) (pF)');
title('(b) Registered no-vibration waveforms');

nexttile; hold on;
for ig = 1:nGap
    plot(uGrid, Hmat(ig, :).', '-', 'Color', colors(ig, :));
end
xlabel('Common coordinate u (mm)');
ylabel('Normalized shape H(u,g)');
title('(c) Shape collapse after registration');

nexttile; hold on;
for ig = 1:nGap
    plot(uGrid, CrefRecon(ig, :).', '-', 'Color', colors(ig, :));
end
plot(uGrid, Cref, 'k--', 'LineWidth', 1.0, 'DisplayName', 'True reference');
xlabel('Common coordinate u (mm)');
ylabel('Reconstructed reference waveform (pF)');
title('(d) Reference-gap reconstruction');

title(tl1, sprintf('No-vibration reconstruction to reference gap g_0 = %.1f mm', gRef), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_reconstruction_pipeline.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_reconstruction_pipeline.pdf'), 'ContentType', 'vector');

%% Figure 2: error summary
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
gapCats = categorical(string(gapList));

nexttile;
bar(gapCats, [summaryTable.registered_rmse_to_ref, summaryTable.reconstructed_rmse_to_ref], 'grouped');
xlabel('Gap g (mm)');
ylabel('RMSE to reference waveform (pF)');
title('(a) Before and after reference-gap reconstruction');
legend({'Registered waveform', 'Reconstructed waveform'}, 'Location', 'northwest', 'FontSize', 7.0);

nexttile;
yyaxis left;
bar(gapCats, summaryTable.shape_rmse_to_ref, 0.55, 'FaceColor', [0.55 0.69 0.89], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Shape RMSE to H_{ref}');
yyaxis right;
plot(gapCats, Ameas / Aref, '-o', 'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10]);
ylabel('Amplitude ratio A(g) / A(g_0)');
title('(b) Residual shape spread and amplitude trend');

exportgraphics(fig2, fullfile(outDir, 'fig2_reconstruction_summary.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_reconstruction_summary.pdf'), 'ContentType', 'vector');

%% Figure 3: representative residuals to reference
repTargetGaps = [gapList(1), gapList(round((nGap+1)/2)), gapList(end)];
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 15]);
tl3 = tiledlayout(fig3, numel(repTargetGaps), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for ir = 1:numel(repTargetGaps)
    gNow = repTargetGaps(ir);
    ig = find(abs(gapList - gNow) == min(abs(gapList - gNow)), 1, 'first');
    regWave = Cmat(ig, :).';
    recWave = CrefRecon(ig, :).';

    nexttile; hold on;
    plot(uGrid, regWave, 'Color', [0.10 0.45 0.75], 'DisplayName', sprintf('Registered g = %.1f', gNow));
    plot(uGrid, recWave, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Reconstructed to g_0');
    plot(uGrid, Cref, 'k--', 'DisplayName', sprintf('True reference g_0 = %.1f', gRef));
    xlabel('Common coordinate u (mm)');
    ylabel('Capacitance (pF)');
    title(sprintf('(%c) Gap %.1f mm: waveform comparison', char('a' + 2*(ir-1)), gNow));
    legend('Location', 'best', 'FontSize', 6.8);

    nexttile; hold on;
    plot(uGrid, regWave - Cref, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Registered - reference');
    plot(uGrid, recWave - Cref, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Reconstructed - reference');
    yline(0, 'k--');
    xlabel('Common coordinate u (mm)');
    ylabel('Residual (pF)');
    title(sprintf('(%c) Gap %.1f mm: residuals', char('b' + 2*(ir-1)), gNow));
    legend('Location', 'best', 'FontSize', 6.8);
end

exportgraphics(fig3, fullfile(outDir, 'fig3_representative_residuals.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_representative_residuals.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== No-vibration reference-gap reconstruction ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Reference gap g0 = %.1f mm (largest gap)\n\n', gRef);
disp(globalSummary);
fprintf('\nGap-by-gap summary:\n');
disp(summaryTable);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
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

function [uOfX, yOnU] = map_waveform_to_reference_domain(curve, refCurve, uGrid)
    valid = isfinite(curve.xLeftMap) & isfinite(curve.xRightMap) & ...
            isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);

    xLeftAnchor = [curve.xLeftMap(valid); 0];
    uLeftAnchor = [refCurve.xLeftMap(valid); 0];
    [xLeftAnchor, leftIdx] = unique(xLeftAnchor, 'stable');
    uLeftAnchor = uLeftAnchor(leftIdx);

    xRightAnchor = [0; flipud(curve.xRightMap(valid))];
    uRightAnchor = [0; flipud(refCurve.xRightMap(valid))];
    [xRightAnchor, rightIdx] = unique(xRightAnchor, 'stable');
    uRightAnchor = uRightAnchor(rightIdx);

    xCenter = curve.xCenter(:);
    yBc = curve.yBc(:);
    [~, iPeak] = max(yBc);

    % Dense anchors are sufficiently informative, so linear interpolation
    % avoids the small overshoot seen with spline-like branch mappings.
    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), 'linear', 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), 'linear', 'extrap');

    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];
    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, 'linear', 'extrap');
end

function dataFile = resolve_data_file(scriptDir)
    candidateList = { ...
        fullfile(scriptDir, '直叶片2mm_不同间隙.txt'), ...
        fullfile(scriptDir, '鐩村彾鐗?mm_涓嶅悓闂撮殭.txt')};
    for i = 1:numel(candidateList)
        if isfile(candidateList{i})
            dataFile = candidateList{i};
            return;
        end
    end
    dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
end
