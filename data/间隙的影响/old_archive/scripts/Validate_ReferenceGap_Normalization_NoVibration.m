%% Validate reference-gap normalization using no-vibration waveforms
% We build:
%   C(u,g)   : baseline-corrected no-vibration waveform surface
%   C_u(u,g) : derivative of C with respect to the common coordinate u
%
% The common coordinate u is constructed by equal-level branch registration
% to a chosen reference gap g0. Then, for each target gap g_i, we test
% whether the no-vibration waveform can be normalized to the reference gap
% using the proposed first-order formula:
%
%   y_ref(u) = C(u,g0) + [C_u(u,g0) / C_u(u,g_i)] * ( y(u) - C(u,g_i) )
%
% For an exact surface and no-vibration input, the residual term vanishes
% and y_ref(u) should collapse to the reference template C(u,g0). To avoid
% a trivial identity test, a leave-one-gap-out interpolation in g is used
% when predicting C(u,g_i) and C_u(u,g_i).

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'reference_gap_normalization_results');
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
refGapTarget = 0.8;            % mm
levelMap = (0.02:0.02:0.98)';  % dense equal-level mapping anchors
nU = 1401;                     % common coordinate samples
epsRel = 0.05;                 % derivative regularization strength
slopeGate = 0.10;              % diagnostic gate for reporting Delta u

%% Load stacked no-vibration waveforms
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = min(abs(gapList - refGapTarget));
g0 = gapList(refIdx);

%% Preprocess waveforms and extract equal-level branches
curve = repmat(struct( ...
    'gap', 0, 'x', [], 'yRaw', [], 'yBc', [], ...
    'base', 0, 'amp', 0, 'xPeak', 0, ...
    'xCenter', [], 'yNorm', [], ...
    'xLeftMap', [], 'xRightMap', []), nGap, 1);

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

%% Build common reference coordinate u from the reference gap
refCurve = curve(refIdx);
validRef = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
uMin = refCurve.xLeftMap(find(validRef, 1, 'first'));
uMax = refCurve.xRightMap(find(validRef, 1, 'first'));
uGrid = linspace(uMin, uMax, nU)';

%% Map every no-vibration waveform into the common u-domain
Cmat = nan(nGap, nU);
phiCell = cell(nGap, 1);
mapRows(nGap, 1) = struct( ...
    'gap_mm', NaN, ...
    'common_domain_rmse_to_ref', NaN, ...
    'common_domain_corr_to_ref', NaN);

for ig = 1:nGap
    [uOfX, yOnU] = map_waveform_to_reference_domain(curve(ig), refCurve, uGrid);
    phiCell{ig} = uOfX;
    Cmat(ig, :) = yOnU(:)';
end

%% Surface derivative C_u(u,g)
CuMat = nan(size(Cmat));
for ig = 1:nGap
    CuMat(ig, :) = gradient(Cmat(ig, :), uGrid);
end

%% Direct common-domain consistency against the reference template
CrefTrue = Cmat(refIdx, :).';
for ig = 1:nGap
    cRow = Cmat(ig, :).';
    mapRows(ig).gap_mm = gapList(ig);
    mapRows(ig).common_domain_rmse_to_ref = sqrt(mean((cRow - CrefTrue).^2, 'omitnan'));
    mapRows(ig).common_domain_corr_to_ref = corr(cRow, CrefTrue, 'Rows', 'complete');
end
mapSummaryTable = struct2table(mapRows);
writetable(mapSummaryTable, fullfile(outDir, 'common_domain_consistency_summary.csv'));

%% Validate reference-gap normalization with leave-one-gap-out interpolation
normRows(nGap, 1) = struct( ...
    'gap_mm', NaN, ...
    'direct_common_domain_rmse_to_ref', NaN, ...
    'loo_fit_rmse_at_g', NaN, ...
    'loo_normalized_rmse_to_ref', NaN, ...
    'loo_delta_u_rms_mm', NaN, ...
    'loo_delta_u_max_mm', NaN);

repTargetGaps = [gapList(1), g0, gapList(end)];
repMask = ismembertol(gapList, repTargetGaps, 1e-9);
repRows = [];

for ig = 1:nGap
    g = gapList(ig);
    targetWave = Cmat(ig, :).';

    [CpredG, CupredG] = interp_surface_leave_one_out(gapList, Cmat, CuMat, ig, g);
    [CpredRef, CupredRef] = interp_surface_leave_one_out(gapList, Cmat, CuMat, ig, g0);

    residual = targetWave - CpredG;
    slopeScale = max(abs(CupredG), [], 'omitnan');
    denom = CupredG.^2 + (epsRel * slopeScale)^2;
    deltaURaw = (CupredG ./ max(denom, eps)) .* residual;
    validSlope = abs(CupredG) >= slopeGate * slopeScale;
    deltaU = deltaURaw;
    deltaU(~validSlope) = NaN;
    gain = (CupredRef .* CupredG) ./ max(denom, eps);
    yRefNorm = CpredRef + gain .* residual;

    normRows(ig).gap_mm = g;
    normRows(ig).direct_common_domain_rmse_to_ref = mapSummaryTable.common_domain_rmse_to_ref(ig);
    normRows(ig).loo_fit_rmse_at_g = sqrt(mean((targetWave - CpredG).^2, 'omitnan'));
    normRows(ig).loo_normalized_rmse_to_ref = sqrt(mean((yRefNorm - CrefTrue).^2, 'omitnan'));
    normRows(ig).loo_delta_u_rms_mm = sqrt(mean(deltaU.^2, 'omitnan'));
    normRows(ig).loo_delta_u_max_mm = max(abs(deltaU), [], 'omitnan');

    if repMask(ig)
        row.gap_mm = g;
        row.targetWave = targetWave;
        row.predG = CpredG;
        row.refTrue = CrefTrue;
        row.predRef = CpredRef;
        row.yRefNorm = yRefNorm;
        row.deltaU = deltaU;
        repRows = [repRows; row]; %#ok<AGROW>
    end
end

normSummaryTable = struct2table(normRows);
writetable(normSummaryTable, fullfile(outDir, 'reference_gap_normalization_summary.csv'));

globalSummary = table( ...
    g0, ...
    mean(mapSummaryTable.common_domain_rmse_to_ref), ...
    max(mapSummaryTable.common_domain_rmse_to_ref), ...
    mean(normSummaryTable.loo_fit_rmse_at_g), ...
    max(normSummaryTable.loo_fit_rmse_at_g), ...
    mean(normSummaryTable.loo_normalized_rmse_to_ref), ...
    max(normSummaryTable.loo_normalized_rmse_to_ref), ...
    mean(normSummaryTable.loo_delta_u_rms_mm), ...
    max(normSummaryTable.loo_delta_u_max_mm), ...
    'VariableNames', {'reference_gap_mm', ...
    'mean_common_domain_rmse_to_ref', 'max_common_domain_rmse_to_ref', ...
    'mean_loo_fit_rmse_at_g', 'max_loo_fit_rmse_at_g', ...
    'mean_loo_normalized_rmse_to_ref', 'max_loo_normalized_rmse_to_ref', ...
    'mean_loo_delta_u_rms_mm', 'max_loo_delta_u_max_mm'});
writetable(globalSummary, fullfile(outDir, 'reference_gap_normalization_global_summary.csv'));

%% Save the calibrated surfaces
surfaceTable = table();
surfaceTable.u_mm = uGrid;
surfaceTable.C_ref_gap = CrefTrue;
surfaceTable.Cu_ref_gap = CuMat(refIdx, :).';
for ig = 1:nGap
    surfaceTable.(sprintf('C_g_%03d', round(10 * gapList(ig)))) = Cmat(ig, :).';
    surfaceTable.(sprintf('Cu_g_%03d', round(10 * gapList(ig)))) = CuMat(ig, :).';
end
writetable(surfaceTable, fullfile(outDir, 'calibrated_surface_C_and_Cu.csv'));

%% Figure 1: transformed no-vibration waveforms in the common domain
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
    plot(uGrid, Cmat(ig, :).', '-', 'Color', colors(ig, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(ig)));
end
xlabel('Common coordinate u (mm)');
ylabel('C(u,g) (pF)');
title('(b) Registered no-vibration waveforms');

nexttile;
imagesc(gapList, uGrid, Cmat.');
set(gca, 'YDir', 'normal');
xlabel('Gap g (mm)');
ylabel('Common coordinate u (mm)');
title('(c) Calibrated surface C(u,g)');
cb = colorbar; cb.Label.String = 'Capacitance (pF)';

nexttile;
imagesc(gapList, uGrid, CuMat.');
set(gca, 'YDir', 'normal');
xlabel('Gap g (mm)');
ylabel('Common coordinate u (mm)');
title('(d) Derivative surface C_u(u,g)');
cb = colorbar; cb.Label.String = 'dC/du';

title(tl1, sprintf('Reference-domain calibration with g_0 = %.1f mm', g0), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_surface_calibration.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_surface_calibration.pdf'), 'ContentType', 'vector');

%% Figure 2: direct common-domain consistency and normalization summary
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
gapCats = categorical(string(gapList));

nexttile;
bar(gapCats, [mapSummaryTable.common_domain_rmse_to_ref, normSummaryTable.loo_normalized_rmse_to_ref], 'grouped');
xlabel('Gap g (mm)');
ylabel('RMSE to reference template (pF)');
title('(a) Before and after reference-gap normalization');
legend({'Registered waveform', 'LOO normalized waveform'}, 'Location', 'northwest', 'FontSize', 7.0);

nexttile;
yyaxis left;
bar(gapCats, normSummaryTable.loo_fit_rmse_at_g, 0.5, 'FaceColor', [0.55 0.69 0.89], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('LOO fit RMSE at g (pF)');
yyaxis right;
plot(gapCats, normSummaryTable.loo_delta_u_rms_mm, '-o', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10]);
ylabel('Equivalent \Delta u RMS (mm)');
title('(b) Residual and equivalent coordinate disturbance');

exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_summary.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_summary.pdf'), 'ContentType', 'vector');

%% Figure 3: representative leave-one-out normalization examples
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 15]);
nRep = numel(repRows);
tl3 = tiledlayout(fig3, nRep, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for ir = 1:nRep
    nexttile; hold on;
    plot(uGrid, repRows(ir).targetWave, 'Color', [0.10 0.45 0.75], 'DisplayName', sprintf('Target g = %.1f', repRows(ir).gap_mm));
    plot(uGrid, repRows(ir).predG, 'k--', 'DisplayName', 'Predicted C(u,g)');
    plot(uGrid, repRows(ir).refTrue, '-', 'Color', [0.30 0.30 0.30], 'DisplayName', sprintf('True ref g_0 = %.1f', g0));
    plot(uGrid, repRows(ir).yRefNorm, '-', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Normalized to g_0');
    xlabel('Common coordinate u (mm)');
    ylabel('Capacitance (pF)');
    title(sprintf('(%c) Gap %.1f mm: normalization result', char('a' + 2*(ir-1)), repRows(ir).gap_mm));
    legend('Location', 'best', 'FontSize', 6.8);

    nexttile; hold on;
    plot(uGrid, repRows(ir).targetWave - repRows(ir).predG, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Residual at g');
    plot(uGrid, repRows(ir).deltaU, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Equivalent \Delta u');
    yline(0, 'k--');
    xlabel('Common coordinate u (mm)');
    ylabel('Residual / \Delta u');
    title(sprintf('(%c) Gap %.1f mm: residual and \x0394u', char('b' + 2*(ir-1)), repRows(ir).gap_mm));
    legend('Location', 'best', 'FontSize', 6.8);
end

exportgraphics(fig3, fullfile(outDir, 'fig3_reference_gap_examples.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_reference_gap_examples.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Reference-gap normalization validation (no vibration) ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Reference gap g0 = %.1f mm\n\n', g0);
disp(globalSummary);
fprintf('\nCommon-domain consistency:\n');
disp(mapSummaryTable);
fprintf('\nLeave-one-gap-out normalization summary:\n');
disp(normSummaryTable);
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

    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), 'pchip', 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), 'pchip', 'extrap');

    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];

    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, 'pchip', 'extrap');
end

function [Cpred, Cupred] = interp_surface_leave_one_out(gaps, Cmat, CuMat, leaveIdx, gTarget)
    keep = true(size(gaps));
    keep(leaveIdx) = false;
    gapsKeep = gaps(keep);

    nU = size(Cmat, 2);
    Cpred = nan(nU, 1);
    Cupred = nan(nU, 1);
    for iu = 1:nU
        Cpred(iu) = interp1(gapsKeep, Cmat(keep, iu), gTarget, 'pchip', 'extrap');
        Cupred(iu) = interp1(gapsKeep, CuMat(keep, iu), gTarget, 'pchip', 'extrap');
    end
end
