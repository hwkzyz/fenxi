%% Validate reference-gap normalization using amplitude-shape surface model
% Compared with the first attempt, this version is more conservative:
% 1) choose the largest gap as reference, which is the broadest / smoothest;
% 2) construct a common coordinate by equal-level branch registration;
% 3) decompose C(u,g) = A(g) * H(u,g), where A(g) is a scalar amplitude law
%    and H(u,g) is the normalized shape field;
% 4) use leave-one-gap-out prediction for A(g_i) and H(u,g_i) separately;
% 5) then apply the first-order reference-gap normalization formula.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'reference_gap_normalization_v2_results');
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
levelMap = (0.02:0.02:0.98)';   % dense branch anchors
nU = 1401;
epsRel = 0.05;
slopeGate = 0.10;

%% Load data
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);
g0 = gapList(refIdx);

%% Preprocess
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

%% Common coordinate from the largest-gap reference
refCurve = curve(refIdx);
validRef = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
uMin = refCurve.xLeftMap(find(validRef, 1, 'first'));
uMax = refCurve.xRightMap(find(validRef, 1, 'first'));
uGrid = linspace(uMin, uMax, nU)';

%% Register all no-vibration waveforms
Cmat = nan(nGap, nU);
for ig = 1:nGap
    [~, yOnU] = map_waveform_to_reference_domain(curve(ig), refCurve, uGrid);
    Cmat(ig, :) = yOnU(:)';
end

CuMat = nan(size(Cmat));
for ig = 1:nGap
    CuMat(ig, :) = gradient(Cmat(ig, :), uGrid);
end

%% Amplitude-shape decomposition
Ameas = max(Cmat, [], 2);
Hmat = Cmat ./ max(Ameas, eps);
HuMat = nan(size(Hmat));
for ig = 1:nGap
    HuMat(ig, :) = gradient(Hmat(ig, :), uGrid);
end

CrefTrue = Cmat(refIdx, :).';
CuRefTrue = CuMat(refIdx, :).';
HrefTrue = Hmat(refIdx, :).';

%% Direct common-domain consistency
mapRows(nGap, 1) = struct('gap_mm', NaN, 'common_domain_rmse_to_ref', NaN, ...
    'shape_rmse_to_ref', NaN, 'shape_corr_to_ref', NaN);
for ig = 1:nGap
    mapRows(ig).gap_mm = gapList(ig);
    mapRows(ig).common_domain_rmse_to_ref = sqrt(mean((Cmat(ig, :).'-CrefTrue).^2, 'omitnan'));
    mapRows(ig).shape_rmse_to_ref = sqrt(mean((Hmat(ig, :).'-HrefTrue).^2, 'omitnan'));
    mapRows(ig).shape_corr_to_ref = corr(Hmat(ig, :).', HrefTrue, 'Rows', 'complete');
end
mapSummaryTable = struct2table(mapRows);
writetable(mapSummaryTable, fullfile(outDir, 'common_domain_consistency_summary.csv'));

%% Leave-one-gap-out validation with amplitude-shape modeling
normRows(nGap, 1) = struct( ...
    'gap_mm', NaN, ...
    'common_domain_rmse_to_ref', NaN, ...
    'loo_amp_fit_abs_err', NaN, ...
    'loo_shape_fit_rmse', NaN, ...
    'loo_fit_rmse_at_g', NaN, ...
    'loo_normalized_rmse_to_ref', NaN, ...
    'loo_delta_u_rms_mm', NaN, ...
    'loo_delta_u_max_mm', NaN);

repTargetGaps = [gapList(1), gapList(round((nGap+1)/2)), gapList(end)];
repMask = ismembertol(gapList, repTargetGaps, 1e-9);
repRows = [];

for ig = 1:nGap
    g = gapList(ig);
    targetWave = Cmat(ig, :).';

    keep = true(nGap, 1);
    keep(ig) = false;
    gapsKeep = gapList(keep);

    Ahat = predict_amplitude_leave_one_out(gapsKeep, Ameas(keep), g);
    Hhat = nan(nU, 1);
    Huhat = nan(nU, 1);
    for iu = 1:nU
        Hhat(iu) = interp1(gapsKeep, Hmat(keep, iu), g, 'pchip', 'extrap');
        Huhat(iu) = interp1(gapsKeep, HuMat(keep, iu), g, 'pchip', 'extrap');
    end

    CpredG = Ahat * Hhat;
    CupredG = Ahat * Huhat;

    residual = targetWave - CpredG;
    slopeScale = max(abs(CupredG), [], 'omitnan');
    denom = CupredG.^2 + (epsRel * slopeScale)^2;
    deltaURaw = (CupredG ./ max(denom, eps)) .* residual;
    validSlope = abs(CupredG) >= slopeGate * slopeScale;
    deltaU = deltaURaw;
    deltaU(~validSlope) = NaN;
    gain = (CuRefTrue .* CupredG) ./ max(denom, eps);
    yRefNorm = CrefTrue + gain .* residual;

    normRows(ig).gap_mm = g;
    normRows(ig).common_domain_rmse_to_ref = mapSummaryTable.common_domain_rmse_to_ref(ig);
    normRows(ig).loo_amp_fit_abs_err = abs(Ahat - Ameas(ig));
    normRows(ig).loo_shape_fit_rmse = sqrt(mean((Hhat - Hmat(ig, :).').^2, 'omitnan'));
    normRows(ig).loo_fit_rmse_at_g = sqrt(mean((targetWave - CpredG).^2, 'omitnan'));
    normRows(ig).loo_normalized_rmse_to_ref = sqrt(mean((yRefNorm - CrefTrue).^2, 'omitnan'));
    normRows(ig).loo_delta_u_rms_mm = sqrt(mean(deltaU.^2, 'omitnan'));
    normRows(ig).loo_delta_u_max_mm = max(abs(deltaU), [], 'omitnan');

    if repMask(ig)
        row.gap_mm = g;
        row.targetWave = targetWave;
        row.predG = CpredG;
        row.refTrue = CrefTrue;
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
    mean(mapSummaryTable.shape_rmse_to_ref), ...
    max(mapSummaryTable.shape_rmse_to_ref), ...
    mean(normSummaryTable.loo_amp_fit_abs_err), ...
    max(normSummaryTable.loo_amp_fit_abs_err), ...
    mean(normSummaryTable.loo_fit_rmse_at_g), ...
    max(normSummaryTable.loo_fit_rmse_at_g), ...
    mean(normSummaryTable.loo_normalized_rmse_to_ref), ...
    max(normSummaryTable.loo_normalized_rmse_to_ref), ...
    'VariableNames', {'reference_gap_mm', ...
    'mean_common_domain_rmse_to_ref', 'max_common_domain_rmse_to_ref', ...
    'mean_shape_rmse_to_ref', 'max_shape_rmse_to_ref', ...
    'mean_amp_fit_abs_err', 'max_amp_fit_abs_err', ...
    'mean_loo_fit_rmse_at_g', 'max_loo_fit_rmse_at_g', ...
    'mean_loo_normalized_rmse_to_ref', 'max_loo_normalized_rmse_to_ref'});
writetable(globalSummary, fullfile(outDir, 'reference_gap_normalization_global_summary.csv'));

%% Export calibrated surfaces
surfaceTable = table();
surfaceTable.u_mm = uGrid;
surfaceTable.C_ref_gap = CrefTrue;
surfaceTable.Cu_ref_gap = CuRefTrue;
surfaceTable.H_ref_gap = HrefTrue;
for ig = 1:nGap
    surfaceTable.(sprintf('C_g_%03d', round(10 * gapList(ig)))) = Cmat(ig, :).';
    surfaceTable.(sprintf('H_g_%03d', round(10 * gapList(ig)))) = Hmat(ig, :).';
end
writetable(surfaceTable, fullfile(outDir, 'calibrated_surface_C_H_Cu.csv'));

%% Figure 1
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
    plot(uGrid, Hmat(ig, :).', '-', 'Color', colors(ig, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(ig)));
end
xlabel('Common coordinate u (mm)');
ylabel('Normalized shape H(u,g)');
title('(b) Registered normalized shapes');

nexttile;
imagesc(gapList, uGrid, Cmat.');
set(gca, 'YDir', 'normal');
xlabel('Gap g (mm)');
ylabel('Common coordinate u (mm)');
title('(c) Surface C(u,g)');
cb = colorbar; cb.Label.String = 'Capacitance (pF)';

nexttile;
imagesc(gapList, uGrid, Hmat.');
set(gca, 'YDir', 'normal');
xlabel('Gap g (mm)');
ylabel('Common coordinate u (mm)');
title('(d) Shape field H(u,g)');
cb = colorbar; cb.Label.String = 'Normalized shape';

title(tl1, sprintf('Amplitude-shape calibration with g_0 = %.1f mm', g0), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_surface_calibration.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_surface_calibration.pdf'), 'ContentType', 'vector');

%% Figure 2
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
gapCats = categorical(string(gapList));

nexttile;
bar(gapCats, [mapSummaryTable.common_domain_rmse_to_ref, normSummaryTable.loo_normalized_rmse_to_ref], 'grouped');
xlabel('Gap g (mm)');
ylabel('RMSE to reference template (pF)');
title('(a) Before and after normalization');
legend({'Registered waveform', 'LOO normalized waveform'}, 'Location', 'northwest', 'FontSize', 7.0);

nexttile;
yyaxis left;
bar(gapCats, normSummaryTable.loo_fit_rmse_at_g, 0.5, 'FaceColor', [0.55 0.69 0.89], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('LOO fit RMSE at g (pF)');
yyaxis right;
plot(gapCats, normSummaryTable.loo_amp_fit_abs_err, '-o', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10]);
ylabel('Amplitude fit absolute error (pF)');
title('(b) Leave-one-gap-out fit quality');

exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_summary.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_summary.pdf'), 'ContentType', 'vector');

%% Figure 3
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
    title(sprintf('(%c) Gap %.1f mm: residual and \Delta u', char('b' + 2*(ir-1)), repRows(ir).gap_mm));
    legend('Location', 'best', 'FontSize', 6.8);
end

exportgraphics(fig3, fullfile(outDir, 'fig3_reference_gap_examples.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_reference_gap_examples.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Reference-gap normalization validation v2 (no vibration) ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Reference gap g0 = %.1f mm (largest gap)\n\n', g0);
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

function Ahat = predict_amplitude_leave_one_out(gKeep, AKeep, gTarget)
    gKeep = gKeep(:);
    AKeep = AKeep(:);
    pHat = fit_powerlaw_amplitude(gKeep, AKeep);
    Ahat = amplitude_model(pHat, gTarget);
end

function pHat = fit_powerlaw_amplitude(g, A)
    g = g(:);
    A = A(:);
    p0 = [max(A) * mean(g), 0.05, 1.2, min(A) * 0.5];
    obj = @(p) mean((amplitude_model(p, g) - A).^2);
    pHat = fminsearch(obj, p0, optimset('Display', 'off'));
end

function A = amplitude_model(p, g)
    a = abs(p(1));
    b = max(p(2), -min(g) + 1e-6);
    n = abs(p(3)) + 1e-6;
    c = max(p(4), 0);
    A = a ./ ((g + b) .^ n) + c;
end
