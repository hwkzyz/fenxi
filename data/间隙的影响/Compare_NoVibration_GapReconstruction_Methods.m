%% Compare candidate no-vibration gap-reconstruction operators
% This script treats V(x;g) = F(x,g) as the governing static response
% manifold and compares several practical implementations of the operator
%
%   T_{g->g0}: V(x;g) -> V_eq(x;g0)
%
% under a unified no-vibration validation protocol.
%
% Compared operators:
%   1) Amplitude-only baseline:
%        peak-centered identity coordinate + reference amplitude rescaling
%   2) Feature-landmark registration:
%        sparse equal-level landmarks on left/right branches
%   3) Area-CDF registration:
%        shape coordinate from normalized cumulative positive area
%   4) Dense equal-level registration:
%        current main method based on dense equal-level branch anchors
%   5) Regularized warp optimization:
%        monotone continuous warping with smoothness regularization
%
% All methods are evaluated using:
%   a) registered waveform consistency to the largest-gap reference;
%   b) normalized-shape consistency;
%   c) reconstructed reference-gap waveform consistency.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = resolve_data_file(scriptDir);

outDir = fullfile(scriptDir, 'compare_no_vibration_gap_reconstruction_results');
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
nU = 1401;
levelDense = (0.02:0.02:0.98)';
levelSparse = [0.10 0.30 0.50 0.70 0.90]';
sGridArea = linspace(0, 1, nU)';
opt.nLeftSeg = 5;
opt.nRightSeg = 5;
opt.lambda = 2e-3;
opt.maxIter = 2500;
opt.nDenseBranch = 240;

%% Load curves
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);
gRef = gapList(refIdx);

curve = repmat(struct( ...
    'gap', 0, 'x', [], 'yRaw', [], 'yBc', [], 'yNorm', [], ...
    'base', 0, 'amp', 0, 'xPeak', 0, 'xCenter', [], ...
    'xLeftDense', [], 'xRightDense', [], ...
    'xLeftSparse', [], 'xRightSparse', [], ...
    'cdfX', [], 'cdfS', []), nGap, 1);

for ig = 1:nGap
    x = xCell{ig}(:);
    y = yCell{ig}(:);
    [base, amp, yNorm] = normalize_waveform(y);
    yBc = y - base;
    [~, iPeak] = max(yBc);
    xPeak = x(iPeak);
    xCenter = x - xPeak;

    [xLeftDense, xRightDense] = level_crossings_normalized(xCenter, yNorm, levelDense);
    [xLeftSparse, xRightSparse] = level_crossings_normalized(xCenter, yNorm, levelSparse);
    [cdfX, cdfS] = build_area_cdf(xCenter, yBc);

    curve(ig).gap = gapList(ig);
    curve(ig).x = x;
    curve(ig).yRaw = y;
    curve(ig).yBc = yBc;
    curve(ig).yNorm = yNorm;
    curve(ig).base = base;
    curve(ig).amp = amp;
    curve(ig).xPeak = xPeak;
    curve(ig).xCenter = xCenter;
    curve(ig).xLeftDense = xLeftDense;
    curve(ig).xRightDense = xRightDense;
    curve(ig).xLeftSparse = xLeftSparse;
    curve(ig).xRightSparse = xRightSparse;
    curve(ig).cdfX = cdfX;
    curve(ig).cdfS = cdfS;
end

refCurve = curve(refIdx);
validRefDense = isfinite(refCurve.xLeftDense) & isfinite(refCurve.xRightDense);
uMin = refCurve.xLeftDense(find(validRefDense, 1, 'first'));
uMax = refCurve.xRightDense(find(validRefDense, 1, 'first'));
uGrid = linspace(uMin, uMax, nU)';

%% Candidate methods
methodDefs = struct( ...
    'key', { ...
        'amplitude_only', ...
        'feature_landmark', ...
        'area_cdf', ...
        'dense_equal_level', ...
        'regularized_warp'}, ...
    'label', { ...
        'Amplitude only', ...
        'Feature-landmark', ...
        'Area-CDF', ...
        'Dense equal-level', ...
        'Regularized warp'});
nMethod = numel(methodDefs);

result = repmat(struct( ...
    'key', "", 'label', "", 'Cmat', [], 'Ameas', [], 'Hmat', [], ...
    'Cref', [], 'Href', [], 'CrefRecon', [], 'summaryTable', table(), ...
    'globalTable', table()), nMethod, 1);

%% Run all methods
for im = 1:nMethod
    method = methodDefs(im);
    Cmat = nan(nGap, nU);

    for ig = 1:nGap
        Cmat(ig, :) = register_curve(curve(ig), refCurve, uGrid, ...
            method.key, levelDense, levelSparse, sGridArea, opt).';
    end

    Ameas = max(Cmat, [], 2);
    Hmat = Cmat ./ max(Ameas, eps);
    Cref = Cmat(refIdx, :).';
    Href = Hmat(refIdx, :).';
    Aref = Ameas(refIdx);

    CrefRecon = nan(size(Cmat));
    for ig = 1:nGap
        CrefRecon(ig, :) = (Aref * Hmat(ig, :).').';
    end

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
    globalTable = table( ...
        string(method.key), string(method.label), gRef, ...
        mean(summaryTable.registered_rmse_to_ref), ...
        max(summaryTable.registered_rmse_to_ref), ...
        mean(summaryTable.shape_rmse_to_ref), ...
        max(summaryTable.shape_rmse_to_ref), ...
        mean(summaryTable.reconstructed_rmse_to_ref), ...
        max(summaryTable.reconstructed_rmse_to_ref), ...
        'VariableNames', {'method_key', 'method_label', 'reference_gap_mm', ...
        'mean_registered_rmse_to_ref', 'max_registered_rmse_to_ref', ...
        'mean_shape_rmse_to_ref', 'max_shape_rmse_to_ref', ...
        'mean_reconstructed_rmse_to_ref', 'max_reconstructed_rmse_to_ref'});

    result(im).key = string(method.key);
    result(im).label = string(method.label);
    result(im).Cmat = Cmat;
    result(im).Ameas = Ameas;
    result(im).Hmat = Hmat;
    result(im).Cref = Cref;
    result(im).Href = Href;
    result(im).CrefRecon = CrefRecon;
    result(im).summaryTable = summaryTable;
    result(im).globalTable = globalTable;

    methodDir = fullfile(outDir, char(method.key));
    if ~isfolder(methodDir)
        mkdir(methodDir);
    end
    writetable(summaryTable, fullfile(methodDir, 'summary_by_gap.csv'));
    writetable(globalTable, fullfile(methodDir, 'global_summary.csv'));

    surfaceTable = table();
    surfaceTable.u_mm = uGrid;
    surfaceTable.C_ref = Cref;
    surfaceTable.H_ref = Href;
    for ig = 1:nGap
        surfaceTable.(sprintf('C_g_%03d', round(10 * gapList(ig)))) = Cmat(ig, :).';
        surfaceTable.(sprintf('H_g_%03d', round(10 * gapList(ig)))) = Hmat(ig, :).';
        surfaceTable.(sprintf('Cref_from_g_%03d', round(10 * gapList(ig)))) = CrefRecon(ig, :).';
    end
    writetable(surfaceTable, fullfile(methodDir, 'registered_and_reconstructed_surfaces.csv'));
end

%% Aggregate global comparison
globalAll = vertcat(result.globalTable);
writetable(globalAll, fullfile(outDir, 'method_global_comparison.csv'));

methodCats = categorical(globalAll.method_label, globalAll.method_label, 'Ordinal', true);

%% Figure 1: global metric comparison
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 11]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(methodCats, globalAll.mean_registered_rmse_to_ref, 0.65, ...
    'FaceColor', [0.58 0.75 0.88], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean registered RMSE (pF)');
title('(a) Common-domain mismatch');

nexttile;
bar(methodCats, globalAll.mean_shape_rmse_to_ref, 0.65, ...
    'FaceColor', [0.74 0.86 0.66], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean shape RMSE');
title('(b) Normalized-shape mismatch');

nexttile;
bar(methodCats, globalAll.mean_reconstructed_rmse_to_ref, 0.65, ...
    'FaceColor', [0.95 0.78 0.50], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean reconstructed RMSE (pF)');
title('(c) Reference-gap reconstruction error');

nexttile;
bar(methodCats, [globalAll.max_registered_rmse_to_ref, globalAll.max_reconstructed_rmse_to_ref], 'grouped');
ylabel('Max RMSE (pF)');
legend({'Registered', 'Reconstructed'}, 'Location', 'northwest', 'FontSize', 7.0);
title('(d) Worst-gap error');

title(tl1, sprintf('No-vibration operator comparison to reference gap g_0 = %.1f mm', gRef), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_method_global_comparison.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_method_global_comparison.pdf'), 'ContentType', 'vector');

%% Figure 2: gap-wise reconstruction RMSE
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 10]);
tl2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(nMethod);

nexttile; hold on;
for im = 1:nMethod
    plot(gapList, result(im).summaryTable.registered_rmse_to_ref, '-o', ...
        'Color', colors(im, :), 'MarkerFaceColor', colors(im, :), ...
        'DisplayName', char(result(im).label));
end
xlabel('Gap g (mm)');
ylabel('Registered RMSE to reference (pF)');
title('(a) Common-domain mismatch by gap');
legend('Location', 'eastoutside', 'FontSize', 7.0);

nexttile; hold on;
for im = 1:nMethod
    plot(gapList, result(im).summaryTable.reconstructed_rmse_to_ref, '-o', ...
        'Color', colors(im, :), 'MarkerFaceColor', colors(im, :));
end
xlabel('Gap g (mm)');
ylabel('Reconstructed RMSE to reference (pF)');
title('(b) Reference-gap reconstruction error by gap');

title(tl2, 'Gap-wise comparison of candidate operators', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_gapwise_comparison.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_gapwise_comparison.pdf'), 'ContentType', 'vector');

%% Figure 3: representative waveform comparison at the hardest non-reference gap
[~, hardGapIdx] = min(gapList);
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 12]);
tl3 = tiledlayout(fig3, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for im = 1:min(nMethod, 4)
    nexttile; hold on;
    regWave = result(im).Cmat(hardGapIdx, :).';
    recWave = result(im).CrefRecon(hardGapIdx, :).';
    plot(uGrid, regWave, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Registered');
    plot(uGrid, recWave, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Reconstructed');
    plot(uGrid, result(im).Cref, 'k--', 'DisplayName', 'True reference');
    xlabel('Reference coordinate u (mm)');
    ylabel('Capacitance (pF)');
    title(sprintf('(%c) %s', char('a' + im - 1), result(im).label));
    if im == 1
        legend('Location', 'best', 'FontSize', 6.8);
    end
end

title(tl3, sprintf('Representative gap %.1f mm -> reference gap %.1f mm', ...
    gapList(hardGapIdx), gRef), 'FontWeight', 'normal');
exportgraphics(fig3, fullfile(outDir, 'fig3_representative_waveforms.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_representative_waveforms.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== No-vibration gap-reconstruction operator comparison ===\n');
fprintf('Reference gap g0 = %.1f mm (largest gap)\n', gRef);
disp(globalAll(:, {'method_label', ...
    'mean_registered_rmse_to_ref', 'mean_shape_rmse_to_ref', ...
    'mean_reconstructed_rmse_to_ref'}));
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
    error('No supported stacked-curve data file was found in: %s', scriptDir);
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

function [cdfX, cdfS] = build_area_cdf(xCenter, yBc)
    xCenter = xCenter(:);
    yBc = yBc(:);
    yPos = max(yBc, 0);
    dx = diff(xCenter);
    w = zeros(size(xCenter));
    if numel(xCenter) >= 2
        w(1) = 0.5 * yPos(1) * dx(1);
        w(end) = 0.5 * yPos(end) * dx(end);
    end
    for i = 2:numel(xCenter)-1
        w(i) = 0.5 * yPos(i) * (xCenter(i+1) - xCenter(i-1));
    end
    w = max(w, 0);
    cum = cumsum(w);
    if cum(end) <= eps
        cdfS = linspace(0, 1, numel(xCenter))';
    else
        cdfS = cum / cum(end);
        cdfS(1) = 0;
        cdfS(end) = 1;
    end
    [cdfS, ia] = unique(cdfS, 'stable');
    cdfX = xCenter(ia);
end

function yOnU = register_curve(cur, ref, uGrid, methodKey, levelDense, levelSparse, sGridArea, opt)
    switch char(methodKey)
        case 'amplitude_only'
            yOnU = interp1(cur.xCenter, cur.yBc, uGrid, 'pchip', 'extrap');

        case 'feature_landmark'
            yOnU = register_by_landmarks(cur, ref, uGrid, levelSparse);

        case 'area_cdf'
            yOnU = register_by_area_cdf(cur, ref, uGrid, sGridArea);

        case 'dense_equal_level'
            yOnU = register_by_landmarks(cur, ref, uGrid, levelDense);

        case 'regularized_warp'
            yOnU = register_by_regularized_warp(cur, ref, uGrid, opt);

        otherwise
            error('Unknown method: %s', methodKey);
    end
end

function yOnU = register_by_landmarks(cur, ref, uGrid, levels)
    [xLeftCur, xRightCur] = level_crossings_normalized(cur.xCenter, cur.yNorm, levels);
    [xLeftRef, xRightRef] = level_crossings_normalized(ref.xCenter, ref.yNorm, levels);
    valid = isfinite(xLeftCur) & isfinite(xRightCur) & isfinite(xLeftRef) & isfinite(xRightRef);
    if nnz(valid) < max(2, ceil(numel(levels) / 2))
        yOnU = interp1(cur.xCenter, cur.yBc, uGrid, 'pchip', 'extrap');
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

    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), 'pchip', 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), 'pchip', 'extrap');
    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];
    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, 'pchip', 'extrap');
end

function yOnU = register_by_area_cdf(cur, ref, uGrid, sGrid)
    xCur = interp1(cur.cdfS, cur.cdfX, sGrid, 'linear', 'extrap');
    xRef = interp1(ref.cdfS, ref.cdfX, sGrid, 'linear', 'extrap');
    yS = interp1(cur.xCenter, cur.yBc, xCur, 'pchip', 'extrap');
    [xRefSorted, idx] = sort(xRef, 'ascend');
    ySorted = yS(idx);
    [xRefSorted, uniqIdx] = unique(xRefSorted, 'stable');
    ySorted = ySorted(uniqIdx);
    yOnU = interp1(xRefSorted, ySorted, uGrid, 'pchip', 'extrap');
    yOnU = yOnU(:);
end

function yOnU = register_by_regularized_warp(cur, ref, uGrid, opt)
    xMinCur = min(cur.xCenter);
    xMaxCur = max(cur.xCenter);
    uMin = min(uGrid);
    uMax = max(uGrid);

    uLeftDense = linspace(uMin, 0, opt.nDenseBranch)';
    uRightDense = linspace(0, uMax, opt.nDenseBranch)';
    refLeftDense = interp1(ref.xCenter, ref.yNorm, uLeftDense, 'pchip', 'extrap');
    refRightDense = interp1(ref.xCenter, ref.yNorm, uRightDense, 'pchip', 'extrap');

    thetaLeft0 = zeros(opt.nLeftSeg, 1);
    thetaRight0 = zeros(opt.nRightSeg, 1);

    leftObj = @(th) branch_warp_objective(th, xMinCur, 0, uMin, 0, ...
        uLeftDense, refLeftDense, cur.xCenter, cur.yNorm, opt.lambda);
    rightObj = @(th) branch_warp_objective(th, 0, xMaxCur, 0, uMax, ...
        uRightDense, refRightDense, cur.xCenter, cur.yNorm, opt.lambda);

    optsFs = optimset('Display', 'off', 'MaxIter', opt.maxIter, ...
        'MaxFunEvals', 10 * opt.maxIter, 'TolX', 1e-8, 'TolFun', 1e-10);
    thetaLeft = fminsearch(leftObj, thetaLeft0, optsFs);
    thetaRight = fminsearch(rightObj, thetaRight0, optsFs);

    [uLeftCtrl, xLeftCtrl] = build_monotone_ctrl(thetaLeft, xMinCur, 0, uMin, 0);
    [uRightCtrl, xRightCtrl] = build_monotone_ctrl(thetaRight, 0, xMaxCur, 0, uMax);

    xWarpLeft = interp1(uLeftCtrl, xLeftCtrl, uGrid(uGrid <= 0), 'pchip', 'extrap');
    xWarpRight = interp1(uRightCtrl, xRightCtrl, uGrid(uGrid > 0), 'pchip', 'extrap');
    xWarp = [xWarpLeft; xWarpRight];
    yOnU = interp1(cur.xCenter, cur.yBc, xWarp, 'pchip', 'extrap');
end

function obj = branch_warp_objective(theta, xStart, xEnd, uStart, uEnd, ...
        uDense, refDense, xCur, yCur, lambda)
    [uCtrl, xCtrl] = build_monotone_ctrl(theta, xStart, xEnd, uStart, uEnd);
    xWarp = interp1(uCtrl, xCtrl, uDense, 'pchip', 'extrap');
    yWarp = interp1(xCur, yCur, xWarp, 'pchip', 'extrap');
    dataTerm = mean((yWarp - refDense).^2, 'omitnan');

    nSeg = numel(theta);
    w = exp(theta(:));
    w = w / sum(w);
    wId = ones(nSeg, 1) / nSeg;
    regTerm = mean((w - wId).^2);
    obj = dataTerm + lambda * regTerm;
end

function [uCtrl, xCtrl] = build_monotone_ctrl(theta, xStart, xEnd, uStart, uEnd)
    nSeg = numel(theta);
    uCtrl = linspace(uStart, uEnd, nSeg + 1)';
    w = exp(theta(:));
    w = w / sum(w);
    xCtrl = xStart + [0; cumsum(w)] * (xEnd - xStart);
    xCtrl(end) = xEnd;
end
