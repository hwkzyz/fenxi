%% Improved unknown-gap validation using physical amplitude and level-set geometry
% This script tests a stricter decoupling model:
%   1) amplitude follows a clearance-dependent physical law;
%   2) normalized shape is represented by equal-level branch widths;
%   3) clearance and lateral shift are estimated jointly;
%   4) local observability is diagnosed using sensitivity-vector angles.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'improved_unknown_gap_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Paper-style figure defaults
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

colors = [0.85 0.35 0.10;
          0.85 0.55 0.10;
          0.84 0.73 0.18;
          0.20 0.60 0.60;
          0.16 0.48 0.66;
          0.10 0.25 0.50;
          0.32 0.32 0.32];

%% Load no-vibration clearance waveforms
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);
xGrid = xCell{1}(:);
S = zeros(nGap, numel(xGrid));
allAmp = zeros(nGap, 1);
for k = 1:nGap
    S(k, :) = yCell{k}(:)';
    [~, allAmp(k)] = normalize_with_edge_baseline(yCell{k}(:));
end

%% Sparse calibration setup
% Only these calibration gaps are treated as known.
calIdx = [1, 3, 5, 7];            % 0.2, 0.6, 1.0, 1.4 mm
testIdx = setdiff(1:nGap, calIdx); % 0.4, 0.8, 1.2 mm
levelGrid = linspace(0.05, 0.97, 50)';

oldModel = build_empirical_levelset_model(gapList(calIdx), xCell(calIdx), yCell(calIdx), levelGrid);
newModel = build_physical_levelset_model(gapList(calIdx), xCell(calIdx), yCell(calIdx), levelGrid);

%% 1) Held-out no-vibration waveform prediction
oldPred = zeros(numel(testIdx), numel(xGrid));
newPred = zeros(numel(testIdx), numel(xGrid));
oldRmse = zeros(numel(testIdx), 1);
newRmse = zeros(numel(testIdx), 1);
oldShapeRmse = zeros(numel(testIdx), 1);
newShapeRmse = zeros(numel(testIdx), 1);

for ii = 1:numel(testIdx)
    k = testIdx(ii);
    oldPred(ii, :) = synthesize_empirical(oldModel, gapList(k), xGrid);
    newPred(ii, :) = synthesize_physical(newModel, gapList(k), xGrid);

    oldRmse(ii) = rms_vec(S(k, :) - oldPred(ii, :));
    newRmse(ii) = rms_vec(S(k, :) - newPred(ii, :));
    oldShapeRmse(ii) = rms_vec(normalize_wave(S(k, :)) - normalize_wave(oldPred(ii, :)));
    newShapeRmse(ii) = rms_vec(normalize_wave(S(k, :)) - normalize_wave(newPred(ii, :)));
end

reconTable = table(gapList(testIdx), oldRmse, newRmse, oldShapeRmse, newShapeRmse, ...
    100 * (oldRmse - newRmse) ./ max(oldRmse, eps), ...
    'VariableNames', {'held_out_gap_mm', 'old_abs_rmse', 'new_abs_rmse', ...
    'old_shape_rmse', 'new_shape_rmse', 'abs_rmse_reduction_percent'});

%% 2) Online estimation when both gap and lateral shift are unknown
trueShiftSet = [-0.20, 0, 0.18]; % mm
modeNames = {'calibrated', 'shape_only'};
rows(numel(testIdx) * numel(trueShiftSet) * numel(modeNames), 1) = struct( ...
    'mode', "", ...
    'true_gap_mm', NaN, ...
    'true_shift_mm', NaN, ...
    'estimated_gap_mm', NaN, ...
    'estimated_shift_mm', NaN, ...
    'gap_error_mm', NaN, ...
    'shift_error_mm', NaN, ...
    'fit_rmse', NaN);
rowId = 0;

for ii = 1:numel(testIdx)
    k = testIdx(ii);
    trueGap = gapList(k);
    for jj = 1:numel(trueShiftSet)
        trueShift = trueShiftSet(jj);
        yOnline = interp1(xGrid - trueShift, S(k, :), xGrid, 'pchip', 'extrap');

        for mm = 1:numel(modeNames)
            mode = modeNames{mm};
            [estGap, estShift, bestRmse] = estimate_gap_shift(newModel, xGrid, yOnline, ...
                [min(gapList(calIdx)), max(gapList(calIdx))], [-0.45, 0.45], mode);

            rowId = rowId + 1;
            rows(rowId).mode = string(mode);
            rows(rowId).true_gap_mm = trueGap;
            rows(rowId).true_shift_mm = trueShift;
            rows(rowId).estimated_gap_mm = estGap;
            rows(rowId).estimated_shift_mm = estShift;
            rows(rowId).gap_error_mm = estGap - trueGap;
            rows(rowId).shift_error_mm = estShift - trueShift;
            rows(rowId).fit_rmse = bestRmse;
        end
    end
end
onlineTable = struct2table(rows);

[modeGroup, modeList] = findgroups(onlineTable.mode);
summary = table(modeList, ...
    splitapply(@(z) mean(abs(z)), onlineTable.gap_error_mm, modeGroup), ...
    splitapply(@(z) max(abs(z)), onlineTable.gap_error_mm, modeGroup), ...
    splitapply(@(z) mean(abs(z)), onlineTable.shift_error_mm, modeGroup), ...
    splitapply(@(z) max(abs(z)), onlineTable.shift_error_mm, modeGroup), ...
    splitapply(@mean, onlineTable.fit_rmse, modeGroup), ...
    'VariableNames', {'mode', 'mean_abs_gap_error_mm', 'max_abs_gap_error_mm', ...
    'mean_abs_shift_error_mm', 'max_abs_shift_error_mm', 'mean_fit_rmse'});

%% 3) Observability diagnostic
obsGap = linspace(min(gapList(calIdx)), max(gapList(calIdx)), 61)';
rhoRaw = zeros(numel(obsGap), 1);
rhoProjected = zeros(numel(obsGap), 1);
sinProjected = zeros(numel(obsGap), 1);

for i = 1:numel(obsGap)
    [rhoRaw(i), rhoProjected(i), sinProjected(i)] = observability_index(newModel, xGrid, obsGap(i));
end
obsTable = table(obsGap, rhoRaw, rhoProjected, sinProjected, ...
    'VariableNames', {'gap_mm', 'raw_abs_corr', 'projected_abs_corr', 'projected_sin_angle'});

%% Export results
writetable(reconTable, fullfile(outDir, 'held_out_reconstruction_comparison.csv'));
writetable(onlineTable, fullfile(outDir, 'online_gap_shift_estimation.csv'));
writetable(summary, fullfile(outDir, 'online_estimation_summary.csv'));
writetable(obsTable, fullfile(outDir, 'observability_index.csv'));

%% Figure 1: improved model validation
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(gapList, allAmp, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'all gaps');
plot(gapList(calIdx), allAmp(calIdx), 's', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], 'DisplayName', 'calibration');
gFine = linspace(min(gapList), max(gapList), 200)';
plot(gFine, amplitude_power_model(newModel.ampParam, gFine), '-', ...
    'Color', [0.10 0.25 0.50], 'DisplayName', 'physical fit');
plot(gFine, interp1(gapList(calIdx), allAmp(calIdx), gFine, 'pchip'), '--', ...
    'Color', [0.35 0.35 0.35], 'DisplayName', 'pchip');
xlabel('Gap (mm)');
ylabel('Amplitude (pF)');
title('(a) Amplitude model');
legend('Location', 'northeast', 'FontSize', 7.5);

nexttile; hold on;
for ii = 1:numel(testIdx)
    k = testIdx(ii);
    plot(xGrid, S(k, :), '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('data g=%.1f', gapList(k)));
    plot(xGrid, newPred(ii, :), '--', 'Color', colors(k, :), ...
        'DisplayName', sprintf('pred g=%.1f', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
title('(b) Held-out prediction');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
bar(gapList(testIdx) - 0.025, oldRmse, 0.05, ...
    'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'old');
bar(gapList(testIdx) + 0.025, newRmse, 0.05, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'new');
xlabel('Held-out gap (mm)');
ylabel('Absolute RMSE');
title('(c) Absolute prediction error');
legend('Location', 'northeast', 'FontSize', 7.5);

nexttile; hold on;
plot(obsTable.gap_mm, obsTable.projected_abs_corr, '-', 'Color', [0.85 0.35 0.10]);
plot(obsTable.gap_mm, obsTable.projected_sin_angle, '-', 'Color', [0.10 0.25 0.50]);
yline(0.2, 'k--', 'LineWidth', 0.8);
xlabel('Gap (mm)');
ylabel('Index');
title('(d) Local observability');
legend({'|corr|', 'sin(angle)', 'weak limit'}, 'Location', 'east', 'FontSize', 7.5);

title(tl, 'Improved unknown-gap waveform manifold', 'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_improved_validation.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_improved_validation.pdf'), 'ContentType', 'vector');

%% Figure 2: online identification
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for mm = 1:numel(modeNames)
    mode = string(modeNames{mm});
    T = onlineTable(onlineTable.mode == mode, :);

    nexttile; hold on;
    plot(T.true_gap_mm, T.estimated_gap_mm, 'o', ...
        'Color', colors(2 * mm, :), 'MarkerFaceColor', colors(2 * mm, :));
    plot([min(gapList), max(gapList)], [min(gapList), max(gapList)], 'k--');
    xlabel('True gap (mm)');
    ylabel('Estimated gap (mm)');
    title(sprintf('(%c) Gap, %s', 'a' + 2 * (mm - 1), replace(mode, '_', ' ')));

    nexttile; hold on;
    plot(T.true_shift_mm, T.estimated_shift_mm, 's', ...
        'Color', colors(2 * mm + 1, :), 'MarkerFaceColor', colors(2 * mm + 1, :));
    plot([-0.25, 0.25], [-0.25, 0.25], 'k--');
    xlabel('True shift (mm)');
    ylabel('Estimated shift (mm)');
    title(sprintf('(%c) Shift, %s', 'b' + 2 * (mm - 1), replace(mode, '_', ' ')));
end
title(tl2, 'Online joint estimation under synthetic lateral shifts', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_online_estimation.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_online_estimation.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Improved unknown-gap validation ===\n');
fprintf('Known calibration gaps: %s mm\n', num2str(gapList(calIdx)', '%.1f '));
fprintf('Held-out gaps: %s mm\n\n', num2str(gapList(testIdx)', '%.1f '));
disp(reconTable);
fprintf('\nAmplitude model: A(g) = alpha / (g + g0)^n\n');
fprintf('alpha = %.6g, g0 = %.6g mm, n = %.6g\n', ...
    newModel.ampParam.alpha, newModel.ampParam.g0, newModel.ampParam.n);
fprintf('\nOnline estimation summary:\n');
disp(summary(:, {'mode', 'mean_abs_gap_error_mm', 'max_abs_gap_error_mm', ...
    'mean_abs_shift_error_mm', 'max_abs_shift_error_mm'}));
fprintf('\nProjected observability: min sin(angle) = %.4f, max |corr| = %.4f\n', ...
    min(obsTable.projected_sin_angle), max(obsTable.projected_abs_corr));
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function model = build_empirical_levelset_model(gaps, xCell, yCell, levels)
    n = numel(gaps);
    base = zeros(n, 1);
    amp = zeros(n, 1);
    xPeak = zeros(n, 1);
    xMin = zeros(n, 1);
    xMax = zeros(n, 1);
    xLeft = zeros(numel(levels), n);
    xRight = zeros(numel(levels), n);

    for i = 1:n
        x = xCell{i}(:);
        y = yCell{i}(:);
        [base(i), amp(i), yn] = normalize_with_edge_baseline(y);
        [~, iPeak] = max(yn);
        xPeak(i) = x(iPeak);
        xMin(i) = min(x);
        xMax(i) = max(x);
        [xLeft(:, i), xRight(:, i)] = branch_level_positions(x, yn, iPeak, levels);
    end

    model.gaps = gaps(:);
    model.levels = levels(:);
    model.base = base;
    model.amp = amp;
    model.xPeak = xPeak;
    model.xMin = xMin;
    model.xMax = xMax;
    model.xLeft = xLeft;
    model.xRight = xRight;
end

function model = build_physical_levelset_model(allGapsForOutput, xCell, yCell, levels)
    gaps = allGapsForOutput(:);
    n = numel(gaps);
    base = zeros(n, 1);
    amp = zeros(n, 1);
    xPeak = zeros(n, 1);
    xMin = zeros(n, 1);
    xMax = zeros(n, 1);
    xLeft = zeros(numel(levels), n);
    xRight = zeros(numel(levels), n);

    for i = 1:n
        x = xCell{i}(:);
        y = yCell{i}(:);
        [base(i), amp(i), yn] = normalize_with_edge_baseline(y);
        [~, iPeak] = max(yn);
        xPeak(i) = x(iPeak);
        xMin(i) = min(x);
        xMax(i) = max(x);
        [xLeft(:, i), xRight(:, i)] = branch_level_positions(x, yn, iPeak, levels);
    end

    width = xRight - xLeft;
    center = 0.5 * (xRight + xLeft);

    model.gaps = gaps;
    model.levels = levels(:);
    model.base = base;
    model.amp = amp;
    model.ampParam = fit_power_amplitude(gaps, amp);
    model.center = center;
    model.width = width;
    model.xPeak = xPeak;
    model.xMin = min(xMin);
    model.xMax = max(xMax);

    % Store all measured amplitudes for plotting. If the caller only uses
    % sparse calibration gaps, this vector is later overwritten by the main script.
    model.measuredAmp = amp;
end

function y = synthesize_empirical(model, gap, xGrid)
    base = interp1(model.gaps, model.base, gap, 'pchip');
    amp = interp1(model.gaps, model.amp, gap, 'pchip');
    xPeak = interp1(model.gaps, model.xPeak, gap, 'pchip');
    xMin = interp1(model.gaps, model.xMin, gap, 'linear');
    xMax = interp1(model.gaps, model.xMax, gap, 'linear');
    xLeft = interp1(model.gaps, model.xLeft', gap, 'pchip')';
    xRight = interp1(model.gaps, model.xRight', gap, 'pchip')';
    y = reconstruct_from_branches(xGrid, base, amp, xMin, xMax, xPeak, model.levels, xLeft, xRight);
end

function y = synthesize_physical(model, gap, xGrid)
    base = interp1(model.gaps, model.base, gap, 'pchip');
    amp = amplitude_power_model(model.ampParam, gap);
    width = interp1(model.gaps, model.width', gap, 'pchip')';
    center = interp1(model.gaps, model.center', gap, 'pchip')';
    xLeft = center - 0.5 * width;
    xRight = center + 0.5 * width;
    xPeak = interp1(model.gaps, model.xPeak, gap, 'pchip');
    y = reconstruct_from_branches(xGrid, base, amp, model.xMin, model.xMax, ...
        xPeak, model.levels, xLeft, xRight);
end

function y = reconstruct_from_branches(xGrid, base, amp, xMin, xMax, xPeak, levels, xLeft, xRight)
    xLeftPts = [xMin; xLeft(:); xPeak];
    yLeftPts = [0; levels(:); 1];
    [xLeftPts, ia] = unique(xLeftPts, 'stable');
    yLeftPts = yLeftPts(ia);

    xRightPts = [xPeak; flipud(xRight(:)); xMax];
    yRightPts = [1; flipud(levels(:)); 0];
    [xRightPts, ia] = unique(xRightPts, 'stable');
    yRightPts = yRightPts(ia);

    yn = zeros(size(xGrid(:)));
    leftMask = xGrid(:) <= xPeak;
    yn(leftMask) = interp1(xLeftPts, yLeftPts, xGrid(leftMask), 'pchip', 'extrap');
    yn(~leftMask) = interp1(xRightPts, yRightPts, xGrid(~leftMask), 'pchip', 'extrap');
    yn = min(max(yn, 0), 1);
    y = (base + amp * yn)';
end

function [estGap, estShift, bestRmse] = estimate_gap_shift(model, xGrid, yOnline, gapRange, shiftRange, mode)
    gapGrid = linspace(gapRange(1), gapRange(2), 61);
    shiftGrid = linspace(shiftRange(1), shiftRange(2), 91);
    bestRmse = inf;
    estGap = NaN;
    estShift = NaN;

    for ig = 1:numel(gapGrid)
        yBase = synthesize_physical(model, gapGrid(ig), xGrid);
        for is = 1:numel(shiftGrid)
            yTpl = interp1(xGrid - shiftGrid(is), yBase, xGrid, 'pchip', 'extrap');
            rmse = objective_rmse(yOnline, yTpl, mode);
            if rmse < bestRmse
                bestRmse = rmse;
                estGap = gapGrid(ig);
                estShift = shiftGrid(is);
            end
        end
    end

    p0 = [estGap, estShift];
    obj = @(p) estimate_objective(model, xGrid, yOnline, p, gapRange, shiftRange, mode);
    pHat = fminsearch(obj, p0, optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10));
    estGap = min(max(pHat(1), gapRange(1)), gapRange(2));
    estShift = min(max(pHat(2), shiftRange(1)), shiftRange(2));
    bestRmse = obj([estGap, estShift]);
end

function val = estimate_objective(model, xGrid, yOnline, p, gapRange, shiftRange, mode)
    gap = p(1);
    shift = p(2);
    if gap < gapRange(1) || gap > gapRange(2) || shift < shiftRange(1) || shift > shiftRange(2)
        val = 1e3 + 1e2 * sum(abs(p));
        return;
    end
    yBase = synthesize_physical(model, gap, xGrid);
    yTpl = interp1(xGrid - shift, yBase, xGrid, 'pchip', 'extrap');
    val = objective_rmse(yOnline, yTpl, mode);
end

function rmse = objective_rmse(yObs, yTpl, mode)
    switch string(mode)
        case "calibrated"
            rmse = rms_vec(yObs(:) - yTpl(:));
        case "shape_only"
            A = [ones(numel(yTpl), 1), yTpl(:)];
            theta = A \ yObs(:);
            resid = yObs(:) - A * theta;
            rmse = rms_vec(resid);
        otherwise
            error('Unknown objective mode: %s', mode);
    end
end

function param = fit_power_amplitude(g, amp)
    g = g(:);
    amp = amp(:);
    obj = @(p) sum((log(max(amp, eps)) - ...
        (p(1) - exp(p(3)) .* log(g + exp(p(2))))).^2);
    p0 = [log(max(amp) * (min(g) + 0.1)), log(0.1), log(1)];
    p = fminsearch(obj, p0, optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12));
    param.alpha = exp(p(1));
    param.g0 = exp(p(2));
    param.n = exp(p(3));
end

function A = amplitude_power_model(param, g)
    A = param.alpha ./ (g + param.g0) .^ param.n;
end

function [rhoRaw, rhoProjected, sinProjected] = observability_index(model, xGrid, gap)
    dg = 1e-3;
    dx = 1e-3;
    y0 = synthesize_physical(model, gap, xGrid);
    yg1 = synthesize_physical(model, gap + dg, xGrid);
    yg0 = synthesize_physical(model, gap - dg, xGrid);
    sx = (interp1(xGrid - dx, y0, xGrid, 'pchip', 'extrap') - ...
          interp1(xGrid + dx, y0, xGrid, 'pchip', 'extrap')) / (2 * dx);
    sg = (yg1 - yg0) / (2 * dg);

    rhoRaw = abs(dot_norm(sg(:), sx(:)));

    nuisance = [ones(numel(y0), 1), y0(:)];
    sgp = project_out(sg(:), nuisance);
    sxp = project_out(sx(:), nuisance);
    rhoProjected = abs(dot_norm(sgp, sxp));
    sinProjected = sqrt(max(0, 1 - rhoProjected ^ 2));
end

function z = project_out(v, basis)
    z = v - basis * (basis \ v);
end

function c = dot_norm(a, b)
    c = (a(:)' * b(:)) / max(norm(a(:)) * norm(b(:)), eps);
end

function [base, amp, yn] = normalize_with_edge_baseline(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yn = y0 / max(amp, eps);
end

function [xLeft, xRight] = branch_level_positions(x, yn, iPeak, levels)
    x = x(:);
    yn = yn(:);
    levels = levels(:);

    xL = x(1:iPeak);
    yL = yn(1:iPeak);
    [yL, ia] = unique(yL, 'stable');
    xL = xL(ia);

    xR = x(iPeak:end);
    yR = yn(iPeak:end);
    [yRAsc, order] = sort(yR, 'ascend');
    xRAsc = xR(order);
    [yRAsc, ia] = unique(yRAsc, 'stable');
    xRAsc = xRAsc(ia);

    xLeft = interp1(yL, xL, levels, 'pchip', 'extrap');
    xRight = interp1(yRAsc, xRAsc, levels, 'pchip', 'extrap');
end

function y = normalize_wave(y)
    y = y(:)';
    y = y - min(y);
    y = y / max(max(y), eps);
end

function r = rms_vec(v)
    r = sqrt(mean(v(:) .^ 2));
end

function [gapList, xCell, yCell] = load_stacked_comsol_curves(filePath)
    txt = fileread(filePath);
    tok = regexp(txt, '%\s*([0-9.]+)mm', 'tokens');
    gapList = cellfun(@(c) str2double(c{1}), tok(:));
    if isempty(gapList)
        gapList = (0.2:0.2:1.4)';
    end

    data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
    x = data(:, 1);
    y = data(:, 2);

    breakIdx = [find(diff(x) < 0); numel(x)];
    startIdx = [1; breakIdx(1:end-1) + 1];
    nCurve = numel(breakIdx);
    if nCurve == 7
        gapList = (0.2:0.2:1.4)';
    elseif numel(gapList) ~= nCurve
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
