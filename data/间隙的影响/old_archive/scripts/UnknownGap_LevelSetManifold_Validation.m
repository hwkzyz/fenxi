%% Unknown-gap validation with a level-set waveform manifold
% Goal:
% Use only several no-vibration calibration waveforms to predict unmeasured
% clearance waveforms, then estimate both clearance and shift for online
% waveforms that may contain vibration-induced lateral displacement.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'unknown_gap_levelset_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Figure defaults
fontName = 'Times New Roman';
set(groot, 'DefaultAxesFontName', fontName, ...
    'DefaultTextFontName', fontName, ...
    'DefaultLegendFontName', fontName, ...
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

%% Load data
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);
xGrid = xCell{1}(:);
S = zeros(nGap, numel(xGrid));
for k = 1:nGap
    S(k, :) = yCell{k}(:)';
end

%% Sparse calibration setup
% Pretend only these gap waveforms are measured during calibration.
calIdx = [1, 3, 5, 7];       % 0.2, 0.6, 1.0, 1.4 mm
testIdx = setdiff(1:nGap, calIdx); % 0.4, 0.8, 1.2 mm

levelGrid = linspace(0.04, 0.98, 48)';
model = build_levelset_model(gapList(calIdx), xCell(calIdx), yCell(calIdx), levelGrid);

%% 1) Predict unknown no-vibration clearance waveforms
pred = zeros(numel(testIdx), numel(xGrid));
reconRmse = zeros(numel(testIdx), 1);
reconShapeRmse = zeros(numel(testIdx), 1);

for ii = 1:numel(testIdx)
    k = testIdx(ii);
    pred(ii, :) = synthesize_waveform(model, gapList(k), xGrid);
    reconRmse(ii) = sqrt(mean((S(k, :) - pred(ii, :)) .^ 2));
    reconShapeRmse(ii) = sqrt(mean((normalize_wave(S(k, :)) - normalize_wave(pred(ii, :))) .^ 2));
end

%% 2) Online estimation when gap is unknown and shift may exist
trueShiftSet = [-0.20, 0, 0.18]; % mm, synthetic vibration-like lateral shifts
numCases = numel(testIdx) * numel(trueShiftSet);
rows(numCases, 1) = struct( ...
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

        [estGap, estShift, bestRmse] = estimate_gap_shift(model, xGrid, yOnline, ...
            [min(gapList(calIdx)), max(gapList(calIdx))], [-0.45, 0.45]);

        rowId = rowId + 1;
        rows(rowId).true_gap_mm = trueGap;
        rows(rowId).true_shift_mm = trueShift;
        rows(rowId).estimated_gap_mm = estGap;
        rows(rowId).estimated_shift_mm = estShift;
        rows(rowId).gap_error_mm = estGap - trueGap;
        rows(rowId).shift_error_mm = estShift - trueShift;
        rows(rowId).fit_rmse = bestRmse;
    end
end
onlineTable = struct2table(rows);

%% Export tables
reconTable = table(gapList(testIdx), reconRmse, reconShapeRmse, ...
    'VariableNames', {'held_out_gap_mm', 'reconstruction_rmse', 'shape_rmse'});
writetable(reconTable, fullfile(outDir, 'held_out_gap_reconstruction.csv'));
writetable(onlineTable, fullfile(outDir, 'online_unknown_gap_shift_estimation.csv'));

summary = table( ...
    mean(abs(onlineTable.gap_error_mm)), max(abs(onlineTable.gap_error_mm)), ...
    mean(abs(onlineTable.shift_error_mm)), max(abs(onlineTable.shift_error_mm)), ...
    mean(onlineTable.fit_rmse), ...
    'VariableNames', {'mean_abs_gap_error_mm', 'max_abs_gap_error_mm', ...
    'mean_abs_shift_error_mm', 'max_abs_shift_error_mm', 'mean_fit_rmse'});
writetable(summary, fullfile(outDir, 'online_estimation_summary.csv'));

%% Figure 1: held-out waveform prediction
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 9]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for kk = 1:numel(calIdx)
    k = calIdx(kk);
    plot(xGrid, S(k, :), '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('cal g = %.1f', gapList(k)));
end
for ii = 1:numel(testIdx)
    k = testIdx(ii);
    plot(xGrid, S(k, :), '--', 'Color', colors(k, :), ...
        'DisplayName', sprintf('test g = %.1f', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'eastoutside', 'Box', 'off');
title('(a) Sparse calibration and held-out gaps');

nexttile; hold on;
for ii = 1:numel(testIdx)
    k = testIdx(ii);
    plot(xGrid, S(k, :), 'o', 'Color', colors(k, :), 'MarkerSize', 3, ...
        'HandleVisibility', 'off');
    plot(xGrid, pred(ii, :), '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('pred g = %.1f', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'best', 'Box', 'off');
title('(b) Held-out waveform prediction');

nexttile; hold on;
for ii = 1:numel(testIdx)
    k = testIdx(ii);
    plot(xGrid, pred(ii, :) - S(k, :), '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('g = %.1f', gapList(k)));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Position x (mm)');
ylabel('Prediction residual');
legend('Location', 'best', 'Box', 'off');
title('(c) Prediction residuals');

nexttile;
bar(gapList(testIdx), reconShapeRmse, 0.5, ...
    'FaceColor', [0.5529, 0.6941, 0.8863], 'EdgeColor', 'k', 'LineWidth', 0.5);
xlabel('Held-out clearance g (mm)');
ylabel('Shape RMSE');
title('(d) Shape prediction error');

exportgraphics(fig1, fullfile(outDir, 'fig1_held_out_gap_prediction.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_held_out_gap_prediction.pdf'), 'ContentType', 'vector');

%% Figure 2: online unknown gap and shift estimation
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 9]);
tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(onlineTable.true_gap_mm, onlineTable.estimated_gap_mm, 'o', ...
    'Color', [0.10, 0.25, 0.50], 'MarkerFaceColor', [0.10, 0.25, 0.50]);
plot([min(gapList), max(gapList)], [min(gapList), max(gapList)], 'k--');
xlabel('True gap (mm)');
ylabel('Estimated gap (mm)');
title('(a) Gap estimation');

nexttile; hold on;
plot(onlineTable.true_shift_mm, onlineTable.estimated_shift_mm, 's', ...
    'Color', [0.85, 0.35, 0.10], 'MarkerFaceColor', [0.85, 0.35, 0.10]);
plot([-0.25, 0.25], [-0.25, 0.25], 'k--');
xlabel('True shift (mm)');
ylabel('Estimated shift (mm)');
title('(b) Shift estimation');

nexttile;
bar(1:height(onlineTable), onlineTable.gap_error_mm, 0.6, ...
    'FaceColor', [0.5529, 0.6941, 0.8863], 'EdgeColor', 'k', 'LineWidth', 0.5);
xticks(1:height(onlineTable));
xticklabels(compose('g=%.1f, dx=%.2f', onlineTable.true_gap_mm, onlineTable.true_shift_mm));
xtickangle(35);
ylabel('Gap error (mm)');
title('(c) Gap errors');

nexttile;
bar(1:height(onlineTable), onlineTable.shift_error_mm, 0.6, ...
    'FaceColor', [0.8431, 0.8902, 0.7490], 'EdgeColor', 'k', 'LineWidth', 0.5);
xticks(1:height(onlineTable));
xticklabels(compose('g=%.1f, dx=%.2f', onlineTable.true_gap_mm, onlineTable.true_shift_mm));
xtickangle(35);
ylabel('Shift error (mm)');
title('(d) Shift errors');

exportgraphics(fig2, fullfile(outDir, 'fig2_online_gap_shift_estimation.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_online_gap_shift_estimation.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Unknown-gap level-set manifold validation ===\n');
fprintf('Known calibration gaps: %s mm\n', num2str(gapList(calIdx)', '%.1f '));
fprintf('Held-out gaps: %s mm\n\n', num2str(gapList(testIdx)', '%.1f '));
disp(reconTable);
fprintf('\nOnline estimation summary:\n');
disp(summary);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function model = build_levelset_model(gaps, xCell, yCell, levels)
    n = numel(gaps);
    xLeft = zeros(numel(levels), n);
    xRight = zeros(numel(levels), n);
    base = zeros(n, 1);
    amp = zeros(n, 1);
    xPeak = zeros(n, 1);
    xMin = zeros(n, 1);
    xMax = zeros(n, 1);

    for i = 1:n
        x = xCell{i}(:);
        y = yCell{i}(:);
        base(i) = min(y);
        y0 = y - base(i);
        amp(i) = max(y0);
        yn = y0 / max(amp(i), eps);
        [~, iPeak] = max(yn);
        xPeak(i) = x(iPeak);
        xMin(i) = min(x);
        xMax(i) = max(x);

        xLeft(:, i) = branch_level_positions(x(1:iPeak), yn(1:iPeak), levels);
        xRight(:, i) = branch_level_positions(flipud(x(iPeak:end)), flipud(yn(iPeak:end)), levels);
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

function xPos = branch_level_positions(xBranch, yBranch, levels)
    [yUniq, ia] = unique(yBranch, 'stable');
    xUniq = xBranch(ia);
    xPos = interp1(yUniq, xUniq, levels, 'pchip', 'extrap');
end

function y = synthesize_waveform(model, gap, xGrid)
    base = interp1(model.gaps, model.base, gap, 'pchip');
    amp = interp1(model.gaps, model.amp, gap, 'pchip');
    xPeak = interp1(model.gaps, model.xPeak, gap, 'pchip');
    xMin = interp1(model.gaps, model.xMin, gap, 'linear');
    xMax = interp1(model.gaps, model.xMax, gap, 'linear');
    xLeft = interp1(model.gaps, model.xLeft', gap, 'pchip')';
    xRight = interp1(model.gaps, model.xRight', gap, 'pchip')';

    levels = model.levels;
    xLeftPts = [xMin; xLeft; xPeak];
    yLeftPts = [0; levels; 1];
    xRightPts = [xPeak; flipud(xRight); xMax];
    yRightPts = [1; flipud(levels); 0];

    [xLeftPts, ia] = unique(xLeftPts, 'stable');
    yLeftPts = yLeftPts(ia);
    [xRightPts, ib] = unique(xRightPts, 'stable');
    yRightPts = yRightPts(ib);

    yn = zeros(size(xGrid));
    leftMask = xGrid <= xPeak;
    rightMask = ~leftMask;
    yn(leftMask) = interp1(xLeftPts, yLeftPts, xGrid(leftMask), 'pchip', 'extrap');
    yn(rightMask) = interp1(xRightPts, yRightPts, xGrid(rightMask), 'pchip', 'extrap');
    yn = min(max(yn, 0), 1);
    y = base + amp * yn;
end

function [estGap, estShift, bestRmse] = estimate_gap_shift(model, xGrid, yOnline, gapRange, shiftRange)
    gapGrid = linspace(gapRange(1), gapRange(2), 121);
    shiftGrid = linspace(shiftRange(1), shiftRange(2), 121);
    bestRmse = inf;
    estGap = gapGrid(1);
    estShift = 0;

    for ig = 1:numel(gapGrid)
        yBase = synthesize_waveform(model, gapGrid(ig), xGrid);
        for is = 1:numel(shiftGrid)
            yTpl = interp1(xGrid - shiftGrid(is), yBase, xGrid, 'pchip', 'extrap');
            rmse = scaled_rmse(yOnline, yTpl);
            if rmse < bestRmse
                bestRmse = rmse;
                estGap = gapGrid(ig);
                estShift = shiftGrid(is);
            end
        end
    end

    p0 = [estGap, estShift];
    obj = @(p) estimate_objective(model, xGrid, yOnline, p, gapRange, shiftRange);
    pHat = fminsearch(obj, p0, optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10));
    estGap = min(max(pHat(1), gapRange(1)), gapRange(2));
    estShift = min(max(pHat(2), shiftRange(1)), shiftRange(2));
    bestRmse = obj([estGap, estShift]);
end

function val = estimate_objective(model, xGrid, yOnline, p, gapRange, shiftRange)
    gap = p(1);
    shift = p(2);
    if gap < gapRange(1) || gap > gapRange(2) || shift < shiftRange(1) || shift > shiftRange(2)
        val = 1e3 + 1e3 * (max(0, gapRange(1) - gap)^2 + max(0, gap - gapRange(2))^2 + ...
            max(0, shiftRange(1) - shift)^2 + max(0, shift - shiftRange(2))^2);
        return;
    end
    yBase = synthesize_waveform(model, gap, xGrid);
    yTpl = interp1(xGrid - shift, yBase, xGrid, 'pchip', 'extrap');
    val = scaled_rmse(yOnline, yTpl);
end

function rmse = scaled_rmse(yObs, yTpl)
    A = [ones(numel(yTpl), 1), yTpl(:)];
    theta = A \ yObs(:);
    yFit = A * theta;
    rmse = sqrt(mean((yObs(:) - yFit) .^ 2));
end

function yn = normalize_wave(y)
    yn = (y - min(y)) ./ max(range(y), eps);
end

function [gapList, xCell, yCell] = load_stacked_comsol_curves(filePath)
    txt = fileread(filePath);
    gapLine = regexp(txt, '%\s*([0-9.]+mm\s*,\s*[0-9.]+mm.*)', 'tokens', 'once');
    if isempty(gapLine)
        error('Gap header was not found in %s', filePath);
    end

    gapText = gapLine{1};
    gapStr = regexp(gapText, '[0-9.]+(?=mm)', 'match');
    gapVals = str2double(gapStr(:));
    if contains(gapText, '...') && numel(gapVals) == 3
        gapList = (gapVals(1):gapVals(2) - gapVals(1):gapVals(3))';
    else
        gapList = gapVals;
    end

    data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data), 2), :);
    x = data(:, 1);
    y = data(:, 2);

    cutIdx = find(diff(x) < 0);
    startIdx = [1; cutIdx + 1];
    endIdx = [cutIdx; numel(x)];
    if numel(startIdx) ~= numel(gapList)
        error('Segment count (%d) does not match gap count (%d).', ...
            numel(startIdx), numel(gapList));
    end

    xCell = cell(numel(gapList), 1);
    yCell = cell(numel(gapList), 1);
    for k = 1:numel(gapList)
        idx = startIdx(k):endIdx(k);
        xCell{k} = x(idx);
        yCell{k} = y(idx);
    end
end
