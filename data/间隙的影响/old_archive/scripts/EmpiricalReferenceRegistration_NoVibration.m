%% Empirical reference registration for no-vibration waveforms
% Theory behind the validation:
% 1) Under no-vibration conditions, all waveforms should share one common
%    latent geometry-passing process.
% 2) Use one reference waveform as the empirical latent basis instead of
%    assuming a Gaussian or another parametric kernel.
% 3) Learn branch-wise monotone coordinate maps from sparse level sets, then
%    test whether the rest of the waveform collapses automatically.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'empirical_reference_registration_results');
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
xRefGrid = xCell{1}(:);
S = zeros(nGap, numel(xRefGrid));
for k = 1:nGap
    S(k, :) = yCell{k}(:)';
end

% Use the middle clearance as the empirical reference basis.
[~, iRef] = min(abs(gapList - median(gapList)));
gRef = gapList(iRef);

%% Preprocess into baseline-removed normalized waveforms
proc = repmat(struct('gap', [], 'x', [], 'y', [], 'yn', [], 'base', [], ...
    'amp', [], 'xPeak', [], 'iPeak', []), nGap, 1);

for k = 1:nGap
    x = xCell{k}(:);
    y = yCell{k}(:);
    y0 = y - min(y);
    amp = max(y0);
    yn = y0 / max(amp, eps);
    [~, iPeak] = max(yn);

    proc(k).gap = gapList(k);
    proc(k).x = x;
    proc(k).y = y;
    proc(k).yn = yn;
    proc(k).base = min(y);
    proc(k).amp = amp;
    proc(k).xPeak = x(iPeak);
    proc(k).iPeak = iPeak;
end

ref = proc(iRef);

%% Sparse control levels for learning the warping maps
levelTrain = [0.08 0.16 0.24 0.36 0.50 0.64 0.78 0.90 0.97]';

simpleNorm = zeros(nGap, numel(xRefGrid));
registeredNorm = zeros(nGap, numel(xRefGrid));
simpleRms = zeros(nGap, 1);
registeredRms = zeros(nGap, 1);
simpleShift = zeros(nGap, 1);
registeredShift = zeros(nGap, 1);
branchMapError = zeros(nGap, 1);

for k = 1:nGap
    simpleNorm(k, :) = interp1(proc(k).x - proc(k).xPeak, proc(k).yn, ...
        xRefGrid - ref.xPeak, 'pchip', 'extrap');
end
simpleRef = simpleNorm(iRef, :);

fitRange = [-0.5, 0.5];
fitMask = xRefGrid >= min(xRefGrid) + 0.4 & xRefGrid <= max(xRefGrid) - 0.4;

for k = 1:nGap
    if k == iRef
        registeredNorm(k, :) = ref.yn;
        continue;
    end

    cur = proc(k);
    [xRefL, xCurL] = level_position_pairs(ref.x, ref.yn, ref.iPeak, cur.x, cur.yn, cur.iPeak, levelTrain, 'left');
    [xRefR, xCurR] = level_position_pairs(ref.x, ref.yn, ref.iPeak, cur.x, cur.yn, cur.iPeak, levelTrain, 'right');

    % Add anchor points near outer tails and at the peak.
    xRefLeftEdge = ref.x(1);
    xCurLeftEdge = cur.x(1);
    xRefRightEdge = ref.x(end);
    xCurRightEdge = cur.x(end);
    xRefPeak = ref.x(ref.iPeak);
    xCurPeak = cur.x(cur.iPeak);

    xRefCtrlL = [xRefLeftEdge; xRefL; xRefPeak];
    xCurCtrlL = [xCurLeftEdge; xCurL; xCurPeak];
    xRefCtrlR = [xRefPeak; xRefR; xRefRightEdge];
    xCurCtrlR = [xCurPeak; xCurR; xCurRightEdge];

    xWarp = zeros(size(xRefGrid));
    leftMask = xRefGrid <= xRefPeak;
    rightMask = xRefGrid > xRefPeak;
    xWarp(leftMask) = interp1(xRefCtrlL, xCurCtrlL, xRefGrid(leftMask), 'pchip', 'extrap');
    xWarp(rightMask) = interp1(xRefCtrlR, xCurCtrlR, xRefGrid(rightMask), 'pchip', 'extrap');

    yWarp = interp1(cur.x, cur.yn, xWarp, 'pchip', 'extrap');
    registeredNorm(k, :) = yWarp;

    branchMapError(k) = mean(abs([xRefL - xCurL; xRefR - xCurR]));
end

regRef = registeredNorm(iRef, :);

for k = 1:nGap
    simpleRms(k) = sqrt(mean((simpleNorm(k, :) - simpleRef) .^ 2));
    registeredRms(k) = sqrt(mean((registeredNorm(k, :) - regRef) .^ 2));
    simpleShift(k) = estimate_shift(xRefGrid, simpleNorm(k, :), simpleRef, fitRange, fitMask);
    registeredShift(k) = estimate_shift(xRefGrid, registeredNorm(k, :), regRef, fitRange, fitMask);
end

%% Export metrics
metrics = table(gapList, simpleShift, registeredShift, ...
    simpleRms, registeredRms, branchMapError, ...
    'VariableNames', {'gap_mm', 'simple_apparent_shift_mm', ...
    'registered_apparent_shift_mm', 'simple_rms_to_ref', ...
    'registered_rms_to_ref', 'mean_control_point_offset_mm'});
writetable(metrics, fullfile(outDir, 'empirical_registration_metrics.csv'));

summary = table( ...
    ["Simple normalization"; "Empirical reference registration"], ...
    [mean(abs(simpleShift)); mean(abs(registeredShift))], ...
    [max(abs(simpleShift)); max(abs(registeredShift))], ...
    [mean(simpleRms); mean(registeredRms)], ...
    'VariableNames', {'method', 'mean_abs_apparent_shift_mm', ...
    'max_abs_apparent_shift_mm', 'mean_rms_to_ref'});
writetable(summary, fullfile(outDir, 'empirical_registration_summary.csv'));

%% Figure 1: waveform collapse
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(xRefGrid - ref.xPeak, simpleNorm(k, :), 'Color', colors(k, :));
end
xlabel('Aligned position x - x_c^{ref} (mm)');
ylabel('Normalized response');
title('(a) Simple normalization');

nexttile; hold on;
for k = 1:nGap
    plot(xRefGrid - ref.xPeak, simpleNorm(k, :) - simpleRef, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Aligned position x - x_c^{ref} (mm)');
ylabel('Residual');
title('(b) Simple residuals');

nexttile; hold on;
for k = 1:nGap
    plot(xRefGrid - ref.xPeak, registeredNorm(k, :), 'Color', colors(k, :));
end
xlabel('Reference-coordinate x - x_c^{ref} (mm)');
ylabel('Registered response');
title('(c) Empirical reference registration');

nexttile; hold on;
for k = 1:nGap
    plot(xRefGrid - ref.xPeak, registeredNorm(k, :) - regRef, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Reference-coordinate x - x_c^{ref} (mm)');
ylabel('Residual');
title('(d) Registration residuals');

exportgraphics(fig1, fullfile(outDir, 'fig1_empirical_registration.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_empirical_registration.pdf'), 'ContentType', 'vector');

%% Figure 2: quantitative comparison
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 9]);
tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(gapList, simpleShift, '-o', 'Color', [0.85, 0.35, 0.10], 'MarkerFaceColor', [0.85, 0.35, 0.10]);
plot(gapList, registeredShift, '-d', 'Color', [0.10, 0.25, 0.50], 'MarkerFaceColor', [0.10, 0.25, 0.50]);
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Clearance g (mm)');
ylabel('Apparent shift (mm)');
legend({'Simple', 'Reference registration'}, 'Location', 'best', 'Box', 'off');
title('(a) False shift');

nexttile; hold on;
plot(gapList, simpleRms, '-o', 'Color', [0.85, 0.35, 0.10], 'MarkerFaceColor', [0.85, 0.35, 0.10]);
plot(gapList, registeredRms, '-d', 'Color', [0.10, 0.25, 0.50], 'MarkerFaceColor', [0.10, 0.25, 0.50]);
xlabel('Clearance g (mm)');
ylabel('RMS to reference');
legend({'Simple', 'Reference registration'}, 'Location', 'best', 'Box', 'off');
title('(b) Collapse residual');

nexttile;
bar(categorical(summary.method), summary.mean_abs_apparent_shift_mm, 0.6, ...
    'FaceColor', [0.5529, 0.6941, 0.8863], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean absolute false shift (mm)');
title('(c) Mean false shift');

nexttile;
bar(categorical(summary.method), summary.mean_rms_to_ref, 0.6, ...
    'FaceColor', [0.8431, 0.8902, 0.7490], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean RMS to reference');
title('(d) Mean collapse residual');

exportgraphics(fig2, fullfile(outDir, 'fig2_empirical_registration_metrics.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_empirical_registration_metrics.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Empirical reference registration validation ===\n');
fprintf('Reference gap: %.1f mm\n', gRef);
fprintf('Control levels used for learning branch-wise maps: %s\n\n', num2str(levelTrain', '%.2f '));
disp(summary);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function [xRefPos, xCurPos] = level_position_pairs(xRef, yRef, iRefPeak, xCur, yCur, iCurPeak, levels, side)
    if strcmpi(side, 'left')
        xRefBranch = xRef(1:iRefPeak);
        yRefBranch = yRef(1:iRefPeak);
        xCurBranch = xCur(1:iCurPeak);
        yCurBranch = yCur(1:iCurPeak);
    else
        xRefBranch = flipud(xRef(iRefPeak:end));
        yRefBranch = flipud(yRef(iRefPeak:end));
        xCurBranch = flipud(xCur(iCurPeak:end));
        yCurBranch = flipud(yCur(iCurPeak:end));
    end

    [yRefUniq, i1] = unique(yRefBranch, 'stable');
    [yCurUniq, i2] = unique(yCurBranch, 'stable');
    xRefUniq = xRefBranch(i1);
    xCurUniq = xCurBranch(i2);

    xRefPos = interp1(yRefUniq, xRefUniq, levels, 'pchip', 'extrap');
    xCurPos = interp1(yCurUniq, xCurUniq, levels, 'pchip', 'extrap');
end

function shift = estimate_shift(x, y, yRef, fitRange, fitMask)
    cost = @(dx) shifted_rmse(x, y, yRef, dx, fitMask);
    shift = fminbnd(cost, fitRange(1), fitRange(2));
end

function val = shifted_rmse(x, y, yRef, dx, fitMask)
    x = x(:);
    y = y(:);
    yRef = yRef(:);
    fitMask = fitMask(:);
    yShift = interp1(x - dx, y, x, 'pchip', 'extrap');
    diffVal = yShift(fitMask) - yRef(fitMask);
    val = sqrt(mean(diffVal .^ 2));
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
