%% Validate no-vibration geometric decoupling idea
% Validation logic:
% 1) All waveforms are known to be no-vibration cases, so the true shift is 0.
% 2) Compare several preprocessing strategies:
%    a) simple amplitude normalization
%    b) width-aligned normalization
%    c) physical-kernel normalization
% 3) The best decoupling should minimize apparent shift and waveform spread.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'validate_no_vibration_results');
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

%% Load stacked curves
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);
x = xCell{1}(:);
dx = mean(diff(x));
S = zeros(nGap, numel(x));
for k = 1:nGap
    S(k, :) = yCell{k}(:)';
end

% Use largest clearance as the reference because its kernel is widest.
[~, iRef] = max(gapList);
gRef = gapList(iRef);
refRaw = S(iRef, :);

%% Step 1: simple amplitude normalization
simpleNorm = zeros(size(S));
for k = 1:nGap
    simpleNorm(k, :) = normalize_wave(S(k, :));
end
refSimple = simpleNorm(iRef, :);

%% Step 2: width-aligned normalization
fwhmList = zeros(nGap, 1);
peakPos = zeros(nGap, 1);
widthNorm = zeros(size(S));

for k = 1:nGap
    fwhmList(k) = calc_fwhm(x, S(k, :));
    [~, idxPeak] = max(S(k, :));
    peakPos(k) = x(idxPeak);
end

fwhmRef = fwhmList(iRef);
peakRef = peakPos(iRef);

for k = 1:nGap
    y0 = normalize_wave(S(k, :));
    xWarp = (x - peakPos(k)) * (fwhmList(k) / fwhmRef) + peakRef;
    widthNorm(k, :) = interp1(xWarp, y0, x, 'pchip', 'extrap');
end
refWidth = widthNorm(iRef, :);

%% Step 3: physical-kernel normalization
W_geom = 2.0; % mm
opts = optimoptions('lsqcurvefit', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 2e4);

sigmaList = zeros(nGap, 1);
AList = zeros(nGap, 1);
C0List = zeros(nGap, 1);
xcList = zeros(nGap, 1);
fitRmse = zeros(nGap, 1);
physNorm = zeros(size(S));

for k = 1:nGap
    y = S(k, :)';
    yMin = min(y);
    yMax = max(y);
    [~, idxPeak] = max(y);
    p0 = [yMin, yMax - yMin, 1.2, x(idxPeak)];
    lb = [yMin - 0.1 * range(y), 0, 0.05, -1.0];
    ub = [yMax, 5 * (yMax - yMin + eps), 5.0, 1.0];

    modelFun = @(p, xx) blurred_rect_model(xx, p(1), p(2), p(3), p(4), W_geom);
    pHat = lsqcurvefit(modelFun, p0, x, y, lb, ub, opts);

    C0List(k) = pHat(1);
    AList(k) = pHat(2);
    sigmaList(k) = pHat(3);
    xcList(k) = pHat(4);
    fitRmse(k) = sqrt(mean((y - modelFun(pHat, x)) .^ 2));
end

sigmaRef = sigmaList(iRef);
for k = 1:nGap
    y0 = (S(k, :) - C0List(k)) / max(AList(k), eps);
    y0 = y0(:)';
    sigmaAdd = sqrt(max(0, sigmaRef^2 - sigmaList(k)^2));
    yEq = gaussian_smooth(y0, dx, sigmaAdd);
    physNorm(k, :) = yEq;
end
refPhys = physNorm(iRef, :);

%% Quantitative validation: apparent shift and residual spread
simpleShift = zeros(nGap, 1);
widthShift = zeros(nGap, 1);
physShift = zeros(nGap, 1);

simpleRms = zeros(nGap, 1);
widthRms = zeros(nGap, 1);
physRms = zeros(nGap, 1);

fitRange = [-0.5, 0.5];
fitMask = x >= min(x) + 0.4 & x <= max(x) - 0.4;

for k = 1:nGap
    simpleShift(k) = estimate_shift(x, simpleNorm(k, :), refSimple, fitRange, fitMask);
    widthShift(k) = estimate_shift(x, widthNorm(k, :), refWidth, fitRange, fitMask);
    physShift(k) = estimate_shift(x, physNorm(k, :), refPhys, fitRange, fitMask);

    simpleRms(k) = sqrt(mean((simpleNorm(k, :) - refSimple) .^ 2));
    widthRms(k) = sqrt(mean((widthNorm(k, :) - refWidth) .^ 2));
    physRms(k) = sqrt(mean((physNorm(k, :) - refPhys) .^ 2));
end

%% Export metrics
metrics = table(gapList, fwhmList, peakPos, sigmaList, fitRmse, ...
    simpleShift, widthShift, physShift, simpleRms, widthRms, physRms, ...
    'VariableNames', {'gap_mm', 'fwhm_mm', 'peak_position_mm', 'sigma_mm', ...
    'physical_fit_rmse', 'simple_apparent_shift_mm', ...
    'width_aligned_apparent_shift_mm', 'physical_apparent_shift_mm', ...
    'simple_rms_to_ref', 'width_aligned_rms_to_ref', 'physical_rms_to_ref'});
writetable(metrics, fullfile(outDir, 'no_vibration_validation_metrics.csv'));

summary = table( ...
    ["Simple amplitude normalization"; "Width-aligned normalization"; "Physical kernel normalization"], ...
    [mean(abs(simpleShift)); mean(abs(widthShift)); mean(abs(physShift))], ...
    [max(abs(simpleShift)); max(abs(widthShift)); max(abs(physShift))], ...
    [mean(simpleRms); mean(widthRms); mean(physRms)], ...
    'VariableNames', {'method', 'mean_abs_apparent_shift_mm', ...
    'max_abs_apparent_shift_mm', 'mean_rms_to_ref'});
writetable(summary, fullfile(outDir, 'no_vibration_validation_summary.csv'));

%% Figure 1: waveform collapse comparison
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tiledlayout(fig1, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(x, simpleNorm(k, :), 'Color', colors(k, :));
end
xlabel('Position x (mm)');
ylabel('Normalized response');
title('(a) Simple amplitude normalization');

nexttile; hold on;
for k = 1:nGap
    plot(x, simpleNorm(k, :) - refSimple, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Position x (mm)');
ylabel('Residual');
title('(b) Simple residuals');

nexttile; hold on;
for k = 1:nGap
    plot(x, widthNorm(k, :), 'Color', colors(k, :));
end
xlabel('Position x (mm)');
ylabel('Width-aligned response');
title('(c) Width-aligned normalization');

nexttile; hold on;
for k = 1:nGap
    plot(x, widthNorm(k, :) - refWidth, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Position x (mm)');
ylabel('Residual');
title('(d) Width-aligned residuals');

nexttile; hold on;
for k = 1:nGap
    plot(x - xcList(k), physNorm(k, :), 'Color', colors(k, :));
end
xlabel('Aligned position x - x_c (mm)');
ylabel('Kernel-normalized response');
title('(e) Physical kernel normalization');

nexttile; hold on;
for k = 1:nGap
    plot(x - xcList(k), physNorm(k, :) - refPhys, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Aligned position x - x_c (mm)');
ylabel('Residual');
title('(f) Physical residuals');

exportgraphics(fig1, fullfile(outDir, 'fig1_waveform_collapse_comparison.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_waveform_collapse_comparison.pdf'), 'ContentType', 'vector');

%% Figure 2: false shift comparison
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 9]);
tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(gapList, simpleShift, '-o', 'Color', [0.85, 0.35, 0.10], 'MarkerFaceColor', [0.85, 0.35, 0.10]);
plot(gapList, widthShift, '-s', 'Color', [0.20, 0.60, 0.60], 'MarkerFaceColor', [0.20, 0.60, 0.60]);
plot(gapList, physShift, '-d', 'Color', [0.10, 0.25, 0.50], 'MarkerFaceColor', [0.10, 0.25, 0.50]);
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Clearance g (mm)');
ylabel('Apparent shift (mm)');
legend({'Simple', 'Width-aligned', 'Physical kernel'}, 'Location', 'best', 'Box', 'off');
title('(a) False shift under no-vibration condition');

nexttile; hold on;
plot(gapList, simpleRms, '-o', 'Color', [0.85, 0.35, 0.10], 'MarkerFaceColor', [0.85, 0.35, 0.10]);
plot(gapList, widthRms, '-s', 'Color', [0.20, 0.60, 0.60], 'MarkerFaceColor', [0.20, 0.60, 0.60]);
plot(gapList, physRms, '-d', 'Color', [0.10, 0.25, 0.50], 'MarkerFaceColor', [0.10, 0.25, 0.50]);
xlabel('Clearance g (mm)');
ylabel('RMS to reference');
legend({'Simple', 'Width-aligned', 'Physical kernel'}, 'Location', 'best', 'Box', 'off');
title('(b) Waveform spread after compensation');

nexttile;
bar(categorical(summary.method), summary.mean_abs_apparent_shift_mm, 0.6, ...
    'FaceColor', [0.5529, 0.6941, 0.8863], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean absolute false shift (mm)');
title('(c) Mean false shift by method');

nexttile;
bar(categorical(summary.method), summary.mean_rms_to_ref, 0.6, ...
    'FaceColor', [0.8431, 0.8902, 0.7490], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean RMS to reference');
title('(d) Mean residual spread by method');

exportgraphics(fig2, fullfile(outDir, 'fig2_false_shift_comparison.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_false_shift_comparison.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== No-vibration validation of geometric decoupling ===\n');
fprintf('Ground truth: all waveforms have zero vibration-induced shift.\n');
fprintf('Reference clearance: %.1f mm\n\n', gRef);
disp(summary);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function y = normalize_wave(x)
    y = (x - min(x)) ./ max(range(x), eps);
end

function width = calc_fwhm(x, y)
    y0 = min(y);
    halfLevel = y0 + 0.5 * range(y);
    above = y >= halfLevel;
    idx = find(diff(above) ~= 0);
    if numel(idx) < 2
        width = NaN;
        return;
    end
    xCross = zeros(numel(idx), 1);
    for ii = 1:numel(idx)
        i = idx(ii);
        xCross(ii) = interp1(y(i:i+1), x(i:i+1), halfLevel, 'linear', 'extrap');
    end
    width = xCross(end) - xCross(1);
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

function y = blurred_rect_model(x, C0, A, sigma, xc, W)
    u = x - xc;
    shape = 0.5 * (erf((u + W / 2) ./ (sqrt(2) * sigma)) - ...
                   erf((u - W / 2) ./ (sqrt(2) * sigma)));
    y = C0 + A * shape;
end

function yOut = gaussian_smooth(yIn, dx, sigma)
    if sigma < 1e-9
        yOut = yIn;
        return;
    end
    halfWidth = max(3, ceil(4 * sigma / dx));
    xx = (-halfWidth:halfWidth) * dx;
    kernel = exp(-0.5 * (xx / sigma) .^ 2);
    kernel = kernel / sum(kernel);
    yOut = conv(yIn, kernel, 'same');
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
