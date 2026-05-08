%% Physical kernel decoupling for straight-blade 2 mm waveforms
% Idea:
% 1) The blade geometry is fixed and represented by a rectangular occupancy.
% 2) Different clearances mainly change the sensor field kernel width sigma(g),
%    plus baseline and gain.
% 3) After fitting sigma(g), every waveform can be mapped to a common
%    reference-gap kernel instead of using amplitude normalization alone.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'physical_kernel_decoupling_results');
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
x = xCell{1}(:);
S = zeros(nGap, numel(x));
for k = 1:nGap
    S(k, :) = yCell{k}(:)';
end
dx = mean(diff(x));

%% Fixed geometry width from file meaning: straight blade 2 mm
W_geom = 2.0; % mm

%% Fit each clearance with a physically motivated kernel model
fitRes = repmat(struct('gap', [], 'C0', [], 'A', [], 'sigma', [], 'xc', [], ...
    'yhat', [], 'rmse', [], 'r2', []), nGap, 1);

opts = optimoptions('lsqcurvefit', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 2e4);

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
    yHat = modelFun(pHat, x);
    rmse = sqrt(mean((y - yHat) .^ 2));
    r2 = 1 - sum((y - yHat) .^ 2) / sum((y - mean(y)) .^ 2);

    fitRes(k).gap = gapList(k);
    fitRes(k).C0 = pHat(1);
    fitRes(k).A = pHat(2);
    fitRes(k).sigma = pHat(3);
    fitRes(k).xc = pHat(4);
    fitRes(k).yhat = yHat;
    fitRes(k).rmse = rmse;
    fitRes(k).r2 = r2;
end

sigmaList = [fitRes.sigma]';
AList = [fitRes.A]';
C0List = [fitRes.C0]';
xcList = [fitRes.xc]';
rmseList = [fitRes.rmse]';
r2List = [fitRes.r2]';

%% Build two normalizations for comparison
% 1) Simple amplitude normalization
simpleNorm = zeros(size(S));
for k = 1:nGap
    simpleNorm(k, :) = normalize_wave(S(k, :));
end

% 2) Physical kernel normalization to the widest kernel
[sigmaRef, iRef] = max(sigmaList);
gRef = gapList(iRef);
physNorm = zeros(size(S));

for k = 1:nGap
    u = x - xcList(k);
    y0 = (S(k, :)' - C0List(k)) / max(AList(k), eps);
    y0 = y0(:)';

    if sigmaRef >= sigmaList(k)
        sigmaAdd = sqrt(max(0, sigmaRef^2 - sigmaList(k)^2));
        yEq = gaussian_smooth(y0, dx, sigmaAdd);
    else
        yEq = y0;
    end
    physNorm(k, :) = yEq;
end

refSimple = simpleNorm(iRef, :);
refPhys = physNorm(iRef, :);
simpleSpread = sqrt(mean((simpleNorm - refSimple) .^ 2, 2));
physSpread = sqrt(mean((physNorm - refPhys) .^ 2, 2));

%% Export metrics
metrics = table(gapList, C0List, AList, sigmaList, xcList, rmseList, r2List, ...
    simpleSpread, physSpread, ...
    'VariableNames', {'gap_mm', 'baseline', 'gain', 'sigma_mm', ...
    'center_shift_mm', 'fit_rmse', 'fit_r2', ...
    'simple_norm_rms_to_ref', 'physical_norm_rms_to_ref'});
writetable(metrics, fullfile(outDir, 'physical_kernel_fit_metrics.csv'));

%% Figure 1: model fit and sigma(g)
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(x, S(k, :), 'o', 'Color', colors(k, :), 'MarkerSize', 3.0, ...
        'HandleVisibility', 'off');
    plot(x, fitRes(k).yhat, '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'eastoutside', 'Box', 'off');
title('(a) Physical kernel model fit');

nexttile; hold on;
plot(gapList, sigmaList, '-o', 'Color', [0.16, 0.48, 0.66], ...
    'MarkerFaceColor', [0.16, 0.48, 0.66]);
xlabel('Clearance g (mm)');
ylabel('\sigma_g (mm)');
title('(b) Fitted field-kernel width');

nexttile; hold on;
plot(gapList, AList, '-o', 'Color', [0.85, 0.35, 0.10], ...
    'MarkerFaceColor', [0.85, 0.35, 0.10]);
xlabel('Clearance g (mm)');
ylabel('Gain A(g)');
title('(c) Fitted gain');

nexttile; hold on;
plot(gapList, rmseList, '-o', 'Color', [0.20, 0.60, 0.60], ...
    'MarkerFaceColor', [0.20, 0.60, 0.60]);
xlabel('Clearance g (mm)');
ylabel('RMSE (pF)');
title('(d) Fit error');

exportgraphics(fig1, fullfile(outDir, 'fig1_physical_fit.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_physical_fit.pdf'), 'ContentType', 'vector');

%% Figure 2: simple normalization vs physical-kernel normalization
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(x, simpleNorm(k, :), 'Color', colors(k, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
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
title(sprintf('(b) Residual to %.1f mm template', gRef));

nexttile; hold on;
for k = 1:nGap
    plot(x - xcList(k), physNorm(k, :), 'Color', colors(k, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
end
xlabel('Aligned position x - x_c (mm)');
ylabel('Kernel-normalized response');
title('(c) Physical kernel normalization');

nexttile; hold on;
for k = 1:nGap
    plot(x - xcList(k), physNorm(k, :) - refPhys, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Aligned position x - x_c (mm)');
ylabel('Residual');
title(sprintf('(d) Residual after kernel mapping to %.1f mm', gRef));

exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_comparison.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_normalization_comparison.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Physical kernel decoupling: straight blade 2 mm ===\n');
fprintf('Reference gap for kernel mapping: %.1f mm\n', gRef);
fprintf('Geometry width fixed at W = %.2f mm\n\n', W_geom);
fprintf('%-8s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'g(mm)', 'sigma', 'gain', 'x_c', 'RMSE', 'simpleRMS', 'physRMS');
for k = 1:nGap
    fprintf('%-8.1f %-10.4f %-10.4f %-10.4f %-10.4g %-10.4f %-10.4f\n', ...
        gapList(k), sigmaList(k), AList(k), xcList(k), rmseList(k), ...
        simpleSpread(k), physSpread(k));
end

fprintf('\nAverage RMS spread:\n');
fprintf('Simple amplitude normalization: %.5f\n', mean(simpleSpread));
fprintf('Physical kernel normalization: %.5f\n', mean(physSpread));
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function y = blurred_rect_model(x, C0, A, sigma, xc, W)
    u = x - xc;
    shape = 0.5 * (erf((u + W / 2) ./ (sqrt(2) * sigma)) - ...
                   erf((u - W / 2) ./ (sqrt(2) * sigma)));
    y = C0 + A * shape;
end

function y = normalize_wave(x)
    y = (x - min(x)) ./ max(range(x), eps);
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
