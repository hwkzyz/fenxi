%% Straight blade 2 mm: clearance effect analysis
% Analyze COMSOL-exported waveforms stacked in a single text file.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'straight_blade_2mm_results');
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

%% Load and split stacked curves
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);

commonX = xCell{1};
S = zeros(nGap, numel(commonX));
for k = 1:nGap
    S(k, :) = yCell{k};
end

%% Reference and normalized shapes
gRef = 0.8;
[~, iRef] = min(abs(gapList - gRef));
gRef = gapList(iRef);

Sn = zeros(size(S));
for k = 1:nGap
    Sn(k, :) = normalize_wave(S(k, :));
end
Sref = Sn(iRef, :);

%% Feature extraction
peakVal = zeros(nGap, 1);
baseVal = zeros(nGap, 1);
peakX = zeros(nGap, 1);
fwhmVal = zeros(nGap, 1);
areaVal = zeros(nGap, 1);
shapeRms = zeros(nGap, 1);
dxBias = zeros(nGap, 1);
biasRmse = zeros(nGap, 1);
rho = zeros(nGap, 1);
centerBySym = zeros(nGap, 1);

fitRange = [-0.5, 0.5];
dX = mean(diff(commonX));

for k = 1:nGap
    [peakVal(k), idxPeak] = max(S(k, :));
    baseVal(k) = min(S(k, :));
    peakX(k) = commonX(idxPeak);
    fwhmVal(k) = calc_fwhm(commonX, S(k, :));
    areaVal(k) = trapz(commonX, S(k, :) - baseVal(k));
    shapeRms(k) = sqrt(mean((Sn(k, :) - Sref).^2));
    centerBySym(k) = estimate_center_by_symmetry(commonX, S(k, :));

    cost = @(dx) scaled_template_rmse(S(k, :), interp1(commonX - dx, S(iRef, :), commonX, 'pchip', 'extrap'));
    [dxBias(k), biasRmse(k)] = fminbnd(cost, fitRange(1), fitRange(2));
end

%% Sensitivity independence
dSdx = gradient(S, dX);
dSdg = zeros(size(S));
for k = 1:nGap
    if k == 1
        dSdg(k, :) = (S(k+1, :) - S(k, :)) ./ (gapList(k+1) - gapList(k));
    elseif k == nGap
        dSdg(k, :) = (S(k, :) - S(k-1, :)) ./ (gapList(k) - gapList(k-1));
    else
        dSdg(k, :) = (S(k+1, :) - S(k-1, :)) ./ (gapList(k+1) - gapList(k-1));
    end

    active = S(k, :) > baseVal(k) + 0.10 * range(S(k, :));
    jx = dSdx(k, active);
    jg = dSdg(k, active);
    rho(k) = dot(jx, jg) / max(norm(jx) * norm(jg), eps);
end

%% Export metrics
metrics = table(gapList, peakVal, baseVal, peakX, centerBySym, ...
    fwhmVal, areaVal, shapeRms, dxBias, biasRmse, rho, ...
    'VariableNames', {'gap_mm', 'peak_value', 'baseline', 'peak_x', ...
    'symmetry_center_x', 'fwhm_x', 'area_above_baseline', ...
    'normalized_shape_rms', 'fixed_template_bias_x', 'fixed_template_rmse', ...
    'sensitivity_correlation'});
writetable(metrics, fullfile(outDir, 'straight_blade_2mm_metrics.csv'));

%% Figure 1: waveforms and normalized residuals
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(commonX, S(k, :), 'Color', colors(k, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'eastoutside', 'Box', 'off');
title('(a) Raw waveforms');

nexttile; hold on;
for k = 1:nGap
    plot(commonX, Sn(k, :), 'Color', colors(k, :));
end
xlabel('Position x (mm)');
ylabel('Normalized response');
title('(b) Amplitude-normalized shapes');

nexttile; hold on;
for k = 1:nGap
    plot(commonX, Sn(k, :) - Sref, 'Color', colors(k, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Position x (mm)');
ylabel('Residual');
title(sprintf('(c) Residual to %.1f mm template', gRef));

nexttile; hold on;
bar(gapList, dxBias, 0.65, 'FaceColor', [0.5529, 0.6941, 0.8863], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Clearance g (mm)');
ylabel('Apparent \Deltax (mm)');
title('(d) Fixed-template bias');

exportgraphics(fig1, fullfile(outDir, 'fig1_waveforms_and_bias.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_waveforms_and_bias.pdf'), 'ContentType', 'vector');

%% Figure 2: feature trends
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tiledlayout(fig2, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(gapList, peakVal, '-o', 'Color', colors(2, :), 'MarkerFaceColor', colors(2, :));
xlabel('Clearance g (mm)');
ylabel('Peak value');
title('(a) Peak trend');

nexttile;
plot(gapList, fwhmVal, '-o', 'Color', colors(4, :), 'MarkerFaceColor', colors(4, :));
xlabel('Clearance g (mm)');
ylabel('FWHM (mm)');
title('(b) Width trend');

nexttile;
plot(gapList, areaVal, '-o', 'Color', colors(5, :), 'MarkerFaceColor', colors(5, :));
xlabel('Clearance g (mm)');
ylabel('Area');
title('(c) Area trend');

nexttile;
plot(gapList, peakX, '-o', 'Color', colors(1, :), 'MarkerFaceColor', colors(1, :)); hold on;
plot(gapList, centerBySym, '--s', 'Color', colors(6, :), 'MarkerFaceColor', 'w');
xlabel('Clearance g (mm)');
ylabel('Characteristic position x (mm)');
legend({'Peak position', 'Symmetry center'}, 'Location', 'best', 'Box', 'off');
title('(d) Position shift');

nexttile;
plot(gapList, shapeRms, '-o', 'Color', colors(3, :), 'MarkerFaceColor', colors(3, :));
xlabel('Clearance g (mm)');
ylabel('Shape RMS');
title('(e) Normalized shape residual');

nexttile;
plot(gapList, rho, '-o', 'Color', colors(7, :), 'MarkerFaceColor', colors(7, :));
yline(0, 'k-', 'LineWidth', 0.6);
ylim([-1, 1]);
xlabel('Clearance g (mm)');
ylabel('Correlation');
title('(f) Sensitivity independence');

exportgraphics(fig2, fullfile(outDir, 'fig2_feature_trends.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_feature_trends.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Straight blade 2 mm: clearance analysis ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Loaded %d clearance cases: %s mm\n', nGap, num2str(gapList', '%.1f '));
fprintf('Reference template clearance: %.1f mm\n\n', gRef);
fprintf('%-8s %-10s %-10s %-10s %-10s %-10s\n', ...
    'g(mm)', 'peak', 'FWHM', 'shapeRMS', 'biasDx', 'rho');
for k = 1:nGap
    fprintf('%-8.1f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n', ...
        gapList(k), peakVal(k), fwhmVal(k), shapeRms(k), dxBias(k), rho(k));
end
fprintf('\nResults saved to: %s\n', outDir);

%% Local functions
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
        gapList = (gapVals(1):gapVals(2)-gapVals(1):gapVals(3))';
    else
        gapList = gapVals;
    end

    data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data), 2), :);
    if size(data, 2) < 2
        error('Expected two numeric columns in %s', filePath);
    end
    x = data(:, 1);
    y = data(:, 2);

    cutIdx = find(diff(x) < 0);
    startIdx = [1; cutIdx + 1];
    endIdx = [cutIdx; numel(x)];
    nSeg = numel(startIdx);
    if nSeg ~= numel(gapList)
        error('Segment count (%d) does not match gap count (%d).', nSeg, numel(gapList));
    end

    xCell = cell(nSeg, 1);
    yCell = cell(nSeg, 1);
    for k = 1:nSeg
        idx = startIdx(k):endIdx(k);
        xCell{k} = x(idx);
        yCell{k} = y(idx);
    end
end

function yn = normalize_wave(y)
    yn = (y - min(y)) ./ max(range(y), eps);
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

function center = estimate_center_by_symmetry(x, y)
    [xu, iu] = unique(x);
    yu = y(iu);
    f = @(xx) interp1(xu, yu, xx, 'pchip', 'extrap');

    [~, idxPeak] = max(yu);
    peakGuess = xu(idxPeak);
    cHalfSpan = min(1.0, 0.2 * (max(xu) - min(xu)));
    cGrid = linspace(peakGuess - cHalfSpan, peakGuess + cHalfSpan, 120);
    best = inf;
    center = peakGuess;
    for c = cGrid
        xa = max(min(xu), 2 * c - max(xu));
        xb = min(max(xu), 2 * c - min(xu));
        if xb <= xa
            continue;
        end
        xx = linspace(xa, xb, 150);
        res = f(xx) - f(2 * c - xx);
        val = sum(res .^ 2);
        if val < best
            best = val;
            center = c;
        end
    end
end

function rmse = scaled_template_rmse(yObs, yTpl)
    A = [ones(numel(yTpl), 1), yTpl(:)];
    theta = A \ yObs(:);
    yFit = A * theta;
    rmse = sqrt(mean((yObs(:) - yFit) .^ 2));
end
