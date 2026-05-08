%% Validate peak-centered symmetrization under circumferential vibration
% We test whether a vibration-distorted waveform can be collapsed back onto
% the no-vibration waveform by:
%   1) centering both waveforms at their peaks;
%   2) symmetrizing the vibrating waveform by left-right mirror averaging.
% The imposed vibration only acts through the horizontal coordinate:
%   y_obs(x0) = y_ref(x0 - u(t))
% with no vertical amplitude drift.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'peak_symmetrization_results');
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
    'DefaultLineLineWidth', 1.2, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% User-specified circumferential vibration
u0 = 0.0;                   % mm
A_vib = [0.5, 0.4];         % mm
f_vib = [500, 1300];        % Hz
phi0 = [pi/4, -pi/3];       % rad
V_tip = 3.0e5;              % mm/s
levels = (0.10:0.10:0.90)';

%% Load one no-vibration waveform as the baseline
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
[~, refIdx] = min(abs(gapList - 0.8));
xRef = xCell{refIdx}(:);
yRef = yCell{refIdx}(:);

% Use a reduced rigid-coordinate domain so the vibrating coordinate stays
% within the support of the no-vibration waveform.
uMax = abs(u0) + sum(abs(A_vib));
x0Grid = linspace(min(xRef) + uMax, max(xRef) - uMax, numel(xRef))';
tGrid = (x0Grid - mean(x0Grid)) / V_tip;
uGrid = u0 + multimode_vibration(tGrid, A_vib, f_vib, phi0);
xiGrid = x0Grid - uGrid;
if any(diff(xiGrid) <= 0)
    error('xi(x0) is non-monotone. Increase V_tip or reduce vibration amplitude.');
end

yObs = interp1(xRef, yRef, xiGrid, 'pchip', 'extrap');

%% Peak-centered coordinates
[~, idxPeakRef] = max(yRef);
[~, idxPeakObs] = max(yObs);
xPeakRef = xRef(idxPeakRef);
xPeakObs = x0Grid(idxPeakObs);

zRef = xRef - xPeakRef;
zObs = x0Grid - xPeakObs;

zMin = max(min(zRef), min(zObs));
zMax = min(max(zRef), max(zObs));
zCommon = linspace(zMin, zMax, 2001)';

yRefPk = interp1(zRef, yRef, zCommon, 'pchip');
yObsPk = interp1(zObs, yObs, zCommon, 'pchip');

%% Symmetrize the vibrating waveform around its peak
[zSym, yObsSym] = symmetric_average(zCommon, yObsPk);
[~, yRefSym] = symmetric_average(zCommon, yRefPk);

%% Compare overlap before and after symmetrization
valid0 = isfinite(yRefPk) & isfinite(yObsPk);
validS = isfinite(yRefSym) & isfinite(yObsSym);
rmse_peak_centered = sqrt(mean((yObsPk(valid0) - yRefPk(valid0)).^2));
rmse_after_sym = sqrt(mean((yObsSym(validS) - yRefSym(validS)).^2));

%% Compare equal-capacitance widths
[baseRef, ampRef, ~] = normalize_waveform(yRef);
yTargets = baseRef + ampRef * levels;
[wRef, xLRef, xRRef] = level_widths_abs(zCommon, yRefPk, yTargets);
[wObsPk, xLObsPk, xRObsPk] = level_widths_abs(zCommon, yObsPk, yTargets);
[wObsSym, xLObsSym, xRObsSym] = level_widths_abs(zSym, yObsSym, yTargets);

summaryTable = table(levels, yTargets, wRef, wObsPk, wObsSym, ...
    wObsPk - wRef, wObsSym - wRef, ...
    'VariableNames', {'level_q', 'capacitance_level', 'width_ref_mm', ...
    'width_peakcentered_mm', 'width_symmetrized_mm', ...
    'peakcentered_width_error_mm', 'symmetrized_width_error_mm'});
writetable(summaryTable, fullfile(outDir, 'peak_symmetrization_width_table.csv'));

globalSummary = table( ...
    rmse_peak_centered, ...
    rmse_after_sym, ...
    mean(abs(summaryTable.peakcentered_width_error_mm), 'omitnan'), ...
    max(abs(summaryTable.peakcentered_width_error_mm), [], 'omitnan'), ...
    mean(abs(summaryTable.symmetrized_width_error_mm), 'omitnan'), ...
    max(abs(summaryTable.symmetrized_width_error_mm), [], 'omitnan'), ...
    'VariableNames', {'rmse_peak_centered', 'rmse_after_symmetrization', ...
    'mean_abs_width_err_peakcentered_mm', 'max_abs_width_err_peakcentered_mm', ...
    'mean_abs_width_err_symmetrized_mm', 'max_abs_width_err_symmetrized_mm'});
writetable(globalSummary, fullfile(outDir, 'peak_symmetrization_summary.csv'));

%% Figure 1: overlap before and after symmetrization
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(zCommon, yRefPk, 'k-', 'DisplayName', 'No vibration, peak-centered');
plot(zCommon, yObsPk, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Vibration, peak-centered');
xlabel('Peak-centered coordinate z (mm)');
ylabel('Capacitance (pF)');
title('(a) Before symmetrization');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(zSym, yRefSym, 'k-', 'DisplayName', 'No vibration, symmetric mean');
plot(zSym, yObsSym, '-', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Vibration, symmetric mean');
xlabel('Peak-centered coordinate z (mm)');
ylabel('Capacitance (pF)');
title('(b) After symmetrization');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(zCommon, yObsPk - yRefPk, '-', 'Color', [0.10 0.45 0.75], ...
    'DisplayName', sprintf('Before, RMSE = %.4g', rmse_peak_centered));
yline(0, 'k--');
xlabel('Peak-centered coordinate z (mm)');
ylabel('Residual (pF)');
title('(c) Residual before symmetrization');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(zSym, yObsSym - yRefSym, '-', 'Color', [0.85 0.35 0.10], ...
    'DisplayName', sprintf('After, RMSE = %.4g', rmse_after_sym));
yline(0, 'k--');
xlabel('Peak-centered coordinate z (mm)');
ylabel('Residual (pF)');
title('(d) Residual after symmetrization');
legend('Location', 'best', 'FontSize', 7.2);

title(tl, sprintf('Peak-centered symmetrization validation at g = %.1f mm', gapList(refIdx)), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_peak_symmetrization_overlap.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_peak_symmetrization_overlap.pdf'), 'ContentType', 'vector');

%% Figure 2: width comparison
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(levels, wRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference');
plot(levels, wObsPk, '-s', 'MarkerSize', 4, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Peak-centered only');
plot(levels, wObsSym, '-d', 'MarkerSize', 4, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Peak + symmetrized');
xlabel('Normalized level q');
ylabel('Equal-capacitance width (mm)');
title('(a) Width comparison');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, summaryTable.peakcentered_width_error_mm, '-s', 'MarkerSize', 4, ...
    'Color', [0.10 0.45 0.75], 'DisplayName', 'Peak-centered only');
plot(levels, summaryTable.symmetrized_width_error_mm, '-d', 'MarkerSize', 4, ...
    'Color', [0.85 0.35 0.10], 'DisplayName', 'Peak + symmetrized');
yline(0, 'k--');
xlabel('Normalized level q');
ylabel('Width error relative to reference (mm)');
title('(b) Width error');
legend('Location', 'best', 'FontSize', 7.2);

exportgraphics(fig2, fullfile(outDir, 'fig2_peak_symmetrization_widths.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_peak_symmetrization_widths.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Peak-centered symmetrization validation ===\n');
fprintf('Reference gap: %.1f mm\n', gapList(refIdx));
fprintf('u(t) = %.3f + sum A_m sin(2*pi*f_m*t + phi_m)\n', u0);
fprintf('A = [%s] mm\n', num2str(A_vib, '%.3f '));
fprintf('f = [%s] Hz\n', num2str(f_vib, '%.1f '));
fprintf('phi = [%s] rad\n', num2str(phi0, '%.3f '));
fprintf('Tip speed V = %.3e mm/s\n\n', V_tip);
disp(globalSummary);
fprintf('\nLevel-by-level width comparison:\n');
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

function u = multimode_vibration(t, A, f, phi)
    u = zeros(size(t));
    for k = 1:numel(A)
        u = u + A(k) * sin(2 * pi * f(k) * t + phi(k));
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function [zHalf, ySym] = symmetric_average(z, y)
    z = z(:);
    y = y(:);
    zHalfMax = min(max(z), -min(z));
    zHalf = linspace(0, zHalfMax, floor(numel(z) / 2))';
    yPlus = interp1(z, y, zHalf, 'pchip', NaN);
    yMinus = interp1(z, y, -zHalf, 'pchip', NaN);
    yMeanHalf = 0.5 * (yPlus + yMinus);
    zHalf = zHalf(isfinite(yMeanHalf));
    yMeanHalf = yMeanHalf(isfinite(yMeanHalf));

    zHalf = zHalf(:);
    yMeanHalf = yMeanHalf(:);
    if zHalf(1) == 0
        zHalfPos = zHalf(2:end);
        yHalfPos = yMeanHalf(2:end);
    else
        zHalfPos = zHalf;
        yHalfPos = yMeanHalf;
    end

    zHalfNeg = -flipud(zHalfPos);
    yHalfNeg = flipud(yHalfPos);
    zHalf = [zHalfNeg; 0; zHalfPos];
    ySym = [yHalfNeg; yMeanHalf(1); yHalfPos];
end

function [widths, xLeft, xRight] = level_widths_abs(x, y, yTargets)
    x = x(:);
    y = y(:);
    yTargets = yTargets(:);
    [~, iPeak] = max(y);

    xLeft = nan(size(yTargets));
    xRight = nan(size(yTargets));
    widths = nan(size(yTargets));

    xRise = x(1:iPeak);
    yRise = y(1:iPeak);
    xFall = x(iPeak:end);
    yFall = y(iPeak:end);

    [yRiseU, idxRiseU] = unique(yRise, 'stable');
    xRiseU = xRise(idxRiseU);
    [yFallU, idxFallU] = unique(yFall, 'stable');
    xFallU = xFall(idxFallU);

    for k = 1:numel(yTargets)
        yt = yTargets(k);
        if yt >= min(yRiseU) && yt <= max(yRiseU) && yt >= min(yFallU) && yt <= max(yFallU)
            xLeft(k) = interp1(yRiseU, xRiseU, yt, 'linear');
            xRight(k) = interp1(flipud(yFallU), flipud(xFallU), yt, 'linear');
            widths(k) = xRight(k) - xLeft(k);
        end
    end
end
