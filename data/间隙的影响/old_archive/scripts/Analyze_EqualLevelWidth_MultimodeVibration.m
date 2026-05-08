%% Analyze equal-level width under user-specified multimode vibration
% The waveform amplitude itself is not drifted vertically. We only impose
% circumferential vibration through the coordinate mapping xi = x0 - u(t).
% The script compares:
%   1) no-vibration equal-capacitance widths;
%   2) apparent widths measured in rigid coordinate x0;
%   3) true widths measured in relative coordinate xi.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'multimode_width_analysis_results');
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

%% User-style vibration settings
u0 = 0.0;                     % mm
A_vib = [1, 1];           % mm
f_vib = [500, 1300];          % Hz
phi0 = [pi/4, -pi/3];         % rad
V_tip = 3.0e5;                % mm/s, representative circumferential speed
levels = (0.10:0.10:0.90)';   % absolute equal-capacitance levels via no-vibration normalization

%% Load baseline waveforms and pick one representative no-vibration pulse
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
[~, refIdx] = min(abs(gapList - 0.8));
xRef = xCell{refIdx}(:);
yRef = yCell{refIdx}(:);
[baseRef, ampRef, ~] = normalize_waveform(yRef);
yTargets = baseRef + ampRef * levels;
[wRef, xLeftRef, xRightRef] = level_widths_abs(xRef, yRef, yTargets);

%% Build a rigid-coordinate grid and impose multimode vibration
% Keep xi inside the original no-vibration support.
uMax = abs(u0) + sum(abs(A_vib));
x0Grid = linspace(min(xRef) + uMax, max(xRef) - uMax, numel(xRef))';
tGrid = (x0Grid - mean(x0Grid)) / V_tip;
uGrid = u0 + multimode_vibration(tGrid, A_vib, f_vib, phi0);
xiGrid = x0Grid - uGrid;

if any(diff(xiGrid) <= 0)
    error('The chosen vibration and tip speed make xi(x0) non-monotone. Increase V_tip or reduce amplitude.');
end

yObs = interp1(xRef, yRef, xiGrid, 'pchip', 'extrap');
[wX0, xLeftX0, xRightX0] = level_widths_abs(x0Grid, yObs, yTargets);
[xiSorted, sortIdx] = sort(xiGrid, 'ascend');
yXi = yObs(sortIdx);
[wXi, xLeftXi, xRightXi] = level_widths_abs(xiSorted, yXi, yTargets);

%% Analyze why the apparent width changes
deltaU = nan(size(levels));
for k = 1:numel(levels)
    if all(isfinite([xLeftX0(k), xRightX0(k)]))
        tLeft = (xLeftX0(k) - mean(x0Grid)) / V_tip;
        tRight = (xRightX0(k) - mean(x0Grid)) / V_tip;
        deltaU(k) = multimode_vibration(tRight, A_vib, f_vib, phi0) - ...
                    multimode_vibration(tLeft, A_vib, f_vib, phi0);
    end
end

checkResidual = wX0 - (wRef + deltaU);

summaryTable = table(levels, yTargets, wRef, wX0, wXi, deltaU, ...
    wX0 - wRef, wXi - wRef, checkResidual, ...
    'VariableNames', {'level_q', 'capacitance_level', 'width_ref_mm', ...
    'width_x0_mm', 'width_xi_mm', 'delta_u_mm', ...
    'apparent_width_error_mm', 'true_width_error_mm', 'identity_residual_mm'});
writetable(summaryTable, fullfile(outDir, 'multimode_equal_level_width_table.csv'));

globalSummary = table( ...
    mean(abs(summaryTable.apparent_width_error_mm), 'omitnan'), ...
    max(abs(summaryTable.apparent_width_error_mm), [], 'omitnan'), ...
    mean(abs(summaryTable.true_width_error_mm), 'omitnan'), ...
    max(abs(summaryTable.true_width_error_mm), [], 'omitnan'), ...
    mean(abs(summaryTable.identity_residual_mm), 'omitnan'), ...
    max(abs(summaryTable.identity_residual_mm), [], 'omitnan'), ...
    'VariableNames', {'mean_abs_err_x0_mm', 'max_abs_err_x0_mm', ...
    'mean_abs_err_xi_mm', 'max_abs_err_xi_mm', ...
    'mean_abs_identity_residual_mm', 'max_abs_identity_residual_mm'});
writetable(globalSummary, fullfile(outDir, 'multimode_equal_level_width_summary.csv'));

%% Figure 1: waveform and widths
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(xRef, yRef, 'k-', 'DisplayName', 'No vibration');
plot(x0Grid, yObs, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Observed in x_0');
xlabel('Rigid coordinate x_0 (mm)');
ylabel('Capacitance (pF)');
title('(a) Waveform in rigid coordinate');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
plot(xRef, yRef, 'k-', 'DisplayName', 'No vibration');
plot(xiSorted, yXi, '--', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Remapped to \xi');
xlabel('Relative coordinate \xi (mm)');
ylabel('Capacitance (pF)');
title('(b) Waveform in relative coordinate');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, wRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference');
plot(levels, wX0, '-s', 'MarkerSize', 4, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Measured in x_0');
xlabel('Normalized level q');
ylabel('Equal-capacitance width (mm)');
title('(c) Apparent width change');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, wRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference');
plot(levels, wXi, '-s', 'MarkerSize', 4, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Measured in \xi');
xlabel('Normalized level q');
ylabel('Equal-capacitance width (mm)');
title('(d) True width invariance');
legend('Location', 'best', 'FontSize', 7.2);

title(tl, sprintf('Equal-level width under multimode vibration at g = %.1f mm', gapList(refIdx)), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_multimode_width_comparison.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_multimode_width_comparison.pdf'), 'ContentType', 'vector');

%% Figure 2: why width changes in x0
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(levels, summaryTable.apparent_width_error_mm, '-o', 'MarkerSize', 4, ...
    'Color', [0.10 0.45 0.75], 'DisplayName', 'w_{x_0} - w_{ref}');
plot(levels, summaryTable.delta_u_mm, '--s', 'MarkerSize', 4, ...
    'Color', [0.85 0.35 0.10], 'DisplayName', 'u(t^+) - u(t^-)');
xlabel('Normalized level q');
ylabel('Width change / displacement difference (mm)');
title('(a) Source of apparent width change');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, summaryTable.identity_residual_mm, '-d', 'MarkerSize', 4, ...
    'Color', [0.20 0.60 0.40]);
yline(0, 'k--');
xlabel('Normalized level q');
ylabel('Residual (mm)');
title('(b) Check: w_{x_0} - w_{ref} - \Delta u');
exportgraphics(fig2, fullfile(outDir, 'fig2_multimode_width_identity.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_multimode_width_identity.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Multimode equal-level width analysis ===\n');
fprintf('Reference gap: %.1f mm\n', gapList(refIdx));
fprintf('Tip speed V = %.3e mm/s\n', V_tip);
fprintf('u(t) = %.3f + sum A_m sin(2*pi*f_m*t + phi_m)\n', u0);
fprintf('A = [%s] mm\n', num2str(A_vib, '%.3f '));
fprintf('f = [%s] Hz\n', num2str(f_vib, '%.1f '));
fprintf('phi = [%s] rad\n\n', num2str(phi0, '%.3f '));
disp(globalSummary);
fprintf('\nLevel-by-level comparison:\n');
disp(summaryTable(:, {'level_q', 'width_ref_mm', 'width_x0_mm', 'width_xi_mm', ...
    'delta_u_mm', 'apparent_width_error_mm', 'true_width_error_mm'}));
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function [gapList, xCell, yCell] = load_stacked_comsol_curves(filePath)
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

function [widths, xLeft, xRight] = level_widths_abs(x, y, yTargets)
    x = x(:);
    y = y(:);
    yTargets = yTargets(:);

    [~, iPeak] = max(y);
    xLeft = nan(size(yTargets));
    xRight = nan(size(yTargets));

    for k = 1:numel(yTargets)
        q = yTargets(k);
        xLeft(k) = find_crossing_left(x(1:iPeak), y(1:iPeak), q);
        xRight(k) = find_crossing_right(x(iPeak:end), y(iPeak:end), q);
    end
    widths = xRight - xLeft;
end

function xq = find_crossing_left(x, y, q)
    idx = find(y >= q, 1, 'first');
    if isempty(idx)
        xq = NaN;
    elseif idx == 1
        xq = x(1);
    else
        xq = interp1(y(idx-1:idx), x(idx-1:idx), q, 'linear');
    end
end

function xq = find_crossing_right(x, y, q)
    idx = find(y < q, 1, 'first');
    if isempty(idx)
        xq = x(end);
    elseif idx == 1
        xq = x(1);
    else
        xq = interp1(y(idx-1:idx), x(idx-1:idx), q, 'linear');
    end
end
