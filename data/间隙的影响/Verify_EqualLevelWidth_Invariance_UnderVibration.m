%% Verify equal-level width invariance under synthetic vibration
% This script uses no-vibration calibration waveforms and imposes several
% synthetic vibration fields u(x0). It compares:
%   1) true equal-level widths in the relative coordinate xi = x0 - u;
%   2) apparent equal-level widths measured directly in the rigid coordinate x0.
%
% Conclusion target:
% Equal-level width is strictly invariant in xi, but only approximately
% invariant in x0 when vibration is not a constant shift.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'equal_level_width_vibration_results');
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

%% Load no-vibration waveforms
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);
xRef = xCell{1}(:);
levels = (0.10:0.10:0.90)';
nLevel = numel(levels);

% Representative gap for detailed plots.
[~, repIdx] = min(abs(gapList - 0.8));

%% Define synthetic vibration fields u(x0)
span = max(xRef) - min(xRef);
xMid = 0.5 * (min(xRef) + max(xRef));

cases = struct([]);
cases(1).name = "No vibration";
cases(1).tag = "no_vibration";
cases(1).u_fun = @(xNorm) zeros(size(xNorm));

cases(2).name = "Constant shift";
cases(2).tag = "constant_shift";
cases(2).u_fun = @(xNorm) 0.18 * ones(size(xNorm));

cases(3).name = "Slow linear drift";
cases(3).tag = "slow_linear";
cases(3).u_fun = @(xNorm) 0.18 + 0.06 * xNorm;

cases(4).name = "Moderate sine drift";
cases(4).tag = "moderate_sine";
cases(4).u_fun = @(xNorm) 0.18 + 0.08 * sin(pi * xNorm);

cases(5).name = "Fast sine drift";
cases(5).tag = "fast_sine";
cases(5).u_fun = @(xNorm) 0.18 + 0.10 * sin(2.5 * pi * xNorm);

nCase = numel(cases);

% Build a common rigid-coordinate grid whose mapped relative coordinate
% always stays inside the original no-vibration waveform support.
uMinGlobal = inf;
uMaxGlobal = -inf;
for ic = 1:nCase
    xNormRef = (xRef - xMid) / max(span, eps);
    uRef = cases(ic).u_fun(xNormRef);
    uMinGlobal = min(uMinGlobal, min(uRef));
    uMaxGlobal = max(uMaxGlobal, max(uRef));
end
x0Grid = linspace(min(xRef) + uMaxGlobal, max(xRef) + uMinGlobal, numel(xRef))';
xNorm = (x0Grid - xMid) / max(span, eps);

%% Evaluate widths for each gap and each vibration case
rows(nGap * nCase, 1) = struct( ...
    'gap_mm', NaN, ...
    'case_name', "", ...
    'mean_abs_err_x0_mm', NaN, ...
    'max_abs_err_x0_mm', NaN, ...
    'mean_abs_err_xi_mm', NaN, ...
    'max_abs_err_xi_mm', NaN, ...
    'mean_rel_err_x0_percent', NaN);
rowId = 0;

repData = struct();

for ig = 1:nGap
    xTrue = xCell{ig}(:);
        yTrue = yCell{ig}(:);
    [baseTrue, ampTrue, ~] = normalize_waveform(yTrue);
    yTargets = baseTrue + ampTrue * levels;
    [wTrue, xLeftTrue, xRightTrue] = level_widths_abs(xTrue, yTrue, yTargets);

    for ic = 1:nCase
        u = cases(ic).u_fun(xNorm);
        u = u(:);
        xiGrid = x0Grid - u;
        if any(diff(xiGrid) <= 0)
            error('Synthetic vibration case "%s" breaks monotonicity of xi(x0).', cases(ic).name);
        end

        % Observed waveform sampled on rigid coordinate x0.
        yObs = interp1(xTrue, yTrue, xiGrid, 'pchip', 'extrap');

        % Apparent widths directly measured in rigid coordinate x0.
        [wX0, xLeftX0, xRightX0] = level_widths_abs(x0Grid, yObs, yTargets);

        % Widths measured after mapping back to true relative coordinate xi.
        [~, sortIdx] = sort(xiGrid, 'ascend');
        yObsSorted = yObs(sortIdx);
        xiSorted = xiGrid(sortIdx);
        [wXi, xLeftXi, xRightXi] = level_widths_abs(xiSorted, yObsSorted, yTargets);

        errX0 = wX0 - wTrue;
        errXi = wXi - wTrue;

        rowId = rowId + 1;
        rows(rowId).gap_mm = gapList(ig);
        rows(rowId).case_name = cases(ic).name;
        rows(rowId).mean_abs_err_x0_mm = mean(abs(errX0), 'omitnan');
        rows(rowId).max_abs_err_x0_mm = max(abs(errX0), [], 'omitnan');
        rows(rowId).mean_abs_err_xi_mm = mean(abs(errXi), 'omitnan');
        rows(rowId).max_abs_err_xi_mm = max(abs(errXi), [], 'omitnan');
        rows(rowId).mean_rel_err_x0_percent = 100 * mean(abs(errX0 ./ max(wTrue, eps)), 'omitnan');

        if ig == repIdx
            repData(ic).name = cases(ic).name;
            repData(ic).x0 = x0Grid;
            repData(ic).xi = xiGrid;
            repData(ic).yObs = yObs;
            repData(ic).wTrue = wTrue;
            repData(ic).wX0 = wX0;
            repData(ic).wXi = wXi;
            repData(ic).xLeftTrue = xLeftTrue;
            repData(ic).xRightTrue = xRightTrue;
            repData(ic).xLeftX0 = xLeftX0;
            repData(ic).xRightX0 = xRightX0;
            repData(ic).xLeftXi = xLeftXi;
            repData(ic).xRightXi = xRightXi;
            repData(ic).u = u;
        end
    end
end

resultTable = struct2table(rows);
writetable(resultTable, fullfile(outDir, 'equal_level_width_invariance_summary.csv'));

summaryTable = groupsummary(resultTable, "case_name", "mean", ...
    ["mean_abs_err_x0_mm", "max_abs_err_x0_mm", "mean_abs_err_xi_mm", ...
     "max_abs_err_xi_mm", "mean_rel_err_x0_percent"]);
writetable(summaryTable, fullfile(outDir, 'equal_level_width_invariance_grouped.csv'));

%% Figure 1: representative waveform in x0 and xi
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(xTrue, yTrue, 'k-', 'DisplayName', 'No vibration');
for ic = 2:nCase
    plot(repData(ic).x0, repData(ic).yObs, '-', 'DisplayName', char(repData(ic).name));
end
xlabel('Rigid coordinate x_0 (mm)');
ylabel('Capacitance (pF)');
title('(a) Observed waveforms in x_0');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
plot(xTrue, yTrue, 'k-', 'DisplayName', 'Reference');
for ic = 2:nCase
    [xiSorted, sortIdx] = sort(repData(ic).xi, 'ascend');
    plot(xiSorted, repData(ic).yObs(sortIdx), '--', 'DisplayName', char(repData(ic).name));
end
xlabel('True relative coordinate \xi (mm)');
ylabel('Capacitance (pF)');
title('(b) Waveforms remapped to \xi');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
for ic = 1:nCase
    plot(levels, repData(ic).wX0, '-o', 'MarkerSize', 4, 'DisplayName', char(repData(ic).name));
end
plot(levels, repData(1).wTrue, 'k--', 'DisplayName', 'Reference width');
xlabel('Normalized level q');
ylabel('Apparent width in x_0 (mm)');
title('(c) Equal-level widths in x_0');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
for ic = 1:nCase
    plot(levels, repData(ic).wXi, '-o', 'MarkerSize', 4, 'DisplayName', char(repData(ic).name));
end
plot(levels, repData(1).wTrue, 'k--', 'DisplayName', 'Reference width');
xlabel('Normalized level q');
ylabel('True width in \xi (mm)');
title('(d) Equal-level widths in \xi');
legend('Location', 'best', 'FontSize', 7.2);

title(tl, sprintf('Width invariance check at representative gap g = %.1f mm', gapList(repIdx)), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_width_invariance_representative.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_width_invariance_representative.pdf'), 'ContentType', 'vector');

%% Figure 2: grouped error summary
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
caseCats = categorical(summaryTable.case_name);

nexttile;
bar(caseCats, summaryTable.mean_mean_abs_err_x0_mm, 0.65, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean |width error| in x_0 (mm)');
title('(a) Apparent width error');

nexttile;
bar(caseCats, summaryTable.mean_mean_abs_err_xi_mm, 0.65, ...
    'FaceColor', [0.84 0.89 0.75], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean |width error| in \xi (mm)');
title('(b) True width error');

title(tl2, 'Grouped equal-level width errors across all gaps', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_width_invariance_grouped.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_width_invariance_grouped.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Equal-level width invariance verification ===\n');
fprintf('Representative gap for detailed plots: %.1f mm\n\n', gapList(repIdx));
disp(summaryTable(:, {'case_name', 'mean_mean_abs_err_x0_mm', 'mean_mean_abs_err_xi_mm', ...
    'mean_mean_rel_err_x0_percent'}));
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
