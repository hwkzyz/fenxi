%% Visualize waveform differences in time, rigid-angle, and relative coordinates
% Goal:
%   1) show why a vibrating waveform differs from the no-vibration waveform
%      in time domain t-V and rigid coordinate x0-V;
%   2) show why they collapse in the true relative coordinate xi-V;
%   3) compare several vibration parameter sets to build intuition.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'time_angle_coordinate_relationship_results');
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

%% Load one baseline no-vibration waveform
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
[~, refIdx] = min(abs(gapList - 0.8));
xRef = xCell{refIdx}(:);
yRef = yCell{refIdx}(:);

V_tip = 3.0e5;  % mm/s
x0Center = mean(xRef);
tRef = (xRef - x0Center) / V_tip;
yNoVib = yRef;

%% Define several vibration cases for intuition
cases = struct( ...
    'name', {'No vibration', 'Constant shift', 'Low-frequency sine', 'User multimode', 'Fast varying vibration'}, ...
    'short', {'No vib', 'Const shift', 'Low-f sine', 'User multimode', 'Fast varying'}, ...
    'u0', {0.0, 0.6, 0.0, 0.0, 0.0}, ...
    'A', {0.0, 0.0, 0.8, [0.5, 0.4], [0.7, 0.45]}, ...
    'f', {0.0, 0.0, 250, [500, 1300], [900, 2200]}, ...
    'phi', {0.0, 0.0, pi/6, [pi/4, -pi/3], [pi/5, -pi/2]} ...
    );

nCase = numel(cases);
caseColors = lines(nCase);

results = repmat(struct( ...
    'name', '', 'short', '', 'x0', [], 't', [], 'u', [], 'xi', [], ...
    'yObs', [], 'yXi', [], 'xiSorted', [], 'peakShiftX0', 0, ...
    'peakShiftT_us', 0, 'rmseX0', 0, 'rmseXi', 0, ...
    'xiSlopeMin', 0, 'xiSlopeMax', 0), nCase, 1);

%% Generate waveforms under each vibration case
for i = 1:nCase
    results(i) = build_case_result(cases(i), xRef, yRef, V_tip, x0Center);
end

%% Summary table
summaryTable = table( ...
    string({results.name})', ...
    [results.peakShiftX0]', ...
    [results.peakShiftT_us]', ...
    [results.rmseX0]', ...
    [results.rmseXi]', ...
    [results.xiSlopeMin]', ...
    [results.xiSlopeMax]', ...
    'VariableNames', {'case_name', 'peak_shift_x0_mm', 'peak_shift_t_us', ...
    'rmse_in_x0', 'rmse_in_xi', 'min_dxi_dx0', 'max_dxi_dx0'});
writetable(summaryTable, fullfile(outDir, 'coordinate_relationship_summary.csv'));

%% Figure 1: displacement functions and coordinate mapping
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 21]);
tl1 = tiledlayout(fig1, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 2:nCase
    nexttile(2*i-3); hold on;
    plot(results(i).t * 1e6, results(i).u, '-', 'Color', caseColors(i, :));
    yline(0, 'k--');
    xlabel('Time t (\mus)');
    ylabel('Displacement u(t) (mm)');
    title(sprintf('(%c) %s: u(t)', char('a' + (i-2)*2), results(i).short));

    nexttile(2*i-2); hold on;
    plot(results(i).x0, results(i).xi, '-', 'Color', caseColors(i, :), 'DisplayName', '\xi(x_0)');
    plot(results(i).x0, results(i).x0, 'k--', 'DisplayName', '\xi = x_0');
    xlabel('Rigid coordinate x_0 (mm)');
    ylabel('Relative coordinate \xi (mm)');
    title(sprintf('(%c) %s: \\xi-x_0 mapping', char('b' + (i-2)*2), results(i).short));
    legend('Location', 'northwest', 'FontSize', 7.0);
end
exportgraphics(fig1, fullfile(outDir, 'fig1_displacement_and_mapping.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_displacement_and_mapping.pdf'), 'ContentType', 'vector');

%% Figure 2: waveform comparison in time domain t-V
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 14]);
tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plotCases = 2:nCase;
for k = 1:numel(plotCases)
    i = plotCases(k);
    nexttile; hold on;
    plot(tRef * 1e6, yNoVib, 'k-', 'DisplayName', 'No vibration');
    plot(results(i).t * 1e6, results(i).yObs, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
    xlabel('Time t (\mus)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) t-V: %s', char('a' + k - 1), results(i).short));
    legend('Location', 'best', 'FontSize', 7.0);
end
exportgraphics(fig2, fullfile(outDir, 'fig2_time_domain_waveforms.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_time_domain_waveforms.pdf'), 'ContentType', 'vector');

%% Figure 3: waveform comparison in rigid coordinate x0-V
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 14]);
tl3 = tiledlayout(fig3, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(plotCases)
    i = plotCases(k);
    nexttile; hold on;
    plot(xRef, yNoVib, 'k-', 'DisplayName', 'No vibration');
    plot(results(i).x0, results(i).yObs, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
    xlabel('Rigid coordinate x_0 (mm)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) x_0-V: %s', char('a' + k - 1), results(i).short));
    legend('Location', 'best', 'FontSize', 7.0);
end
exportgraphics(fig3, fullfile(outDir, 'fig3_rigid_coordinate_waveforms.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_rigid_coordinate_waveforms.pdf'), 'ContentType', 'vector');

%% Figure 4: waveform comparison in true relative coordinate xi-V
fig4 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 14]);
tl4 = tiledlayout(fig4, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(plotCases)
    i = plotCases(k);
    nexttile; hold on;
    plot(xRef, yNoVib, 'k-', 'DisplayName', 'No vibration');
    plot(results(i).xiSorted, results(i).yXi, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
    xlabel('Relative coordinate \xi (mm)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) \\xi-V: %s', char('a' + k - 1), results(i).short));
    legend('Location', 'best', 'FontSize', 7.0);
end
exportgraphics(fig4, fullfile(outDir, 'fig4_relative_coordinate_waveforms.png'), 'Resolution', 300);
exportgraphics(fig4, fullfile(outDir, 'fig4_relative_coordinate_waveforms.pdf'), 'ContentType', 'vector');

%% Figure 5: residual comparison in x0 and xi
fig5 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 14]);
tl5 = tiledlayout(fig5, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(plotCases)
    i = plotCases(k);
    nexttile; hold on;
    yRefOnX0 = interp1(xRef, yRef, results(i).x0, 'pchip', NaN);
    resX0 = results(i).yObs - yRefOnX0;
    plot(results(i).x0, resX0, '-', 'Color', caseColors(i, :), 'DisplayName', 'Residual in x_0');
    yline(0, 'k--');
    yyaxis right;
    plot(results(i).xiSorted, results(i).yXi - interp1(xRef, yRef, results(i).xiSorted, 'pchip', NaN), ...
        '--', 'Color', [0.4 0.4 0.4], 'DisplayName', 'Residual in \xi');
    ylabel('Residual in \xi (pF)');
    yyaxis left;
    xlabel('Coordinate (mm)');
    ylabel('Residual in x_0 (pF)');
    title(sprintf('(%c) Residuals: %s', char('a' + k - 1), results(i).short));
end
exportgraphics(fig5, fullfile(outDir, 'fig5_residual_comparison.png'), 'Resolution', 300);
exportgraphics(fig5, fullfile(outDir, 'fig5_residual_comparison.pdf'), 'ContentType', 'vector');

%% Figure 6: all cases overlaid in one view for each coordinate
fig6 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 18]);
tl6 = tiledlayout(fig6, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(tRef * 1e6, yNoVib, 'k-', 'DisplayName', 'No vibration');
for i = 2:nCase
    plot(results(i).t * 1e6, results(i).yObs, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
end
xlabel('Time t (\mus)');
ylabel('Capacitance V (pF)');
title('(a) Time-domain overlay');
legend('Location', 'eastoutside', 'FontSize', 7.0);

nexttile; hold on;
plot(xRef, yNoVib, 'k-', 'DisplayName', 'No vibration');
for i = 2:nCase
    plot(results(i).x0, results(i).yObs, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
end
xlabel('Rigid coordinate x_0 (mm)');
ylabel('Capacitance V (pF)');
title('(b) Rigid-coordinate overlay');

nexttile; hold on;
plot(xRef, yNoVib, 'k-', 'DisplayName', 'No vibration');
for i = 2:nCase
    plot(results(i).xiSorted, results(i).yXi, '-', 'Color', caseColors(i, :), ...
        'DisplayName', results(i).short);
end
xlabel('Relative coordinate \xi (mm)');
ylabel('Capacitance V (pF)');
title('(c) Relative-coordinate overlay');
exportgraphics(fig6, fullfile(outDir, 'fig6_overlay_all_coordinates.png'), 'Resolution', 300);
exportgraphics(fig6, fullfile(outDir, 'fig6_overlay_all_coordinates.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Time / angle / relative coordinate relationship visualization ===\n');
fprintf('Reference gap: %.1f mm\n', gapList(refIdx));
fprintf('Tip speed V = %.3e mm/s\n\n', V_tip);
disp(summaryTable);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function result = build_case_result(caseDef, xRef, yRef, V_tip, x0Center)
    result.name = caseDef.name;
    result.short = caseDef.short;

    u0 = caseDef.u0;
    A = caseDef.A;
    f = caseDef.f;
    phi = caseDef.phi;

    if isscalar(A)
        AAbsMax = abs(A);
    else
        AAbsMax = sum(abs(A));
    end
    uMax = abs(u0) + AAbsMax;

    x0 = linspace(min(xRef) + uMax, max(xRef) - uMax, numel(xRef))';
    t = (x0 - x0Center) / V_tip;
    u = u0 + multimode_vibration(t, A, f, phi);
    xi = x0 - u;

    if any(diff(xi) <= 0)
        error('Case "%s" makes xi(x0) non-monotone. Increase V_tip or reduce vibration.', caseDef.name);
    end

    yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    [xiSorted, idx] = sort(xi, 'ascend');
    yXi = yObs(idx);

    yRefOnX0 = interp1(xRef, yRef, x0, 'pchip', NaN);
    yRefOnXi = interp1(xRef, yRef, xiSorted, 'pchip', NaN);
    rmseX0 = sqrt(mean((yObs - yRefOnX0).^2, 'omitnan'));
    rmseXi = sqrt(mean((yXi - yRefOnXi).^2, 'omitnan'));

    [~, iPeakObs] = max(yObs);
    [~, iPeakRef] = max(yRefOnX0);
    peakShiftX0 = x0(iPeakObs) - x0(iPeakRef);
    peakShiftT_us = (peakShiftX0 / V_tip) * 1e6;

    dxi_dx0 = gradient(xi, x0);

    result.x0 = x0;
    result.t = t;
    result.u = u;
    result.xi = xi;
    result.yObs = yObs;
    result.yXi = yXi;
    result.xiSorted = xiSorted;
    result.peakShiftX0 = peakShiftX0;
    result.peakShiftT_us = peakShiftT_us;
    result.rmseX0 = rmseX0;
    result.rmseXi = rmseXi;
    result.xiSlopeMin = min(dxi_dx0);
    result.xiSlopeMax = max(dxi_dx0);
end

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
    if isscalar(A)
        if A == 0
            u = zeros(size(t));
        else
            u = A * sin(2 * pi * f * t + phi);
        end
        return;
    end

    u = zeros(size(t));
    for k = 1:numel(A)
        u = u + A(k) * sin(2 * pi * f(k) * t + phi(k));
    end
end
