%% Plot 03: single-case waveform fit and residual comparison.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
matPath = fullfile(ctx.outDir, 'inv_log_2_key_cases.mat');
if ~exist(matPath, 'file')
    run(fullfile(ctx.thisDir, 'Step01_Run_InvLog2_Key_Cases.m'));
end

S = load(matPath, 'caseResults');
caseResults = S.caseResults;

if ~exist('caseName', 'var') || isempty(caseName)
    caseName = "default_10dB";
end

targetModels = ["baseline_1_over_g", "field_basis_inv_log_2"];
pick = cell(numel(targetModels), 1);
for im = 1:numel(targetModels)
    mask = arrayfun(@(s) s.case_name == caseName && s.model_label == targetModels(im), caseResults);
    idx = find(mask, 1);
    if isempty(idx)
        error('Case %s with model %s was not found in cached results.', caseName, targetModels(im));
    end
    pick{im} = caseResults(idx);
end

fig = figure('Name', 'inv_log_2 single-case waveform fit', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
highMap = pick{1}.highMap;
plot(ax1, highMap.x_v(:), highMap.V_a(:), 'k.', 'MarkerSize', 6, ...
    'DisplayName', 'measured waveform');
hold(ax1, 'on');
colors = [0.35 0.35 0.35; 0.00 0.45 0.74];
for im = 1:numel(pick)
    xFit = highMap.x_v(:) - pick{im}.result.dx_used - fit_u(pick{im}.result.fit.p, highMap.t_v(:));
    [xSorted, sortIdx] = sort(highMap.x_v(:), 'ascend');
    plot(ax1, xSorted, pick{im}.result.VFit(sortIdx), '-', ...
        'Color', colors(im, :), 'LineWidth', 1.3, ...
        'DisplayName', char(pick{im}.model_label));
    pick{im}.residual = highMap.V_a(:) - pick{im}.result.VFit(:);
    pick{im}.xFit = xFit;
end
title(ax1, '(a) Waveform fit', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
for im = 1:numel(pick)
    plot(ax2, highMap.t_v(:), pick{im}.residual(:), '.', ...
        'Color', colors(im, :), 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s residual', pick{im}.model_label));
    hold(ax2, 'on');
end
yline(ax2, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
title(ax2, '(b) Residual comparison', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Time t (s)', 'FontSize', 9);
ylabel(ax2, 'Residual (V)', 'FontSize', 9);
legend(ax2, 'Location', 'best');

summaryText = sprintf([ ...
    '%s\n', ...
    'baseline: gap err %.4g mm, freq err %.4g Hz, rmse %.4g\n', ...
    'inv\\_log\\_2: gap err %.4g mm, freq err %.4g Hz, rmse %.4g'], ...
    caseName, ...
    abs(pick{1}.result.g_used - pick{1}.truthCase.g_high), ...
    mean(abs(sort(pick{1}.result.f_id(:)) - sort(pick{1}.truthCase.f(:)))), ...
    pick{1}.result.rmse, ...
    abs(pick{2}.result.g_used - pick{2}.truthCase.g_high), ...
    mean(abs(sort(pick{2}.result.f_id(:)) - sort(pick{2}.truthCase.f(:)))), ...
    pick{2}.result.rmse);
annotation(fig, 'textbox', [0.16, 0.78, 0.52, 0.12], ...
    'String', summaryText, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'Interpreter', 'tex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8);

apply_inv_log_2_figure_style(ax1, 170, 110);
