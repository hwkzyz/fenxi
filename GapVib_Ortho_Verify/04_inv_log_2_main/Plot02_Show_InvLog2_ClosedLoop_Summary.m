%% Plot 02: closed-loop quick summary for inv_log_2.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
matPath = fullfile(ctx.outDir, 'inv_log_2_closed_loop_quick.mat');
if ~exist(matPath, 'file')
    run(fullfile(ctx.thisDir, 'Step02_Run_InvLog2_ClosedLoop_Quick.m'));
end

S = load(matPath, 'summary');
summary = S.summary;

metricNames = { ...
    'mean_gap_error_mm', ...
    'mean_freq_error_Hz', ...
    'mean_amp_error_mm', ...
    'mean_wave_rmse'};
metricTitles = { ...
    '(a) Mean gap error', ...
    '(b) Mean frequency error', ...
    '(c) Mean amplitude error', ...
    '(d) Mean waveform RMSE'};
metricUnits = {'mm', 'Hz', 'mm', ''};

labels = cellstr(summary.model_label);
colors = [0.35 0.35 0.35; 0.00 0.45 0.74; 0.85 0.33 0.10];

fig = figure('Name', 'inv_log_2 closed-loop summary', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for imetric = 1:numel(metricNames)
    ax = nexttile;
    y = summary.(metricNames{imetric});
    b = bar(ax, y, 0.65, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = colors(1:numel(y), :);
    xticks(ax, 1:numel(labels));
    xticklabels(ax, labels);
    xtickangle(ax, 18);
    ylabel(ax, metricUnits{imetric}, 'FontSize', 9);
    title(ax, metricTitles{imetric}, 'FontWeight', 'normal', 'FontSize', 9);

    yMax = max(y);
    if yMax <= 0
        yMax = 1;
    end
    ylim(ax, [0, 1.18 * yMax]);
    for ii = 1:numel(y)
        text(ax, ii, y(ii) + 0.03 * yMax, sprintf('%.3g', y(ii)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 8);
    end
end

apply_inv_log_2_figure_style(gca, 170, 120);
