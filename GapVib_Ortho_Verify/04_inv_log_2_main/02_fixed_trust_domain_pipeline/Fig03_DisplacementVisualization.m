%% Plot 03: fixed-trust displacement post-processing visualization.
%
% This script is a standalone visualization entry for Step03. It runs
% Step03 if the displacement result file is missing, then displays the saved
% displacement reconstruction and summary metrics.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
csvPath = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_displacement.csv');
summaryPath = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_displacement_summary.csv');
if ~exist(csvPath, 'file') || ~exist(summaryPath, 'file')
    run(fullfile(ctx.thisDir, '02_fixed_trust_domain_pipeline', ...
        'Step03_Displacement_FixedTrust.m'));
end
D = readtable(csvPath);
S = readtable(summaryPath);

t = D.t_s(:);
uTruth = D.u_truth_mm(:);
uIdentified = D.u_identified_mm(:);
uError = D.u_error_mm(:);

fig = figure('Name', 'fixed trust displacement visualization', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
preview = t <= min(t) + 0.018;
plot(ax1, t(preview), uTruth(preview), 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'truth');
hold(ax1, 'on');
plot(ax1, t(preview), uIdentified(preview), '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.2, 'DisplayName', 'identified');
ylabel(ax1, 'Displacement (mm)', 'FontSize', 9);
title(ax1, '(a) Displacement preview', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
plot(ax2, uTruth, uIdentified, 'k.', 'MarkerSize', 5);
hold(ax2, 'on');
lims = [min([uTruth; uIdentified]), max([uTruth; uIdentified])];
plot(ax2, lims, lims, 'k--', 'LineWidth', 0.8);
xlabel(ax2, 'Truth displacement (mm)', 'FontSize', 9);
ylabel(ax2, 'Identified displacement (mm)', 'FontSize', 9);
title(ax2, '(b) Sample agreement', 'FontWeight', 'normal', 'FontSize', 9);

ax3 = nexttile;
histogram(ax3, uError, 28, 'FaceColor', [0.00 0.45 0.74], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.85);
xline(ax3, 0, 'k--', 'LineWidth', 0.8);
xlabel(ax3, 'Error (mm)', 'FontSize', 9);
ylabel(ax3, 'Count', 'FontSize', 9);
title(ax3, '(c) Error distribution', 'FontWeight', 'normal', 'FontSize', 9);

ax4 = nexttile;
bar(ax4, categorical({'gap','u RMSE','u max'}), ...
    [S.gap_error_mm(1), S.u_rmse_mm(1), S.u_max_abs_error_mm(1)]);
ylabel(ax4, 'Error (mm)', 'FontSize', 9);
title(ax4, '(d) Key displacement metrics', 'FontWeight', 'normal', 'FontSize', 9);

apply_inv_log_2_figure_style(ax1, 170, 125);
