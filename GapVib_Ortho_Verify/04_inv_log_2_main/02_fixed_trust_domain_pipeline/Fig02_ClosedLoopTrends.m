%% Plot 02: fixed-trust closed-loop trends.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
csvPath = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_closed_loop.csv');
if ~exist(csvPath, 'file')
    runMode = "smoke";
    run(fullfile(ctx.thisDir, '02_fixed_trust_domain_pipeline', ...
        'Step02_ClosedLoop_FixedTrust.m'));
end
T = readtable(csvPath);

modelOrder = ["baseline_1_over_g", "field_basis_inv_log_2", "power_law_p025_order2"];
modelLabels = ["1/g", "inv log", "power"];
colors = [0.55 0.55 0.55; 0.00 0.45 0.74; 0.85 0.33 0.10];

fig = figure('Name', 'fixed trust closed-loop trends', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
firstAx = [];

metricNames = ["gap_error_mm", "mean_freq_error_Hz", "mean_amp_error_mm", "rmse"];
yLabels = ["Gap error (mm)", "Frequency error (Hz)", "Amplitude error (mm)", "Voltage RMSE (V)"];
titles = ["(a) Gap", "(b) Frequency", "(c) Amplitude", "(d) Voltage fit"];

for ia = 1:4
    ax = nexttile;
    if isempty(firstAx)
        firstAx = ax;
    end
    hold(ax, 'on');
    for im = 1:numel(modelOrder)
        idx = string(T.model_label) == modelOrder(im);
        if ~any(idx)
            continue;
        end
        [gGroup, gVals] = findgroups(T.g_true_mm(idx));
        yMean = splitapply(@mean, T.(metricNames(ia))(idx), gGroup);
        plot(ax, gVals, yMean, 'o-', 'Color', colors(im, :), ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', modelLabels(im));
    end
    xlabel(ax, 'True gap g (mm)', 'FontSize', 9);
    ylabel(ax, yLabels(ia), 'FontSize', 9);
    title(ax, titles(ia), 'FontWeight', 'normal', 'FontSize', 9);
    if ia == 1
        legend(ax, 'Location', 'best');
    end
end
apply_inv_log_2_figure_style(firstAx, 180, 130);
