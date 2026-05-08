%% Plot 07: fixed-trust model performance comparison.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('runMode', 'var') || isempty(runMode)
    runMode = "smoke";
end

summaryPath = fullfile(ctx.outDir, "inv_log_2_fixed_trust_closed_loop_" + runMode + "_summary.csv");
if runMode == "smoke"
    summaryPath = fullfile(ctx.outDir, "inv_log_2_fixed_trust_closed_loop_summary.csv");
end
if ~exist(summaryPath, 'file')
    run(fullfile(ctx.thisDir, '02_fixed_trust_domain_pipeline', ...
        'Step02_ClosedLoop_FixedTrust.m'));
end
summary = readtable(summaryPath);
idx = strcmp(summary.model_label, 'field_basis_inv_log_2') | ...
    strcmp(summary.model_label, 'baseline_1_over_g') | ...
    strcmp(summary.model_label, 'power_law_p025_order2');
summary = summary(idx, :);
order = ["baseline_1_over_g", "field_basis_inv_log_2", "power_law_p025_order2"];
[~, orderIdx] = ismember(order, string(summary.model_label));
orderIdx = orderIdx(orderIdx > 0);
summary = summary(orderIdx, :);

fig = figure('Name', 'fixed trust performance comparison', 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
modelCats = categorical(short_model_labels(summary.model_label));
modelCats = reordercats(modelCats, short_model_labels(summary.model_label));

ax1 = nexttile;
bar(ax1, modelCats, summary.mean_gap_error_mm);
ylabel(ax1, 'Gap error (mm)', 'FontSize', 9);
title(ax1, '(a) Gap', 'FontWeight', 'normal', 'FontSize', 9);
xtickangle(ax1, 18);

ax2 = nexttile;
bar(ax2, modelCats, summary.mean_freq_error_Hz);
ylabel(ax2, 'Frequency error (Hz)', 'FontSize', 9);
title(ax2, '(b) Frequency', 'FontWeight', 'normal', 'FontSize', 9);
xtickangle(ax2, 18);

ax3 = nexttile;
bar(ax3, modelCats, summary.mean_amp_error_mm);
ylabel(ax3, 'Amplitude error (mm)', 'FontSize', 9);
title(ax3, '(c) Amplitude', 'FontWeight', 'normal', 'FontSize', 9);
xtickangle(ax3, 18);

apply_inv_log_2_figure_style(ax1, 180, 72);

function labels = short_model_labels(modelLabel)
labels = strings(size(modelLabel));
for ii = 1:numel(modelLabel)
    switch string(modelLabel(ii))
        case "baseline_1_over_g"
            labels(ii) = "1/g";
        case "field_basis_inv_log_2"
            labels(ii) = "inv log";
        case "power_law_p025_order2"
            labels(ii) = "power";
        otherwise
            labels(ii) = string(modelLabel(ii));
    end
end
end
