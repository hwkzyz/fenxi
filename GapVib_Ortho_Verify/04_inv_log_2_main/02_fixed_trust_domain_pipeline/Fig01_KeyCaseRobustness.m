%% Plot 01: fixed-trust key-case robustness.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
csvPath = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_key_cases.csv');
if ~exist(csvPath, 'file')
    run(fullfile(ctx.thisDir, '02_fixed_trust_domain_pipeline', ...
        'Step01_KeyCases_FixedTrust.m'));
end
T = readtable(csvPath);

methodKeep = strcmp(T.model_label, 'baseline_1_over_g') | ...
    strcmp(T.model_label, 'field_basis_inv_log_2');
T = T(methodKeep, :);
[caseNames, caseOrder] = unique(string(T.case_name), 'stable');
methodNames = ["baseline_1_over_g", "field_basis_inv_log_2"];
methodLabels = ["1/g", "inv log"];

gapErr = nan(numel(caseNames), numel(methodNames));
freqErr = nan(numel(caseNames), numel(methodNames));
ampErr = nan(numel(caseNames), numel(methodNames));
for ic = 1:numel(caseNames)
    for im = 1:numel(methodNames)
        idx = string(T.case_name) == caseNames(ic) & string(T.model_label) == methodNames(im);
        gapErr(ic, im) = mean(T.gap_error_mm(idx), 'omitnan');
        freqErr(ic, im) = mean(T.mean_freq_error_Hz(idx), 'omitnan');
        ampErr(ic, im) = mean(T.mean_amp_error_mm(idx), 'omitnan');
    end
end

caseLabels = shorten_case_labels(caseNames);
caseCats = categorical(caseLabels);
caseCats = reordercats(caseCats, caseLabels);
fig = figure('Name', 'fixed trust key-case robustness', 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
bar(ax1, caseCats, gapErr);
ylabel(ax1, 'Gap error (mm)', 'FontSize', 9);
title(ax1, '(a) Gap robustness', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, methodLabels, 'Location', 'northwest');

ax2 = nexttile;
bar(ax2, caseCats, freqErr);
ylabel(ax2, 'Frequency error (Hz)', 'FontSize', 9);
title(ax2, '(b) Frequency robustness', 'FontWeight', 'normal', 'FontSize', 9);

ax3 = nexttile;
bar(ax3, caseCats, ampErr);
ylabel(ax3, 'Amplitude error (mm)', 'FontSize', 9);
xlabel(ax3, 'Key case', 'FontSize', 9);
title(ax3, '(c) Amplitude robustness', 'FontWeight', 'normal', 'FontSize', 9);

apply_inv_log_2_figure_style(ax1, 180, 140);

function labels = shorten_case_labels(caseNames)
labels = strings(size(caseNames));
for ii = 1:numel(caseNames)
    switch caseNames(ii)
        case "default_10dB"
            labels(ii) = "default";
        case "phase_sensitive_20dB"
            labels(ii) = "phase";
        case "weak_mode2_5dB"
            labels(ii) = "weak mode2";
        case "gap_bias_minus005_10dB"
            labels(ii) = "gap bias";
        otherwise
            labels(ii) = caseNames(ii);
    end
end
end
