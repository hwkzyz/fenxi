%% Visualize Step07J method comparison for 20250527
% Reads Step07J trend CSV files and plots method-level EO, frequency,
% amplitude, and RMSE comparisons. This script does not rerun identification.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_method_comparison');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cases = {
    'main18_legacy_basexc', fullfile(outDir, 'Step07J_NestedStaticWarp_VPFullWave_Trend_20250527_B1_S136.csv')
    'main80_basexc_gradienttpl', fullfile(outDir, 'Step07J_NestedStaticWarp_VPFullWave_Trend_20250527_B1_S136_Diag80L_W3S1_BaseXc_GradientTpl.csv')
    };

style = paper_style_local();

for ic = 1:size(cases, 1)
    caseName = cases{ic, 1};
    trendFile = cases{ic, 2};
    if ~isfile(trendFile)
        warning('Trend file not found, skipping %s: %s', caseName, trendFile);
        continue;
    end

    T = readtable(trendFile);
    Summary = build_method_summary_local(T);
    summaryFile = fullfile(outDir, sprintf('Step07J_MethodComparison_%s_20250527_B1_S136.csv', caseName));
    writetable(Summary, summaryFile);

    plot_method_bars_local(Summary, caseName, figDir, style);
    plot_window_trends_local(T, caseName, figDir, style);
    fprintf('Visualized %s (%d windows). Summary: %s\n', caseName, height(T), summaryFile);
end

fprintf('Saved figures under: %s\n', figDir);

%% Local functions
function Summary = build_method_summary_local(T)
methods = {'direct','fixed','gap','tilt','shift'};
labels = {'Direct low template','Fixed','Gap only (main)','Gap + tilt ablation','Gap + tilt + shift ablation'};
rows = cell(numel(methods), 1);
for i = 1:numel(methods)
    m = methods{i};
    eo = T.([m '_EO']);
    f = T.([m '_frequency_hz']);
    A = T.([m '_amplitude_mm']);
    rmse = T.([m '_rmse_mV']);
    is14 = eo == 14;
    wrong = eo ~= 14;
    wrongText = format_wrong_eo_text_local(eo(wrong));

    rows{i} = table(string(labels{i}), height(T), nnz(is14), nnz(wrong), wrongText, ...
        mean(f, 'omitnan'), std(f, 'omitnan'), median(f, 'omitnan'), ...
        mean(A, 'omitnan'), median(A, 'omitnan'), ...
        mean(rmse, 'omitnan'), median(rmse, 'omitnan'), ...
        mean(f(is14), 'omitnan'), std(f(is14), 'omitnan'), ...
        mean(A(is14), 'omitnan'), median(A(is14), 'omitnan'), ...
        mean(rmse(is14), 'omitnan'), median(rmse(is14), 'omitnan'), ...
        'VariableNames', {'method','n_windows','eo14_count','wrong_count','wrong_eo_breakdown', ...
        'mean_freq_all_hz','std_freq_all_hz','median_freq_all_hz', ...
        'mean_amp_all_mm','median_amp_all_mm','mean_rmse_all_mV','median_rmse_all_mV', ...
        'mean_freq_eo14_hz','std_freq_eo14_hz','mean_amp_eo14_mm','median_amp_eo14_mm', ...
        'mean_rmse_eo14_mV','median_rmse_eo14_mV'});
end
Summary = vertcat(rows{:});
Summary.eo14_rate = Summary.eo14_count ./ Summary.n_windows;
end

function textOut = format_wrong_eo_text_local(eoWrong)
if isempty(eoWrong)
    textOut = "none";
    return;
end
u = unique(eoWrong(:)).';
parts = strings(1, numel(u));
for i = 1:numel(u)
    parts(i) = sprintf('EO%d=%d', u(i), nnz(eoWrong == u(i)));
end
textOut = strjoin(parts, '; ');
end

function plot_method_bars_local(Summary, caseName, figDir, style)
methodLabels = cellstr(Summary.method);
x = 1:height(Summary);
colors = method_colors_local(height(Summary));

fig = figure('Name', ['Step07J method bars ' caseName], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 13]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar_with_colors_local(x, Summary.eo14_rate, colors);
ylim([0, 1.05]);
ylabel('EO14 rate');
title('Frequency lock');
set_method_xticks_local(gca, x, methodLabels, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_freq_eo14_hz, colors);
ylabel('Frequency (Hz)');
title('EO14-only frequency');
set_method_xticks_local(gca, x, methodLabels, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_amp_eo14_mm, colors);
ylabel('Amplitude (mm)');
title('EO14-only amplitude');
set_method_xticks_local(gca, x, methodLabels, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_rmse_eo14_mV, colors);
ylabel('RMSE (mV)');
title('EO14-only waveform RMSE');
set_method_xticks_local(gca, x, methodLabels, style);
format_axes_local(gca, style);

export_figure_local(fig, fullfile(figDir, sprintf('Step07J_MethodBars_%s_20250527_B1_S136.png', caseName)));
end

function plot_window_trends_local(T, caseName, figDir, style)
methods = {'direct','fixed','gap','tilt','shift'};
labels = {'Direct','Fixed','Gap only (main)','Gap + tilt ablation','Gap + tilt + shift ablation'};
colors = method_colors_local(numel(methods));
markers = {'o','s','^','d','v'};

fig = figure('Name', ['Step07J window trends ' caseName], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 15]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
for i = 1:numel(methods)
    y = T.([methods{i} '_EO']);
    plot(T.window_id, y, '-', 'Color', colors(i, :), 'LineWidth', 1.0, ...
        'Marker', markers{i}, 'MarkerSize', 3.4, 'DisplayName', labels{i});
end
yline(14, '--', 'Color', style.gray, 'LineWidth', 0.9, 'DisplayName', 'EO14');
ylabel('EO');
title('EO trend');
format_axes_local(gca, style);
legend('Location', 'eastoutside');

nexttile;
hold on;
for i = 1:numel(methods)
    y = T.([methods{i} '_amplitude_mm']);
    plot(T.window_id, y, '-', 'Color', colors(i, :), 'LineWidth', 1.0, ...
        'Marker', markers{i}, 'MarkerSize', 3.4, 'DisplayName', labels{i});
end
ylabel('Amplitude (mm)');
title('Amplitude trend');
format_axes_local(gca, style);

nexttile;
hold on;
for i = 1:numel(methods)
    y = T.([methods{i} '_rmse_mV']);
    plot(T.window_id, y, '-', 'Color', colors(i, :), 'LineWidth', 1.0, ...
        'Marker', markers{i}, 'MarkerSize', 3.4, 'DisplayName', labels{i});
end
xlabel('Window index');
ylabel('RMSE (mV)');
title('Waveform RMSE trend');
format_axes_local(gca, style);

export_figure_local(fig, fullfile(figDir, sprintf('Step07J_WindowTrends_%s_20250527_B1_S136.png', caseName)));
end

function bar_with_colors_local(x, y, colors)
b = bar(x, y, 0.72, 'FaceColor', 'flat', 'EdgeColor', [0.15, 0.15, 0.15], 'LineWidth', 0.6);
b.CData = colors;
end

function set_method_xticks_local(ax, x, labels, style)
set(ax, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 28);
ax.FontName = style.fontName;
ax.FontSize = style.fontSize;
end

function colors = method_colors_local(n)
base = [
    0.08, 0.08, 0.08
    0.36, 0.36, 0.36
    0.00, 0.32, 0.72
    0.00, 0.45, 0.28
    0.84, 0.10, 0.10
    ];
colors = base(1:n, :);
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.gray = [0.55, 0.55, 0.55];
end

function format_axes_local(ax, style)
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.TickDir = 'in';
ax.FontName = style.fontName;
ax.FontSize = style.fontSize;
end

function export_figure_local(fig, pngFile)
[folder, name] = fileparts(pngFile);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
exportgraphics(fig, pngFile, 'Resolution', 300);
pdfFile = fullfile(folder, [name, '.pdf']);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
end
