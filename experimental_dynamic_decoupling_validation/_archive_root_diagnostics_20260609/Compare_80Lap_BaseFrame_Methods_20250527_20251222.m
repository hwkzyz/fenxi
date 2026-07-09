%% Compare_80Lap_OPRAnchored_Methods_20250527_20251222
% Compare 80-lap / 3-lap sliding-window Step07J results under the agreed
% OPR-anchored route:
%   x_abs is determined from OPR; the low-speed no-vibration template center
%   defines the per-sensor x_rel frame. GradientXRange030 changes the
%   waveform trust region/weights, not the sensor center. gap_only is the
%   paper-facing model; tilt/shift are ablations.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
routeTag = 'OPRAnchoredDirectDrivenGapOnly';
outDir = fullfile(rootDir, 'diagnostics_80lap_opranchored_direct_driven_method_compare');
figDir = fullfile(outDir, 'figures');
if exist(outDir, 'dir') ~= 7; mkdir(outDir); end
if exist(figDir, 'dir') ~= 7; mkdir(figDir); end

cases = {
    '20250527', 'B1_S136', fullfile(rootDir, '20250527_low_speed_gap_prior_decoupling', 'outputs', ...
        'Step07J_NestedStaticWarp_VPFullWave_Trend_20250527_B1_S136_OPRA80D.csv')
    '20251222', 'B1_S123', fullfile(rootDir, '20251222_low_speed_gap_prior_decoupling', 'outputs', ...
        'Step07J_NestedStaticWarp_VPFullWave_Trend_20251222_B1_S123_OPRA80D.csv')
    };

methods = {'direct','fixed','gap','tilt','shift'};
methodLabels = {'Direct low template','Fixed','Gap only','Gap + tilt','Gap + tilt + shift'};
style = paper_style_local();

allSummary = {};
for ic = 1:size(cases, 1)
    dataset = cases{ic, 1};
    sensorTag = cases{ic, 2};
    trendFile = cases{ic, 3};
    if ~isfile(trendFile)
        error('Missing OPR-anchored 80-lap trend file: %s', trendFile);
    end

    T = readtable(trendFile);
    Summary = build_method_summary_local(T, methods, methodLabels);
    Summary.dataset = repmat(string(dataset), height(Summary), 1);
    Summary.sensor_tag = repmat(string(sensorTag), height(Summary), 1);
    Summary = movevars(Summary, {'dataset','sensor_tag'}, 'Before', 1);
    allSummary{ic} = Summary; %#ok<AGROW>

    writetable(Summary, fullfile(outDir, sprintf('MethodSummary_80Lap_%s_%s_%s.csv', routeTag, dataset, sensorTag)));
    writetable(T, fullfile(outDir, sprintf('WindowTrend_80Lap_%s_%s_%s.csv', routeTag, dataset, sensorTag)));
    plot_method_bars_local(Summary, dataset, sensorTag, figDir, style, routeTag);
    plot_window_trends_local(T, methods, methodLabels, dataset, sensorTag, figDir, style, routeTag);
end

CombinedSummary = vertcat(allSummary{:});
writetable(CombinedSummary, fullfile(outDir, sprintf('Combined_MethodSummary_80Lap_%s_20250527_20251222.csv', routeTag)));
plot_cross_dataset_bars_local(CombinedSummary, methodLabels, figDir, style, routeTag);

fprintf('\nSaved OPR-anchored 80-lap comparison under:\n  %s\n', outDir);
disp(CombinedSummary(:, {'dataset','method','n_windows','eo14_count','wrong_count', ...
    'wrong_eo_breakdown','mean_freq_eo14_hz','mean_amp_eo14_mm','mean_rmse_eo14_mV','mean_rmse_all_mV'}));

%% Local functions
function Summary = build_method_summary_local(T, methods, labels)
rows = cell(numel(methods), 1);
for i = 1:numel(methods)
    m = methods{i};
    eo = T.([m '_EO']);
    f = T.([m '_frequency_hz']);
    A = T.([m '_amplitude_mm']);
    rmse = T.([m '_rmse_mV']);
    is14 = eo == 14;
    wrong = ~is14;

    rows{i} = table(string(labels{i}), height(T), nnz(is14), nnz(wrong), ...
        format_wrong_eo_text_local(eo(wrong)), ...
        mean(f, 'omitnan'), std(f, 'omitnan'), median(f, 'omitnan'), ...
        mean(A, 'omitnan'), median(A, 'omitnan'), ...
        mean(rmse, 'omitnan'), median(rmse, 'omitnan'), ...
        mean(f(is14), 'omitnan'), std(f(is14), 'omitnan'), median(f(is14), 'omitnan'), ...
        mean(A(is14), 'omitnan'), median(A(is14), 'omitnan'), ...
        mean(rmse(is14), 'omitnan'), median(rmse(is14), 'omitnan'), ...
        'VariableNames', {'method','n_windows','eo14_count','wrong_count','wrong_eo_breakdown', ...
        'mean_freq_all_hz','std_freq_all_hz','median_freq_all_hz', ...
        'mean_amp_all_mm','median_amp_all_mm','mean_rmse_all_mV','median_rmse_all_mV', ...
        'mean_freq_eo14_hz','std_freq_eo14_hz','median_freq_eo14_hz', ...
        'mean_amp_eo14_mm','median_amp_eo14_mm','mean_rmse_eo14_mV','median_rmse_eo14_mV'});
end
Summary = vertcat(rows{:});
Summary.eo14_rate = Summary.eo14_count ./ Summary.n_windows;
end

function textOut = format_wrong_eo_text_local(eoWrong)
missingCount = nnz(~isfinite(eoWrong));
eoWrong = eoWrong(isfinite(eoWrong));
if isempty(eoWrong) && missingCount == 0
    textOut = "none";
    return;
end
u = unique(eoWrong(:)).';
parts = strings(1, numel(u) + double(missingCount > 0));
for i = 1:numel(u)
    parts(i) = sprintf('EO%d=%d', u(i), nnz(eoWrong == u(i)));
end
if missingCount > 0
    parts(end) = sprintf('missing=%d', missingCount);
end
textOut = strjoin(parts, '; ');
end

function plot_method_bars_local(Summary, dataset, sensorTag, figDir, style, routeTag)
x = 1:height(Summary);
colors = method_colors_local(height(Summary));

fig = figure('Name', [routeTag ' 80-lap method bars ' dataset], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 14]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar_with_colors_local(x, Summary.eo14_rate, colors);
ylim([0, 1.05]);
ylabel('EO14 rate');
title('Frequency identification');
set_method_ticks_local(gca, x, Summary.method, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_freq_eo14_hz, colors);
ylabel('Frequency (Hz)');
title('EO14-only frequency');
set_method_ticks_local(gca, x, Summary.method, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_amp_eo14_mm, colors);
ylabel('Amplitude (mm)');
title('EO14-only amplitude');
set_method_ticks_local(gca, x, Summary.method, style);
format_axes_local(gca, style);

nexttile;
bar_with_colors_local(x, Summary.mean_rmse_eo14_mV, colors);
ylabel('RMSE (mV)');
title('EO14-only waveform RMSE');
set_method_ticks_local(gca, x, Summary.method, style);
format_axes_local(gca, style);

export_figure_local(fig, fullfile(figDir, sprintf('MethodBars_80Lap_%s_%s_%s.png', routeTag, dataset, sensorTag)));
end

function plot_cross_dataset_bars_local(Summary, methodLabels, figDir, style, routeTag)
datasets = unique(Summary.dataset, 'stable');
methods = string(methodLabels(:));
metricNames = {'eo14_rate','mean_amp_eo14_mm','mean_rmse_eo14_mV'};
titles = {'EO14 rate','EO14-only amplitude','EO14-only RMSE'};
ylabels = {'Rate','Amplitude (mm)','RMSE (mV)'};

fig = figure('Name', [routeTag ' 80-lap cross-dataset bars'], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 13]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for imetric = 1:numel(metricNames)
    nexttile;
    Y = nan(numel(methods), numel(datasets));
    for id = 1:numel(datasets)
        for im = 1:numel(methods)
            row = Summary.dataset == datasets(id) & Summary.method == methods(im);
            if any(row)
                Y(im, id) = Summary.(metricNames{imetric})(find(row, 1));
            end
        end
    end
    bar(Y, 'grouped');
    title(titles{imetric});
    ylabel(ylabels{imetric});
    set(gca, 'XTick', 1:numel(methods), 'XTickLabel', methods, 'XTickLabelRotation', 25);
    legend(cellstr(datasets), 'Location', 'best');
    format_axes_local(gca, style);
end
export_figure_local(fig, fullfile(figDir, sprintf('CrossDataset_MethodBars_80Lap_%s_20250527_20251222.png', routeTag)));
end

function plot_window_trends_local(T, methods, labels, dataset, sensorTag, figDir, style, routeTag)
colors = method_colors_local(numel(methods));
markers = {'o','s','^','d','v'};

fig = figure('Name', [routeTag ' 80-lap trends ' dataset], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 16]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
for i = 1:numel(methods)
    plot(T.window_id, T.([methods{i} '_EO']), '-', 'Color', colors(i, :), ...
        'LineWidth', 0.9, 'Marker', markers{i}, 'MarkerSize', 3.0, 'DisplayName', labels{i});
end
yline(14, '--', 'Color', style.gray, 'LineWidth', 0.9, 'DisplayName', 'EO14');
ylabel('EO');
title('EO trend');
format_axes_local(gca, style);
legend('Location', 'eastoutside');

nexttile;
hold on;
for i = 1:numel(methods)
    plot(T.window_id, T.([methods{i} '_amplitude_mm']), '-', 'Color', colors(i, :), ...
        'LineWidth', 0.9, 'Marker', markers{i}, 'MarkerSize', 3.0, 'DisplayName', labels{i});
end
ylabel('Amplitude (mm)');
title('Amplitude trend');
format_axes_local(gca, style);

nexttile;
hold on;
for i = 1:numel(methods)
    plot(T.window_id, T.([methods{i} '_rmse_mV']), '-', 'Color', colors(i, :), ...
        'LineWidth', 0.9, 'Marker', markers{i}, 'MarkerSize', 3.0, 'DisplayName', labels{i});
end
xlabel('Window index');
ylabel('RMSE (mV)');
title('Waveform RMSE trend');
format_axes_local(gca, style);

export_figure_local(fig, fullfile(figDir, sprintf('WindowTrends_80Lap_%s_%s_%s.png', routeTag, dataset, sensorTag)));
end

function bar_with_colors_local(x, y, colors)
b = bar(x, y, 0.72, 'FaceColor', 'flat', 'EdgeColor', [0.15, 0.15, 0.15], 'LineWidth', 0.6);
b.CData = colors;
end

function set_method_ticks_local(ax, x, labels, style)
set(ax, 'XTick', x, 'XTickLabel', cellstr(labels), 'XTickLabelRotation', 25);
ax.FontName = style.fontName;
ax.FontSize = style.tickFontSize;
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
style.tickFontSize = 8.2;
style.gray = [0.45, 0.45, 0.45];
end

function format_axes_local(ax, style)
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.TickDir = 'in';
ax.FontName = style.fontName;
ax.FontSize = style.tickFontSize;
ax.LineWidth = 0.8;
end

function export_figure_local(fig, pngFile)
[folder, name] = fileparts(pngFile);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
end
