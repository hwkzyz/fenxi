%% Step03_Visualize_Identification_Results_20250527
% Visualization script for the stabilized Step03 identification results.
%
% Inputs:
%   output/identification/Result_ProposedFixedTemplateSliding_B1_S16_20250527.mat
%   output/identification/Result_ProposedFixedTemplateSliding_B1_S136_20250527.mat
%
% Outputs:
%   output/figures/step03_visualization/*.png
%   output/figures/step03_visualization/*.pdf
%   output/figures/step03_visualization/*.emf

clear; clc; close all;

%% Settings
route_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(route_dir, 'output', 'identification');
figure_dir = fullfile(route_dir, 'output', 'figures', 'step03_visualization');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

sensor_tags = {'S13', 'S16', 'S136'};
target_blade = 1;
example_tag = 'S136';
example_window_id = 9;

font_name = 'Times New Roman';
font_size_axis = 8.5;
font_size_label = 9;
line_width = 1.15;
marker_size = 4.5;
colors = [0.8500 0.3250 0.0980;
          0.0000 0.4470 0.7410;
          0.4660 0.6740 0.1880];

export_png = true;
export_pdf = true;
export_emf = true;

%% Load results
Results = {};
for k = 1:numel(sensor_tags)
    result_file = fullfile(result_dir, sprintf( ...
        'Result_ProposedFixedTemplateSliding_B%d_%s_20250527.mat', ...
        target_blade, sensor_tags{k}));
    if isfile(result_file)
        loaded = load(result_file, 'Result');
        Results{end+1} = loaded.Result; %#ok<SAGROW>
        fprintf('Loaded %s\n', result_file);
    else
        warning('Result file not found: %s', result_file);
    end
end

if isempty(Results)
    error('No Step03 result files were found. Run Step03 first.');
end

%% Figure 1: sliding-window identification trends
fig1 = figure('Name', 'Step03 trend comparison', 'Color', 'w');
set(fig1, 'Units', 'centimeters', 'Position', [2 2 17 15]);
tl = tiledlayout(fig1, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
hold(ax, 'on');
for k = 1:numel(Results)
    R = Results{k};
    plot(ax, R.Trend.window_id, R.Trend.frequency_hz, 'o-', ...
        'Color', colors(k, :), 'LineWidth', line_width, 'MarkerSize', marker_size, ...
        'DisplayName', R.SensorTag);
end
for k = 1:numel(Results)
    R = Results{k};
    if isfield(R, 'SG') && isfield(R.SG, 'frequency_hz')
        yline(ax, R.SG.frequency_hz, ':', 'Color', [0.25 0.25 0.25], ...
            'LineWidth', 0.9, 'DisplayName', 'old SG');
        break;
    end
end
ylabel(ax, '{\it f} (Hz)', 'Interpreter', 'tex');
legend(ax, 'Location', 'eastoutside', 'Box', 'off');
title(ax, 'Frequency');

ax = nexttile(tl);
hold(ax, 'on');
for k = 1:numel(Results)
    R = Results{k};
    plot(ax, R.Trend.window_id, R.Trend.EO, 'o-', ...
        'Color', colors(k, :), 'LineWidth', line_width, 'MarkerSize', marker_size, ...
        'DisplayName', R.SensorTag);
end
yline(ax, 14, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
ylabel(ax, 'EO');
title(ax, 'Engine order');

ax = nexttile(tl);
hold(ax, 'on');
for k = 1:numel(Results)
    R = Results{k};
    plot(ax, R.Trend.window_id, R.Trend.amplitude_mm, 'o-', ...
        'Color', colors(k, :), 'LineWidth', line_width, 'MarkerSize', marker_size, ...
        'DisplayName', R.SensorTag);
end
for k = 1:numel(Results)
    R = Results{k};
    if isfield(R, 'SG') && isfield(R.SG, 'amplitude_mm')
        yline(ax, R.SG.amplitude_mm, ':', 'Color', [0.25 0.25 0.25], ...
            'LineWidth', 0.9);
        break;
    end
end
ylabel(ax, '{\it A} (mm)', 'Interpreter', 'tex');
title(ax, 'Amplitude');

ax = nexttile(tl);
hold(ax, 'on');
for k = 1:numel(Results)
    R = Results{k};
    plot(ax, R.Trend.window_id, R.Trend.weighted_rmse_V, 'o-', ...
        'Color', colors(k, :), 'LineWidth', line_width, 'MarkerSize', marker_size, ...
        'DisplayName', R.SensorTag);
end
xlabel(ax, 'Sliding-window index');
ylabel(ax, 'RMSE (V)');
title(ax, 'Weighted residual');

all_axes = findall(fig1, 'Type', 'axes');
set(all_axes, 'FontName', font_name, 'FontSize', font_size_axis, ...
    'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.75);
set(findall(fig1, 'Type', 'text'), 'FontName', font_name, 'FontSize', font_size_label);

base = fullfile(figure_dir, 'Step03_Fig1_TrendComparison');
if export_png; exportgraphics(fig1, [base, '.png'], 'Resolution', 300); end
if export_pdf; exportgraphics(fig1, [base, '.pdf'], 'ContentType', 'vector'); end
if export_emf; print(fig1, [base, '.emf'], '-dmeta'); end

%% Figure 2: EO lock diagnostic
fig2 = figure('Name', 'Step03 EO lock diagnostic', 'Color', 'w');
set(fig2, 'Units', 'centimeters', 'Position', [2.5 2.5 17 9.2]);
tl = tiledlayout(fig2, 1, numel(Results), 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(Results)
    R = Results{k};
    T = sortrows(R.GlobalEOSelectionTable, 'EO');
    ax = nexttile(tl);
    hold(ax, 'on');
    yyaxis(ax, 'left');
    plot(ax, T.EO, T.median_rel_rmse, 'o-', ...
        'Color', colors(k, :), 'LineWidth', line_width, 'MarkerSize', marker_size);
    ylabel(ax, 'Median relative RMSE');
    yyaxis(ax, 'right');
    stem(ax, T.EO, T.support_count, 'filled', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 0.8, 'MarkerSize', 3.5);
    ylabel(ax, 'Support count');
    xline(ax, R.GlobalEO, ':', 'Color', [0.1 0.1 0.1], 'LineWidth', 0.9);
    xlabel(ax, 'EO candidate');
    title(ax, sprintf('%s: locked EO %d', R.SensorTag, R.GlobalEO));
end

all_axes = findall(fig2, 'Type', 'axes');
set(all_axes, 'FontName', font_name, 'FontSize', font_size_axis, ...
    'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.75);
set(findall(fig2, 'Type', 'text'), 'FontName', font_name, 'FontSize', font_size_label);

base = fullfile(figure_dir, 'Step03_Fig2_EOLockDiagnostic');
if export_png; exportgraphics(fig2, [base, '.png'], 'Resolution', 300); end
if export_pdf; exportgraphics(fig2, [base, '.pdf'], 'ContentType', 'vector'); end
if export_emf; print(fig2, [base, '.emf'], '-dmeta'); end

%% Figure 3: measured and fitted waveforms in one representative window
example_idx = [];
for k = 1:numel(Results)
    if strcmpi(Results{k}.SensorTag, example_tag)
        example_idx = k;
        break;
    end
end
if isempty(example_idx)
    example_idx = 1;
end

R = Results{example_idx};
example_window_id = min(example_window_id, numel(R.WindowResult));
Rw = R.WindowResult(example_window_id);
sensor_ids = unique(Rw.sensor_id(:)).';

fig3 = figure('Name', 'Step03 representative waveform fit', 'Color', 'w');
set(fig3, 'Units', 'centimeters', 'Position', [3 3 17 10.5]);
tl = tiledlayout(fig3, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensor_ids)
    sid = sensor_ids(is);
    ax = nexttile(tl);
    idx = Rw.sensor_id == sid & isfinite(Rw.V_meas) & isfinite(Rw.V_pred);
    [t_sort, order] = sort((Rw.t(idx) - Rw.t0) * 1000);
    V_meas = Rw.V_meas(idx);
    V_pred = Rw.V_pred(idx);
    V_meas = V_meas(order);
    V_pred = V_pred(order);
    plot(ax, t_sort, V_meas, '.', 'Color', [0.65 0.65 0.65], ...
        'MarkerSize', 3.0, 'DisplayName', 'measured');
    hold(ax, 'on');
    plot(ax, t_sort, V_pred, '-', 'Color', colors(min(is, size(colors, 1)), :), ...
        'LineWidth', 1.0, 'DisplayName', 'fit');
    ylabel(ax, sprintf('CH%d (V)', sid));
    if is == 1
        title(ax, sprintf('%s, window %d, EO %.4f, f %.3f Hz', ...
            R.SensorTag, Rw.window_id, Rw.EO, Rw.frequency_hz));
        legend(ax, 'Location', 'eastoutside', 'Box', 'off');
    end
    if is == numel(sensor_ids)
        xlabel(ax, 'Time in window (ms)');
    end
end

all_axes = findall(fig3, 'Type', 'axes');
set(all_axes, 'FontName', font_name, 'FontSize', font_size_axis, ...
    'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.75);
set(findall(fig3, 'Type', 'text'), 'FontName', font_name, 'FontSize', font_size_label);

base = fullfile(figure_dir, sprintf('Step03_Fig3_WindowFit_%s_W%02d', R.SensorTag, Rw.window_id));
if export_png; exportgraphics(fig3, [base, '.png'], 'Resolution', 300); end
if export_pdf; exportgraphics(fig3, [base, '.pdf'], 'ContentType', 'vector'); end
if export_emf; print(fig3, [base, '.emf'], '-dmeta'); end

%% Console summary
fprintf('\nSaved Step03 visualization figures to:\n  %s\n', figure_dir);
for k = 1:numel(Results)
    R = Results{k};
    fprintf('%s: EO %d, f std %.6f Hz, A std %.6f mm, mean RMSE %.6f V\n', ...
        R.SensorTag, R.GlobalEO, std(R.Trend.frequency_hz), ...
        std(R.Trend.amplitude_mm), mean(R.Trend.weighted_rmse_V));
end
