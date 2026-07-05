%% Step04_Main_Visualize_GradientXRange030_20250527
% Visualize the main ratio=0.30 gradient-xrange identification result.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(route_dir, 'output', 'identification');
figure_dir = fullfile(route_dir, 'output', 'figures', 'main_gradient_xrange030');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

result_suffix = strtrim(getenv('STEP03_CURRENT_RESULT_SUFFIX'));
if isempty(result_suffix)
    result_suffix = 'Main_GradientXRange030';
end
result_candidates = {
    fullfile(result_dir, sprintf( ...
        'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_%s_20250527.mat', result_suffix))
    fullfile(result_dir, sprintf( ...
        'Result_Step03_04_FirstOrderVPAdaptive_B1_S136_%s_20250527.mat', result_suffix))
    fullfile(result_dir, sprintf( ...
        'Result_Step03_01_DirectLowSpeedWaveform_B1_S136_%s_20250527.mat', result_suffix))
    };
result_file = '';
for i = 1:numel(result_candidates)
    if exist(result_candidates{i}, 'file') == 2
        result_file = result_candidates{i};
        break;
    end
end
if exist(result_file, 'file') ~= 2
    error(['Main Step03 result not found. Expected first candidate:\n  %s\n' ...
        'Run Step03_Main_Run_GradientXRange030_Identification_20250527 first.'], result_candidates{1});
end

loaded = load(result_file, 'Result');
Result = loaded.Result;
T = Result.Trend;

fig = figure('Name', '20250527 main identification trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(T.window_id, T.EO_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5); hold on;
yline(mode(T.EO_id), ':', 'LineWidth', 1.0);
ylabel('EO');
title('Engine order');
style_axes_local();

nexttile;
plot(T.window_id, T.fn_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('Frequency (Hz)');
title('Identified frequency');
style_axes_local();

nexttile;
yyaxis left;
plot(T.window_id, T.weighted_voltage_rmse, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('RMSE (V)');
yyaxis right;
plot(T.window_id, T.A_id, 's-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('A (mm)');
xlabel('Window');
title('Fit quality and amplitude');
style_axes_local();

png_file = fullfile(figure_dir, 'Step04_Main_GradientXRange030_Trend_20250527.png');
pdf_file = fullfile(figure_dir, 'Step04_Main_GradientXRange030_Trend_20250527.pdf');
exportgraphics(fig, png_file, 'Resolution', 300);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');

summary = table();
summary.ResultFile = string(result_file);
summary.DominantEO = mode(T.EO_id);
summary.EOConsistency = mean(T.EO_id == summary.DominantEO, 'omitnan');
summary.MeanFrequencyHz = mean(T.fn_id, 'omitnan');
summary.MedianFrequencyHz = median(T.fn_id, 'omitnan');
summary.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
summary.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
summary.MeanAmplitudeMM = mean(T.A_id, 'omitnan');
summary_csv = fullfile(result_dir, 'Step04_Main_GradientXRange030_Summary_20250527.csv');
writetable(summary, summary_csv);

fprintf('Loaded result: %s\n', result_file);
fprintf('Saved figure: %s\n', png_file);
fprintf('Saved summary: %s\n', summary_csv);
disp(summary);

function style_axes_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end
