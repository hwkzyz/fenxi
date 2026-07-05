%% Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20250527
% Visualize the fixed OPRCenterStd direct low-speed-template result.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(routeDir, 'output', 'identification');
figureDir = fullfile(routeDir, 'output', 'figures', 'main_direct_template_oprcenterstd');
if exist(figureDir, 'dir') ~= 7; mkdir(figureDir); end

resultSuffix = 'Main_DirectTemplate_OPRCenterStd';
resultFile = fullfile(resultDir, sprintf( ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_%s_20250527.mat', resultSuffix));
if exist(resultFile, 'file') ~= 2
    error(['Direct OPRCenterStd Step03 result not found:\n  %s\n' ...
        'Run Step03_Main_Run_OPRCenterStd_DirectTemplate_20250527 first.'], resultFile);
end

loaded = load(resultFile, 'Result');
Result = loaded.Result;
T = Result.Trend;

fig = figure('Name', '20250527 direct OPRCenterStd identification trend', 'Color', 'w', ...
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

pngFile = fullfile(figureDir, 'Step04_DirectTemplate_OPRCenterStd_Trend_20250527.png');
pdfFile = fullfile(figureDir, 'Step04_DirectTemplate_OPRCenterStd_Trend_20250527.pdf');
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

summary = table();
summary.ResultFile = string(resultFile);
summary.DominantEO = mode(T.EO_id);
summary.EOConsistency = mean(T.EO_id == summary.DominantEO, 'omitnan');
summary.MeanFrequencyHz = mean(T.fn_id, 'omitnan');
summary.MedianFrequencyHz = median(T.fn_id, 'omitnan');
summary.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
summary.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
summary.MeanAmplitudeMM = mean(T.A_id, 'omitnan');
summaryCsv = fullfile(resultDir, 'Step04_DirectTemplate_OPRCenterStd_Summary_20250527.csv');
writetable(summary, summaryCsv);

fprintf('Loaded result: %s\n', resultFile);
fprintf('Saved figure: %s\n', pngFile);
fprintf('Saved summary: %s\n', summaryCsv);
disp(summary);

function style_axes_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end
