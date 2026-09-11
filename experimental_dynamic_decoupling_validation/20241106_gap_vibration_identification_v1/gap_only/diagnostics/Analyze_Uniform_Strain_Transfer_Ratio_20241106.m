function Analyze_Uniform_Strain_Transfer_Ratio_20241106()
%ANALYZE_UNIFORM_STRAIN_TRANSFER_RATIO_20241106 Test one common scale factor.
% This diagnostic reads prior audit outputs only. Formal programs and
% identification results are not changed.

clc;
close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(rootDir, 'diagnostics', ...
    'strain_single_tone_amplitude_20260715');
inputFile = fullfile(inputDir, ...
    'StrainSingleToneAmplitude_WindowDetail_20241106.csv');
assert(isfile(inputFile), 'Run Diagnose_Strain_SingleTone_Amplitude_20241106 first.');

T = readtable(inputFile);
T = T(isfinite(T.btt_amplitude_mm) & isfinite(T.main12_amplitude_mm) & ...
    T.main12_amplitude_mm > 0, :);
oldTransferRatio = mean([1.1308777609, 1.13085644867]);

allFit = fit_scale_local(T.main12_amplitude_mm, T.btt_amplitude_mm, ...
    oldTransferRatio, "all_windows_one_ratio");
caseTimes = unique(T.case_time_s, 'stable');
caseFits = repmat(allFit, numel(caseTimes), 1);
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    caseFits(i) = fit_scale_local(T.main12_amplitude_mm(mask), ...
        T.btt_amplitude_mm(mask), oldTransferRatio, ...
        sprintf('case_%.1f_s', caseTimes(i)));
end

fitStruct = [allFit; caseFits];
fitTable = struct2table(fitStruct);
writetable(fitTable, fullfile(inputDir, ...
    'UniformTransferRatio_FitSummary_20241106.csv'));

T.global_scale = repmat(allFit.scaleFactor, height(T), 1);
T.global_scaled_strain_amplitude_mm = allFit.scaleFactor * T.main12_amplitude_mm;
T.global_residual_mm = T.btt_amplitude_mm - ...
    T.global_scaled_strain_amplitude_mm;
T.case_scale = nan(height(T), 1);
T.case_scaled_strain_amplitude_mm = nan(height(T), 1);
T.case_residual_mm = nan(height(T), 1);
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    T.case_scale(mask) = caseFits(i).scaleFactor;
    T.case_scaled_strain_amplitude_mm(mask) = ...
        caseFits(i).scaleFactor * T.main12_amplitude_mm(mask);
    T.case_residual_mm(mask) = T.btt_amplitude_mm(mask) - ...
        T.case_scaled_strain_amplitude_mm(mask);
end
writetable(T, fullfile(inputDir, ...
    'UniformTransferRatio_WindowResiduals_20241106.csv'));

plot_result_local(T, allFit, caseFits, caseTimes, inputDir);
disp(fitTable);
fprintf('\nUniform-ratio outputs: %s\n', inputDir);
end


function fit = fit_scale_local(strainAmp, bttAmp, oldRatio, label)
strainAmp = strainAmp(:);
bttAmp = bttAmp(:);
scaleFactor = (strainAmp' * bttAmp) / max(strainAmp' * strainAmp, eps);
prediction = scaleFactor * strainAmp;
residual = bttAmp - prediction;
relativeError = abs(residual) ./ max(abs(bttAmp), eps);
amplitudeRatio = bttAmp ./ strainAmp;
fit = struct();
fit.fit_scope = string(label);
fit.window_count = numel(strainAmp);
fit.scaleFactor = scaleFactor;
fit.original_transfer_ratio_1_per_m = oldRatio;
fit.fitted_transfer_ratio_1_per_m = oldRatio / scaleFactor;
fit.mean_btt_to_strain_ratio = mean(amplitudeRatio);
fit.std_btt_to_strain_ratio = std(amplitudeRatio);
fit.cv_btt_to_strain_ratio_percent = 100 * std(amplitudeRatio) / ...
    max(mean(amplitudeRatio), eps);
fit.min_btt_to_strain_ratio = min(amplitudeRatio);
fit.max_btt_to_strain_ratio = max(amplitudeRatio);
fit.rmse_mm = sqrt(mean(residual .^ 2));
fit.mae_mm = mean(abs(residual));
fit.mape_vs_btt_percent = 100 * mean(relativeError);
fit.max_error_vs_btt_percent = 100 * max(relativeError);
fit.trend_correlation = corr(strainAmp, bttAmp, 'Rows', 'complete');
fit.r_squared_about_btt_mean = 1 - sum(residual .^ 2) / ...
    max(sum((bttAmp - mean(bttAmp)) .^ 2), eps);
end


function plot_result_local(T, allFit, caseFits, caseTimes, outDir)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 22, 15], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(numel(caseTimes));

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    plot(ax, T.window_id(mask), ...
        T.btt_amplitude_mm(mask) ./ T.main12_amplitude_mm(mask), '-o', ...
        'Color', colors(i, :), 'DisplayName', sprintf('%.1f s', caseTimes(i)));
end
hGlobal = yline(ax, allFit.scaleFactor, '--k');
hGlobal.DisplayName = 'Global scale';
xlabel(ax, 'Window index');
ylabel(ax, 'BTT/strain amplitude ratio');
title(ax, 'Ratio must be constant for one transfer ratio');
legend(ax, 'Location', 'best', 'Box', 'off');
grid(ax, 'on');

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    plot(ax, T.window_id(mask), T.btt_amplitude_mm(mask), '-o', ...
        'Color', colors(i, :), 'DisplayName', sprintf('BTT %.1f s', caseTimes(i)));
    plot(ax, T.window_id(mask), ...
        T.global_scaled_strain_amplitude_mm(mask), '--s', ...
        'Color', colors(i, :), 'HandleVisibility', 'off');
end
xlabel(ax, 'Window index');
ylabel(ax, 'Amplitude (mm)');
title(ax, 'Global scaling: solid BTT, dashed strain');
legend(ax, 'Location', 'best', 'Box', 'off');
grid(ax, 'on');

ax = nexttile;
hold(ax, 'on');
maximum = 1.08 * max([T.main12_amplitude_mm; T.btt_amplitude_mm]);
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    scatter(ax, T.main12_amplitude_mm(mask), T.btt_amplitude_mm(mask), ...
        28, colors(i, :), 'filled', 'DisplayName', sprintf('%.1f s', caseTimes(i)));
    xLine = [0, maximum];
    plot(ax, xLine, caseFits(i).scaleFactor * xLine, ':', ...
        'Color', colors(i, :), 'HandleVisibility', 'off');
end
plot(ax, [0, maximum], allFit.scaleFactor * [0, maximum], '--k', ...
    'DisplayName', 'Global scale');
xlim(ax, [0, maximum]);
ylim(ax, [0, maximum]);
axis(ax, 'square');
xlabel(ax, 'Strain-derived amplitude (mm)');
ylabel(ax, 'BTT amplitude (mm)');
title(ax, 'One slope cannot represent both cases');
legend(ax, 'Location', 'best', 'Box', 'off');
grid(ax, 'on');

ax = nexttile;
globalMape = zeros(numel(caseTimes), 1);
caseMape = zeros(numel(caseTimes), 1);
for i = 1:numel(caseTimes)
    mask = T.case_time_s == caseTimes(i);
    globalMape(i) = 100 * mean(abs(T.global_residual_mm(mask)) ./ ...
        T.btt_amplitude_mm(mask));
    caseMape(i) = caseFits(i).mape_vs_btt_percent;
end
bar(ax, [globalMape, caseMape]);
set(ax, 'XTickLabel', arrayfun(@(x) sprintf('%.1f s', x), caseTimes, ...
    'UniformOutput', false));
ylabel(ax, 'MAPE relative to BTT (%)');
title(ax, 'Global versus case-specific scaling');
legend(ax, {'One global ratio', 'Separate case ratios'}, ...
    'Location', 'best', 'Box', 'off');
grid(ax, 'on');

exportgraphics(fig, fullfile(outDir, ...
    'UniformTransferRatio_Diagnostic_20241106.png'), 'Resolution', 300);
close(fig);
end
