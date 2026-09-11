function Analyze_Resonance_TimePoints_StrainRatio_20241106()
%ANALYZE_RESONANCE_TIMEPOINTS_STRAINRATIO_20241106 Resonance-region audit.
% Reads consistently configured diagnostic BTT results and
% the frozen strain evidence. Formal programs and results are not changed.

clc;
close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
validationRoot = fileparts(rootDir);
oldRoute = fullfile(validationRoot, ...
    '20241106_low_speed_gap_prior_decoupling');
sourceRoot = fullfile(oldRoute, ...
    'diagnostics_weight_objective_ablation_20260714', 'results');
outDir = fullfile(rootDir, 'diagnostics', ...
    'resonance_timepoint_strain_ratio_20260715');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

timesSec = [75.0; 76.0; 78.0; 80.0; 81.0; 82.0; 83.0; 84.0; 84.5; 85.5; 86.6];
resultFiles = {
    fullfile(sourceRoot, 'unlocked_topk_weighted_T075p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_unlocked_topk_weighted_T075p000.mat')
    fullfile(sourceRoot, 'unlocked_topk_weighted_T076p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_unlocked_topk_weighted_T076p000.mat')
    fullfile(sourceRoot, 'unlocked_topk_weighted_T078p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_unlocked_topk_weighted_T078p000.mat')
    fullfile(sourceRoot, 'unlocked_topk_weighted_T080p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_T080p000_unlockedTopK_weighted.mat')
    fullfile(outDir, 'btt_weighted_T081p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_resonance_T081p000.mat')
    fullfile(outDir, 'btt_weighted_T082p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_resonance_T082p000.mat')
    fullfile(outDir, 'btt_weighted_T083p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_resonance_T083p000.mat')
    fullfile(outDir, 'btt_weighted_T084p000', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_resonance_T084p000.mat')
    fullfile(outDir, 'btt_weighted_T084p500', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_resonance_T084p500.mat')
    fullfile(sourceRoot, 'unlocked_topk_weighted_T085p500', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_T085p500_unlockedTopK_weighted.mat')
    fullfile(sourceRoot, 'unlocked_topk_weighted', ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_T086p600_unlockedTopK_weighted.mat')
    };

cfg = Config_20241106();
step04File = fullfile(cfg.paths.strainEvidence, ...
    'Step04_BTT_STE_Resonance_Regions_20241106.mat');
alignmentFile = fullfile(cfg.paths.strainValidation, ...
    'Step12_GlobalTimeAlignment_20241106_B4_S257.mat');
S4 = load(step04File, 'strain');
SA = load(alignmentFile, 'Alignment');
transferRatioPerM = mean([1.1308777609, 1.13085644867]);
strainToMm = 1 / transferRatioPerM / 1000;
strainTime = S4.strain.time(:) + SA.Alignment.strainToBttOffsetSec;
strainMm = S4.strain.value(:) * strainToMm;

windowTables = cell(numel(timesSec), 1);
for iCase = 1:numel(timesSec)
    resultFile = resultFiles{iCase};
    assert(isfile(resultFile), 'Missing diagnostic BTT result: %s', resultFile);
    S = load(resultFile, 'Result');
    W = S.Result.WindowResult(:);
    rows = cell(numel(W), 1);
    bounds = vertcat(W.timeWindow);
    spanMask = strainTime >= min(bounds(:, 1)) & strainTime <= max(bounds(:, 2));
    globalFit = search_tone_local(strainTime(spanMask), strainMm(spanMask), ...
        (600:0.02:650).');
    for iWin = 1:numel(W)
        B = W(iWin).modelFits.gap_only;
        mask = strainTime >= W(iWin).timeWindow(1) & ...
            strainTime <= W(iWin).timeWindow(2);
        localFit = search_tone_local(strainTime(mask), strainMm(mask), ...
            (globalFit.frequencyHz - 3:0.02:globalFit.frequencyHz + 3).');
        rows{iWin} = table(timesSec(iCase), W(iWin).windowId, B.EO, ...
            B.freqHz, B.amplitudeMm, localFit.frequencyHz, ...
            localFit.amplitudeMm, localFit.rSquared, ...
            B.amplitudeMm / localFit.amplitudeMm, ...
            'VariableNames', {'case_time_s', 'window_id', 'EO', ...
            'btt_frequency_hz', 'btt_amplitude_mm', ...
            'strain_frequency_hz', 'strain_amplitude_mm', ...
            'strain_single_tone_r_squared', 'btt_to_strain_ratio'});
    end
    windowTables{iCase} = vertcat(rows{:});
end

T = vertcat(windowTables{:});
T.used_for_amplitude_ratio = T.EO == 12;
Tfit = T(T.used_for_amplitude_ratio, :);
globalScale = (Tfit.strain_amplitude_mm' * Tfit.btt_amplitude_mm) / ...
    (Tfit.strain_amplitude_mm' * Tfit.strain_amplitude_mm);
mask75 = Tfit.case_time_s == 75.0;
scale75 = (Tfit.strain_amplitude_mm(mask75)' * Tfit.btt_amplitude_mm(mask75)) / ...
    (Tfit.strain_amplitude_mm(mask75)' * Tfit.strain_amplitude_mm(mask75));

summaryRows = cell(numel(timesSec), 1);
for iCase = 1:numel(timesSec)
    C = T(T.case_time_s == timesSec(iCase), :);
    Ca = C(C.used_for_amplitude_ratio, :);
    caseScale = (Ca.strain_amplitude_mm' * Ca.btt_amplitude_mm) / ...
        (Ca.strain_amplitude_mm' * Ca.strain_amplitude_mm);
    casePrediction = caseScale * Ca.strain_amplitude_mm;
    globalPrediction = globalScale * Ca.strain_amplitude_mm;
    prediction75 = scale75 * Ca.strain_amplitude_mm;
    summaryRows{iCase} = table(timesSec(iCase), height(C), nnz(C.EO == 12), ...
        mean(Ca.btt_frequency_hz), mean(Ca.strain_frequency_hz), ...
        mean(Ca.btt_amplitude_mm), mean(Ca.strain_amplitude_mm), ...
        mean(Ca.btt_to_strain_ratio), std(Ca.btt_to_strain_ratio), ...
        100 * std(Ca.btt_to_strain_ratio) / mean(Ca.btt_to_strain_ratio), ...
        min(Ca.btt_to_strain_ratio), max(Ca.btt_to_strain_ratio), ...
        corr(Ca.strain_amplitude_mm, Ca.btt_amplitude_mm), ...
        mean(Ca.strain_single_tone_r_squared), caseScale, ...
        transferRatioPerM / caseScale, ...
        mape_local(casePrediction, Ca.btt_amplitude_mm), ...
        mape_local(globalPrediction, Ca.btt_amplitude_mm), ...
        mape_local(prediction75, Ca.btt_amplitude_mm), ...
        'VariableNames', {'case_time_s', 'window_count', 'EO12_count', ...
        'mean_btt_frequency_hz', 'mean_strain_frequency_hz', ...
        'mean_btt_amplitude_mm', 'mean_strain_amplitude_mm', ...
        'mean_btt_to_strain_ratio', 'std_btt_to_strain_ratio', ...
        'cv_btt_to_strain_ratio_percent', 'min_btt_to_strain_ratio', ...
        'max_btt_to_strain_ratio', 'within_case_trend_correlation', ...
        'mean_strain_single_tone_r_squared', 'case_best_scale', ...
        'case_implied_transfer_ratio_1_per_m', ...
        'case_specific_scale_mape_percent', ...
        'all_case_global_scale_mape_percent', ...
        'scale_fixed_from_75s_mape_percent'});
end
Summary = vertcat(summaryRows{:});
T.global_scale = repmat(globalScale, height(T), 1);
T.scale_fixed_from_75s = repmat(scale75, height(T), 1);
T.global_scaled_strain_mm = globalScale * T.strain_amplitude_mm;
T.scaled_from_75s_strain_mm = scale75 * T.strain_amplitude_mm;

writetable(T, fullfile(outDir, ...
    'ResonanceTimePoints_WindowDetail_20241106.csv'));
writetable(Summary, fullfile(outDir, ...
    'ResonanceTimePoints_Summary_20241106.csv'));
save(fullfile(outDir, 'ResonanceTimePoints_Audit_20241106.mat'), ...
    'T', 'Summary', 'globalScale', 'scale75', 'transferRatioPerM');
plot_summary_local(T, Summary, globalScale, scale75, outDir);

fprintf('%d-case global strain scale = %.6f; implied transfer ratio = %.6f 1/m.\n', ...
    numel(timesSec), globalScale, transferRatioPerM / globalScale);
fprintf('75 s anchored strain scale = %.6f; implied transfer ratio = %.6f 1/m.\n', ...
    scale75, transferRatioPerM / scale75);
disp(Summary);
fprintf('\nOutputs: %s\n', outDir);
end


function fit = search_tone_local(t, y, frequencyGrid)
t = t(:);
y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
assert(numel(t) >= 20, 'Fewer than 20 valid strain samples.');
t = t - mean(t);
y = detrend(y);
bestSse = inf;
fit = struct();
for i = 1:numel(frequencyGrid)
    phase = 2 * pi * frequencyGrid(i) * t;
    design = [cos(phase), sin(phase), ones(size(t))];
    coef = design \ y;
    residual = y - design * coef;
    sse = sum(residual .^ 2);
    if sse < bestSse
        bestSse = sse;
        totalSse = sum((y - mean(y)) .^ 2);
        fit.frequencyHz = frequencyGrid(i);
        fit.amplitudeMm = hypot(coef(1), coef(2));
        fit.rSquared = 1 - sse / max(totalSse, eps);
    end
end
end


function value = mape_local(prediction, reference)
value = 100 * mean(abs(prediction - reference) ./ max(abs(reference), eps));
end


function plot_summary_local(T, S, globalScale, scale75, outDir)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 22, 15], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
plot(ax, S.case_time_s, S.mean_btt_amplitude_mm, '-o', 'LineWidth', 1.2, ...
    'DisplayName', 'BTT');
hold(ax, 'on');
plot(ax, S.case_time_s, S.mean_strain_amplitude_mm, '-s', 'LineWidth', 1.2, ...
    'DisplayName', 'Strain-derived');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Mean amplitude (mm)');
title(ax, 'Resonance-region amplitude');
legend(ax, 'Location', 'best', 'Box', 'off'); grid(ax, 'on');

ax = nexttile;
errorbar(ax, S.case_time_s, S.mean_btt_to_strain_ratio, ...
    S.std_btt_to_strain_ratio, '-o', 'LineWidth', 1.2);
hold(ax, 'on');
yline(ax, globalScale, '--', 'All-case scale');
yline(ax, scale75, ':', '75 s scale');
xlabel(ax, 'Time (s)'); ylabel(ax, 'BTT/strain amplitude ratio');
title(ax, 'Required scale is not constant'); grid(ax, 'on');

ax = nexttile;
plot(ax, S.case_time_s, S.case_implied_transfer_ratio_1_per_m, ...
    '-o', 'LineWidth', 1.2);
yline(ax, 1.1308671048, '--', 'FE transfer ratio');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Implied transfer ratio (1/m)');
title(ax, 'Transfer ratio required to force agreement'); grid(ax, 'on');

ax = nexttile;
bar(ax, [S.case_specific_scale_mape_percent, ...
    S.all_case_global_scale_mape_percent, ...
    S.scale_fixed_from_75s_mape_percent]);
set(ax, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), S.case_time_s, ...
    'UniformOutput', false));
xlabel(ax, 'Time (s)'); ylabel(ax, 'MAPE relative to BTT (%)');
title(ax, 'Scaling error by time point');
legend(ax, {'Case-specific', 'All-case global', 'Fixed from 75 s'}, ...
    'Location', 'best', 'Box', 'off'); grid(ax, 'on');

exportgraphics(fig, fullfile(outDir, ...
    'ResonanceTimePoints_StrainRatio_20241106.png'), 'Resolution', 300);
close(fig);
end
