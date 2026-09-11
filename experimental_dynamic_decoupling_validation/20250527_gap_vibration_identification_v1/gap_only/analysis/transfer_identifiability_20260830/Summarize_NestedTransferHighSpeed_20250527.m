%% Summarize P0-P3 high-speed propagation results.
clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
resultDir = fullfile(caseDir, 'results', 'gap_aware');
outputDir = fullfile(thisDir, 'outputs');
modelNames = {'P0_g0','P1_g0_gain','P2_registration','P3_full_transfer'};
suffixes = {'p0_g0','p1_g0_gain','p2_registration','p3_full_transfer'};
rows = cell(numel(modelNames), 1);
windowRows = cell(numel(modelNames), 1);

for im = 1:numel(modelNames)
    matFile = fullfile(resultDir, sprintf( ...
        'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_transfer_%s_20260830.mat', suffixes{im}));
    S = load(matFile, 'Result', 'Trend');
    T = S.Trend;
    n = height(T);
    dg = nan(n, 3);
    for iw = 1:n
        fit = S.Result.WindowResult(iw).modelFits.gap_only;
        dg(iw, :) = fit.deltaGapMm(:).';
    end
    rows{im} = table(string(modelNames{im}), mode(T.gap_EO), nnz(T.gap_EO == 14), ...
        mean(T.gap_frequency_hz), std(T.gap_frequency_hz), ...
        mean(T.gap_amplitude_mm), std(T.gap_amplitude_mm), ...
        mean(T.gap_waveform_rmse_mV), median(T.gap_waveform_rmse_mV), ...
        mean(T.gap_dx_mm), mean(T.gap_mean_delta_gap_mm), ...
        mean(dg(:,1)), mean(dg(:,2)), mean(dg(:,3)), ...
        std(dg(:,1)), std(dg(:,2)), std(dg(:,3)), ...
        all(logical(T.gap_domain_feasible)), max(T.gap_domain_max_violation_mm), ...
        'VariableNames', {'model','dominantEO','eo14WindowCount','meanFrequencyHz','stdFrequencyHz', ...
        'meanAmplitudeMm','stdAmplitudeMm','meanWaveformRmseMv','medianWaveformRmseMv', ...
        'meanDxMm','meanDeltaGapMm','meanDeltaGapCH1Mm','meanDeltaGapCH3Mm','meanDeltaGapCH6Mm', ...
        'stdDeltaGapCH1Mm','stdDeltaGapCH3Mm','stdDeltaGapCH6Mm', ...
        'allDomainFeasible','maxDomainViolationMm'});
    windowRows{im} = table(repmat(string(modelNames{im}), n, 1), T.window_id, T.gap_EO, ...
        T.gap_frequency_hz, T.gap_amplitude_mm, T.gap_waveform_rmse_mV, T.gap_dx_mm, ...
        dg(:,1), dg(:,2), dg(:,3), T.gap_domain_feasible, T.gap_domain_max_violation_mm, ...
        'VariableNames', {'model','windowId','EO','frequencyHz','amplitudeMm','waveformRmseMv', ...
        'dxMm','deltaGapCH1Mm','deltaGapCH3Mm','deltaGapCH6Mm','domainFeasible','domainMaxViolationMm'});
end

Summary = vertcat(rows{:});
Window = vertcat(windowRows{:});
ref = Summary(Summary.model == "P3_full_transfer", :);
Summary.frequencyDifferenceFromP3Hz = Summary.meanFrequencyHz - ref.meanFrequencyHz;
Summary.amplitudeDifferenceFromP3Mm = Summary.meanAmplitudeMm - ref.meanAmplitudeMm;
Summary.rmseDifferenceFromP3Mv = Summary.meanWaveformRmseMv - ref.meanWaveformRmseMv;
Summary.meanDeltaGapDifferenceFromP3Mm = Summary.meanDeltaGapMm - ref.meanDeltaGapMm;

writetable(Summary, fullfile(outputDir, 'NestedTransferHighSpeedSummary_20250527.csv'));
writetable(Window, fullfile(outputDir, 'NestedTransferHighSpeedByWindow_20250527.csv'));
save(fullfile(outputDir, 'NestedTransferHighSpeedSummary_20250527.mat'), 'Summary', 'Window');

fig = figure('Color', 'w', 'Position', [80 80 1280 780]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = strrep(cellstr(Summary.model), '_', '\_');
nexttile; bar(Summary.meanFrequencyHz); grid on; box on;
xticks(1:height(Summary)); xticklabels(labels); ylabel('Frequency (Hz)'); title('High-speed frequency');
nexttile; bar(Summary.meanAmplitudeMm); grid on; box on;
xticks(1:height(Summary)); xticklabels(labels); ylabel('Amplitude (mm)'); title('Bending-vibration amplitude');
nexttile; bar([Summary.meanDeltaGapCH1Mm, Summary.meanDeltaGapCH3Mm, Summary.meanDeltaGapCH6Mm]);
grid on; box on; xticks(1:height(Summary)); xticklabels(labels); ylabel('\Delta g (mm)');
title('Per-sensor gap increment'); legend({'CH1','CH3','CH6'}, 'Location', 'best');
nexttile; bar(Summary.meanWaveformRmseMv); grid on; box on;
xticks(1:height(Summary)); xticklabels(labels); ylabel('Waveform RMSE (mV)'); title('Pure high-speed residual');
exportgraphics(fig, fullfile(outputDir, 'NestedTransferHighSpeedComparison_20250527.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outputDir, 'NestedTransferHighSpeedComparison_20250527.pdf'), 'ContentType', 'vector');
close(fig);

disp(Summary);
fprintf('Outputs: %s\n', outputDir);
