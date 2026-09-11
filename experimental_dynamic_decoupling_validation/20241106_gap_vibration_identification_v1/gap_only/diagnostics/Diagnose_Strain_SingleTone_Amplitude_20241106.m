function Diagnose_Strain_SingleTone_Amplitude_20241106()
%DIAGNOSE_STRAIN_SINGLETONE_AMPLITUDE_20241106 Independent amplitude audit.
% Reads the frozen identification results and strain evidence. It does not
% modify Main01-Main13, Config_20241106, or any formal result file.

clc;
close all;

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
cfg = Config_20241106();

caseTimesSec = [75.0, 86.6];
caseTags = {'T075p000', 'T086p600'};
transferRatioPerM = mean([1.1308777609, 1.13085644867]);
strainToMm = 1 / transferRatioPerM / 1000;

outDir = fullfile(cfg.paths.diagnostics, ...
    'strain_single_tone_amplitude_20260715');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

step04File = fullfile(cfg.paths.strainEvidence, ...
    'Step04_BTT_STE_Resonance_Regions_20241106.mat');
alignmentFile = fullfile(cfg.paths.strainValidation, ...
    'Step12_GlobalTimeAlignment_20241106_B4_S257.mat');
assert(isfile(step04File), 'Missing strain evidence: %s', step04File);
assert(isfile(alignmentFile), 'Missing time alignment: %s', alignmentFile);

S4 = load(step04File, 'strain');
SA = load(alignmentFile, 'Alignment');
assert(isfield(S4, 'strain') && isfield(SA, 'Alignment'), ...
    'Required strain or alignment structure is missing.');
offsetSec = SA.Alignment.strainToBttOffsetSec;
strainTime = S4.strain.time(:) + offsetSec;
strainValue = S4.strain.value(:);
strainRaw = S4.strain.rawValue(:);
valid = isfinite(strainTime) & isfinite(strainValue) & isfinite(strainRaw);
strainTime = strainTime(valid);
strainValue = strainValue(valid);
strainRaw = strainRaw(valid);

detailTables = cell(numel(caseTimesSec), 1);
summaryRows = cell(numel(caseTimesSec), 1);
for iCase = 1:numel(caseTimesSec)
    resultFile = resolve_result_file_local(cfg.paths.gapResults, caseTags{iCase});
    S7 = load(resultFile, 'Result');
    R = S7.Result;
    W = R.WindowResult(:);

    bounds = vertcat(W.timeWindow);
    caseSpan = [min(bounds(:, 1)), max(bounds(:, 2))];
    caseMask = strainTime >= caseSpan(1) & strainTime <= caseSpan(2);
    globalFit = search_sse_tone_local(strainTime(caseMask), ...
        strainValue(caseMask), (600:0.02:650).', false);

    rows = cell(numel(W), 1);
    for iWin = 1:numel(W)
        fitBtt = W(iWin).modelFits.gap_only;
        timeWindow = W(iWin).timeWindow;
        mask = strainTime >= timeWindow(1) & strainTime <= timeWindow(2);
        t = strainTime(mask);
        y = strainValue(mask) * strainToMm;
        yRaw = strainRaw(mask) * strainToMm;
        assert(numel(t) >= 20, 'Too few strain samples in case %.1f window %d.', ...
            caseTimesSec(iCase), iWin);

        localGrid = (globalFit.frequencyHz - 3:0.02:globalFit.frequencyHz + 3).';
        main12Fit = search_sse_tone_local(t, y, localGrid, false);
        main13Grid = (fitBtt.freqHz - 80:0.1:fitBtt.freqHz + 80).';
        main13Fit = search_max_amplitude_local(t, y, main13Grid);
        fixedFit = fit_tone_local(t, y, fitBtt.freqHz, false, 0);
        hannFit = search_sse_tone_local(t, y, localGrid, true);
        rawFit = search_sse_tone_local(t, yRaw, localGrid, false);
        chirpFit = search_chirp_local(t, y, main12Fit.frequencyHz);

        centerTime = mean(timeWindow);
        halfDuration = diff(timeWindow) / 4;
        centerMask = strainTime >= centerTime - halfDuration & ...
            strainTime <= centerTime + halfDuration;
        centerFit = search_sse_tone_local(strainTime(centerMask), ...
            strainValue(centerMask) * strainToMm, localGrid, false);

        extendedHalfDuration = 0.75 * diff(timeWindow);
        extendedMask = strainTime >= centerTime - extendedHalfDuration & ...
            strainTime <= centerTime + extendedHalfDuration;
        extendedFit = search_sse_tone_local(strainTime(extendedMask), ...
            strainValue(extendedMask) * strainToMm, localGrid, false);

        [residualPeakHz, residualPeakAmpMm] = residual_peak_local( ...
            main12Fit.timeLocal, main12Fit.residual, [600, 650]);
        methodAmps = [main12Fit.amplitudeMm, main13Fit.amplitudeMm, ...
            fixedFit.amplitudeMm, hannFit.amplitudeMm, rawFit.amplitudeMm, ...
            chirpFit.amplitudeMm, centerFit.amplitudeMm, extendedFit.amplitudeMm];
        spreadPercent = 100 * (max(methodAmps) - min(methodAmps)) / ...
            max(main12Fit.amplitudeMm, eps);

        rows{iWin} = table(caseTimesSec(iCase), W(iWin).windowId, ...
            timeWindow(1), timeWindow(2), numel(t), fitBtt.EO, ...
            fitBtt.freqHz, fitBtt.amplitudeMm, globalFit.frequencyHz, ...
            main12Fit.frequencyHz, main12Fit.amplitudeMm, main12Fit.rSquared, ...
            main13Fit.frequencyHz, main13Fit.amplitudeMm, ...
            fixedFit.amplitudeMm, hannFit.frequencyHz, hannFit.amplitudeMm, ...
            rawFit.amplitudeMm, chirpFit.frequencyHz, chirpFit.chirpRateHzPerSec, ...
            chirpFit.amplitudeMm, chirpFit.rSquared, centerFit.amplitudeMm, ...
            extendedFit.amplitudeMm, residualPeakHz, residualPeakAmpMm, ...
            residualPeakAmpMm / max(main12Fit.amplitudeMm, eps), spreadPercent, ...
            'VariableNames', {'case_time_s', 'window_id', 'time_start_s', ...
            'time_end_s', 'sample_count', 'EO', 'btt_frequency_hz', ...
            'btt_amplitude_mm', 'strain_global_frequency_hz', ...
            'main12_frequency_hz', 'main12_amplitude_mm', 'main12_r_squared', ...
            'main13_frequency_hz', 'main13_amplitude_mm', ...
            'fixed_btt_frequency_amplitude_mm', 'hann_frequency_hz', ...
            'hann_amplitude_mm', 'raw_strain_amplitude_mm', ...
            'chirp_frequency_hz', 'chirp_rate_hz_per_s', ...
            'chirp_amplitude_mm', 'chirp_r_squared', ...
            'center_half_window_amplitude_mm', ...
            'extended_1p5_window_amplitude_mm', 'residual_peak_frequency_hz', ...
            'residual_peak_amplitude_mm', 'residual_to_main12_amp_ratio', ...
            'method_spread_percent'});
    end

    T = vertcat(rows{:});
    detailTables{iCase} = T;
    requiredRatio = transferRatioPerM * mean(T.main12_amplitude_mm, 'omitnan') / ...
        mean(T.btt_amplitude_mm, 'omitnan');
    alternativeMax = max([mean(T.main13_amplitude_mm, 'omitnan'), ...
        mean(T.fixed_btt_frequency_amplitude_mm, 'omitnan'), ...
        mean(T.hann_amplitude_mm, 'omitnan'), ...
        mean(T.raw_strain_amplitude_mm, 'omitnan'), ...
        mean(T.chirp_amplitude_mm, 'omitnan'), ...
        mean(T.center_half_window_amplitude_mm, 'omitnan'), ...
        mean(T.extended_1p5_window_amplitude_mm, 'omitnan')]);
    summaryRows{iCase} = table(caseTimesSec(iCase), height(T), ...
        mode(T.EO), mean(T.btt_frequency_hz, 'omitnan'), ...
        mean(T.btt_amplitude_mm, 'omitnan'), ...
        mean(T.main12_frequency_hz, 'omitnan'), ...
        mean(T.main12_amplitude_mm, 'omitnan'), ...
        mean(T.main13_amplitude_mm, 'omitnan'), ...
        mean(T.fixed_btt_frequency_amplitude_mm, 'omitnan'), ...
        mean(T.hann_amplitude_mm, 'omitnan'), ...
        mean(T.raw_strain_amplitude_mm, 'omitnan'), ...
        mean(T.chirp_amplitude_mm, 'omitnan'), ...
        mean(T.center_half_window_amplitude_mm, 'omitnan'), ...
        mean(T.extended_1p5_window_amplitude_mm, 'omitnan'), ...
        mean(T.main12_r_squared, 'omitnan'), ...
        mean(T.residual_to_main12_amp_ratio, 'omitnan'), ...
        mean(T.method_spread_percent, 'omitnan'), ...
        alternativeMax / mean(T.main12_amplitude_mm, 'omitnan'), ...
        mean(T.btt_amplitude_mm ./ T.main12_amplitude_mm, 'omitnan'), ...
        requiredRatio, ...
        'VariableNames', {'case_time_s', 'window_count', 'dominant_EO', ...
        'mean_btt_frequency_hz', 'mean_btt_amplitude_mm', ...
        'mean_main12_frequency_hz', 'mean_main12_amplitude_mm', ...
        'mean_main13_amplitude_mm', 'mean_fixed_btt_frequency_amplitude_mm', ...
        'mean_hann_amplitude_mm', 'mean_raw_strain_amplitude_mm', ...
        'mean_chirp_amplitude_mm', 'mean_center_half_window_amplitude_mm', ...
        'mean_extended_1p5_window_amplitude_mm', 'mean_main12_r_squared', ...
        'mean_residual_to_main12_amp_ratio', 'mean_method_spread_percent', ...
        'max_alternative_to_main12_ratio', 'mean_btt_to_main12_ratio', ...
        'transfer_ratio_required_to_match_btt_1_per_m'});
end

detailTable = vertcat(detailTables{:});
summaryTable = vertcat(summaryRows{:});
writetable(detailTable, fullfile(outDir, ...
    'StrainSingleToneAmplitude_WindowDetail_20241106.csv'));
writetable(summaryTable, fullfile(outDir, ...
    'StrainSingleToneAmplitude_CaseSummary_20241106.csv'));
save(fullfile(outDir, 'StrainSingleToneAmplitude_Audit_20241106.mat'), ...
    'detailTable', 'summaryTable', 'offsetSec', 'transferRatioPerM', ...
    'strainToMm');
plot_diagnostic_local(detailTable, summaryTable, outDir);

disp(summaryTable);
fprintf('\nDiagnostic outputs: %s\n', outDir);
end


function resultFile = resolve_result_file_local(resultDir, timeTag)
pattern = sprintf(['Step07J_NestedStaticWarp_VPFullWave_20241106_' ...
    'B4_S257_%s_main_gapfixedtilt_W3S1_F300to1000Hz_DFpm2Hz.mat'], timeTag);
files = dir(fullfile(resultDir, pattern));
assert(numel(files) == 1, 'Expected one formal result for %s; found %d.', ...
    timeTag, numel(files));
resultFile = fullfile(files(1).folder, files(1).name);
end


function fit = search_sse_tone_local(t, y, frequencyGrid, useHann)
bestSse = inf;
fit = struct();
for i = 1:numel(frequencyGrid)
    candidate = fit_tone_local(t, y, frequencyGrid(i), useHann, 0);
    if candidate.weightedSse < bestSse
        bestSse = candidate.weightedSse;
        fit = candidate;
    end
end
end


function fit = search_max_amplitude_local(t, y, frequencyGrid)
bestAmp = -inf;
fit = struct();
for i = 1:numel(frequencyGrid)
    candidate = fit_tone_local(t, y, frequencyGrid(i), false, 0);
    if candidate.amplitudeMm > bestAmp
        bestAmp = candidate.amplitudeMm;
        fit = candidate;
    end
end
end


function fit = search_chirp_local(t, y, centerFrequencyHz)
bestSse = inf;
fit = struct();
frequencyGrid = (centerFrequencyHz - 1:0.05:centerFrequencyHz + 1).';
chirpGrid = (-100:10:100).';
for i = 1:numel(frequencyGrid)
    for j = 1:numel(chirpGrid)
        candidate = fit_tone_local(t, y, frequencyGrid(i), false, chirpGrid(j));
        if candidate.weightedSse < bestSse
            bestSse = candidate.weightedSse;
            fit = candidate;
        end
    end
end
end


function fit = fit_tone_local(t, y, frequencyHz, useHann, chirpRateHzPerSec)
t = t(:);
y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
tLocal = t - mean(t);
yDetrended = detrend(y);
phase = 2 * pi * (frequencyHz * tLocal + 0.5 * chirpRateHzPerSec * tLocal .^ 2);
design = [cos(phase), sin(phase), ones(size(tLocal))];
if useHann
    n = numel(tLocal);
    weight = 0.5 - 0.5 * cos(2 * pi * (0:n-1).' / max(n - 1, 1));
else
    weight = ones(size(tLocal));
end
sqrtWeight = sqrt(weight);
coef = (design .* sqrtWeight) \ (yDetrended .* sqrtWeight);
prediction = design * coef;
residual = yDetrended - prediction;
weightedSse = sum(weight .* residual .^ 2);
totalSse = sum(weight .* (yDetrended - ...
    sum(weight .* yDetrended) / max(sum(weight), eps)) .^ 2);
fit = struct('frequencyHz', frequencyHz, ...
    'chirpRateHzPerSec', chirpRateHzPerSec, ...
    'amplitudeMm', hypot(coef(1), coef(2)), ...
    'rSquared', 1 - weightedSse / max(totalSse, eps), ...
    'weightedSse', weightedSse, 'timeLocal', tLocal, ...
    'residual', residual);
end


function [peakHz, peakAmp] = residual_peak_local(t, residual, frequencyBandHz)
n = numel(residual);
window = 0.5 - 0.5 * cos(2 * pi * (0:n-1).' / max(n - 1, 1));
dt = median(diff(t));
nfft = 2 ^ nextpow2(n);
spectrum = fft(residual(:) .* window, nfft);
frequency = (0:nfft/2).' / (dt * nfft);
amplitude = abs(spectrum(1:nfft/2 + 1)) / max(sum(window) / 2, eps);
band = frequency >= frequencyBandHz(1) & frequency <= frequencyBandHz(2);
bandFrequency = frequency(band);
bandAmplitude = amplitude(band);
[peakAmp, index] = max(bandAmplitude);
peakHz = bandFrequency(index);
end


function plot_diagnostic_local(T, S, outDir)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 22, 15], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(5);
for iCase = 1:height(S)
    caseMask = T.case_time_s == S.case_time_s(iCase);
    Tc = T(caseMask, :);
    ax = nexttile(iCase);
    hold(ax, 'on');
    plot(ax, Tc.window_id, Tc.btt_amplitude_mm, '-o', ...
        'Color', colors(1, :), 'DisplayName', 'BTT');
    plot(ax, Tc.window_id, Tc.main12_amplitude_mm, '-s', ...
        'Color', colors(2, :), 'DisplayName', 'Main12 SSE');
    plot(ax, Tc.window_id, Tc.hann_amplitude_mm, '-^', ...
        'Color', colors(3, :), 'DisplayName', 'Hann SSE');
    plot(ax, Tc.window_id, Tc.chirp_amplitude_mm, '-d', ...
        'Color', colors(4, :), 'DisplayName', 'Chirp');
    xlabel(ax, 'Window index');
    ylabel(ax, 'Amplitude (mm)');
    title(ax, sprintf('%.1f s', S.case_time_s(iCase)));
    grid(ax, 'on');
    legend(ax, 'Location', 'best', 'Box', 'off');
end

ax = nexttile(3);
methodMeans = [S.mean_btt_amplitude_mm, S.mean_main12_amplitude_mm, ...
    S.mean_main13_amplitude_mm, S.mean_fixed_btt_frequency_amplitude_mm, ...
    S.mean_hann_amplitude_mm, S.mean_chirp_amplitude_mm, ...
    S.mean_center_half_window_amplitude_mm, ...
    S.mean_extended_1p5_window_amplitude_mm];
bar(ax, methodMeans.');
set(ax, 'XTick', 1:size(methodMeans, 2), 'XTickLabel', ...
    {'BTT', 'Main12', 'Main13', 'Fixed f', 'Hann', 'Chirp', 'Half', '1.5x'}, ...
    'XTickLabelRotation', 30);
ylabel(ax, 'Mean amplitude (mm)');
legend(ax, arrayfun(@(x) sprintf('%.1f s', x), S.case_time_s, ...
    'UniformOutput', false), 'Location', 'best', 'Box', 'off');
grid(ax, 'on');

ax = nexttile(4);
yyaxis(ax, 'left');
bar(ax, S.case_time_s, S.mean_btt_to_main12_ratio, 0.45);
ylabel(ax, 'BTT/Main12 amplitude ratio');
yyaxis(ax, 'right');
plot(ax, S.case_time_s, S.mean_main12_r_squared, '-o', 'LineWidth', 1.2);
ylabel(ax, 'Main12 single-tone R^2');
xlabel(ax, 'Case time (s)');
grid(ax, 'on');
title(ax, 'Amplitude discrepancy and fit quality');

exportgraphics(fig, fullfile(outDir, ...
    'StrainSingleToneAmplitude_Diagnostic_20241106.png'), 'Resolution', 300);
close(fig);
end
