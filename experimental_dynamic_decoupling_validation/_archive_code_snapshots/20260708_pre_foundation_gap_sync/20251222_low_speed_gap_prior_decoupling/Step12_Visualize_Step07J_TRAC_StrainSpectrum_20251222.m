%% Compare main Step07J result with strain-gauge waveform and spectrum
% Current validation route:
%   1. use only the main method, gap_tilt;
%   2. compare identified frequency with the 50-55 s strain local FFT;
%   3. compute TRAC with the recorded time alignment and with a small
%      local strain-window shift search, to diagnose possible time offset.
%
% Important:
%   The strain gauge signal is converted to equivalent displacement in mm
%   using the same transfer ratio as Step7_true. The waveform comparison
%   and TRAC calculation use mm directly; no amplitude normalization is
%   applied before TRAC.

clear; clc; close all;

%% 1. User-adjustable settings
thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_trac_strain_spectrum');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

datasetLabel = sprintf('%s %s', C0.dataset, C0.caseTag);
selectedTimeSec = [50, 55];
methodName = 'gap_tilt';
methodLabel = 'gap + tilt';
timeShiftSearchSec = -0.003:0.00005:0.003;
K_Strain2MM = 1 / 1.266 / 1000;  % Strain-to-displacement transfer ratio, same as Step7_true.
plotFreqHz = [0, 1000];
peakSearchHz = [450, 700];
primaryStrainChannelFile = 'AI1-01_20251222204737.mat';
primaryStrainChannelLabel = 'AI1-01';
legacyStrainTimeOffsetSec = -30.1765;  % Offset used by the original Step7_true script.
legacyResultFiles = {'Vib_Ident_B1_Start39p2s.mat', ...
    'Vib_Ident_B1_Start41p7s.mat', 'Vib_Ident_B1_Start50p0s.mat'};
legacyAnalysisWindowRange = 2:16;
legacyOprPinlvHz = 5e6;
legacyOprFileList = 4000:1000:57000;

%% 2. Load Step07J and Step04 strain FFT evidence
step07jFile = resolve_step07j_result_local(outDir, C0);
step04File = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.mat');
if ~isfile(step07jFile)
    error('Missing Step07J result: %s', step07jFile);
end
if ~isfile(step04File)
    error('Run Step04 first. Missing file: %s', step04File);
end

S7 = load(step07jFile, 'Result');
R7 = S7.Result;
windowResult = R7.WindowResult(:);

fprintf('Using Step07J result: %s\n', step07jFile);

S4 = load(step04File, 'regionTable', 'localFftTable', 'localFftResult', 'strainSource');
if isempty(S4.regionTable)
    error('Step04 selected region table is empty.');
end
strainOrderFreqHz = S4.localFftTable.syncFreqAtRegionOrderHz(1);
primaryStrainSource = override_strain_channel_source_local( ...
    S4.strainSource, primaryStrainChannelFile);

%% 3. Load aligned strain-gauge time series
[strainTime, strainMm, strainFile] = load_aligned_strain_mm_local(primaryStrainSource, K_Strain2MM);
[strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, strainFftPeakAmpMm] = ...
    compute_selected_strain_fft_mm_local(strainTime, strainMm, selectedTimeSec, peakSearchHz);

%% 4. Compute raw and shift-optimized TRAC
rows = {};
detailData = {};
detailPickWindow = [1, round(numel(windowResult) / 2), numel(windowResult)];

for w = 1:numel(windowResult)
    tWin = windowResult(w).timeWindow;
    if tWin(2) < selectedTimeSec(1) || tWin(1) > selectedTimeSec(2)
        continue;
    end

    fit = windowResult(w).modelFits.(methodName);
    [rawTrac, rawPhase, tRaw, yRaw, yReconRaw] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, 0, fit.amplitudeMm, fit.freqHz);
    [bestTrac, bestPhase, bestShift, tBest, yBest, yReconBest] = ...
        search_best_trac_shift_local(strainTime, strainMm, tWin, ...
        timeShiftSearchSec, fit.amplitudeMm, fit.freqHz);

    rows(end + 1, :) = {windowResult(w).windowId, tWin(1), tWin(2), ...
        methodLabel, fit.EO, fit.freqHz, fit.amplitudeMm, fit.weightedRmseMv, ...
        rawPhase, rawTrac, bestPhase, bestTrac, bestShift, ...
        numel(tRaw), rms(yRaw, 'omitnan'), strainFftPeakHz, ...
        strainFftPeakAmpMm, fit.freqHz - strainFftPeakHz}; %#ok<AGROW>

    if ismember(windowResult(w).windowId, detailPickWindow)
        s = struct();
        s.windowId = windowResult(w).windowId;
        s.methodLabel = methodLabel;
        s.timeRaw = tRaw;
        s.strainRaw = yRaw;
        s.reconRaw = yReconRaw;
        s.rawTrac = rawTrac;
        s.timeBest = tBest;
        s.strainBest = yBest;
        s.reconBest = yReconBest;
        s.bestTrac = bestTrac;
        s.bestShift = bestShift;
        s.freqHz = fit.freqHz;
        s.ampMm = fit.amplitudeMm;
        detailData{end + 1} = s; %#ok<SAGROW>
    end
end

tracTable = cell2table(rows, 'VariableNames', {'window_id', ...
    'time_start_s', 'time_end_s', 'method', 'EO', 'frequency_hz', ...
    'amplitude_mm', 'rmse_mV', 'raw_phase_rad', 'raw_TRAC', ...
    'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm', 'strain_fft_peak_hz', ...
    'strain_fft_peak_amp_mm', 'freq_error_vs_strain_fft_hz'});

csvFile = fullfile(outDir, sprintf('Step07J_TRAC_StrainSpectrum_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(tracTable, csvFile);

%% 5. Audit strain-channel and legacy-window sensitivity
strainSources = make_strain_source_options_local(S4.strainSource, legacyStrainTimeOffsetSec);
channelOffsetTable = compute_step07j_channel_offset_table_local(windowResult, ...
    methodName, methodLabel, strainSources, timeShiftSearchSec, selectedTimeSec, ...
    K_Strain2MM);
channelOffsetCsvFile = fullfile(outDir, sprintf( ...
    'Step07J_TRAC_ChannelOffsetSensitivity_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(channelOffsetTable, channelOffsetCsvFile);

legacyDataDir = infer_legacy_dynamic_data_dir_local(S4.strainSource, C0);
legacyTracTable = compute_legacy_start_trac_audit_local(legacyDataDir, ...
    legacyResultFiles, strainSources, legacyAnalysisWindowRange, ...
    legacyOprPinlvHz, legacyOprFileList, K_Strain2MM);
legacyCsvFile = fullfile(outDir, sprintf( ...
    'Step07J_LegacyStart_TRAC_Audit_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(legacyTracTable, legacyCsvFile);
natureFigFile = plot_ai101_nature_validation_local(tracTable, channelOffsetTable, ...
    detailData, strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, ...
    strainOrderFreqHz, figDir, C0, primaryStrainChannelLabel, selectedTimeSec);

%% 6. Plot spectrum, frequency, TRAC and shift
fig = figure('Name', '20251222 Step07J gap_tilt TRAC and strain spectrum', ...
    'Color', 'w', 'Position', [80, 80, 1360, 900], 'NumberTitle', 'off');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
mainColor = [0.10, 0.55, 0.25];

nexttile;
hold on; grid on; box on;
plot(strainFftFreqHz, strainFftAmpMm, ...
    'Color', [0.10, 0.35, 0.75], 'LineWidth', 1.1, ...
    'DisplayName', 'strain local FFT, converted to mm');
xline(strainFftPeakHz, 'r--', sprintf('strain FFT peak %.3f Hz', strainFftPeakHz), ...
    'LineWidth', 1.2, 'DisplayName', 'strain FFT peak');
xline(strainOrderFreqHz, 'k:', sprintf('14X order line %.3f Hz', strainOrderFreqHz), ...
    'LineWidth', 1.0, 'DisplayName', '14X order line');
xlim(plotFreqHz);
xlabel('Frequency (Hz)');
ylabel('Equivalent displacement FFT amp. (mm)');
title(sprintf('%s strain spectrum, %.1f-%.1f s, K=1/1.266/1000, source %s', ...
    datasetLabel, selectedTimeSec(1), selectedTimeSec(2), short_name_local(strainFile)));
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(tracTable.window_id, tracTable.frequency_hz, '-o', ...
    'Color', mainColor, 'LineWidth', 1.2, 'MarkerSize', 4, ...
    'DisplayName', methodLabel);
yline(strainFftPeakHz, 'r--', 'DisplayName', 'strain FFT peak');
xlabel('Step07J window');
ylabel('Frequency (Hz)');
title('Main-method frequency compared with strain-gauge FFT');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(tracTable.window_id, tracTable.raw_TRAC, '-o', ...
    'Color', [0.45, 0.45, 0.45], 'LineWidth', 1.0, 'MarkerSize', 4, ...
    'DisplayName', 'raw time alignment');
plot(tracTable.window_id, tracTable.best_TRAC, '-o', ...
    'Color', mainColor, 'LineWidth', 1.3, 'MarkerSize', 4, ...
    'DisplayName', 'best local time shift');
xlabel('Step07J window');
ylabel('TRAC');
ylim([0, 1]);
title(sprintf('TRAC for %s, shift search %.1f to %.1f ms', ...
    methodLabel, 1000 * min(timeShiftSearchSec), 1000 * max(timeShiftSearchSec)));
legend('Location', 'best');

nexttile;
yyaxis left;
plot(tracTable.window_id, tracTable.rmse_mV, '-o', ...
    'Color', mainColor, 'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Step07J waveform RMSE (mV)');
yyaxis right;
plot(tracTable.window_id, 1000 * tracTable.best_time_shift_s, '-s', ...
    'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.1, 'MarkerSize', 4);
ylabel('Best time shift (ms)');
grid on; box on;
xlabel('Step07J window');
title('Waveform objective and local strain-time shift');

figFile = fullfile(figDir, sprintf('Step07J_TRAC_StrainSpectrum_%s_%s.png', C0.dataset, C0.caseTag));
saveas(fig, figFile);

%% 7. Plot channel/offset sensitivity and legacy audit
[channelOffsetFigFile, legacyFigFile] = plot_trac_audit_figures_local( ...
    channelOffsetTable, legacyTracTable, figDir, C0, selectedTimeSec);

%% 8. Plot representative waveform details
detailFig = figure('Name', '20251222 Step07J gap_tilt TRAC waveform details', ...
    'Color', 'w', 'Position', [100, 100, 1320, 760], 'NumberTitle', 'off');
nDetail = numel(detailData);
tiledlayout(max(1, nDetail), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nDetail
    s = detailData{i};
    nexttile;
    hold on; grid on; box on;
    plot(s.timeRaw, s.strainRaw, ...
        'Color', [0.45, 0.45, 0.45], 'LineWidth', 1.0, ...
        'DisplayName', 'strain gauge displacement (mm)');
    plot(s.timeRaw, s.reconRaw, 'r--', ...
        'LineWidth', 0.9, 'DisplayName', 'BTT sine, raw time, amp in mm');
    plot(s.timeBest, s.reconBest, '-', ...
        'Color', mainColor, 'LineWidth', 0.9, ...
        'DisplayName', 'BTT sine, shifted strain window, amp in mm');
    xlabel('BTT-aligned time (s)');
    ylabel('Equivalent displacement (mm)');
    title(sprintf('W%d %s, f=%.3f Hz, raw TRAC=%.4f, shifted TRAC=%.4f, dt=%.2f ms', ...
        s.windowId, s.methodLabel, s.freqHz, s.rawTrac, ...
        s.bestTrac, 1000 * s.bestShift));
    legend('Location', 'best');
end
detailFigFile = fullfile(figDir, ...
    sprintf('Step07J_TRAC_WaveformDetails_%s_%s.png', C0.dataset, C0.caseTag));
saveas(detailFig, detailFigFile);

fprintf('TRAC and strain spectrum comparison complete.\n');
fprintf('Main method: %s. Selected %.3f-%.3f s.\n', methodLabel, ...
    selectedTimeSec(1), selectedTimeSec(2));
fprintf('Strain FFT peak %.6f Hz, peak amplitude %.6g mm, order line %.6f Hz.\n', ...
    strainFftPeakHz, strainFftPeakAmpMm, strainOrderFreqHz);
disp(groupsummary(tracTable, 'method', 'mean', ...
    {'raw_TRAC', 'best_TRAC', 'best_time_shift_s', ...
    'freq_error_vs_strain_fft_hz', 'rmse_mV'}));
fprintf('Saved:\n  %s\n  %s\n  %s\n', csvFile, figFile, detailFigFile);
fprintf('Additional audit outputs:\n  %s\n  %s\n  %s\n  %s\n', ...
    channelOffsetCsvFile, legacyCsvFile, channelOffsetFigFile, legacyFigFile);
fprintf('Primary AI1-01 publication figure:\n  %s\n', natureFigFile);

function [strainTime, strainMm, strainFile] = load_aligned_strain_mm_local(strainSource, strainToMm)
strainFile = strainSource.strain_file;
if ~isfile(strainFile)
    error('Missing strain file recorded by Step04/Step03: %s', strainFile);
end
S = load(strainFile, 'Datas');
if ~isfield(S, 'Datas') || size(S.Datas, 2) < 2
    error('Invalid strain Datas in file: %s', strainFile);
end
strainTime = S.Datas(:, 1) + strainSource.strain_time_offset_total_sec;
strainMm = (S.Datas(:, 2) - mean(S.Datas(:, 2), 'omitnan')) * strainToMm;
valid = isfinite(strainTime) & isfinite(strainMm);
strainTime = strainTime(valid);
strainMm = strainMm(valid);
[strainTime, orderIdx] = sort(strainTime);
strainMm = strainMm(orderIdx);
end

function strainSource = override_strain_channel_source_local(strainSource, channelFileName)
strainDir = fileparts(strainSource.strain_file);
candidate = fullfile(strainDir, channelFileName);
if ~isfile(candidate)
    error('Requested primary strain channel is missing: %s', candidate);
end
strainSource.strain_file = candidate;
end

function [freqHz, ampMm, peakFreqHz, peakAmpMm] = compute_selected_strain_fft_mm_local( ...
    strainTime, strainMm, selectedTimeSec, peakSearchHz)
mask = strainTime >= selectedTimeSec(1) & strainTime <= selectedTimeSec(2);
tSeg = strainTime(mask);
ySeg = strainMm(mask);
if numel(tSeg) < 16
    error('Selected strain segment %.3f-%.3f s has too few samples.', ...
        selectedTimeSec(1), selectedTimeSec(2));
end
ySeg = detrend(ySeg(:) - mean(ySeg, 'omitnan'));
dt = median(diff(tSeg), 'omitnan');
fs = 1 / dt;
n = numel(ySeg);
win = hann_window_local(n);
yWin = ySeg .* win;
nfft = 2^nextpow2(n);
Y = fft(yWin, nfft);
ampMm = abs(Y(1:nfft/2 + 1)) / max(sum(win) / 2, eps);
freqHz = (0:nfft/2).' * fs / nfft;

peakMask = freqHz >= peakSearchHz(1) & freqHz <= peakSearchHz(2);
if ~any(peakMask)
    error('No FFT bins inside peak search band %.1f-%.1f Hz.', ...
        peakSearchHz(1), peakSearchHz(2));
end
freqBand = freqHz(peakMask);
ampBand = ampMm(peakMask);
[peakAmpMm, peakIdx] = max(ampBand);
peakFreqHz = freqBand(peakIdx);
end

function win = hann_window_local(n)
if n <= 1
    win = ones(n, 1);
else
    k = (0:n-1).';
    win = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
end
end

function [bestTrac, bestPhase, bestShift, bestTime, bestY, bestRecon] = ...
    search_best_trac_shift_local(strainTime, strainMm, tWin, shiftList, ampMm, freqHz)
bestTrac = -inf;
bestPhase = NaN;
bestShift = NaN;
bestTime = [];
bestY = [];
bestRecon = [];
for shiftNow = shiftList
    [tracNow, phaseNow, tNow, yNow, reconNow] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, shiftNow, ampMm, freqHz);
    if isfinite(tracNow) && tracNow > bestTrac
        bestTrac = tracNow;
        bestPhase = phaseNow;
        bestShift = shiftNow;
        bestTime = tNow;
        bestY = yNow;
        bestRecon = reconNow;
    end
end
if isinf(bestTrac)
    bestTrac = NaN;
end
end

function [tracValue, phaseRad, tSeg, ySeg, yRecon] = calc_trac_at_shift_local( ...
    strainTime, strainMm, tWin, shiftSec, ampMm, freqHz)
mask = strainTime >= (tWin(1) + shiftSec) & strainTime <= (tWin(2) + shiftSec);
tSeg = strainTime(mask);
ySeg = strainMm(mask);
if numel(tSeg) < 16
    tracValue = NaN;
    phaseRad = NaN;
    yRecon = [];
    return;
end
ySeg = detrend(ySeg(:) - mean(ySeg, 'omitnan'));
tSeg = tSeg(:);

sinBase = sin(2 * pi * freqHz * tSeg);
cosBase = cos(2 * pi * freqHz * tSeg);
a = ySeg' * sinBase;
b = ySeg' * cosBase;
phaseRad = atan2(b, a);
yRecon = ampMm * sin(2 * pi * freqHz * tSeg + phaseRad);

num = (ySeg' * yRecon)^2;
den = (ySeg' * ySeg) * (yRecon' * yRecon);
if den <= 0 || ~isfinite(den)
    tracValue = NaN;
else
    tracValue = num / den;
end
end

function name = short_name_local(pathText)
[~, name, ext] = fileparts(char(pathText));
name = [name, ext];
end

function pngFile = plot_ai101_nature_validation_local(tracTable, channelOffsetTable, ...
    detailData, fftFreqHz, fftAmpMm, strainPeakHz, orderLineHz, figDir, C0, ...
    primaryChannelLabel, selectedTimeSec)
baseName = sprintf('Step07J_AI101_NatureValidation_%s_%s', C0.dataset, C0.caseTag);
pngFile = fullfile(figDir, [baseName, '.png']);

signalColor = [0.000, 0.447, 0.741];
accentColor = [0.835, 0.369, 0.000];
mutedColor = [0.62, 0.62, 0.62];
darkGray = [0.22, 0.22, 0.22];

fig = figure('Name', 'Step07J AI1-01 validation figure', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18.3, 13.2], ...
    'NumberTitle', 'off');
set(fig, 'Renderer', 'painters');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
plot(ax, fftFreqHz, fftAmpMm, 'Color', darkGray, 'LineWidth', 0.9);
xline(ax, strainPeakHz, '-', sprintf('%.2f Hz', strainPeakHz), ...
    'Color', signalColor, 'LineWidth', 1.0, 'LabelVerticalAlignment', 'bottom');
xline(ax, orderLineHz, '--', '14X', 'Color', mutedColor, 'LineWidth', 0.8, ...
    'LabelVerticalAlignment', 'middle');
xlim(ax, [540, 620]);
ylim(ax, [0, max(fftAmpMm(fftFreqHz >= 540 & fftFreqHz <= 620)) * 1.18]);
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'FFT amplitude (mm)');
title(ax, sprintf('%s strain spectrum, %.0f-%.0f s', ...
    primaryChannelLabel, selectedTimeSec(1), selectedTimeSec(2)));
panel_label_local(ax, 'a');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
freqError = tracTable.frequency_hz - strainPeakHz;
plot(ax, tracTable.window_id, freqError, '-o', 'Color', signalColor, ...
    'LineWidth', 1.0, 'MarkerSize', 3.8, 'MarkerFaceColor', 'w');
yline(ax, 0, '-', 'Color', mutedColor, 'LineWidth', 0.7);
yline(ax, orderLineHz - strainPeakHz, '--', 'Color', mutedColor, 'LineWidth', 0.7);
xlim(ax, [min(tracTable.window_id) - 0.5, max(tracTable.window_id) + 0.5]);
ylim(ax, [-0.13, 0.07]);
xlabel(ax, 'Step07J window');
ylabel(ax, '\Deltaf to strain peak (Hz)');
title(ax, 'Identified frequency error');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
hAi101 = plot(ax, tracTable.window_id, tracTable.raw_TRAC, '-o', 'Color', signalColor, ...
    'LineWidth', 1.1, 'MarkerSize', 3.8, 'MarkerFaceColor', 'w', ...
    'DisplayName', 'AI1-01');
hAi101Shift = plot(ax, tracTable.window_id, tracTable.best_TRAC, '-', 'Color', accentColor, ...
    'LineWidth', 1.0, 'DisplayName', 'AI1-01, local shift');
ai103Mask = strcmp(channelOffsetTable.strain_channel, 'AI1-03') & ...
    strcmp(channelOffsetTable.offset_mode, 'Step04 offset');
if any(ai103Mask)
    hAi103 = plot(ax, channelOffsetTable.window_id(ai103Mask), ...
        channelOffsetTable.raw_TRAC(ai103Mask), '-', 'Color', mutedColor, ...
        'LineWidth', 0.9, 'DisplayName', 'AI1-03 reference');
else
    hAi103 = [];
end
yline(ax, 0.90, ':', 'Color', mutedColor, 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
xlim(ax, [min(tracTable.window_id) - 0.5, max(tracTable.window_id) + 2.0]);
ylim(ax, [0.65, 0.95]);
xlabel(ax, 'Step07J window');
ylabel(ax, 'TRAC');
title(ax, 'Waveform-shape agreement');
text(ax, 18.45, 0.898, 'AI1-01', 'Color', hAi101.Color, ...
    'FontName', 'Arial', 'FontSize', 6.5, 'HorizontalAlignment', 'left');
text(ax, 18.45, 0.916, 'local shift', ...
    'Color', hAi101Shift.Color, 'FontName', 'Arial', 'FontSize', 6.5, ...
    'HorizontalAlignment', 'left');
if ~isempty(hAi103)
    text(ax, 18.45, 0.733, 'AI1-03', 'Color', hAi103.Color, ...
        'FontName', 'Arial', 'FontSize', 6.5, 'HorizontalAlignment', 'left');
end
panel_label_local(ax, 'c');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
if isempty(detailData)
    text(ax, 0.5, 0.5, 'No waveform detail', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
else
    s = detailData{max(1, round(numel(detailData) / 2))};
    t0 = min(s.timeBest);
    tMs = 1000 * (s.timeBest - t0);
    yComponent = project_single_frequency_component_local( ...
        s.timeBest, s.strainBest, s.freqHz);
    plot(ax, tMs, s.strainBest, '-', 'Color', [0.78, 0.78, 0.78], ...
        'LineWidth', 0.6);
    plot(ax, tMs, yComponent, '-', 'Color', darkGray, 'LineWidth', 0.9);
    plot(ax, tMs, s.reconBest, '-', 'Color', signalColor, 'LineWidth', 1.0);
    xlabel(ax, 'Time from window start (ms)');
    ylabel(ax, 'Displacement (mm)');
    title(ax, sprintf('Representative window W%d, TRAC %.3f', ...
        s.windowId, s.bestTrac));
    text(ax, 0.03, 0.92, 'raw AI1-01', 'Units', 'normalized', ...
        'Color', [0.55, 0.55, 0.55], 'FontName', 'Arial', 'FontSize', 6.5);
    text(ax, 0.03, 0.84, 'synchronous component', 'Units', 'normalized', ...
        'Color', darkGray, 'FontName', 'Arial', 'FontSize', 6.5);
    text(ax, 0.03, 0.76, 'BTT reconstruction', 'Units', 'normalized', ...
        'Color', signalColor, 'FontName', 'Arial', 'FontSize', 6.5);
end
panel_label_local(ax, 'd');
apply_nature_axes_local(ax);

export_nature_figure_local(fig, fullfile(figDir, baseName));
end

function apply_nature_axes_local(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 7, 'LineWidth', 0.6, ...
    'TickDir', 'out', 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off');
ax.Title.FontSize = 7.5;
ax.Title.FontWeight = 'normal';
ax.XLabel.FontSize = 7.5;
ax.YLabel.FontSize = 7.5;
end

function panel_label_local(ax, labelText)
text(ax, -0.12, 1.08, labelText, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function export_nature_figure_local(fig, baseFile)
exportgraphics(fig, [baseFile, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [baseFile, '.png'], 'Resolution', 600);
print(fig, [baseFile, '.tiff'], '-dtiff', '-r600');
try
    saveas(fig, [baseFile, '.emf']);
catch ME
    warning('EMF export failed for %s: %s', baseFile, ME.message);
end
end

function yComponent = project_single_frequency_component_local(t, y, freqHz)
t = t(:);
y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
if isempty(t)
    yComponent = [];
    return;
end
y = detrend(y - mean(y, 'omitnan'));
sinBase = sin(2 * pi * freqHz * t);
cosBase = cos(2 * pi * freqHz * t);
design = [sinBase, cosBase];
coef = design \ y;
yComponent = design * coef;
end

function strainSources = make_strain_source_options_local(strainSource, legacyOffsetSec)
strainDir = fileparts(strainSource.strain_file);
channelFiles = {'AI1-01_20251222204737.mat', 'AI1-03_20251222204737.mat'};
channelLabels = {'AI1-01', 'AI1-03'};
offsetLabels = {'Step04 offset', 'Step7 true offset'};
offsetValues = [strainSource.strain_time_offset_total_sec, legacyOffsetSec];

strainSources = struct('channel', {}, 'offsetLabel', {}, 'file', {}, ...
    'offsetSec', {}, 'sourceLabel', {});
for ci = 1:numel(channelFiles)
    strainFile = fullfile(strainDir, channelFiles{ci});
    if ~isfile(strainFile)
        warning('Skipping missing strain channel file: %s', strainFile);
        continue;
    end
    for oi = 1:numel(offsetValues)
        s = struct();
        s.channel = channelLabels{ci};
        s.offsetLabel = offsetLabels{oi};
        s.file = strainFile;
        s.offsetSec = offsetValues(oi);
        s.sourceLabel = sprintf('%s, %s', s.channel, s.offsetLabel);
        strainSources(end + 1) = s; %#ok<AGROW>
    end
end
if isempty(strainSources)
    error('No usable strain channel files were found under %s.', strainDir);
end
end

function T = compute_step07j_channel_offset_table_local(windowResult, methodName, ...
    methodLabel, strainSources, shiftList, selectedTimeSec, strainToMm)
rows = {};
for si = 1:numel(strainSources)
    [strainTime, strainMm] = load_strain_file_mm_local( ...
        strainSources(si).file, strainSources(si).offsetSec, strainToMm);
    for w = 1:numel(windowResult)
        tWin = windowResult(w).timeWindow;
        if tWin(2) < selectedTimeSec(1) || tWin(1) > selectedTimeSec(2)
            continue;
        end
        fit = windowResult(w).modelFits.(methodName);
        [rawTrac, rawPhase, tRaw, yRaw] = calc_trac_at_shift_local( ...
            strainTime, strainMm, tWin, 0, fit.amplitudeMm, fit.freqHz);
        [bestTrac, bestPhase, bestShift] = search_best_trac_shift_local( ...
            strainTime, strainMm, tWin, shiftList, fit.amplitudeMm, fit.freqHz);
        rows(end + 1, :) = {windowResult(w).windowId, ...
            tWin(1), tWin(2), methodLabel, fit.EO, fit.freqHz, ...
            fit.amplitudeMm, fit.weightedRmseMv, strainSources(si).channel, ...
            strainSources(si).offsetLabel, strainSources(si).sourceLabel, ...
            strainSources(si).file, strainSources(si).offsetSec, rawPhase, ...
            rawTrac, bestPhase, bestTrac, bestShift, numel(tRaw), ...
            rms(yRaw, 'omitnan')}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'window_id', 'time_start_s', ...
    'time_end_s', 'method', 'EO', 'frequency_hz', 'amplitude_mm', ...
    'rmse_mV', 'strain_channel', 'offset_mode', 'source_label', ...
    'strain_file', 'strain_time_offset_s', 'raw_phase_rad', 'raw_TRAC', ...
    'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm'});
end

function legacyDataDir = infer_legacy_dynamic_data_dir_local(strainSource, C0)
strainDir = fileparts(strainSource.strain_file);
experimentDir = fileparts(strainDir);
dynamicCaseName = C0.flowConfig.identification.dynamicCaseName;
legacyDataDir = fullfile(experimentDir, '传感器数据', dynamicCaseName);
if ~isfolder(legacyDataDir)
    warning('Legacy dynamic data directory was not found: %s', legacyDataDir);
end
end

function T = compute_legacy_start_trac_audit_local(legacyDataDir, legacyFiles, ...
    strainSources, analysisWindows, pinlvHz, oprFileList, strainToMm)
rows = {};
if ~isfolder(legacyDataDir)
    T = empty_legacy_trac_table_local();
    return;
end

availableFiles = legacyFiles(cellfun(@(f) isfile(fullfile(legacyDataDir, f)), legacyFiles));
if isempty(availableFiles)
    warning('No legacy Vib_Ident_B1_Start*.mat files found in %s.', legacyDataDir);
    T = empty_legacy_trac_table_local();
    return;
end

maxStartSec = max(cellfun(@parse_legacy_start_time_sec_local, availableFiles));
maxReadTime = maxStartSec + (max(analysisWindows) * 3 * 0.1) + 8;
oprTimes = reconstruct_legacy_opr_local(legacyDataDir, oprFileList, pinlvHz, maxReadTime);
if isempty(oprTimes)
    warning('Could not reconstruct legacy OPR times from %s.', legacyDataDir);
    T = empty_legacy_trac_table_local();
    return;
end

legacyStrainSources = strainSources(strcmp({strainSources.offsetLabel}, 'Step7 true offset'));
for fi = 1:numel(availableFiles)
    legacyFile = fullfile(legacyDataDir, availableFiles{fi});
    S = load(legacyFile, 'Result_Struct');
    if ~isfield(S, 'Result_Struct') || ~isfield(S.Result_Struct, 'Trends')
        warning('Skipping invalid legacy result: %s', legacyFile);
        continue;
    end
    R = S.Result_Struct;
    startSec = parse_legacy_start_time_sec_local(availableFiles{fi});
    startTag = sprintf('Start%.1fs', startSec);
    baseIdx = find(oprTimes >= startSec, 1, 'first');
    if isempty(baseIdx)
        warning('No OPR pulse after %.3f s for %s.', startSec, availableFiles{fi});
        continue;
    end
    winSize = R.Config.WinSize;
    for si = 1:numel(legacyStrainSources)
        [strainTime, strainMm] = load_strain_file_mm_local( ...
            legacyStrainSources(si).file, legacyStrainSources(si).offsetSec, strainToMm);
        for kk = 1:numel(analysisWindows)
            w = analysisWindows(kk);
            if w > numel(R.Trends.Amp) || w > numel(R.Trends.Freq) || ...
                    ~isfinite(R.Trends.Amp(w)) || ~isfinite(R.Trends.Freq(w))
                continue;
            end
            idxStart = baseIdx + (w - 1);
            idxEnd = baseIdx + (w - 1) + winSize;
            if idxEnd > numel(oprTimes)
                continue;
            end
            tWin = [oprTimes(idxStart), oprTimes(idxEnd)];
            [tracVal, phaseRad, tSeg, ySeg] = calc_trac_at_shift_local( ...
                strainTime, strainMm, tWin, 0, R.Trends.Amp(w), R.Trends.Freq(w));
            rows(end + 1, :) = {availableFiles{fi}, startTag, startSec, ...
                w, tWin(1), tWin(2), legacyStrainSources(si).channel, ...
                legacyStrainSources(si).file, R.Trends.Freq(w), ...
                R.Trends.Amp(w), phaseRad, tracVal, numel(tSeg), ...
                rms(ySeg, 'omitnan')}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    T = empty_legacy_trac_table_local();
else
    T = cell2table(rows, 'VariableNames', {'legacy_file', 'start_tag', ...
        'start_time_s', 'window_id', 'time_start_s', 'time_end_s', ...
        'strain_channel', 'strain_file', 'frequency_hz', 'amplitude_mm', ...
        'phase_rad', 'TRAC', 'strain_sample_count', 'strain_rms_mm'});
end
end

function T = empty_legacy_trac_table_local()
T = cell2table(cell(0, 14), 'VariableNames', {'legacy_file', 'start_tag', ...
    'start_time_s', 'window_id', 'time_start_s', 'time_end_s', ...
    'strain_channel', 'strain_file', 'frequency_hz', 'amplitude_mm', ...
    'phase_rad', 'TRAC', 'strain_sample_count', 'strain_rms_mm'});
end

function [channelOffsetFigFile, legacyFigFile] = plot_trac_audit_figures_local( ...
    channelOffsetTable, legacyTracTable, figDir, C0, selectedTimeSec)
channelOffsetFigFile = fullfile(figDir, sprintf( ...
    'Step07J_TRAC_ChannelOffsetSensitivity_%s_%s.png', C0.dataset, C0.caseTag));
legacyFigFile = fullfile(figDir, sprintf( ...
    'Step07J_LegacyStart_TRAC_Audit_%s_%s.png', C0.dataset, C0.caseTag));

sourceLabels = unique(channelOffsetTable.source_label, 'stable');
colors = lines(max(1, numel(sourceLabels)));
fig = figure('Name', 'Step07J TRAC channel and offset sensitivity', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10], ...
    'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; box on;
for i = 1:numel(sourceLabels)
    mask = strcmp(channelOffsetTable.source_label, sourceLabels{i});
    plot(channelOffsetTable.window_id(mask), channelOffsetTable.raw_TRAC(mask), ...
        '-o', 'LineWidth', 1.1, 'MarkerSize', 4, 'Color', colors(i, :), ...
        'DisplayName', sourceLabels{i});
end
ylabel('TRAC');
ylim([0, 1]);
title(sprintf('Current Step07J windows, %.1f-%.1f s', ...
    selectedTimeSec(1), selectedTimeSec(2)));
lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 8, 'Interpreter', 'none');
apply_paper_axes_local(gca);

nexttile;
hold on; box on;
for i = 1:numel(sourceLabels)
    mask = strcmp(channelOffsetTable.source_label, sourceLabels{i});
    plot(channelOffsetTable.window_id(mask), 1000 * channelOffsetTable.best_time_shift_s(mask), ...
        '-s', 'LineWidth', 1.1, 'MarkerSize', 4, 'Color', colors(i, :), ...
        'DisplayName', sourceLabels{i});
end
xlabel('Step07J window');
ylabel('Best shift (ms)');
title('Local shift selected by TRAC search');
apply_paper_axes_local(gca);
export_figure_multi_local(fig, channelOffsetFigFile);

fig2 = figure('Name', 'Legacy Step7_true start-window TRAC audit', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8], ...
    'NumberTitle', 'off');
hold on; box on;
if ~isempty(legacyTracTable)
    startTags = unique(legacyTracTable.start_tag, 'stable');
    channels = unique(legacyTracTable.strain_channel, 'stable');
    colors = lines(max(1, numel(channels)));
    x = 1:numel(startTags);
    for ci = 1:numel(channels)
        y = nan(size(x));
        for xi = 1:numel(startTags)
            mask = strcmp(legacyTracTable.start_tag, startTags{xi}) & ...
                strcmp(legacyTracTable.strain_channel, channels{ci});
            y(xi) = mean(legacyTracTable.TRAC(mask), 'omitnan');
        end
        plot(x, y, '-o', 'LineWidth', 1.2, 'MarkerSize', 5, ...
            'Color', colors(ci, :), 'DisplayName', channels{ci});
    end
    xlim([0.8, numel(startTags) + 0.2]);
    xticks(x);
    xticklabels(startTags);
    lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal');
    set(lgd, 'FontName', 'Times New Roman', 'FontSize', 8, 'Interpreter', 'none');
else
    text(0.5, 0.5, 'No legacy audit rows', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
end
ylabel('Mean TRAC');
ylim([0, 1]);
title('Legacy Step7 true windows, W2-W16');
apply_paper_axes_local(gca);
export_figure_multi_local(fig2, legacyFigFile);
end

function apply_paper_axes_local(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'TickDir', 'in', 'LineWidth', 0.8, 'Box', 'on', 'XGrid', 'off', ...
    'YGrid', 'off');
ax.Title.FontSize = 9;
ax.XLabel.FontSize = 9;
ax.YLabel.FontSize = 9;
end

function export_figure_multi_local(fig, pngFile)
[folder, name] = fileparts(pngFile);
if ~exist(folder, 'dir')
    mkdir(folder);
end
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
try
    saveas(fig, fullfile(folder, [name, '.emf']));
catch ME
    warning('EMF export failed for %s: %s', name, ME.message);
end
end

function startSec = parse_legacy_start_time_sec_local(fileName)
tokens = regexp(fileName, 'Start(\d+p\d+)s', 'tokens', 'once');
if isempty(tokens)
    error('Cannot parse start time from legacy file name: %s', fileName);
end
startSec = str2double(strrep(tokens{1}, 'p', '.'));
end

function [strainTime, strainMm] = load_strain_file_mm_local(strainFile, offsetSec, strainToMm)
S = load(strainFile, 'Datas');
if ~isfield(S, 'Datas') || size(S.Datas, 2) < 2
    error('Invalid strain Datas in file: %s', strainFile);
end
strainTime = S.Datas(:, 1) + offsetSec;
strainMm = (S.Datas(:, 2) - mean(S.Datas(:, 2), 'omitnan')) * strainToMm;
valid = isfinite(strainTime) & isfinite(strainMm);
strainTime = strainTime(valid);
strainMm = strainMm(valid);
[strainTime, orderIdx] = sort(strainTime);
strainMm = strainMm(orderIdx);
end

function oprTimes = reconstruct_legacy_opr_local(dataDir, fileList, pinlvHz, maxReadTime)
bigT = [];
bigV = [];
lastEnd = 0;
offset = 0;
for fId = fileList
    fn = fullfile(dataDir, sprintf('4-4-%d.mat', fId));
    if ~isfile(fn)
        continue;
    end
    tmp = load(fn);
    if isfield(tmp, 'jilu04')
        ro = tmp.jilu04;
    elseif isfield(tmp, 'jilu4')
        ro = tmp.jilu4;
    else
        continue;
    end
    ro(ro(:, 1) == 0, :) = [];
    to = ro(:, 1) / pinlvHz;
    vo = ro(:, 2);
    if lastEnd > 0
        offset = lastEnd + 1 / pinlvHz - to(1);
    end
    currT = to + offset;
    bigT = [bigT; currT]; %#ok<AGROW>
    bigV = [bigV; vo]; %#ok<AGROW>
    lastEnd = currT(end);
    if lastEnd > maxReadTime
        break;
    end
end
oprTimes = extract_rising_edge_cluster_local(bigT, bigV, 1.5, 1000);
end

function edgeTimes = extract_rising_edge_cluster_local(t, v, threshold, gapPts)
edgeTimes = [];
above = find(v > threshold);
if isempty(above)
    return;
end
jumps = find(diff(above) > gapPts);
starts = [above(1); above(jumps + 1)];
for k = 1:numel(starts)
    idx = starts(k);
    if idx <= 1
        tc = t(idx);
    else
        t1 = t(idx - 1);
        v1 = v(idx - 1);
        t2 = t(idx);
        v2 = v(idx);
        if v2 ~= v1
            tc = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
        else
            tc = t1;
        end
    end
    edgeTimes(end + 1, 1) = tc; %#ok<AGROW>
end
end

function step07jFile = resolve_step07j_result_local(outDir, C0)
regionTag = '';
if isfield(C0, 'flowConfig') && isfield(C0.flowConfig, 'identification') && ...
        isfield(C0.flowConfig.identification, 'resonanceRegionShortTag')
    regionTag = C0.flowConfig.identification.resonanceRegionShortTag;
end

candidates = {
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_%s_main_gaptilt.mat', C0.dataset, C0.bladeId, C0.sensorTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaptilt.mat', C0.dataset, C0.bladeId, C0.sensorTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_method_compare.mat', C0.dataset, C0.bladeId, C0.sensorTag)
    };
candidates = candidates(~cellfun(@(s) contains(s, '__'), candidates));
for i = 1:numel(candidates)
    candidate = fullfile(outDir, candidates{i});
    if isfile(candidate)
        step07jFile = candidate;
        return;
    end
end

files = dir(fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s*.mat', C0.dataset, C0.caseTag)));
if isempty(files)
    files = dir(fullfile(outDir, sprintf( ...
        'Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s*.mat', ...
        C0.dataset, C0.bladeId, C0.sensorTag)));
end
if ~isempty(files)
    [~, idx] = max([files.datenum]);
    step07jFile = fullfile(files(idx).folder, files(idx).name);
    warning('Using latest matching Step07J result because canonical files were not found:\n  %s', step07jFile);
    return;
end

error(['Missing Step07J result for %s. Expected a main_gaptilt or method_compare file under:\n' ...
    '  %s\nRun Step07J/Step00 first.'], C0.caseTag, outDir);
end
