%% Step12: validate Step07J against strain-gauge spectrum and waveform, 20250527
% Main figure logic:
%   a) strain-gauge spectrum around the selected resonance band;
%   b) Step07J frequency error relative to the strain FFT peak;
%   c) TRAC trend using the primary strain channel, with the sibling channel
%      as a quiet reference-source sensitivity check;
%   d) representative-window waveform comparison.

clear; clc; close all;

%% 1. Settings
thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20250527();
C0 = CaseConfig();
formalOutDir = packageCfg.paths.results;
outDir = fullfile(thisDir, 'strain_validation_results');
figDir = fullfile(outDir, 'figures');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

methodName = packageCfg.model.mainGapModel;
methodLabel = 'gap only';
timeShiftSearchSec = -0.003:0.00005:0.003;
K_Strain2MM = 1 / packageCfg.strain.transferRatio1PerM / 1000;
peakSearchHz = packageCfg.strain.frequencyBandHz;

%% 2. Load Step07J and Step04 evidence
step07jFile = fullfile(formalOutDir, 'gap_aware', ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_diag_opr6_unequal_production_fullbundle.mat');
step04File = fullfile(formalOutDir, 'Step04_BTT_STE_Resonance_Regions_20250527.mat');
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

S4 = load(step04File, 'regionTable', 'localFftTable', 'strainSource');
% Validate exactly the time span identified by Step07J. Step04 regions are
% diagnostics and may rank a different resonance elsewhere in the run.
windowStartSec = arrayfun(@(x) x.timeWindow(1), windowResult);
windowEndSec = arrayfun(@(x) x.timeWindow(2), windowResult);
selectedTimeSec = [min(windowStartSec), max(windowEndSec)];
if isfield(R7, 'Trend') && ismember('gap_integer_frequency_hz', ...
        R7.Trend.Properties.VariableNames)
    strainOrderFreqHz = mean(R7.Trend.gap_integer_frequency_hz, 'omitnan');
else
    strainOrderFreqHz = mean(arrayfun(@(x) ...
        x.modelFits.(methodName).EO .* x.rotFreqHz, windowResult), 'omitnan');
end

strainSources = make_strain_source_options_local(S4.strainSource);
primaryIdx = find(strcmpi(string({strainSources.channel}), ...
    string(packageCfg.strain.primaryChannel)), 1, 'first');
if isempty(primaryIdx)
    error('Configured primary strain channel %s was not found.', ...
        packageCfg.strain.primaryChannel);
end
primarySource = S4.strainSource;
primarySource.strain_file = strainSources(primaryIdx).file;
[strainTime, strainMm, strainFile] = load_aligned_strain_mm_local( ...
    primarySource, K_Strain2MM);
[strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, strainFftPeakAmpMm] = ...
    compute_selected_strain_fft_mm_local(strainTime, strainMm, selectedTimeSec, peakSearchHz);
primaryChannel = channel_label_from_file_local(strainFile);

%% 3. Compute primary TRAC table
rows = {};
detailData = {};
detailPickWindow = unique([1, round(numel(windowResult) / 2), R7.BestWindowIndex]);
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
        primaryChannel, rawPhase, rawTrac, bestPhase, bestTrac, bestShift, ...
        numel(tRaw), rms(yRaw, 'omitnan'), strainFftPeakHz, ...
        strainFftPeakAmpMm, fit.freqHz - strainFftPeakHz}; %#ok<AGROW>

    if ismember(windowResult(w).windowId, detailPickWindow)
        s = struct();
        s.windowId = windowResult(w).windowId;
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
    'amplitude_mm', 'rmse_mV', 'strain_channel', 'raw_phase_rad', ...
    'raw_TRAC', 'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm', 'strain_fft_peak_hz', ...
    'strain_fft_peak_amp_mm', 'freq_error_vs_strain_fft_hz'});

csvFile = fullfile(outDir, sprintf('Step12_TRAC_StrainSpectrum_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(tracTable, csvFile);

channelTable = compute_channel_sensitivity_table_local(windowResult, methodName, ...
    methodLabel, strainSources, timeShiftSearchSec, selectedTimeSec, K_Strain2MM);
channelCsvFile = fullfile(outDir, sprintf( ...
    'Step12_TRAC_ChannelSensitivity_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(channelTable, channelCsvFile);

%% 4. Publication-style main figure
natureFigFile = plot_step12_nature_validation_local(tracTable, channelTable, ...
    detailData, strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, ...
    strainOrderFreqHz, figDir, C0, primaryChannel, selectedTimeSec);
[channelTimeFigFile, channelTimeTable] = plot_channel_time_domain_compare_local( ...
    windowResult, methodName, strainSources, channelTable, figDir, C0, ...
    K_Strain2MM, R7.BestWindowIndex);
channelTimeCsvFile = fullfile(outDir, sprintf( ...
    'Step12_ChannelTimeDomainCompare_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(channelTimeTable, channelTimeCsvFile);

fprintf('Step12 20250527 validation complete.\n');
fprintf('Primary strain channel: %s (%s)\n', primaryChannel, strainFile);
fprintf('Selected time %.3f-%.3f s. Strain FFT peak %.6f Hz, 14X line %.6f Hz.\n', ...
    selectedTimeSec(1), selectedTimeSec(2), strainFftPeakHz, strainOrderFreqHz);
disp(groupsummary(tracTable, 'strain_channel', 'mean', ...
    {'raw_TRAC', 'best_TRAC', 'best_time_shift_s', ...
    'freq_error_vs_strain_fft_hz', 'rmse_mV'}));
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    csvFile, channelCsvFile, channelTimeCsvFile, natureFigFile, channelTimeFigFile);

%% Helpers
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

function strainSources = make_strain_source_options_local(primarySource)
strainDir = fileparts(primarySource.strain_file);
files = dir(fullfile(strainDir, 'AI1-*.mat'));
strainSources = struct('channel', {}, 'file', {}, 'offsetSec', {});
for i = 1:numel(files)
    s = struct();
    s.file = fullfile(files(i).folder, files(i).name);
    s.channel = channel_label_from_file_local(s.file);
    s.offsetSec = primarySource.strain_time_offset_total_sec;
    strainSources(end + 1) = s; %#ok<AGROW>
end
if isempty(strainSources)
    error('No AI1-*.mat strain channel files found in %s.', strainDir);
end
end

function T = compute_channel_sensitivity_table_local(windowResult, methodName, ...
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
        rows(end + 1, :) = {windowResult(w).windowId, tWin(1), tWin(2), ...
            methodLabel, fit.EO, fit.freqHz, fit.amplitudeMm, fit.weightedRmseMv, ...
            strainSources(si).channel, strainSources(si).file, rawPhase, rawTrac, ...
            bestPhase, bestTrac, bestShift, numel(tRaw), rms(yRaw, 'omitnan')}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'window_id', 'time_start_s', ...
    'time_end_s', 'method', 'EO', 'frequency_hz', 'amplitude_mm', ...
    'rmse_mV', 'strain_channel', 'strain_file', 'raw_phase_rad', ...
    'raw_TRAC', 'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm'});
end

function [strainTime, strainMm] = load_strain_file_mm_local(strainFile, offsetSec, strainToMm)
S = load(strainFile, 'Datas');
strainTime = S.Datas(:, 1) + offsetSec;
strainMm = (S.Datas(:, 2) - mean(S.Datas(:, 2), 'omitnan')) * strainToMm;
valid = isfinite(strainTime) & isfinite(strainMm);
strainTime = strainTime(valid);
strainMm = strainMm(valid);
[strainTime, orderIdx] = sort(strainTime);
strainMm = strainMm(orderIdx);
end

function pngFile = plot_step12_nature_validation_local(tracTable, channelTable, ...
    detailData, fftFreqHz, fftAmpMm, strainPeakHz, orderLineHz, figDir, C0, ...
    primaryChannel, selectedTimeSec)
baseName = sprintf('Step12_NatureValidation_%s_%s', C0.dataset, C0.caseTag);
pngFile = fullfile(figDir, [baseName, '.png']);

signalColor = [0.000, 0.447, 0.741];
accentColor = [0.835, 0.369, 0.000];
mutedColor = [0.62, 0.62, 0.62];
darkGray = [0.22, 0.22, 0.22];

fig = figure('Name', 'Step12 20250527 strain validation figure', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 13.2], ...
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
fftMask = fftFreqHz >= 540 & fftFreqHz <= 620;
ylim(ax, [0, max(fftAmpMm(fftMask)) * 1.18]);
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'FFT amplitude (mm)');
title(ax, sprintf('%s strain spectrum, %.1f-%.1f s', ...
    primaryChannel, selectedTimeSec(1), selectedTimeSec(2)));
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
ylim(ax, [min(freqError) - 0.06, max(freqError) + 0.06]);
xlabel(ax, 'Step07J window');
ylabel(ax, '\Deltaf to strain peak (Hz)');
title(ax, 'Identified frequency error');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
plot(ax, tracTable.window_id, tracTable.raw_TRAC, '-o', 'Color', signalColor, ...
    'LineWidth', 1.1, 'MarkerSize', 3.8, 'MarkerFaceColor', 'w');
plot(ax, tracTable.window_id, tracTable.best_TRAC, '-', 'Color', accentColor, ...
    'LineWidth', 1.0);
otherChannels = setdiff(unique(channelTable.strain_channel, 'stable'), {primaryChannel});
if ~isempty(otherChannels)
    maskOther = strcmp(channelTable.strain_channel, otherChannels{1});
    plot(ax, channelTable.window_id(maskOther), channelTable.raw_TRAC(maskOther), ...
        '-', 'Color', mutedColor, 'LineWidth', 0.9);
    text(ax, max(tracTable.window_id) + 0.45, ...
        channelTable.raw_TRAC(find(maskOther, 1, 'last')), otherChannels{1}, ...
        'Color', mutedColor, 'FontName', 'Arial', 'FontSize', 6.5);
end
yline(ax, 0.90, ':', 'Color', mutedColor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(ax, [min(tracTable.window_id) - 0.5, max(tracTable.window_id) + 2.0]);
ylim(ax, [max(0, min([tracTable.raw_TRAC; channelTable.raw_TRAC]) - 0.08), 1.0]);
xlabel(ax, 'Step07J window');
ylabel(ax, 'TRAC');
title(ax, 'Waveform-shape agreement');
text(ax, max(tracTable.window_id) + 0.45, tracTable.raw_TRAC(end), ...
    primaryChannel, 'Color', signalColor, 'FontName', 'Arial', 'FontSize', 6.5);
text(ax, max(tracTable.window_id) + 0.45, min(0.98, tracTable.best_TRAC(end) + 0.025), ...
    'local shift', 'Color', accentColor, 'FontName', 'Arial', 'FontSize', 6.5);
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
    plot(ax, tMs, s.strainBest, '-', 'Color', [0.78, 0.78, 0.78], 'LineWidth', 0.6);
    plot(ax, tMs, yComponent, '-', 'Color', darkGray, 'LineWidth', 0.9);
    plot(ax, tMs, s.reconBest, '-', 'Color', signalColor, 'LineWidth', 1.0);
    xlabel(ax, 'Time from window start (ms)');
    ylabel(ax, 'Displacement (mm)');
    title(ax, sprintf('Representative window W%d, TRAC %.3f', ...
        s.windowId, s.bestTrac));
    text(ax, 0.64, 0.92, ['raw ', primaryChannel], 'Units', 'normalized', ...
        'Color', [0.55, 0.55, 0.55], 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.84, 'synchronous component', 'Units', 'normalized', ...
        'Color', darkGray, 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.76, 'BTT reconstruction', 'Units', 'normalized', ...
        'Color', signalColor, 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
end
panel_label_local(ax, 'd');
apply_nature_axes_local(ax);

export_nature_figure_local(fig, fullfile(figDir, baseName));
end

function [pngFile, T] = plot_channel_time_domain_compare_local(windowResult, ...
    methodName, strainSources, channelTable, figDir, C0, strainToMm, representativeIndex)
if nargin < 8 || ~isfinite(representativeIndex) || representativeIndex < 1 || ...
        representativeIndex > numel(windowResult)
    representativeIndex = max(1, round(numel(windowResult) / 2));
end
wr = windowResult(representativeIndex);
fit = wr.modelFits.(methodName);
tWin = wr.timeWindow;

rows = {};
for si = 1:numel(strainSources)
    [strainTime, strainMm] = load_strain_file_mm_local( ...
        strainSources(si).file, strainSources(si).offsetSec, strainToMm);
    rowMask = strcmp(channelTable.strain_channel, strainSources(si).channel) & ...
        channelTable.window_id == wr.windowId;
    if any(rowMask)
        bestShift = channelTable.best_time_shift_s(find(rowMask, 1, 'first'));
    else
        bestShift = 0;
    end
    [tracVal, ~, tSeg, ySeg, yRecon] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, bestShift, fit.amplitudeMm, fit.freqHz);
    yComponent = project_single_frequency_component_local(tSeg, ySeg, fit.freqHz);
    tMs = 1000 * (tSeg - min(tSeg));
    for ii = 1:numel(tMs)
        rows(end + 1, :) = {wr.windowId, strainSources(si).channel, ...
            strainSources(si).file, bestShift, fit.freqHz, fit.amplitudeMm, ...
            tracVal, tMs(ii), ySeg(ii), yComponent(ii), yRecon(ii)}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'window_id', 'strain_channel', ...
    'strain_file', 'best_time_shift_s', 'frequency_hz', 'amplitude_mm', ...
    'TRAC', 'time_ms', 'raw_mm', 'sync_component_mm', 'btt_reconstruction_mm'});

baseName = sprintf('Step12_ChannelTimeDomainCompare_%s_%s_W%02d', ...
    C0.dataset, C0.caseTag, wr.windowId);
pngFile = fullfile(figDir, [baseName, '.png']);

channels = unique(T.strain_channel, 'stable');
palette = [0.000, 0.447, 0.741; 0.835, 0.369, 0.000; ...
    0.25, 0.25, 0.25; 0.45, 0.45, 0.45];
darkGray = [0.22, 0.22, 0.22];
mutedGray = [0.72, 0.72, 0.72];

fig = figure('Name', 'Step12 channel time-domain comparison', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 13.0], ...
    'NumberTitle', 'off');
set(fig, 'Renderer', 'painters');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
offsetStep = max(0.12, 1.25 * max(abs(T.raw_mm), [], 'omitnan'));
for ci = 1:numel(channels)
    mask = strcmp(T.strain_channel, channels{ci});
    yOffset = (numel(channels) - ci) * offsetStep;
    plot(ax, T.time_ms(mask), T.raw_mm(mask) + yOffset, '-', ...
        'Color', palette(ci, :), 'LineWidth', 0.8);
    text(ax, max(T.time_ms(mask)) * 1.01, yOffset, channels{ci}, ...
        'Color', palette(ci, :), 'FontName', 'Arial', 'FontSize', 6.5, ...
        'VerticalAlignment', 'middle');
end
xlabel(ax, 'Time from window start (ms)');
ylabel(ax, 'Raw displacement + offset (mm)');
title(ax, sprintf('Raw channel waveforms, W%d', wr.windowId));
panel_label_local(ax, 'a');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
for ci = 1:numel(channels)
    mask = strcmp(T.strain_channel, channels{ci});
    plot(ax, T.time_ms(mask), T.sync_component_mm(mask), '-', ...
        'Color', palette(ci, :), 'LineWidth', 1.0);
    text(ax, max(T.time_ms(mask)) * 1.01, ...
        T.sync_component_mm(find(mask, 1, 'last')), channels{ci}, ...
        'Color', palette(ci, :), 'FontName', 'Arial', 'FontSize', 6.5);
end
xlabel(ax, 'Time from window start (ms)');
ylabel(ax, 'Synchronous component (mm)');
title(ax, 'Single-frequency components');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

for ci = 1:min(2, numel(channels))
    ax = nexttile;
    hold(ax, 'on');
    mask = strcmp(T.strain_channel, channels{ci});
    plot(ax, T.time_ms(mask), T.raw_mm(mask), '-', 'Color', mutedGray, 'LineWidth', 0.55);
    plot(ax, T.time_ms(mask), T.sync_component_mm(mask), '-', ...
        'Color', darkGray, 'LineWidth', 0.9);
    plot(ax, T.time_ms(mask), T.btt_reconstruction_mm(mask), '-', ...
        'Color', palette(ci, :), 'LineWidth', 1.0);
    tracVal = T.TRAC(find(mask, 1, 'first'));
    title(ax, sprintf('%s vs BTT, TRAC %.3f', channels{ci}, tracVal));
    xlabel(ax, 'Time from window start (ms)');
    ylabel(ax, 'Displacement (mm)');
    text(ax, 0.64, 0.92, 'raw', 'Units', 'normalized', 'Color', [0.55, 0.55, 0.55], ...
        'FontName', 'Arial', 'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.84, 'sync component', 'Units', 'normalized', 'Color', darkGray, ...
        'FontName', 'Arial', 'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.76, 'BTT reconstruction', 'Units', 'normalized', ...
        'Color', palette(ci, :), 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    panel_label_local(ax, char('b' + ci));
    apply_nature_axes_local(ax);
end

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
design = [sin(2 * pi * freqHz * t), cos(2 * pi * freqHz * t)];
coef = design \ y;
yComponent = design * coef;
end

function label = channel_label_from_file_local(pathText)
[~, name] = fileparts(pathText);
token = regexp(name, 'AI1-\d+', 'match', 'once');
if isempty(token)
    label = name;
else
    label = token;
end
end

function step07jFile = resolve_step07j_result_local(outDir, C0)
candidates = {
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaponly.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaponly.mat', C0.dataset, C0.bladeId, C0.sensorTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag)
    };
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
error('Missing Step07J result for %s under %s.', C0.caseTag, outDir);
end
