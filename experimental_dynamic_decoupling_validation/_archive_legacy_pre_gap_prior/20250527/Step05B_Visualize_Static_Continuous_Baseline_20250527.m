%% Step05B: Visualize static calibration continuous-waveform baseline
% This diagnostic checks whether the static-gap baseline used in Step05 is
% taken from stable no-pulse portions of the original continuous waveform.
%
% It does not change Step05 outputs.  Use the figures and CSV summary to
% decide whether Step05 should use full-waveform p05, p10, median, or a
% pulse-masked background baseline.

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
libraryDir = fullfile(rootDir, 'static_gap_waveform_library_20260524_2000Hz');
rawBinaryDir = fullfile(libraryDir, 'raw_binary_source');
manifestFile = fullfile(libraryDir, 'gap_waveform_manifest.csv');
outDir = fullfile(thisDir, 'outputs', 'Step05B_static_continuous_baseline');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

selectedChannel = 1;  % MATLAB 1-based; original plots label this as CH0
sourceNamesForDetailedPlot = "ACTS1000_data_1_5mm";
plotAllDetailed = false;
maxPlotPoints = 30000;
zoomDurationS = 20;
pulseThresholdFraction = 0.12;
pulsePadS = 0.35;
backgroundSegmentTrimFraction = 0.12;  % Trim each no-pulse segment inward from both edges
minBackgroundSegmentS = 0.08;
maxBackgroundPlotPoints = 12000;
baselineRevCount = 2;        % Use only the first 1-2 stable revolutions for baseline
bladeCount = 6;
baselineWindowPadS = 0.60;

%% 2. Discover static-gap raw sources
manifest = readtable(manifestFile, 'TextType', 'string', 'Delimiter', ',', ...
    'NumHeaderLines', 0, 'ReadVariableNames', true);
if width(manifest) == 6 && strcmp(manifest.Properties.VariableNames{1}, 'Var1')
    manifest.Properties.VariableNames = {'recordedGapMm', 'trueReferenceGapMm', ...
        'rowCount', 'bladeIds', 'relativeCsv', 'ExtraVar'};
elseif width(manifest) == 5 && strcmp(manifest.Properties.VariableNames{1}, 'Var1')
    manifest.Properties.VariableNames = {'recordedGapMm', 'trueReferenceGapMm', ...
        'rowCount', 'bladeIds', 'relativeCsv'};
end
manifest = sortrows(manifest, 'trueReferenceGapMm');

summaryRows = cell(height(manifest), 1);
detailCache = struct();
for ig = 1:height(manifest)
    T = readtable(fullfile(libraryDir, manifest.relativeCsv(ig)), 'TextType', 'string');
    sourceName = string(T.sourceName(find(T.sourceName ~= "", 1, 'first')));
    binPath = fullfile(rawBinaryDir, sourceName + ".bin");
    headerPath = fullfile(rawBinaryDir, sourceName + "_header.txt");
    [t, yMv, fs] = read_acts1000_binary_channel_local(binPath, headerPath, selectedChannel);

    [baselineWindowMask, baselineWindowInfo] = select_initial_revolution_window_local( ...
        yMv, fs, bladeCount, baselineRevCount, baselineWindowPadS);
    diag = compute_static_background_baselines_local(yMv, fs, baselineWindowMask, ...
        pulseThresholdFraction, pulsePadS, backgroundSegmentTrimFraction, ...
        minBackgroundSegmentS);

    summaryRows{ig} = table(sourceName, manifest.recordedGapMm(ig), ...
        manifest.trueReferenceGapMm(ig), fs, numel(yMv), ...
        baselineWindowInfo.startTimeS, baselineWindowInfo.endTimeS, ...
        baselineWindowInfo.peakCount, diag.fullP05Mv, diag.fullP10Mv, diag.fullMedianMv, ...
        diag.backgroundMeanMv, diag.backgroundMedianMv, diag.backgroundP05Mv, ...
        diag.backgroundStdMv, diag.backgroundFraction, diag.thresholdMv, ...
        'VariableNames', {'sourceName', 'recordedGapMm', 'trueReferenceGapMm', ...
        'fsHz', 'sampleCount', 'baselineWindowStartS', 'baselineWindowEndS', ...
        'baselineWindowPeakCount', 'fullP05Mv', 'fullP10Mv', 'fullMedianMv', ...
        'backgroundMeanMv', 'backgroundMedianMv', 'backgroundP05Mv', ...
        'backgroundStdMv', 'backgroundFraction', 'pulseThresholdMv'});

    makeDetail = plotAllDetailed || any(sourceName == string(sourceNamesForDetailedPlot(:)));
    if makeDetail
        detailCache.(matlab.lang.makeValidName(sourceName)) = struct( ...
            'sourceName', sourceName, 'recordedGapMm', manifest.recordedGapMm(ig), ...
            'trueReferenceGapMm', manifest.trueReferenceGapMm(ig), ...
            't', t, 'yMv', yMv, 'fs', fs, 'diag', diag, ...
            'baselineWindowInfo', baselineWindowInfo);
    end
end

baselineSummary = vertcat(summaryRows{:});
summaryFile = fullfile(outDir, 'Step05B_Static_Continuous_Baseline_Summary.csv');
writetable(baselineSummary, summaryFile);

%% 3. Detailed waveform figures for selected gaps
detailNames = fieldnames(detailCache);
for id = 1:numel(detailNames)
    D = detailCache.(detailNames{id});
    fig = figure('Name', char(D.sourceName), 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 26, 17], 'NumberTitle', 'off');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    idxPlot = decimate_indices_local(numel(D.t), maxPlotPoints);
    bgIdxPlot = idxPlot(D.diag.backgroundMask(idxPlot));
    pulseIdxPlot = idxPlot(~D.diag.backgroundMask(idxPlot));
    bgIdxPlot = thin_indices_local(bgIdxPlot, maxBackgroundPlotPoints);

    nexttile; hold on; box on;
    plot(D.t(idxPlot), D.yMv(idxPlot), '-', 'Color', [0.15, 0.15, 0.15], ...
        'LineWidth', 0.6, 'DisplayName', 'raw continuous waveform');
    scatter(D.t(bgIdxPlot), D.yMv(bgIdxPlot), 4, [0.10, 0.45, 0.80], ...
        'filled', 'DisplayName', 'auto no-pulse samples');
    scatter(D.t(pulseIdxPlot), D.yMv(pulseIdxPlot), 4, [0.85, 0.33, 0.10], ...
        'filled', 'DisplayName', 'excluded pulse vicinity');
    xline(D.baselineWindowInfo.startTimeS, '--', 'baseline window start', ...
        'Color', [0.30, 0.30, 0.30], 'HandleVisibility', 'off');
    xline(D.baselineWindowInfo.endTimeS, '--', 'baseline window end', ...
        'Color', [0.30, 0.30, 0.30], 'HandleVisibility', 'off');
    add_baseline_lines_local(D.diag);
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    title(sprintf('%s, recorded gap %.2f mm', D.sourceName, D.recordedGapMm), ...
        'Interpreter', 'none');
    legend('Location', 'best');

    nexttile; hold on; box on;
    zoomMask = D.t >= D.baselineWindowInfo.startTimeS & ...
        D.t <= min(max(D.t), D.baselineWindowInfo.startTimeS + zoomDurationS);
    plot(D.t(zoomMask), D.yMv(zoomMask), '-', 'Color', [0.15, 0.15, 0.15], ...
        'LineWidth', 0.8, 'DisplayName', 'raw');
    bgZoomIdx = find(zoomMask & D.diag.backgroundMask);
    bgZoomIdx = thin_indices_local(bgZoomIdx, maxBackgroundPlotPoints);
    scatter(D.t(bgZoomIdx), D.yMv(bgZoomIdx), 8, [0.10, 0.45, 0.80], ...
        'filled', 'DisplayName', 'no-pulse');
    add_baseline_lines_local(D.diag);
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    title(sprintf('Baseline window: first %d revolutions, trimmed no-pulse selection', ...
        baselineRevCount));
    legend('Location', 'best');

    nexttile; hold on; box on;
    histogram(D.yMv, 140, 'FaceColor', [0.75, 0.75, 0.75], ...
        'EdgeColor', 'none', 'DisplayName', 'all samples');
    histogram(D.yMv(D.diag.backgroundMask), 80, 'FaceColor', [0.10, 0.45, 0.80], ...
        'FaceAlpha', 0.65, 'EdgeColor', 'none', 'DisplayName', 'no-pulse samples');
    xline(D.diag.fullP05Mv, '--', 'p05 all', 'Color', [0.65, 0.25, 0.20]);
    xline(D.diag.backgroundMedianMv, '-', 'median no-pulse', 'Color', [0.10, 0.35, 0.75], ...
        'LineWidth', 1.2);
    xlabel('Voltage (mV)');
    ylabel('Count');
    title('Voltage distribution');
    legend('Location', 'best');

    nexttile; hold on; box on;
    plot(baselineSummary.recordedGapMm, baselineSummary.fullP05Mv, 'o-', ...
        'LineWidth', 1.0, 'DisplayName', 'full p05');
    plot(baselineSummary.recordedGapMm, baselineSummary.fullP10Mv, 's-', ...
        'LineWidth', 1.0, 'DisplayName', 'full p10');
    plot(baselineSummary.recordedGapMm, baselineSummary.backgroundMedianMv, '^-', ...
        'LineWidth', 1.0, 'DisplayName', 'no-pulse median');
    xline(D.recordedGapMm, '--k', 'current gap', 'LabelOrientation', 'horizontal');
    xlabel('Recorded gap (mm)');
    ylabel('Baseline candidate (mV)');
    title('Baseline candidates across static gaps');
    legend('Location', 'best');

    set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);
    figFile = fullfile(outDir, D.sourceName + "_static_baseline_diagnostic.png");
    exportgraphics(fig, figFile, 'Resolution', 300);
end

fprintf('\nStep05B complete.\n');
fprintf('Baseline summary saved to:\n  %s\n', summaryFile);
disp(baselineSummary(:, {'sourceName', 'recordedGapMm', 'fullP05Mv', ...
    'fullP10Mv', 'fullMedianMv', 'backgroundMedianMv', 'backgroundFraction'}));

%% Local functions
function [windowMask, info] = select_initial_revolution_window_local(yMv, fs, bladeCount, revCount, padS)
yMv = yMv(:);
n = numel(yMv);
t = (0:n-1).' / fs;
baseline0 = prctile(yMv, 10);
peakHeight = prctile(yMv, 99) - baseline0;
minPeakHeight = baseline0 + 0.35 * max(peakHeight, eps);
minPeakDistance = max(1, round(0.7 * median_blade_interval_guess_local(yMv, fs)));
[~, locs] = findpeaks(yMv, 'MinPeakHeight', minPeakHeight, ...
    'MinPeakDistance', minPeakDistance);
if numel(locs) < bladeCount
    warning('Could not find enough peaks for initial revolution window; using first 20 s.');
    startIdx = 1;
    endIdx = min(n, round(20 * fs));
    peakCount = numel(locs);
else
    requiredPeaks = min(numel(locs), revCount * bladeCount);
    useLocs = locs(1:requiredPeaks);
    startIdx = max(1, useLocs(1) - round(padS * fs));
    endIdx = min(n, useLocs(end) + round(padS * fs));
    peakCount = requiredPeaks;
end
windowMask = false(n, 1);
windowMask(startIdx:endIdx) = true;
info = struct();
info.startTimeS = t(startIdx);
info.endTimeS = t(endIdx);
info.startIdx = startIdx;
info.endIdx = endIdx;
info.peakCount = peakCount;
end

function intervalGuess = median_blade_interval_guess_local(yMv, fs)
yMv = yMv(:);
baseline0 = prctile(yMv, 10);
peakHeight = prctile(yMv, 99) - baseline0;
[~, locs] = findpeaks(yMv, 'MinPeakHeight', baseline0 + 0.35 * max(peakHeight, eps), ...
    'MinPeakDistance', round(0.5 * fs));
if numel(locs) >= 2
    intervalGuess = median(diff(locs), 'omitnan');
else
    intervalGuess = 2.7 * fs;
end
end

function diag = compute_static_background_baselines_local(yMv, fs, windowMask, pulseThresholdFraction, pulsePadS, ...
    backgroundSegmentTrimFraction, minBackgroundSegmentS)
yMv = yMv(:);
windowMask = windowMask(:) & isfinite(yMv);
yWin = yMv(windowMask);
fullP05Mv = prctile(yWin, 5);
fullP10Mv = prctile(yWin, 10);
fullMedianMv = median(yWin, 'omitnan');
upperMv = prctile(yWin, 99.5);
thresholdMv = fullMedianMv + pulseThresholdFraction * max(upperMv - fullMedianMv, eps);
pulseCore = yMv > thresholdMv & windowMask;
padN = max(1, round(pulsePadS * fs));
kernel = ones(2 * padN + 1, 1);
pulseVicinity = conv(double(pulseCore), kernel, 'same') > 0;
backgroundMask0 = windowMask & ~pulseVicinity & isfinite(yMv);
backgroundMask = trim_background_segments_local(backgroundMask0, fs, ...
    backgroundSegmentTrimFraction, minBackgroundSegmentS);
if nnz(backgroundMask) < 0.10 * nnz(backgroundMask0)
    warning('Trimmed background mask is too small; falling back to the untrimmed no-pulse mask.');
    backgroundMask = backgroundMask0;
end
bg = yMv(backgroundMask);
diag = struct();
diag.fullP05Mv = fullP05Mv;
diag.fullP10Mv = fullP10Mv;
diag.fullMedianMv = fullMedianMv;
diag.backgroundMeanMv = mean(bg, 'omitnan');
diag.backgroundMedianMv = median(bg, 'omitnan');
diag.backgroundP05Mv = prctile(bg, 5);
diag.backgroundStdMv = std(bg, 'omitnan');
diag.backgroundFraction = nnz(backgroundMask) / nnz(windowMask);
diag.thresholdMv = thresholdMv;
diag.backgroundSegmentTrimFraction = backgroundSegmentTrimFraction;
diag.minBackgroundSegmentS = minBackgroundSegmentS;
diag.backgroundMask = backgroundMask;
diag.untrimmedBackgroundMask = backgroundMask0;
diag.windowMask = windowMask;
end

function trimmedMask = trim_background_segments_local(mask, fs, trimFraction, minSegmentS)
mask = mask(:) & isfinite(mask(:));
trimmedMask = false(size(mask));
idx = find(mask);
if isempty(idx)
    return;
end
breaks = [1; find(diff(idx) > 1) + 1; numel(idx) + 1];
minSegmentN = max(3, round(minSegmentS * fs));
for iseg = 1:numel(breaks)-1
    segIdx = idx(breaks(iseg):breaks(iseg+1)-1);
    n = numel(segIdx);
    if n < minSegmentN
        continue;
    end
    trimN = floor(trimFraction * n);
    if 2 * trimN >= n - 2
        trimN = max(0, floor((n - 2) / 2));
    end
    keepIdx = segIdx((1 + trimN):(n - trimN));
    trimmedMask(keepIdx) = true;
end
end

function add_baseline_lines_local(diag)
yline(diag.fullP05Mv, '--', 'p05 all', 'Color', [0.65, 0.25, 0.20], ...
    'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
yline(diag.fullP10Mv, ':', 'p10 all', 'Color', [0.55, 0.35, 0.20], ...
    'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
yline(diag.backgroundMedianMv, '-', 'median no-pulse', ...
    'Color', [0.10, 0.35, 0.75], 'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left');
end

function idx = decimate_indices_local(n, maxPoints)
if n <= maxPoints
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxPoints))).';
end
end

function idxThin = thin_indices_local(idx, maxPoints)
idx = idx(:);
if numel(idx) <= maxPoints
    idxThin = idx;
else
    pick = unique(round(linspace(1, numel(idx), maxPoints)));
    idxThin = idx(pick);
end
end

function [t, y, fs] = read_acts1000_binary_channel_local(binPath, headerPath, selectedChannel)
headerText = fileread(headerPath);
chanCount = max(1, round(read_header_number_local(headerText, 'chanEnableCount', 1)));
fs = read_header_number_local(headerText, 'SampleRate', 2000);
resolution = round(read_header_number_local(headerText, 'resolution', 12));
if selectedChannel > chanCount
    error('selectedChannel=%d exceeds channel count %d for %s.', selectedChannel, chanCount, binPath);
end

rangeMax = zeros(1, chanCount);
rangeMin = zeros(1, chanCount);
for ch = 1:chanCount
    rangeMax(ch) = read_header_number_local(headerText, sprintf('%drangeMaxValue', ch - 1), 5000);
    rangeMin(ch) = read_header_number_local(headerText, sprintf('%drangeMinValue', ch - 1), -5000);
end

fid = fopen(binPath, 'rb');
if fid < 0
    error('Cannot open binary file: %s', binPath);
end
raw = fread(fid, inf, 'uint16=>uint16');
fclose(fid);

validLen = floor(numel(raw) / chanCount) * chanCount;
raw = reshape(raw(1:validLen), chanCount, []);
codeMask = uint16(2^resolution - 1);
code = double(bitand(raw(selectedChannel, :), codeMask));
y = (rangeMax(selectedChannel) - rangeMin(selectedChannel)) / 2^resolution .* code + ...
    rangeMin(selectedChannel);
t = (0:numel(y)-1).' / fs;
y = y(:);
end

function value = read_header_number_local(headerText, key, defaultValue)
expr = [regexptranslate('escape', key), '\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'];
tok = regexp(headerText, expr, 'tokens', 'once');
if isempty(tok)
    value = defaultValue;
else
    value = str2double(tok{1});
end
end
