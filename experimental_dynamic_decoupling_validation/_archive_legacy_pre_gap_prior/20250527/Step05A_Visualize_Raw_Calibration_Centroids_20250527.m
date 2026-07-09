%% Step05A: Visualize raw calibration time waveforms and centroids
% This diagnostic script reads the archived ARTSCOPE binary files and
% overlays the peak centers / voltage centroids from the staged waveform CSV.
% It is intended for checking the calibration-speed estimate used by Step05.

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
libraryDir = fullfile(rootDir, 'static_gap_waveform_library_20260524_2000Hz');
rawBinaryDir = fullfile(libraryDir, 'raw_binary_source');
csvFile = fullfile(libraryDir, 'raw_source', 'same_blade_waveforms_across_gaps_2000Hz.csv');
outDir = fullfile(thisDir, 'outputs', 'Step05A_single_gap_speed_summary');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sourceNamesToProcess = "all";
sourceNameForWaveformPlot = "ACTS1000_data_1_5mm";
selectedChannel = 1;   % 1-based channel used by the original low-sample visualization script
bladeCount = 6;
rTipMm = 62.0;
edgeFraction = 0.10;
zoomPaddingS = 1.0;

if ~isfile(csvFile)
    error('Missing staged waveform CSV: %s', csvFile);
end
Tcsv = readtable(csvFile, 'TextType', 'string');
if isstring(sourceNamesToProcess) && isscalar(sourceNamesToProcess) && sourceNamesToProcess == "all"
    sourceSummary = unique(Tcsv(:, {'sourceName', 'gapMm'}), 'rows');
    sourceSummary = sortrows(sourceSummary, 'gapMm');
    sourceNamesToProcess = sourceSummary.sourceName;
end

%% 2. Compute centroids for all gaps and visualize one representative waveform
allRows = {};
allCentroidRows = {};
for isrc = 1:numel(sourceNamesToProcess)
    sourceName = sourceNamesToProcess(isrc);
    binPath = fullfile(rawBinaryDir, sourceName + ".bin");
    headerPath = fullfile(rawBinaryDir, sourceName + "_header.txt");
    if ~isfile(binPath) || ~isfile(headerPath)
        warning('Missing raw binary/header for %s. Skipping.', sourceName);
        continue;
    end

    [tRaw, yRaw, fs] = read_acts1000_binary_channel_local(binPath, headerPath, selectedChannel);
    rowsSource = Tcsv(Tcsv.sourceName == sourceName, :);
    if isempty(rowsSource)
        warning('No staged CSV rows for %s. Skipping.', sourceName);
        continue;
    end

    waveKeys = unique(rowsSource(:, {'sourceName', 'gapMm', 'bladeId', 'eventNo', 'peakCenterSample'}), 'rows');
    waveKeys = sortrows(waveKeys, 'peakCenterSample');
    nWave = height(waveKeys);
    centroidSample = nan(nWave, 1);
    centroidTime = nan(nWave, 1);
    peakTime = nan(nWave, 1);
    peakSample = nan(nWave, 1);
    centroidVoltage = nan(nWave, 1);
    peakVoltage = nan(nWave, 1);

    for iw = 1:nWave
        rows = rowsSource.bladeId == waveKeys.bladeId(iw) & ...
            rowsSource.eventNo == waveKeys.eventNo(iw) & ...
            rowsSource.peakCenterSample == waveKeys.peakCenterSample(iw);
        y = rowsSource.voltageMv(rows);
        sampleIndex = rowsSource.sampleIndex(rows);
        edgeCount = max(10, round(edgeFraction * numel(y)));
        baseline = median([y(1:edgeCount); y(end-edgeCount+1:end)], 'omitnan');
        w = max(y - baseline, 0);
        if sum(w, 'omitnan') > eps
            centroidSample(iw) = sum(sampleIndex .* w, 'omitnan') / sum(w, 'omitnan');
            centroidVoltage(iw) = interp1(sampleIndex, y, centroidSample(iw), 'linear', NaN);
        else
            centroidSample(iw) = waveKeys.peakCenterSample(iw);
            centroidVoltage(iw) = NaN;
        end
        peakSample(iw) = waveKeys.peakCenterSample(iw);
        centroidTime(iw) = centroidSample(iw) / fs;
        peakTime(iw) = peakSample(iw) / fs;
        peakVoltage(iw) = interp1(tRaw, yRaw, peakTime(iw), 'linear', NaN);
    end
    centroidVoltageRaw = interp1(tRaw, yRaw, centroidTime, 'linear', NaN);

    if nWave >= bladeCount
        c = centroidSample(1:bladeCount);
        p = peakSample(1:bladeCount);
        bladeIntervalSamplesCentroid = median(diff(c), 'omitnan');
        bladeIntervalSamplesPeak = median(diff(p), 'omitnan');
        rpmCentroid = 60 / (bladeCount * bladeIntervalSamplesCentroid / fs);
        rpmPeak = 60 / (bladeCount * bladeIntervalSamplesPeak / fs);
        tangentialSpeedCentroid = 2 * pi * rTipMm * rpmCentroid / 60;
        tangentialSpeedPeak = 2 * pi * rTipMm * rpmPeak / 60;
    else
        bladeIntervalSamplesCentroid = NaN;
        bladeIntervalSamplesPeak = NaN;
        rpmCentroid = NaN;
        rpmPeak = NaN;
        tangentialSpeedCentroid = NaN;
        tangentialSpeedPeak = NaN;
    end

    oneSummary = table(sourceName, waveKeys.gapMm(1), fs, nWave, ...
        bladeIntervalSamplesCentroid, rpmCentroid, tangentialSpeedCentroid, ...
        bladeIntervalSamplesPeak, rpmPeak, tangentialSpeedPeak, ...
        'VariableNames', {'sourceName', 'recordedGapMm', 'fsHz', 'waveCount', ...
        'bladeIntervalSamplesCentroid', 'rpmCentroid', 'tangentialSpeedCentroidMmS', ...
        'bladeIntervalSamplesPeak', 'rpmPeak', 'tangentialSpeedPeakMmS'});
    allRows{end + 1, 1} = oneSummary; %#ok<SAGROW>

    intervalFromPrevCentroidSamples = [NaN; diff(centroidSample)];
    intervalFromPrevPeakSamples = [NaN; diff(peakSample)];
    oneCentroidTable = table(repmat(sourceName, nWave, 1), ...
        repmat(waveKeys.gapMm(1), nWave, 1), waveKeys.bladeId, waveKeys.eventNo, ...
        peakSample, centroidSample, centroidSample - peakSample, ...
        peakTime, centroidTime, peakVoltage, centroidVoltageRaw, ...
        intervalFromPrevPeakSamples, intervalFromPrevCentroidSamples, ...
        'VariableNames', {'sourceName', 'recordedGapMm', 'bladeId', 'eventNo', ...
        'peakSample', 'centroidSample', 'centroidMinusPeakSamples', ...
        'peakTimeS', 'centroidTimeS', 'peakVoltageMv', 'centroidVoltageMv', ...
        'intervalFromPrevPeakSamples', 'intervalFromPrevCentroidSamples'});
    allCentroidRows{end + 1, 1} = oneCentroidTable; %#ok<SAGROW>

    if sourceName == sourceNameForWaveformPlot
        fig = figure('Name', char(sourceName), 'Color', 'w', ...
            'Units', 'centimeters', 'Position', [2, 2, 24, 16], 'NumberTitle', 'off');
        tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile; hold on; box on;
        plot(tRaw, yRaw, 'k-', 'LineWidth', 0.8, 'DisplayName', 'raw continuous waveform');
        scatter(peakTime, interp1(tRaw, yRaw, peakTime, 'linear', NaN), ...
            36, [0.85, 0.33, 0.10], 'filled', 'DisplayName', 'peak center');
        scatter(centroidTime, centroidVoltageRaw, ...
            46, [0.10, 0.45, 0.80], 'filled', 'DisplayName', 'voltage centroid');
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title(sprintf('%s full continuous waveform, channel 1; centroid rpm %.4f', ...
            sourceName, rpmCentroid), 'Interpreter', 'none');
        legend('Location', 'best');

        nexttile; hold on; box on;
        plot(tRaw, yRaw, 'k-', 'LineWidth', 0.9, 'DisplayName', 'raw continuous waveform');
        scatter(peakTime, peakVoltage, ...
            42, [0.85, 0.33, 0.10], 'filled', 'DisplayName', 'peak center');
        scatter(centroidTime, centroidVoltageRaw, ...
            52, [0.10, 0.45, 0.80], 'filled', 'DisplayName', 'voltage centroid');
        for iw = 1:nWave
            xline(peakTime(iw), '--', 'Color', [0.85, 0.33, 0.10, 0.35], 'HandleVisibility', 'off');
            xline(centroidTime(iw), '-', 'Color', [0.10, 0.45, 0.80, 0.40], 'HandleVisibility', 'off');
            text(centroidTime(iw), centroidVoltageRaw(iw), sprintf(' B%d', waveKeys.bladeId(iw)), ...
                'FontSize', 8, 'Color', [0.10, 0.25, 0.55], 'VerticalAlignment', 'bottom');
        end
        xMin = max(0, min([peakTime; centroidTime]) - zoomPaddingS);
        xMax = min(tRaw(end), max([peakTime; centroidTime]) + zoomPaddingS);
        xlim([xMin, xMax]);
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title('Zoomed segment with peak centers and voltage centroids');
        legend('Location', 'best');

        nexttile; hold on; box on;
        plot(1:nWave, centroidSample, 'o-', 'LineWidth', 1.2, 'DisplayName', 'centroid sample');
        plot(1:nWave, peakSample, 's--', 'LineWidth', 1.0, 'DisplayName', 'peak sample');
        xlabel('Detected waveform order');
        ylabel('Sample index');
        title(sprintf('Median interval: centroid %.3f samples, peak %.3f samples', ...
            bladeIntervalSamplesCentroid, bladeIntervalSamplesPeak));
        legend('Location', 'best');
        set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
        set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);

        figFile = fullfile(outDir, sourceName + "_single_gap_raw_centroid_overlay.png");
        exportgraphics(fig, figFile, 'Resolution', 300);
    end
end

%% 3. Save summary and average calibration speed
if isempty(allRows)
    error('No raw calibration files were visualized.');
end
centroidSummary = vertcat(allRows{:});
summaryFile = fullfile(outDir, 'Step05A_Raw_Centroid_Speed_Summary.csv');
writetable(centroidSummary, summaryFile);
centroidTable = vertcat(allCentroidRows{:});
centroidTableFile = fullfile(outDir, 'Step05A_Raw_Centroid_Table.csv');
writetable(centroidTable, centroidTableFile);

averageRpmCentroid = mean(centroidSummary.rpmCentroid, 'omitnan');
stdRpmCentroid = std(centroidSummary.rpmCentroid, 'omitnan');
averageTangentialSpeedMmS = mean(centroidSummary.tangentialSpeedCentroidMmS, 'omitnan');
stdTangentialSpeedMmS = std(centroidSummary.tangentialSpeedCentroidMmS, 'omitnan');
averageBladeIntervalSamples = mean(centroidSummary.bladeIntervalSamplesCentroid, 'omitnan');

averageSpeedTable = table(selectedChannel, bladeCount, rTipMm, averageBladeIntervalSamples, ...
    averageRpmCentroid, stdRpmCentroid, averageTangentialSpeedMmS, stdTangentialSpeedMmS, ...
    'VariableNames', {'selectedChannel', 'bladeCount', 'rTipMm', ...
    'averageBladeIntervalSamples', 'averageRpmCentroid', 'stdRpmCentroid', ...
    'averageTangentialSpeedMmS', 'stdTangentialSpeedMmS'});
averageSpeedFile = fullfile(outDir, 'Step05A_Average_Calibration_Speed.csv');
writetable(averageSpeedTable, averageSpeedFile);

calibrationSpeed = struct();
calibrationSpeed.selectedChannel = selectedChannel;
calibrationSpeed.bladeCount = bladeCount;
calibrationSpeed.rTipMm = rTipMm;
calibrationSpeed.averageBladeIntervalSamples = averageBladeIntervalSamples;
calibrationSpeed.averageRpmCentroid = averageRpmCentroid;
calibrationSpeed.stdRpmCentroid = stdRpmCentroid;
calibrationSpeed.averageTangentialSpeedMmS = averageTangentialSpeedMmS;
calibrationSpeed.stdTangentialSpeedMmS = stdTangentialSpeedMmS;
calibrationSpeed.sourceSummary = centroidSummary;
save(fullfile(outDir, 'Step05A_Average_Calibration_Speed.mat'), 'calibrationSpeed');

figSpeed = figure('Name', 'Calibration speed summary', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [3, 3, 18, 11], 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on;
plot(centroidSummary.recordedGapMm, centroidSummary.rpmCentroid, 'o-', ...
    'Color', [0.10, 0.45, 0.80], 'MarkerFaceColor', [0.10, 0.45, 0.80], ...
    'LineWidth', 1.2);
yline(averageRpmCentroid, '--', sprintf('mean %.4f rpm', averageRpmCentroid), ...
    'Color', [0.25, 0.25, 0.25], 'LabelHorizontalAlignment', 'left');
xlabel('Recorded gap (mm)');
ylabel('Speed (rpm)');
title('Centroid-based calibration speed');

nexttile; hold on; box on;
plot(centroidSummary.recordedGapMm, centroidSummary.tangentialSpeedCentroidMmS, 'o-', ...
    'Color', [0.10, 0.55, 0.35], 'MarkerFaceColor', [0.10, 0.55, 0.35], ...
    'LineWidth', 1.2);
yline(averageTangentialSpeedMmS, '--', sprintf('mean %.4f mm/s', averageTangentialSpeedMmS), ...
    'Color', [0.25, 0.25, 0.25], 'LabelHorizontalAlignment', 'left');
xlabel('Recorded gap (mm)');
ylabel('Tip speed (mm/s)');
title('Tip circumferential speed for displacement conversion');
set(findall(figSpeed, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(figSpeed, '-property', 'FontSize'), 'FontSize', 9);
exportgraphics(figSpeed, fullfile(outDir, 'Step05A_Calibration_Speed_By_Gap.png'), 'Resolution', 300);

fprintf('\nStep05A complete. Summary saved to:\n  %s\n', summaryFile);
fprintf('Centroid table saved to:\n  %s\n', centroidTableFile);
fprintf('Average calibration speed saved to:\n  %s\n', averageSpeedFile);
fprintf('Mean rpm = %.6f rpm, mean tip speed = %.6f mm/s\n', ...
    averageRpmCentroid, averageTangentialSpeedMmS);
disp(centroidSummary);

%% Local functions
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
