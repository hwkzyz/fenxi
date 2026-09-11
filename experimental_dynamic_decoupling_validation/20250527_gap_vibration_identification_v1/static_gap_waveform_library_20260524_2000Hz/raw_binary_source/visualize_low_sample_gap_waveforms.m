%% Visualize 2000 Hz ACTS1000 clearance waveforms
% Main outputs:
%   1) Blade 1..6: same blade waveform comparison across gaps.
%   2) Blade 1 normalized waveform comparison across gaps.
%   3) Optional selected single-gap raw waveform figures.
%
% The code is intentionally written as a straight script. Only three small
% helper functions are kept at the end for header parsing, label formatting,
% and figure export.

clear;
clc;
close all;

%% Settings
dataFolder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(dataFolder, 'low_sample_gap_visualization');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

selectedChannel = 1;          % 1-based. CH0 in figure labels.
timeRangeSeconds = [];        % [] = whole file; example: [0 3].
minGapToPlotMm = -inf;
maxGapToPlotMm = inf;

makeSelectedSingleGapFigures = false;
singleGapNamesToPlot = "none"; % "none", "all", or ["ACTS1000_data_1_1mm", ...].

bladeCount = 6;


bladeIdsToPlot = 1:6;
normalizedBladeIdToPlot = 1;
firstBladeIdAfterDetection = 1;
sameBladeEventOccurrence = 1; % 1 = first passage of each blade in every gap file.
referenceGapForBladeMatchingMm = NaN; % NaN = first sorted gap.
peakVectorStartEvent = 1;     % Use 6 continuous detected peaks from this event.

prePeakSeconds = 0.50;
postPeakSeconds = 0.50;
detectIgnoreSeconds = 0.50;
expectedBladePeriodSeconds = 2.725;
eventMinDistanceFactor = 0.75;
thresholdSigma = 6;
minProminenceMv = 30;
minWindowPeakAboveBaselineMv = 300;

peakCenterMethod = "sgolay";  % "sgolay", "poly7", or "raw".
sgolayOrder = 3;
sgolayFrameSeconds = 0.025;
polyOrder = 7;

showFigures = true;
savePng = true;
savePdf = false;
saveEmf = false;

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 8.5);
set(groot, 'defaultTextFontSize', 9);
set(groot, 'defaultAxesTickDir', 'in');
set(groot, 'defaultAxesBox', 'on');

%% Find and sort files by clearance
binFiles = dir(fullfile(dataFolder, 'ACTS1000_data*mm.bin'));
if isempty(binFiles)
    error('No ACTS1000_data*mm.bin files found in %s', dataFolder);
end

sourceNames = strings(numel(binFiles), 1);
gapMm = zeros(numel(binFiles), 1);
for k = 1:numel(binFiles)
    [~, sourceNames(k)] = fileparts(binFiles(k).name);
    token = regexp(char(sourceNames(k)), 'ACTS1000_data_(.+?)mm$', 'tokens', 'once');
    if isempty(token)
        error('Cannot parse clearance from file name: %s', sourceNames(k));
    end
    gapMm(k) = str2double(strrep(token{1}, '_', '.'));
end

[gapMm, order] = sort(gapMm);
binFiles = binFiles(order);
sourceNames = sourceNames(order);

keepGap = gapMm >= minGapToPlotMm & gapMm <= maxGapToPlotMm;
binFiles = binFiles(keepGap);
sourceNames = sourceNames(keepGap);
gapMm = gapMm(keepGap);
if isempty(binFiles)
    error('No files remain after applying the selected gap range.');
end

%% Read files, detect peaks, assign blades, and extract aligned raw windows
summaryTables = {};
waveformTables = {};
peakMatchTables = {};
referencePeakVector = [];

for fileNo = 1:numel(binFiles)
    sourceName = sourceNames(fileNo);
    thisGap = gapMm(fileNo);
    binPath = fullfile(binFiles(fileNo).folder, binFiles(fileNo).name);
    headerPath = fullfile(binFiles(fileNo).folder, sourceName + "_header.txt");
    if ~isfile(headerPath)
        error('Missing header file: %s', headerPath);
    end

    fprintf('[%d/%d] Reading %s, gap = %.3f mm\n', ...
        fileNo, numel(binFiles), binFiles(fileNo).name, thisGap);

    headerText = fileread(headerPath);
    chanCount = max(1, round(readHeaderNumber(headerText, 'chanEnableCount', 1)));
    fs = readHeaderNumber(headerText, 'SampleRate', 2000);
    resolution = round(readHeaderNumber(headerText, 'resolution', 12));
    rangeMax = zeros(1, chanCount);
    rangeMin = zeros(1, chanCount);
    for ch = 1:chanCount
        rangeMax(ch) = readHeaderNumber(headerText, sprintf('%drangeMaxValue', ch - 1), 5000);
        rangeMin(ch) = readHeaderNumber(headerText, sprintf('%drangeMinValue', ch - 1), -5000);
    end

    fid = fopen(binPath, 'rb');
    if fid < 0
        error('Cannot open binary file: %s', binPath);
    end
    raw = fread(fid, inf, 'uint16=>uint16');
    fclose(fid);

    validLen = floor(numel(raw) / chanCount) * chanCount;
    raw = reshape(raw(1:validLen), chanCount, []);
    dataMv = zeros(size(raw));
    codeMask = uint16(2^resolution - 1);
    for ch = 1:chanCount
        code = double(bitand(raw(ch, :), codeMask));
        dataMv(ch, :) = (rangeMax(ch) - rangeMin(ch)) / 2^resolution .* code + rangeMin(ch);
    end

    t = (0:size(dataMv, 2) - 1) / fs;
    if ~isempty(timeRangeSeconds)
        keepTime = t >= timeRangeSeconds(1) & t <= timeRangeSeconds(2);
        t = t(keepTime);
        dataMv = dataMv(:, keepTime);
    end

    if selectedChannel > size(dataMv, 1)
        error('selectedChannel=%d exceeds channel count %d for %s.', ...
            selectedChannel, size(dataMv, 1), binFiles(fileNo).name);
    end

    %% Optional single-gap figure
    needSingleGapFigure = makeSelectedSingleGapFigures && ...
        (singleGapNamesToPlot == "all" || any(sourceName == string(singleGapNamesToPlot(:))));
    if needSingleGapFigure
        fig = figure('Color', 'w', 'Visible', figureVisibility(showFigures), ...
            'Units', 'centimeters', 'Position', [2 2 17 8]);
        ax = axes(fig);
        plot(ax, t, dataMv.', 'LineWidth', 0.8);
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Voltage (mV)');
        title(ax, sprintf('%s raw waveform, gap %s', sourceName, gapLabel(thisGap)), ...
            'Interpreter', 'none');
        legend(ax, compose('CH%d', 0:chanCount - 1), 'Location', 'best');
        ax.Box = 'on';
        ax.TickDir = 'in';
        grid(ax, 'off');
        exportFigure(fig, outputFolder, sourceName + "_raw_waveform_2000Hz", ...
            savePng, savePdf, saveEmf);
    end

    %% Summary rows
    channel = (0:chanCount - 1).';
    summaryTables{end + 1, 1} = table( ...
        repmat(sourceName, chanCount, 1), repmat(thisGap, chanCount, 1), channel, ...
        repmat(fs, chanCount, 1), repmat(size(dataMv, 2), chanCount, 1), ...
        repmat(t(end) - t(1), chanCount, 1), mean(dataMv, 2), std(dataMv, 0, 2), ...
        min(dataMv, [], 2), max(dataMv, [], 2), max(dataMv, [], 2) - min(dataMv, [], 2), ...
        'VariableNames', {'sourceName', 'gapMm', 'channel', 'fsHz', ...
        'sampleCount', 'durationSeconds', 'meanMv', 'stdMv', ...
        'minMv', 'maxMv', 'peakToPeakMv'});

    %% Peak detection on selected channel
    yAll = dataMv(selectedChannel, :);
    baselineAll = median(yAll, 'omitnan');
    eventSignal = yAll(:) - baselineAll;
    noiseSigma = 1.4826 * mad(eventSignal, 1);
    minProminence = max(thresholdSigma * noiseSigma, minProminenceMv);
    minPeakDistance = max(1, round(eventMinDistanceFactor * expectedBladePeriodSeconds * fs));
    [peakAmplitudes, peakSamples] = findpeaks(eventSignal, ...
        'MinPeakDistance', minPeakDistance, ...
        'MinPeakProminence', minProminence);
    keepPeaks = peakSamples >= round(detectIgnoreSeconds * fs);
    peakSamples = peakSamples(keepPeaks);
    peakAmplitudes = peakAmplitudes(keepPeaks);

    if isempty(peakSamples)
        warning('No peaks detected in %s.', sourceName);
        continue;
    end

    %% Build a six-peak vector and match its cyclic order to the reference gap
    if numel(peakSamples) < bladeCount
        warning('Less than %d peaks detected in %s.', bladeCount, sourceName);
        continue;
    end

    vectorStart = min(max(1, peakVectorStartEvent), numel(peakSamples) - bladeCount + 1);
    vectorEvents = vectorStart:(vectorStart + bladeCount - 1);
    peakVector = peakAmplitudes(vectorEvents).';
    peakVector = peakVector - mean(peakVector, 'omitnan');
    peakNorm = norm(peakVector);
    if peakNorm > 0
        peakVector = peakVector / peakNorm;
    end

    useAsReference = isempty(referencePeakVector);
    if ~isnan(referenceGapForBladeMatchingMm)
        useAsReference = abs(thisGap - referenceGapForBladeMatchingMm) < 1e-9;
    end
    if useAsReference
        referencePeakVector = peakVector;
        bestShift = 0;
        bestScore = 1;
    else
        scores = zeros(bladeCount, 1);
        for shift = 0:(bladeCount - 1)
            scores(shift + 1) = referencePeakVector(:).' * circshift(peakVector(:), -shift);
        end
        [bestScore, bestIdx] = max(scores);
        bestShift = bestIdx - 1;
    end

    peakMatchTables{end + 1, 1} = table(sourceName, thisGap, ...
        vectorStart, bestShift, bestScore, peakVector(1), peakVector(2), ...
        peakVector(3), peakVector(4), peakVector(5), peakVector(6), ...
        'VariableNames', {'sourceName', 'gapMm', 'vectorStartEvent', ...
        'cyclicShift', 'matchScore', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6'});

    eventNo = (1:numel(peakSamples)).';
    bladeIdOfEvent = mod(eventNo - 1 - bestShift + ...
        firstBladeIdAfterDetection - 1, bladeCount) + 1;

    %% Extract one aligned raw waveform for every blade
    for bladeId = bladeIdsToPlot(:).'
        selectedEvents = find(bladeIdOfEvent == bladeId);
        if numel(selectedEvents) < sameBladeEventOccurrence
            warning('Only %d events for blade %d in %s.', ...
                numel(selectedEvents), bladeId, sourceName);
            continue;
        end

        selectedEventIndex = selectedEvents(sameBladeEventOccurrence);
        roughPeakSample = peakSamples(selectedEventIndex);

        % Refine peak center with raw, SG-smoothed, or polynomial-fitted curve.
        halfRefineSamples = max(2, round(0.10 * fs));
        idxFit = (max(1, roughPeakSample - halfRefineSamples): ...
            min(numel(yAll), roughPeakSample + halfRefineSamples)).';
        yFitRaw = yAll(idxFit).';
        switch lower(string(peakCenterMethod))
            case "raw"
                yFit = yFitRaw;
            case "sgolay"
                frameSamples = max(3, round(sgolayFrameSeconds * fs));
                if mod(frameSamples, 2) == 0
                    frameSamples = frameSamples + 1;
                end
                frameSamples = min(frameSamples, numel(yFitRaw) - mod(numel(yFitRaw) + 1, 2));
                frameSamples = max(frameSamples, sgolayOrder + 2 + mod(sgolayOrder + 2, 2));
                yFit = sgolayfilt(yFitRaw, sgolayOrder, frameSamples);
            case "poly7"
                xFit = linspace(-1, 1, numel(yFitRaw)).';
                orderFit = min(polyOrder, numel(yFitRaw) - 1);
                pFit = polyfit(xFit, yFitRaw, orderFit);
                yFit = polyval(pFit, xFit);
            otherwise
                error('Unknown peakCenterMethod: %s', peakCenterMethod);
        end
        [~, localCenter] = max(yFit);
        peakCenterSample = idxFit(localCenter);
        peakCenterOffsetSamples = peakCenterSample - roughPeakSample;

        idxWin = (peakCenterSample - round(prePeakSeconds * fs)): ...
            (peakCenterSample + round(postPeakSeconds * fs));
        idxWin = idxWin(idxWin >= 1 & idxWin <= numel(yAll));
        timeRelative = (idxWin(:) - peakCenterSample) / fs;
        yWin = yAll(idxWin).';

        edgeCount = max(3, round(0.10 * numel(yWin)));
        baselineWindow = median([yWin(1:edgeCount); yWin(end-edgeCount+1:end)], 'omitnan');
        peakAboveBaseline = max(yWin) - baselineWindow;
        if peakAboveBaseline < minWindowPeakAboveBaselineMv
            warning('Skip %s blade %d event %d: peak above baseline is %.1f mV.', ...
                sourceName, bladeId, eventNo(selectedEventIndex), peakAboveBaseline);
            continue;
        end

        n = numel(yWin);
        waveformTables{end + 1, 1} = table( ...
            repmat(sourceName, n, 1), repmat(thisGap, n, 1), repmat(bladeId, n, 1), ...
            repmat(eventNo(selectedEventIndex), n, 1), repmat(peakCenterSample, n, 1), ...
            repmat(peakCenterOffsetSamples, n, 1), repmat(peakAboveBaseline, n, 1), ...
            idxWin(:), timeRelative(:), yWin(:), ...
            'VariableNames', {'sourceName', 'gapMm', 'bladeId', 'eventNo', ...
            'peakCenterSample', 'peakCenterOffsetSamples', 'peakAboveBaselineMv', ...
            'sampleIndex', 'timeRelative_s', 'voltageMv'});
    end
end

%% Save data tables
summaryTable = vertcat(summaryTables{:});
writetable(summaryTable, fullfile(outputFolder, 'low_sample_gap_summary.csv'));
if ~isempty(peakMatchTables)
    peakMatchTable = vertcat(peakMatchTables{:});
    writetable(peakMatchTable, fullfile(outputFolder, 'peak_vector_matching_summary.csv'));
end

if isempty(waveformTables)
    error('No valid blade waveforms were extracted.');
end
sameBladeTable = vertcat(waveformTables{:});
writetable(sameBladeTable, fullfile(outputFolder, 'same_blade_waveforms_across_gaps_2000Hz.csv'));

%% Plot each blade across gaps
for bladeId = bladeIdsToPlot
    rowsBlade = sameBladeTable.bladeId == bladeId;
    if ~any(rowsBlade)
        warning('No waveform rows for blade %d.', bladeId);
        continue;
    end

    Tb = sameBladeTable(rowsBlade, :);
    gaps = unique(Tb.gapMm, 'stable');
    colors = turbo(numel(gaps));

    fig = figure('Color', 'w', 'Visible', figureVisibility(showFigures), ...
        'Units', 'centimeters', 'Position', [2 2 17 9.5]);
    ax = axes(fig);
    hold(ax, 'on');
    for k = 1:numel(gaps)
        rows = Tb.gapMm == gaps(k);
        plot(ax, Tb.timeRelative_s(rows), Tb.voltageMv(rows), ...
            'LineWidth', 0.85, 'Color', colors(k, :), ...
            'DisplayName', gapLabel(gaps(k)));
    end
    xline(ax, 0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Voltage (mV)');
    title(ax, sprintf('Blade %d raw waveforms at different clearances, CH%d, 2000 Hz', ...
        bladeId, selectedChannel - 1));
    legend(ax, 'Location', 'eastoutside', 'Box', 'off');
    ax.Box = 'on';
    ax.TickDir = 'in';
    grid(ax, 'off');
    exportFigure(fig, outputFolder, ...
        sprintf('blade_%d_raw_waveforms_across_gaps_2000Hz', bladeId), ...
        savePng, savePdf, saveEmf);
end

%% Plot normalized comparison for one blade
rowsBlade = sameBladeTable.bladeId == normalizedBladeIdToPlot;
if any(rowsBlade)
    Tb = sameBladeTable(rowsBlade, :);
    gaps = unique(Tb.gapMm, 'stable');
    colors = turbo(numel(gaps));

    fig = figure('Color', 'w', 'Visible', figureVisibility(showFigures), ...
        'Units', 'centimeters', 'Position', [2 2 17 9.5]);
    ax = axes(fig);
    hold(ax, 'on');
    for k = 1:numel(gaps)
        rows = Tb.gapMm == gaps(k);
        tt = Tb.timeRelative_s(rows);
        yy = Tb.voltageMv(rows);
        edgeCount = max(3, round(0.10 * numel(yy)));
        baseline = median([yy(1:edgeCount); yy(end-edgeCount+1:end)], 'omitnan');
        amplitude = max(yy) - baseline;
        if amplitude <= 0 || isnan(amplitude)
            warning('Skip normalized curve %.3f mm: invalid amplitude.', gaps(k));
            continue;
        end
        plot(ax, tt, (yy - baseline) / amplitude, ...
            'LineWidth', 0.9, 'Color', colors(k, :), ...
            'DisplayName', gapLabel(gaps(k)));
    end
    xline(ax, 0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Normalized voltage');
    title(ax, sprintf('Blade %d normalized waveforms at different clearances, CH%d, 2000 Hz', ...
        normalizedBladeIdToPlot, selectedChannel - 1));
    ylim(ax, [-0.1, 1.1]);
    legend(ax, 'Location', 'eastoutside', 'Box', 'off');
    ax.Box = 'on';
    ax.TickDir = 'in';
    grid(ax, 'off');
    exportFigure(fig, outputFolder, ...
        sprintf('blade_%d_normalized_waveforms_across_gaps_2000Hz', normalizedBladeIdToPlot), ...
        savePng, savePdf, saveEmf);
end

fprintf('\nDone.\n');
fprintf('Figures and summary saved to:\n%s\n', outputFolder);

%% Small helper functions
function value = readHeaderNumber(txt, key, defaultValue)
    token = regexp(txt, [regexptranslate('escape', key), ...
        '\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'], 'tokens', 'once');
    if isempty(token)
        value = defaultValue;
    else
        value = str2double(token{1});
    end
end

function label = gapLabel(gapMm)
    text = regexprep(sprintf('%.2f', gapMm), '0+$', '');
    text = regexprep(text, '\.$', '');
    label = sprintf('%s mm', text);
end

function visibility = figureVisibility(showFigures)
    if showFigures
        visibility = 'on';
    else
        visibility = 'off';
    end
end

function exportFigure(fig, outputFolder, baseName, savePng, savePdf, saveEmf)
    if savePng
        exportgraphics(fig, fullfile(outputFolder, baseName + ".png"), 'Resolution', 300);
    end
    if savePdf
        exportgraphics(fig, fullfile(outputFolder, baseName + ".pdf"), 'ContentType', 'vector');
    end
    if saveEmf
        print(fig, fullfile(outputFolder, baseName + ".emf"), '-dmeta', '-r300');
    end
    drawnow;
end
