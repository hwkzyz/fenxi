%% Step05: Build static gap response surface for 20250527 dynamic validation
% This step reads the staged multi-gap waveform library and fits the
% response-surface model used by the waveform-domain decoupling method:
%
%   F(g,x) = B0(x) + B1(x)/g + B2(x)*log(g/g0)
%
% The source library is copied under:
%   ../static_gap_waveform_library_20260524_2000Hz/
%
% Normal downstream use should load:
%   outputs/Step05_Response_Surface_20250527.mat

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20250527();
libraryDir = packageCfg.paths.staticGapLibrary;
rawBinaryDir = fullfile(libraryDir, 'raw_binary_source');
manifestFile = fullfile(libraryDir, 'gap_waveform_manifest.csv');
outDir = packageCfg.paths.results;
speedSummaryFile = fullfile(outDir, 'Step05A_single_gap_speed_summary', ...
    'Step05A_Average_Calibration_Speed.csv');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

targetBladeMode = 'auto_reference';  % auto_reference, explicit_reference, or average_all_blades
explicitReferenceBladeId = 1;
targetBladeModeOverride = strtrim(getenv('STEP05_TARGET_BLADE_MODE'));
if ~isempty(targetBladeModeOverride)
    targetBladeMode = targetBladeModeOverride;
end
explicitBladeOverride = str2double(strtrim(getenv('STEP05_REFERENCE_BLADE_ID')));
if isfinite(explicitBladeOverride)
    explicitReferenceBladeId = explicitBladeOverride;
end
gapZeroOffsetMm = 0.2;
g0Mm = 1.0;
fitWindowHalfWidthS = 0.20;
pointsPerWaveform = 801;
effectiveWindowThreshold = 0.12;
edgeFraction = 0.10;
staticBaselineMethod = 'initial_2rev_trimmed_background_median';
staticBaselineRevCount = 2;
staticBaselineWindowPadS = 0.60;
staticBaselinePulseThresholdFraction = 0.12;
staticBaselineSegmentTrimFraction = 0.12;
staticBaselineMinSegmentS = 0.08;
minValidGapsForFit = 4;
bladeCount = 6;
minBladesPerAverage = 3;
sampleRateHz = 2000;
rTipMm = 62.0;
sourceChannel = 1;  % MATLAB 1-based channel index used to generate the staged calibration CSV
saveFigures = true;

if ~isfile(manifestFile)
    error('Missing staged gap waveform manifest: %s', manifestFile);
end

%% 2. Load staged multi-gap waveforms and estimate calibration speed
manifest = readtable(manifestFile, 'TextType', 'string', 'Delimiter', ',', ...
    'NumHeaderLines', 0, 'ReadVariableNames', true);
if width(manifest) == 6 && strcmp(manifest.Properties.VariableNames{1}, 'Var1')
    manifest.Properties.VariableNames = {'recordedGapMm', 'trueReferenceGapMm', ...
        'rowCount', 'bladeIds', 'relativeCsv', 'ExtraVar'};
elseif width(manifest) == 5 && strcmp(manifest.Properties.VariableNames{1}, 'Var1')
    manifest.Properties.VariableNames = {'recordedGapMm', 'trueReferenceGapMm', ...
        'rowCount', 'bladeIds', 'relativeCsv'};
end
if ~ismember('recordedGapMm', manifest.Properties.VariableNames)
    error('Manifest column recordedGapMm was not found. Columns are: %s', ...
        strjoin(manifest.Properties.VariableNames, ', '));
end
recordedGaps = manifest.recordedGapMm(:);
trueGaps = manifest.trueReferenceGapMm(:);
[trueGaps, orderGap] = sort(trueGaps);
recordedGaps = recordedGaps(orderGap);
manifest = manifest(orderGap, :);

allBladeIds = [];
rawRows = cell(height(manifest), 1);
centroidSamplesByGap = cell(height(manifest), 1);
bladeIntervalSamples = nan(height(manifest), 1);
revPeriodS = nan(height(manifest), 1);
rotFreqHz = nan(height(manifest), 1);
rpm = nan(height(manifest), 1);
tangentialSpeedMmS = nan(height(manifest), 1);
staticBaselineMvByGap = nan(height(manifest), 1);
staticBaselineWindowStartS = nan(height(manifest), 1);
staticBaselineWindowEndS = nan(height(manifest), 1);
staticBaselineWindowPeakCount = nan(height(manifest), 1);
staticBaselineBackgroundFraction = nan(height(manifest), 1);
staticBaselineFullP05MvByGap = nan(height(manifest), 1);
staticBaselineFullP10MvByGap = nan(height(manifest), 1);
staticBaselineFullMedianMvByGap = nan(height(manifest), 1);
staticBaselineSourceByGap = strings(height(manifest), 1);

for ig = 1:height(manifest)
    csvPath = fullfile(libraryDir, manifest.relativeCsv(ig));
    T = readtable(csvPath, 'TextType', 'string');
    rawRows{ig} = T;
    allBladeIds = union(allBladeIds, unique(T.bladeId(:)).');
    sourceNameForGap = string(T.sourceName(find(T.sourceName ~= "", 1, 'first')));
    [baselineMv, baselineDiag] = estimate_static_gap_baseline_local( ...
        rawBinaryDir, sourceNameForGap, sourceChannel, staticBaselineMethod, ...
        bladeCount, staticBaselineRevCount, staticBaselineWindowPadS, ...
        staticBaselinePulseThresholdFraction, staticBaselineSegmentTrimFraction, ...
        staticBaselineMinSegmentS);
    staticBaselineMvByGap(ig) = baselineMv;
    staticBaselineWindowStartS(ig) = baselineDiag.windowStartS;
    staticBaselineWindowEndS(ig) = baselineDiag.windowEndS;
    staticBaselineWindowPeakCount(ig) = baselineDiag.windowPeakCount;
    staticBaselineBackgroundFraction(ig) = baselineDiag.backgroundFraction;
    staticBaselineFullP05MvByGap(ig) = baselineDiag.fullP05Mv;
    staticBaselineFullP10MvByGap(ig) = baselineDiag.fullP10Mv;
    staticBaselineFullMedianMvByGap(ig) = baselineDiag.fullMedianMv;
    staticBaselineSourceByGap(ig) = sourceNameForGap;

    waveMeta = table();
    groupKeys = unique(T(:, {'sourceName', 'gapMm', 'bladeId', 'eventNo', 'peakCenterSample'}), 'rows');
    centroidSamples = nan(height(groupKeys), 1);
    for iw = 1:height(groupKeys)
        rows = T.sourceName == groupKeys.sourceName(iw) & ...
            T.gapMm == groupKeys.gapMm(iw) & ...
            T.bladeId == groupKeys.bladeId(iw) & ...
            T.eventNo == groupKeys.eventNo(iw) & ...
            T.peakCenterSample == groupKeys.peakCenterSample(iw);
        y = T.voltageMv(rows);
        sampleIndex = T.sampleIndex(rows);
        if isempty(y)
            continue;
        end
        weight = max(y - staticBaselineMvByGap(ig), 0);
        if sum(weight, 'omitnan') > eps
            centroidSamples(iw) = sum(sampleIndex .* weight, 'omitnan') / sum(weight, 'omitnan');
        else
            centroidSamples(iw) = groupKeys.peakCenterSample(iw);
        end
    end
    waveMeta = groupKeys;
    waveMeta.centroidSample = centroidSamples;
    waveMeta = sortrows(waveMeta, 'centroidSample');
    centroidSamplesByGap{ig} = waveMeta;

    if height(waveMeta) >= bladeCount
        c = waveMeta.centroidSample(1:bladeCount);
        c = c(isfinite(c));
        if numel(c) >= 2
            bladeIntervalSamples(ig) = median(diff(sort(c)), 'omitnan');
            revPeriodS(ig) = bladeCount * bladeIntervalSamples(ig) / sampleRateHz;
            rotFreqHz(ig) = 1 / revPeriodS(ig);
            rpm(ig) = 60 * rotFreqHz(ig);
            tangentialSpeedMmS(ig) = 2 * pi * rTipMm * rotFreqHz(ig);
        end
    end
end
bladeIds = sort(allBladeIds(:).');

if any(~isfinite(tangentialSpeedMmS))
    bad = find(~isfinite(tangentialSpeedMmS));
    error('Failed to estimate calibration tangential speed for gap rows: %s', mat2str(bad(:).'));
end

if isfile(speedSummaryFile)
    speedSummary = readtable(speedSummaryFile);
    if ismember('averageTangentialSpeedMmS', speedSummary.Properties.VariableNames)
        referenceTangentialSpeedMmS = speedSummary.averageTangentialSpeedMmS(1);
        referenceRpm = speedSummary.averageRpmCentroid(1);
        speedSource = speedSummaryFile;
    else
        error('Missing averageTangentialSpeedMmS in speed summary: %s', speedSummaryFile);
    end
else
    referenceTangentialSpeedMmS = mean(tangentialSpeedMmS, 'omitnan');
    referenceRpm = mean(rpm, 'omitnan');
    speedSource = 'computed inside Step05 from all staged calibration gaps';
end

commonHalfWidthMm = fitWindowHalfWidthS * referenceTangentialSpeedMmS;
xGrid = linspace(-commonHalfWidthMm, commonHalfWidthMm, pointsPerWaveform).';

waveforms = nan(numel(xGrid), numel(bladeIds), numel(trueGaps));
sourceCounts = zeros(numel(bladeIds), numel(trueGaps));
for ig = 1:numel(trueGaps)
    T = rawRows{ig};
    speedNow = referenceTangentialSpeedMmS;
    baselineNow = staticBaselineMvByGap(ig);
    for ib = 1:numel(bladeIds)
        bid = bladeIds(ib);
        rows = T.bladeId == bid;
        if ~any(rows)
            continue;
        end
        tt = T.timeRelative_s(rows);
        yy = T.voltageMv(rows);
        [tt, idxSort] = sort(tt(:));
        yy = yy(idxSort);
        [tt, idxUnique] = unique(tt, 'stable');
        yy = yy(idxUnique);
        xDisp = tt * speedNow;
        [xDisp, idxUniqueX] = unique(xDisp, 'stable');
        yy = yy(idxUniqueX);
        yInterp = interp1(xDisp, yy, xGrid, 'linear', NaN);

        waveforms(:, ib, ig) = yInterp - baselineNow;
        sourceCounts(ib, ig) = nnz(rows);
    end
end

%% 3. Select static library aggregation and effective response window
peakByBladeGap = squeeze(max(waveforms, [], 1, 'omitnan'));
meanPeakByBlade = mean(peakByBladeGap, 2, 'omitnan');
waveformCount = sum(isfinite(waveforms), 2);
waveformsZero = waveforms;
waveformsZero(~isfinite(waveformsZero)) = 0;
waveformSum = sum(waveformsZero, 2);
YmeanAll = squeeze(waveformSum ./ max(waveformCount, 1));
YmeanAll(squeeze(waveformCount) < minBladesPerAverage) = NaN;

if strcmpi(targetBladeMode, 'average_all_blades')
    refBladeIndex = NaN;
    referenceBladeId = NaN;
    aggregationLabel = sprintf('average of %d blades', numel(bladeIds));
    YrefAll = YmeanAll;
elseif strcmpi(targetBladeMode, 'auto_reference')
    [~, refBladeIndex] = max(meanPeakByBlade);
    referenceBladeId = bladeIds(refBladeIndex);
    aggregationLabel = sprintf('reference blade %d', referenceBladeId);
    YrefAll = squeeze(waveforms(:, refBladeIndex, :));
else
    refBladeIndex = find(bladeIds == explicitReferenceBladeId, 1, 'first');
    if isempty(refBladeIndex)
        error('Explicit reference blade %d was not found in the staged library.', explicitReferenceBladeId);
    end
    referenceBladeId = bladeIds(refBladeIndex);
    aggregationLabel = sprintf('reference blade %d', referenceBladeId);
    YrefAll = squeeze(waveforms(:, refBladeIndex, :));
end
validGapMask = all(isfinite(YrefAll), 1).';
gTrain = trueGaps(validGapMask);
recordedGapTrain = recordedGaps(validGapMask);
Yref = YrefAll(:, validGapMask);

if numel(gTrain) < minValidGapsForFit
    error('Too few valid gaps for response-surface fitting: %d.', numel(gTrain));
end

peakEnvelope = max(Yref, [], 2, 'omitnan');
effectiveWindow = peakEnvelope >= effectiveWindowThreshold * max(peakEnvelope, [], 'omitnan');
idxWin = find(effectiveWindow);
if ~isempty(idxWin)
    effectiveWindow(idxWin(1):idxWin(end)) = true;
end

%% 4. Fit F(g,x) = B0 + B1/g + B2 log(g/g0)
basis = [ones(numel(gTrain), 1), 1 ./ gTrain(:), log(gTrain(:) ./ g0Mm)];
coeff = nan(numel(xGrid), 3);
Yfit = nan(size(Yref));
for ix = 1:numel(xGrid)
    if ~effectiveWindow(ix)
        continue;
    end
    y = Yref(ix, :).';
    if nnz(isfinite(y)) < minValidGapsForFit
        continue;
    end
    coeff(ix, :) = (basis \ y).';
    Yfit(ix, :) = (basis * coeff(ix, :).').';
end

fitRows = cell(numel(gTrain), 1);
for ig = 1:numel(gTrain)
    keep = effectiveWindow & isfinite(Yfit(:, ig)) & isfinite(Yref(:, ig));
    residual = Yref(keep, ig) - Yfit(keep, ig);
    rmseMv = sqrt(mean(residual.^2, 'omitnan'));
    relRmse = rmseMv / max(abs(Yref(keep, ig)), [], 'omitnan');
    fitRows{ig} = table(recordedGapTrain(ig), gTrain(ig), rmseMv, relRmse, ...
        'VariableNames', {'recordedGapMm', 'trueReferenceGapMm', 'rmseMv', 'relativeRmse'});
end
fitTable = vertcat(fitRows{:});

dFdgGrid = nan(numel(xGrid), numel(gTrain));
dFdxGrid = nan(numel(xGrid), numel(gTrain));
for ig = 1:numel(gTrain)
    g = gTrain(ig);
    dFdgGrid(:, ig) = -coeff(:, 2) ./ (g .^ 2) + coeff(:, 3) ./ g;
    dFdxGrid(:, ig) = gradient(Yfit(:, ig), xGrid);
end

responseSurface = struct();
responseSurface.dataset = '20250527';
responseSurface.libraryDir = libraryDir;
responseSurface.manifestFile = manifestFile;
responseSurface.sourceNote = 'Built from static_gap_waveform_library_20260524_2000Hz staged multi-gap waveforms.';
responseSurface.gapZeroOffsetMm = gapZeroOffsetMm;
responseSurface.sampleRateHz = sampleRateHz;
responseSurface.rTipMm = rTipMm;
responseSurface.sourceChannel = sourceChannel;
responseSurface.sourceChannelIndexing = 'MATLAB 1-based channel index';
responseSurface.sourceChannelNote = 'The staged static-gap CSV was generated with selectedChannel=1; original visualization labels this as CH0.';
responseSurface.staticBaselineMethod = staticBaselineMethod;
responseSurface.staticBaselineByGap = table(recordedGaps(:), trueGaps(:), ...
    staticBaselineMvByGap(:), staticBaselineWindowStartS(:), ...
    staticBaselineWindowEndS(:), staticBaselineWindowPeakCount(:), ...
    staticBaselineBackgroundFraction(:), staticBaselineFullP05MvByGap(:), ...
    staticBaselineFullP10MvByGap(:), staticBaselineFullMedianMvByGap(:), ...
    staticBaselineSourceByGap(:), ...
    'VariableNames', {'recordedGapMm', 'trueReferenceGapMm', ...
    'baselineMv', 'baselineWindowStartS', 'baselineWindowEndS', ...
    'baselineWindowPeakCount', 'backgroundFraction', 'fullP05Mv', ...
    'fullP10Mv', 'fullMedianMv', 'sourceName'});
responseSurface.referenceRpm = referenceRpm;
responseSurface.referenceTangentialSpeedMmS = referenceTangentialSpeedMmS;
responseSurface.referenceTangentialSpeedSource = speedSource;
responseSurface.recordedGapMm = recordedGaps(:);
responseSurface.trueGapMm = trueGaps(:);
responseSurface.gTrainMm = gTrain(:);
responseSurface.recordedGapTrainMm = recordedGapTrain(:);
responseSurface.g0Mm = g0Mm;
responseSurface.xGrid = xGrid;
responseSurface.xUnit = 'mm circumferential displacement from waveform centroid/peak coordinate';
responseSurface.fitWindowHalfWidthS = fitWindowHalfWidthS;
responseSurface.commonHalfWidthMm = commonHalfWidthMm;
responseSurface.voltageUnit = 'mV baseline-corrected using one raw continuous-waveform baseline per recorded static gap';
responseSurface.bladeIds = bladeIds;
responseSurface.targetBladeMode = targetBladeMode;
responseSurface.waveformAggregationLabel = aggregationLabel;
responseSurface.minBladesPerAverage = minBladesPerAverage;
responseSurface.referenceBladeId = referenceBladeId;
responseSurface.referenceBladeIndex = refBladeIndex;
responseSurface.waveforms = waveforms;
responseSurface.waveformCount = waveformCount;
responseSurface.YmeanAll = YmeanAll;
responseSurface.Yref = Yref;
responseSurface.Yfit = Yfit;
responseSurface.coeff = coeff;
responseSurface.effectiveWindow = effectiveWindow;
responseSurface.dFdgGrid = dFdgGrid;
responseSurface.dFdxGrid = dFdxGrid;
responseSurface.fitTable = fitTable;
responseSurface.sourceCounts = sourceCounts;
responseSurface.calibrationMotion = table(recordedGaps(:), trueGaps(:), ...
    bladeIntervalSamples(:), revPeriodS(:), rotFreqHz(:), rpm(:), ...
    tangentialSpeedMmS(:), ...
    'VariableNames', {'recordedGapMm', 'trueReferenceGapMm', ...
    'bladeIntervalSamples', 'revPeriodS', 'rotFreqHz', 'rpm', ...
    'tangentialSpeedMmS'});
responseSurface.centroidSamplesByGap = centroidSamplesByGap;
responseSurface.note = ['Step05 response surface uses the staged static gap library. ' ...
    'Under the current assumption that blade-thickness differences have only a small effect on waveform shape, ' ...
    'the default aggregation uses the multi-blade average waveform across static gaps. ' ...
    'The response coordinate is circumferential displacement computed using ' ...
    'the average calibration tip speed from Step05A. ' ...
    'Dynamic highMap coordinates should be mapped to the same circumferential-displacement coordinate in Step06; ' ...
    'dynamic clearance estimates should be treated as effective until a same-rig multi-gap calibration is available.'];

%% 5. Save outputs
matFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
csvFile = fullfile(outDir, 'Step05_Response_Surface_Fit_20250527.csv');
save(matFile, 'responseSurface', '-v7.3');
motionCsvFile = fullfile(outDir, 'Step05_Calibration_Motion_20250527.csv');
baselineCsvFile = fullfile(outDir, 'Step05_Static_Gap_Continuous_Baseline_20250527.csv');
writetable(fitTable, csvFile);
writetable(responseSurface.calibrationMotion, motionCsvFile);
writetable(responseSurface.staticBaselineByGap, baselineCsvFile);

%% 6. Diagnostics
if saveFigures
    fig = figure('Name', '20250527 Step05 response surface', ...
        'Color', 'w', 'Position', [80, 80, 1350, 880], 'NumberTitle', 'off');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    colors = turbo(numel(gTrain));
    nexttile; hold on; grid on; box on;
    for ig = 1:numel(gTrain)
        plot(xGrid, Yref(:, ig), '-', 'Color', colors(ig, :), 'LineWidth', 0.8);
        plot(xGrid, Yfit(:, ig), '--', 'Color', colors(ig, :), 'LineWidth', 1.0);
    end
    yl = ylim;
    patch([xGrid(find(effectiveWindow, 1, 'first')), xGrid(find(effectiveWindow, 1, 'last')), ...
        xGrid(find(effectiveWindow, 1, 'last')), xGrid(find(effectiveWindow, 1, 'first'))], ...
        [yl(1), yl(1), yl(2), yl(2)], [1.0, 0.90, 0.45], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(findobj(gca, 'Type', 'line'), 'top');
    xlabel('Circumferential displacement x (mm)');
    ylabel('Voltage (mV)');
    title(sprintf('%s: measured solid, fitted dashed', aggregationLabel));

    nexttile; hold on; grid on; box on;
    plot(fitTable.trueReferenceGapMm, fitTable.rmseMv, 'o-', 'LineWidth', 1.2);
    xlabel('True reference gap (mm)');
    ylabel('RMSE (mV)');
    title('Fit error in effective window');

    nexttile; hold on; grid on; box on;
    imagesc(gTrain, xGrid, Yfit);
    axis xy; colorbar;
    xlabel('Gap (mm)');
    ylabel('Circumferential displacement x (mm)');
    title('Fitted response surface F(g,x)');

    nexttile; hold on; grid on; box on;
    plot(xGrid, dFdgGrid(:, max(1, round(numel(gTrain)/2))), 'LineWidth', 1.2);
    plot(xGrid, dFdxGrid(:, max(1, round(numel(gTrain)/2))), 'LineWidth', 1.2);
    xlabel('Circumferential displacement x (mm)');
    ylabel('Derivative');
    title('Representative sensitivities');
    legend('dF/dg', 'dF/dx', 'Location', 'best');

    figFile = fullfile(outDir, 'Step05_Response_Surface_20250527.png');
    saveas(fig, figFile);
end

fprintf('\nStep05 complete.\n');
fprintf('Response surface saved to:\n  %s\n', matFile);
fprintf('Fit summary saved to:\n  %s\n', csvFile);
fprintf('Calibration speed summary saved to:\n  %s\n', motionCsvFile);
fprintf('Static continuous-baseline summary saved to:\n  %s\n', baselineCsvFile);

%% Local functions
function [baselineMv, diagInfo] = estimate_static_gap_baseline_local(rawBinaryDir, sourceName, selectedChannel, method, ...
    bladeCount, revCount, windowPadS, pulseThresholdFraction, segmentTrimFraction, minSegmentS)
binPath = fullfile(rawBinaryDir, sourceName + ".bin");
headerPath = fullfile(rawBinaryDir, sourceName + "_header.txt");
if ~isfile(binPath) || ~isfile(headerPath)
    error('Missing raw static calibration binary/header for %s.', sourceName);
end
[~, yMv, fs] = read_acts1000_binary_channel_local(binPath, headerPath, selectedChannel);
diagInfo = struct();
switch lower(method)
    case 'initial_2rev_trimmed_background_median'
        [windowMask, windowInfo] = select_initial_revolution_window_local( ...
            yMv, fs, bladeCount, revCount, windowPadS);
        bgDiag = compute_trimmed_background_baseline_local(yMv, fs, windowMask, ...
            pulseThresholdFraction, windowPadS, segmentTrimFraction, minSegmentS);
        baselineMv = bgDiag.backgroundMedianMv;
        diagInfo = bgDiag;
        diagInfo.windowStartS = windowInfo.startTimeS;
        diagInfo.windowEndS = windowInfo.endTimeS;
        diagInfo.windowPeakCount = windowInfo.peakCount;
    case 'raw_continuous_p05'
        diagInfo.fullP05Mv = prctile(yMv, 5);
        diagInfo.fullP10Mv = prctile(yMv, 10);
        diagInfo.fullMedianMv = median(yMv, 'omitnan');
        baselineMv = diagInfo.fullP05Mv;
    case 'raw_continuous_p10'
        diagInfo.fullP05Mv = prctile(yMv, 5);
        diagInfo.fullP10Mv = prctile(yMv, 10);
        diagInfo.fullMedianMv = median(yMv, 'omitnan');
        baselineMv = diagInfo.fullP10Mv;
    case 'raw_continuous_median'
        diagInfo.fullP05Mv = prctile(yMv, 5);
        diagInfo.fullP10Mv = prctile(yMv, 10);
        diagInfo.fullMedianMv = median(yMv, 'omitnan');
        baselineMv = diagInfo.fullMedianMv;
    otherwise
        error('Unknown static baseline method: %s', method);
end
end

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
info = struct('startTimeS', t(startIdx), 'endTimeS', t(endIdx), ...
    'startIdx', startIdx, 'endIdx', endIdx, 'peakCount', peakCount);
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

function diagInfo = compute_trimmed_background_baseline_local(yMv, fs, windowMask, pulseThresholdFraction, pulsePadS, segmentTrimFraction, minSegmentS)
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
backgroundMask = trim_background_segments_local(backgroundMask0, fs, segmentTrimFraction, minSegmentS);
if nnz(backgroundMask) < 0.10 * nnz(backgroundMask0)
    warning('Trimmed background mask is too small; falling back to the untrimmed no-pulse mask.');
    backgroundMask = backgroundMask0;
end
bg = yMv(backgroundMask);
diagInfo = struct();
diagInfo.fullP05Mv = fullP05Mv;
diagInfo.fullP10Mv = fullP10Mv;
diagInfo.fullMedianMv = fullMedianMv;
diagInfo.backgroundMeanMv = mean(bg, 'omitnan');
diagInfo.backgroundMedianMv = median(bg, 'omitnan');
diagInfo.backgroundP05Mv = prctile(bg, 5);
diagInfo.backgroundStdMv = std(bg, 'omitnan');
diagInfo.backgroundFraction = nnz(backgroundMask) / nnz(windowMask);
diagInfo.thresholdMv = thresholdMv;
end

function trimmedMask = trim_background_segments_local(mask, fs, trimFraction, minSegmentS)
mask = mask(:);
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
