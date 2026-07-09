%% Step04: resonance detection and visualization dashboard for 20241106
% Unified Step04 workflow:
%   1. aligned strain time response and OPR speed variation
%   2. raw segmented FFT with direct 100-microstrain resonance regions
%   3. resonance-region summary
%
% Figures are opened directly. The selected-channel regions are saved using
% the Step04 filenames and the legacy Step05 CSV name so downstream scripts
% can keep reading the same file.

clear; clc; close all;

%% 1. Simple settings
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

dataRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
highDir = fullfile(dataRoot, '3000_3150');

strainFiles = { ...
    'AI1-01_20241106170148.mat'
    'AI1-02_20241106170148.mat'
    'AI1-03_20241106170148.mat'
    'AI1-04_20241106170148.mat'
    };
strainLabels = {'AI1-01','AI1-02','AI1-03','AI1-04'};

strainChannelsToAnalyze = 1:numel(strainFiles);
selectedStrainChannel = 4;      % Use AI1-04 to match Untitled2.m and fft1.m.
detrendMethod = 'linear';       % remove slow strain drift before STFT.
sampleRateOprHz = 5e6;
strainSampleRateHz = 10000;
strainToBttOffsetSec = -102.6;

stftWindowSec = 2.0;
stftStepSec = 0.25;
stftFreqBandHz = [0 1000];

orders = 6:18;
orderHalfBandHz = 3.0;
regionMergeGapSec = 0.00;
minRegionWidthSec = 0.00;
maxRegionWidthSec = inf;
resonanceFftAmpThresholdMicrostrain = 100.0;
resonanceTimeAmpThresholdMicrostrain = 100.0;
resonanceLowThresholdRatio = 0.55;
resonanceCombineMode = "fft_surface_only";
preferredRegionTimeSec = [75, 90];

dashboardFreqBandHz = [0 1000];
make3DFigure = true;

% Raw segmented FFT used as the primary strain check. This matches
% fft1.m/Untitled2.m: raw AI1 data, rectangular segments, microstrain
% amplitude, no detrend, no Hann window, no dB conversion.
makeFft1Style3DFigure = true;
fft1StyleSampleRange = [1 2700000];
fft1StyleNumSegments = 100;
fft1StyleFreqBandHz = [0 800];
fft1StylePeakBandHz = [20 800];
fft1StyleZLim = [0 250];

% Local direct FFT check. Enter the window in strain raw time so it matches
% the original AI1 scripts directly. Untitled2.m uses raw strain time
% [170 180] s, and fft1.m uses [150 210] s. The dashboard axis is BTT time,
% so the window is converted below by strainToBttOffsetSec.
makeLocal2DFFT = true;
localFftStrainTimeWindowSec = [170 180];
localFftBttTimeWindowSec = localFftStrainTimeWindowSec + strainToBttOffsetSec;
localFftFreqBandHz = [0 1000];

%% 2. Analyze each strain channel separately
opr = load_opr_local(highDir, sampleRateOprHz);
channelResults = repmat(empty_channel_result_local(), numel(strainChannelsToAnalyze), 1);
scoreTables = cell(numel(strainChannelsToAnalyze), 1);
peakTables = cell(numel(strainChannelsToAnalyze), 1);
regionTables = cell(numel(strainChannelsToAnalyze), 1);
localFftTables = cell(numel(strainChannelsToAnalyze), 1);
summaryTables = cell(numel(strainChannelsToAnalyze), 1);

for ii = 1:numel(strainChannelsToAnalyze)
    channelIdx = strainChannelsToAnalyze(ii);
    strainFile = fullfile(dataRoot, strainFiles{channelIdx});
    strainNow = load_strain_local(strainFile, strainLabels{channelIdx}, ...
        strainSampleRateHz, detrendMethod);
    rawFftNow = raw_segmented_fft_local(strainNow, strainSampleRateHz, ...
        fft1StyleSampleRange, fft1StyleNumSegments, fft1StylePeakBandHz);
    stftNow = strain_stft_local(strainNow.time, strainNow.value, ...
        strainSampleRateHz, stftWindowSec, stftStepSec, stftFreqBandHz);

    tBttNow = stftNow.timeSec(:).' + strainToBttOffsetSec;
    validTimeRangeNow = [max(min(tBttNow), min(opr.rpmTime)), ...
        min(max(tBttNow), max(opr.rpmTime))];

    orderEnergyNow = sample_order_energy_local(stftNow, tBttNow, opr, ...
        orders, orderHalfBandHz);
    [scoreTableNow, peakTableNow, regionTableNow, localFftTableNow, localFftResultNow] = detect_resonance_regions_from_raw_fft_local( ...
        rawFftNow, orderEnergyNow, orders, strainToBttOffsetSec, regionMergeGapSec, ...
        resonanceFftAmpThresholdMicrostrain, resonanceTimeAmpThresholdMicrostrain, ...
        resonanceLowThresholdRatio, resonanceCombineMode, ...
        minRegionWidthSec, maxRegionWidthSec, preferredRegionTimeSec);
    peakTableNow = sortrows(peakTableNow, 'ampDb', 'descend');
    peakTableNow = add_strain_labels_local(peakTableNow, channelIdx, ...
        strainLabels{channelIdx});
    regionTableNow = sortrows(regionTableNow, {'preferredOverlapSec','regionPeakScore'}, {'descend','descend'});
    regionTableNow = add_strain_labels_local(regionTableNow, channelIdx, ...
        strainLabels{channelIdx});
    scoreTableNow = add_strain_labels_local(scoreTableNow, channelIdx, ...
        strainLabels{channelIdx});
    localFftTableNow = add_strain_labels_local(localFftTableNow, channelIdx, ...
        strainLabels{channelIdx});

    channelResults(ii).channelIdx = channelIdx;
    channelResults(ii).strain = strainNow;
    channelResults(ii).rawFft = rawFftNow;
    channelResults(ii).stft = stftNow;
    channelResults(ii).tBtt = tBttNow;
    channelResults(ii).validTimeRange = validTimeRangeNow;
    channelResults(ii).orderEnergy = orderEnergyNow;
    channelResults(ii).scoreTable = scoreTableNow;
    channelResults(ii).peakTable = peakTableNow;
    channelResults(ii).regionTable = regionTableNow;
    channelResults(ii).localFftTable = localFftTableNow;
    channelResults(ii).localFftResult = localFftResultNow;

    scoreTables{ii} = scoreTableNow;
    peakTables{ii} = peakTableNow;
    regionTables{ii} = regionTableNow;
    localFftTables{ii} = localFftTableNow;
    summaryTables{ii} = summarize_channel_local(channelIdx, ...
        strainLabels{channelIdx}, peakTableNow, regionTableNow);
end

allScoreTable = vertcat(scoreTables{:});
allPeakTable = vertcat(peakTables{:});
allRegionTable = vertcat(regionTables{:});
allLocalFftTable = vertcat(localFftTables{:});
channelSummary = vertcat(summaryTables{:});

if isnumeric(selectedStrainChannel)
    selectedIdx = find(strainChannelsToAnalyze == selectedStrainChannel, 1);
    if isempty(selectedIdx)
        error('selectedStrainChannel is not in strainChannelsToAnalyze.');
    end
else
    [~, selectedIdx] = max(channelSummary.maxPeakAmpDb);
end

strain = channelResults(selectedIdx).strain;
rawFft = channelResults(selectedIdx).rawFft;
stft = channelResults(selectedIdx).stft;
tBtt = channelResults(selectedIdx).tBtt;
validTimeRange = channelResults(selectedIdx).validTimeRange;
orderEnergy = channelResults(selectedIdx).orderEnergy;
scoreTable = channelResults(selectedIdx).scoreTable;
peakTable = channelResults(selectedIdx).peakTable;
regionTable = channelResults(selectedIdx).regionTable;
localFftTable = channelResults(selectedIdx).localFftTable;
localFftResult = channelResults(selectedIdx).localFftResult;
resonanceParams = localFftResult.params;

%% 3. Print multi-channel summary
fprintf('\n=== Step04 multi-strain resonance dashboard ===\n');
fprintf('Detrend method: %s\n', detrendMethod);
fprintf('strainToBttOffsetSec = %+0.3f s\n', strainToBttOffsetSec);
fprintf('Local FFT strain raw time window = %.3f-%.3f s, BTT window = %.3f-%.3f s\n', ...
    localFftStrainTimeWindowSec(1), localFftStrainTimeWindowSec(2), ...
    localFftBttTimeWindowSec(1), localFftBttTimeWindowSec(2));
fprintf('Orders: %s\n', mat2str(orders));
fprintf(['Resonance rule: raw segmented FFT surface max >= %.1f microstrain; ', ...
    'continuous above-threshold intervals are resonance regions.\n'], ...
    resonanceFftAmpThresholdMicrostrain);
fprintf(['Raw FFT settings: peak band %.0f-%.0f Hz, segments %d, ', ...
    'merge gap %.1f s, max width %s s\n\n'], ...
    fft1StylePeakBandHz(1), fft1StylePeakBandHz(2), ...
    rawFft.numSegments, regionMergeGapSec, string_or_inf_local(maxRegionWidthSec));

fprintf('Channel summary:\n');
disp(channelSummary);

fprintf('\nAll-channel region summary:\n');
disp(allRegionTable);

fprintf('\nSelected detail channel: %s\n', strain.label);
fprintf('Selected-channel region summary:\n');
disp(regionTable);

fprintf('\nSelected-channel peak table used for regions:\n');
disp(peakTable);

%% 4. Channel comparison and selected-channel dashboard
plot_channel_comparison_local(channelSummary, allPeakTable, dashboardFreqBandHz);
plot_dashboard_local(strain, rawFft, strainToBttOffsetSec, opr, stft, tBtt, ...
    orderEnergy, orders, peakTable, regionTable, validTimeRange, dashboardFreqBandHz);

%% 5. Optional 3D auxiliary view
if make3DFigure
    plot_3d_auxiliary_local(rawFft, strainToBttOffsetSec, validTimeRange, ...
        dashboardFreqBandHz, regionTable);
end

%% 6. Optional original-script 3D FFT view
if makeFft1Style3DFigure
    plot_fft1_style_3d_fft_local(strain, rawFft, strainSampleRateHz, fft1StyleSampleRange, ...
        fft1StyleNumSegments, fft1StyleFreqBandHz, fft1StyleZLim);
end

%% 7. Optional local 2D FFT view for a user-selected steady region
if makeLocal2DFFT
    plot_local_2d_fft_local(strain, strainToBttOffsetSec, opr, ...
        localFftBttTimeWindowSec, strainSampleRateHz, localFftFreqBandHz, orders);
end

matFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20241106.mat');
csvFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20241106.csv');
legacyCsvFile = fullfile(outDir, 'Step05_ResonanceRegions_20241106.csv');
scoreCsvFile = fullfile(outDir, 'Step04_DirectStrainResonance_Score_20241106.csv');
localFftCsvFile = fullfile(outDir, 'Step04_LocalStrainFFT_20241106.csv');
save(matFile, 'scoreTable', 'peakTable', 'regionTable', 'localFftTable', ...
    'localFftResult', 'rawFft', 'strain', 'resonanceParams', ...
    'allScoreTable', 'allPeakTable', 'allRegionTable', 'allLocalFftTable', ...
    'channelSummary', 'strainToBttOffsetSec', 'selectedStrainChannel', '-v7');
writetable(regionTable, csvFile);
writetable(regionTable, legacyCsvFile);
writetable(scoreTable, scoreCsvFile);
writetable(localFftTable, localFftCsvFile);

fprintf('\nSaved Step04 resonance outputs:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    matFile, csvFile, legacyCsvFile, scoreCsvFile, localFftCsvFile);

%% Local functions
function R = empty_channel_result_local()
R = struct('channelIdx', NaN, 'strain', [], 'rawFft', [], 'stft', [], 'tBtt', [], ...
    'validTimeRange', [], 'orderEnergy', [], 'scoreTable', table(), ...
    'peakTable', table(), 'regionTable', table(), 'localFftTable', table(), ...
    'localFftResult', struct());
end

function opr = load_opr_local(caseDir, sampleRateHz)
oprFile = fullfile(caseDir, 'jiluOPR.mat');
if exist(oprFile, 'file') ~= 2
    error('Missing OPR file: %s', oprFile);
end
S = load(oprFile, 'jiluOPR');
raw = S.jiluOPR;
opr.tCenter = mean(raw(:, 1:2), 2) / sampleRateHz;
opr.rpmTime = opr.tCenter(1:end-1);
opr.rpm = 60 ./ diff(opr.tCenter);
end

function strain = load_strain_local(filePath, label, Fs, detrendMethod)
if exist(filePath, 'file') ~= 2
    error('Missing strain file: %s', filePath);
end
S = load(filePath, 'Datas');
strain.label = label;
strain.file = filePath;
strain.time = S.Datas(:, 1);
strain.rawValue = S.Datas(:, 2);
strain.value = preprocess_strain_local(strain.rawValue, Fs, detrendMethod);
strain.detrendMethod = detrendMethod;
end

function y = preprocess_strain_local(x, Fs, detrendMethod)
x = x(:);
x(~isfinite(x)) = median(x(isfinite(x)), 'omitnan');
switch lower(string(detrendMethod))
    case "none"
        y = x - median(x, 'omitnan');
    case "linear"
        y = detrend(x, 'linear');
    case "movingmedian"
        win = max(5, round(2.0 * Fs));
        y = x - movmedian(x, win, 'omitnan');
    otherwise
        error('Unknown detrendMethod: %s', detrendMethod);
end
y = y - median(y, 'omitnan');
end

function rawFft = raw_segmented_fft_local(strain, Fs, sampleRange, numSegments, peakBandHz)
xRaw = strain.rawValue(:);
tRaw = strain.time(:);
nTotal = numel(xRaw);
i1 = max(1, round(sampleRange(1)));
i2 = min(nTotal, round(sampleRange(2)));
if i2 <= i1
    error('Raw FFT sample range is invalid: %d:%d.', i1, i2);
end

idx = (i1:i2).';
x = xRaw(idx);
t = tRaw(idx);
N = numel(x);
segmentLength = floor(N / numSegments);
if segmentLength < 2
    error('Raw FFT segment length is too short: %d samples.', segmentLength);
end

Nused = segmentLength * numSegments;
x = x(1:Nused);
t = t(1:Nused);
freq = (0:floor(segmentLength / 2)).' * (Fs / segmentLength);
amplitudes = zeros(numSegments, numel(freq));
segmentCenterTimeSec = zeros(numSegments, 1);
segmentPeakFreqHz = zeros(numSegments, 1);
segmentPeakAmpMicrostrain = zeros(numSegments, 1);
segmentDetrendedMaxAbsMicrostrain = zeros(numSegments, 1);
segmentDetrendedP95AbsMicrostrain = zeros(numSegments, 1);
peakMask = freq >= peakBandHz(1) & freq <= peakBandHz(2);
freqForPeak = freq(peakMask);

for i = 1:numSegments
    localIdx = (i - 1) * segmentLength + (1:segmentLength);
    segmentData = x(localIdx);
    segmentDetrended = detrend(segmentData, 'linear');
    Y = fft(segmentData);
    P2 = abs(Y / segmentLength);
    P1 = P2(1:(floor(segmentLength / 2) + 1));
    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end
    P1(1) = 0;

    amplitudes(i, :) = P1(:).';
    segmentCenterTimeSec(i) = mean(t(localIdx), 'omitnan');
    [segmentPeakAmpMicrostrain(i), peakIdx] = max(P1(peakMask));
    segmentPeakFreqHz(i) = freqForPeak(peakIdx);
    absDetrended = abs(segmentDetrended);
    segmentDetrendedMaxAbsMicrostrain(i) = max(absDetrended, [], 'omitnan');
    segmentDetrendedP95AbsMicrostrain(i) = percentile_local(absDetrended, 95);
end

rawFft.sampleStart = i1;
rawFft.sampleEnd = i1 + Nused - 1;
rawFft.segmentLength = segmentLength;
rawFft.numSegments = numSegments;
rawFft.Nused = Nused;
rawFft.freqHz = freq;
rawFft.segmentCenterTimeSec = segmentCenterTimeSec;
rawFft.amplitudesMicrostrain = amplitudes;
rawFft.peakBandHz = peakBandHz;
rawFft.segmentPeakFreqHz = segmentPeakFreqHz;
rawFft.segmentPeakAmpMicrostrain = segmentPeakAmpMicrostrain;
rawFft.segmentDetrendedMaxAbsMicrostrain = segmentDetrendedMaxAbsMicrostrain;
rawFft.segmentDetrendedP95AbsMicrostrain = segmentDetrendedP95AbsMicrostrain;
end

function stft = strain_stft_local(t, x, Fs, windowSec, stepSec, freqBandHz)
t = t(:);
x = x(:);
winN = max(8, round(windowSec * Fs));
stepN = max(1, round(stepSec * Fs));
nfft = 2 ^ nextpow2(winN);
freq = (0:(nfft/2)).' * Fs / nfft;
freqMask = freq >= freqBandHz(1) & freq <= freqBandHz(2);
freqKeep = freq(freqMask);

starts = 1:stepN:(numel(x) - winN + 1);
amp = NaN(numel(freqKeep), numel(starts));
timeSec = NaN(1, numel(starts));
w = hann_local(winN);
gain = mean(w);

for i = 1:numel(starts)
    idx = starts(i):(starts(i) + winN - 1);
    xw = x(idx);
    xw = xw - mean(xw, 'omitnan');
    y = fft(xw .* w, nfft);
    a2 = abs(y / (winN * gain));
    a1 = a2(1:(nfft/2 + 1));
    if numel(a1) > 2
        a1(2:end-1) = 2 * a1(2:end-1);
    end
    amp(:, i) = a1(freqMask);
    timeSec(i) = mean(t(idx([1, end])));
end

stft.timeSec = timeSec;
stft.freqHz = freqKeep;
stft.amp = amp;
stft.windowSec = windowSec;
stft.stepSec = stepSec;
end

function orderEnergy = sample_order_energy_local(stft, tBtt, opr, orders, halfBandHz)
nT = numel(stft.timeSec);
nO = numel(orders);
amp = NaN(nT, nO);
ampDb = NaN(nT, nO);
orderFreqHz = NaN(nT, nO);
rpm = interp1(opr.rpmTime, opr.rpm, tBtt, 'linear', NaN);

for j = 1:nO
    fLine = orders(j) * rpm / 60;
    orderFreqHz(:, j) = fLine(:);
    for i = 1:nT
        f0 = fLine(i);
        if ~isfinite(f0) || f0 < min(stft.freqHz) || f0 > max(stft.freqHz)
            continue;
        end
        fMask = stft.freqHz >= f0 - halfBandHz & stft.freqHz <= f0 + halfBandHz;
        if ~any(fMask)
            [~, idx] = min(abs(stft.freqHz - f0));
            amp(i, j) = stft.amp(idx, i);
        else
            amp(i, j) = max(stft.amp(fMask, i), [], 'omitnan');
        end
    end
    ampDb(:, j) = 20 * log10(amp(:, j) + eps);
end

orderEnergy.tBtt = tBtt(:);
orderEnergy.strainTimeSec = stft.timeSec(:);
orderEnergy.rpm = rpm(:);
orderEnergy.amp = amp;
orderEnergy.ampDb = ampDb;
orderEnergy.orderFreqHz = orderFreqHz;
end

function [S, P, R, L, localFftResult] = detect_resonance_regions_from_raw_fft_local(rawFft, E, orders, ...
    offsetSec, mergeGapSec, fftAmpThr, timeAmpThr, lowRatio, combineMode, ...
    minRegionWidthSec, maxRegionWidthSec, preferredRegionTimeSec)
t = rawFft.segmentCenterTimeSec(:) + offsetSec;
fftAmp = rawFft.segmentPeakAmpMicrostrain(:);
timeAmp = rawFft.segmentDetrendedMaxAbsMicrostrain(:);
timeP95 = rawFft.segmentDetrendedP95AbsMicrostrain(:);
freq = rawFft.segmentPeakFreqHz(:);
rpmRaw = interp1(E.tBtt, E.rpm, t, 'linear', NaN);
valid = isfinite(t) & isfinite(fftAmp) & isfinite(timeAmp) & isfinite(freq) & isfinite(rpmRaw);
surfaceAmp = fftAmp;
timeAmpOut = timeAmp;
score = surfaceAmp;
highMask = valid & surfaceAmp >= fftAmpThr;
lowMask = highMask;

shaftFreqHz = rpmRaw / 60;
orderApprox = freq ./ shaftFreqHz;
nearestOrder = NaN(size(orderApprox));
syncFreqAtOrderHz = NaN(size(orderApprox));
for i = 1:numel(orderApprox)
    if ~isfinite(orderApprox(i))
        continue;
    end
    [~, idxOrder] = min(abs(orders - orderApprox(i)));
    nearestOrder(i) = orders(idxOrder);
    syncFreqAtOrderHz(i) = nearestOrder(i) * shaftFreqHz(i);
end

S = table(t, rpmRaw, freq, surfaceAmp, timeAmpOut, timeP95, ...
    orderApprox, nearestOrder, syncFreqAtOrderHz, score, highMask, lowMask, ...
    'VariableNames', {'bttTimeSec','rpm','peakFreqHz', ...
    'peakFftAmpMicrostrain','peakTimeAmpMicrostrain', ...
    'peakTimeP95AbsMicrostrain','fftPeakOrderApprox','nearestOrder', ...
    'syncFreqAtNearestOrderHz','resonanceScore','isHighMask','isLowMask'});

edge = diff([false; S.isHighMask(:); false]);
startIdx = find(edge == 1);
endIdx = find(edge == -1) - 1;
rows = [];
for k = 1:numel(startIdx)
    segIdx = startIdx(k):endIdx(k);
    [~, rel] = max(S.peakFftAmpMicrostrain(segIdx));
    ip = segIdx(rel);
    i0 = startIdx(k);
    i1 = endIdx(k);
    dt = median(diff(t(valid)), 'omitnan');
    if ~isfinite(dt) || dt <= 0
        dt = rawFft.segmentLength / max(rawFft.Nused / range(rawFft.segmentCenterTimeSec), eps);
    end
    t0 = t(i0) - 0.5 * dt;
    t1 = t(i1) + 0.5 * dt;
    widthSec = t1 - t0;
    if widthSec < minRegionWidthSec
        continue;
    end
    if widthSec > maxRegionWidthSec
        halfW = 0.5 * maxRegionWidthSec;
        t0 = t(ip) - halfW;
        t1 = t(ip) + halfW;
        widthSec = maxRegionWidthSec;
    end

    peakOrder = S.nearestOrder(ip);
    if ~isfinite(peakOrder)
        peakOrder = nearest_finite_order_local(S.fftPeakOrderApprox(ip), orders);
    end
    peakAmpDb = 20 * log10(max(S.peakFftAmpMicrostrain(ip), eps));
    idx = find(S.bttTimeSec >= t0 & S.bttTimeSec <= t1);
    rows = [rows; k, t0, t1, t0 - offsetSec, t1 - offsetSec, widthSec, ...
        S.bttTimeSec(ip), S.bttTimeSec(ip) - offsetSec, ...
        S.peakFreqHz(ip), peakOrder, peakAmpDb, S.peakFftAmpMicrostrain(ip), ...
        S.peakTimeAmpMicrostrain(ip), S.rpm(ip), S.syncFreqAtNearestOrderHz(ip), ...
        S.fftPeakOrderApprox(ip), S.resonanceScore(ip), numel(idx), ...
        double(overlap_duration_local([t0, t1], preferredRegionTimeSec))]; %#ok<AGROW>
end

if isempty(rows)
    R = empty_region_table_local();
else
    R = array2table(rows, 'VariableNames', {'regionId','bttStartSec','bttEndSec', ...
        'strainStartSec','strainEndSec','regionDuration','peakBttTimeSec', ...
        'peakStrainTimeSec','dominantFreqHz','dominantOrder','peakAmpDb', ...
        'peakFftAmpMicrostrain','peakTimeAmpMicrostrain','peakRpm', ...
        'syncFreqAtRegionOrderHz','fftPeakOrderApprox','regionPeakScore', ...
        'peakCount','preferredOverlapSec'});
    R = merge_score_regions_local(R, mergeGapSec);
    R.ordersIncluded = strings(height(R), 1);
    for i = 1:height(R)
        inRegion = S.bttTimeSec >= R.bttStartSec(i) & S.bttTimeSec <= R.bttEndSec(i) & ...
            isfinite(S.nearestOrder);
        ord = unique(S.nearestOrder(inRegion)).';
        if isempty(ord)
            R.ordersIncluded(i) = "";
        else
            R.ordersIncluded(i) = join(string(ord), ",");
        end
    end
    R.regionStart = R.bttStartSec;
    R.regionEnd = R.bttEndSec;
    R.regionPeakTime = R.peakBttTimeSec;
    R.regionPeakFreqHz = R.dominantFreqHz;
    R.regionPeakOrder = R.dominantOrder;
    R.regionWidthSec = R.bttEndSec - R.bttStartSec;
    R = sortrows(R, {'preferredOverlapSec','regionPeakScore'}, {'descend','descend'});
    R.regionScanRank = (1:height(R)).';
end

P = build_peak_table_from_regions_local(R);
L = build_local_fft_table_from_regions_local(R);
localFftResult = struct();
localFftResult.method = 'raw_segmented_fft_surface_threshold';
localFftResult.description = ['Raw strain segmented FFT amplitude and ', ...
    'continuous intervals above the surface-amplitude threshold. ', ...
    'Detrended time amplitude is reported only as a check.'];
localFftResult.scoreTable = S;
localFftResult.params = struct('fftAmpThresholdMicrostrain', fftAmpThr, ...
    'timeAmpThresholdMicrostrain', timeAmpThr, 'lowThresholdRatio', lowRatio, ...
    'combineMode', combineMode, 'regionMergeGapSec', mergeGapSec, ...
    'minRegionDurationSec', minRegionWidthSec, ...
    'maxRegionWidthSec', maxRegionWidthSec, ...
    'preferredRegionTimeSec', preferredRegionTimeSec);
end

function R2 = merge_score_regions_local(R, mergeGapSec)
if isempty(R)
    R2 = R;
    return;
end
R = sortrows(R, 'bttStartSec', 'ascend');
R2 = R([],:);
newId = 0;
i = 1;
while i <= height(R)
    newId = newId + 1;
    t1 = R.bttEndSec(i);
    rows = i;
    j = i + 1;
    while j <= height(R)
        if R.bttStartSec(j) - t1 <= mergeGapSec
            rows(end + 1) = j; %#ok<AGROW>
            t1 = max(t1, R.bttEndSec(j));
            j = j + 1;
        else
            break;
        end
    end

    sub = R(rows, :);
    [~, idxBest] = max(sub.peakFftAmpMicrostrain);
    newRow = sub(idxBest, :);
    newRow.regionId = newId;
    newRow.bttStartSec = min(sub.bttStartSec);
    newRow.bttEndSec = max(sub.bttEndSec);
    newRow.strainStartSec = min(sub.strainStartSec);
    newRow.strainEndSec = max(sub.strainEndSec);
    newRow.peakCount = height(sub);
    newRow.regionDuration = newRow.bttEndSec - newRow.bttStartSec;
    if ismember('regionWidthSec', R.Properties.VariableNames)
        newRow.regionWidthSec = newRow.bttEndSec - newRow.bttStartSec;
    end
    if ismember('preferredOverlapSec', R.Properties.VariableNames)
        newRow.preferredOverlapSec = max(sub.preferredOverlapSec);
    end
    R2 = [R2; newRow]; %#ok<AGROW>
    i = j;
end
end

function P = build_peak_table_from_regions_local(R)
if isempty(R)
    P = empty_peak_table_local();
    return;
end
P = table(R.dominantOrder, R.regionScanRank, R.peakStrainTimeSec, ...
    R.peakBttTimeSec, R.dominantFreqHz, R.peakFftAmpMicrostrain, ...
    R.peakAmpDb, R.peakRpm, R.peakFftAmpMicrostrain, ...
    R.peakTimeAmpMicrostrain, R.regionDuration, ...
    'VariableNames', {'order','rankInOrder','strainTimeSec','bttTimeSec', ...
    'orderFreqHz','amp','ampDb','rpm','fftAmpMicrostrain', ...
    'timeAmpMicrostrain','peakWidthSec'});
end

function L = build_local_fft_table_from_regions_local(R)
if isempty(R)
    L = table(NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
        NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
        'VariableNames', {'regionId','regionStart','regionEnd', ...
        'regionPeakTime','regionPeakFreqHz','regionPeakOrder', ...
        'fftPeakAmpMicrostrain','peakTimeAmpMicrostrain'});
    return;
end
L = table(R.regionId, R.regionStart, R.regionEnd, R.regionPeakTime, ...
    R.regionPeakFreqHz, R.regionPeakOrder, R.peakFftAmpMicrostrain, ...
    R.peakTimeAmpMicrostrain, ...
    'VariableNames', {'regionId','regionStart','regionEnd', ...
    'regionPeakTime','regionPeakFreqHz','regionPeakOrder', ...
    'fftPeakAmpMicrostrain','peakTimeAmpMicrostrain'});
end

function order = nearest_finite_order_local(orderApprox, orders)
if ~isfinite(orderApprox)
    order = NaN;
    return;
end
[~, idx] = min(abs(orders - orderApprox));
order = orders(idx);
end

function d = overlap_duration_local(a, b)
d = max(0, min(a(2), b(2)) - max(a(1), b(1)));
end

function s = string_or_inf_local(x)
if isinf(x)
    s = "inf";
else
    s = string(x);
end
end

function T = empty_peak_table_local()
T = table(NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    'VariableNames', {'order','rankInOrder','strainTimeSec','bttTimeSec', ...
    'orderFreqHz','amp','ampDb','rpm','fftAmpMicrostrain','timeAmpMicrostrain','peakWidthSec'});
end

function R = empty_region_table_local()
R = table(NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), strings(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    NaN(0, 1), NaN(0, 1), NaN(0, 1), ...
    'VariableNames', {'regionId','bttStartSec','bttEndSec', ...
    'strainStartSec','strainEndSec','regionDuration','peakBttTimeSec', ...
    'peakStrainTimeSec','dominantFreqHz','dominantOrder','peakAmpDb', ...
    'peakFftAmpMicrostrain','peakTimeAmpMicrostrain','peakRpm', ...
    'syncFreqAtRegionOrderHz','fftPeakOrderApprox','regionPeakScore', ...
    'peakCount','preferredOverlapSec','ordersIncluded','regionStart', ...
    'regionEnd','regionPeakTime','regionPeakFreqHz','regionPeakOrder', ...
    'regionWidthSec','regionScanRank'});
end

function T = add_strain_labels_local(T, channelIdx, channelLabel)
T.strainChannel = repmat(channelIdx, height(T), 1);
T.strainLabel = repmat(string(channelLabel), height(T), 1);
if width(T) > 2
    T = movevars(T, {'strainChannel','strainLabel'}, 'Before', 1);
end
end

function S = summarize_channel_local(channelIdx, channelLabel, peakTable, regionTable)
if isempty(peakTable)
    maxPeakAmpDb = NaN;
    strongestOrder = NaN;
    strongestFreqHz = NaN;
    peakBttTimeSec = NaN;
else
    [maxPeakAmpDb, idx] = max(peakTable.ampDb);
    strongestOrder = peakTable.order(idx);
    strongestFreqHz = peakTable.orderFreqHz(idx);
    peakBttTimeSec = peakTable.bttTimeSec(idx);
end
S = table(channelIdx, string(channelLabel), height(peakTable), ...
    height(regionTable), maxPeakAmpDb, strongestOrder, strongestFreqHz, ...
    peakBttTimeSec, ...
    'VariableNames', {'strainChannel','strainLabel','peakCount', ...
    'regionCount','maxPeakAmpDb','strongestOrder','strongestFreqHz', ...
    'peakBttTimeSec'});
end

function plot_channel_comparison_local(channelSummary, allPeakTable, freqBand)
figure('Name', 'Step04 strain-channel comparison', ...
    'Color', 'w', 'Position', [80, 70, 1180, 680], 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
bar(ax1, channelSummary.strainChannel, channelSummary.maxPeakAmpDb, 0.65, ...
    'FaceColor', [0.20 0.50 0.85]);
grid(ax1, 'on');
box(ax1, 'on');
xticks(ax1, channelSummary.strainChannel);
xticklabels(ax1, cellstr(channelSummary.strainLabel));
ylabel(ax1, 'Max large-peak level (dB)');
title(ax1, '1) Strongest large resonance peak of each strain channel');
for i = 1:height(channelSummary)
    if isfinite(channelSummary.maxPeakAmpDb(i))
        text(ax1, channelSummary.strainChannel(i), channelSummary.maxPeakAmpDb(i), ...
            sprintf('%dX %.0f Hz', channelSummary.strongestOrder(i), ...
            channelSummary.strongestFreqHz(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Interpreter', 'none');
    end
end

ax2 = nexttile;
hold(ax2, 'on');
colors = lines(height(channelSummary));
for i = 1:height(channelSummary)
    ch = channelSummary.strainChannel(i);
    rows = allPeakTable.strainChannel == ch & ...
        allPeakTable.orderFreqHz >= freqBand(1) & ...
        allPeakTable.orderFreqHz <= freqBand(2);
    scatter(ax2, allPeakTable.bttTimeSec(rows), allPeakTable.orderFreqHz(rows), ...
        28 + max(0, allPeakTable.ampDb(rows) - 30) * 4, colors(i, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, ...
        'DisplayName', char(channelSummary.strainLabel(i)));
end
grid(ax2, 'on');
box(ax2, 'on');
xlabel(ax2, 'BTT-aligned time (s)');
ylabel(ax2, 'Order frequency (Hz)');
ylim(ax2, freqBand);
title(ax2, '2) Large peaks from each detrended strain channel');
legend(ax2, 'Location', 'bestoutside');
end

function plot_dashboard_local(strain, rawFft, offsetSec, opr, ~, ~, ~, ~, ...
    peakTable, regionTable, timeRange, freqBand)
tSignal = strain.time(:) + offsetSec;
x = strain.rawValue(:);
valid = tSignal >= timeRange(1) & tSignal <= timeRange(2) & isfinite(x);
tPlot = tSignal(valid);
xPlot = x(valid);
idxThin = decimate_index_local(numel(tPlot), 60000);

rawTimeBtt = rawFft.segmentCenterTimeSec(:).' + offsetSec;
freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);
rawAmp = rawFft.amplitudesMicrostrain(:, freqMask).';
peakMask = rawTimeBtt >= timeRange(1) & rawTimeBtt <= timeRange(2);

figure('Name', 'Step04 resonance visualization dashboard', ...
    'Color', 'w', 'Position', [50, 35, 1420, 940], 'NumberTitle', 'off');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, tPlot(idxThin), xPlot(idxThin), 'Color', [0.16 0.16 0.16], ...
    'LineWidth', 0.5);
hold(ax1, 'on');
draw_region_patches_local(ax1, regionTable);
box(ax1, 'on');
set(ax1, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8.5);
xlim(ax1, timeRange);
ylabel(ax1, 'Raw strain (\muepsilon)');
title(ax1, sprintf('1) %s raw strain time history', strain.label), ...
    'Interpreter', 'none');

ax2 = nexttile;
plot_speed_trace_local(ax2, opr, peakTable, regionTable, timeRange);

ax3 = nexttile;
imagesc(ax3, rawTimeBtt, rawFft.freqHz(freqMask), rawAmp);
axis(ax3, 'xy');
box(ax3, 'on');
set(ax3, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8.5);
xlim(ax3, timeRange);
ylim(ax3, freqBand);
colormap(ax3, high_contrast_colormap_local(256));
clim(ax3, raw_amplitude_color_limits_local(rawAmp, 99.5));
cb = colorbar(ax3);
cb.Label.String = 'FFT amplitude (\muepsilon)';
ylabel(ax3, 'Frequency (Hz)');
title(ax3, sprintf('3) Raw segmented FFT amplitude, %.0f-%.0f Hz', ...
    freqBand(1), freqBand(2)));
hold(ax3, 'on');
draw_resonance_ellipses_local(ax3, regionTable, freqBand);
plot_major_peak_points_local(ax3, peakTable, regionTable, freqBand);

ax4 = nexttile;
hold(ax4, 'on');
box(ax4, 'on');
set(ax4, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8.5);
yyaxis(ax4, 'left');
plot(ax4, rawTimeBtt(peakMask), rawFft.segmentPeakFreqHz(peakMask), ...
    'o-', 'Color', [0.05 0.22 0.65], 'MarkerFaceColor', [0.20 0.45 0.95], ...
    'MarkerSize', 3.8, 'LineWidth', 1.0);
ylabel(ax4, 'Peak frequency (Hz)');
ylim(ax4, freqBand);
draw_region_patches_local(ax4, regionTable);
yyaxis(ax4, 'right');
plot(ax4, rawTimeBtt(peakMask), rawFft.segmentPeakAmpMicrostrain(peakMask), ...
    's-', 'Color', [0.75 0.18 0.08], 'MarkerFaceColor', [0.95 0.38 0.18], ...
    'MarkerSize', 3.5, 'LineWidth', 1.0);
ylabel(ax4, 'Peak amplitude (\muepsilon)');
xlim(ax4, timeRange);
xlabel(ax4, 'BTT-aligned time (s)');
title(ax4, sprintf('4) Segment peak trajectory in %.0f-%.0f Hz raw FFT band', ...
    rawFft.peakBandHz(1), rawFft.peakBandHz(2)));
if ~isempty(regionTable)
    for i = 1:height(regionTable)
        xline(ax4, regionTable.peakBttTimeSec(i), '--', ...
            sprintf('R%d %dX', regionTable.regionId(i), regionTable.dominantOrder(i)), ...
            'Color', [0.20 0.20 0.20], 'LineWidth', 0.8, ...
            'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
    end
end

linkaxes([ax1, ax2, ax3, ax4], 'x');
end

function plot_speed_trace_local(ax, opr, peakTable, regionTable, timeRange)
valid = opr.rpmTime >= timeRange(1) & opr.rpmTime <= timeRange(2) & ...
    isfinite(opr.rpm);
t = opr.rpmTime(valid);
rpm = opr.rpm(valid);
plot(ax, t, rpm, 'k-', 'LineWidth', 1.0);
hold(ax, 'on');
if numel(rpm) > 9
    smoothPts = max(5, round(numel(rpm) / 350));
    plot(ax, t, movmean(rpm, smoothPts, 'omitnan'), ...
        'Color', [0.85 0.20 0.10], 'LineWidth', 1.3);
end
draw_region_patches_local(ax, regionTable);
plot_peak_rpm_points_local(ax, peakTable, regionTable);
grid(ax, 'on');
box(ax, 'on');
xlim(ax, timeRange);
ylabel(ax, 'RPM');
title(ax, '2) OPR speed variation and large resonance peak times');
end

function plot_peak_rpm_points_local(ax, T, R)
if isempty(T)
    return;
end
if isempty(R)
    colors = repmat([1.00 0.80 0.10], height(T), 1);
else
    colors = colors_for_peak_rows_local(T, R);
end
ampMin = min(T.ampDb);
ampMax = max(T.ampDb);
if ampMax <= ampMin
    markerSize = 28 * ones(height(T), 1);
else
    markerSize = 18 + 52 * (T.ampDb - ampMin) / (ampMax - ampMin);
end
for k = 1:height(T)
    if ~isfinite(T.rpm(k))
        continue;
    end
    scatter(ax, T.bttTimeSec(k), T.rpm(k), markerSize(k), colors(k, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7, ...
        'MarkerFaceAlpha', 0.85, 'HandleVisibility', 'off');
end
end

function plot_3d_auxiliary_local(rawFft, offsetSec, timeRange, freqBand, regionTable)
tBtt = rawFft.segmentCenterTimeSec(:).' + offsetSec;
freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);
timeMask = tBtt >= timeRange(1) & tBtt <= timeRange(2);
freq = rawFft.freqHz(freqMask);
t = tBtt(timeMask);
A = rawFft.amplitudesMicrostrain(timeMask, freqMask).';

[T, F] = meshgrid(t, freq);
figure('Name', 'Step04 auxiliary raw segmented FFT', ...
    'Color', 'w', 'Position', [90, 90, 1280, 760], 'NumberTitle', 'off');
surf(T, F, A, 'EdgeColor', 'none', 'FaceColor', 'interp');
grid on; box on;
colormap(high_contrast_colormap_local(256));
ampLimits = raw_amplitude_color_limits_local(A, 100);
clim(ampLimits);
zlim(ampLimits);
cb = colorbar;
cb.Label.String = 'FFT amplitude (\muepsilon)';
view(42, 34);
axis tight;
xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
zlabel('FFT amplitude (\muepsilon)');
title(sprintf('Auxiliary raw segmented FFT, %.0f-%.0f Hz', ...
    freqBand(1), freqBand(2)));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'TickDir', 'in', 'LineWidth', 0.8);
hold on;
if ~isempty(regionTable)
    zTop = ampLimits(2);
    for i = 1:height(regionTable)
        plot3([regionTable.bttStartSec(i), regionTable.bttEndSec(i)], ...
            [regionTable.dominantFreqHz(i), regionTable.dominantFreqHz(i)], ...
            [zTop, zTop], 'w-', 'LineWidth', 2.0);
        text(regionTable.peakBttTimeSec(i), regionTable.dominantFreqHz(i), zTop, ...
            sprintf('R%d %dX', regionTable.regionId(i), regionTable.dominantOrder(i)), ...
            'Color', 'w', 'FontWeight', 'bold');
    end
end
end

function plot_fft1_style_3d_fft_local(strain, rawFft, ~, sampleRange, numSegments, freqBandHz, zLim)
freq = rawFft.freqHz(:);
segmentTime = rawFft.segmentCenterTimeSec(:) - rawFft.segmentCenterTimeSec(1);
amplitudes = rawFft.amplitudesMicrostrain;
figure('Name', 'Step04 fft1-style 3D segmented FFT', ...
    'Color', 'w', 'NumberTitle', 'off');
surf(segmentTime(:), freq, amplitudes.');
xlim([segmentTime(1), segmentTime(end)]);
ylim(freqBandHz);
zlim(zLim);
shading interp;
colorbar;
colormap(jet);
view(30, 30);
set(gca, 'FontSize', 9);
set(gca, 'FontName', 'Times');
set(gca, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [5 10 9 6.5]);
xlabel('Raw strain time from selected start (s)');
ylabel('Frequency (Hz)');
zlabel('FFT amplitude (\muepsilon)');
title(sprintf('%s fft1-style 3D FFT, samples %d:%d, %d segments', ...
    strain.label, sampleRange(1), rawFft.sampleEnd, numSegments));
end

function plot_local_2d_fft_local(strain, offsetSec, opr, bttWindowSec, Fs, ...
    freqBandHz, orders)
tBttSignal = strain.time(:) + offsetSec;
x = strain.rawValue(:);
keep = tBttSignal >= bttWindowSec(1) & tBttSignal <= bttWindowSec(2) & isfinite(x);
if nnz(keep) < 8
    warning('Local FFT window %.2f-%.2f s has too few strain samples.', ...
        bttWindowSec(1), bttWindowSec(2));
    return;
end

tLocal = tBttSignal(keep);
xLocal = x(keep);

N = numel(xLocal);
Y = fft(xLocal);
P2 = abs(Y / N);
P1 = P2(1:floor(N / 2) + 1);
if numel(P1) > 2
    P1(2:end-1) = 2 * P1(2:end-1);
end
freq = (0:(numel(P1) - 1)).' * Fs / N;
P1(1) = 0;
freqMask = freq >= freqBandHz(1) & freq <= freqBandHz(2);

rpmLocal = interp1(opr.rpmTime, opr.rpm, tLocal, 'linear', NaN);
rpmMean = mean(rpmLocal, 'omitnan');
rotFreqMeanHz = rpmMean / 60;
orderFreqMeanHz = orders(:) * rotFreqMeanHz;

figure('Name', sprintf('Step04 local direct FFT %.1f-%.1f s', ...
    bttWindowSec(1), bttWindowSec(2)), ...
    'Color', 'w', 'Position', [80, 80, 1250, 620], 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
idxThin = decimate_index_local(numel(tLocal), 30000);
plot(ax1, tLocal(idxThin), xLocal(idxThin), 'k-', 'LineWidth', 0.8);
xlabel(ax1, 'BTT-aligned time (s)');
ylabel(ax1, 'Raw strain (\muepsilon)');
title(ax1, sprintf('%s local time signal, %.2f-%.2f s', ...
    strain.label, bttWindowSec(1), bttWindowSec(2)));

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
plot(ax2, freq(freqMask), P1(freqMask), 'k-', 'LineWidth', 1.2);
for i = 1:numel(orderFreqMeanHz)
    f = orderFreqMeanHz(i);
    if f >= freqBandHz(1) && f <= freqBandHz(2)
        xline(ax2, f, '--', sprintf('%dX', orders(i)), ...
            'Color', [0.85 0.25 0.10], 'LineWidth', 0.8, ...
            'LabelVerticalAlignment', 'bottom');
    end
end

xlim(ax2, freqBandHz);
xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'FFT amplitude (\muepsilon)');
title(ax2, sprintf('Raw direct FFT of full local segment, N=%d, df=%.3f Hz, mean speed %.1f rpm', ...
    N, Fs / N, rpmMean));
end

function draw_region_patches_local(ax, R)
if isempty(R)
    return;
end
yl = ylim(ax);
colors = region_palette_local(max(1, height(R)));
for i = 1:height(R)
    c = colors(i, :);
    patch(ax, [R.bttStartSec(i), R.bttEndSec(i), R.bttEndSec(i), R.bttStartSec(i)], ...
        [yl(1), yl(1), yl(2), yl(2)], c, ...
        'FaceAlpha', 0.055, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(ax, R.bttStartSec(i), ':', 'Color', c, ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, R.bttEndSec(i), ':', 'Color', c, ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, R.peakBttTimeSec(i), '--', sprintf('R%d', R.regionId(i)), ...
        'Color', c, 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
ylim(ax, yl);
lineObj = findobj(ax, 'Type', 'line');
if ~isempty(lineObj)
    uistack(lineObj, 'top');
end
end

function draw_resonance_ellipses_local(ax, R, freqBand)
if isempty(R)
    return;
end
colors = region_palette_local(height(R));
theta = linspace(0, 2*pi, 240);
freqSpan = diff(freqBand);
for i = 1:height(R)
    c = colors(i, :);
    cx = R.peakBttTimeSec(i);
    cy = R.dominantFreqHz(i);
    rx = max(4.0, 0.58 * (R.bttEndSec(i) - R.bttStartSec(i)));
    ry = max(18.0, 0.075 * freqSpan);
    x = cx + rx * cos(theta);
    y = cy + ry * sin(theta);
    outside = y < freqBand(1) | y > freqBand(2);
    x(outside) = NaN;
    y(outside) = NaN;

    plot(ax, x, y, 'w-', 'LineWidth', 2.6, ...
        'HandleVisibility', 'off');
    plot(ax, x, y, '--', 'Color', c, 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
    plot(ax, cx, cy, 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');

    labelY = min(freqBand(2) - 6, cy + ry + 9);
    text(ax, cx, labelY, ...
        sprintf('R%d  %dX  %.0f Hz', R.regionId(i), ...
        R.dominantOrder(i), R.dominantFreqHz(i)), ...
        'Color', c, 'FontWeight', 'bold', 'FontSize', 8, ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0], ...
        'Margin', 1, 'Interpreter', 'none');
end
end

function plot_major_peak_points_local(ax, T, R, freqBand)
if isempty(T)
    return;
end
if isempty(R)
    colors = repmat([1.00 0.80 0.10], height(T), 1);
else
    colors = colors_for_peak_rows_local(T, R);
end

ampMin = min(T.ampDb);
ampMax = max(T.ampDb);
if ampMax <= ampMin
    markerSize = 42 * ones(height(T), 1);
else
    markerSize = 30 + 70 * (T.ampDb - ampMin) / (ampMax - ampMin);
end

for k = 1:height(T)
    y = T.orderFreqHz(k);
    if y < freqBand(1) || y > freqBand(2)
        continue;
    end
    scatter(ax, T.bttTimeSec(k), y, markerSize(k), colors(k, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.82, 'HandleVisibility', 'off');

    labelOffset = 7 + 5 * mod(k, 2);
    labelY = min(freqBand(2) - 5, y + labelOffset);
    text(ax, T.bttTimeSec(k), labelY, sprintf('%dX', T.order(k)), ...
        'Color', colors(k, :), 'FontWeight', 'bold', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0], ...
        'Margin', 0.5, 'Interpreter', 'none');
end
end

function colors = colors_for_peak_rows_local(T, R)
colors = zeros(height(T), 3);
regionColors = region_palette_local(height(R));
for k = 1:height(T)
    regionIdx = find(T.bttTimeSec(k) >= R.bttStartSec & ...
        T.bttTimeSec(k) <= R.bttEndSec, 1, 'first');
    if isempty(regionIdx)
        colors(k, :) = [1.00 1.00 1.00];
    else
        colors(k, :) = regionColors(regionIdx, :);
    end
end
end

function colors = region_palette_local(n)
base = [ ...
    1.000 0.720 0.080
    0.250 0.900 0.450
    1.000 0.280 0.520
    0.980 0.520 0.100
    0.680 0.920 0.220
    0.950 0.720 0.720];
idx = mod((1:n) - 1, size(base, 1)) + 1;
colors = base(idx, :);
end

function colorLimits = raw_amplitude_color_limits_local(x, highPct)
x = x(isfinite(x));
if isempty(x)
    colorLimits = [0, 1];
    return;
end
lo = 0;
hi = percentile_local(x, highPct);
if ~isfinite(hi) || hi <= lo
    hi = max(x);
end
if ~isfinite(hi) || hi <= lo
    hi = 1;
end
colorLimits = [lo, hi];
end

function p = percentile_local(x, q)
x = sort(x(isfinite(x)));
if isempty(x)
    p = NaN;
    return;
end
idx = 1 + (numel(x) - 1) * q / 100;
i0 = floor(idx);
i1 = ceil(idx);
if i0 == i1
    p = x(i0);
else
    p = x(i0) + (idx - i0) * (x(i1) - x(i0));
end
end

function cmap = high_contrast_colormap_local(n)
base = [ ...
    0.000 0.000 0.000
    0.020 0.040 0.160
    0.120 0.030 0.320
    0.380 0.030 0.360
    0.700 0.060 0.180
    0.950 0.300 0.050
    1.000 0.720 0.120
    1.000 1.000 0.780];
xBase = linspace(0, 1, size(base, 1));
xq = linspace(0, 1, n);
cmap = interp1(xBase, base, xq, 'pchip');
cmap = min(max(cmap, 0), 1);
end

function idx = decimate_index_local(n, maxN)
if n <= maxN
    idx = 1:n;
else
    idx = unique(round(linspace(1, n, maxN)));
end
end

function w = hann_local(n)
if n <= 1
    w = 1;
else
    k = (0:n-1).';
    w = 0.5 - 0.5 * cos(2*pi*k/(n-1));
end
end
