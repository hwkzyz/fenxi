%% Step04: resonance detection and visualization dashboard for 20251222
% Unified Step04 workflow:
%   1. aligned strain time response and target-blade BTT vibration
%   2. OPR speed variation
%   3. raw segmented FFT with direct 100-microstrain resonance regions
%   4. resonance-region summary
%
% Figures are opened directly. The selected regions are also saved using the
% legacy Step04 filenames so downstream scripts keep reading the same files.

clear; clc; close all;

%% 1. Simple settings
thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20251222();
legacyDir = packageCfg.paths.preparation;
addpath(legacyDir);

cfg = Get_20251222_BTT_Config();
caseName = cfg.dynamic_cases{1};
caseOutputDir = fullfile(cfg.output_root, caseName);
outDir = packageCfg.paths.prepared;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sensorIds = cfg.step3_default_btt_sensor_ids;
targetBladeId = 1;              % Must match the blade with the strain gauge.
baselineTimeSec = 3.0;
detrendMethod = 'linear';       % 'linear', 'movingmedian', or 'none' for strain time plot.

orders = cfg.step3_align_order_candidates(:).';
resonanceFftAmpThresholdMicrostrain = 100.0;
resonanceTimeAmpThresholdMicrostrain = 100.0;
resonanceLowThresholdRatio = 0.55;
resonanceCombineMode = "both";
regionMergeGapSec = 0.00;
minRegionDurationSec = 0.00;
maxRegionWidthSec = inf;
preferredRegionTimeSec = [50, 55];  % Later identification uses this constant-speed segment.
rawFftNumSegments = 200;
rawFftPeakBandHz = [20, 800];

dashboardFreqBandHz = cfg.step3_align_freq_plot_hz;
make3DFigure = true;

%% 2. Load Step03 alignment and BTT displacement
alignFile = fullfile(caseOutputDir, 'Step3_Spectrum_RPM_Alignment_20251222.mat');
if exist(alignFile, 'file') ~= 2
    error('Run Step03 first. Missing file: %s', alignFile);
end
A = load(alignFile, 'result');
alignResult = A.result;
check_stft_preprocess_local(alignResult, detrendMethod, alignFile);

stft.freqHz = alignResult.stft_freq_hz(:);
stft.timeSec = alignResult.stft_time_s(:).';
stft.amp = alignResult.stft_amp;
tBtt = stft.timeSec - alignResult.best_tau_sec;

rpm.timeSec = alignResult.rpm_time_s(:);
rpm.value = alignResult.rpm_values(:);
validTimeRange = [max(min(tBtt), min(rpm.timeSec)), ...
    min(max(tBtt), max(rpm.timeSec))];

step3ResultFile = fullfile(caseOutputDir, 'Step3_Result_20251222.mat');
[strainTimeBtt, strainValue, strainSource, strainRaw] = load_aligned_strain_time_series_local( ...
    step3ResultFile, detrendMethod, alignResult.best_tau_sec);
[bttTime, bttVib] = load_btt_vibration_local(caseOutputDir, sensorIds, ...
    targetBladeId, baselineTimeSec);

%% 3. Direct resonance detection from raw segmented FFT and detrended strain
rawFft = raw_segmented_fft_local(strainRaw, validTimeRange, ...
    rawFftNumSegments, rawFftPeakBandHz);
resonanceParams = struct();
resonanceParams.freqBandHz = dashboardFreqBandHz;
resonanceParams.peakBandHz = rawFftPeakBandHz;
resonanceParams.fftAmpThresholdMicrostrain = resonanceFftAmpThresholdMicrostrain;
resonanceParams.timeAmpThresholdMicrostrain = resonanceTimeAmpThresholdMicrostrain;
resonanceParams.lowThresholdRatio = resonanceLowThresholdRatio;
resonanceParams.combineMode = resonanceCombineMode;
resonanceParams.regionMergeGapSec = regionMergeGapSec;
resonanceParams.minRegionDurationSec = minRegionDurationSec;
resonanceParams.maxRegionWidthSec = maxRegionWidthSec;
resonanceParams.preferredRegionTimeSec = preferredRegionTimeSec;

[scoreTable, peakTable, regionTable, localFftTable, localFftResult] = ...
    detect_resonance_regions_from_raw_fft_local(rawFft, rpm, orders, ...
    resonanceParams);

%% 4. Print summary
fprintf('\n=== 20251222 Step04 resonance visualization dashboard ===\n');
fprintf('Case: %s\n', caseName);
fprintf('Target blade for BTT/strain time comparison: blade %d\n', targetBladeId);
fprintf('Strain time source: %s\n', string_or_empty_local(strainSource, 'strain_file'));
fprintf('Strain time detrend method: %s\n', detrendMethod);
fprintf('best_tau_sec = %+0.3f s, BTT_time = strain_STFT_time - best_tau_sec\n', ...
    alignResult.best_tau_sec);
fprintf('Orders: %s\n', mat2str(orders));
fprintf(['Resonance rule: raw segmented FFT surface max >= %.1f microstrain; ', ...
    'continuous above-threshold intervals are resonance regions.\n'], ...
    resonanceFftAmpThresholdMicrostrain);
fprintf('Preferred downstream segment priority: %.2f-%.2f s\n\n', preferredRegionTimeSec);

fprintf('Region summary:\n');
disp(regionTable);
fprintf('\nPeak table used for regions:\n');
disp(peakTable);

%% 5. Main dashboard
plot_dashboard_local(strainTimeBtt, strainValue, bttTime, bttVib, targetBladeId, ...
    rpm, rawFft, peakTable, ...
    regionTable, validTimeRange, dashboardFreqBandHz);

%% 6. Optional 3D auxiliary view
if make3DFigure
    plot_3d_auxiliary_local(rawFft, validTimeRange, dashboardFreqBandHz, regionTable);
end

matFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.mat');
csvFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.csv');
scoreCsvFile = fullfile(outDir, 'Step04_DirectStrainResonance_Score_20251222.csv');
localFftCsvFile = fullfile(outDir, 'Step04_LocalStrainFFT_20251222.csv');
save(matFile, 'scoreTable', 'peakTable', 'regionTable', 'localFftTable', ...
    'localFftResult', 'rawFft', 'strainSource', 'resonanceParams', ...
    'alignResult', 'caseName', '-v7');
writetable(regionTable, csvFile);
writetable(scoreTable, scoreCsvFile);
writetable(localFftTable, localFftCsvFile);

fprintf('\nSaved Step04 resonance outputs:\n  %s\n  %s\n  %s\n  %s\n', ...
    matFile, csvFile, scoreCsvFile, localFftCsvFile);
fprintf('Done. Figures are open in MATLAB.\n');

%% Local functions
function check_stft_preprocess_local(alignResult, expectedMethod, alignFile)
if ~isfield(alignResult, 'strain_preprocess_method')
    error(['The Step03 STFT result has no strain_preprocess_method field.\n', ...
        'Run Step03 again so STFT identification uses detrended strain:\n  %s'], alignFile);
end
if ~strcmpi(string(alignResult.strain_preprocess_method), string(expectedMethod))
    error(['The Step03 STFT was built with strain_preprocess_method = "%s", ', ...
    'but Step04 expects "%s".\nRun Step03 again or set the same method in both files.'], ...
        string(alignResult.strain_preprocess_method), string(expectedMethod));
end
end

function [strainTimeBtt, strainValue, strainSource, strainRaw] = load_aligned_strain_time_series_local( ...
    step3ResultFile, detrendMethod, bestTauSec)
strainTimeBtt = [];
strainValue = [];
strainSource = struct('step3_file', step3ResultFile, 'strain_file', '', ...
    'strain_time_offset_total_sec', NaN, 'best_tau_sec', bestTauSec, ...
    'time_coordinate', 'BTT = raw strain + offset - best tau');
strainRaw = struct('timeBtt', [], 'rawValue', [], 'detrendedValue', [], ...
    'sampleRateHz', NaN, 'file', '', 'timeOffsetSec', NaN);
if exist(step3ResultFile, 'file') ~= 2
    warning('Missing Step3 result for strain time-domain plot: %s', step3ResultFile);
    return;
end
D = load(step3ResultFile, 'result');
if ~isfield(D, 'result') || ~isfield(D.result, 'strain_file') || ...
        exist(D.result.strain_file, 'file') ~= 2
    warning('Invalid strain source recorded in Step3 result: %s', step3ResultFile);
    return;
end
S = load(D.result.strain_file, 'Datas');
if ~isfield(S, 'Datas') || size(S.Datas, 2) < 2
    warning('Invalid strain Datas in file: %s', D.result.strain_file);
    return;
end
timeOffsetSec = D.result.strain_time_offset_total_sec;
% Step03 defines best_tau_sec by BTT_time = strain_time - best_tau_sec.
% Apply that same conversion to the raw strain trace.  Previously only the
% STFT axis was shifted, while the raw segmented FFT retained strain-clock
% time and was incorrectly labelled as BTT time.
strainTimeBtt = S.Datas(:, 1) + timeOffsetSec - bestTauSec;
strainValue = preprocess_strain_local(S.Datas(:, 2), strainTimeBtt, detrendMethod);
valid = isfinite(strainTimeBtt) & isfinite(strainValue);
rawValue = S.Datas(:, 2);
strainTimeBtt = strainTimeBtt(valid);
strainValue = strainValue(valid);
rawValue = rawValue(valid);
[strainTimeBtt, orderIdx] = sort(strainTimeBtt);
strainValue = strainValue(orderIdx);
rawValue = rawValue(orderIdx);
strainSource.strain_file = D.result.strain_file;
strainSource.strain_time_offset_total_sec = timeOffsetSec;
strainRaw.timeBtt = strainTimeBtt;
strainRaw.rawValue = rawValue;
strainRaw.detrendedValue = strainValue;
strainRaw.sampleRateHz = estimate_sample_rate_local(strainTimeBtt);
strainRaw.file = D.result.strain_file;
strainRaw.timeOffsetSec = timeOffsetSec;
end

function y = preprocess_strain_local(x, t, detrendMethod)
x = x(:);
t = t(:);
finiteX = isfinite(x);
if any(finiteX)
    x(~finiteX) = median(x(finiteX), 'omitnan');
else
    x(:) = 0;
end
switch lower(string(detrendMethod))
    case "none"
        y = x - median(x, 'omitnan');
    case "linear"
        y = detrend(x, 'linear');
    case "movingmedian"
        Fs = estimate_sample_rate_local(t);
        win = max(5, round(2.0 * Fs));
        if mod(win, 2) == 0
            win = win + 1;
        end
        y = x - movmedian(x, win, 'omitnan');
    otherwise
        error('Unknown detrendMethod: %s', detrendMethod);
end
y = y - median(y, 'omitnan');
end

function [allTime, allVib] = load_btt_vibration_local(caseOutputDir, sensorIds, targetBladeId, baselineTimeSec)
allTime = [];
allVib = [];
for sid = sensorIds
    vibFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d_vib_final.mat', sid));
    if exist(vibFile, 'file') ~= 2
        warning('Missing BTT displacement file: %s', vibFile);
        continue;
    end
    S = load(vibFile, 'jilublade');
    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 6
        warning('Invalid BTT displacement file: %s', vibFile);
        continue;
    end
    for bid = targetBladeId
        bladeMask = S.jilublade(:, 4) == bid;
        t = S.jilublade(bladeMask, 3);
        y = S.jilublade(bladeMask, 6);
        valid = isfinite(t) & isfinite(y);
        t = t(valid);
        y = y(valid);
        if numel(t) < 10
            continue;
        end
        [t, orderIdx] = sort(t);
        y = y(orderIdx);
        dtBlade = median(diff(t), 'omitnan');
        if ~isfinite(dtBlade) || dtBlade <= 0
            dtBlade = 0.05;
        end
        baselinePts = max(5, round(baselineTimeSec / dtBlade));
        if mod(baselinePts, 2) == 0
            baselinePts = baselinePts + 1;
        end
        yVib = y - movmedian(y, baselinePts, 'omitnan');
        allTime = [allTime; t(:)]; %#ok<AGROW>
        allVib = [allVib; yVib(:)]; %#ok<AGROW>
    end
end
if isempty(allTime)
    error('No valid BTT displacement points were loaded.');
end
[allTime, idx] = sort(allTime);
allVib = allVib(idx);
end

function rawFft = raw_segmented_fft_local(strainRaw, timeRange, numSegments, peakBandHz)
tRaw = strainRaw.timeBtt(:);
xRaw = strainRaw.rawValue(:);
xDetrended = strainRaw.detrendedValue(:);
Fs = strainRaw.sampleRateHz;
validRange = isfinite(tRaw) & isfinite(xRaw) & tRaw >= timeRange(1) & tRaw <= timeRange(2);
idxAll = find(validRange);
if numel(idxAll) < 16
    error('Not enough raw strain samples inside BTT time range %.3f-%.3f s.', ...
        timeRange(1), timeRange(2));
end
i1 = idxAll(1);
i2 = idxAll(end);
x = xRaw(i1:i2);
xd = xDetrended(i1:i2);
t = tRaw(i1:i2);
N = numel(x);
numSegments = min(numSegments, floor(N / 16));
segmentLength = floor(N / numSegments);
if segmentLength < 2
    error('Raw FFT segment length is too short: %d samples.', segmentLength);
end

Nused = segmentLength * numSegments;
x = x(1:Nused);
xd = xd(1:Nused);
t = t(1:Nused);
freq = (0:floor(segmentLength / 2)).' * (Fs / segmentLength);
amplitudes = zeros(numSegments, numel(freq));
segmentCenterTimeSec = zeros(numSegments, 1);
segmentPeakFreqHz = zeros(numSegments, 1);
segmentPeakAmpMicrostrain = zeros(numSegments, 1);
segmentDetrendedMaxAbsMicrostrain = zeros(numSegments, 1);
segmentDetrendedHalfPeakToPeakMicrostrain = zeros(numSegments, 1);
segmentDetrendedP95AbsMicrostrain = zeros(numSegments, 1);
peakMask = freq >= peakBandHz(1) & freq <= peakBandHz(2);
freqForPeak = freq(peakMask);

for i = 1:numSegments
    localIdx = (i - 1) * segmentLength + (1:segmentLength);
    segmentData = x(localIdx);
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
    segmentDetrended = xd(localIdx);
    halfPeakToPeak = 0.5 * (max(segmentDetrended, [], 'omitnan') - ...
        min(segmentDetrended, [], 'omitnan'));
    absDetrended = abs(segmentDetrended);
    segmentDetrendedHalfPeakToPeakMicrostrain(i) = halfPeakToPeak;
    segmentDetrendedMaxAbsMicrostrain(i) = halfPeakToPeak;
    segmentDetrendedP95AbsMicrostrain(i) = percentile_local(absDetrended, 95);
end

rawFft.sampleStart = i1;
rawFft.sampleEnd = i1 + Nused - 1;
rawFft.segmentLength = segmentLength;
rawFft.numSegments = numSegments;
rawFft.Nused = Nused;
rawFft.sampleRateHz = Fs;
rawFft.freqHz = freq;
rawFft.segmentCenterTimeSec = segmentCenterTimeSec;
rawFft.amplitudesMicrostrain = amplitudes;
rawFft.peakBandHz = peakBandHz;
rawFft.segmentPeakFreqHz = segmentPeakFreqHz;
rawFft.segmentPeakAmpMicrostrain = segmentPeakAmpMicrostrain;
rawFft.segmentDetrendedMaxAbsMicrostrain = segmentDetrendedMaxAbsMicrostrain;
rawFft.segmentDetrendedHalfPeakToPeakMicrostrain = segmentDetrendedHalfPeakToPeakMicrostrain;
rawFft.segmentDetrendedP95AbsMicrostrain = segmentDetrendedP95AbsMicrostrain;
end

function [scoreTable, peakTable, regionTable, localFftTable, localFftResult] = ...
    detect_resonance_regions_from_raw_fft_local(rawFft, rpm, orders, P)
t = rawFft.segmentCenterTimeSec(:);
fftAmp = rawFft.segmentPeakAmpMicrostrain(:);
if isfield(rawFft, 'segmentDetrendedHalfPeakToPeakMicrostrain')
    timeAmp = rawFft.segmentDetrendedHalfPeakToPeakMicrostrain(:);
else
    timeAmp = rawFft.segmentDetrendedMaxAbsMicrostrain(:);
end
timeP95 = rawFft.segmentDetrendedP95AbsMicrostrain(:);
freq = rawFft.segmentPeakFreqHz(:);
rpmAt = interp1(rpm.timeSec, rpm.value, t, 'linear', NaN);
valid = isfinite(t) & isfinite(fftAmp) & isfinite(timeAmp) & ...
    isfinite(freq) & isfinite(rpmAt);

surfaceAmp = fftAmp;
timeAmpOut = timeAmp;
highMask = valid & surfaceAmp >= P.fftAmpThresholdMicrostrain;
lowMask = highMask;
score = surfaceAmp;

shaftFreqHz = rpmAt / 60;
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

scoreTable = table(t, rpmAt, freq, surfaceAmp, timeAmpOut, timeP95, ...
    orderApprox, nearestOrder, syncFreqAtOrderHz, score, highMask, lowMask, ...
    'VariableNames', {'bttTimeSec','rpm','peakFreqHz', ...
    'peakFftAmpMicrostrain','peakTimeAmpMicrostrain', ...
    'peakTimeP95AbsMicrostrain','fftPeakOrderApprox','nearestOrder', ...
    'syncFreqAtNearestOrderHz','resonanceScore','isHighMask','isLowMask'});

regionTable = build_raw_region_table_local(scoreTable, rawFft, orders, P);
peakTable = build_peak_table_from_regions_local(regionTable);
localFftTable = build_local_fft_table_from_regions_local(regionTable);
localFftResult = struct();
localFftResult.method = 'raw_segmented_fft_surface_threshold';
localFftResult.description = ['Raw strain segmented FFT amplitude and ', ...
    'continuous intervals above the surface-amplitude threshold. ', ...
    'Detrended half peak-to-peak amplitude is reported only as a check.'];
localFftResult.scoreTable = scoreTable;
localFftResult.params = P;
end

function R = build_raw_region_table_local(S, rawFft, orders, P)
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

    dt = median(diff(S.bttTimeSec(isfinite(S.bttTimeSec))), 'omitnan');
    if ~isfinite(dt) || dt <= 0
        dt = rawFft.segmentLength / rawFft.sampleRateHz;
    end
    t0 = S.bttTimeSec(i0) - 0.5 * dt;
    t1 = S.bttTimeSec(i1) + 0.5 * dt;
    if t1 - t0 < P.minRegionDurationSec
        continue;
    end
    peakOrder = S.nearestOrder(ip);
    if ~isfinite(peakOrder)
        peakOrder = nearest_finite_order_local(S.fftPeakOrderApprox(ip), orders);
    end
    peakAmpDb = 20 * log10(max(S.peakFftAmpMicrostrain(ip), eps));
    idx = find(S.bttTimeSec >= t0 & S.bttTimeSec <= t1);
    rows = [rows; k, t0, t1, t1 - t0, S.bttTimeSec(ip), ...
        S.peakFreqHz(ip), peakOrder, peakAmpDb, ...
        S.peakFftAmpMicrostrain(ip), S.peakTimeAmpMicrostrain(ip), ...
        S.rpm(ip), S.syncFreqAtNearestOrderHz(ip), ...
        S.fftPeakOrderApprox(ip), S.resonanceScore(ip), ...
        numel(idx), double(overlap_duration_local([t0, t1], P.preferredRegionTimeSec))]; %#ok<AGROW>
end

if isempty(rows)
    R = empty_region_table_local();
    return;
end
R = array2table(rows, 'VariableNames', {'regionId','bttStartSec','bttEndSec', ...
    'regionDuration','peakBttTimeSec','dominantFreqHz','dominantOrder', ...
    'peakAmpDb','peakFftAmpMicrostrain','peakTimeAmpMicrostrain', ...
    'peakRpm','syncFreqAtRegionOrderHz','fftPeakOrderApprox', ...
    'regionPeakScore','peakCount','preferredOverlapSec'});
R = merge_raw_regions_local(R, P.regionMergeGapSec, P.maxRegionWidthSec);
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
R = sortrows(R, {'preferredOverlapSec','regionPeakScore'}, {'descend','descend'});
R.regionScanRank = (1:height(R)).';
R.regionId = R.regionScanRank;
end

function R2 = merge_raw_regions_local(R, mergeGapSec, maxRegionWidthSec)
if isempty(R)
    R2 = R;
    return;
end
R = sortrows(R, 'bttStartSec', 'ascend');
R2 = R([],:);
newId = 0;
i = 1;
while i <= height(R)
    rows = i;
    tEnd = R.bttEndSec(i);
    j = i + 1;
    while j <= height(R)
        mergedStart = min(R.bttStartSec(rows(1)), R.bttStartSec(j));
        mergedEnd = max(tEnd, R.bttEndSec(j));
        if R.bttStartSec(j) - tEnd <= mergeGapSec && ...
                mergedEnd - mergedStart <= maxRegionWidthSec
            rows(end + 1) = j; %#ok<AGROW>
            tEnd = mergedEnd;
            j = j + 1;
        else
            break;
        end
    end
    sub = R(rows, :);
    [~, idxBest] = max(sub.regionPeakScore);
    newId = newId + 1;
    newRow = sub(idxBest, :);
    newRow.regionId = newId;
    newRow.bttStartSec = min(sub.bttStartSec);
    newRow.bttEndSec = max(sub.bttEndSec);
    newRow.regionDuration = newRow.bttEndSec - newRow.bttStartSec;
    newRow.peakCount = sum(sub.peakCount);
    newRow.preferredOverlapSec = max(sub.preferredOverlapSec);
    R2 = [R2; newRow]; %#ok<AGROW>
    i = j;
end
end

function T = build_peak_table_from_regions_local(R)
if isempty(R)
    T = empty_peak_table_local();
    return;
end
T = table(R.dominantOrder, R.regionScanRank, R.peakBttTimeSec, ...
    R.dominantFreqHz, R.peakFftAmpMicrostrain, R.peakAmpDb, R.peakRpm, ...
    R.peakTimeAmpMicrostrain, R.regionId, ...
    'VariableNames', {'order','rankInOrder','bttTimeSec','orderFreqHz', ...
    'amp','ampDb','rpm','timeAmpMicrostrain','regionId'});
end

function T = build_local_fft_table_from_regions_local(R)
if isempty(R)
    T = table();
    return;
end
T = table(R.regionId, R.regionStart, R.regionEnd, R.regionPeakTime, ...
    R.regionPeakOrder, R.regionPeakFreqHz, R.regionPeakFreqHz, ...
    R.peakFftAmpMicrostrain, R.fftPeakOrderApprox, R.syncFreqAtRegionOrderHz, ...
    R.peakTimeAmpMicrostrain, ...
    'VariableNames', {'regionId','regionStart','regionEnd','regionPeakTime', ...
    'regionPeakOrder','regionPeakFreqHz','fftPeakFreqHz', ...
    'fftPeakAmpMicrostrain','fftPeakOrderApprox','syncFreqAtRegionOrderHz', ...
    'peakTimeAmpMicrostrain'});
end

function order = nearest_finite_order_local(orderApprox, orders)
if ~isfinite(orderApprox)
    order = NaN;
    return;
end
[~, idx] = min(abs(orders - orderApprox));
order = orders(idx);
end

function dt = overlap_duration_local(a, b)
if isempty(b) || numel(b) ~= 2 || any(~isfinite(b))
    dt = 0;
    return;
end
dt = max(0, min(a(2), b(2)) - max(a(1), b(1)));
end

function T = empty_region_table_local()
T = table([], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
    strings(0, 1), [], [], [], [], [], [], ...
    'VariableNames', {'regionId','bttStartSec','bttEndSec','regionDuration', ...
    'peakBttTimeSec','dominantFreqHz','dominantOrder','peakAmpDb', ...
    'peakFftAmpMicrostrain','peakTimeAmpMicrostrain','peakRpm', ...
    'syncFreqAtRegionOrderHz','fftPeakOrderApprox','regionPeakScore', ...
    'peakCount','preferredOverlapSec','ordersIncluded','regionStart', ...
    'regionEnd','regionPeakTime','regionPeakFreqHz','regionPeakOrder', ...
    'regionScanRank'});
end

function T = empty_peak_table_local()
T = table([], [], [], [], [], [], [], [], [], ...
    'VariableNames', {'order','rankInOrder','bttTimeSec','orderFreqHz', ...
    'amp','ampDb','rpm','timeAmpMicrostrain','regionId'});
end

function plot_dashboard_local(strainTimeBtt, strainValue, bttTime, bttVib, targetBladeId, ...
    rpm, rawFft, peakTable, regionTable, timeRange, freqBand)
rawTime = rawFft.segmentCenterTimeSec(:).';
freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);
rawAmp = rawFft.amplitudesMicrostrain(:, freqMask).';

figure('Name', '20251222 Step04 resonance visualization dashboard', ...
    'Color', 'w', 'Position', [50, 35, 1420, 940], 'NumberTitle', 'off');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot_target_blade_btt_and_strain_local(ax1, strainTimeBtt, strainValue, ...
    bttTime, bttVib, targetBladeId, regionTable, timeRange);

ax2 = nexttile;
plot_speed_trace_local(ax2, rpm, peakTable, regionTable, timeRange);

ax3 = nexttile;
imagesc(ax3, rawTime, rawFft.freqHz(freqMask), rawAmp);
axis(ax3, 'xy'); grid(ax3, 'on'); box(ax3, 'on');
xlim(ax3, timeRange); ylim(ax3, freqBand);
colormap(ax3, high_contrast_colormap_local(256));
clim(ax3, raw_amplitude_color_limits_local(rawAmp, 99.5));
cb = colorbar(ax3); cb.Label.String = 'Amplitude (microstrain)';
ylabel(ax3, 'Frequency (Hz)');
title(ax3, sprintf('3) Raw segmented FFT %.0f-%.0f Hz with direct 100-microstrain regions', ...
    freqBand(1), freqBand(2)));
hold(ax3, 'on');
draw_resonance_ellipses_local(ax3, regionTable, freqBand);
plot_major_peak_points_local(ax3, peakTable, regionTable, freqBand);

ax4 = nexttile;
if ~isempty(regionTable)
    bh = bar(ax4, regionTable.regionId, regionTable.peakAmpDb, 0.65, 'FaceColor', 'flat');
    bh.CData = region_palette_local(height(regionTable));
    ylim(ax4, [0, max(regionTable.peakAmpDb) * 1.22]);
    for i = 1:height(regionTable)
        txt = sprintf('R%d: %dX\n%.1f-%.1fs', regionTable.regionId(i), ...
            regionTable.dominantOrder(i), regionTable.bttStartSec(i), regionTable.bttEndSec(i));
        text(ax4, regionTable.regionId(i), regionTable.peakAmpDb(i), txt, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Interpreter', 'none');
    end
end
grid(ax4, 'on'); box(ax4, 'on');
xlabel(ax4, 'Region ID'); ylabel(ax4, 'Peak raw FFT amplitude (dB)');
title(ax4, '4) Region summary: dominant order and peak level');
linkaxes([ax1, ax2, ax3], 'x');
end

function plot_target_blade_btt_and_strain_local(ax, strainTimeBtt, strainValue, ...
    bttTime, bttVib, targetBladeId, regionTable, timeRange)
validBtt = bttTime >= timeRange(1) & bttTime <= timeRange(2) & isfinite(bttVib);
tBtt = bttTime(validBtt);
yBtt = bttVib(validBtt);
idxThin = decimate_index_local(numel(tBtt), 60000);
if isempty(tBtt)
    warning('No BTT samples for target blade %d inside display range.', targetBladeId);
    tBtt = NaN;
    yBtt = NaN;
    idxThin = 1;
end

yyaxis(ax, 'left');
plot(ax, tBtt(idxThin), yBtt(idxThin), '.', ...
    'Color', [0.18 0.18 0.18 0.35], 'MarkerSize', 3);
hold(ax, 'on');
bttEnv = moving_rms_local(yBtt, 0.15, estimate_sample_rate_local(tBtt));
idxEnv = decimate_index_local(numel(tBtt), 12000);
plot(ax, tBtt(idxEnv), bttEnv(idxEnv), 'r-', 'LineWidth', 1.0);
plot(ax, tBtt(idxEnv), -bttEnv(idxEnv), 'r-', 'LineWidth', 1.0);
ylabel(ax, sprintf('Blade %d BTT vib. (mm)', targetBladeId));
ax.YColor = [0.10 0.10 0.10];

yyaxis(ax, 'right');
validStrain = strainTimeBtt >= timeRange(1) & strainTimeBtt <= timeRange(2) & ...
    isfinite(strainValue);
tStrain = strainTimeBtt(validStrain);
yStrain = strainValue(validStrain);
if ~isempty(tStrain)
    yStrain = yStrain - median(yStrain, 'omitnan');
    strainScale = percentile_local(abs(yStrain), 98);
    if ~isfinite(strainScale) || strainScale <= 0
        strainScale = max(abs(yStrain), [], 'omitnan');
    end
    if ~isfinite(strainScale) || strainScale <= 0
        strainScale = 1;
    end
    yStrainNorm = yStrain / strainScale;
    idxStrain = decimate_index_local(numel(tStrain), 50000);
    plot(ax, tStrain(idxStrain), yStrainNorm(idxStrain), ...
        'Color', [0.05 0.35 0.85 0.55], 'LineWidth', 0.6);
    strainEnv = moving_rms_local(yStrainNorm, 0.15, estimate_sample_rate_local(tStrain));
    idxStrainEnv = decimate_index_local(numel(tStrain), 12000);
    plot(ax, tStrain(idxStrainEnv), strainEnv(idxStrainEnv), ...
        'Color', [0.05 0.35 0.85], 'LineWidth', 1.0);
end
ylabel(ax, 'Strain norm.');
ax.YColor = [0.05 0.35 0.85];

yyaxis(ax, 'left');
draw_region_patches_local(ax, regionTable);
grid(ax, 'on'); box(ax, 'on'); xlim(ax, timeRange);
title(ax, sprintf('1) Strain response and BTT vibration of target blade %d', targetBladeId));
end

function plot_speed_trace_local(ax, rpm, peakTable, regionTable, timeRange)
valid = rpm.timeSec >= timeRange(1) & rpm.timeSec <= timeRange(2) & isfinite(rpm.value);
t = rpm.timeSec(valid);
y = rpm.value(valid);
plot(ax, t, y, 'k-', 'LineWidth', 1.0);
hold(ax, 'on');
if numel(y) > 9
    smoothPts = max(5, round(numel(y) / 350));
    plot(ax, t, movmean(y, smoothPts, 'omitnan'), 'Color', [0.85 0.20 0.10], 'LineWidth', 1.3);
end
draw_region_patches_local(ax, regionTable);
plot_peak_rpm_points_local(ax, peakTable, regionTable);
grid(ax, 'on'); box(ax, 'on'); xlim(ax, timeRange);
ylabel(ax, 'RPM');
title(ax, '2) OPR speed variation and selected resonance peak times');
end

function plot_peak_rpm_points_local(ax, T, R)
if isempty(T)
    return;
end
colors = colors_for_peak_rows_local(T, R);
sizes = scale_marker_sizes_local(T.ampDb, 18, 70);
for k = 1:height(T)
    if ~isfinite(T.rpm(k))
        continue;
    end
    scatter(ax, T.bttTimeSec(k), T.rpm(k), sizes(k), colors(k, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7, ...
        'MarkerFaceAlpha', 0.85, 'HandleVisibility', 'off');
end
end

function plot_major_peak_points_local(ax, T, R, freqBand)
if isempty(T)
    return;
end
colors = colors_for_peak_rows_local(T, R);
sizes = scale_marker_sizes_local(T.ampDb, 30, 100);
for k = 1:height(T)
    y = T.orderFreqHz(k);
    if y < freqBand(1) || y > freqBand(2)
        continue;
    end
    scatter(ax, T.bttTimeSec(k), y, sizes(k), colors(k, :), ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.82, 'HandleVisibility', 'off');
    text(ax, T.bttTimeSec(k), min(freqBand(2)-5, y+7), sprintf('%dX', T.order(k)), ...
        'Color', colors(k, :), 'FontWeight', 'bold', 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0], ...
        'Margin', 0.5, 'Interpreter', 'none');
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
    rx = max(2.0, 0.58 * (R.bttEndSec(i) - R.bttStartSec(i) + 1));
    ry = max(12.0, 0.070 * freqSpan);
    x = cx + rx * cos(theta);
    y = cy + ry * sin(theta);
    outside = y < freqBand(1) | y > freqBand(2);
    x(outside) = NaN; y(outside) = NaN;
    plot(ax, x, y, 'w-', 'LineWidth', 2.4, 'HandleVisibility', 'off');
    plot(ax, x, y, '--', 'Color', c, 'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax, cx, cy, 'o', 'MarkerSize', 6, 'MarkerFaceColor', c, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    text(ax, cx, min(freqBand(2)-6, cy+ry+8), ...
        sprintf('R%d  %dX  %.0f Hz', R.regionId(i), R.dominantOrder(i), R.dominantFreqHz(i)), ...
        'Color', c, 'FontWeight', 'bold', 'FontSize', 8, ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0], ...
        'Margin', 1, 'Interpreter', 'none');
end
end

function draw_region_patches_local(ax, R)
if isempty(R)
    return;
end
yl = ylim(ax);
colors = region_palette_local(height(R));
for i = 1:height(R)
    c = colors(i, :);
    patch(ax, [R.bttStartSec(i), R.bttEndSec(i), R.bttEndSec(i), R.bttStartSec(i)], ...
        [yl(1), yl(1), yl(2), yl(2)], c, ...
        'FaceAlpha', 0.055, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(ax, R.bttStartSec(i), ':', 'Color', c, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, R.bttEndSec(i), ':', 'Color', c, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax, R.peakBttTimeSec(i), '--', sprintf('R%d', R.regionId(i)), ...
        'Color', c, 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
ylim(ax, yl);
end

function plot_3d_auxiliary_local(rawFft, timeRange, freqBand, regionTable)
tRaw = rawFft.segmentCenterTimeSec(:).';
freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);
timeMask = tRaw >= timeRange(1) & tRaw <= timeRange(2);
freq = rawFft.freqHz(freqMask);
t = tRaw(timeMask);
A = rawFft.amplitudesMicrostrain(timeMask, freqMask).';
[T, F] = meshgrid(t, freq);
figure('Name', '20251222 Step04 auxiliary raw segmented FFT', ...
    'Color', 'w', 'Position', [90, 90, 1280, 760], 'NumberTitle', 'off');
surf(T, F, A, 'EdgeColor', 'none', 'FaceColor', 'interp');
grid on; box on;
colormap(high_contrast_colormap_local(256));
clim(raw_amplitude_color_limits_local(A, 100));
colorbar; view(42, 34); axis tight;
xlabel('BTT-aligned time (s)'); ylabel('Frequency (Hz)'); zlabel('Amplitude (microstrain)');
title(sprintf('Auxiliary raw segmented FFT, %.0f-%.0f Hz', freqBand(1), freqBand(2)));
hold on;
if ~isempty(regionTable)
    zTop = max(A(:), [], 'omitnan') + 1;
    for i = 1:height(regionTable)
        plot3(regionTable.peakBttTimeSec(i), regionTable.dominantFreqHz(i), zTop, ...
            'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        text(regionTable.peakBttTimeSec(i), regionTable.dominantFreqHz(i), zTop, ...
            sprintf('R%d %dX', regionTable.regionId(i), regionTable.dominantOrder(i)), ...
            'Color', 'w', 'FontWeight', 'bold');
    end
end
end

function sizes = scale_marker_sizes_local(x, smallSize, largeSize)
xMin = min(x);
xMax = max(x);
if xMax <= xMin
    sizes = smallSize * ones(size(x));
else
    sizes = smallSize + (largeSize - smallSize) * (x - xMin) / (xMax - xMin);
end
end

function colors = colors_for_peak_rows_local(T, R)
if isempty(R)
    colors = repmat([1.00 0.80 0.10], height(T), 1);
    return;
end
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
hi = percentile_local(x, highPct);
if ~isfinite(hi) || hi <= 0
    hi = max(x);
end
if ~isfinite(hi) || hi <= 0
    hi = 1;
end
colorLimits = [0, hi];
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

function idx = decimate_index_local(n, maxPoints)
if n <= maxPoints
    idx = 1:n;
else
    idx = unique(round(linspace(1, n, maxPoints)));
end
end

function Fs = estimate_sample_rate_local(t)
dt = median(diff(t), 'omitnan');
if ~isfinite(dt) || dt <= 0
    Fs = 1000;
else
    Fs = 1 / dt;
end
end

function y = moving_rms_local(x, windowSec, Fs)
win = max(5, round(windowSec * Fs));
y = sqrt(movmean(x(:).^2, win, 'omitnan'));
end

function txt = string_or_empty_local(S, fieldName)
if isstruct(S) && isfield(S, fieldName)
    txt = S.(fieldName);
else
    txt = '';
end
if isstring(txt)
    txt = char(txt);
end
if isempty(txt)
    txt = '(empty)';
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
    1.000 0.980 0.550];
x = linspace(0, 1, size(base, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, base, xi, 'pchip');
cmap = min(max(cmap, 0), 1);
end
