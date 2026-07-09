%% Step07K: check 20241106 blade-slot and waveform selection near 75 s
% Purpose:
%   Before trusting Step07J, check whether the selected blade slot is really
%   the resonant / instrumented blade and whether raw pulse waveforms are
%   centered correctly.
%
% This script only uses the small steady resonance window near 75 s:
%   time window = 72-78 s, probes = P5/P7, blade slots = 1-6.
%
% Diagnostics:
%   1. xAtPeakMedianMm: pulse peak should stay close to x=0 if waveform
%      cutting is centered correctly.
%   2. peakStdMv / changeRmsMv: large values indicate strong dynamic
%      waveform variation for that blade slot.
%   3. bttStdMm: Step02 timing displacement scatter in the same window.
%   4. bestEo / bestEoAmpMv: EO content in peak-amplitude modulation.

clear; clc; close all;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

step02File = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
if ~isfile(step02File)
    error('Run Step02 first. Missing file: %s', step02File);
end
S02 = load(step02File);

timeWindowSec = [72 78];
analysisSensors = [5 7];
bladeSlots = 1:S02.bladeCount;
pulsesPerSensorBlade = 24;
pulseHalfWindowSec = 2.5e-4;
pointsPerPulse = 151;
edgeFractionForBaseline = 0.18;
blockStepSamples = 1e7;
blockLabelStep = 1000;
rawCacheMaxEntries = 8;
eoCandidates = 8:16;

analysisSensors = parse_int_env_local('STEP07K_ANALYSIS_SENSORS', analysisSensors);
bladeSlots = parse_int_env_local('STEP07K_BLADE_SLOTS', bladeSlots);
timeWindowSec = parse_range_env_local('STEP07K_TIME_WINDOW_SEC', timeWindowSec);
pulsesPerSensorBlade = round(parse_scalar_env_local('STEP07K_PULSES_PER_SLOT', pulsesPerSensorBlade));
pointsPerPulse = round(parse_scalar_env_local('STEP07K_POINTS_PER_PULSE', pointsPerPulse));

fprintf('\n=== Step07K 20241106 blade-slot / waveform-selection check ===\n');
fprintf('High-speed folder: %s\n', S02.highDir);
fprintf('Time window: %.2f-%.2f s\n', timeWindowSec(1), timeWindowSec(2));
fprintf('Sensors: %s, blade slots: %s\n', mat2str(analysisSensors), mat2str(bladeSlots));
fprintf('Pulses per sensor/blade: %d, points per pulse: %d\n', ...
    pulsesPerSensorBlade, pointsPerPulse);

opr = load_opr_local(S02.highDir, S02.sampleRateHz);
rawCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

rows = {};
slotScores = table();
examples = struct('sensorId', {}, 'bladeSlot', {}, 'eventTime', {}, ...
    'xMm', {}, 'vMv', {});

%% 2. Extract small-window pulses for every slot
for sid = analysisSensors
    sensorIdx = find([S02.highDisp.sensorId] == sid, 1, 'first');
    if isempty(sensorIdx)
        warning('No Step02 data for sensor %d. Skipped.', sid);
        continue;
    end

    for bladeSlot = bladeSlots
        D = S02.highDisp(sensorIdx).slot(bladeSlot);
        keepEvent = D.time(:) >= timeWindowSec(1) & D.time(:) <= timeWindowSec(2);
        eventTimes = D.time(keepEvent);
        bttX = D.xMm(keepEvent);
        if isempty(eventTimes)
            continue;
        end

        pick = pick_evenly_local(numel(eventTimes), pulsesPerSensorBlade);
        eventTimes = eventTimes(pick);
        bttX = bttX(pick);

        pulseGrid = [];
        peakMv = NaN(numel(eventTimes), 1);
        xAtPeak = NaN(numel(eventTimes), 1);
        acceptedTime = NaN(numel(eventTimes), 1);
        rotFreqHz = NaN(numel(eventTimes), 1);

        for ie = 1:numel(eventTimes)
            tEvent = eventTimes(ie);
            rotFreqHz(ie) = local_rot_freq_local(opr.tCenter, tEvent);
            tipSpeedMmS = 2 * pi * S02.rTipMm * rotFreqHz(ie);

            [tq, vq, ok] = read_raw_pulse_local(S02.highDir, sid, tEvent, ...
                pulseHalfWindowSec, pointsPerPulse, S02.sampleRateHz, ...
                blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries);
            if ~ok
                continue;
            end

            x = (tq - tEvent) * tipSpeedMmS;
            edgeCount = max(5, round(edgeFractionForBaseline * numel(vq)));
            baselineV = median([vq(1:edgeCount); vq(end-edgeCount+1:end)], 'omitnan');
            vMv = 1000 * (vq - baselineV);

            if isempty(pulseGrid)
                pulseGrid = NaN(numel(x), numel(eventTimes));
                xGrid = x(:);
            end
            pulseGrid(:, ie) = interp1(x(:), vMv(:), xGrid, 'linear', NaN);
            [peakMv(ie), ip] = max(vMv, [], 'omitnan');
            if isfinite(peakMv(ie))
                xAtPeak(ie) = x(ip);
            end
            acceptedTime(ie) = tEvent;

            if numel(examples) < 18 && ismember(bladeSlot, [1 2 3 4 5 6])
                examples(end + 1).sensorId = sid; %#ok<SAGROW>
                examples(end).bladeSlot = bladeSlot;
                examples(end).eventTime = tEvent;
                examples(end).xMm = x(:);
                examples(end).vMv = vMv(:);
            end
        end

        goodPulse = isfinite(peakMv) & isfinite(xAtPeak);
        if nnz(goodPulse) < 6
            continue;
        end

        template = median(pulseGrid(:, goodPulse), 2, 'omitnan');
        diffPulse = pulseGrid(:, goodPulse) - template;
        changeRmsMv = sqrt(mean(diffPulse(:).^2, 'omitnan'));

        y = peakMv(goodPulse);
        t = acceptedTime(goodPulse);
        fRot = mean(rotFreqHz(goodPulse), 'omitnan');
        eoAmp = NaN(size(eoCandidates));
        for ieo = 1:numel(eoCandidates)
            eoAmp(ieo) = sinusoid_amplitude_local(t, y, eoCandidates(ieo) * fRot);
        end
        [bestEoAmpMv, ibest] = max(eoAmp);
        bestEo = eoCandidates(ibest);

        rows{end + 1, 1} = table(sid, bladeSlot, nnz(goodPulse), ...
            median(peakMv(goodPulse), 'omitnan'), std(peakMv(goodPulse), 'omitnan'), ...
            changeRmsMv, median(xAtPeak(goodPulse), 'omitnan'), ...
            std(xAtPeak(goodPulse), 'omitnan'), std(bttX(goodPulse), 'omitnan'), ...
            bestEo, bestEoAmpMv, fRot, ...
            'VariableNames', {'sensorId','bladeSlot','pulseCount', ...
            'peakMedianMv','peakStdMv','changeRmsMv','xAtPeakMedianMm', ...
            'xAtPeakStdMm','bttStdMm','bestEo','bestEoAmpMv','rotFreqHz'}); %#ok<SAGROW>
    end
end

if isempty(rows)
    error('No usable pulses found.');
end
diagnosticTable = vertcat(rows{:});

%% 3. Combine P5/P7 into one score per blade slot
slotRows = {};
for bladeSlot = bladeSlots
    T = diagnosticTable(diagnosticTable.bladeSlot == bladeSlot, :);
    if isempty(T)
        continue;
    end
    peakScore = mean(T.peakStdMv ./ max(T.peakMedianMv, eps), 'omitnan');
    changeScore = mean(T.changeRmsMv ./ max(T.peakMedianMv, eps), 'omitnan');
    bttScore = mean(T.bttStdMm, 'omitnan');
    eoScore = mean(T.bestEoAmpMv ./ max(T.peakMedianMv, eps), 'omitnan');
    centerError = mean(abs(T.xAtPeakMedianMm), 'omitnan');
    totalScore = 0.35 * normalize_positive_local(peakScore) + ...
        0.30 * normalize_positive_local(changeScore) + ...
        0.20 * normalize_positive_local(bttScore) + ...
        0.15 * normalize_positive_local(eoScore);
    slotRows{end + 1, 1} = table(bladeSlot, peakScore, changeScore, bttScore, ...
        eoScore, centerError, totalScore, ...
        'VariableNames', {'bladeSlot','relativePeakStd','relativeChangeRms', ...
        'bttStdMm','relativeEoAmp','meanAbsXAtPeakMm','roughScore'}); %#ok<SAGROW>
end
slotScores = vertcat(slotRows{:});

% Normalize roughScore after all rows are available.
scoreRaw = 0.35 * rescale_safe_local(slotScores.relativePeakStd) + ...
    0.30 * rescale_safe_local(slotScores.relativeChangeRms) + ...
    0.20 * rescale_safe_local(slotScores.bttStdMm) + ...
    0.15 * rescale_safe_local(slotScores.relativeEoAmp);
slotScores.roughScore = scoreRaw;
[~, bestIdx] = max(slotScores.roughScore);
recommendedBladeSlot = slotScores.bladeSlot(bestIdx);

%% 4. Save and plot
csvFile = fullfile(outDir, 'Step07K_BladeSlot_WaveformSelection_20241106.csv');
slotCsvFile = fullfile(outDir, 'Step07K_BladeSlot_RoughScore_20241106.csv');
matFile = fullfile(outDir, 'Step07K_BladeSlot_WaveformSelection_20241106.mat');
writetable(diagnosticTable, csvFile);
writetable(slotScores, slotCsvFile);
save(matFile, 'diagnosticTable', 'slotScores', 'recommendedBladeSlot', ...
    'examples', 'timeWindowSec', 'analysisSensors', 'bladeSlots', '-v7.3');

plot_diagnostics_local(diagnosticTable, slotScores, examples, recommendedBladeSlot);

fprintf('\nStep07K complete.\n');
fprintf('Recommended blade slot by rough waveform/BTT score: B%d\n', recommendedBladeSlot);
disp(slotScores);
fprintf('Saved:\n  %s\n  %s\n', csvFile, slotCsvFile);

%% Local functions
function values = parse_int_env_local(name, defaultValues)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValues;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    values = defaultValues;
end
end

function value = parse_scalar_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
tmp = str2double(raw);
if isfinite(tmp)
    value = tmp;
else
    value = defaultValue;
end
end

function values = parse_range_env_local(name, defaultValues)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValues;
    return;
end
values = sscanf(raw, '%f').';
if numel(values) < 2
    values = defaultValues;
else
    values = values(1:2);
end
end

function opr = load_opr_local(caseDir, sampleRateHz)
f = fullfile(caseDir, 'jiluOPR.mat');
S = load(f, 'jiluOPR');
opr.tCenter = mean(S.jiluOPR(:, 1:2), 2) / sampleRateHz;
end

function idx = pick_evenly_local(n, maxCount)
if n <= maxCount
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxCount))).';
end
end

function f = local_rot_freq_local(oprTime, t)
idx = find(oprTime <= t, 1, 'last');
if isempty(idx)
    idx = 1;
end
idx = min(idx, numel(oprTime) - 1);
dt = oprTime(idx + 1) - oprTime(idx);
if ~(isfinite(dt) && dt > 0)
    dt = median(diff(oprTime), 'omitnan');
end
f = 1 / dt;
end

function amp = sinusoid_amplitude_local(t, y, fHz)
t = t(:);
y = y(:);
mask = isfinite(t) & isfinite(y);
t = t(mask);
y = y(mask);
if numel(t) < 4
    amp = NaN;
    return;
end
y = y - mean(y, 'omitnan');
A = [sin(2*pi*fHz*t), cos(2*pi*fHz*t), ones(size(t))];
coef = A \ y;
amp = hypot(coef(1), coef(2));
end

function y = normalize_positive_local(x)
if ~isfinite(x)
    y = 0;
else
    y = max(0, x);
end
end

function y = rescale_safe_local(x)
x = x(:);
if all(~isfinite(x)) || max(x, [], 'omitnan') <= min(x, [], 'omitnan')
    y = zeros(size(x));
else
    y = (x - min(x, [], 'omitnan')) ./ ...
        (max(x, [], 'omitnan') - min(x, [], 'omitnan'));
end
end

function [tq, vq, ok] = read_raw_pulse_local(caseDir, sensorId, tCenter, halfWindowSec, ...
    pointsPerPulse, sampleRateHz, blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries)
tq = linspace(tCenter - halfWindowSec, tCenter + halfWindowSec, pointsPerPulse).';
sampleStart = floor((tCenter - halfWindowSec) * sampleRateHz);
sampleEnd = ceil((tCenter + halfWindowSec) * sampleRateHz);
blockStart = sample_to_block_label_local(sampleStart, blockStepSamples, blockLabelStep);
blockEnd = sample_to_block_label_local(sampleEnd, blockStepSamples, blockLabelStep);
rawAll = zeros(0, 2);
for blockLabel = blockStart:blockLabelStep:blockEnd
    raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries);
    rawAll = [rawAll; raw]; %#ok<AGROW>
end
if isempty(rawAll)
    vq = NaN(size(tq));
    ok = false;
    return;
end
rawAll(rawAll(:, 1) == 0, :) = [];
rawAll = sortrows(rawAll, 1);
rawT = rawAll(:, 1) / sampleRateHz;
rawV = rawAll(:, 2);
keep = rawT >= tq(1) & rawT <= tq(end) & isfinite(rawV);
if nnz(keep) < 4
    vq = NaN(size(tq));
    ok = false;
    return;
end
vq = interp1(rawT(keep), rawV(keep), tq, 'linear', NaN);
ok = nnz(isfinite(vq)) >= 0.8 * numel(vq);
end

function blockLabel = sample_to_block_label_local(sample, blockStepSamples, blockLabelStep)
blockIndex = ceil((sample - 10000) / blockStepSamples);
blockIndex = max(1, blockIndex);
blockLabel = blockIndex * blockLabelStep;
end

function raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries)
key = sprintf('S%d_B%d', sensorId, blockLabel);
if isKey(rawCache, key)
    raw = rawCache(key);
    return;
end
filePath = fullfile(caseDir, sprintf('4-%d-%d.mat', sensorId, blockLabel));
if ~isfile(filePath)
    raw = zeros(0, 2);
    rawCache(key) = raw;
    return;
end
varName = sprintf('jilu0%d', sensorId);
S = load(filePath, varName);
raw = S.(varName);
rawCache(key) = raw;
if rawCache.Count > rawCacheMaxEntries
    cacheKeys = keys(rawCache);
    remove(rawCache, cacheKeys(1:max(1, numel(cacheKeys) - rawCacheMaxEntries)));
end
end

function plot_diagnostics_local(T, slotScores, examples, recommendedBladeSlot)
figure('Name', '20241106 Step07K blade-slot waveform check', 'Color', 'w', ...
    'Position', [80, 60, 1450, 850], 'NumberTitle', 'off');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for sid = unique(T.sensorId).'
    M = T(T.sensorId == sid, :);
    plot(M.bladeSlot, M.xAtPeakMedianMm, 'o-', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('P%d', sid));
end
yline(0, '--k');
xline(recommendedBladeSlot, 'r--', 'recommended');
xlabel('Blade slot');
ylabel('Median x at pulse peak (mm)');
title('Waveform centering check');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
for sid = unique(T.sensorId).'
    M = T(T.sensorId == sid, :);
    plot(M.bladeSlot, M.peakStdMv, 'o-', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('P%d', sid));
end
xline(recommendedBladeSlot, 'r--');
xlabel('Blade slot');
ylabel('Peak std (mV)');
title('Pulse peak modulation');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
for sid = unique(T.sensorId).'
    M = T(T.sensorId == sid, :);
    plot(M.bladeSlot, M.changeRmsMv, 'o-', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('P%d', sid));
end
xline(recommendedBladeSlot, 'r--');
xlabel('Blade slot');
ylabel('Pulse-to-template RMS (mV)');
title('Whole-waveform variation');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
bar(slotScores.bladeSlot, slotScores.bttStdMm);
xline(recommendedBladeSlot, 'r--');
xlabel('Blade slot');
ylabel('BTT std (mm)');
title('Step02 timing displacement scatter');

nexttile; hold on; grid on; box on;
bar(slotScores.bladeSlot, slotScores.roughScore);
xline(recommendedBladeSlot, 'r--');
xlabel('Blade slot');
ylabel('Rough score');
title('Combined rough blade-slot score');

nexttile; hold on; grid on; box on;
for i = 1:min(numel(examples), 18)
    E = examples(i);
    if E.bladeSlot == recommendedBladeSlot
        lw = 1.4;
    else
        lw = 0.6;
    end
    plot(E.xMm, E.vMv, 'LineWidth', lw, ...
        'DisplayName', sprintf('P%d-B%d', E.sensorId, E.bladeSlot));
end
xlabel('Pulse-local x (mm)');
ylabel('Voltage (mV)');
title('Example extracted pulses');
legend('Location', 'bestoutside');
end
