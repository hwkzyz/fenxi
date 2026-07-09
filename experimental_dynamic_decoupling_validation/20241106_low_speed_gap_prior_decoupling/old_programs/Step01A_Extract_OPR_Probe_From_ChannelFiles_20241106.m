%% Step01A: extract OPR and probe pulses from 20241106 channel MAT files
% Clean, readable replacement for the OPR/probe branches of maichongcompute.m.
%
% Input:
%   D:\...\20241106_2\900\4-1-1000.mat, 4-2-1000.mat, ...
%   D:\...\20241106_2\3000_3150\4-1-1000.mat, 4-2-1000.mat, ...
%
% Output:
%   By default, this script only prints diagnostics and opens figures.
%   It DOES NOT overwrite existing jiluOPR.mat or jilublade_probe*.mat.

clear; clc; close all;

%% 1. Simple settings: edit here when debugging
dataRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
sampleRateHz = 5e6;

caseList = struct( ...
    'name', {'Low speed 900 rpm', 'High speed 3000-3150 rpm'}, ...
    'dir', {fullfile(dataRoot, '900'), fullfile(dataRoot, '3000_3150')}, ...
    'lastBlock', {20000, 81000});

% Channel 1 is the 20241106 OPR channel in these extracted 4-x files.
oprChannel = 1;

% Probe channels available in this experiment.
probeChannels = [2 3 4 5 6 7];

% Debug selection.
% selectedCaseIndex: [] = both cases, 1 = low speed, 2 = high speed.
% selectedProbeChannels: choose any subset, for example [2 5].
% selectedBlockRange: [] = all blocks, [36000 36000] = only 4-x-36000.mat.
% selectedTimeRangeSec: [] = no time filter, [80 90] = only 80-90 s data.
processOpr = true;
selectedCaseIndex = [];  % []=低速+高速；1=只低速；2=只高速
selectedProbeChannels = probeChannels;
selectedBlockRange = []; % []=全部；[36000 36000]=只处理 4-x-36000.mat
selectedTimeRangeSec = []; % []=不截取；[70 72]=只看 70-72 s

% Thresholds copied from the 20241106 maichongcompute.m idea:
% 2/3: eddy-current-like negative baseline, 4/6: fiber-like, 5/7: capacitive.
thresholdByChannel = containers.Map('KeyType', 'double', 'ValueType', 'double');
thresholdByChannel(1) = 2.0;
thresholdByChannel(2) = -0.2;
thresholdByChannel(3) = -0.2;
thresholdByChannel(4) = 2.0;
thresholdByChannel(5) = 1.0;
thresholdByChannel(6) = 2.0;
thresholdByChannel(7) = 1.0;

% Pulse splitting. Larger values merge nearby above-threshold points into
% one pulse. Old scripts used roughly 1e4 as the initial value.
minPulseGapSamples = 10000;
smoothWindow = 16;

% Trim only the first loaded block, so the first pulse has a clean left edge.
firstBlockTrimSamples = 20000;

% Debug controls.
% Set maxBlocksToRead = 3 to quickly test the first 3 files per channel.
% Set [] to process all available files up to lastBlock.
maxBlocksToRead = [];
maxBlocksEnv = str2double(strtrim(getenv('STEP01A_MAX_BLOCKS')));
if isfinite(maxBlocksEnv) && maxBlocksEnv > 0
    maxBlocksToRead = floor(maxBlocksEnv);
end
caseEnv = sscanf(strtrim(getenv('STEP01A_CASE_INDEX')), '%d').';
if ~isempty(caseEnv)
    selectedCaseIndex = caseEnv;
end
channelsEnv = sscanf(strtrim(getenv('STEP01A_PROBE_CHANNELS')), '%d').';
if ~isempty(channelsEnv)
    selectedProbeChannels = channelsEnv;
end
blockEnv = sscanf(strtrim(getenv('STEP01A_BLOCK_RANGE')), '%d').';
if numel(blockEnv) >= 2
    selectedBlockRange = blockEnv(1:2);
elseif numel(blockEnv) == 1
    selectedBlockRange = [blockEnv(1), blockEnv(1)];
end
timeEnv = sscanf(strtrim(getenv('STEP01A_TIME_RANGE_SEC')), '%f').';
if numel(timeEnv) >= 2
    selectedTimeRangeSec = timeEnv(1:2);
end
plotExamplePulseCount = 3;

% Saving controls. Keep false while tuning thresholds.
saveRecomputedFiles = false;
outputSuffix = '_recomputed_step01a';

%% 2. Extract OPR and all probes
fprintf('\n=== Step01A 20241106 OPR/probe extraction from channel files ===\n');
fprintf('Data root: %s\n', dataRoot);
fprintf('saveRecomputedFiles = %d. Existing files are not overwritten.\n', saveRecomputedFiles);
fprintf('maxBlocksToRead = %s\n', mat2str(maxBlocksToRead));
fprintf('selectedCaseIndex = %s, selectedProbeChannels = %s\n', ...
    mat2str(selectedCaseIndex), mat2str(selectedProbeChannels));
fprintf('selectedBlockRange = %s, selectedTimeRangeSec = %s\n', ...
    mat2str(selectedBlockRange), mat2str(selectedTimeRangeSec));

if isempty(selectedCaseIndex)
    caseIndexList = 1:numel(caseList);
else
    caseIndexList = selectedCaseIndex;
end

for ic = caseIndexList
    C = caseList(ic);
    fprintf('\n================ %s ================\n', C.name);
    fprintf('Folder: %s\n', C.dir);

    if processOpr
        oprThreshold = thresholdByChannel(oprChannel);
        oprResult = extract_signal_local(C.dir, oprChannel, C.lastBlock, maxBlocksToRead, ...
            selectedBlockRange, selectedTimeRangeSec, oprThreshold, minPulseGapSamples, ...
            smoothWindow, firstBlockTrimSamples, sampleRateHz, true);
        existingOpr = load_existing_opr_local(C.dir);
        print_opr_summary_local(oprResult, existingOpr, sampleRateHz);
        plot_opr_summary_local(C.name, oprResult, existingOpr, sampleRateHz, plotExamplePulseCount);

        if saveRecomputedFiles
            jiluOPR = oprResult.sampleStartEnd; %#ok<NASGU>
            save(fullfile(C.dir, ['jiluOPR', outputSuffix, '.mat']), 'jiluOPR');
        end
    end

    for ch = selectedProbeChannels
        threshold = thresholdByChannel(ch);
        probeResult = extract_signal_local(C.dir, ch, C.lastBlock, maxBlocksToRead, ...
            selectedBlockRange, selectedTimeRangeSec, threshold, minPulseGapSamples, ...
            smoothWindow, firstBlockTrimSamples, sampleRateHz, false);
        existingProbe = load_existing_probe_local(C.dir, ch);
        print_probe_summary_local(probeResult, existingProbe, sampleRateHz);
        plot_probe_summary_local(C.name, probeResult, existingProbe, sampleRateHz, plotExamplePulseCount);

        if saveRecomputedFiles
            jilublade = build_jilublade_from_probe_local(probeResult, sampleRateHz); %#ok<NASGU>
            save(fullfile(C.dir, sprintf('jilublade_probe%d%s.mat', ch, outputSuffix)), 'jilublade');
        end
    end
end

fprintf('\nDone. Figures are open in MATLAB. No file was saved unless saveRecomputedFiles=true.\n');

%% Local functions
function result = extract_signal_local(caseDir, channelId, lastBlock, maxBlocksToRead, ...
    selectedBlockRange, selectedTimeRangeSec, thresholdV, minGapSamples, smoothWindow, ...
    firstBlockTrimSamples, sampleRateHz, keepRawExamples)
files = build_channel_file_list_local(caseDir, channelId, lastBlock, maxBlocksToRead, selectedBlockRange);
tail = zeros(0, 2);
sampleStartEnd = zeros(0, 2);
arrivalSample = zeros(0, 1);
examplePulse = struct('sample', {}, 'voltage', {}, 'startSample', {}, 'endSample', {}, 'arrivalSample', {});

for iFile = 1:numel(files)
    raw = load_channel_file_local(files(iFile).path, channelId);
    raw(raw(:, 1) == 0, :) = [];

    if files(iFile).block == 1000 && firstBlockTrimSamples > 0 && size(raw, 1) > firstBlockTrimSamples
        raw(1:firstBlockTrimSamples, :) = [];
    end

    if ~isempty(selectedTimeRangeSec)
        tRaw = raw(:, 1) / sampleRateHz;
        keep = tRaw >= selectedTimeRangeSec(1) & tRaw <= selectedTimeRangeSec(2);
        raw = raw(keep, :);
        if isempty(raw)
            continue;
        end
    end

    raw = [tail; raw]; %#ok<AGROW>
    [pulseStartEnd, pulseArrival, pulseExamples, tail] = detect_pulses_local( ...
        raw, thresholdV, minGapSamples, smoothWindow, keepRawExamples);

    sampleStartEnd = [sampleStartEnd; pulseStartEnd]; %#ok<AGROW>
    arrivalSample = [arrivalSample; pulseArrival]; %#ok<AGROW>
    examplePulse = [examplePulse, pulseExamples]; %#ok<AGROW>
end

result.caseDir = caseDir;
result.channelId = channelId;
result.thresholdV = thresholdV;
result.fileCount = numel(files);
result.sampleStartEnd = sampleStartEnd;
result.arrivalSample = arrivalSample;
result.arrivalTime = arrivalSample / sampleRateHz;
result.examplePulse = examplePulse;
end

function files = build_channel_file_list_local(caseDir, channelId, lastBlock, maxBlocksToRead, selectedBlockRange)
files = struct('path', {}, 'block', {});
for block = 1000:1000:lastBlock
    if ~isempty(selectedBlockRange)
        if block < selectedBlockRange(1) || block > selectedBlockRange(2)
            continue;
        end
    end
    f = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, block));
    if exist(f, 'file') == 2
        files(end + 1).path = f; %#ok<AGROW>
        files(end).block = block;
    end
end

if ~isempty(maxBlocksToRead)
    files = files(1:min(numel(files), maxBlocksToRead));
end

if isempty(files)
    error('No channel files found for channel %d in %s.', channelId, caseDir);
end
end

function raw = load_channel_file_local(filePath, channelId)
varName = sprintf('jilu0%d', channelId);
S = load(filePath, varName);
if ~isfield(S, varName)
    error('File %s does not contain %s.', filePath, varName);
end
raw = S.(varName);
end

function [pulseStartEnd, pulseArrival, pulseExamples, tail] = detect_pulses_local(raw, thresholdV, ...
    minGapSamples, smoothWindow, keepRawExamples)
if isempty(raw)
    pulseStartEnd = zeros(0, 2);
    pulseArrival = zeros(0, 1);
    pulseExamples = struct('sample', {}, 'voltage', {}, 'startSample', {}, 'endSample', {}, 'arrivalSample', {});
    tail = zeros(0, 2);
    return;
end

vSmooth = smooth(raw(:, 2), smoothWindow);
aboveIdx = find(vSmooth > thresholdV);
if isempty(aboveIdx)
    pulseStartEnd = zeros(0, 2);
    pulseArrival = zeros(0, 1);
    pulseExamples = struct('sample', {}, 'voltage', {}, 'startSample', {}, 'endSample', {}, 'arrivalSample', {});
    tail = raw(max(1, end - minGapSamples):end, :);
    return;
end

sampleAbove = raw(aboveIdx, 1);
gap = diff(sampleAbove);
splitAt = find(gap > minGapSamples);
segStart = [1; splitAt + 1];
segEnd = [splitAt; numel(aboveIdx)];

% Keep the last segment as tail because it may continue in the next file.
completeCount = max(0, numel(segStart) - 1);
pulseStartEnd = zeros(completeCount, 2);
pulseArrival = zeros(completeCount, 1);
pulseExamples = struct('sample', {}, 'voltage', {}, 'startSample', {}, 'endSample', {}, 'arrivalSample', {});

for i = 1:completeCount
    aAbove = aboveIdx(segStart(i));
    bAbove = aboveIdx(segEnd(i));
    pad = max(1, floor(0.05 * (bAbove - aAbove + 1)));
    a = max(1, aAbove - pad);
    b = min(size(raw, 1), bAbove + pad);

    pulseStartEnd(i, :) = [sampleAbove(segStart(i)), sampleAbove(segEnd(i))];
    pulseArrival(i) = half_area_arrival_local(raw(a:b, 1), raw(a:b, 2), thresholdV);

    if keepRawExamples && numel(pulseExamples) < 5
        pulseExamples(end + 1).sample = raw(a:b, 1); %#ok<AGROW>
        pulseExamples(end).voltage = raw(a:b, 2);
        pulseExamples(end).startSample = pulseStartEnd(i, 1);
        pulseExamples(end).endSample = pulseStartEnd(i, 2);
        pulseExamples(end).arrivalSample = pulseArrival(i);
    end
end

tailStart = max(1, aboveIdx(segStart(end)) - floor(minGapSamples / 2));
tail = raw(tailStart:end, :);
end

function arrival = half_area_arrival_local(sample, voltage, thresholdV)
sample = sample(:);
voltage = voltage(:);
if numel(sample) < 2
    arrival = sample(1);
    return;
end

% Use positive area above threshold. If this is empty, fall back to midpoint.
weight = voltage - thresholdV;
weight(weight < 0) = 0;
if sum(weight) <= 0
    arrival = mean(sample([1, end]));
    return;
end

cumArea = cumsum(weight);
target = 0.5 * cumArea(end);
idx = find(cumArea >= target, 1, 'first');
if isempty(idx) || idx == 1
    arrival = sample(1);
else
    arrival = 0.5 * (sample(idx - 1) + sample(idx));
end
end

function existing = load_existing_opr_local(caseDir)
f = fullfile(caseDir, 'jiluOPR.mat');
if exist(f, 'file') ~= 2
    existing = [];
    return;
end
S = load(f, 'jiluOPR');
existing.sampleStartEnd = S.jiluOPR(:, 1:2);
existing.arrivalSample = mean(existing.sampleStartEnd, 2);
end

function existing = load_existing_probe_local(caseDir, channelId)
f1 = fullfile(caseDir, sprintf('jilublade_probe%d_nihe.mat', channelId));
f2 = fullfile(caseDir, sprintf('jilublade_probe%d.mat', channelId));
if exist(f1, 'file') == 2
    f = f1;
elseif exist(f2, 'file') == 2
    f = f2;
else
    existing = [];
    return;
end
S = load(f, 'jilublade');
existing.file = f;
existing.jilublade = S.jilublade;
existing.arrivalTime = S.jilublade(:, 3);
end

function print_opr_summary_local(result, existing, sampleRateHz)
fprintf('\nOPR CH%d, threshold %.3f V, files %d\n', ...
    result.channelId, result.thresholdV, result.fileCount);
fprintf('  Detected pulses = %d\n', numel(result.arrivalSample));
if numel(result.arrivalSample) > 1
    t = result.arrivalSample / sampleRateHz;
    rpm = 60 ./ diff(t);
    fprintf('  Time = %.6f to %.6f s, rpm median/min/max = %.2f / %.2f / %.2f\n', ...
        t(1), t(end), median(rpm, 'omitnan'), min(rpm, [], 'omitnan'), max(rpm, [], 'omitnan'));
end
if ~isempty(existing)
    print_arrival_compare_local('  Existing jiluOPR', result.arrivalSample / sampleRateHz, ...
        existing.arrivalSample / sampleRateHz);
end
end

function print_probe_summary_local(result, existing, sampleRateHz)
fprintf('\nProbe CH%d, threshold %.3f V, files %d\n', ...
    result.channelId, result.thresholdV, result.fileCount);
fprintf('  Detected pulses = %d\n', numel(result.arrivalSample));
if ~isempty(result.arrivalSample)
    fprintf('  Time = %.6f to %.6f s\n', result.arrivalSample(1) / sampleRateHz, ...
        result.arrivalSample(end) / sampleRateHz);
end
if ~isempty(existing)
    print_arrival_compare_local(sprintf('  Existing %s', strip_path_local(existing.file)), ...
        result.arrivalSample / sampleRateHz, existing.arrivalTime);
end
end

function print_arrival_compare_local(label, tNew, tOld)
tNew = tNew(:);
tOld = tOld(:);
if isempty(tNew) || isempty(tOld)
    fprintf('%s: no overlap for comparison.\n', label);
    return;
end

[dtUs, nMatched] = nearest_time_difference_us_local(tNew, tOld);
fprintf('%s count = %d, nearest-time compare n = %d, dt us median/maxabs = %.3f / %.3f\n', ...
    label, numel(tOld), nMatched, median(dtUs, 'omitnan'), max(abs(dtUs), [], 'omitnan'));
end

function [dtUs, nMatched] = nearest_time_difference_us_local(tNew, tOld)
% Compare local-window recomputation with the existing full-length file.
% Direct index-by-index comparison is only valid for full runs; nearest-time
% matching also works when selectedBlockRange or selectedTimeRangeSec is used.
tOld = sort(tOld(:));
tNew = tNew(:);
idxNearest = interp1(tOld, 1:numel(tOld), tNew, 'nearest', 'extrap');
idxNearest = max(1, min(numel(tOld), round(idxNearest)));
dtSec = tNew - tOld(idxNearest);

if numel(tOld) > 1
    maxAbsDt = 0.5 * median(diff(tOld), 'omitnan');
else
    maxAbsDt = inf;
end
keep = abs(dtSec) <= maxAbsDt;
dtUs = dtSec(keep) * 1e6;
nMatched = nnz(keep);
end

function plot_opr_summary_local(caseName, result, existing, sampleRateHz, exampleCount)
figure('Name', sprintf('Step01A OPR CH%d %s', result.channelId, caseName), 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

t = result.arrivalSample / sampleRateHz;
nexttile;
if numel(t) > 1
    plot(t(1:end-1), 60 ./ diff(t), 'k-', 'LineWidth', 1.0);
end
grid on; box on;
xlabel('Time (s)');
ylabel('Speed (rpm)');
title([caseName, ': recomputed OPR speed']);

nexttile;
if ~isempty(existing)
    [dtUs, n] = nearest_time_difference_us_local(result.arrivalSample / sampleRateHz, ...
        existing.arrivalSample / sampleRateHz);
    plot((1:n).', dtUs, 'b.');
    ylabel('New - old (us)');
else
    plot(t, result.arrivalSample, 'b.');
    ylabel('Sample');
end
grid on; box on;
xlabel('Pulse index');
title('Comparison with existing jiluOPR.mat');

nexttile;
plot_examples_local(result, sampleRateHz, exampleCount);
title('Example detected OPR pulses');
end

function plot_probe_summary_local(caseName, result, existing, sampleRateHz, exampleCount)
figure('Name', sprintf('Step01A Probe CH%d %s', result.channelId, caseName), 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
t = result.arrivalSample / sampleRateHz;
if numel(t) > 1
    plot(t(1:end-1), diff(t) * 1000, 'k.', 'MarkerSize', 3);
end
grid on; box on;
xlabel('Time (s)');
ylabel('Adjacent interval (ms)');
title(sprintf('%s: probe CH%d interval check', caseName, result.channelId));

nexttile;
if ~isempty(existing)
    [dtUs, n] = nearest_time_difference_us_local(result.arrivalSample / sampleRateHz, ...
        existing.arrivalTime);
    plot((1:n).', dtUs, 'b.');
    grid on; box on;
    xlabel('Pulse index');
    ylabel('New - old (us)');
    title('Comparison with existing jilublade file');
else
    plot_examples_local(result, sampleRateHz, exampleCount);
    title('Example detected probe pulses');
end
end

function plot_examples_local(result, sampleRateHz, exampleCount)
hold on;
n = min([numel(result.examplePulse), exampleCount]);
for i = 1:n
    E = result.examplePulse(i);
    plot((E.sample - E.sample(1)) / sampleRateHz * 1e3, E.voltage, 'LineWidth', 1.0);
    xline((E.arrivalSample - E.sample(1)) / sampleRateHz * 1e3, '--');
end
grid on; box on;
xlabel('Local time (ms)');
ylabel('Voltage (V)');
yline(result.thresholdV, 'r--', 'threshold');
end

function jilublade = build_jilublade_from_probe_local(result, sampleRateHz)
n = size(result.sampleStartEnd, 1);
jilublade = zeros(n, 4);
jilublade(:, 1) = result.sampleStartEnd(:, 1) / sampleRateHz;
jilublade(:, 2) = result.sampleStartEnd(:, 2) / sampleRateHz;
jilublade(:, 3) = result.arrivalSample(:) / sampleRateHz;
end

function s = strip_path_local(filePath)
[~, name, ext] = fileparts(filePath);
s = [name, ext];
end
