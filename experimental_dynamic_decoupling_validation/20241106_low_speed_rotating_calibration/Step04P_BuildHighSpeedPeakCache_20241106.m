%% Step04P_BuildHighSpeedPeakCache_20241106
% Build reusable per-sensor high-speed peak caches for the whole case.
% This step does the expensive 7th-order polynomial peak fitting only once.
% Later Step04 / Step06 can reuse the fitted six-peak values when scanning
% different startTimeSec regions.

clear; close all; clc;

P = NewFlow_Config_20241106();

%% Parameters to tune
S04P = struct();
S04P.analysisSensors = [2 3 5 7];
S04P.matchPolyDegree = P.numbering.matchPolyDegree;
S04P.sampleRateHz = P.machine.sampleRateHz;
S04P.pulsePadSec = 0;
S04P.forceRebuild = false;
S04P.progressEveryRows = 1000;

P = apply_step04p_local_options_local(P, S04P);
analysisSensors = P.sensors.analysis;

fprintf('\n=== Step04P: build high-speed peak cache ===\n');
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Fit degree: %d\n', S04P.matchPolyDegree);
fprintf('Pulse pad: %.9f s\n', S04P.pulsePadSec);
drawnow;

oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
require_file_local(oprFile, 'high-speed OPR pulse file');
loadedOpr = load(oprFile, 'jiluOPR');
oprTimes = loadedOpr.jiluOPR(:, 1);

fileRanges = build_case_file_ranges_local(P.data.highSpeedDir, P.machine.oprChannel, S04P.sampleRateHz);

for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    outFile = peak_cache_file_local(P, sid);
    if exist(outFile, 'file') == 2 && ~S04P.forceRebuild
        loaded = load(outFile, 'HighSpeedPeakCacheSensor');
        if isfield(loaded, 'HighSpeedPeakCacheSensor') && ...
                cache_matches_request_local(loaded.HighSpeedPeakCacheSensor, sid, P, S04P)
            fprintf('CH%d cache already matches current settings. Skip: %s\n', sid, outFile);
            drawnow;
            continue;
        end
    end

    fprintf('Building CH%d peak cache (%d/%d) ...\n', sid, i, numel(analysisSensors));
    drawnow;

    probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
    require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
    loadedProbe = load(probeFile, 'jilublade');
    jilublade = loadedProbe.jilublade;
    revolutionIndex = map_pulses_to_revolutions_local(jilublade(:, 3), oprTimes);

    [tRaw, vRaw] = load_full_raw_channel_local(P.data.highSpeedDir, fileRanges, sid, S04P.sampleRateHz);
    cache = build_sensor_peak_cache_local( ...
        sid, jilublade, revolutionIndex, tRaw, vRaw, P, S04P);

    outDir = fileparts(outFile);
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    HighSpeedPeakCacheSensor = cache; %#ok<NASGU>
    save(outFile, 'HighSpeedPeakCacheSensor', '-v7.3');
    fprintf('Saved CH%d cache: %s\n', sid, outFile);
    drawnow;
end

function P = apply_step04p_local_options_local(P, S04P)
P.sensors.analysis = S04P.analysisSensors(:).';
end

function outFile = peak_cache_file_local(P, sid)
outFile = fullfile(P.outputDir, '04P_high_speed_peak_cache', ...
    sprintf('HighSpeedPeakCache_CH%d_20241106.mat', sid));
end

function tf = cache_matches_request_local(cache, sid, P, S04P)
tf = isstruct(cache) && ...
    isfield(cache, 'sensor_id') && isequal(cache.sensor_id, sid) && ...
    isfield(cache, 'dataset') && strcmpi(string(cache.dataset), "20241106") && ...
    isfield(cache, 'high_speed_case') && strcmpi(string(cache.high_speed_case), string(P.data.highSpeedCase)) && ...
    isfield(cache, 'sample_rate_hz') && isequal(cache.sample_rate_hz, S04P.sampleRateHz) && ...
    isfield(cache, 'match_poly_degree') && isequal(cache.match_poly_degree, S04P.matchPolyDegree) && ...
    isfield(cache, 'pulse_pad_sec') && abs(cache.pulse_pad_sec - S04P.pulsePadSec) <= eps;
end

function cache = build_sensor_peak_cache_local(sensorId, jilublade, revolutionIndex, tRaw, vRaw, P, S04P)
nRows = size(jilublade, 1);
cache = struct();
cache.dataset = '20241106';
cache.high_speed_case = P.data.highSpeedCase;
cache.sensor_id = sensorId;
cache.sample_rate_hz = S04P.sampleRateHz;
cache.match_poly_degree = S04P.matchPolyDegree;
cache.pulse_pad_sec = S04P.pulsePadSec;
cache.row_count = nRows;
cache.row_id = (1:nRows).';
cache.t_start = jilublade(:, 1) - S04P.pulsePadSec;
cache.t_end = jilublade(:, 2) + S04P.pulsePadSec;
cache.t_peak_nominal = jilublade(:, 3);
cache.revolution_index = revolutionIndex(:);
cache.raw_peak_value = nan(nRows, 1);
cache.raw_peak_time = nan(nRows, 1);
cache.fit_peak_value = nan(nRows, 1);
cache.fit_peak_time = nan(nRows, 1);
cache.fit_succeeded = false(nRows, 1);
cache.degree_used = nan(nRows, 1);
cache.sample_count = zeros(nRows, 1);

    [~, order] = sort(cache.t_peak_nominal, 'ascend');
    sortedRowIds = cache.row_id(order);
    sortedStarts = cache.t_start(order);
    sortedEnds = cache.t_end(order);

    idxStart = 1;
    idxEnd = 1;
    nRaw = numel(tRaw);
    for k = 1:nRows
        rowId = sortedRowIds(k);
        t0 = sortedStarts(k);
        t1 = sortedEnds(k);

        while idxStart <= nRaw && tRaw(idxStart) < t0
            idxStart = idxStart + 1;
        end
        if idxEnd < idxStart
            idxEnd = idxStart;
        end
        while idxEnd <= nRaw && tRaw(idxEnd) <= t1
            idxEnd = idxEnd + 1;
        end
        i0 = idxStart;
        i1 = idxEnd - 1;
        if i0 > i1 || i0 > nRaw || i1 < 1
            continue;
        end

        tSeg = tRaw(i0:i1);
        vSeg = vRaw(i0:i1);
        cache.sample_count(rowId) = numel(tSeg);
        if numel(tSeg) < max(8, S04P.matchPolyDegree + 2)
            continue;
        end
        fitInfo = fit_polynomial_peak_trace_local(tSeg, vSeg, S04P.matchPolyDegree);
        cache.raw_peak_value(rowId) = fitInfo.rawPeakValue;
        cache.raw_peak_time(rowId) = fitInfo.rawPeakTime;
        cache.fit_peak_value(rowId) = fitInfo.peakValue;
        cache.fit_peak_time(rowId) = fitInfo.peakTime;
        cache.fit_succeeded(rowId) = fitInfo.fitSucceeded;
        cache.degree_used(rowId) = fitInfo.degreeUsed;

        if mod(k, max(1, S04P.progressEveryRows)) == 0 || k == nRows
            fprintf('  CH%d rows: %d / %d\n', sensorId, k, nRows);
            drawnow;
        end
    end
end

function [tRaw, vRaw] = load_full_raw_channel_local(caseDir, fileRanges, sensorId, sampleRateHz)
tRaw = [];
vRaw = [];
for k = 1:numel(fileRanges)
    fileId = fileRanges(k).file_id;
    offset = fileRanges(k).offset;
    [tLocal, vLocal] = load_raw_channel_local(caseDir, sensorId, fileId, sampleRateHz);
    tRaw = [tRaw; tLocal(:) + offset]; %#ok<AGROW>
    vRaw = [vRaw; vLocal(:)]; %#ok<AGROW>
end
end

function fitInfo = fit_polynomial_peak_trace_local(t, v, degree)
fitInfo = struct( ...
    'peakValue', NaN, ...
    'peakTime', NaN, ...
    'rawPeakValue', NaN, ...
    'rawPeakTime', NaN, ...
    'fitSucceeded', false, ...
    'degreeUsed', NaN);
t = t(:);
v = v(:);
valid = isfinite(t) & isfinite(v);
t = t(valid);
v = v(valid);
if numel(t) < max(3, degree + 1) || range(t) <= eps
    return;
end
[t, order] = sort(t);
v = v(order);
[rawPeakValue, rawIdx] = max(v);
fitInfo.rawPeakValue = rawPeakValue;
fitInfo.rawPeakTime = t(rawIdx);
realDegree = min(degree, numel(t) - 1);
fitInfo.degreeUsed = realDegree;
[coef, ~, mu] = polyfit(t, v, realDegree);
tFine = linspace(min(t), max(t), 200).';
vFine = polyval(coef, tFine, [], mu);
[peakValue, peakIdx] = max(vFine);
fitInfo.peakValue = peakValue;
fitInfo.peakTime = tFine(peakIdx);
fitInfo.fitSucceeded = true;
end

function fileRanges = build_case_file_ranges_local(caseDir, oprChannel, sampleRateHz)
files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.dat', oprChannel)));
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.txt', oprChannel)));
end
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('4-%d-*.mat', oprChannel)));
end
if isempty(files)
    error('No OPR raw files found in %s.', caseDir);
end
fileIds = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, 'Data_(\d+)', 'tokens', 'once');
    if isempty(token)
        token = regexp(files(i).name, '4-\d+-(\d+)\.mat', 'tokens', 'once');
    end
    if isempty(token)
        error('Could not parse file id from %s.', files(i).name);
    end
    fileIds(i) = str2double(token{1});
end
[fileIds, order] = sort(fileIds);
files = files(order); %#ok<NASGU>
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = 0;
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_channel_local(caseDir, oprChannel, fileIds(i), sampleRateHz);
    if i == 1
        offset = 0;
    else
        offset = lastEnd + 1 / sampleRateHz - tOpr(1);
    end
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
    lastEnd = fileRanges(i).t_end;
end
end

function [tSec, v] = load_raw_channel_local(caseDir, channelId, fileId, sampleRateHz)
patterns = { ...
    sprintf('Probe%d_Data_%d.dat', channelId, fileId), ...
    sprintf('Probe%d_Data_%d.txt', channelId, fileId), ...
    sprintf('4-%d-%d.mat', channelId, fileId)};
filePath = '';
for i = 1:numel(patterns)
    candidate = fullfile(caseDir, patterns{i});
    if exist(candidate, 'file') == 2
        filePath = candidate;
        break;
    end
end
if isempty(filePath)
    error('Missing raw file for CH%d file %d in %s.', channelId, fileId, caseDir);
end
varName = sprintf('jilu%02d', channelId);
vars = whos('-file', filePath);
if ~ismember(varName, {vars.name})
    error('Raw MAT file does not contain %s: %s', varName, filePath);
end
loaded = load(filePath, varName);
raw = loaded.(varName);
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / sampleRateHz;
v = raw(:, 2);
end

function revIndex = map_pulses_to_revolutions_local(pulseTimes, oprTimes)
revIndex = nan(size(pulseTimes));
for i = 1:numel(pulseTimes)
    idx = find(oprTimes <= pulseTimes(i), 1, 'last');
    if ~isempty(idx)
        revIndex(i) = idx;
    end
end
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end
