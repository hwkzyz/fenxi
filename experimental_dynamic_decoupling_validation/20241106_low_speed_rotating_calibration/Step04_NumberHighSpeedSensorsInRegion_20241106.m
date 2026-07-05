%% Step04_NumberHighSpeedSensorsInRegion_20241106
% Number high-speed sensors in the selected region using the same logic as
% the low-speed chain:
% 1. OPR defines each revolution.
% 2. Each valid revolution contributes exactly six fitted peaks.
% 3. Every revolution is cyclically matched to the low-speed reference.

clear; close all; clc;

%% Parameters to tune
P = NewFlow_Config_20241106();
S04 = struct();
S04.analysisSensors = [2 3 5 7];
S04.startTimeSec = 75.0;
S04.bladeCount = P.machine.bladeCount;
S04.targetLaps = P.region.targetLaps;
S04.matchPolyDegree = P.numbering.matchPolyDegree;
S04.sampleRateHz = P.machine.sampleRateHz;
S04.pulsePadSec = 0;
S04.usePeakCache = true;
S04.autoBuildPeakCacheIfMissing = true;
P = apply_step04_local_options_local(P, S04);

analysisSensors = P.sensors.analysis;
bladeCount = S04.bladeCount;
matchPolyDegree = S04.matchPolyDegree;
sampleRateHz = S04.sampleRateHz;
pulsePadSec = S04.pulsePadSec;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
lowSpeedFingerprintFile = P.files.lowSpeedFingerprint;
outFile = P.files.highSpeedNumbering;
outDir = fileparts(outFile);

highSpeedDir = P.data.highSpeedDir;
highSpeedPulseDir = P.data.highSpeedPulseDir;
jiluFile = fullfile(highSpeedPulseDir, 'jiluOPR.mat');

if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

%% Load upstream artifacts
require_file_local(lowSpeedFingerprintFile, 'Step01 low-speed fingerprint');
require_file_local(jiluFile, 'high-speed OPR pulse file');

loadedRef = load(lowSpeedFingerprintFile, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;
RegionSelection = build_step04_region_selection_local(S04);
validate_upstream_artifacts_local(LowSpeedReference, RegionSelection, analysisSensors);

loadedOpr = load(jiluFile, 'jiluOPR');
oprTimes = loadedOpr.jiluOPR(:, 1);

%% Build per-sensor revolution libraries inside the selected region
probeCache = load_probe_cache_local(highSpeedPulseDir, analysisSensors, oprTimes);
selectedRevIds = select_region_revolutions_local(probeCache, RegionSelection, P);
timeWindow = build_selected_region_time_window_local(probeCache, selectedRevIds, pulsePadSec);
peakCacheSet = [];
if S04.usePeakCache
    peakCacheSet = ensure_peak_cache_set_local(P, analysisSensors, matchPolyDegree, sampleRateHz, pulsePadSec, S04.autoBuildPeakCacheIfMissing);
end

raw = [];
if isempty(peakCacheSet)
    fileRanges = build_case_file_ranges_local(highSpeedDir, P.machine.oprChannel, sampleRateHz);
    raw = load_raw_subset_local(highSpeedDir, fileRanges, analysisSensors, sampleRateHz, timeWindow);
end

sensorResults = repmat(empty_sensor_result_local(bladeCount), numel(analysisSensors), 1);
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ValidRevolutionCount', NaN, ...
    'MeanBestScore', NaN, ...
    'MinBestScore', NaN, ...
    'MeanScoreMargin', NaN), numel(analysisSensors), 1);

for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    probe = probeCache(i);
    refFingerprint = read_low_speed_fingerprint_local(LowSpeedReference, sid, bladeCount);
    if ~isempty(peakCacheSet)
        [revLib, numberingTable, summaryStruct] = build_high_speed_revolution_library_from_cache_local( ...
            peakCacheSet{i}, probe, selectedRevIds, refFingerprint, bladeCount);
    else
        [revLib, numberingTable, summaryStruct] = build_high_speed_revolution_library_local( ...
            raw(sid), probe, selectedRevIds, refFingerprint, bladeCount, pulsePadSec, matchPolyDegree);
    end

    sensorResults(i).sensor_id = sid;
    sensorResults(i).dominant_shift = summaryStruct.dominant_shift;
    sensorResults(i).local_to_physical_blade_ids = summaryStruct.local_to_physical;
    sensorResults(i).physical_to_local_blade_ids = summaryStruct.physical_to_local;
    sensorResults(i).best_score = summaryStruct.mean_best_score;
    sensorResults(i).score_margin = summaryStruct.mean_score_margin;
    sensorResults(i).selected_revolution_ids = revLib.revolution_ids(:);
    sensorResults(i).selected_rows_by_physical = revLib.selected_rows_by_physical;
    sensorResults(i).revolution_local_rows = revLib.local_rows;
    sensorResults(i).revolution_local_peak_values = revLib.local_peak_values;
    sensorResults(i).revolution_physical_peak_values = revLib.physical_peak_values;
    sensorResults(i).revolution_best_shifts = revLib.best_shifts(:);
    sensorResults(i).revolution_best_scores = revLib.best_scores(:);
    sensorResults(i).revolution_score_margins = revLib.score_margins(:);
    sensorResults(i).revolution_shift_scores = revLib.shift_scores;
    sensorResults(i).revolution_peak_times = revLib.peak_times;
    sensorResults(i).revolution_row_windows = revLib.row_windows;
    sensorResults(i).numbering_table = numberingTable;

    summaryRows(i).SensorID = sid;
    summaryRows(i).DominantShift = summaryStruct.dominant_shift;
    summaryRows(i).LocalSlotOfB1 = summaryStruct.physical_to_local(1);
    summaryRows(i).ValidRevolutionCount = numel(revLib.revolution_ids);
    summaryRows(i).MeanBestScore = summaryStruct.mean_best_score;
    summaryRows(i).MinBestScore = summaryStruct.min_best_score;
    summaryRows(i).MeanScoreMargin = summaryStruct.mean_score_margin;
end

HighSpeedNumbering = struct();
HighSpeedNumbering.mode = 'opr_defined_revolutions_with_six_peak_cyclic_matching';
HighSpeedNumbering.analysis_sensors = analysisSensors;
HighSpeedNumbering.blade_count = bladeCount;
HighSpeedNumbering.region_selection = RegionSelection;
HighSpeedNumbering.selected_revolution_ids = selectedRevIds(:);
HighSpeedNumbering.time_window = timeWindow(:).';
HighSpeedNumbering.low_speed_reference_file = lowSpeedFingerprintFile;
HighSpeedNumbering.match_poly_degree = matchPolyDegree;
HighSpeedNumbering.summary = struct2table(summaryRows);
HighSpeedNumbering.sensor = sensorResults;
HighSpeedNumbering.note = ['Each sensor contains a per-revolution six-peak table in physical blade order. ' ...
    'selected_rows_by_physical is retained for downstream waveform extraction compatibility.'];

save(outFile, 'HighSpeedNumbering', '-v7.3');
writetable(HighSpeedNumbering.summary, strrep(outFile, '.mat', '.csv'));

fprintf('\n=== Step04: number high-speed sensors in region ===\n');
fprintf('Selected revolutions: %d\n', numel(selectedRevIds));
fprintf('Fit degree: %d\n', matchPolyDegree);
disp(HighSpeedNumbering.summary);
fprintf('Saved: %s\n', outFile);

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function RegionSelection = build_step04_region_selection_local(S04)
RegionSelection = struct( ...
    'selection_source', 'step04_inline_manual_start', ...
    'region_id', NaN, ...
    'start_mode', 'manual', ...
    'analysis_sensors', S04.analysisSensors(:).', ...
    'target_laps', S04.targetLaps, ...
    'start_time_sec', S04.startTimeSec, ...
    'region_start_sec', S04.startTimeSec, ...
    'end_time_sec', NaN, ...
    'peak_time_sec', NaN, ...
    'dominant_order', NaN, ...
    'dominant_freq_hz', NaN, ...
    'region_plan_file', '');
end

function validate_upstream_artifacts_local(LowSpeedReference, RegionSelection, analysisSensors)
if isfield(LowSpeedReference, 'analysis_sensors')
    available = unique(LowSpeedReference.analysis_sensors(:).', 'stable');
    requested = unique(analysisSensors(:).', 'stable');
    if ~all(ismember(requested, available))
        error(['Step04 analysisSensors %s are not covered by Step01 artifact sensors %s.'], ...
            mat2str(analysisSensors), mat2str(LowSpeedReference.analysis_sensors));
    end
end
if isfield(RegionSelection, 'analysis_sensors')
    available = unique(RegionSelection.analysis_sensors(:).', 'stable');
    requested = unique(analysisSensors(:).', 'stable');
    if ~all(ismember(requested, available))
        warning(['Step04 analysisSensors %s are not fully listed in Step03 region artifact sensors %s. ', ...
            'Using the region time definition only.'], ...
            mat2str(analysisSensors), mat2str(RegionSelection.analysis_sensors));
    end
end
end

function probeCache = load_probe_cache_local(highSpeedPulseDir, sensorIds, oprTimes)
probeCache = repmat(struct('sensor_id', NaN, 'jilublade', [], 'revolution_index', []), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    probeFile = fullfile(highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
    require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
    loaded = load(probeFile, 'jilublade');
    jilublade = loaded.jilublade;
    probeCache(i).sensor_id = sid;
    probeCache(i).jilublade = jilublade;
    probeCache(i).revolution_index = map_pulses_to_revolutions_local(jilublade(:, 3), oprTimes);
end
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

function selectedRevIds = select_region_revolutions_local(probeCache, RegionSelection, P)
anchorSensor = P.sensors.anchorPreference(1);
if ~ismember(anchorSensor, [probeCache.sensor_id])
    anchorSensor = probeCache(1).sensor_id;
end
anchorProbe = probeCache([probeCache.sensor_id] == anchorSensor);
rows = find(anchorProbe.jilublade(:, 3) >= RegionSelection.start_time_sec);
if isempty(rows)
    error('No anchor pulses after %.6f s.', RegionSelection.start_time_sec);
end
revIds = unique(anchorProbe.revolution_index(rows), 'stable');
revIds = revIds(isfinite(revIds));
if numel(revIds) < RegionSelection.target_laps
    error('Only %d revolutions found after %.6f s; need %d.', ...
        numel(revIds), RegionSelection.start_time_sec, RegionSelection.target_laps);
end
selectedRevIds = revIds(1:RegionSelection.target_laps);
end

function timeWindow = build_selected_region_time_window_local(probeCache, selectedRevIds, pulsePadSec)
timeWindow = [inf, -inf];
for i = 1:numel(probeCache)
    probe = probeCache(i);
    keep = ismember(probe.revolution_index, selectedRevIds);
    rows = find(keep);
    if isempty(rows)
        continue;
    end
    timeWindow(1) = min(timeWindow(1), min(probe.jilublade(rows, 1)) - pulsePadSec);
    timeWindow(2) = max(timeWindow(2), max(probe.jilublade(rows, 2)) + pulsePadSec);
end
if ~all(isfinite(timeWindow))
    error('Could not determine selected high-speed region time window.');
end
end

function peakCacheSet = ensure_peak_cache_set_local(P, sensorIds, matchPolyDegree, sampleRateHz, pulsePadSec, autoBuildIfMissing)
peakCacheSet = cell(numel(sensorIds), 1);
fileRanges = [];
oprTimes = [];
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    cacheFile = peak_cache_file_local(P, sid);
    cache = struct();
    if exist(cacheFile, 'file') == 2
        loaded = load(cacheFile, 'HighSpeedPeakCacheSensor');
        if isfield(loaded, 'HighSpeedPeakCacheSensor') && ...
                peak_cache_matches_request_local(loaded.HighSpeedPeakCacheSensor, sid, P, matchPolyDegree, sampleRateHz, pulsePadSec)
            cache = loaded.HighSpeedPeakCacheSensor;
        end
    end
    if isempty(fieldnames(cache))
        if ~autoBuildIfMissing
            peakCacheSet = [];
            return;
        end
        fprintf('Step04: building missing CH%d peak cache ...\n', sid);
        drawnow;
        if isempty(oprTimes)
            oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
            require_file_local(oprFile, 'high-speed OPR pulse file');
            loadedOpr = load(oprFile, 'jiluOPR');
            oprTimes = loadedOpr.jiluOPR(:, 1);
        end
        if isempty(fileRanges)
            fileRanges = build_case_file_ranges_local(P.data.highSpeedDir, P.machine.oprChannel, sampleRateHz);
        end
        probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
        require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
        loadedProbe = load(probeFile, 'jilublade');
        jilublade = loadedProbe.jilublade;
        revolutionIndex = map_pulses_to_revolutions_local(jilublade(:, 3), oprTimes);
        [tRaw, vRaw] = load_full_raw_channel_local(P.data.highSpeedDir, fileRanges, sid, sampleRateHz);
        cache = build_sensor_peak_cache_local(sid, jilublade, revolutionIndex, tRaw, vRaw, P, pulsePadSec, matchPolyDegree, sampleRateHz);
        cacheDir = fileparts(cacheFile);
        if exist(cacheDir, 'dir') ~= 7
            mkdir(cacheDir);
        end
        HighSpeedPeakCacheSensor = cache; %#ok<NASGU>
        save(cacheFile, 'HighSpeedPeakCacheSensor', '-v7.3');
    end
    peakCacheSet{i} = cache;
end
end

function cacheFile = peak_cache_file_local(P, sid)
cacheFile = fullfile(P.outputDir, '04P_high_speed_peak_cache', ...
    sprintf('HighSpeedPeakCache_CH%d_20241106.mat', sid));
end

function tf = peak_cache_matches_request_local(cache, sid, P, matchPolyDegree, sampleRateHz, pulsePadSec)
tf = isstruct(cache) && ...
    isfield(cache, 'sensor_id') && isequal(cache.sensor_id, sid) && ...
    isfield(cache, 'dataset') && strcmpi(string(cache.dataset), "20241106") && ...
    isfield(cache, 'high_speed_case') && strcmpi(string(cache.high_speed_case), string(P.data.highSpeedCase)) && ...
    isfield(cache, 'sample_rate_hz') && isequal(cache.sample_rate_hz, sampleRateHz) && ...
    isfield(cache, 'match_poly_degree') && isequal(cache.match_poly_degree, matchPolyDegree) && ...
    isfield(cache, 'pulse_pad_sec') && abs(cache.pulse_pad_sec - pulsePadSec) <= eps;
end

function cache = build_sensor_peak_cache_local(sensorId, jilublade, revolutionIndex, tRaw, vRaw, P, pulsePadSec, matchPolyDegree, sampleRateHz)
nRows = size(jilublade, 1);
cache = struct();
cache.dataset = '20241106';
cache.high_speed_case = P.data.highSpeedCase;
cache.sensor_id = sensorId;
cache.sample_rate_hz = sampleRateHz;
cache.match_poly_degree = matchPolyDegree;
cache.pulse_pad_sec = pulsePadSec;
cache.row_count = nRows;
cache.row_id = (1:nRows).';
cache.t_start = jilublade(:, 1) - pulsePadSec;
cache.t_end = jilublade(:, 2) + pulsePadSec;
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
    while idxStart <= nRaw && tRaw(idxStart) < sortedStarts(k)
        idxStart = idxStart + 1;
    end
    if idxEnd < idxStart
        idxEnd = idxStart;
    end
    while idxEnd <= nRaw && tRaw(idxEnd) <= sortedEnds(k)
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
    if numel(tSeg) < max(8, matchPolyDegree + 2)
        continue;
    end
    fitInfo = fit_polynomial_peak_trace_local(tSeg, vSeg, matchPolyDegree);
    cache.raw_peak_value(rowId) = fitInfo.rawPeakValue;
    cache.raw_peak_time(rowId) = fitInfo.rawPeakTime;
    cache.fit_peak_value(rowId) = fitInfo.peakValue;
    cache.fit_peak_time(rowId) = fitInfo.peakTime;
    cache.fit_succeeded(rowId) = fitInfo.fitSucceeded;
    cache.degree_used(rowId) = fitInfo.degreeUsed;
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

function [revLib, numberingTable, summaryStruct] = build_high_speed_revolution_library_from_cache_local( ...
        peakCache, probe, selectedRevIds, refFingerprint, bladeCount)
nRev = numel(selectedRevIds);
localRows = nan(nRev, bladeCount);
peakTimes = nan(nRev, bladeCount);
localPeakValues = nan(nRev, bladeCount);
shiftScores = nan(nRev, bladeCount);
bestShifts = nan(nRev, 1);
bestScores = nan(nRev, 1);
scoreMargins = nan(nRev, 1);
rowWindows = nan(nRev, bladeCount, 2);

for k = 1:nRev
    revId = selectedRevIds(k);
    rows = find(probe.revolution_index == revId);
    [~, order] = sort(probe.jilublade(rows, 3), 'ascend');
    rows = rows(order);
    if numel(rows) ~= bladeCount
        continue;
    end
    localRows(k, :) = rows(:).';
    peakTimes(k, :) = probe.jilublade(rows, 3).';
    rowWindows(k, :, 1) = peakCache.t_start(rows);
    rowWindows(k, :, 2) = peakCache.t_end(rows);
    localPeakValues(k, :) = peakCache.fit_peak_value(rows);
    peakVec = normalize_vector_max_local(localPeakValues(k, :));
    if ~all(isfinite(peakVec))
        continue;
    end
    for shift = 0:(bladeCount - 1)
        shifted = circshift(peakVec(:), -shift);
        shiftScores(k, shift + 1) = corr_scalar_local(shifted, refFingerprint(:));
    end
    [bestScores(k), pos] = max(shiftScores(k, :));
    if ~isempty(pos) && isfinite(bestScores(k))
        bestShifts(k) = pos - 1;
        scoreMargins(k) = score_margin_local(shiftScores(k, :));
    end
end

dominantShift = choose_dominant_shift_local(bestShifts, bestScores, scoreMargins);
localToPhysical = mod((1:bladeCount) - 1 - dominantShift, bladeCount) + 1;
physicalToLocal = invert_mapping_local(localToPhysical);
physicalPeakValues = reorder_peak_matrix_local(localPeakValues, physicalToLocal);
selectedRowsByPhysical = reorder_peak_matrix_local(localRows, physicalToLocal);

revRows = struct([]);
for k = 1:nRev
    row = struct();
    row.SensorID = probe.sensor_id;
    row.RevolutionID = selectedRevIds(k);
    row.DominantShift = dominantShift;
    row.BestShift = bestShifts(k);
    row.BestScore = bestScores(k);
    row.ScoreMargin = scoreMargins(k);
    row.LocalSlotOfB1 = physicalToLocal(1);
    for b = 1:bladeCount
        row.(sprintf('B%d', b)) = physicalPeakValues(k, b);
    end
    revRows = [revRows; row]; %#ok<AGROW>
end

summaryStruct = struct();
summaryStruct.dominant_shift = dominantShift;
summaryStruct.local_to_physical = localToPhysical;
summaryStruct.physical_to_local = physicalToLocal;
summaryStruct.mean_best_score = mean(bestScores, 'omitnan');
summaryStruct.min_best_score = min(bestScores, [], 'omitnan');
summaryStruct.mean_score_margin = mean(scoreMargins, 'omitnan');

revLib = struct();
revLib.revolution_ids = selectedRevIds(:).';
revLib.local_rows = localRows;
revLib.peak_times = peakTimes;
revLib.local_peak_values = localPeakValues;
revLib.physical_peak_values = physicalPeakValues;
revLib.best_shifts = bestShifts;
revLib.best_scores = bestScores;
revLib.score_margins = scoreMargins;
revLib.shift_scores = shiftScores;
revLib.selected_rows_by_physical = selectedRowsByPhysical;
revLib.row_windows = rowWindows;
numberingTable = struct2table(revRows);
end

function result = empty_sensor_result_local(bladeCount)
result = struct( ...
    'sensor_id', NaN, ...
    'dominant_shift', NaN, ...
    'local_to_physical_blade_ids', nan(1, bladeCount), ...
    'physical_to_local_blade_ids', nan(1, bladeCount), ...
    'best_score', NaN, ...
    'score_margin', NaN, ...
    'selected_revolution_ids', [], ...
    'selected_rows_by_physical', [], ...
    'revolution_local_rows', [], ...
    'revolution_local_peak_values', [], ...
    'revolution_physical_peak_values', [], ...
    'revolution_best_shifts', [], ...
    'revolution_best_scores', [], ...
    'revolution_score_margins', [], ...
    'revolution_shift_scores', [], ...
    'revolution_peak_times', [], ...
    'revolution_row_windows', [], ...
    'numbering_table', table());
end

function fp = read_low_speed_fingerprint_local(LowSpeedReference, sid, bladeCount)
if ~isKey(LowSpeedReference.fingerprints, sid)
    error('Low-speed fingerprint missing CH%d.', sid);
end
raw = LowSpeedReference.fingerprints(sid);
if numel(raw) < bladeCount
    error('Low-speed fingerprint for CH%d has only %d values.', sid, numel(raw));
end
fp = normalize_vector_max_local(raw(1:bladeCount));
end

function [revLib, numberingTable, summaryStruct] = build_high_speed_revolution_library_local( ...
        rawSensor, probe, selectedRevIds, refFingerprint, bladeCount, pulsePadSec, matchPolyDegree)
nRev = numel(selectedRevIds);
localRows = nan(nRev, bladeCount);
peakTimes = nan(nRev, bladeCount);
localPeakValues = nan(nRev, bladeCount);
shiftScores = nan(nRev, bladeCount);
bestShifts = nan(nRev, 1);
bestScores = nan(nRev, 1);
scoreMargins = nan(nRev, 1);
rowWindows = nan(nRev, bladeCount, 2);

for k = 1:nRev
    revId = selectedRevIds(k);
    rows = find(probe.revolution_index == revId);
    [~, order] = sort(probe.jilublade(rows, 3), 'ascend');
    rows = rows(order);
    if numel(rows) ~= bladeCount
        continue;
    end
    localRows(k, :) = rows(:).';
    peakTimes(k, :) = probe.jilublade(rows, 3).';
    for slot = 1:bladeCount
        rowId = rows(slot);
        t0 = probe.jilublade(rowId, 1) - pulsePadSec;
        t1 = probe.jilublade(rowId, 2) + pulsePadSec;
        rowWindows(k, slot, :) = [t0, t1];
        keep = rawSensor.T >= t0 & rawSensor.T <= t1;
        if nnz(keep) < max(8, matchPolyDegree + 2)
            continue;
        end
        fitInfo = fit_polynomial_peak_trace_local(rawSensor.T(keep), rawSensor.V(keep), matchPolyDegree);
        if fitInfo.fitSucceeded
            localPeakValues(k, slot) = fitInfo.peakValue;
        end
    end
    peakVec = normalize_vector_max_local(localPeakValues(k, :));
    if ~all(isfinite(peakVec))
        continue;
    end
    for shift = 0:(bladeCount - 1)
        shifted = circshift(peakVec(:), -shift);
        shiftScores(k, shift + 1) = corr_scalar_local(shifted, refFingerprint(:));
    end
    [bestScores(k), pos] = max(shiftScores(k, :));
    if ~isempty(pos) && isfinite(bestScores(k))
        bestShifts(k) = pos - 1;
        scoreMargins(k) = score_margin_local(shiftScores(k, :));
    end
end

dominantShift = choose_dominant_shift_local(bestShifts, bestScores, scoreMargins);
localToPhysical = mod((1:bladeCount) - 1 - dominantShift, bladeCount) + 1;
physicalToLocal = invert_mapping_local(localToPhysical);
physicalPeakValues = reorder_peak_matrix_local(localPeakValues, physicalToLocal);
selectedRowsByPhysical = reorder_peak_matrix_local(localRows, physicalToLocal);

revRows = struct([]);
for k = 1:nRev
    row = struct();
    row.SensorID = probe.sensor_id;
    row.RevolutionID = selectedRevIds(k);
    row.DominantShift = dominantShift;
    row.BestShift = bestShifts(k);
    row.BestScore = bestScores(k);
    row.ScoreMargin = scoreMargins(k);
    row.LocalSlotOfB1 = physicalToLocal(1);
    for b = 1:bladeCount
        row.(sprintf('B%d', b)) = physicalPeakValues(k, b);
    end
    revRows = [revRows; row]; %#ok<AGROW>
end

summaryStruct = struct();
summaryStruct.dominant_shift = dominantShift;
summaryStruct.local_to_physical = localToPhysical;
summaryStruct.physical_to_local = physicalToLocal;
summaryStruct.mean_best_score = mean(bestScores, 'omitnan');
summaryStruct.min_best_score = min(bestScores, [], 'omitnan');
summaryStruct.mean_score_margin = mean(scoreMargins, 'omitnan');

revLib = struct();
revLib.revolution_ids = selectedRevIds(:).';
revLib.local_rows = localRows;
revLib.peak_times = peakTimes;
revLib.local_peak_values = localPeakValues;
revLib.physical_peak_values = physicalPeakValues;
revLib.best_shifts = bestShifts;
revLib.best_scores = bestScores;
revLib.score_margins = scoreMargins;
revLib.shift_scores = shiftScores;
revLib.selected_rows_by_physical = selectedRowsByPhysical;
revLib.row_windows = rowWindows;
numberingTable = struct2table(revRows);
end

function physicalMat = reorder_peak_matrix_local(localMat, physicalToLocal)
physicalMat = nan(size(localMat));
for bladeId = 1:numel(physicalToLocal)
    localSlot = physicalToLocal(bladeId);
    physicalMat(:, bladeId) = localMat(:, localSlot);
end
end

function dominantShift = choose_dominant_shift_local(bestShiftByRev, bestScoreByRev, scoreMarginByRev)
valid = isfinite(bestShiftByRev) & isfinite(bestScoreByRev);
if ~any(valid)
    dominantShift = 0;
    return;
end
shiftValues = bestShiftByRev(valid);
scoreValues = bestScoreByRev(valid);
marginValues = scoreMarginByRev(valid);
uniqueShifts = unique(shiftValues(:).');
rankRows = nan(numel(uniqueShifts), 4);
for i = 1:numel(uniqueShifts)
    keep = shiftValues == uniqueShifts(i);
    rankRows(i, :) = [nnz(keep), mean(scoreValues(keep), 'omitnan'), ...
        mean(marginValues(keep), 'omitnan'), -uniqueShifts(i)];
end
[~, order] = sortrows(rankRows, [-1 -2 -3 -4]);
dominantShift = uniqueShifts(order(1));
end

function fitInfo = fit_polynomial_peak_trace_local(t, v, degree)
fitInfo = struct( ...
    'peakValue', NaN, ...
    'peakTime', NaN, ...
    'rawPeakValue', NaN, ...
    'rawPeakTime', NaN, ...
    't', [], ...
    'v', [], ...
    'tFine', [], ...
    'vFine', [], ...
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
fitInfo.t = t;
fitInfo.v = v;
realDegree = min(degree, numel(t) - 1);
fitInfo.degreeUsed = realDegree;
[coef, ~, mu] = polyfit(t, v, realDegree);
tFine = linspace(min(t), max(t), 200).';
vFine = polyval(coef, tFine, [], mu);
[peakValue, peakIdx] = max(vFine);
fitInfo.peakValue = peakValue;
fitInfo.peakTime = tFine(peakIdx);
fitInfo.tFine = tFine;
fitInfo.vFine = vFine;
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

function raw = load_raw_subset_local(caseDir, fileRanges, sensorIds, sampleRateHz, timeWindow)
raw(max(sensorIds)) = struct('T', [], 'V', []);
useFiles = fileRanges([fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2));
for k = 1:numel(useFiles)
    fileId = useFiles(k).file_id;
    offset = useFiles(k).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, sampleRateHz);
        tGlobal = tLocal(:) + offset;
        keep = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        raw(sid).T = [raw(sid).T; tGlobal(keep)]; %#ok<AGROW>
        raw(sid).V = [raw(sid).V; vLocal(keep)]; %#ok<AGROW>
    end
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

function y = normalize_vector_max_local(x)
x = x(:);
mx = max(x);
if ~isfinite(mx) || mx <= eps
    y = nan(size(x));
else
    y = x / mx;
end
end

function value = corr_scalar_local(a, b)
good = isfinite(a) & isfinite(b);
if nnz(good) < 3
    value = NaN;
    return;
end
a = a(good);
b = b(good);
if std(a) <= eps || std(b) <= eps
    value = NaN;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function margin = score_margin_local(scoreRow)
sortedScores = sort(scoreRow, 'descend', 'MissingPlacement', 'last');
if numel(sortedScores) >= 2 && isfinite(sortedScores(1)) && isfinite(sortedScores(2))
    margin = sortedScores(1) - sortedScores(2);
else
    margin = NaN;
end
end

function invMap = invert_mapping_local(map)
invMap = nan(size(map));
for i = 1:numel(map)
    v = map(i);
    if isfinite(v) && v >= 1 && v <= numel(map)
        invMap(v) = i;
    end
end
end

function P = apply_step04_local_options_local(P, S04)
P.sensors.analysis = S04.analysisSensors(:).';
P.region.startTimeSec = S04.startTimeSec;
P.waveform.pulsePadSec = S04.pulsePadSec;
sensorTag = step04_sensor_tag_local(P.sensors.analysis);
timeLabel = step04_time_label_local(P.region.startTimeSec);
P.files.regionSelection = fullfile(P.outputDir, '03_region_selection', ...
    sprintf('HighSpeedRegion_%s_20241106.mat', timeLabel));
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
end

function tag = step04_sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = step04_time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end
