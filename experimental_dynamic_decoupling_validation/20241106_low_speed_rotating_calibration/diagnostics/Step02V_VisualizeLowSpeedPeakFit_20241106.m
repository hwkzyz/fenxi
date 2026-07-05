%% Step02V_VisualizeLowSpeedPeakFit_20241106
% Visualize the low-speed seven-order peak fitting for one representative
% OPR-defined revolution per sensor. This script is intentionally separate
% from Step02 so the main numbering path stays lightweight.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

%% Parameters to tune
analysisSensors = [2 3 5 7];
bladeCount = 6;
representativeMode = 'reference_revolution';   % 'reference_revolution' | 'best_score'
peakFitDegreeOverride = [];
peakFitRefinePointsOverride = [];
saveFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
numberingFile = fullfile(routeDir, 'output', 'new_flow', '02_low_speed_numbering', ...
    ['S', sprintf('%d', analysisSensors)], 'LowSpeedNumbering_20241106.mat');
figureDir = fullfile(routeDir, 'output', 'new_flow', 'figures', '02_low_speed_numbering');

if exist(numberingFile, 'file') ~= 2
    error(['Missing Step02 artifact:\n  %s\n\n' ...
        'Run Step02_NumberLowSpeedSensors_20241106 first.'], numberingFile);
end

loaded = load(numberingFile, 'LowSpeedNumbering', 'LowSpeedSummary');
if ~isfield(loaded, 'LowSpeedNumbering') || ~isfield(loaded, 'LowSpeedSummary')
    error('Step02 file must contain LowSpeedNumbering and LowSpeedSummary:\n  %s', numberingFile);
end
LowSpeedNumbering = loaded.LowSpeedNumbering;
LowSpeedSummary = loaded.LowSpeedSummary;

lowSpeedFingerprintFile = fullfile(routeDir, 'output', 'new_flow', '01_low_speed_reference', ...
    'LowSpeedReferenceFingerprint_20241106.mat');
loadedRef = load(lowSpeedFingerprintFile, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;
Sensor_Config = load_source_sensor_config_local(LowSpeedReference, lowSpeedFingerprintFile);

P = NewFlow_Config_20241106();
datasetRoot = P.data.datasetRoot;
initialTrimPoints = 20000;

if isempty(peakFitDegreeOverride)
    peakFitDegree = read_or_default_local(Sensor_Config, 'Peak_Fit_Degree', 7);
else
    peakFitDegree = peakFitDegreeOverride;
end
if isempty(peakFitRefinePointsOverride)
    peakFitRefinePoints = read_or_default_local(Sensor_Config, 'Peak_Fit_Refine_Points', 200);
else
    peakFitRefinePoints = peakFitRefinePointsOverride;
end

validate_step02v_inputs_local(LowSpeedNumbering, LowSpeedSummary, LowSpeedReference, analysisSensors, bladeCount);

visualize_low_speed_peak_fit_audit_local( ...
    LowSpeedNumbering, LowSpeedSummary, LowSpeedReference, Sensor_Config, ...
    datasetRoot, analysisSensors, bladeCount, initialTrimPoints, ...
    peakFitDegree, peakFitRefinePoints, representativeMode, saveFigures, figureDir);

fprintf('\n=== Step02V: visualize low-speed peak fit ===\n');
fprintf('Source numbering: %s\n', numberingFile);
fprintf('Representative mode: %s\n', representativeMode);
fprintf('Fit degree: %d\n', peakFitDegree);
fprintf('Figures: %s\n', figureDir);

function value = read_or_default_local(S, fieldName, defaultValue)
if isfield(S, fieldName) && isfinite(S.(fieldName))
    value = double(S.(fieldName));
else
    value = defaultValue;
end
end

function validate_step02v_inputs_local( ...
        LowSpeedNumbering, LowSpeedSummary, LowSpeedReference, sensorIds, bladeCount)
if isempty(LowSpeedNumbering) || isempty(LowSpeedSummary)
    error('Step02 outputs are empty.');
end
if ~isa(LowSpeedReference.revolution_pulse_times, 'containers.Map') || ...
        ~isa(LowSpeedReference.sensor_physical_to_local, 'containers.Map')
    error('LowSpeedReference is missing revolution pulse timing maps needed by Step02V.');
end
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    if ~any(LowSpeedSummary.SensorID == sid)
        error('LowSpeedSummary does not contain CH%d.', sid);
    end
    rows = LowSpeedNumbering(LowSpeedNumbering.SensorID == sid, :);
    if height(rows) < 1
        error('LowSpeedNumbering does not contain any revolution rows for CH%d.', sid);
    end
    for b = 1:bladeCount
        if ~ismember(sprintf('B%d', b), LowSpeedNumbering.Properties.VariableNames)
            error('LowSpeedNumbering is missing column B%d.', b);
        end
    end
end
end

function Sensor_Config = load_source_sensor_config_local(LowSpeedReference, lowSpeedFingerprintFile)
if ~isfield(LowSpeedReference, 'source_file') || isempty(LowSpeedReference.source_file)
    error(['LowSpeedReference.source_file is missing.\n' ...
        'Cannot trace low-speed fit audit source from:\n  %s'], lowSpeedFingerprintFile);
end
sourceFile = LowSpeedReference.source_file;
if exist(sourceFile, 'file') ~= 2
    error('Low-speed source Sensor_Config file not found:\n  %s', sourceFile);
end
loaded = load(sourceFile, 'Sensor_Config');
if ~isfield(loaded, 'Sensor_Config')
    error('Source file does not contain Sensor_Config:\n  %s', sourceFile);
end
Sensor_Config = loaded.Sensor_Config;
end

function visualize_low_speed_peak_fit_audit_local( ...
        LowSpeedNumbering, LowSpeedSummary, LowSpeedReference, Sensor_Config, ...
        datasetRoot, analysisSensors, bladeCount, initialTrimPoints, ...
        peakFitDegree, peakFitRefinePoints, representativeMode, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

caseDir = fullfile(datasetRoot, Sensor_Config.ReferenceCase);
if exist(caseDir, 'dir') ~= 7
    error('Low-speed case directory not found:\n  %s', caseDir);
end

thresholdMap = build_sensor_threshold_map_local(Sensor_Config);
bladeColors = lines(bladeCount);
pulseCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    pulseCache(sid) = extract_low_speed_sensor_pulses_local( ...
        caseDir, sid, thresholdMap(sid), Sensor_Config.Gap_Points, Sensor_Config.Pinlv, ...
        initialTrimPoints, peakFitDegree, peakFitRefinePoints);
end

fig = figure('Name', 'Step02V low-speed peak-fit audit', 'Color', 'w', ...
    'Position', [80, 60, 1800, 920], 'NumberTitle', 'off');
tiledlayout(fig, numel(analysisSensors), bladeCount, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    pulseInfo = pulseCache(sid);
    [repRow, repRevId] = choose_representative_revolution_local( ...
        LowSpeedNumbering, LowSpeedSummary, sid, representativeMode, LowSpeedReference.reference_revolution_id);
    pulseTimeMat = LowSpeedReference.revolution_pulse_times(sid);
    physicalToLocal = force_row_vector_local(LowSpeedReference.sensor_physical_to_local(sid));
    localSlots = nan(1, bladeCount);
    pulseIndices = nan(1, bladeCount);
    for bladeId = 1:bladeCount
        localSlot = physicalToLocal(bladeId);
        localSlots(bladeId) = localSlot;
        tPeak = pulseTimeMat(repRow, localSlot);
        pulseIndices(bladeId) = find_closest_pulse_index_local(pulseInfo.arrival_times, tPeak);
    end

    for bladeId = 1:bladeCount
        nexttile;
        render_low_speed_single_pulse_fit_local( ...
            pulseInfo, pulseIndices(bladeId), sid, bladeId, repRevId, ...
            localSlots(bladeId), bladeColors(bladeId, :));
    end
end

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step02V_LowSpeedPeakFitAudit_20241106.png'), ...
        'Resolution', 300);
end
end

function [rowIdx, revId] = choose_representative_revolution_local( ...
        LowSpeedNumbering, LowSpeedSummary, sid, representativeMode, referenceRevId)
rows = LowSpeedNumbering(LowSpeedNumbering.SensorID == sid, :);
modeText = lower(strtrim(representativeMode));
switch modeText
    case 'reference_revolution'
        rowIdx = find(rows.RevolutionID == referenceRevId, 1, 'first');
        if isempty(rowIdx)
            [~, rowIdx] = max(rows.BestScore);
        end
    case 'best_score'
        [~, rowIdx] = max(rows.BestScore);
    otherwise
        error('Unsupported representativeMode: %s', representativeMode);
end
if isempty(rowIdx)
    error('Could not choose representative revolution for CH%d.', sid);
end
revId = rows.RevolutionID(rowIdx);

summaryRow = LowSpeedSummary(LowSpeedSummary.SensorID == sid, :);
if height(summaryRow) == 1 && isfinite(summaryRow.ReferenceRevolutionID) && ...
        strcmpi(modeText, 'reference_revolution') && revId ~= summaryRow.ReferenceRevolutionID
    warning('CH%d does not contain reference revolution %d; Step02V falls back to best-score revolution %d.', ...
        sid, summaryRow.ReferenceRevolutionID, revId);
end
end

function idx = find_closest_pulse_index_local(arrivalTimes, tPeak)
arrivalTimes = arrivalTimes(:);
[~, idx] = min(abs(arrivalTimes - tPeak));
end

function thresholdMap = build_sensor_threshold_map_local(Sensor_Config)
thresholdMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
ids = Sensor_Config.Sensor_IDs(:).';
if isa(Sensor_Config.Sensor_Thresholds, 'containers.Map')
    sourceMap = Sensor_Config.Sensor_Thresholds;
    for i = 1:numel(ids)
        sid = double(ids(i));
        if isKey(sourceMap, sid)
            thresholdMap(sid) = sourceMap(sid);
        end
    end
else
    vals = Sensor_Config.Sensor_Thresholds(:).';
    for i = 1:min(numel(ids), numel(vals))
        thresholdMap(double(ids(i))) = vals(i);
    end
end
end

function pulseInfo = extract_low_speed_sensor_pulses_local( ...
        caseDir, sid, threshold, gapPoints, sampleRateHz, initialTrimPoints, fitDegree, refinePoints)
fileIds = list_channel_file_ids_local(caseDir, sid);
tail = struct('T', [], 'V', []);
timeOffset = 0;
pulseInfo = struct( ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'start_times', [], ...
    'end_times', [], ...
    'pulses', struct('t', {}, 'v', {}, 'start_time', {}, 'end_time', {}, 'fit', {}));

for i = 1:numel(fileIds)
    [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileIds(i), sampleRateHz);
    if isempty(tLocal)
        continue;
    end
    if i == 1 && initialTrimPoints > 0 && numel(tLocal) > initialTrimPoints
        tLocal = tLocal((initialTrimPoints + 1):end);
        vLocal = vLocal((initialTrimPoints + 1):end);
    end
    tGlobal = tLocal + timeOffset;
    [pulsePart, tail] = extract_sensor_pulses_with_fit_local( ...
        tGlobal, vLocal, threshold, gapPoints, fitDegree, refinePoints, tail);
    pulseInfo.arrival_times = [pulseInfo.arrival_times; pulsePart.arrival_times]; %#ok<AGROW>
    pulseInfo.peak_features = [pulseInfo.peak_features; pulsePart.peak_features]; %#ok<AGROW>
    pulseInfo.start_times = [pulseInfo.start_times; pulsePart.start_times]; %#ok<AGROW>
    pulseInfo.end_times = [pulseInfo.end_times; pulsePart.end_times]; %#ok<AGROW>
    pulseInfo.pulses = [pulseInfo.pulses; pulsePart.pulses]; %#ok<AGROW>
    timeOffset = tGlobal(end) + 1 / sampleRateHz;
end
end

function [pulseInfo, nextTail] = extract_sensor_pulses_with_fit_local( ...
        t, v, threshold, gapPoints, fitDegree, refinePoints, tail)
if isempty(tail) || ~isfield(tail, 'T') || isempty(tail.T)
    tAll = t(:);
    vAll = v(:);
else
    tAll = [tail.T(:); t(:)];
    vAll = [tail.V(:); v(:)];
end

smoothV = smooth_signal_local(vAll);
idxAbove = find(smoothV > threshold);
pulseInfo = struct( ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'start_times', [], ...
    'end_times', [], ...
    'pulses', struct('t', {}, 'v', {}, 'start_time', {}, 'end_time', {}, 'fit', {}));

if isempty(idxAbove)
    nextTail = prepare_tail_local(tAll, vAll, []);
    return;
end

[segStart, segEnd] = split_segments_local(idxAbove, gapPoints);
if isempty(segStart)
    nextTail = prepare_tail_local(tAll, vAll, []);
    return;
end

keepSeg = true(numel(segStart), 1);
if segEnd(end) >= numel(tAll) - 1
    keepSeg(end) = false;
end

for i = find(keepSeg).'
    idxRange = segStart(i):segEnd(i);
    tSeg = tAll(idxRange);
    vSeg = vAll(idxRange);
    smoothSeg = smoothV(idxRange);
    fitInfo = fit_pulse_peak_trace_local(tSeg, vSeg, smoothSeg, threshold, fitDegree, refinePoints);
    k = numel(pulseInfo.arrival_times) + 1;
    pulseInfo.arrival_times(k, 1) = fitInfo.arrivalTime;
    pulseInfo.peak_features(k, 1) = fitInfo.peakValue;
    pulseInfo.start_times(k, 1) = tSeg(1);
    pulseInfo.end_times(k, 1) = tSeg(end);
    pulseInfo.pulses(k, 1) = struct( ...
        't', tSeg(:), ...
        'v', vSeg(:), ...
        'start_time', tSeg(1), ...
        'end_time', tSeg(end), ...
        'fit', fitInfo);
end

nextTail = prepare_tail_local(tAll, vAll, segStart);
end

function nextTail = prepare_tail_local(tAll, vAll, segStart)
nextTail = struct('T', [], 'V', []);
if isempty(tAll)
    return;
end
if isempty(segStart)
    keepFrom = max(1, numel(tAll) - 4 * max(1, round(numel(tAll) / 50)));
else
    keepFrom = max(1, segStart(end) - 4 * max(1, round(numel(tAll) / 50)));
end
nextTail.T = tAll(keepFrom:end);
nextTail.V = vAll(keepFrom:end);
end

function fitInfo = fit_pulse_peak_trace_local(tSeg, vSeg, smoothSeg, threshold, fitDegree, refinePoints)
fitInfo = struct( ...
    'arrivalTime', NaN, ...
    'peakValue', NaN, ...
    'peakTime', NaN, ...
    'rawPeakValue', NaN, ...
    'rawPeakTime', NaN, ...
    't', tSeg(:), ...
    'v', vSeg(:), ...
    'tFine', [], ...
    'vFine', [], ...
    'degreeUsed', NaN, ...
    'fitSucceeded', false);

validMask = smoothSeg > threshold;
if sum(validMask) < 4
    [fitInfo.rawPeakValue, rawIdx] = max(vSeg);
    fitInfo.rawPeakTime = tSeg(rawIdx);
    fitInfo.peakValue = fitInfo.rawPeakValue;
    fitInfo.peakTime = fitInfo.rawPeakTime;
    fitInfo.arrivalTime = fitInfo.rawPeakTime;
    return;
end

tFit = tSeg(validMask);
vFit = vSeg(validMask);
[fitInfo.rawPeakValue, rawIdx] = max(vSeg);
fitInfo.rawPeakTime = tSeg(rawIdx);
mu = mean(tFit);
uniqueCount = numel(unique(tFit));
order = min(fitDegree, uniqueCount - 1);
fitInfo.degreeUsed = order;
if order < 1
    fitInfo.peakValue = fitInfo.rawPeakValue;
    fitInfo.peakTime = fitInfo.rawPeakTime;
    fitInfo.arrivalTime = fitInfo.rawPeakTime;
    return;
end

try
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    p = polyfit(tFit - mu, vFit, order);
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    dt = linspace(min(tFit) - mu, max(tFit) - mu, refinePoints).';
    vf = polyval(p, dt);
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;
    if sum(weights) > 0
        fitInfo.arrivalTime = sum(tf .* weights) / sum(weights);
    else
        fitInfo.arrivalTime = mean(tf);
    end
    [fitInfo.peakValue, idxMax] = max(vf);
    fitInfo.peakTime = tf(idxMax);
    fitInfo.tFine = tf;
    fitInfo.vFine = vf;
    fitInfo.fitSucceeded = true;
catch
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    fitInfo.peakValue = fitInfo.rawPeakValue;
    fitInfo.peakTime = fitInfo.rawPeakTime;
    fitInfo.arrivalTime = fitInfo.rawPeakTime;
end
end

function render_low_speed_single_pulse_fit_local( ...
        pulseInfo, pulseIdx, sid, bladeId, revId, localSlot, bladeColor)
if pulseIdx < 1 || pulseIdx > numel(pulseInfo.pulses)
    axis off;
    title(sprintf('CH%d R%d B%d missing', sid, revId, bladeId));
    return;
end

P = pulseInfo.pulses(pulseIdx);
fitInfo = P.fit;
t0 = P.start_time;
tMs = (P.t - t0) * 1000;
plot(tMs, P.v, 'k-', 'LineWidth', 0.9);
hold on; box on; grid on;
if fitInfo.fitSucceeded
    plot((fitInfo.tFine - t0) * 1000, fitInfo.vFine, '-', ...
        'Color', [0.10 0.45 0.90], 'LineWidth', 1.1);
    plot((fitInfo.peakTime - t0) * 1000, fitInfo.peakValue, 'o', ...
        'MarkerSize', 5.5, 'MarkerFaceColor', bladeColor, 'MarkerEdgeColor', bladeColor);
end
plot((fitInfo.rawPeakTime - t0) * 1000, fitInfo.rawPeakValue, 'x', ...
    'MarkerSize', 5.5, 'LineWidth', 1.0, 'Color', [0.55 0.55 0.55]);
xline(0, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
xline((P.end_time - t0) * 1000, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
title(sprintf('CH%d R%d B%d(L%d) fit %.3f raw %.3f', ...
    sid, revId, bladeId, localSlot, fitInfo.peakValue, fitInfo.rawPeakValue), 'FontSize', 8);
xlabel('Pulse window time (ms)');
ylabel('Voltage');
yl = ylim;
text(0.02 * max(tMs), yl(2) - 0.08 * range(yl), ...
    sprintf('pulse%d deg%d', pulseIdx, fitInfo.degreeUsed), ...
    'FontSize', 7, 'Color', [0.25 0.25 0.25], 'VerticalAlignment', 'top');
end

function smoothV = smooth_signal_local(v)
if numel(v) >= 21
    smoothV = sgolayfilt(v, 3, 21);
elseif numel(v) >= 5
    span = max(3, 2 * floor(numel(v) / 4) + 1);
    smoothV = smoothdata(v, 'movmean', span);
else
    smoothV = v;
end
end

function [segStart, segEnd] = split_segments_local(idxAbove, gapPoints)
if isempty(idxAbove)
    segStart = zeros(0, 1);
    segEnd = zeros(0, 1);
    return;
end
jumps = find(diff(idxAbove) > gapPoints);
segStart = [idxAbove(1); idxAbove(jumps + 1)];
segEnd = [idxAbove(jumps); idxAbove(end)];
end

function fileIds = list_channel_file_ids_local(caseDir, sid)
files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.dat', sid)));
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.txt', sid)));
end
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('4-%d-*.mat', sid)));
end
if isempty(files)
    error('No raw files found for CH%d in %s.', sid, caseDir);
end
fileIds = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, 'Data_(\d+)', 'tokens', 'once');
    if isempty(token)
        token = regexp(files(i).name, '4-\d+-(\d+)\.mat', 'tokens', 'once');
    end
    fileIds(i) = str2double(token{1});
end
fileIds = sort(fileIds);
end

function [t, v] = load_raw_channel_local(caseDir, sid, fileId, sampleRateHz)
datFile = fullfile(caseDir, sprintf('Probe%d_Data_%d.dat', sid, fileId));
txtFile = fullfile(caseDir, sprintf('Probe%d_Data_%d.txt', sid, fileId));
matFile = fullfile(caseDir, sprintf('4-%d-%d.mat', sid, fileId));

if exist(datFile, 'file') == 2
    A = load(datFile);
    v = A(:, 2);
elseif exist(txtFile, 'file') == 2
    A = load(txtFile);
    v = A(:, 2);
elseif exist(matFile, 'file') == 2
    loaded = load(matFile);
    fields = fieldnames(loaded);
    raw = loaded.(fields{1});
    raw(raw(:, 1) == 0, :) = [];
    t = raw(:, 1) / sampleRateHz;
    v = raw(:, 2);
    return;
else
    error('Could not find raw file for CH%d file %d in %s.', sid, fileId, caseDir);
end
t = (0:(numel(v) - 1)).' / sampleRateHz;
end

function x = force_row_vector_local(x)
x = x(:).';
end

