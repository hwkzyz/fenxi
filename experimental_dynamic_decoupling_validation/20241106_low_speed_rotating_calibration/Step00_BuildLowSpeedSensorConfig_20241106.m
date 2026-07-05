%% Step00_BuildLowSpeedSensorConfig_20241106
% Build the low-speed Sensor_Config artifact inside the current main folder.
% This is the source step for the new Step01 low-speed reference fingerprint.

clear; close all; clc;

%% Parameters to tune
datasetRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
lowSpeedCase = '900';

sensorIds = [2 3 4 5 6 7];
capacitanceIds = [5 7];
eddyCurrentIds = [2 3];
unusedIds = [];
oprId = 1;

bladesNum = 6;
oprPulsesPerRev = 1;
pinlv = 5e6;
rTipMM = 65.0;
initialTrimPoints = 20000;
gapPoints = 10000;
oprThreshold = 2.0;

referenceSensorId = 5;
referenceSearchRevs = 3;
maxAlignmentSearchStarts = 24;
maxLapsProcess = 5000;

sensorThresholdIds = [1 2 3 4 5 6 7];
sensorThresholdVals = [2.0 -0.2 -0.2 2.0 1.0 2.0 1.0];
sensorPeakFitDegree = 7;
sensorPeakFitRefinePoints = 200;

showPlots = true;
saveFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
caseDir = fullfile(datasetRoot, lowSpeedCase);
referenceOutputDir = fullfile(routeDir, 'output', 'reference', [lowSpeedCase, '_reference']);
sensorConfigFile = fullfile(referenceOutputDir, 'Sensor_Config_20241106.mat');
sensorConfigSummaryFile = fullfile(referenceOutputDir, 'Sensor_Config_Summary_20241106.csv');
figureDir = fullfile(routeDir, 'output', 'figures', '00_low_speed_sensor_config');

if exist(referenceOutputDir, 'dir') ~= 7
    mkdir(referenceOutputDir);
end
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end
if exist(caseDir, 'dir') ~= 7
    error('Low-speed data folder not found:\n  %s', caseDir);
end

sensorThresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:numel(sensorThresholdIds)
    sensorThresholds(sensorThresholdIds(i)) = sensorThresholdVals(i);
end

cfg = struct();
cfg.dataset_root = datasetRoot;
cfg.low_speed_case = lowSpeedCase;
cfg.reference_output_dir = referenceOutputDir;
cfg.sensor_ids = sensorIds;
cfg.capacitance_ids = capacitanceIds;
cfg.eddy_current_ids = eddyCurrentIds;
cfg.unused_ids = unusedIds;
cfg.opr_id = oprId;
cfg.blades_num = bladesNum;
cfg.opr_pulses_per_rev = oprPulsesPerRev;
cfg.pinlv = pinlv;
cfg.r_tip_mm = rTipMM;
cfg.initial_trim_points = initialTrimPoints;
cfg.gap_points = gapPoints;
cfg.opr_threshold = oprThreshold;
cfg.reference_sensor_id = referenceSensorId;
cfg.reference_search_revs = referenceSearchRevs;
cfg.max_alignment_search_starts = maxAlignmentSearchStarts;
cfg.max_laps_process = maxLapsProcess;
cfg.sensor_thresholds = sensorThresholds;
cfg.sensor_peak_fit_degree = sensorPeakFitDegree;
cfg.sensor_peak_fit_refine_points = sensorPeakFitRefinePoints;

fprintf('\n=== Step00: build low-speed Sensor_Config ===\n');
fprintf('Low-speed data folder: %s\n', caseDir);
fprintf('Sensors: %s\n', mat2str(sensorIds));
fprintf('Output: %s\n', sensorConfigFile);

caseData = extract_btt_features_local(caseDir, cfg);

revolutionIdsMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionLocalPeaksMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionPulseIndicesMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionPulseTimesMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionBestShiftMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionBestScoreMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionScoreMarginMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionShiftScoresMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
revolutionPhysicalPeaksMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
sensorDominantShiftMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
sensorLocalToPhysicalMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
sensorPhysicalToLocalMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
fingerprints = containers.Map('KeyType', 'double', 'ValueType', 'any');
targetIndices = containers.Map('KeyType', 'double', 'ValueType', 'double');

revolutionLibraries = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    revolutionLibraries(sid) = build_sensor_revolution_library_local(caseData, sid, cfg);
end

refLibrary = revolutionLibraries(referenceSensorId);
if isempty(refLibrary.revolution_ids)
    error('Reference sensor CH%d does not have any valid OPR-defined six-pulse revolutions.', referenceSensorId);
end
[referenceRowIdx, referenceRevolutionId, refFp] = choose_reference_revolution_local(refLibrary, cfg);

summaryRows = repmat(struct( ...
    'sensor_id', NaN, ...
    'dominant_shift', NaN, ...
    'valid_revolution_count', NaN, ...
    'mean_best_score', NaN, ...
    'reference_revolution_id', NaN), numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    lib = revolutionLibraries(sid);
    if isempty(lib.revolution_ids)
        error('CH%d does not have valid OPR-defined six-pulse revolutions.', sid);
    end

    match = match_revolution_library_local(lib.local_peak_vectors, refFp, bladesNum);
    physicalPeaks = reorder_peak_matrix_local(lib.local_peak_vectors, match.physical_to_local);
    fingerprints(sid) = mean(physicalPeaks, 1, 'omitnan');

    rowForTarget = find(lib.revolution_ids == referenceRevolutionId, 1, 'first');
    if isempty(rowForTarget)
        rowForTarget = referenceRowIdx;
        rowForTarget = min(max(rowForTarget, 1), size(lib.pulse_indices, 1));
    end
    b1LocalSlot = match.physical_to_local(1);
    targetIndices(sid) = lib.pulse_indices(rowForTarget, b1LocalSlot);

    revolutionIdsMap(sid) = lib.revolution_ids(:).';
    revolutionLocalPeaksMap(sid) = lib.local_peak_vectors;
    revolutionPulseIndicesMap(sid) = lib.pulse_indices;
    revolutionPulseTimesMap(sid) = lib.pulse_times;
    revolutionBestShiftMap(sid) = match.best_shift_by_revolution(:).';
    revolutionBestScoreMap(sid) = match.best_score_by_revolution(:).';
    revolutionScoreMarginMap(sid) = match.score_margin_by_revolution(:).';
    revolutionShiftScoresMap(sid) = match.shift_scores_by_revolution;
    revolutionPhysicalPeaksMap(sid) = physicalPeaks;
    sensorDominantShiftMap(sid) = match.dominant_shift;
    sensorLocalToPhysicalMap(sid) = match.local_to_physical;
    sensorPhysicalToLocalMap(sid) = match.physical_to_local;

    summaryRows(i).sensor_id = sid;
    summaryRows(i).dominant_shift = match.dominant_shift;
    summaryRows(i).valid_revolution_count = numel(lib.revolution_ids);
    summaryRows(i).mean_best_score = mean(match.best_score_by_revolution, 'omitnan');
    summaryRows(i).reference_revolution_id = referenceRevolutionId;
end

stdAnglesStartEdge = nan(max([sensorIds, oprId]), bladesNum);
if numel(caseData.opr_times) > oprPulsesPerRev
    spdT = caseData.opr_times(1:(end - oprPulsesPerRev));
    spdV = 360 ./ (caseData.opr_times((oprPulsesPerRev + 1):end) - caseData.opr_times(1:(end - oprPulsesPerRev)));
    FomegaDeg = griddedInterpolant(spdT, spdV, 'linear', 'nearest');

    for sid = sensorIds
        measTimes = caseData.channels(sid).arrival_times(:);
        lib = revolutionLibraries(sid);
        physicalToLocal = sensorPhysicalToLocalMap(sid);
        pulseIdxMat = lib.pulse_indices;
        for bladeId = 1:bladesNum
            localSlot = physicalToLocal(bladeId);
            pulseIdxList = pulseIdxMat(:, localSlot);
            pulseIdxList = pulseIdxList(isfinite(pulseIdxList));
            relAngles = [];
            for j = 1:min(numel(pulseIdxList), maxLapsProcess)
                pulseIdx = pulseIdxList(j);
                tMeas = measTimes(pulseIdx);
                idxPrevOpr = find(caseData.opr_times < tMeas, 1, 'last');
                if isempty(idxPrevOpr)
                    continue;
                end
                tRef = caseData.opr_times(idxPrevOpr);
                tGrid = linspace(tRef, tMeas, 10);
                relAngles(end + 1) = trapz(tGrid, FomegaDeg(tGrid)); %#ok<AGROW>
            end
            if ~isempty(relAngles)
                stdAnglesStartEdge(sid, bladeId) = mean(relAngles);
            end
        end
    end
end

oprReference = build_opr_center_reference_for_case_local(caseDir, cfg);
stdAnglesOprCenter = stdAnglesStartEdge - oprReference.phase_shift_deg;

Sensor_Config = struct();
Sensor_Config.ReferenceCase = lowSpeedCase;
Sensor_Config.Sensor_IDs = sensorIds;
Sensor_Config.Capacitance_IDs = capacitanceIds;
Sensor_Config.Eddy_Current_IDs = eddyCurrentIds;
Sensor_Config.OPR_ID = oprId;
Sensor_Config.Blades_Num = bladesNum;
Sensor_Config.Pinlv = pinlv;
Sensor_Config.Gap_Points = gapPoints;
Sensor_Config.OPR_Threshold = oprThreshold;
Sensor_Config.Sensor_Thresholds = sensorThresholds;
Sensor_Config.Fingerprints = fingerprints;
Sensor_Config.Target_Indices = targetIndices;
Sensor_Config.ReferenceSensorID = referenceSensorId;
Sensor_Config.ReferenceRevolutionID = referenceRevolutionId;
Sensor_Config.ReferenceFingerprint = refFp(:).';
Sensor_Config.Revolution_IDs = revolutionIdsMap;
Sensor_Config.Revolution_Local_Peaks = revolutionLocalPeaksMap;
Sensor_Config.Revolution_Local_Pulse_Indices = revolutionPulseIndicesMap;
Sensor_Config.Revolution_Pulse_Times = revolutionPulseTimesMap;
Sensor_Config.Revolution_Best_Shifts = revolutionBestShiftMap;
Sensor_Config.Revolution_Best_Scores = revolutionBestScoreMap;
Sensor_Config.Revolution_Score_Margins = revolutionScoreMarginMap;
Sensor_Config.Revolution_Shift_Scores = revolutionShiftScoresMap;
Sensor_Config.Revolution_Physical_Peaks = revolutionPhysicalPeaksMap;
Sensor_Config.Sensor_Dominant_Shifts = sensorDominantShiftMap;
Sensor_Config.Sensor_Local_To_Physical = sensorLocalToPhysicalMap;
Sensor_Config.Sensor_Physical_To_Local = sensorPhysicalToLocalMap;
Sensor_Config.Standard_Relative_Angles = stdAnglesOprCenter;
Sensor_Config.Standard_Relative_Angles_OPRCenter = stdAnglesOprCenter;
Sensor_Config.Standard_Relative_Angles_StartEdge = stdAnglesStartEdge;
Sensor_Config.Standard_Relative_Angles_Reference = 'opr_pulse_center';
Sensor_Config.Standard_Relative_Angles_StartEdge_Reference = 'opr_threshold_start_edge';
Sensor_Config.OPRReference = oprReference;
Sensor_Config.Peak_Fit_Degree = cfg.sensor_peak_fit_degree;
Sensor_Config.Peak_Fit_Refine_Points = cfg.sensor_peak_fit_refine_points;

summaryTable = struct2table(summaryRows);
save(sensorConfigFile, 'Sensor_Config');
writetable(summaryTable, sensorConfigSummaryFile);

fprintf('Saved Sensor_Config to:\n  %s\n', sensorConfigFile);
disp(summaryTable);
print_sensor_angle_summary_local(Sensor_Config, cfg);

if showPlots
    visualize_sensor_config_local(caseData, summaryTable, Sensor_Config, cfg, saveFigures, figureDir);
end

function caseData = extract_btt_features_local(caseDir, cfg)
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR files found in %s.', caseDir);
end

maxChannelId = max([cfg.sensor_ids, cfg.opr_id]);
channels(maxChannelId) = struct( ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'start_times', [], ...
    'end_times', []);
tailData(maxChannelId) = struct('t', [], 'v', []);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);

oprTimes = [];
lastEndTime = [];

for iFile = 1:numel(fileIds)
    fileId = fileIds(iFile);
    [tOprLocal, vOprLocal] = load_raw_channel_local(caseDir, cfg.opr_id, fileId, cfg.pinlv);
    if isempty(lastEndTime)
        offset = 0;
    else
        offset = lastEndTime + 1 / cfg.pinlv - tOprLocal(1);
    end

    tOprGlobal = tOprLocal + offset;
    lastEndTime = tOprGlobal(end);
    fileRanges(iFile).file_id = fileId;
    fileRanges(iFile).offset = offset;
    fileRanges(iFile).t_start = tOprGlobal(1);
    fileRanges(iFile).t_end = tOprGlobal(end);

    if iFile == 1 && cfg.initial_trim_points > 0
        keepIdx = (cfg.initial_trim_points + 1):numel(tOprGlobal);
        tOprGlobal = tOprGlobal(keepIdx);
        vOprLocal = vOprLocal(keepIdx);
    end

    [oprArrivals, tailData(cfg.opr_id)] = extract_opr_edges_local( ...
        tOprGlobal, vOprLocal, cfg.opr_threshold, cfg.gap_points, tailData(cfg.opr_id));
    oprTimes = [oprTimes; oprArrivals]; %#ok<AGROW>

    for sid = cfg.sensor_ids
        if ~channel_file_exists_local(caseDir, sid, fileId)
            continue;
        end
        [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, cfg.pinlv);
        tGlobal = tLocal + offset;
        if iFile == 1 && cfg.initial_trim_points > 0
            keepIdx = (cfg.initial_trim_points + 1):numel(tGlobal);
            tGlobal = tGlobal(keepIdx);
            vLocal = vLocal(keepIdx);
        end

        [pulseInfo, nextTail] = extract_sensor_pulses_local( ...
            tGlobal, vLocal, cfg.sensor_thresholds(sid), cfg.gap_points, ...
            cfg.sensor_peak_fit_degree, cfg.sensor_peak_fit_refine_points, tailData(sid));
        tailData(sid) = nextTail;
        if isempty(pulseInfo.arrival_times)
            continue;
        end

        channels(sid).arrival_times = [channels(sid).arrival_times; pulseInfo.arrival_times]; %#ok<AGROW>
        channels(sid).peak_features = [channels(sid).peak_features; pulseInfo.peak_features]; %#ok<AGROW>
        channels(sid).start_times = [channels(sid).start_times; pulseInfo.start_times]; %#ok<AGROW>
        channels(sid).end_times = [channels(sid).end_times; pulseInfo.end_times]; %#ok<AGROW>
    end
end

[omegaTimeS, omegaRadS, omegaRpm] = compute_speed_local(oprTimes, cfg.opr_pulses_per_rev);
caseData = struct();
caseData.case_name = cfg.low_speed_case;
caseData.case_dir = caseDir;
caseData.file_ids = fileIds;
caseData.sensor_ids = cfg.sensor_ids(:).';
caseData.channels = channels;
caseData.opr_times = oprTimes(:);
caseData.omega_time_s = omegaTimeS(:);
caseData.omega_rad_s = omegaRadS(:);
caseData.omega_rpm = omegaRpm(:);
caseData.file_ranges = fileRanges;
end

function lib = build_sensor_revolution_library_local(caseData, sid, cfg)
arrivals = caseData.channels(sid).arrival_times(:);
peaks = caseData.channels(sid).peak_features(:);
revIds = [];
localPeakRows = zeros(0, cfg.blades_num);
pulseIdxRows = zeros(0, cfg.blades_num);
pulseTimeRows = zeros(0, cfg.blades_num);
numRevs = numel(caseData.opr_times) - cfg.opr_pulses_per_rev;
for revId = 1:numRevs
    t0 = caseData.opr_times(revId);
    t1 = caseData.opr_times(revId + cfg.opr_pulses_per_rev);
    idx = find(arrivals >= t0 & arrivals < t1);
    if numel(idx) ~= cfg.blades_num
        continue;
    end
    [~, order] = sort(arrivals(idx), 'ascend');
    idx = idx(order);
    revIds(end + 1, 1) = revId; %#ok<AGROW>
    localPeakRows(end + 1, :) = peaks(idx).'; %#ok<AGROW>
    pulseIdxRows(end + 1, :) = idx(:).'; %#ok<AGROW>
    pulseTimeRows(end + 1, :) = arrivals(idx).'; %#ok<AGROW>
end
lib = struct( ...
    'sensor_id', sid, ...
    'revolution_ids', revIds, ...
    'local_peak_vectors', localPeakRows, ...
    'pulse_indices', pulseIdxRows, ...
    'pulse_times', pulseTimeRows);
end

function [rowIdx, revId, refFp] = choose_reference_revolution_local(refLibrary, cfg)
candidateCount = min(cfg.reference_search_revs, size(refLibrary.local_peak_vectors, 1));
if candidateCount < 1
    error('Reference sensor does not contain enough candidate revolutions.');
end
candidateMat = refLibrary.local_peak_vectors(1:candidateCount, :);
[~, rowIdx] = max(max(candidateMat, [], 2));
revId = refLibrary.revolution_ids(rowIdx);
refFp = candidateMat(rowIdx, :);
end

function match = match_revolution_library_local(localPeakMat, refFp, bladeCount)
nRevs = size(localPeakMat, 1);
shiftScoresByRev = nan(nRevs, bladeCount);
bestShiftByRev = nan(nRevs, 1);
bestScoreByRev = nan(nRevs, 1);
scoreMarginByRev = nan(nRevs, 1);
refNorm = normalize_vector_max_local(refFp(:));
for i = 1:nRevs
    peakVec = normalize_vector_max_local(localPeakMat(i, :).');
    if ~all(isfinite(peakVec))
        continue;
    end
    for shift = 0:(bladeCount - 1)
        shifted = circshift(peakVec, -shift);
        shiftScoresByRev(i, shift + 1) = corr_scalar_local(shifted, refNorm);
    end
    [bestScoreByRev(i), pos] = max(shiftScoresByRev(i, :));
    if ~isempty(pos) && isfinite(bestScoreByRev(i))
        bestShiftByRev(i) = pos - 1;
        scoreMarginByRev(i) = score_margin_local(shiftScoresByRev(i, :));
    end
end
dominantShift = choose_dominant_shift_local(bestShiftByRev, bestScoreByRev, scoreMarginByRev);
localToPhysical = mod((1:bladeCount) - 1 - dominantShift, bladeCount) + 1;
physicalToLocal = invert_mapping_local(localToPhysical);
match = struct( ...
    'shift_scores_by_revolution', shiftScoresByRev, ...
    'best_shift_by_revolution', bestShiftByRev, ...
    'best_score_by_revolution', bestScoreByRev, ...
    'score_margin_by_revolution', scoreMarginByRev, ...
    'dominant_shift', dominantShift, ...
    'local_to_physical', localToPhysical, ...
    'physical_to_local', physicalToLocal);
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

function physicalMat = reorder_peak_matrix_local(localPeakMat, physicalToLocal)
physicalMat = nan(size(localPeakMat));
for bladeId = 1:numel(physicalToLocal)
    localSlot = physicalToLocal(bladeId);
    physicalMat(:, bladeId) = localPeakMat(:, localSlot);
end
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
C = corrcoef(a, b);
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

function visualize_sensor_config_local(caseData, summaryTable, Sensor_Config, cfg, saveFigures, figureDir)
fig = figure('Name', 'Step00 low-speed Sensor_Config', 'Color', 'w', ...
    'Position', [120, 120, 1200, 760], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(caseData.omega_time_s, caseData.omega_rpm, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('RPM');
title('Low-speed OPR-derived speed');
grid on; box on;

nexttile;
fpMat = nan(numel(cfg.sensor_ids), cfg.blades_num);
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    fpMat(i, :) = Sensor_Config.Fingerprints(sid);
end
plot(1:cfg.blades_num, fpMat.', '-o', 'LineWidth', 1.1);
xlabel('Physical blade ID');
ylabel('Fingerprint peak');
title('Low-speed six-peak fingerprints');
legend(compose('CH%d', cfg.sensor_ids), 'Location', 'bestoutside');
grid on; box on;

nexttile;
yyaxis left;
bar(categorical(compose('CH%d', summaryTable.sensor_id)), summaryTable.dominant_shift);
ylabel('Dominant cyclic shift');
yyaxis right;
plot(categorical(compose('CH%d', summaryTable.sensor_id)), summaryTable.valid_revolution_count, 'o-k', 'LineWidth', 1.0);
ylabel('Valid revolutions');
title('OPR-defined revolution numbering summary');
grid on; box on;

nexttile;
axis off;
summaryText = evalc('disp(summaryTable)');
text(0, 1, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Consolas', 'FontSize', 9, 'Interpreter', 'none');
title('Summary table');

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step00_LowSpeedSensorConfig_20241106.png'), 'Resolution', 300);
end
end

function print_sensor_angle_summary_local(Sensor_Config, cfg)
sensorIds = cfg.sensor_ids(:).';
angleMat = Sensor_Config.Standard_Relative_Angles(sensorIds, :);
bladeIds = (1:cfg.blades_num).';

fprintf('\nSensor standard relative angles by blade (deg, referenced to previous OPR pulse)\n');
angleTable = array2table(angleMat.', 'VariableNames', compose('CH%d', sensorIds));
angleTable = addvars(angleTable, bladeIds, 'Before', 1, 'NewVariableNames', 'BladeID');
disp(angleTable);
end

function oprReference = build_opr_center_reference_for_case_local(caseDir, cfg)
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
centerTimes = [];
startTimes = [];
lastEndTime = [];
tail = struct('t', [], 'v', []);

for iFile = 1:numel(fileIds)
    fileId = fileIds(iFile);
    [tLocal, vLocal] = load_raw_channel_local(caseDir, cfg.opr_id, fileId, cfg.pinlv);
    if isempty(lastEndTime)
        offset = 0;
    else
        offset = lastEndTime + 1 / cfg.pinlv - tLocal(1);
    end
    tGlobal = tLocal + offset;
    lastEndTime = tGlobal(end);
    if iFile == 1 && cfg.initial_trim_points > 0
        keepIdx = (cfg.initial_trim_points + 1):numel(tGlobal);
        tGlobal = tGlobal(keepIdx);
        vLocal = vLocal(keepIdx);
    end
    [segCenter, segStart, tail] = extract_opr_centers_local( ...
        tGlobal, vLocal, cfg.opr_threshold, cfg.gap_points, tail);
    centerTimes = [centerTimes; segCenter(:)]; %#ok<AGROW>
    startTimes = [startTimes; segStart(:)]; %#ok<AGROW>
end

oprReference = build_opr_reference_from_center_start_local( ...
    centerTimes, startTimes, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
oprReference.mode = 'multi_threshold_center';
oprReference.standard_angle_reference = 'opr_pulse_center';
oprReference.source_case = cfg.low_speed_case;
end

function fileIds = list_case_file_ids_local(caseDir, channelId)
d = dir(fullfile(caseDir, sprintf('4-%d-*.mat', channelId)));
fileIds = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channelId) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        fileIds(i) = str2double(tok{1});
    end
end
fileIds = sort(unique(fileIds(~isnan(fileIds))));
end

function [tSec, v] = load_raw_channel_local(caseDir, channelId, fileId, pinlv)
filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
if ~isfile(filepath)
    error('Missing channel file: %s', filepath);
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function tf = channel_file_exists_local(caseDir, channelId, fileId)
tf = isfile(fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId)));
end

function [pulseInfo, nextTail] = extract_sensor_pulses_local(t, v, threshold, gapPoints, fitDegree, refinePoints, tail)
if nargin < 7 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

smoothV = smooth_signal_local(v);
idxAbove = find(smoothV > threshold);
pulseInfo = struct('arrival_times', [], 'peak_features', [], 'start_times', [], 'end_times', []);
nextTail = struct('t', [], 'v', []);
if isempty(idxAbove)
    return;
end

[segStart, segEnd] = split_segments_local(idxAbove, gapPoints);
if numel(v) - idxAbove(end) < gapPoints
    carryStart = max(1, segStart(end) - gapPoints);
    nextTail.t = t(carryStart:end);
    nextTail.v = v(carryStart:end);
    segStart(end) = [];
    segEnd(end) = [];
end
if isempty(segStart)
    return;
end

nSeg = numel(segStart);
arrivalTimes = nan(nSeg, 1);
peakFeatures = nan(nSeg, 1);
startTimes = nan(nSeg, 1);
endTimes = nan(nSeg, 1);
for i = 1:nSeg
    idxRange = segStart(i):segEnd(i);
    tSeg = t(idxRange);
    vSeg = v(idxRange);
    smoothSeg = smoothV(idxRange);
    [arrivalTimes(i), peakFeatures(i)] = fit_arrival_time_local( ...
        tSeg, vSeg, smoothSeg, threshold, fitDegree, refinePoints);
    startTimes(i) = tSeg(1);
    endTimes(i) = tSeg(end);
end

pulseInfo.arrival_times = arrivalTimes;
pulseInfo.peak_features = peakFeatures;
pulseInfo.start_times = startTimes;
pulseInfo.end_times = endTimes;
end

function [oprArrivals, nextTail] = extract_opr_edges_local(t, v, threshold, gapPoints, tail)
if nargin < 5 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

idxAbove = find(v > threshold);
oprArrivals = [];
nextTail = struct('t', [], 'v', []);
if isempty(idxAbove)
    return;
end

[segStart, ~] = split_segments_local(idxAbove, gapPoints);
if numel(v) - idxAbove(end) < gapPoints
    carryStart = max(1, segStart(end) - gapPoints);
    nextTail.t = t(carryStart:end);
    nextTail.v = v(carryStart:end);
    segStart(end) = [];
end
if isempty(segStart)
    return;
end

oprArrivals = nan(numel(segStart), 1);
for i = 1:numel(segStart)
    idx = segStart(i);
    oprArrivals(i) = interpolate_threshold_time_local(t, v, idx, threshold);
end
end

function [segStart, segEnd] = split_segments_local(idxAbove, gapPoints)
jumps = find(diff(idxAbove) > gapPoints);
segStart = [idxAbove(1); idxAbove(jumps + 1)];
segEnd = [idxAbove(jumps); idxAbove(end)];
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

function [arrivalTime, peakFeature] = fit_arrival_time_local(tSeg, vSeg, smoothSeg, threshold, fitDegree, refinePoints)
validMask = smoothSeg > threshold;
if sum(validMask) < 4
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
    return;
end

tFit = tSeg(validMask);
vFit = vSeg(validMask);
mu = mean(tFit);
uniqueCount = numel(unique(tFit));
order = min(fitDegree, uniqueCount - 1);
if order < 1
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
    return;
end

try
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    p = polyfit(tFit - mu, vFit, order);
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    dt = linspace(min(tFit) - mu, max(tFit) - mu, refinePoints);
    vf = polyval(p, dt);
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;
    if sum(weights) > 0
        arrivalTime = sum(tf .* weights) / sum(weights);
        peakFeature = max(vf);
    else
        [peakFeature, idxMax] = max(vf);
        arrivalTime = tf(idxMax);
    end
catch
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
end
end

function [omegaTimeS, omegaRadS, omegaRpm] = compute_speed_local(oprTimes, oprPulsesPerRev)
if numel(oprTimes) <= oprPulsesPerRev
    omegaTimeS = [];
    omegaRadS = [];
    omegaRpm = [];
    return;
end
revPeriod = oprTimes((oprPulsesPerRev + 1):end) - oprTimes(1:(end - oprPulsesPerRev));
omegaTimeS = oprTimes(1:(end - oprPulsesPerRev));
omegaRadS = 2 * pi ./ revPeriod;
omegaRpm = 60 ./ revPeriod;
end

function [centerTimes, startTimes, nextTail] = extract_opr_centers_local(t, v, threshold, gapPoints, tail)
if nargin < 5 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

idxAbove = find(v > threshold);
centerTimes = [];
startTimes = [];
nextTail = struct('t', [], 'v', []);
if isempty(idxAbove)
    return;
end

[segStart, segEnd] = split_segments_local(idxAbove, gapPoints);
if numel(v) - idxAbove(end) < gapPoints
    carryStart = max(1, segStart(end) - gapPoints);
    nextTail.t = t(carryStart:end);
    nextTail.v = v(carryStart:end);
    segStart(end) = [];
    segEnd(end) = [];
end

for iSeg = 1:numel(segStart)
    a = segStart(iSeg);
    b = segEnd(iSeg);
    startTimes(end + 1, 1) = interpolate_threshold_time_local(t, v, a, threshold); %#ok<AGROW>
    centerTimes(end + 1, 1) = compute_multithreshold_center_local(t, v, a, b); %#ok<AGROW>
end
end

function t0 = interpolate_threshold_time_local(t, v, idx, threshold)
if idx <= 1
    t0 = t(idx);
    return;
end
t1 = t(idx - 1);
t2 = t(idx);
v1 = v(idx - 1);
v2 = v(idx);
if abs(v2 - v1) < eps
    t0 = t2;
else
    t0 = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
end
end

function tCenter = compute_multithreshold_center_local(t, v, a, b)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(numel(v), b + pad);
tSeg = t(a0:b0);
vSeg = v(a0:b0);
try
    vSmooth = smooth(vSeg, 16);
catch
    vSmooth = smoothdata(vSeg, 'movmean', 16);
end
[peakVal, iPeak] = max(vSmooth);
baseVal = median(vSmooth(vSmooth <= prctile(vSmooth, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vSmooth);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] * (peakVal - baseVal);
centers = nan(numel(levels), 1);
for i = 1:numel(levels)
    level = levels(i);
    iRise = find(vSmooth(1:iPeak) >= level, 1, 'first');
    iFallRel = find(vSmooth(iPeak:end) <= level, 1, 'first');
    if isempty(iRise) || isempty(iFallRel)
        continue;
    end
    iFall = iPeak + iFallRel - 1;
    tRise = interpolate_level_crossing_local(tSeg, vSmooth, iRise, level);
    tFall = interpolate_level_crossing_local(tSeg, vSmooth, iFall, level);
    centers(i) = 0.5 * (tRise + tFall);
end
tCenter = median(centers, 'omitnan');
if ~isfinite(tCenter)
    tCenter = tSeg(iPeak);
end
end

function tx = interpolate_level_crossing_local(t, v, idx, level)
i1 = max(1, idx - 1);
i2 = idx;
if i1 == i2 || abs(v(i2) - v(i1)) < eps
    tx = t(idx);
else
    tx = t(i1) + (level - v(i1)) * (t(i2) - t(i1)) / (v(i2) - v(i1));
end
end

function oprReference = build_opr_reference_from_center_start_local(centerTimes, startTimes, oprPulsesPerRev, rTipMM)
oprReference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0, 'median_center_minus_start_s', 0);
n = min(numel(centerTimes) - oprPulsesPerRev, numel(startTimes));
if n < 1
    return;
end
dtCenter = centerTimes(1:n) - startTimes(1:n);
dtRev = centerTimes((1:n) + oprPulsesPerRev) - centerTimes(1:n);
valid = isfinite(dtCenter) & isfinite(dtRev) & dtRev > eps;
if ~any(valid)
    return;
end
shiftDeg = 360 * dtCenter(valid) ./ dtRev(valid);
oprReference.phase_shift_deg = median(shiftDeg, 'omitnan');
oprReference.phase_shift_mm = oprReference.phase_shift_deg * (pi / 180) * rTipMM;
oprReference.median_center_minus_start_s = median(dtCenter(valid), 'omitnan');
end
