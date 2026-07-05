%% Step01_Build_LowSpeed_Reference_20241106
% Build the method-neutral low-speed reference layer for the 20241106 BTT
% data-foundation route.

clear; close all; clc;

%% Runtime switches
datasetRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
showPlots = true;

dataCfg = BTTDataConfig_20241106();
dataCfg = fill_missing_step01_defaults_local(dataCfg);

datasetRoot = dataCfg.dataset_root;
lowSpeedCase = dataCfg.low_speed_case;
extractSensorIds = dataCfg.low_speed_extract_sensor_ids;
sensorIds = extractSensorIds;
analysisSensorIds = dataCfg.sensor_ids;
capacitanceIds = dataCfg.capacitance_ids;
eddyCurrentIds = dataCfg.eddy_current_ids;
unusedIds = dataCfg.unused_ids;
oprId = dataCfg.opr_id;

bladesNum = dataCfg.blades_num;
oprPulsesPerRev = dataCfg.opr_pulses_per_rev;
pinlv = dataCfg.sample_rate_hz;
rTipMM = dataCfg.r_tip_mm;
initialTrimPoints = dataCfg.initial_trim_points;
gapPoints = dataCfg.gap_points;
oprThreshold = dataCfg.opr_threshold;

referenceSensorId = dataCfg.reference_sensor_id;
referenceSearchRevs = dataCfg.reference_search_revs;
maxAlignmentSearchStarts = dataCfg.max_alignment_search_starts;
maxLapsProcess = dataCfg.max_laps_process;

sensorThresholds = dataCfg.sensor_thresholds;
sensorThresholdIds = dataCfg.sensor_threshold_ids;
sensorThresholdVals = dataCfg.sensor_threshold_values;
sensorPeakFitDegree = dataCfg.sensor_peak_fit_degree;
sensorPeakFitRefinePoints = dataCfg.sensor_peak_fit_refine_points;
saveFigures = dataCfg.save_figures;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
caseDir = fullfile(datasetRoot, lowSpeedCase);
referenceOutputDir = dataCfg.step01_output_dir;
sensorConfigFile = fullfile(referenceOutputDir, 'Sensor_Config_20241106.mat');
sensorConfigSummaryFile = fullfile(referenceOutputDir, 'Sensor_Config_Summary_20241106.csv');
lowSpeedFeatureFile = fullfile(referenceOutputDir, 'LowSpeed_Features_20241106.mat');
lowSpeedReferenceFile = fullfile(referenceOutputDir, 'LowSpeedReference_20241106.mat');
lowSpeedFingerprintFile = fullfile(referenceOutputDir, 'LowSpeedReferenceFingerprint_20241106.mat');
lowSpeedNumberingFile = fullfile(referenceOutputDir, 'LowSpeedNumbering_20241106.mat');
lowSpeedNumberingSummaryFile = fullfile(referenceOutputDir, 'LowSpeedNumbering_Summary_20241106.csv');
channelStatsFile = fullfile(referenceOutputDir, 'LowSpeed_Channel_Stats_20241106.csv');
angleQualityFile = fullfile(referenceOutputDir, 'Standard_Relative_Angles_Quality_20241106.csv');
figureDir = dataCfg.step01_figure_dir;

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
cfg.dataset = dataCfg.dataset;
cfg.route_dir = routeDir;
cfg.output_root = dataCfg.output_root;
cfg.step01_output_dir = dataCfg.step01_output_dir;
cfg.step01_figure_dir = dataCfg.step01_figure_dir;
cfg.dataset_root = datasetRoot;
cfg.low_speed_case = lowSpeedCase;
cfg.reference_output_dir = referenceOutputDir;
cfg.extract_sensor_ids = extractSensorIds;
cfg.analysis_sensor_ids = analysisSensorIds;
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

fprintf('\n=== Step01: build low-speed reference layer ===\n');
fprintf('Low-speed data folder: %s\n', caseDir);
fprintf('Extract sensors: %s\n', mat2str(sensorIds));
fprintf('Analysis sensors: %s\n', mat2str(analysisSensorIds));
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
Sensor_Config.Dataset = dataCfg.dataset;
Sensor_Config.ReferenceCase = lowSpeedCase;
Sensor_Config.Sensor_IDs = sensorIds;
Sensor_Config.Extract_Sensor_IDs = sensorIds;
Sensor_Config.Analysis_Sensor_IDs = analysisSensorIds;
Sensor_Config.Candidate_Sensor_IDs = dataCfg.candidate_sensor_ids;
Sensor_Config.Optional_Sensor_IDs = dataCfg.optional_sensor_ids;
Sensor_Config.Capacitance_IDs = capacitanceIds;
Sensor_Config.Eddy_Current_IDs = eddyCurrentIds;
Sensor_Config.Unused_IDs = unusedIds;
Sensor_Config.OPR_ID = oprId;
Sensor_Config.Blades_Num = bladesNum;
Sensor_Config.Sample_Rate_Hz = pinlv;
Sensor_Config.Pinlv = pinlv;
Sensor_Config.R_Tip_mm = rTipMM;
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
Sensor_Config.Standard_Relative_Angles_SelectedStatistic = dataCfg.standard_angle_value;
Sensor_Config.OPR_Timing_Method = dataCfg.opr_timing_method;
Sensor_Config.OPR_Center_Correction_Method = dataCfg.opr_center_correction_method;
Sensor_Config.OPR_Start_Times = caseData.opr_start_times;
Sensor_Config.OPR_Center_Times = oprReference.center_times;
Sensor_Config.OPR_Center_Phase_Shift_Deg = oprReference.phase_shift_deg;
Sensor_Config.OPR_Center_Phase_Shift_mm = oprReference.phase_shift_mm;
Sensor_Config.OPR_Center_Phase_Shift_Deg_IQR = oprReference.phase_shift_deg_iqr;
Sensor_Config.OPR_Center_Minus_Start_s_Median = oprReference.median_center_minus_start_s;
Sensor_Config.OPR_Center_Minus_Start_s_IQR = oprReference.center_minus_start_s_iqr;
Sensor_Config.OPR_Center_Valid_Count = oprReference.valid_count;
Sensor_Config.OPR_Start_Time_Count = oprReference.start_time_count;
Sensor_Config.OPR_Center_Time_Count = oprReference.center_time_count;
Sensor_Config.Probe_Arrival_Method = 'polynomial_peak_fit';
Sensor_Config.Blade_ID_Rule = dataCfg.low_speed_blade1_rule;
Sensor_Config.OPRReference = oprReference;
Sensor_Config.Peak_Fit_Degree = cfg.sensor_peak_fit_degree;
Sensor_Config.Peak_Fit_Refine_Points = cfg.sensor_peak_fit_refine_points;
Sensor_Config.CreatedBy = mfilename;
Sensor_Config.CreatedOn = datestr(now, 31);

summaryTable = struct2table(summaryRows);
channelStats = build_channel_stats_table_local(caseData, sensorIds);
angleQualityTable = build_angle_quality_table_local(Sensor_Config, sensorIds, bladesNum);
[LowSpeedReference, fingerprintRows] = build_low_speed_reference_local( ...
    Sensor_Config, analysisSensorIds, bladesNum, sensorConfigFile);
[LowSpeedNumbering, LowSpeedSummary] = build_low_speed_numbering_tables_local( ...
    LowSpeedReference, analysisSensorIds, bladesNum);

case_data = caseData; %#ok<NASGU>
save(sensorConfigFile, 'Sensor_Config');
save(lowSpeedFeatureFile, 'caseData', 'case_data', 'cfg', '-v7.3');
save(lowSpeedReferenceFile, 'LowSpeedReference', 'fingerprintRows');
save(lowSpeedFingerprintFile, 'LowSpeedReference', 'fingerprintRows');
save(lowSpeedNumberingFile, 'LowSpeedNumbering', 'LowSpeedSummary');
writetable(summaryTable, sensorConfigSummaryFile);
writetable(channelStats, channelStatsFile);
writetable(angleQualityTable, angleQualityFile);
writetable(LowSpeedNumbering, strrep(lowSpeedNumberingFile, '.mat', '.csv'));
writetable(LowSpeedSummary, lowSpeedNumberingSummaryFile);

fprintf('Saved Step01 low-speed reference products to:\n  %s\n', referenceOutputDir);
disp(summaryTable);
disp(fingerprintRows);
disp(LowSpeedSummary);
print_sensor_angle_summary_local(Sensor_Config, cfg);

if showPlots
    visualize_sensor_config_local(caseData, summaryTable, Sensor_Config, cfg, saveFigures, figureDir);
end

function cfg = fill_missing_step01_defaults_local(cfg)
if ~isfield(cfg, 'low_speed_extract_sensor_ids') || isempty(cfg.low_speed_extract_sensor_ids)
    if isfield(cfg, 'sensor_ids') && ~isempty(cfg.sensor_ids)
        cfg.low_speed_extract_sensor_ids = cfg.sensor_ids;
    else
        cfg.low_speed_extract_sensor_ids = cfg.candidate_sensor_ids;
    end
end
if ~isfield(cfg, 'pinlv') || isempty(cfg.pinlv)
    cfg.pinlv = cfg.sample_rate_hz;
end
if ~isfield(cfg, 'sample_rate_hz') || isempty(cfg.sample_rate_hz)
    cfg.sample_rate_hz = cfg.pinlv;
end
if ~isfield(cfg, 'candidate_sensor_ids') || isempty(cfg.candidate_sensor_ids)
    cfg.candidate_sensor_ids = cfg.low_speed_extract_sensor_ids;
end
if ~isfield(cfg, 'sensor_thresholds') || isempty(cfg.sensor_thresholds)
    cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
    ids = unique([cfg.candidate_sensor_ids(:); cfg.sensor_ids(:); cfg.opr_id]);
    for i = 1:numel(ids)
        cfg.sensor_thresholds(ids(i)) = cfg.sensor_threshold_default;
    end
end
end

function T = build_channel_stats_table_local(caseData, sensorIds)
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'pulse_count', NaN, ...
    'first_arrival_s', NaN, ...
    'last_arrival_s', NaN, ...
    'peak_mean', NaN, ...
    'peak_std', NaN, ...
    'peak_min', NaN, ...
    'peak_max', NaN), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    arrivals = caseData.channels(sid).arrival_times(:);
    peaks = caseData.channels(sid).peak_features(:);
    rows(i).sensor_id = sid;
    rows(i).pulse_count = numel(arrivals);
    if ~isempty(arrivals)
        rows(i).first_arrival_s = arrivals(1);
        rows(i).last_arrival_s = arrivals(end);
    end
    rows(i).peak_mean = mean(peaks, 'omitnan');
    rows(i).peak_std = std(peaks, 'omitnan');
    rows(i).peak_min = min(peaks, [], 'omitnan');
    rows(i).peak_max = max(peaks, [], 'omitnan');
end
T = struct2table(rows);
end

function T = build_angle_quality_table_local(Sensor_Config, sensorIds, bladeCount)
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'selected_angle_deg', NaN, ...
    'start_edge_angle_deg', NaN, ...
    'angle_reference', "", ...
    'is_finite', false), numel(sensorIds) * bladeCount, 1);
k = 0;
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    for bladeId = 1:bladeCount
        k = k + 1;
        rows(k).sensor_id = sid;
        rows(k).blade_id = bladeId;
        rows(k).selected_angle_deg = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, bladeId);
        rows(k).start_edge_angle_deg = Sensor_Config.Standard_Relative_Angles_StartEdge(sid, bladeId);
        rows(k).angle_reference = string(Sensor_Config.Standard_Relative_Angles_Reference);
        rows(k).is_finite = isfinite(rows(k).selected_angle_deg);
    end
end
T = struct2table(rows);
end

function [LowSpeedReference, fingerprintRows] = build_low_speed_reference_local( ...
        Sensor_Config, analysisSensors, bladeCount, sensorConfigFile)
validate_sensor_config_for_reference_local(Sensor_Config, analysisSensors, bladeCount, sensorConfigFile);
fingerprintRows = build_fingerprint_table_local(Sensor_Config, analysisSensors, bladeCount);

LowSpeedReference = struct();
LowSpeedReference.source_file = sensorConfigFile;
LowSpeedReference.dataset = Sensor_Config.Dataset;
LowSpeedReference.analysis_sensors = analysisSensors;
LowSpeedReference.blade_count = bladeCount;
LowSpeedReference.reference_sensor_id = Sensor_Config.ReferenceSensorID;
LowSpeedReference.reference_revolution_id = Sensor_Config.ReferenceRevolutionID;
LowSpeedReference.reference_fingerprint = Sensor_Config.ReferenceFingerprint;
LowSpeedReference.fingerprint_table = fingerprintRows;
LowSpeedReference.fingerprints = Sensor_Config.Fingerprints;
LowSpeedReference.target_indices = Sensor_Config.Target_Indices;
LowSpeedReference.revolution_ids = Sensor_Config.Revolution_IDs;
LowSpeedReference.revolution_local_peaks = Sensor_Config.Revolution_Local_Peaks;
LowSpeedReference.revolution_pulse_times = Sensor_Config.Revolution_Pulse_Times;
LowSpeedReference.revolution_physical_peaks = Sensor_Config.Revolution_Physical_Peaks;
LowSpeedReference.revolution_best_shifts = Sensor_Config.Revolution_Best_Shifts;
LowSpeedReference.revolution_best_scores = Sensor_Config.Revolution_Best_Scores;
LowSpeedReference.revolution_score_margins = Sensor_Config.Revolution_Score_Margins;
LowSpeedReference.sensor_dominant_shifts = Sensor_Config.Sensor_Dominant_Shifts;
LowSpeedReference.sensor_local_to_physical = Sensor_Config.Sensor_Local_To_Physical;
LowSpeedReference.sensor_physical_to_local = Sensor_Config.Sensor_Physical_To_Local;
LowSpeedReference.standard_angles_opr_center = Sensor_Config.Standard_Relative_Angles_OPRCenter;
LowSpeedReference.standard_angles_start_edge = Sensor_Config.Standard_Relative_Angles_StartEdge;
LowSpeedReference.opr_timing_method = Sensor_Config.OPR_Timing_Method;
LowSpeedReference.opr_center_correction_method = Sensor_Config.OPR_Center_Correction_Method;
LowSpeedReference.opr_start_times = Sensor_Config.OPR_Start_Times;
LowSpeedReference.opr_center_times = Sensor_Config.OPR_Center_Times;
LowSpeedReference.opr_center_phase_shift_deg = Sensor_Config.OPR_Center_Phase_Shift_Deg;
LowSpeedReference.opr_center_phase_shift_mm = Sensor_Config.OPR_Center_Phase_Shift_mm;
LowSpeedReference.opr_center_phase_shift_deg_iqr = Sensor_Config.OPR_Center_Phase_Shift_Deg_IQR;
LowSpeedReference.opr_center_minus_start_s_median = Sensor_Config.OPR_Center_Minus_Start_s_Median;
LowSpeedReference.opr_center_minus_start_s_iqr = Sensor_Config.OPR_Center_Minus_Start_s_IQR;
LowSpeedReference.opr_reference = Sensor_Config.OPRReference;
LowSpeedReference.probe_arrival_method = Sensor_Config.Probe_Arrival_Method;
LowSpeedReference.blade_id_rule = Sensor_Config.Blade_ID_Rule;
LowSpeedReference.created_by = mfilename;
LowSpeedReference.created_on = datestr(now, 31);
end

function validate_sensor_config_for_reference_local(Sensor_Config, sensorIds, bladeCount, sensorConfigFile)
requiredFields = {'Fingerprints', 'Target_Indices', 'ReferenceSensorID', 'ReferenceRevolutionID', ...
    'Revolution_IDs', 'Revolution_Pulse_Times', 'Revolution_Physical_Peaks', ...
    'Sensor_Dominant_Shifts', 'Sensor_Physical_To_Local'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(Sensor_Config, name)
        error('Sensor_Config is missing field %s:\n  %s', name, sensorConfigFile);
    end
end
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Sensor_Config.Fingerprints does not contain CH%d.', sid);
    end
    fp = Sensor_Config.Fingerprints(sid);
    if numel(fp) < bladeCount
        error('Fingerprint for CH%d has only %d values, but bladeCount=%d.', ...
            sid, numel(fp), bladeCount);
    end
end
end

function T = build_fingerprint_table_local(Sensor_Config, sensorIds, bladeCount)
rows = repmat(struct('SensorID', NaN, 'ReferenceRevolutionID', NaN, 'DominantShift', NaN, ...
    'ValidRevolutionCount', NaN, 'Fingerprint', ''), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = Sensor_Config.Fingerprints(sid);
    rows(i).SensorID = sid;
    rows(i).ReferenceRevolutionID = Sensor_Config.ReferenceRevolutionID;
    rows(i).DominantShift = Sensor_Config.Sensor_Dominant_Shifts(sid);
    rows(i).ValidRevolutionCount = numel(Sensor_Config.Revolution_IDs(sid));
    rows(i).Fingerprint = mat2str(fp(1:bladeCount), 6);
end
T = struct2table(rows);
end

function [numberingTable, summaryTable] = build_low_speed_numbering_tables_local( ...
        LowSpeedReference, sensorIds, bladeCount)
numberingRows = struct([]);
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'ReferenceRevolutionID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ValidRevolutionCount', NaN, ...
    'MeanBestScore', NaN, ...
    'MinBestScore', NaN, ...
    'MeanScoreMargin', NaN, ...
    'FingerprintB1', NaN, ...
    'FingerprintB2', NaN, ...
    'FingerprintB3', NaN, ...
    'FingerprintB4', NaN, ...
    'FingerprintB5', NaN, ...
    'FingerprintB6', NaN), numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    revIds = force_row_vector_local(LowSpeedReference.revolution_ids(sid));
    peakMat = LowSpeedReference.revolution_physical_peaks(sid);
    bestShift = force_row_vector_local(LowSpeedReference.revolution_best_shifts(sid));
    bestScore = force_row_vector_local(LowSpeedReference.revolution_best_scores(sid));
    scoreMargin = force_row_vector_local(LowSpeedReference.revolution_score_margins(sid));
    physicalToLocal = force_row_vector_local(LowSpeedReference.sensor_physical_to_local(sid));
    dominantShift = LowSpeedReference.sensor_dominant_shifts(sid);
    fingerprint = force_row_vector_local(LowSpeedReference.fingerprints(sid));

    nRows = size(peakMat, 1);
    for k = 1:nRows
        row = struct();
        row.SensorID = sid;
        row.RevolutionID = revIds(k);
        row.ReferenceRevolutionID = LowSpeedReference.reference_revolution_id;
        row.DominantShift = dominantShift;
        row.BestShift = bestShift(k);
        row.BestScore = bestScore(k);
        row.ScoreMargin = scoreMargin(k);
        row.LocalSlotOfB1 = physicalToLocal(1);
        for b = 1:bladeCount
            row.(sprintf('B%d', b)) = peakMat(k, b);
        end
        numberingRows = [numberingRows; row]; %#ok<AGROW>
    end

    summaryRows(i).SensorID = sid;
    summaryRows(i).ReferenceRevolutionID = LowSpeedReference.reference_revolution_id;
    summaryRows(i).DominantShift = dominantShift;
    summaryRows(i).LocalSlotOfB1 = physicalToLocal(1);
    summaryRows(i).ValidRevolutionCount = nRows;
    summaryRows(i).MeanBestScore = mean(bestScore, 'omitnan');
    summaryRows(i).MinBestScore = min(bestScore, [], 'omitnan');
    summaryRows(i).MeanScoreMargin = mean(scoreMargin, 'omitnan');
    for b = 1:bladeCount
        summaryRows(i).(sprintf('FingerprintB%d', b)) = fingerprint(b);
    end
end

if isempty(numberingRows)
    numberingTable = table();
else
    numberingTable = struct2table(numberingRows);
end
summaryTable = struct2table(summaryRows);
end

function row = force_row_vector_local(x)
row = x(:).';
end

function caseData = extract_btt_features_local(caseDir, cfg)
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR files found in %s.', caseDir);
end

maxChannelId = max([cfg.extract_sensor_ids, cfg.opr_id]);
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

    for sid = cfg.extract_sensor_ids
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
caseData.sensor_ids = cfg.extract_sensor_ids(:).';
caseData.channels = channels;
caseData.opr_times = oprTimes(:);
caseData.opr_start_times = oprTimes(:);
caseData.opr_timing_method = 'threshold_start_edge';
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
fig = figure('Name', 'Step01 low-speed reference', 'Color', 'w', ...
    'Position', [120, 120, 1200, 760], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(caseData.omega_time_s, caseData.omega_rpm, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('RPM');
title('Low-speed OPR-derived speed');
grid on; box on;

nexttile;
plotSensorIds = cfg.analysis_sensor_ids(:).';
fpMat = nan(numel(plotSensorIds), cfg.blades_num);
for i = 1:numel(plotSensorIds)
    sid = plotSensorIds(i);
    fpMat(i, :) = Sensor_Config.Fingerprints(sid);
end
plot(1:cfg.blades_num, fpMat.', '-o', 'LineWidth', 1.1);
xlabel('Physical blade ID');
ylabel('Fingerprint peak');
title('Low-speed six-peak fingerprints');
legend(compose('CH%d', plotSensorIds), 'Location', 'bestoutside');
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
    exportgraphics(fig, fullfile(figureDir, 'Step01_LowSpeedReference_20241106.png'), 'Resolution', 300);
end
end

function print_sensor_angle_summary_local(Sensor_Config, cfg)
sensorIds = cfg.analysis_sensor_ids(:).';
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
oprReference = struct( ...
    'phase_shift_deg', 0, ...
    'phase_shift_deg_iqr', NaN, ...
    'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0, ...
    'center_minus_start_s_iqr', NaN, ...
    'valid_count', 0, ...
    'start_time_count', numel(startTimes), ...
    'center_time_count', numel(centerTimes), ...
    'start_times', startTimes(:), ...
    'center_times', centerTimes(:));
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
oprReference.phase_shift_deg_iqr = iqr(shiftDeg(isfinite(shiftDeg)));
oprReference.phase_shift_mm = oprReference.phase_shift_deg * (pi / 180) * rTipMM;
oprReference.median_center_minus_start_s = median(dtCenter(valid), 'omitnan');
oprReference.center_minus_start_s_iqr = iqr(dtCenter(valid));
oprReference.valid_count = nnz(valid);
end
