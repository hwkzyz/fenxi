%% Step05_BuildLowSpeedTemplateLibrary_20241106
% Build low-speed templates directly from the new OPR-based numbering chain.
% Each blade uses OPR-defined revolutions, extracts the matching six-pulse
% waveform segments, restores them to the x-domain, then applies a
% gradient-based trust-region rule before fitting the final template.
%
% Semantic role in the 20241106 main chain:
%   This is the low-speed waveform calibration step.
%   It corresponds to the "build low-speed template" role that is Step01
%   in the 20251222 OPRCenterStd route.

clear; close all; clc;

P = NewFlow_Config_20241106();

% Step05 local parameter block.
% Edit here directly when you want to adjust the low-speed template domain
% for selected sensors. Each row is [sensor_id, delta_mm_per_side].
% Positive = wider, negative = narrower.
S05 = struct();
S05.sensorThresholdOverrideTable = [2 0.20; 3 0.20];
S05.sensorDomainExpandMMTable = [2 -0.35; 3 -0.35];
S05.etaLimitMM = 0.20;
S05.etaAdaptiveQuantile = 85;
S05.etaAdaptiveSafetyFactor = 1.20;
S05.etaAdaptiveMinMM = 0.02;
S05.etaAdaptiveMaxMM = S05.etaLimitMM;
P.template.sensorThresholdOverrideTable = S05.sensorThresholdOverrideTable;
P.template.sensorDomainExpandMMTable = S05.sensorDomainExpandMMTable;
P.template.etaLimitMM = S05.etaLimitMM;
P.template.etaAdaptiveQuantile = S05.etaAdaptiveQuantile;
P.template.etaAdaptiveSafetyFactor = S05.etaAdaptiveSafetyFactor;
P.template.etaAdaptiveMinMM = S05.etaAdaptiveMinMM;
P.template.etaAdaptiveMaxMM = S05.etaAdaptiveMaxMM;

require_file_local(P.files.sensorConfig, 'Step00 low-speed sensor config');
require_file_local(P.files.lowSpeedFingerprint, 'Step01 low-speed fingerprint');
require_file_local(P.files.lowSpeedNumbering, 'Step02 low-speed numbering');

loadedCfg = load(P.files.sensorConfig, 'Sensor_Config');
Sensor_Config = loadedCfg.Sensor_Config;
loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;
loadedNum = load(P.files.lowSpeedNumbering, 'LowSpeedNumbering', 'LowSpeedSummary');
LowSpeedNumbering = loadedNum.LowSpeedNumbering;
LowSpeedSummary = loadedNum.LowSpeedSummary;

outDir = fileparts(P.files.lowSpeedTemplateLibrary);
templateDir = fullfile(P.outputDir, 'templates', sensor_tag_local(P.sensors.analysis));
singleTemplateDir = fullfile(P.outputDir, 'template_library');
figDir = fullfile(P.view.figureDir, '05A_low_speed_template_library');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
if exist(templateDir, 'dir') ~= 7
    mkdir(templateDir);
end
if exist(singleTemplateDir, 'dir') ~= 7
    mkdir(singleTemplateDir);
end
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

legacyCfg = build_legacy_cfg_local(P, Sensor_Config);
[rawStream, oprTimes, oprReference] = load_low_speed_streams_local( ...
    P.data.lowSpeedDir, legacyCfg, P.sensors.analysis);
F_omega_deg = build_phase_speed_local(oprTimes, P.machine.oprPulsesPerRev);

templateEntries = repmat(struct('blade_id', NaN, 'template_file', "", 'exists', false), ...
    numel(P.template.blades), 1);
summaryRows = repmat(struct( ...
    'BladeID', NaN, ...
    'SensorID', NaN, ...
    'SelectedPulseCount', NaN, ...
    'SelectedRevolutionCount', NaN, ...
    'WidePointCount', NaN, ...
    'TrustPointCount', NaN, ...
    'XDomainLeftMM', NaN, ...
    'XDomainRightMM', NaN, ...
    'XcMM', NaN, ...
    'DetectedXcMM', NaN, ...
    'ThresholdV', NaN, ...
    'BaselineV', NaN, ...
    'EtaMedianMM', NaN, ...
    'EtaIQRMM', NaN, ...
    'EtaLimitMM', NaN, ...
    'GradientMode', "", ...
    'StablePlanMode', "", ...
    'TemplateFile', ""), 0, 1);
auditBlade = repmat(struct( ...
    'blade_id', NaN, ...
    'bundle_file', "", ...
    'sensor', []), numel(P.template.blades), 1);

fprintf('\n=== Step05: low-speed template library from new numbering ===\n');
fprintf('Low-speed raw-data folder: %s\n', P.data.lowSpeedDir);
fprintf('Template xrange mode: %s, center mode: %s\n', ...
    P.template.xrangeMode, P.template.centerMode);

for ib = 1:numel(P.template.blades)
    bladeId = P.template.blades(ib);
    stablePlan = load_stable_window_plan_local(P, bladeId);
    [TemplateBundle, bladeSummaryRows, bladeAudit] = build_single_blade_template_bundle_local( ...
        P, bladeId, Sensor_Config, LowSpeedReference, LowSpeedNumbering, LowSpeedSummary, ...
        rawStream, oprTimes, oprReference, F_omega_deg, stablePlan, templateDir, singleTemplateDir);

    bundleFile = fullfile(templateDir, sprintf( ...
        'TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
        bladeId, sensor_tag_local(P.sensors.analysis), P.template.suffix));
    save(bundleFile, 'TemplateBundle', '-v7.3');

    templateEntries(ib).blade_id = bladeId;
    templateEntries(ib).template_file = string(bundleFile);
    templateEntries(ib).exists = true;
    summaryRows = [summaryRows; bladeSummaryRows]; %#ok<AGROW>
    auditBlade(ib).blade_id = bladeId;
    auditBlade(ib).bundle_file = string(bundleFile);
    auditBlade(ib).sensor = bladeAudit;

    fprintf('B%d template bundle saved: %s\n', bladeId, bundleFile);
end

LowSpeedTemplateLibrary = struct();
LowSpeedTemplateLibrary.sensor_ids = P.sensors.analysis(:).';
LowSpeedTemplateLibrary.sensor_tag = sensor_tag_local(P.sensors.analysis);
LowSpeedTemplateLibrary.template_suffix = P.template.suffix;
LowSpeedTemplateLibrary.entry = templateEntries;
LowSpeedTemplateLibrary.note = ['Built from Step02 low-speed OPR-defined numbering. ' ...
    'Each blade template is reconstructed from selected revolutions, gradient-trimmed in x, ' ...
    'and exported as TemplateBundle-compatible per-sensor files.'];

LowSpeedTemplateSummary = struct2table(summaryRows);
LowSpeedTemplateAudit = struct();
LowSpeedTemplateAudit.source_sensor_config = P.files.sensorConfig;
LowSpeedTemplateAudit.source_low_speed_fingerprint = P.files.lowSpeedFingerprint;
LowSpeedTemplateAudit.source_low_speed_numbering = P.files.lowSpeedNumbering;
LowSpeedTemplateAudit.blade = auditBlade;
LowSpeedTemplateAudit.template_settings = P.template;
LowSpeedTemplateAudit.opr_reference = oprReference;

save(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary', '-v7.3');
save(P.files.lowSpeedTemplateAudit, 'LowSpeedTemplateAudit', '-v7.3');
writetable(LowSpeedTemplateSummary, P.files.lowSpeedTemplateSummary);

fprintf('Saved low-speed template library: %s\n', P.files.lowSpeedTemplateLibrary);
fprintf('Saved low-speed template summary: %s\n', P.files.lowSpeedTemplateSummary);
fprintf('Saved low-speed template audit: %s\n', P.files.lowSpeedTemplateAudit);

if P.view.enable
    visualize_low_speed_template_summary_local(LowSpeedTemplateSummary, P, figDir);
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function cfg = build_legacy_cfg_local(P, Sensor_Config)
cfg = struct();
cfg.low_speed_case = P.data.lowSpeedCase;
cfg.dataset_root = P.data.datasetRoot;
cfg.pinlv = P.machine.sampleRateHz;
cfg.opr_id = P.machine.oprChannel;
cfg.opr_threshold = Sensor_Config.OPR_Threshold;
cfg.opr_pulses_per_rev = P.machine.oprPulsesPerRev;
cfg.blades_num = P.machine.bladeCount;
cfg.gap_points = Sensor_Config.Gap_Points;
cfg.r_tip_mm = P.machine.tipRadiusMM;
cfg.initial_trim_points = 0;
thresholdMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
sensorThresholds = Sensor_Config.Sensor_Thresholds;
for sid = P.sensors.analysis(:).'
    if isa(sensorThresholds, 'containers.Map')
        if ~isKey(sensorThresholds, sid)
            error('Sensor_Config.Sensor_Thresholds does not contain CH%d.', sid);
        end
        thresholdMap(sid) = sensorThresholds(sid);
    elseif isnumeric(sensorThresholds) && numel(sensorThresholds) >= sid
        thresholdMap(sid) = sensorThresholds(sid);
    else
        error('Unsupported Sensor_Thresholds format for CH%d.', sid);
    end
end
cfg.sensor_thresholds = thresholdMap;
end

function [TemplateBundle, bladeSummaryRows, bladeAudit] = build_single_blade_template_bundle_local( ...
        P, bladeId, Sensor_Config, LowSpeedReference, LowSpeedNumbering, LowSpeedSummary, ...
        rawStream, oprTimes, oprReference, F_omega_deg, stablePlan, templateDir, singleTemplateDir)

Template = struct();
Template.Route = 'low_speed_numbering_to_gradient_template';
Template.TargetBlade = bladeId;
Template.SensorIDs = P.sensors.analysis(:).';
Template.SensorTag = sensor_tag_local(P.sensors.analysis);
Template.TemplateSuffix = P.template.suffix;
Template.LowSpeedCase = P.data.lowSpeedCase;
Template.LowSpeedDataDir = P.data.lowSpeedDir;
Template.SensorConfigFile = P.files.sensorConfig;
Template.LowSpeedFingerprintFile = P.files.lowSpeedFingerprint;
Template.LowSpeedNumberingFile = P.files.lowSpeedNumbering;
Template.StableWindowPlanFile = stable_plan_file_text_local(stablePlan);
Template.OPRReference = oprReference;
Template.CreatedBy = mfilename;
Template.Settings = P.template;
Template.Sensor = repmat(struct( ...
    'sensor_id', NaN, ...
    'baseline', NaN, ...
    'threshold', NaN, ...
    'xc', NaN, ...
    'xc_detected', NaN, ...
    'xc_reference_source', '', ...
    'x_grid', [], ...
    'v_grid', [], ...
    'dv_dx', [], ...
    'bin_weight', [], ...
    'x_domain', [], ...
    'selection_mode', '', ...
    'point_count', NaN, ...
    'wide_point_count', NaN, ...
    'lap_count', NaN, ...
    'selected_pulse_indices', [], ...
    'selected_revolution_ids', [], ...
    'fit_method', '', ...
    'spline_smoothing', NaN, ...
    'eta_app_mm', [], ...
    'eta_median_mm', NaN, ...
    'eta_iqr_mm', NaN, ...
    'eta_limit_mm', NaN), numel(P.sensors.analysis), 1);

    bladeSummaryRows = repmat(struct( ...
        'BladeID', NaN, ...
        'SensorID', NaN, ...
        'SelectedPulseCount', NaN, ...
        'SelectedRevolutionCount', NaN, ...
        'WidePointCount', NaN, ...
        'TrustPointCount', NaN, ...
        'XDomainLeftMM', NaN, ...
        'XDomainRightMM', NaN, ...
        'XcMM', NaN, ...
        'DetectedXcMM', NaN, ...
    'ThresholdV', NaN, ...
    'BaselineV', NaN, ...
    'EtaMedianMM', NaN, ...
    'EtaIQRMM', NaN, ...
    'EtaLimitMM', NaN, ...
    'GradientMode', "", ...
    'StablePlanMode', "", ...
    'TemplateFile', ""), numel(P.sensors.analysis), 1);

bladeAudit = repmat(struct( ...
    'sensor_id', NaN, ...
    'selected_revolution_ids', [], ...
    'selected_pulse_indices', [], ...
    'x_wide', [], ...
    'v_wide', [], ...
    'x_selected', [], ...
    'v_selected', [], ...
    'weight_selected', [], ...
    'x_profile', [], ...
    'v_profile', [], ...
    'gradient_abs', [], ...
    'gradient_threshold', NaN, ...
    'effective_mask', [], ...
    'xc_detected_rel', NaN, ...
    'x_domain_rel', [NaN NaN], ...
    'eta_app_mm', [], ...
    'eta_median_mm', NaN, ...
    'eta_iqr_mm', NaN, ...
    'eta_limit_mm', NaN, ...
    'stable_plan_mode', "", ...
    'selection_mode', "", ...
    'baseline', NaN, ...
    'threshold', NaN, ...
    'template_file', ""), numel(P.sensors.analysis), 1);

sensorTemplateFiles = cell(numel(P.sensors.analysis), 1);
referenceCenter = build_reference_center_table_local(P, bladeId, stablePlan, Sensor_Config, LowSpeedReference);

for is = 1:numel(P.sensors.analysis)
    sid = P.sensors.analysis(is);
    threshold = sensor_threshold_local(Sensor_Config, sid, P);
    thetaStd = read_opr_center_standard_angle_local(LowSpeedReference, sid, bladeId, Sensor_Config, oprReference);
    R = rawStream(sid);
    segments = extract_sensor_segments_local(R.T, R.V, threshold, legacy_gap_points_local(Sensor_Config));
    baseline = estimate_background_baseline_local(R.V, segments);

    [selectedRevIds, selectedPulseTimes, stablePlanMode] = choose_low_speed_template_pulses_local( ...
        P, bladeId, sid, Sensor_Config, LowSpeedReference, LowSpeedNumbering, LowSpeedSummary, stablePlan);
    selectedPulseIdx = match_selected_pulse_times_to_segments_local( ...
        segments.arrival_time(:), selectedPulseTimes(:), sid, bladeId);
    if numel(selectedPulseIdx) < P.template.minSelectedPulseCount
        error('B%d CH%d has only %d selected pulses.', bladeId, sid, numel(selectedPulseIdx));
    end

    [xWide, vWide] = restore_low_speed_point_cloud_local( ...
        R, segments, selectedPulseIdx, oprTimes, F_omega_deg, thetaStd, ...
        P.machine.tipRadiusMM, P.template.segmentExpandFactor);
    selectCfg = build_select_cfg_local(P, threshold, baseline, sid);
    selected = select_sg_like_template_points_local(xWide, vWide, selectCfg);

    xcForCoordinate = selected.xc;
    xcReferenceSource = 'detected_from_current_template_points';
    iref = find(referenceCenter.sensor_id == sid, 1);
    if ~isempty(iref) && isfinite(referenceCenter.xc_mm(iref))
        xcForCoordinate = referenceCenter.xc_mm(iref);
        xcReferenceSource = 'reference_center_from_same_blade';
    end

    xRel = selected.x_selected(:) - xcForCoordinate;
    vSel = selected.v_selected(:);
    wSel = selected.weight_selected(:);
    finiteMask = isfinite(xRel) & isfinite(vSel) & isfinite(wSel);
    xRel = xRel(finiteMask);
    vSel = vSel(finiteMask);
    wSel = wSel(finiteMask);

    [xGrid, vGrid, dvDx, binWeight, fitMethod] = fit_template_curve_local(P, xRel, vSel, wSel);
    xDomain = [min(xRel), max(xRel)];
    [etaApp, etaMedian, etaIqr, etaLimit] = estimate_low_speed_eta_scatter_local( ...
        xRel, vSel, xGrid, vGrid, P.template);

    Template.Sensor(is).sensor_id = sid;
    Template.Sensor(is).baseline = baseline;
    Template.Sensor(is).threshold = threshold;
    Template.Sensor(is).xc = xcForCoordinate;
    Template.Sensor(is).xc_detected = selected.xc;
    Template.Sensor(is).xc_reference_source = xcReferenceSource;
    Template.Sensor(is).x_grid = xGrid;
    Template.Sensor(is).v_grid = vGrid;
    Template.Sensor(is).dv_dx = dvDx;
    Template.Sensor(is).bin_weight = binWeight;
    Template.Sensor(is).x_domain = xDomain;
    Template.Sensor(is).selection_mode = selected.mode;
    Template.Sensor(is).point_count = numel(xRel);
    Template.Sensor(is).wide_point_count = numel(xWide);
    Template.Sensor(is).lap_count = numel(selectedRevIds);
    Template.Sensor(is).selected_pulse_indices = selectedPulseIdx(:);
    Template.Sensor(is).selected_revolution_ids = selectedRevIds(:);
    Template.Sensor(is).fit_method = fitMethod;
    Template.Sensor(is).spline_smoothing = P.template.splineSmoothing;
    Template.Sensor(is).eta_app_mm = etaApp(:);
    Template.Sensor(is).eta_median_mm = etaMedian;
    Template.Sensor(is).eta_iqr_mm = etaIqr;
    Template.Sensor(is).eta_limit_mm = etaLimit;

    sensorFile = fullfile(singleTemplateDir, sprintf( ...
        'Template_20241106_B%d_S%d.mat', bladeId, sid));
    save_single_sensor_template_local(sensorFile, Template, Template.Sensor(is));
    legacySensorFile = fullfile(templateDir, sprintf( ...
        'Template_LowSpeedRotating_B%d_CH%d_%s_20241106.mat', ...
        bladeId, sid, P.template.suffix));
    save_single_sensor_template_local(legacySensorFile, Template, Template.Sensor(is));
    sensorTemplateFiles{is} = sensorFile;

    [xProfile, vProfile, gradientAbs, gradientThreshold, effectiveMask] = build_gradient_diagnostic_local( ...
        selected.x_wide, selected.v_wide, baseline, selected.xc, P);
    bladeAudit(is).sensor_id = sid;
    bladeAudit(is).selected_revolution_ids = selectedRevIds(:);
    bladeAudit(is).selected_pulse_indices = selectedPulseIdx(:);
    bladeAudit(is).x_wide = selected.x_wide(:) - xcForCoordinate;
    bladeAudit(is).v_wide = selected.v_wide(:);
    bladeAudit(is).x_selected = selected.x_selected(:) - xcForCoordinate;
    bladeAudit(is).v_selected = selected.v_selected(:);
    bladeAudit(is).weight_selected = selected.weight_selected(:);
    bladeAudit(is).x_profile = xProfile(:) - xcForCoordinate;
    bladeAudit(is).v_profile = vProfile(:);
    bladeAudit(is).gradient_abs = gradientAbs(:);
    bladeAudit(is).gradient_threshold = gradientThreshold;
    bladeAudit(is).effective_mask = effectiveMask(:);
    bladeAudit(is).xc_detected_rel = selected.xc - xcForCoordinate;
    bladeAudit(is).x_domain_rel = xDomain(:).';
    bladeAudit(is).eta_app_mm = etaApp(:);
    bladeAudit(is).eta_median_mm = etaMedian;
    bladeAudit(is).eta_iqr_mm = etaIqr;
    bladeAudit(is).eta_limit_mm = etaLimit;
    bladeAudit(is).stable_plan_mode = string(stablePlanMode);
    bladeAudit(is).selection_mode = string(selected.mode);
    bladeAudit(is).baseline = baseline;
    bladeAudit(is).threshold = threshold;
    bladeAudit(is).template_file = string(sensorFile);

    bladeSummaryRows(is).BladeID = bladeId;
    bladeSummaryRows(is).SensorID = sid;
    bladeSummaryRows(is).SelectedPulseCount = numel(selectedPulseIdx);
    bladeSummaryRows(is).SelectedRevolutionCount = numel(selectedRevIds);
    bladeSummaryRows(is).WidePointCount = numel(xWide);
    bladeSummaryRows(is).TrustPointCount = numel(xRel);
    bladeSummaryRows(is).XDomainLeftMM = xDomain(1);
    bladeSummaryRows(is).XDomainRightMM = xDomain(2);
    bladeSummaryRows(is).XcMM = xcForCoordinate;
    bladeSummaryRows(is).DetectedXcMM = selected.xc;
    bladeSummaryRows(is).ThresholdV = threshold;
    bladeSummaryRows(is).BaselineV = baseline;
    bladeSummaryRows(is).EtaMedianMM = etaMedian;
    bladeSummaryRows(is).EtaIQRMM = etaIqr;
    bladeSummaryRows(is).EtaLimitMM = etaLimit;
    bladeSummaryRows(is).GradientMode = string(P.template.xrangeMode);
    bladeSummaryRows(is).StablePlanMode = string(stablePlanMode);
    bladeSummaryRows(is).TemplateFile = string(sensorFile);
end

bundleFile = fullfile(templateDir, sprintf( ...
    'TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    bladeId, sensor_tag_local(P.sensors.analysis), P.template.suffix));
TemplateBundle = build_template_bundle_local(Template, sensorTemplateFiles, bundleFile);
end

function value = sensor_threshold_local(Sensor_Config, sid, P)
if isfield(P, 'template') && isfield(P.template, 'sensorThresholdOverrideTable')
    value = resolve_sensor_table_value_local(P.template.sensorThresholdOverrideTable, sid, NaN);
    if isfinite(value)
        return;
    end
end
thresholds = Sensor_Config.Sensor_Thresholds;
value = thresholds(sid);
end

function gapPoints = legacy_gap_points_local(Sensor_Config)
gapPoints = Sensor_Config.Gap_Points;
end

function stablePlan = load_stable_window_plan_local(P, bladeId)
stablePlan = struct('has_file', false, 'file', '', 'table', table());
if ~P.template.preferStableWindowPlan
    return;
end
planFile = fullfile(P.routeDir, 'output', 'stable_window_plan', sprintf( ...
    'Step01_CoverageFirstStableWindowPlan_B%d_%s_20241106.csv', ...
    bladeId, sensor_tag_local(P.sensors.analysis)));
if exist(planFile, 'file') ~= 2
    return;
end
T = readtable(planFile);
requiredVars = {'BladeID', 'SensorID', 'LapCount', 'StartOffsetLap'};
for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, T.Properties.VariableNames)
        error('Stable-window plan missing column %s: %s', requiredVars{i}, planFile);
    end
end
T = T(T.BladeID == bladeId & ismember(T.SensorID, P.sensors.analysis), :);
stablePlan.has_file = true;
stablePlan.file = planFile;
stablePlan.table = T;
end

function textValue = stable_plan_file_text_local(stablePlan)
textValue = '';
if isstruct(stablePlan) && isfield(stablePlan, 'has_file') && stablePlan.has_file
    textValue = stablePlan.file;
end
end

function referenceCenter = build_reference_center_table_local(P, ~, ~, ~, ~)
referenceCenter = table(P.sensors.analysis(:), nan(numel(P.sensors.analysis), 1), ...
    'VariableNames', {'sensor_id', 'xc_mm'});
% This route builds each blade template from its own low-speed cloud and
% therefore defaults to the detected center. If a future cross-blade
% reference-center file is introduced, load it here and fill xc_mm.
end

function [selectedRevIds, selectedPulseTimes, modeText] = choose_low_speed_template_pulses_local( ...
        P, bladeId, sid, Sensor_Config, ~, LowSpeedNumbering, ~, stablePlan)

revIds = map_get_numeric_row_local(Sensor_Config.Revolution_IDs, sid);
pulseTimeMat = map_get_numeric_matrix_local(Sensor_Config.Revolution_Pulse_Times, sid);
physicalToLocal = map_get_numeric_row_local(Sensor_Config.Sensor_Physical_To_Local, sid);
localSlot = physicalToLocal(bladeId);
if size(pulseTimeMat, 1) ~= numel(revIds)
    error('CH%d revolution ids and pulse-time rows are inconsistent.', sid);
end

selectedMask = false(size(revIds));
modeText = 'auto_from_low_speed_numbering';
if isstruct(stablePlan) && stablePlan.has_file
    row = stablePlan.table(stablePlan.table.SensorID == sid, :);
    if height(row) == 1
        startOffset = round(row.StartOffsetLap(1));
        lapCount = round(row.LapCount(1));
        selectedRange = (1:lapCount) + startOffset;
        selectedRange = selectedRange(selectedRange >= 1 & selectedRange <= numel(revIds));
        selectedMask(selectedRange) = true;
        modeText = 'stable_window_plan';
    end
end
if ~any(selectedMask)
    sensorRows = LowSpeedNumbering(LowSpeedNumbering.SensorID == sid, :);
    sensorRows = sensorRows(sensorRows.RevolutionID >= min(revIds) & sensorRows.RevolutionID <= max(revIds), :);
    nKeep = min(P.template.defaultCalibrationLaps, height(sensorRows));
    if nKeep < P.template.minSelectedPulseCount
        error('B%d CH%d does not have enough low-speed revolutions to build template.', bladeId, sid);
    end
    selectedMask = ismember(revIds, sensorRows.RevolutionID(1:nKeep));
end

selectedRevIds = revIds(selectedMask);
selectedPulseTimes = pulseTimeMat(selectedMask, localSlot);
valid = isfinite(selectedPulseTimes) & selectedPulseTimes > 0;
selectedPulseTimes = selectedPulseTimes(valid);
selectedRevIds = selectedRevIds(valid);
end

function selectedIdx = match_selected_pulse_times_to_segments_local(segmentArrivalTimes, selectedPulseTimes, sid, bladeId)
selectedIdx = nan(numel(selectedPulseTimes), 1);
if numel(segmentArrivalTimes) < 1
    error('No extracted segments found while matching B%d CH%d selected pulse times.', bladeId, sid);
end
diffArrival = diff(segmentArrivalTimes(:));
bladePeriodEstimate = median(diffArrival(diffArrival > 0), 'omitnan');
if ~isfinite(bladePeriodEstimate) || bladePeriodEstimate <= 0
    bladePeriodEstimate = 1e-3;
end
matchTol = max(5e-4, 0.25 * bladePeriodEstimate);
for i = 1:numel(selectedPulseTimes)
    [bestDt, idx] = min(abs(segmentArrivalTimes(:) - selectedPulseTimes(i)));
    if ~isfinite(bestDt) || bestDt > matchTol
        error(['B%d CH%d selected pulse time could not be matched to extracted segment. ' ...
            'Pulse time = %.9f s, best dt = %.9f s, tol = %.9f s.'], ...
            bladeId, sid, selectedPulseTimes(i), bestDt, matchTol);
    end
    selectedIdx(i) = idx;
end
end

function row = map_get_numeric_row_local(mapObj, sid)
row = mapObj(sid);
row = row(:).';
end

function mat = map_get_numeric_matrix_local(mapObj, sid)
mat = mapObj(sid);
end

function thetaStd = read_opr_center_standard_angle_local(LowSpeedReference, sid, bladeId, Sensor_Config, oprReference)
if isfield(LowSpeedReference, 'standard_angles_opr_center')
    thetaStd = LowSpeedReference.standard_angles_opr_center(sid, bladeId);
    return;
end
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    thetaStd = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, bladeId);
elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
    thetaStd = Sensor_Config.Standard_Relative_Angles(sid, bladeId);
else
    error('No OPR-center standard angle found for CH%d B%d.', sid, bladeId);
end
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
if nargin >= 5 && isstruct(oprReference) && isfield(oprReference, 'phase_shift_deg') && isfinite(oprReference.phase_shift_deg)
    thetaStd = thetaStd - oprReference.phase_shift_deg;
end
end

function selectCfg = build_select_cfg_local(P, threshold, baseline, sensorId)
selectCfg = struct();
selectCfg.sensor_id = sensorId;
selectCfg.threshold = threshold;
selectCfg.baseline = baseline;
selectCfg.weight_floor = P.template.weightFloor;
selectCfg.trust_quantile = P.template.trustQuantile;
selectCfg.trust_edge_margin_mm = P.template.trustEdgeMarginMM;
selectCfg.gradient_energy_quantile = P.template.gradientEnergyQuantile;
selectCfg.gradient_edge_margin_mm = P.template.gradientEdgeMarginMM;
selectCfg.xrange_mode = P.template.xrangeMode;
selectCfg.xrange_gradient_min_ratio = P.template.gradientMinRatio;
selectCfg.xrange_amplitude_min_ratio = P.template.amplitudeMinRatio;
selectCfg.xrange_min_half_width_mm = P.template.minHalfWidthMM;
selectCfg.xrange_max_half_width_mm = P.template.maxHalfWidthMM;
selectCfg.sensor_domain_expand_mm = resolve_sensor_table_value_local( ...
    P.template.sensorDomainExpandMMTable, sensorId, 0);
selectCfg.center_mode = P.template.centerMode;
end

function [xGrid, vGrid, dvDx, binWeight, fitMethod] = fit_template_curve_local(P, xRel, vSel, wSel)
xLeft = min(xRel);
xRight = max(xRel);
xGrid = linspace(xLeft, xRight, P.template.templateGridN).';
edges = linspace(xLeft, xRight, P.template.templateGridN + 1).';
binId = discretize(xRel, edges);
validBin = ~isnan(binId);
vMed = accumarray(binId(validBin), vSel(validBin), [P.template.templateGridN, 1], @median, NaN);
wBin = accumarray(binId(validBin), wSel(validBin), [P.template.templateGridN, 1], @mean, NaN);
binCount = accumarray(binId(validBin), 1, [P.template.templateGridN, 1], @sum, 0);
validGrid = isfinite(vMed) & binCount >= 3;
if nnz(validGrid) < P.template.minValidTemplateBins
    error('Too few valid template bins: %d.', nnz(validGrid));
end
xFit = xGrid(validGrid);
vFit = vMed(validGrid);
wFit = wBin(validGrid);
[xUnique, ~, ic] = unique(round(xFit, 6));
vUnique = accumarray(ic, vFit, [], @median);
wUnique = accumarray(ic, wFit, [], @mean);
if exist('csaps', 'file') == 2
    try
        vGrid = csaps(xUnique, vUnique, P.template.splineSmoothing, xGrid, wUnique);
        fitMethod = 'weighted csaps smoothing spline';
    catch
        vGrid = csaps(xUnique, vUnique, P.template.splineSmoothing, xGrid);
        fitMethod = 'csaps smoothing spline';
    end
else
    vFill = fillmissing(vMed, 'linear', 'EndValues', 'nearest');
    vGrid = smoothdata(vFill, 'sgolay', 41);
    fitMethod = 'median bins + sgolay fallback';
end
vGrid = vGrid(:);
dvDx = gradient(vGrid, xGrid);
binWeight = fillmissing(wBin, 'linear', 'EndValues', 'nearest');
end

function sensor_template = build_single_sensor_template_local(template, sensorStruct)
sensor_template = template;
sensor_template.SensorIDs = sensorStruct.sensor_id;
sensor_template.SensorTag = sprintf('S%d', sensorStruct.sensor_id);
sensor_template.Sensor = sensorStruct;
sensor_template.SourceTemplateMode = 'single_blade_sensor_template';
sensor_template.SourceDataset = '20241106';
end

function save_single_sensor_template_local(filePath, template, sensorStruct)
Template = build_single_sensor_template_local(template, sensorStruct); %#ok<NASGU>
save(filePath, 'Template', '-v7.3');
end

function [etaApp, etaMedian, etaIqr, etaLimit] = estimate_low_speed_eta_scatter_local( ...
        xSel, vSel, xGrid, vGrid, settings)
etaApp = [];
etaMedian = NaN;
etaIqr = NaN;
etaLimit = NaN;
if isempty(xSel) || isempty(vSel) || isempty(xGrid) || isempty(vGrid)
    return;
end
keep = isfinite(xSel(:)) & isfinite(vSel(:));
if nnz(keep) < 5
    return;
end
xUse = xSel(keep);
vUse = vSel(keep);
tpl = struct('x_grid', xGrid(:), 'v_grid', vGrid(:));
xTpl = invert_template_voltage_local_refine_local(tpl, vUse, xUse);
etaApp = xUse - xTpl;
etaApp = etaApp(isfinite(etaApp));
if isempty(etaApp)
    return;
end
etaMedian = median(etaApp, 'omitnan');
etaIqr = iqr(etaApp);
if ~isfield(settings, 'etaLimitMM') || ~isfinite(settings.etaLimitMM) || settings.etaLimitMM <= 0
    return;
end
if ~isfield(settings, 'etaAdaptiveQuantile') || ~isfinite(settings.etaAdaptiveQuantile)
    q = 85;
else
    q = min(max(settings.etaAdaptiveQuantile, 0), 100);
end
if ~isfield(settings, 'etaAdaptiveSafetyFactor') || ~isfinite(settings.etaAdaptiveSafetyFactor)
    safety = 1.2;
else
    safety = max(settings.etaAdaptiveSafetyFactor, 1);
end
if ~isfield(settings, 'etaAdaptiveMinMM') || ~isfinite(settings.etaAdaptiveMinMM)
    minMM = 0.02;
else
    minMM = max(settings.etaAdaptiveMinMM, 0);
end
if ~isfield(settings, 'etaAdaptiveMaxMM') || ~isfinite(settings.etaAdaptiveMaxMM)
    maxMM = settings.etaLimitMM;
else
    maxMM = abs(settings.etaAdaptiveMaxMM);
end
etaCand = safety * prctile(abs(etaApp(:)), q);
if ~isfinite(etaCand) || etaCand <= 0
    etaCand = safety * abs(etaIqr);
end
etaLimit = min(max(etaCand, minMM), max(abs(settings.etaLimitMM), maxMM));
end

function xStat = invert_template_voltage_local_refine_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = nan(size(v));
for i = 1:numel(v)
    vv = v(i);
    if ~isfinite(vv)
        continue;
    end
    diffV = vGrid - vv;
    crossingX = [];
    for k = 1:(numel(vGrid) - 1)
        if diffV(k) == 0
            crossingX(end + 1, 1) = xGrid(k); %#ok<AGROW>
        elseif diffV(k) * diffV(k + 1) <= 0
            denom = vGrid(k + 1) - vGrid(k);
            if abs(denom) < eps
                continue;
            end
            alpha = (vv - vGrid(k)) / denom;
            crossingX(end + 1, 1) = xGrid(k) + alpha * (xGrid(k + 1) - xGrid(k)); %#ok<AGROW>
        end
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function bundle = build_template_bundle_local(template, sensorTemplateFiles, bundleFile)
bundle = struct();
bundle.Route = 'low_speed_rotating_template_sensor_bundle_index';
bundle.TargetBlade = template.TargetBlade;
bundle.SensorIDs = template.SensorIDs;
bundle.SensorTag = template.SensorTag;
bundle.TemplateSuffix = template.TemplateSuffix;
bundle.LowSpeedCase = template.LowSpeedCase;
bundle.LowSpeedDataDir = template.LowSpeedDataDir;
bundle.SensorConfigFile = template.SensorConfigFile;
bundle.LowSpeedFingerprintFile = template.LowSpeedFingerprintFile;
bundle.LowSpeedNumberingFile = template.LowSpeedNumberingFile;
bundle.StableWindowPlanFile = template.StableWindowPlanFile;
bundle.OPRReference = template.OPRReference;
bundle.CreatedBy = template.CreatedBy;
bundle.Settings = template.Settings;
bundle.BundleFile = bundleFile;
bundle.SensorFiles = table(template.SensorIDs(:), string(sensorTemplateFiles(:)), ...
    'VariableNames', {'sensor_id', 'template_file'});
end

function [rawStream, oprTimes, oprReference] = load_low_speed_streams_local(caseDir, cfg, sensorIds)
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR files found in %s.', caseDir);
end
maxSid = max([sensorIds(:); cfg.opr_id]);
rawStream(maxSid) = struct('T', [], 'V', []);
oprT = [];
oprV = [];
lastEnd = [];
for i = 1:numel(fileIds)
    fileId = fileIds(i);
    [tOpr, vOpr] = load_raw_channel_local(caseDir, cfg.opr_id, fileId, cfg.pinlv);
    if isempty(tOpr)
        continue;
    end
    if isempty(lastEnd)
        offset = 0;
    else
        offset = lastEnd + 1 / cfg.pinlv - tOpr(1);
    end
    tOpr = tOpr(:) + offset;
    lastEnd = tOpr(end);
    oprT = [oprT; tOpr(:)]; %#ok<AGROW>
    oprV = [oprV; vOpr(:)]; %#ok<AGROW>

    for sid = sensorIds(:).'
        filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', sid, fileId));
        if exist(filepath, 'file') ~= 2
            continue;
        end
        [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, cfg.pinlv);
        rawStream(sid).T = [rawStream(sid).T; tLocal(:) + offset]; %#ok<AGROW>
        rawStream(sid).V = [rawStream(sid).V; vLocal(:)]; %#ok<AGROW>
    end
end
oprSegments = extract_opr_segments_local(oprT, oprV, cfg.opr_threshold, cfg.gap_points);
oprTimes = oprSegments.arrival_time(:);
oprReference = build_opr_reference_local(oprSegments, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
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

function [tSec, v] = load_raw_channel_local(caseDir, channelId, fileId, fs)
filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
if exist(filepath, 'file') ~= 2
    tSec = [];
    v = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / fs;
v = raw(:, 2);
end

function segments = extract_opr_segments_local(t, v, threshold, gapPoints)
idx = find(v(:) > threshold);
segments = split_segments_to_struct_local(t, v, idx, gapPoints, threshold, false);
end

function segments = extract_sensor_segments_local(t, v, threshold, gapPoints)
if numel(v) >= 21
    vSmooth = sgolayfilt(v(:), 3, 21);
elseif numel(v) >= 5
    vSmooth = smoothdata(v(:), 'movmean', max(3, 2 * floor(numel(v) / 4) + 1));
else
    vSmooth = v(:);
end
idx = find(vSmooth > threshold);
segments = split_segments_to_struct_local(t, v, idx, gapPoints, threshold, true);
end

function segments = split_segments_to_struct_local(t, v, idx, gapPoints, threshold, useCentroid)
segments = struct('start_idx', [], 'end_idx', [], 'start_time', [], ...
    'end_time', [], 'arrival_time', [], 'peak_value', []);
if isempty(idx)
    return;
end
jumps = find(diff(idx) > gapPoints);
segStartPos = [1; jumps(:) + 1];
segEndPos = [jumps(:); numel(idx)];
n = numel(segStartPos);
segments.start_idx = zeros(n, 1);
segments.end_idx = zeros(n, 1);
segments.start_time = zeros(n, 1);
segments.end_time = zeros(n, 1);
segments.arrival_time = zeros(n, 1);
segments.peak_value = zeros(n, 1);
for k = 1:n
    a = idx(segStartPos(k));
    b = idx(segEndPos(k));
    segments.start_idx(k) = a;
    segments.end_idx(k) = b;
    segments.start_time(k) = t(a);
    segments.end_time(k) = t(b);
    [pk, ipk] = max(v(a:b));
    segments.peak_value(k) = pk;
    if useCentroid
        tt = t(a:b);
        vv = max(v(a:b) - threshold, 0);
        if sum(vv) > eps
            segments.arrival_time(k) = sum(tt(:) .* vv(:)) / sum(vv);
        else
            segments.arrival_time(k) = t(a + ipk - 1);
        end
    else
        segments.arrival_time(k) = compute_opr_multithreshold_center_local(t, v, a, b);
    end
end
end

function tCenter = compute_opr_multithreshold_center_local(t, v, a, b)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(numel(v), b + pad);
tt = t(a0:b0);
vv = v(a0:b0);
try
    vs = smooth(vv, 16);
catch
    vs = smoothdata(vv, 'movmean', 16);
end
[peakVal, iPeak] = max(vs);
baseVal = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vs);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] .* max(peakVal - baseVal, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(tt(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = crossing_time_local(tt(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
tCenter = median(centers, 'omitnan');
if ~isfinite(tCenter)
    tCenter = 0.5 * (t(a) + t(b));
end
end

function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < level & v(2:end) >= level, 1, 'first');
else
    idx = find(v(1:end-1) >= level & v(2:end) < level, 1, 'last');
end
if isempty(idx)
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (level - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end

function baseline = estimate_background_baseline_local(v, segments)
mask = true(size(v));
for k = 1:numel(segments.start_idx)
    a = segments.start_idx(k);
    b = segments.end_idx(k);
    pad = max(10, round(1.5 * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(v), b + pad);
    mask(aa:bb) = false;
end
if nnz(mask) > 100
    baseline = median(v(mask), 'omitnan');
else
    baseline = median(v(v <= prctile(v, 30)), 'omitnan');
end
end

function F = build_phase_speed_local(oprTimes, pulsesPerRev)
spdT = oprTimes(1:(end - pulsesPerRev));
spdV = 360 ./ max(oprTimes((pulsesPerRev + 1):end) - oprTimes(1:(end - pulsesPerRev)), eps);
F = griddedInterpolant(spdT, spdV, 'linear', 'nearest');
end

function thetaPoints = map_segment_to_relative_angle_local(tRef, tSeg, F_omega_deg)
thetaPoints = zeros(size(tSeg));
for i = 1:numel(tSeg)
    tGrid = linspace(tRef, tSeg(i), 10);
    thetaPoints(i) = trapz(tGrid, F_omega_deg(tGrid));
end
end

function [xWide, vWide] = restore_low_speed_point_cloud_local( ...
        R, segments, selectedIdx, oprTimes, F_omega_deg, thetaStd, rTipMM, expandFactor)
xWide = [];
vWide = [];
for ii = 1:numel(selectedIdx)
    pidx = selectedIdx(ii);
    a = segments.start_idx(pidx);
    b = segments.end_idx(pidx);
    pad = max(8, round(expandFactor * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(R.T), b + pad);
    tSeg = R.T(aa:bb);
    vSeg = R.V(aa:bb);
    tArrival = segments.arrival_time(pidx);
    idxPrev = find(oprTimes < tArrival, 1, 'last');
    if isempty(idxPrev)
        continue;
    end
    thetaPointsDeg = map_segment_to_relative_angle_local(oprTimes(idxPrev), tSeg, F_omega_deg);
    thetaDiffDeg = mod(thetaPointsDeg - thetaStd + 180, 360) - 180;
    xPointsMM = thetaDiffDeg * (pi / 180) * rTipMM;
    xWide = [xWide; xPointsMM(:)]; %#ok<AGROW>
    vWide = [vWide; vSeg(:)]; %#ok<AGROW>
end
end

function selected = select_sg_like_template_points_local(xWide, vWide, cfg)
finite = isfinite(xWide) & isfinite(vWide);
xWide = xWide(finite);
vWide = vWide(finite);
vZero = max(vWide - cfg.baseline, 0);
if numel(xWide) < 50
    error('Too few low-speed cloud points for template selection.');
end
highLevel = prctile(vZero, 85);
highMask = vZero >= highLevel & vZero > 0;
if nnz(highMask) >= 5 && sum(vZero(highMask)) > eps
    xc = sum(xWide(highMask) .* vZero(highMask)) / sum(vZero(highMask));
else
    [~, imax] = max(vZero);
    xc = xWide(imax);
end
strictMask = vWide >= cfg.threshold;
dxStrict = xWide(strictMask) - xc;
if nnz(strictMask) < 30
    strictMask = vZero >= prctile(vZero, 60);
    dxStrict = xWide(strictMask) - xc;
end
[leftL, rightL, xrangeOk, xrangeModeUsed] = estimate_xrange_span_local(xWide, vWide, cfg.baseline, xc, cfg);
if ~xrangeOk
    leftDist = -dxStrict(dxStrict < 0);
    rightDist = dxStrict(dxStrict > 0);
    if isempty(leftDist) || isempty(rightDist)
        L = prctile(abs(dxStrict), 95);
        leftL = L;
        rightL = L;
    else
        leftL = prctile(leftDist, 100 * cfg.trust_quantile);
        rightL = prctile(rightDist, 100 * cfg.trust_quantile);
    end
    leftL = max(leftL - cfg.trust_edge_margin_mm, 0.1);
    rightL = max(rightL - cfg.trust_edge_margin_mm, 0.1);
end
if xrangeOk
    trustMask = xWide >= xc - leftL & xWide <= xc + rightL;
else
    trustMask = strictMask & xWide >= xc - leftL & xWide <= xc + rightL;
end
if nnz(trustMask) < 30
    trustMask = strictMask;
end
if isfield(cfg, 'sensor_domain_expand_mm') && isfinite(cfg.sensor_domain_expand_mm) && ...
        abs(cfg.sensor_domain_expand_mm) > 0
    trustMask = trustMask | (xWide >= xc - leftL - cfg.sensor_domain_expand_mm & ...
        xWide <= xc + rightL + cfg.sensor_domain_expand_mm);
end
xSelected = xWide(trustMask);
vSelected = vWide(trustMask);
if strcmpi(cfg.center_mode, 'sgfit')
    xc = refine_center_by_sg_fit_local(xSelected, vSelected, cfg.baseline, xc);
end
weightSelected = build_gradient_weight_from_wide_cloud_local(xSelected, xWide, vWide, cfg.baseline, cfg.weight_floor);
selected = struct();
if xrangeOk
    selected.mode = ['raw_low_speed_' xrangeModeUsed '_gate_gradient_weight'];
else
    selected.mode = 'raw_low_speed_sg_like_strict_gate_trust_window_gradient_weight';
end
if isfield(cfg, 'sensor_domain_expand_mm') && isfinite(cfg.sensor_domain_expand_mm) && ...
        abs(cfg.sensor_domain_expand_mm) > 0
    if cfg.sensor_domain_expand_mm > 0
        selected.mode = [selected.mode '_plus_domain_expand'];
    else
        selected.mode = [selected.mode '_plus_domain_shrink'];
    end
end
selected.xc = xc;
selected.x_wide = xWide;
selected.v_wide = vWide;
selected.x_selected = xSelected;
selected.v_selected = vSelected;
selected.weight_selected = weightSelected;
end

function value = resolve_sensor_table_value_local(tableNow, sensorId, defaultValue)
value = defaultValue;
if nargin < 3
    defaultValue = 0;
    value = defaultValue;
end
if isempty(tableNow) || size(tableNow, 2) < 2
    return;
end
idx = find(tableNow(:, 1) == sensorId, 1);
if isempty(idx)
    return;
end
candidate = tableNow(idx, 2);
if isfinite(candidate)
    value = candidate;
end
end

function [leftL, rightL, ok, modeUsed] = estimate_xrange_span_local(x, v, baseline, xc, cfg)
if strcmpi(cfg.xrange_mode, 'threshold')
    [leftL, rightL, ok] = estimate_gradient_threshold_span_local( ...
        x, v, baseline, xc, cfg.xrange_gradient_min_ratio, cfg.xrange_amplitude_min_ratio, ...
        cfg.xrange_min_half_width_mm, cfg.xrange_max_half_width_mm);
    modeUsed = 'gradient_threshold';
else
    [leftL, rightL, ok] = estimate_gradient_energy_span_local( ...
        x, v, baseline, xc, cfg.gradient_energy_quantile, cfg.gradient_edge_margin_mm);
    modeUsed = 'gradient_energy';
end
end

function [leftL, rightL, ok] = estimate_gradient_threshold_span_local( ...
        x, v, baseline, xc, gradientMinRatio, amplitudeMinRatio, minHalfWidth, maxHalfWidth)
leftL = NaN;
rightL = NaN;
ok = false;
[xGrid, vSmooth] = build_denoised_static_profile_local(x, v, baseline);
if numel(xGrid) < 20
    return;
end
vAboveBase = max(vSmooth(:), 0);
peakV = max(vAboveBase, [], 'omitnan');
if ~isfinite(peakV) || peakV <= 0
    return;
end
gAbs = abs(gradient(vSmooth(:), xGrid(:)));
gThreshold = estimate_noise_aware_gradient_threshold_local(gAbs, vAboveBase, gradientMinRatio);
effective = gAbs >= gThreshold;
if amplitudeMinRatio > 0
    effective = effective & vAboveBase >= amplitudeMinRatio * peakV;
end
leftCandidates = find(xGrid(:) < xc & effective(:));
rightCandidates = find(xGrid(:) > xc & effective(:));
if isempty(leftCandidates) || isempty(rightCandidates)
    return;
end
leftIdx = leftCandidates(1);
rightIdx = rightCandidates(end);
leftL = max(xc - xGrid(leftIdx), 0);
rightL = max(xGrid(rightIdx) - xc, 0);
leftL = max(leftL, minHalfWidth);
rightL = max(rightL, minHalfWidth);
if isfinite(maxHalfWidth) && maxHalfWidth > 0
    leftL = min(leftL, maxHalfWidth);
    rightL = min(rightL, maxHalfWidth);
end
ok = isfinite(leftL) && isfinite(rightL) && leftL > 0 && rightL > 0;
end

function [leftL, rightL, ok] = estimate_gradient_energy_span_local( ...
        x, v, baseline, xc, energyQuantile, edgeMarginMM)
leftL = NaN;
rightL = NaN;
ok = false;
[xGrid, vSmooth] = build_denoised_static_profile_local(x, v, baseline);
if numel(xGrid) < 20 || max(vSmooth) <= 0
    return;
end
g = abs(gradient(vSmooth, xGrid));
signalGate = vSmooth >= prctile(vSmooth, 40);
vAboveBase = max(vSmooth(:), 0);
peakV = max(vAboveBase, [], 'omitnan');
if ~isfinite(peakV) || peakV <= 0
    return;
end
energy = g(:) .* signalGate(:) + 0.05 * max(g(:)) * (vAboveBase ./ peakV);
if sum(energy) <= eps
    return;
end
tail = (1 - energyQuantile) / 2;
cum = cumsum(energy) ./ sum(energy);
[cumUnique, ia] = unique(cum, 'stable');
xUnique = xGrid(ia);
if numel(cumUnique) < 2
    return;
end
xLo = interp1(cumUnique, xUnique, tail, 'linear', 'extrap');
xHi = interp1(cumUnique, xUnique, 1 - tail, 'linear', 'extrap');
if ~isfinite(xLo) || ~isfinite(xHi) || xHi <= xLo || xc <= xLo || xc >= xHi
    return;
end
leftL = max(xc - xLo - edgeMarginMM, 0.1);
rightL = max(xHi - xc - edgeMarginMM, 0.1);
ok = true;
end

function [xGrid, vSmooth] = build_denoised_static_profile_local(x, v, baseline)
finite = isfinite(x) & isfinite(v);
x = x(finite);
vZero = max(v(finite) - baseline, 0);
xGrid = [];
vSmooth = [];
if numel(x) < 50 || max(vZero) <= 0
    return;
end
[xSort, idx] = sort(x);
vSort = vZero(idx);
gridN = min(700, max(120, round(numel(xSort) / 120)));
edges = linspace(min(xSort), max(xSort), gridN + 1).';
xGrid = 0.5 * (edges(1:end-1) + edges(2:end));
binId = discretize(xSort, edges);
valid = ~isnan(binId);
vMed = accumarray(binId(valid), vSort(valid), [gridN, 1], @median, NaN);
if nnz(isfinite(vMed)) < 15
    xGrid = [];
    vSmooth = [];
    return;
end
vFill = fillmissing(vMed, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(vFill) - 1) / 2) + 1);
if span >= 5
    vSmooth = smoothdata(vFill, 'sgolay', span);
else
    vSmooth = vFill;
end
end

function gThreshold = estimate_noise_aware_gradient_threshold_local(gAbs, vAboveBase, ratio)
gAbs = gAbs(:);
vAboveBase = vAboveBase(:);
finite = isfinite(gAbs) & isfinite(vAboveBase);
gAbs = gAbs(finite);
vAboveBase = vAboveBase(finite);
if isempty(gAbs)
    gThreshold = Inf;
    return;
end
peakV = max(vAboveBase, [], 'omitnan');
if ~isfinite(peakV) || peakV <= 0
    gThreshold = prctile(gAbs, 95);
    return;
end
noiseGate = vAboveBase <= 0.10 * peakV;
signalGate = vAboveBase >= 0.20 * peakV;
if nnz(noiseGate) >= 10
    gNoise = prctile(gAbs(noiseGate), 95);
else
    gNoise = prctile(gAbs, 10);
end
if nnz(signalGate) >= 10
    gSignal = prctile(gAbs(signalGate), 95);
else
    gSignal = prctile(gAbs, 95);
end
if ~isfinite(gNoise)
    gNoise = 0;
end
if ~isfinite(gSignal) || gSignal <= gNoise
    gSignal = max(gAbs);
end
gThreshold = gNoise + ratio * max(gSignal - gNoise, 0);
end

function xc = refine_center_by_sg_fit_local(x, v, baseline, xcInitial)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
if numel(x) < 50
    xc = xcInitial;
    return;
end
[peakVal, idxPeak] = max(v);
B0 = max(peakVal - baseline, 0.1);
w0 = max(std(x), 0.2);
n0 = 3.0;
xc0 = xcInitial;
if ~isfinite(xc0)
    xc0 = x(idxPeak);
end
sgModel = @(p, xx) p(1) .* exp(-abs((xx - p(4)) ./ max(p(2), 1e-6)).^max(p(3), 1e-6)) + baseline;
objFun = @(p) sum((v - sgModel(p, x)).^2);
p0 = [B0, w0, n0, xc0];
lb = [0.05, 0.05, 1.2, min(x) - 0.5];
ub = [max(10, 2 * B0 + 0.5), 6.0, 10.0, max(x) + 0.5];
try
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        pOpt = fmincon(objFun, p0, [], [], [], [], lb, ub, [], opts);
    else
        penaltyObj = @(p) objFun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
        pOpt = fminsearch(penaltyObj, p0, opts);
        pOpt = min(max(pOpt, lb), ub);
    end
    if isfinite(pOpt(4))
        xc = pOpt(4);
    else
        xc = xcInitial;
    end
catch
    xc = xcInitial;
end
end

function W = build_gradient_weight_from_wide_cloud_local(xQuery, xWide, vWide, baseline, weightFloor)
finite = isfinite(xWide) & isfinite(vWide);
xWide = xWide(finite);
vZero = max(vWide(finite) - baseline, 0);
if numel(xWide) < 20 || all(vZero == 0)
    W = ones(size(xQuery));
    return;
end
[xSort, idx] = sort(xWide);
vSort = vZero(idx);
gridN = min(500, max(80, round(numel(xSort) / 200)));
edges = linspace(min(xSort), max(xSort), gridN + 1).';
xGrid = 0.5 * (edges(1:end-1) + edges(2:end));
binId = discretize(xSort, edges);
valid = ~isnan(binId);
vMed = accumarray(binId(valid), vSort(valid), [gridN, 1], @median, NaN);
validGrid = isfinite(vMed);
if nnz(validGrid) < 10
    W = ones(size(xQuery));
    return;
end
vFill = fillmissing(vMed, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(vFill) - 1) / 2) + 1);
if span >= 5
    vSmooth = smoothdata(vFill, 'sgolay', span);
else
    vSmooth = vFill;
end
dVdx = abs(gradient(vSmooth, xGrid));
wGrid = weightFloor + (1 - weightFloor) * dVdx ./ max(dVdx + eps);
W = interp1(xGrid, wGrid, xQuery(:), 'linear', weightFloor);
W(~isfinite(W)) = weightFloor;
W = max(W, weightFloor);
end

function oprReference = build_opr_reference_local(oprSegments, pulsesPerRev, rTipMM)
oprReference = struct('mode', 'multi_threshold_center', ...
    'standard_angle_reference', 'opr_pulse_center', ...
    'phase_shift_deg', 0, ...
    'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if isempty(oprSegments.arrival_time) || isempty(oprSegments.start_time) || pulsesPerRev < 1
    return;
end
centerTime = oprSegments.arrival_time(:);
startTime = oprSegments.start_time(:);
n = min(numel(centerTime) - pulsesPerRev, numel(startTime));
if n < 1
    return;
end
dtCenter = centerTime(1:n) - startTime(1:n);
dtRev = centerTime((1:n) + pulsesPerRev) - centerTime(1:n);
valid = isfinite(dtCenter) & isfinite(dtRev) & dtRev > eps;
if ~any(valid)
    return;
end
shiftDeg = 360 * dtCenter(valid) ./ dtRev(valid);
oprReference.phase_shift_deg = median(shiftDeg, 'omitnan');
oprReference.phase_shift_mm = oprReference.phase_shift_deg * (pi / 180) * rTipMM;
oprReference.median_center_minus_start_s = median(dtCenter(valid), 'omitnan');
end

function [xProfile, vProfile, gradientAbs, gradientThreshold, effectiveMask] = build_gradient_diagnostic_local( ...
        xWide, vWide, baseline, xc, P)
[xProfile, vProfile] = build_denoised_static_profile_local(xWide, vWide, baseline);
if isempty(xProfile)
    gradientAbs = [];
    gradientThreshold = NaN;
    effectiveMask = [];
    return;
end
gradientAbs = abs(gradient(vProfile(:), xProfile(:)));
vAboveBase = max(vProfile(:), 0);
if strcmpi(P.template.xrangeMode, 'threshold')
    gradientThreshold = estimate_noise_aware_gradient_threshold_local( ...
        gradientAbs, vAboveBase, P.template.gradientMinRatio);
    effectiveMask = gradientAbs >= gradientThreshold;
    if P.template.amplitudeMinRatio > 0
        peakV = max(vAboveBase, [], 'omitnan');
        if isfinite(peakV) && peakV > 0
            effectiveMask = effectiveMask & vAboveBase >= P.template.amplitudeMinRatio * peakV;
        end
    end
else
    gradientThreshold = prctile(gradientAbs, 75);
    effectiveMask = gradientAbs >= gradientThreshold;
end
if isfinite(xc)
    effectiveMask = effectiveMask & isfinite(xProfile(:)) & abs(xProfile(:) - xc) <= 1.5 * P.template.maxHalfWidthMM;
end
end

function visualize_low_speed_template_summary_local(T, P, figDir)
if isempty(T)
    return;
end
fig = figure('Name', 'Step05 low-speed template summary', 'Color', 'w', ...
    'Position', [100, 100, 1300, 780], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
gscatter(T.SensorID, T.TrustPointCount, T.BladeID);
xlabel('Sensor ID');
ylabel('Trust-cloud point count');
title('Selected template points by blade and sensor');
grid on; box on;

nexttile;
scatter(T.XDomainLeftMM, T.XDomainRightMM, 40, T.BladeID, 'filled');
xlabel('Domain left (mm)');
ylabel('Domain right (mm)');
title('Template domain span');
grid on; box on; colorbar;

nexttile;
bar(categorical(compose('B%d-CH%d', T.BladeID, T.SensorID)), [T.DetectedXcMM, T.XcMM]);
ylabel('Center location (mm)');
title('Detected center vs final coordinate center');
legend({'Detected xc', 'Final xc'}, 'Location', 'best');
grid on; box on;

nexttile;
axis off;
summaryText = evalc('disp(T(:, {''BladeID'',''SensorID'',''SelectedPulseCount'',''SelectedRevolutionCount'',''XDomainLeftMM'',''XDomainRightMM'',''StablePlanMode''}))');
text(0, 1, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Consolas', 'FontSize', 8.5, 'Interpreter', 'none');
title('Template build table');

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step05_LowSpeedTemplateLibrary_20241106.png'), 'Resolution', 300);
end
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end
