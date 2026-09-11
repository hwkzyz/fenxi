%% Step06: build gap-aware DynamicMap for 20241106
% This script keeps the verified continuous raw waveform window extraction
% from the old direct Step06K route, but repackages it as the main
% gap-aware DynamicMap input for the new Step07J/Step07K split workflow.
%
% DynamicMap.Window(w).Sensor(s) keeps the mature direct-style structure:
%   DynamicMap.Window(w).Sensor(s).t / x_abs / x_rel / V / W / theta
%
% All runtime calibration inputs are resolved from this route's outputs.
% The raw measurement folder is recorded separately through Step02.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
packageCfg = Config_20241106();
outDir = packageCfg.paths.gapRuntime;
figDir = fullfile(packageCfg.paths.figures, ...
    'step06_gapaware_dynamic_map_20241106');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
Flow = ProjectionFlow_Config_20241106();

%% Parameters to tune
P = struct();
P.targetBlade = Flow.identification.targetBlade;
P.analysisSensors = Flow.identification.analysisSensors;
P.gapSensors = Flow.identification.gapSensors;
P.analysisStartTimeSec = Flow.identification.analysisStartTimeSec;
P.targetLaps = Flow.identification.targetBladePasses;
P.windowLaps = Flow.identification.windowBladePasses;
P.slidingStepLaps = Flow.identification.slidingStepBladePasses;
P.pulseHalfWindowSec = 3.0e-4;  % raw waveform window around each blade pass
P.edgeFractionForBaseline = 0.18;
P.minFinitePointsPerPulse = 100;
P.blockStepSamples = 1e7;
P.blockLabelStep = 1000;
P.rawCacheMaxEntries = 10;
P.templateFile = '';
P.correctedLibFile = '';
P.templateBankFile = fullfile(packageCfg.paths.calibrationGap, ...
    Flow.calibration.templateBankFile);
P.gapBankFile = fullfile(packageCfg.paths.calibrationGap, ...
    Flow.calibration.gapBankFile);
P.rotatingDynamicMapFile = '';
P.directReferenceFile = '';

P.targetBlade = parse_scalar_env_local('STEP06G_TARGET_BLADE', P.targetBlade);
P.analysisSensors = parse_int_env_local('STEP06G_ANALYSIS_SENSORS', P.analysisSensors);
P.gapSensors = parse_int_env_local('STEP06G_GAP_SENSORS', P.gapSensors);
P.gapSensors = intersect(P.gapSensors(:).', P.analysisSensors(:).', 'stable');
assert_capacitive_gap_sensors_local(P.gapSensors, 'STEP06G_ALLOW_NONCAP_GAP_SENSORS');
P.analysisStartTimeSec = parse_float_env_local('STEP06G_START_TIME_SEC', P.analysisStartTimeSec);
P.targetLaps = round(parse_scalar_env_local('STEP06G_TARGET_LAPS', P.targetLaps));
P.windowLaps = round(parse_scalar_env_local('STEP06G_WINDOW_LAPS', P.windowLaps));
P.slidingStepLaps = round(parse_scalar_env_local('STEP06G_SLIDING_STEP_LAPS', P.slidingStepLaps));
P.pulseHalfWindowSec = parse_scalar_env_local('STEP06G_PULSE_HALF_WINDOW_SEC', P.pulseHalfWindowSec);
outputSuffix = sanitize_suffix_local(strtrim(getenv('STEP06G_OUTPUT_SUFFIX')));

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
gapSensorTag = ['S', sprintf('%d', P.gapSensors)];
if isempty(P.correctedLibFile)
    P.correctedLibFile = fullfile(packageCfg.paths.calibrationGap, sprintf( ...
        'Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
        P.targetBlade, gapSensorTag));
end

step02File = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
if ~isfile(step02File)
    error('Run Step02 first. Missing file: %s', step02File);
end

S02 = load(step02File);
if ~isfield(S02, 'bladeNumbering') || ...
        ~isfield(S02.bladeNumbering, 'scheme') || ...
        ~strcmpi(S02.bladeNumbering.scheme, 'foundation_physical_blade_id')
    error(['Step02 cache does not use Foundation physical blade IDs. ' ...
        'Rerun Step02_Calc_BTT_Displacement_20241106.m before Step06.']);
end
if ~isfield(S02.bladeNumbering, 'physicalSensorIds') || ...
        ~all(ismember(P.analysisSensors, S02.bladeNumbering.physicalSensorIds))
    error('Step02 cache lacks Foundation physical blade labels for sensors %s.', ...
        mat2str(P.analysisSensors));
end
[Template, templateSourceLabel] = load_template_for_run_local(P);
[CorrectedGapLibrary, correctedSourceLabel] = load_gap_library_for_run_local(P, Template);
if isempty(CorrectedGapLibrary)
    if ~isfile(P.correctedLibFile)
        error(['Run Step06I_Calibrate_AllBladeGapLibrary_20241106 first. ' ...
            'Missing gap bank %s and fallback combined library %s.'], ...
            P.gapBankFile, P.correctedLibFile);
    end
    Sc = load(P.correctedLibFile, 'CorrectedGapLibrary');
    CorrectedGapLibrary = Sc.CorrectedGapLibrary;
    correctedSourceLabel = P.correctedLibFile;
    if isfield(CorrectedGapLibrary, 'lowSpeedTemplate') && isstruct(CorrectedGapLibrary.lowSpeedTemplate)
        gapTemplate = clean_template_local(CorrectedGapLibrary.lowSpeedTemplate, P.gapSensors);
        Template = merge_gap_template_sensors_local(Template, gapTemplate, P.gapSensors);
    end
end
if isempty(Template)
    if isfield(CorrectedGapLibrary, 'lowSpeedTemplate') && isstruct(CorrectedGapLibrary.lowSpeedTemplate)
        Template = clean_template_local(CorrectedGapLibrary.lowSpeedTemplate, P.analysisSensors);
        templateSourceLabel = '[CorrectedGapLibrary.lowSpeedTemplate]';
    else
        if ~isfile(P.templateFile)
            error('Missing low-speed template bank and fallback template: %s', P.templateFile);
        end
        S = load(P.templateFile, 'Template');
        Template = clean_template_local(S.Template, P.analysisSensors);
        templateSourceLabel = P.templateFile;
    end
else
    Template = clean_template_local(Template, P.analysisSensors);
end
RotatingDynamicMap = [];
if isfile(P.rotatingDynamicMapFile)
    tmp = load(P.rotatingDynamicMapFile, 'DynamicMap');
    RotatingDynamicMap = tmp.DynamicMap;
end
DirectReference = [];
if isfile(P.directReferenceFile)
    tmp = load(P.directReferenceFile, 'IdentificationResult');
    if isfield(tmp, 'IdentificationResult')
        DirectReference = tmp.IdentificationResult;
    end
end

fprintf('\n=== Step06 20241106 gap-aware DynamicMap ===\n');
fprintf('Low-speed template: %s\n', templateSourceLabel);
fprintf('Corrected gap library: %s\n', correctedSourceLabel);
fprintf('High-speed raw folder: %s\n', S02.highDir);
fprintf('Physical blade B%d, sensors %s\n', P.targetBlade, mat2str(P.analysisSensors));
fprintf('Gap sensors %s; direct-only sensors %s\n', ...
    mat2str(P.gapSensors), mat2str(setdiff(P.analysisSensors, P.gapSensors, 'stable')));
fprintf('Start time %.6f s, target laps %d, window %d laps, step %d lap(s)\n', ...
    P.analysisStartTimeSec, P.targetLaps, P.windowLaps, P.slidingStepLaps);
fprintf('Planned sliding windows: %d\n', planned_window_count_local(P.targetLaps, ...
    P.windowLaps, P.slidingStepLaps));
fprintf('Raw pulse window: Foundation jilublade(:,1:2) bounds\n');

opr = load_foundation_opr_local(S02.bladeNumbering.oprFile);
rawCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

LapData = repmat(struct('sensor_id', NaN, 'Lap', []), numel(P.analysisSensors), 1);
selectionRows = cell(numel(P.analysisSensors), 1);
globalTime = [inf, -inf];

for is = 1:numel(P.analysisSensors)
    sid = P.analysisSensors(is);
    sensorIdx = find([S02.highDisp.sensorId] == sid, 1, 'first');
    if isempty(sensorIdx)
        error('No Step02 high-speed displacement data for sensor %d.', sid);
    end
    if P.targetBlade > numel(S02.highDisp(sensorIdx).slot)
        error('Physical blade %d is not available for sensor %d.', P.targetBlade, sid);
    end

    slot = S02.highDisp(sensorIdx).slot(P.targetBlade);
    if ~isfield(slot, 'physicalBladeId') || ...
            slot.physicalBladeId ~= P.targetBlade
        error('CH%d slot %d is not verified as physical blade B%d.', ...
            sid, P.targetBlade, P.targetBlade);
    end
    eventTimesAll = slot.time(:);
    if ~isfield(slot, 'pulseStart') || ~isfield(slot, 'pulseEnd') || ...
            numel(slot.pulseStart) ~= numel(eventTimesAll) || ...
            numel(slot.pulseEnd) ~= numel(eventTimesAll)
        error('CH%d physical B%d lacks Foundation pulse boundaries.', sid, P.targetBlade);
    end
    selectedIndex = find(eventTimesAll >= P.analysisStartTimeSec & ...
        isfinite(eventTimesAll) & isfinite(slot.pulseStart(:)) & ...
        isfinite(slot.pulseEnd(:)));
    eventTimes = eventTimesAll(selectedIndex);
    if numel(eventTimes) < P.targetLaps
        error('P%d has only %d blade passes after %.6f s.', ...
            sid, numel(eventTimes), P.analysisStartTimeSec);
    end
    selectedIndex = selectedIndex(1:P.targetLaps);
    eventTimes = eventTimesAll(selectedIndex);
    pulseStart = slot.pulseStart(selectedIndex);
    pulseEnd = slot.pulseEnd(selectedIndex);
    selectionRows{is} = eventTimes(:);
    globalTime(1) = min(globalTime(1), min(pulseStart));
    globalTime(2) = max(globalTime(2), max(pulseEnd));

    Lap = repmat(struct('lap_id', NaN, 'event_time', NaN, ...
        'pulse_start', NaN, 'pulse_end', NaN, 't', [], ...
        'x_abs', [], 'V', [], 'W', [], 'theta', [], 'rot_freq_hz', NaN), ...
        P.targetLaps, 1);

    for lapId = 1:P.targetLaps
        tEvent = eventTimes(lapId);
        tPulseStart = pulseStart(lapId);
        tPulseEnd = pulseEnd(lapId);
        rotFreqHz = local_rot_freq_local(opr.tCenter, tEvent);
        tipSpeedMmS = 2 * pi * S02.rTipMm * rotFreqHz;

        [tRaw, vRaw, ok] = read_raw_time_window_local(S02.highDir, sid, ...
            tPulseStart, tPulseEnd, S02.sampleRateHz, P.blockStepSamples, ...
            P.blockLabelStep, rawCache, P.rawCacheMaxEntries);
        if ~ok
            continue;
        end

        xAbs = (tRaw - tEvent) * tipSpeedMmS;
        edgeCount = max(5, round(P.edgeFractionForBaseline * numel(vRaw)));
        baselineV = median([vRaw(1:edgeCount); vRaw(end-edgeCount+1:end)], 'omitnan');
        vMv = 1000 * (vRaw - baselineV);
        theta = local_rot_phase_local(opr.tCenter, tRaw);
        finite = isfinite(tRaw) & isfinite(xAbs) & isfinite(vMv) & isfinite(theta);
        if nnz(finite) < P.minFinitePointsPerPulse
            continue;
        end

        tRaw = tRaw(finite);
        xAbs = xAbs(finite);
        vMv = vMv(finite);
        theta = theta(finite);
        W = build_weight_local(tRaw, vMv);

        Lap(lapId).lap_id = lapId;
        Lap(lapId).event_time = tEvent;
        Lap(lapId).pulse_start = tPulseStart;
        Lap(lapId).pulse_end = tPulseEnd;
        Lap(lapId).t = tRaw(:);
        Lap(lapId).x_abs = xAbs(:);
        Lap(lapId).V = vMv(:);
        Lap(lapId).W = W(:);
        Lap(lapId).theta = theta(:);
        Lap(lapId).rot_freq_hz = rotFreqHz;
    end

    LapData(is).sensor_id = sid;
    LapData(is).Lap = Lap;
end

numWindows = planned_window_count_local(P.targetLaps, P.windowLaps, P.slidingStepLaps);
Window = repmat(struct('window_id', NaN, 'lap_range', [], ...
    'event_time_window', [], 'time_window', [], 'rot_freq_mean_hz', NaN, ...
    'rot_rpm_mean', NaN, 'Sensor', []), numWindows, 1);

for iw = 1:numWindows
    lapStart = 1 + (iw - 1) * P.slidingStepLaps;
    lapEnd = lapStart + P.windowLaps - 1;
    lapRange = lapStart:lapEnd;

    Sensor = repmat(struct('sensor_id', NaN, 't', [], 'x_abs', [], ...
        'x_rel', [], 'V', [], 'W', [], 'theta', [], 'event_time', [], ...
        'point_count', NaN), numel(P.analysisSensors), 1);
    allT = [];
    allEventT = [];
    allRotFreq = [];

    for is = 1:numel(P.analysisSensors)
        sid = P.analysisSensors(is);
        Tpl = Template.Sensor(is);
        t = [];
        xAbs = [];
        V = [];
        W = [];
        theta = [];
        eventTime = [];
        rotFreq = [];

        for lapId = lapRange
            D = LapData(is).Lap(lapId);
            t = [t; D.t(:)]; %#ok<AGROW>
            xAbs = [xAbs; D.x_abs(:)]; %#ok<AGROW>
            V = [V; D.V(:)]; %#ok<AGROW>
            W = [W; D.W(:)]; %#ok<AGROW>
            theta = [theta; D.theta(:)]; %#ok<AGROW>
            eventTime = [eventTime; repmat(D.event_time, numel(D.t), 1)]; %#ok<AGROW>
            rotFreq = [rotFreq; D.rot_freq_hz]; %#ok<AGROW>
        end

        xRel = xAbs - Tpl.xc;
        Sensor(is).sensor_id = sid;
        Sensor(is).t = t(:);
        Sensor(is).x_abs = xAbs(:);
        Sensor(is).x_rel = xRel(:);
        Sensor(is).V = V(:);
        Sensor(is).W = normalize_weight_local(W(:));
        Sensor(is).theta = theta(:);
        Sensor(is).event_time = eventTime(:);
        Sensor(is).point_count = numel(t);
        allT = [allT; t(:)]; %#ok<AGROW>
        allEventT = [allEventT; eventTime(:)]; %#ok<AGROW>
        allRotFreq = [allRotFreq; rotFreq(:)]; %#ok<AGROW>
    end

    Window(iw).window_id = iw;
    Window(iw).lap_range = lapRange;
    Window(iw).event_time_window = [min(allEventT), max(allEventT)];
    Window(iw).time_window = [min(allT), max(allT)];
    Window(iw).rot_freq_mean_hz = mean(allRotFreq, 'omitnan');
    Window(iw).rot_rpm_mean = 60 * Window(iw).rot_freq_mean_hz;
    Window(iw).Sensor = Sensor;
end

DynamicMap = struct();
DynamicMap.dataset = '20241106';
DynamicMap.method = 'gap_aware_continuous_raw_waveform_window_map';
DynamicMap.description = ['Continuous raw high-speed waveform windows for the 20241106 gap-aware main route. ' ...
    'Runtime template, gap calibration, Step02 and DynamicMap inputs are folder-local.'];
DynamicMap.TemplateFile = templateSourceLabel;
DynamicMap.CorrectedGapLibraryFile = correctedSourceLabel;
DynamicMap.DirectReferenceFile = '';
DynamicMap.RotatingDynamicMapFile = '';
DynamicMap.Step02File = step02File;
DynamicMap.lowDir = S02.lowDir;
DynamicMap.highDir = S02.highDir;
DynamicMap.SourceSettings = P;
DynamicMap.analysisSensors = P.analysisSensors(:).';
DynamicMap.targetBlade = P.targetBlade;
DynamicMap.globalTimeWindow = globalTime;
DynamicMap.selectionEventTimesBySensor = selectionRows;
DynamicMap.Route = '20241106_gap_tilt_main';
DynamicMap.SensorIDs = P.analysisSensors(:).';
DynamicMap.SensorTag = sensorTag;
DynamicMap.TargetBlade = P.targetBlade;
DynamicMap.BladeNumbering = S02.bladeNumbering;
DynamicMap.SourceResultFile = step02File;
DynamicMap.SourceResultSensorTag = sensorTag;
DynamicMap.TemplateFileForXRel = templateSourceLabel;
DynamicMap.TemplateSuffixForXRel = 'LowSpeedTemplateBank_20241106';
DynamicMap.XCenterBySensor = collect_template_xcenter_local(Template, P.analysisSensors);
DynamicMap.SelectionInfo = build_selection_info_local(P, RotatingDynamicMap);
DynamicMap.GlobalTimeWindow = globalTime;
DynamicMap.SelectedRawFileIDs = [];
DynamicMap.CorrectedGapLibrary = rmfield_if_present_local(CorrectedGapLibrary, 'responseSurface');
if ~isempty(RotatingDynamicMap)
    DynamicMap.RotatingReference = RotatingDynamicMap;
end
if ~isempty(DirectReference)
    DynamicMap.DirectReference = DirectReference;
end
DynamicMap.Window = Window;

timeTag = strrep(sprintf('T%07.3f', P.analysisStartTimeSec), '.', 'p');
matFile = fullfile(outDir, sprintf( ...
    'Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s_W%dS%d%s.mat', ...
    P.targetBlade, sensorTag, timeTag, P.windowLaps, P.slidingStepLaps, outputSuffix));
save(matFile, 'DynamicMap', 'Template', 'CorrectedGapLibrary', 'P', '-v7.3');

summary = build_window_summary_local(DynamicMap);
csvFile = fullfile(outDir, sprintf( ...
    'Step06_BuildGapAwareDynamicMap_WindowSummary_20241106_B%d_%s_%s_W%dS%d%s.csv', ...
    P.targetBlade, sensorTag, timeTag, P.windowLaps, P.slidingStepLaps, outputSuffix));
writetable(summary, csvFile);

figFile = fullfile(figDir, sprintf( ...
    'Step06_BuildGapAwareDynamicMap_20241106_B%d_%s_%s_W%dS%d%s.png', ...
    P.targetBlade, sensorTag, timeTag, P.windowLaps, P.slidingStepLaps, outputSuffix));
plot_dynamic_map_local(DynamicMap, Template, figFile);

fprintf('\nStep06 gap-aware DynamicMap complete.\n');
fprintf('Selected raw high-speed time: %.6f to %.6f s\n', globalTime(1), globalTime(2));
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, csvFile, figFile);
disp(summary);

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

function value = parse_float_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end

value = str2double(raw);
if ~isfinite(value) || value <= 0
    error('%s must be a positive number.', name);
end
end

function suffix = sanitize_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_' suffix];
end
end

function Template = clean_template_local(Template, sensorIds)
Template.Sensor = Template.Sensor(ismember([Template.Sensor.sensor_id], sensorIds));
[~, order] = ismember(sensorIds, [Template.Sensor.sensor_id]);
Template.Sensor = Template.Sensor(order);
for is = 1:numel(Template.Sensor)
    Tpl = Template.Sensor(is);
    x = Tpl.x_grid(:);
    v = Tpl.v_grid(:);
    finite = isfinite(x) & isfinite(v);
    x = x(finite);
    v = v(finite);
    [x, ord] = sort(x);
    v = v(ord);
    [x, keep] = unique(x, 'stable');
    v = v(keep);
    if numel(x) < 20
        error('Template P%d has too few finite points.', Tpl.sensor_id);
    end
    Template.Sensor(is).x_grid = x;
    Template.Sensor(is).v_grid = v;
    Template.Sensor(is).dv_dx = gradient(v, x);
    Template.Sensor(is).x_domain = [min(x), max(x)];
    if ~isfield(Template.Sensor(is), 'xc') || ~isfinite(Template.Sensor(is).xc)
        Template.Sensor(is).xc = 0;
    end
end
Template.analysisSensors = sensorIds(:).';
end

function opr = load_foundation_opr_local(oprFile)
if ~isfile(oprFile)
    error('Missing Foundation OPR file: %s', oprFile);
end
S = load(oprFile, 'jiluOPR');
raw = S.jiluOPR;
if isempty(raw) || size(raw, 2) < 1
    error('Invalid Foundation jiluOPR in %s.', oprFile);
end
t = raw(:, 1);
t = t(isfinite(t) & t > 0);
opr.tCenter = t(:);
opr.speedTime = t(1:end-1);
opr.rpm = 60 ./ diff(t);
end

function f = local_rot_freq_local(oprTime, t)
idx = find(oprTime <= t, 1, 'last');
if isempty(idx) || idx >= numel(oprTime)
    f = NaN;
else
    f = 1 / (oprTime(idx + 1) - oprTime(idx));
end
end

function theta = local_rot_phase_local(oprTime, t)
if numel(oprTime) < 2
    theta = NaN(size(t));
    return;
end
phaseAnchor = 2 * pi * (0:numel(oprTime)-1).';
theta = interp1(oprTime(:), phaseAnchor, t(:), 'linear', NaN);
theta = reshape(theta, size(t));
end

function [tRaw, vRaw, ok] = read_raw_time_window_local(caseDir, sensorId, ...
    timeStartSec, timeEndSec, sampleRateHz, blockStepSamples, blockLabelStep, ...
    rawCache, rawCacheMaxEntries)
sampleStart = floor(timeStartSec * sampleRateHz);
sampleEnd = ceil(timeEndSec * sampleRateHz);
blockStart = sample_to_block_label_local(sampleStart, blockStepSamples, blockLabelStep);
blockEnd = sample_to_block_label_local(sampleEnd, blockStepSamples, blockLabelStep);

rawAll = zeros(0, 2);
for blockLabel = blockStart:blockLabelStep:blockEnd
    raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries);
    rawAll = [rawAll; raw]; %#ok<AGROW>
end
if isempty(rawAll)
    tRaw = [];
    vRaw = [];
    ok = false;
    return;
end

rawAll(rawAll(:, 1) == 0, :) = [];
rawAll = sortrows(rawAll, 1);
rawT = rawAll(:, 1) / sampleRateHz;
keep = rawT >= timeStartSec & rawT <= timeEndSec & ...
    isfinite(rawAll(:, 2));
tRaw = rawT(keep);
vRaw = rawAll(keep, 2);
ok = numel(tRaw) >= 10;
end

function xc = collect_template_xcenter_local(Template, sensorIds)
xc = NaN(numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
    if isempty(idx) || ~isfield(Template.Sensor(idx), 'xc')
        continue;
    end
    xc(i) = Template.Sensor(idx).xc;
end
end

function info = build_selection_info_local(P, RotatingDynamicMap)
info = struct();
info.region_mode = 'manual_start';
info.start_time_sec = P.analysisStartTimeSec;
info.target_laps = P.targetLaps;
info.window_laps = P.windowLaps;
info.sliding_step_laps = P.slidingStepLaps;
info.analysis_sensors = P.analysisSensors(:).';
if ~isempty(RotatingDynamicMap) && isfield(RotatingDynamicMap, 'SelectionInfo')
    info.rotating_selection_info = RotatingDynamicMap.SelectionInfo;
else
    info.rotating_selection_info = struct();
end
end

function S = rmfield_if_present_local(S, fieldName)
if isstruct(S) && isfield(S, fieldName)
    S = rmfield(S, fieldName);
end
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
if isfield(S, varName)
    raw = S.(varName);
else
    raw = zeros(0, 2);
end
rawCache(key) = raw;
if rawCache.Count > rawCacheMaxEntries
    cacheKeys = keys(rawCache);
    remove(rawCache, cacheKeys(1:max(1, numel(cacheKeys) - rawCacheMaxEntries)));
end
end

function w = build_weight_local(t, v)
if numel(v) < 3
    w = ones(size(v));
    return;
end
vSmooth = movmedian(v(:), min(11, numel(v)));
dv = abs(gradient(vSmooth, t(:)));
if max(dv) > 0
    w = dv ./ max(dv);
else
    w = ones(size(v));
end
w = max(0.05, w);
end

function w = normalize_weight_local(w)
if isempty(w)
    return;
end
w = max(0.05, w(:));
if max(w) > 0
    w = w ./ max(w);
end
w = max(0.05, w);
end

function T = build_window_summary_local(DynamicMap)
rows = {};
for iw = 1:numel(DynamicMap.Window)
    W = DynamicMap.Window(iw);
    for is = 1:numel(W.Sensor)
        D = W.Sensor(is);
        rows{end + 1, 1} = table(W.window_id, string(mat2str(W.lap_range)), ...
            D.sensor_id, W.time_window(1), W.time_window(2), ...
            W.event_time_window(1), W.event_time_window(2), D.point_count, ...
            min(D.x_rel, [], 'omitnan'), max(D.x_rel, [], 'omitnan'), ...
            max(D.V, [], 'omitnan'), W.rot_freq_mean_hz, ...
            'VariableNames', {'windowId','lapRange','sensorId','timeStartSec', ...
            'timeEndSec','eventStartSec','eventEndSec','pointCount', ...
            'xMinMm','xMaxMm','maxVoltageMv','rotFreqMeanHz'});
    end
end
T = vertcat(rows{:});
end

function plot_dynamic_map_local(DynamicMap, Template, figFile)
fig = figure('Name', '20241106 Step06K direct DynamicMap', 'Color', 'w', ...
    'Position', [80, 60, 1450, 860], 'NumberTitle', 'off');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for is = 1:numel(DynamicMap.analysisSensors)
    sid = DynamicMap.analysisSensors(is);
    eventTimes = DynamicMap.selectionEventTimesBySensor{is};
    plot(eventTimes, is * ones(size(eventTimes)), 'o-', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('P%d', sid));
end
xlabel('Event time (s)');
ylabel('Sensor row');
title(sprintf('Selected blade-pass events, %.3f-%.3f s', ...
    DynamicMap.globalTimeWindow(1), DynamicMap.globalTimeWindow(2)));
legend('Location', 'best');

plotWindow = DynamicMap.Window(1);
colors = lines(numel(plotWindow.Sensor));
nexttile; hold on; grid on; box on;
for is = 1:numel(plotWindow.Sensor)
    D = plotWindow.Sensor(is);
    scatter(D.t, D.V, 5, colors(is, :), 'filled', ...
        'MarkerFaceAlpha', 0.18, 'DisplayName', sprintf('P%d', D.sensor_id));
end
xlabel('Time (s)');
ylabel('Voltage (mV)');
title(sprintf('Window %d raw waveform samples', plotWindow.window_id));
legend('Location', 'best');

nexttile; hold on; grid on; box on;
for is = 1:numel(Template.Sensor)
    Tpl = Template.Sensor(is);
    plot(Tpl.x_grid, Tpl.v_grid, '-', 'LineWidth', 1.3, ...
        'DisplayName', sprintf('P%d template', Tpl.sensor_id));
end
for is = 1:numel(plotWindow.Sensor)
    D = plotWindow.Sensor(is);
    scatter(D.x_rel, D.V, 5, colors(is, :), 'filled', ...
        'MarkerFaceAlpha', 0.12, 'DisplayName', sprintf('P%d dynamic', D.sensor_id));
end
xlabel('x relative to template center (mm)');
ylabel('Voltage (mV)');
title(sprintf('Window %d continuous waveform in x', plotWindow.window_id));
legend('Location', 'bestoutside');

for is = 1:min(3, numel(plotWindow.Sensor))
    D = plotWindow.Sensor(is);
    nexttile; hold on; grid on; box on;
    events = unique(D.event_time(isfinite(D.event_time)));
    for ie = 1:numel(events)
        mask = abs(D.event_time - events(ie)) < 1e-10;
        plot(D.x_rel(mask), D.V(mask), '-', 'LineWidth', 0.8);
    end
    xlabel('x relative to template center (mm)');
    ylabel('Voltage (mV)');
    title(sprintf('P%d window %d, laps %s', D.sensor_id, ...
        plotWindow.window_id, mat2str(plotWindow.lap_range)));
end

exportgraphics(fig, figFile, 'Resolution', 300);
end

function assert_capacitive_gap_sensors_local(sensorIds, overrideEnv)
allowed = [5 7];
if all(ismember(sensorIds(:).', allowed))
    return;
end
if parse_logical_env_local(overrideEnv, false)
    warning('Using non-capacitive sensors %s in the gap-aware DynamicMap because %s is enabled.', ...
        mat2str(sensorIds), overrideEnv);
    return;
end
error(['The 20241106 gap-aware route is calibrated for capacitive sensors [5 7]. ' ...
    'Requested sensors %s include non-capacitive probes. Set %s=1 only for diagnostics.'], ...
    mat2str(sensorIds), overrideEnv);
end

function value = parse_logical_env_local(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    value = logical(defaultValue);
    return;
end
if ismember(raw, {'1','true','yes','on'})
    value = true;
elseif ismember(raw, {'0','false','no','off'})
    value = false;
else
    error('%s must be one of on/off, true/false, yes/no, or 1/0.', name);
end
end

function [Template, label] = load_template_for_run_local(P)
Template = [];
label = '';
if ~isfile(P.templateBankFile)
    return;
end
S = load(P.templateBankFile, 'LowSpeedTemplateBank');
bank = S.LowSpeedTemplateBank;
Template = struct();
Template.dataset = bank.dataset;
Template.targetBlade = P.targetBlade;
Template.analysisSensors = P.analysisSensors(:).';
Template.SensorIDs = P.analysisSensors(:).';
Template.SensorTag = ['S', sprintf('%d', P.analysisSensors)];
Template.Sensor = struct([]);
for sid = P.analysisSensors(:).'
    idx = find([bank.entry.bladeId] == P.targetBlade & [bank.entry.sensorId] == sid, 1, 'first');
    if isempty(idx)
        error('LowSpeedTemplateBank does not contain B%d CH%d.', P.targetBlade, sid);
    end
    Template.Sensor = append_sensor_struct_local(Template.Sensor, bank.entry(idx).sensorTemplate);
end
label = P.templateBankFile;
end

function [CorrectedGapLibrary, label] = load_gap_library_for_run_local(P, Template)
CorrectedGapLibrary = [];
label = '';
if ~isfile(P.gapBankFile)
    return;
end
S = load(P.gapBankFile, 'GapCalibrationBank');
bank = S.GapCalibrationBank;
sensorCorr = struct([]);
for is = 1:numel(P.gapSensors)
    sid = P.gapSensors(is);
    idx = find([bank.entry.bladeId] == P.targetBlade & [bank.entry.sensorId] == sid, 1, 'first');
    if isempty(idx)
        error('GapCalibrationBank does not contain B%d CH%d.', P.targetBlade, sid);
    end
    sensorCorr = append_sensor_struct_local(sensorCorr, bank.entry(idx).sensor);
end
gapTemplate = clean_template_local(Template, P.gapSensors);
CorrectedGapLibrary = struct();
CorrectedGapLibrary.dataset = bank.dataset;
CorrectedGapLibrary.method = 'gap_bank_view_for_dynamic_map';
CorrectedGapLibrary.description = 'Per-run view assembled from all-blade GapCalibrationBank.';
CorrectedGapLibrary.responseFile = bank.responseFile;
CorrectedGapLibrary.templateBankFile = bank.templateBankFile;
CorrectedGapLibrary.templateFile = string(bank.templateBankFile);
CorrectedGapLibrary.templateSourceMode = 'low_speed_template_bank';
CorrectedGapLibrary.templateSourceFile = string(bank.templateBankFile);
CorrectedGapLibrary.targetBlade = P.targetBlade;
CorrectedGapLibrary.analysisSensors = P.gapSensors(:).';
CorrectedGapLibrary.gapSensors = P.gapSensors(:).';
CorrectedGapLibrary.directOnlySensors = setdiff(P.analysisSensors, P.gapSensors, 'stable');
CorrectedGapLibrary.responseSurface = bank.responseSurface;
CorrectedGapLibrary.lowSpeedTemplate = gapTemplate;
CorrectedGapLibrary.sensor = sensorCorr;
CorrectedGapLibrary.formula = 'T_low(x)+a*(F_raw(g+mu*(x-tau),k*(x-tau))-F_raw(g0+mu*(x-tau),k*(x-tau)))';
CorrectedGapLibrary.cfg = bank.cfg;
label = P.gapBankFile;
end

function Template = merge_gap_template_sensors_local(Template, gapTemplate, gapSensors)
for sid = gapSensors(:).'
    idxDst = find([Template.Sensor.sensor_id] == sid, 1, 'first');
    idxSrc = find([gapTemplate.Sensor.sensor_id] == sid, 1, 'first');
    if isempty(idxDst) || isempty(idxSrc)
        continue;
    end
    [Template.Sensor, gapSensor] = align_struct_array_fields_local(Template.Sensor, gapTemplate.Sensor(idxSrc));
    Template.Sensor(idxDst) = gapSensor;
end
end

function out = append_sensor_struct_local(out, sensor)
if isempty(out)
    out = sensor;
    return;
end
[out, sensor] = align_struct_array_fields_local(out, sensor);
out(end + 1) = sensor; %#ok<AGROW>
end

function [arrayOut, itemOut] = align_struct_array_fields_local(arrayOut, itemOut)
fieldsAll = unique([fieldnames(arrayOut); fieldnames(itemOut)], 'stable');
for i = 1:numel(fieldsAll)
    name = fieldsAll{i};
    if ~isfield(arrayOut, name)
        [arrayOut.(name)] = deal([]);
    end
    if ~isfield(itemOut, name)
        itemOut.(name) = [];
    end
end
arrayOut = orderfields(arrayOut, fieldsAll);
itemOut = orderfields(itemOut, fieldsAll);
end

function numWindows = planned_window_count_local(targetLaps, windowLaps, slidingStepLaps)
if ~isfinite(targetLaps) || ~isfinite(windowLaps) || ~isfinite(slidingStepLaps) || ...
        targetLaps <= 0 || windowLaps <= 0 || slidingStepLaps <= 0
    error('target/window/step laps must be positive finite values.');
end
if targetLaps < windowLaps
    error('target laps (%d) must be >= window laps (%d).', targetLaps, windowLaps);
end
numWindows = floor((targetLaps - windowLaps) / slidingStepLaps) + 1;
end
