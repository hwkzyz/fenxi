%% Step06: build dynamic highMap from 20241106 raw probe waveforms
% This is the first waveform-domain dynamic data structure after Step05I.
%
% Main output remains the gap-prior highMap structure. What can change now
% is only the event reference source:
%   1) legacy Step02 blade-slot timing
%   2) rotating-calibration Step03/Step04 region + numbering reference
%
% Output highMap contains:
%   t_v              raw sample time of waveform points
%   x_v              pulse-local circumferential coordinate, mm
%   V_a              edge-baseline-corrected voltage, mV
%   S_v              probe channel id
%   bladeSlot_v      blade slot relative to the 20241106 OPR pulse
%   regionId_v       resonance region id from Step04/Step05 diagnostic table
%   theta_v          continuous rotor phase, rad, from the OPR pulse train
%
% Figures are opened directly. MAT/CSV tables are saved for later steps.

clear; clc; close all;

%% 1. Paths and simple settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

step02File = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
step05File = fullfile(outDir, 'Step05_Response_Surface_20241106.mat');
step05IFile = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_20241106.mat');
regionCsvFile = fullfile(outDir, 'Step05_ResonanceRegions_20241106.csv');
rotDir = fullfile(rootDir, '20241106_low_speed_rotating_calibration');

if ~isfile(step02File)
    error('Run Step02 first. Missing file: %s', step02File);
end
if ~isfile(step05File)
    error('Run Step05 first. Missing file: %s', step05File);
end

S02 = load(step02File);
S05 = load(step05File, 'responseSurface');
responseSurface = S05.responseSurface;

if isfile(step05IFile)
    S05I = load(step05IFile, 'OffsetTiltResponseSurface');
    offsetTiltSurface = S05I.OffsetTiltResponseSurface; %#ok<NASGU>
else
    offsetTiltSurface = []; %#ok<NASGU>
end

eventSourceMode = 'gap_prior_legacy_step02'; % gap_prior_legacy_step02 | rotating_step03_step04
highDir = S02.highDir;
sensorIds = S02.sensorIds(:).';
bladeCount = S02.bladeCount;
sampleRateHz = S02.sampleRateHz;
rTipMm = S02.rTipMm;
highDisp = S02.highDisp;

if isfile(regionCsvFile)
    regionTable = readtable(regionCsvFile);
else
    warning('No resonance region table found. Step06 will use the whole high-speed time range.');
    regionTable = table(1, 0, inf, 'VariableNames', {'regionId','bttStartSec','bttEndSec'});
end

% Debug-first defaults:
%   R2 is the strongest 20241106 resonance region found in Step04.
%   Only the nearly steady resonance near 75 s is extracted first.
%   Probe 5/7 are capacitive probes and are used for gap-library waveform
%   decoupling. Probe 2/3 are eddy-current probes; keep them for direct
%   low-speed calibration and vibration-parameter identification, not this
%   shared gap-library highMap.
%   Blade slot 1 is only a slot label relative to the 20241106 OPR pulse.
% Use the STEP06_* environment variables below to expand this after checking.
analysisSensors = intersect([5 7], sensorIds, 'stable');
bladeSlots = 1;
selectedRegionIds = 2;
if ~any(regionTable.regionId == selectedRegionIds)
    selectedRegionIds = regionTable.regionId(1);
end
manualTimeWindowSec = [72 78];  % around the steady 75 s resonance inside R2.

% Keep this modest first. Increase after the workflow is verified.
pulsesPerRegionSensorBlade = 20;
pulseHalfWindowSec = 2.5e-4;
pointsPerPulse = 101;
edgeFractionForBaseline = 0.18;
minFinitePointsPerPulse = 30;
blockStepSamples = 1e7;
blockLabelStep = 1000;
rawCacheMaxEntries = 10;

analysisSensors = parse_int_env_local('STEP06_ANALYSIS_SENSORS', analysisSensors);
bladeSlots = parse_int_env_local('STEP06_BLADE_SLOTS', bladeSlots);
selectedRegionIds = parse_int_env_local('STEP06_REGION_IDS', selectedRegionIds);
manualTimeWindowSec = parse_range_env_local('STEP06_TIME_WINDOW_SEC', manualTimeWindowSec);
pulsesPerRegionSensorBlade = parse_scalar_env_local('STEP06_PULSES_PER_REGION_SLOT', pulsesPerRegionSensorBlade);
pulseHalfWindowSec = parse_scalar_env_local('STEP06_PULSE_HALF_WINDOW_SEC', pulseHalfWindowSec);
pointsPerPulse = round(parse_scalar_env_local('STEP06_POINTS_PER_PULSE', pointsPerPulse));
eventSourceMode = parse_mode_env_local('STEP06_EVENT_SOURCE_MODE', eventSourceMode, ...
    {'gap_prior_legacy_step02', 'rotating_step03_step04'});

xGrid = responseSurface.xGrid(:);
if isfield(responseSurface, 'effectiveWindow') && any(responseSurface.effectiveWindow)
    xMin = min(xGrid(responseSurface.effectiveWindow));
    xMax = max(xGrid(responseSurface.effectiveWindow));
else
    xMin = min(xGrid);
    xMax = max(xGrid);
end

%% 2. Load OPR speed reference
opr = load_opr_local(highDir, sampleRateHz);

fprintf('\n=== Step06 20241106 dynamic highMap ===\n');
fprintf('High-speed folder: %s\n', highDir);
fprintf('Event source mode: %s\n', eventSourceMode);
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Blade slots: %s\n', mat2str(bladeSlots));
fprintf('Region ids: %s\n', mat2str(selectedRegionIds));
fprintf('Manual time window: %s s\n', mat2str(manualTimeWindowSec));
fprintf('Pulses per region/sensor/blade: %d\n', pulsesPerRegionSensorBlade);
fprintf('Pulse half window: %.1f us, points per pulse: %d\n', ...
    pulseHalfWindowSec * 1e6, pointsPerPulse);

%% 3. Select blade-pass events and cut raw waveform pulses
rawCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

rotEventRef = struct();
if strcmpi(eventSourceMode, 'rotating_step03_step04')
    rotEventRef = load_rotating_event_reference_local(rotDir, analysisSensors, bladeSlots, sampleRateHz, highDir);
    if isfield(rotEventRef, 'timeWindowSec') && all(isfinite(rotEventRef.timeWindowSec))
        manualTimeWindowSec = intersect_time_ranges_local(manualTimeWindowSec, rotEventRef.timeWindowSec);
    end
    if isfield(rotEventRef, 'regionIds') && ~isempty(rotEventRef.regionIds)
        selectedRegionIds = rotEventRef.regionIds(:).';
    end
end

t_v = [];
x_v = [];
V_raw_mV = [];
V_a = [];
baseline_mV_v = [];
S_v = [];
bladeSlot_v = [];
regionId_v = [];
eventTime_v = [];
bttDisplacementMm_v = [];
rotFreqHz_v = [];
theta_v = [];
W_v = [];

windowRows = {};
baselineRows = {};
examplePulse = struct('sensorId', {}, 'bladeSlot', {}, 'regionId', {}, ...
    'eventTime', {}, 'xMm', {}, 'vCorrMv', {});

for ir = 1:numel(selectedRegionIds)
    rid = selectedRegionIds(ir);
    [tStart, tEnd, regionRowsOk] = resolve_region_time_window_local( ...
        eventSourceMode, rid, regionTable, regionCsvFile, manualTimeWindowSec, rotEventRef);
    if ~regionRowsOk
        continue;
    end

    for sid = analysisSensors
        for bladeSlot = bladeSlots
            [eventTimesAll, xBttAll, eventMeta] = load_event_reference_local( ...
                eventSourceMode, sid, bladeSlot, tStart, tEnd, highDisp, rotEventRef);
            eventTimes = eventTimesAll(:);
            xBtt = xBttAll(:);
            if isempty(eventTimes)
                continue;
            end

            pick = pick_evenly_local(numel(eventTimes), pulsesPerRegionSensorBlade);
            eventTimes = eventTimes(pick);
            xBtt = xBtt(pick);

            pulsePeakMv = NaN(numel(eventTimes), 1);
            pulseValidPoints = zeros(numel(eventTimes), 1);

            for ie = 1:numel(eventTimes)
                tEvent = eventTimes(ie);
                rotFreqHz = local_rot_freq_local(opr.tCenter, tEvent);
                tipSpeedMmS = 2 * pi * rTipMm * rotFreqHz;

                [tq, vq, rawOk] = read_raw_pulse_local(highDir, sid, tEvent, ...
                    pulseHalfWindowSec, pointsPerPulse, sampleRateHz, ...
                    blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries);
                if ~rawOk
                    continue;
                end

                xPulse = (tq - tEvent) * tipSpeedMmS;
                thetaPulse = local_rot_phase_local(opr.tCenter, tq);
                edgeCount = max(5, round(edgeFractionForBaseline * numel(vq)));
                edgeValues = [vq(1:edgeCount); vq(end-edgeCount+1:end)];
                baselineV = median(edgeValues, 'omitnan');
                vRawMv = 1000 * vq;
                vCorrMv = 1000 * (vq - baselineV);

                keepPoint = isfinite(xPulse) & isfinite(vCorrMv) & ...
                    xPulse >= xMin & xPulse <= xMax;
                if nnz(keepPoint) < minFinitePointsPerPulse
                    continue;
                end

                w = build_weight_local(tq(keepPoint), vCorrMv(keepPoint));
                n = nnz(keepPoint);

                t_v = [t_v; tq(keepPoint)]; %#ok<AGROW>
                x_v = [x_v; xPulse(keepPoint)]; %#ok<AGROW>
                V_raw_mV = [V_raw_mV; vRawMv(keepPoint)]; %#ok<AGROW>
                V_a = [V_a; vCorrMv(keepPoint)]; %#ok<AGROW>
                baseline_mV_v = [baseline_mV_v; 1000 * baselineV * ones(n, 1)]; %#ok<AGROW>
                S_v = [S_v; sid * ones(n, 1)]; %#ok<AGROW>
                bladeSlot_v = [bladeSlot_v; eventMeta.bladeLabel * ones(n, 1)]; %#ok<AGROW>
                regionId_v = [regionId_v; eventMeta.regionId * ones(n, 1)]; %#ok<AGROW>
                eventTime_v = [eventTime_v; tEvent * ones(n, 1)]; %#ok<AGROW>
                bttDisplacementMm_v = [bttDisplacementMm_v; xBtt(ie) * ones(n, 1)]; %#ok<AGROW>
                rotFreqHz_v = [rotFreqHz_v; rotFreqHz * ones(n, 1)]; %#ok<AGROW>
                theta_v = [theta_v; thetaPulse(keepPoint)]; %#ok<AGROW>
                W_v = [W_v; w(:)]; %#ok<AGROW>

                pulsePeakMv(ie) = max(vCorrMv(keepPoint), [], 'omitnan');
                pulseValidPoints(ie) = n;

                if numel(examplePulse) < 12
                    examplePulse(end + 1).sensorId = sid; %#ok<AGROW>
                    examplePulse(end).bladeSlot = eventMeta.bladeLabel;
                    examplePulse(end).regionId = eventMeta.regionId;
                    examplePulse(end).eventTime = tEvent;
                    examplePulse(end).xMm = xPulse(keepPoint);
                    examplePulse(end).vCorrMv = vCorrMv(keepPoint);
                end
            end

            goodPulse = isfinite(pulsePeakMv) & pulseValidPoints > 0;
            if any(goodPulse)
                windowRows{end + 1, 1} = table(eventMeta.regionId, sid, eventMeta.bladeLabel, tStart, tEnd, ...
                    nnz(goodPulse), median(pulsePeakMv(goodPulse), 'omitnan'), ...
                    max(pulsePeakMv(goodPulse), [], 'omitnan'), ...
                    median(pulseValidPoints(goodPulse), 'omitnan'), ...
                    'VariableNames', {'regionId','sensorId','bladeSlot', ...
                    'regionStartSec','regionEndSec','pulseCount', ...
                    'medianPeakMv','maxPeakMv','medianPointCount'});
            end
        end
    end
end

if isempty(t_v)
    error('No valid highMap samples were constructed. Try increasing pulseHalfWindowSec or check raw channel files.');
end

W_v = W_v ./ max(W_v);
W_v = max(0.05, W_v);

%% 4. Build highMap structure
if isempty(windowRows)
    windowTable = table();
else
    windowTable = vertcat(windowRows{:});
end

for sid = analysisSensors
    mask = S_v == sid;
    if any(mask)
        baselineRows{end + 1, 1} = table(sid, median(baseline_mV_v(mask), 'omitnan'), ...
            min(baseline_mV_v(mask), [], 'omitnan'), max(baseline_mV_v(mask), [], 'omitnan'), ...
            nnz(mask), 'VariableNames', {'sensorId','medianBaselineMv', ...
            'minBaselineMv','maxBaselineMv','pointCount'}); %#ok<SAGROW>
    end
end
if isempty(baselineRows)
    baselineTable = table();
else
    baselineTable = vertcat(baselineRows{:});
end

highMap = struct();
highMap.dataset = '20241106';
highMap.method = 'raw_probe_pulse_highmap_from_step02_timing';
highMap.highDir = highDir;
highMap.analysisSensors = analysisSensors;
highMap.bladeSlots = bladeSlots;
highMap.selectedRegionIds = selectedRegionIds;
highMap.manualTimeWindowSec = manualTimeWindowSec;
highMap.eventSourceMode = eventSourceMode;
highMap.regionTable = regionTable;
highMap.t_v = t_v;
highMap.x_v = x_v;
highMap.x_fit_v = x_v;
highMap.V_raw_mV = V_raw_mV;
highMap.V_a = V_a;
highMap.baseline_mV_v = baseline_mV_v;
highMap.S_v = S_v;
highMap.bladeSlot_v = bladeSlot_v;
highMap.regionId_v = regionId_v;
highMap.eventTime_v = eventTime_v;
highMap.bttDisplacementMm_v = bttDisplacementMm_v;
highMap.rotFreqHz_v = rotFreqHz_v;
highMap.theta_v = theta_v;
highMap.W_v = W_v;
highMap.pointCount = numel(t_v);
highMap.sampleRateHz = sampleRateHz;
highMap.rTipMm = rTipMm;
highMap.pulseHalfWindowSec = pulseHalfWindowSec;
highMap.pointsPerPulse = pointsPerPulse;
highMap.responseSurfaceXMinMm = xMin;
highMap.responseSurfaceXMaxMm = xMax;
highMap.windowTable = windowTable;
highMap.dynamicBaselineBySensor = baselineTable;
highMap.examplePulse = examplePulse;
highMap.coordinateNote = ['x_v is pulse-local x=(t-t-bladepass)*2*pi*r_tip*f_rot. ' ...
    'No old-date low-speed template or sensor geometry is used. ' ...
    'bttDisplacementMm_v stores legacy Step02 displacement when available; rotating Step04 event mode leaves it as NaN context.'];
highMap.voltageNote = 'V_a is mV after subtracting each pulse edge median baseline.';

%% 5. Save reusable outputs
sensorTag = ['S', sprintf('%d', analysisSensors)];
bladeTag = ['B', sprintf('%d', bladeSlots)];
matFile = fullfile(outDir, 'Step06_HighMap_20241106.mat');
csvFile = fullfile(outDir, 'Step06_HighMap_Window_Summary_20241106.csv');
baselineCsvFile = fullfile(outDir, 'Step06_HighMap_Baseline_Summary_20241106.csv');
taggedMatFile = fullfile(outDir, sprintf('Step06_HighMap_20241106_%s_%s.mat', bladeTag, sensorTag));
taggedCsvFile = fullfile(outDir, sprintf('Step06_HighMap_Window_Summary_20241106_%s_%s.csv', bladeTag, sensorTag));
save(matFile, 'highMap', '-v7.3');
save(taggedMatFile, 'highMap', '-v7.3');
writetable(windowTable, csvFile);
writetable(windowTable, taggedCsvFile);
writetable(baselineTable, baselineCsvFile);

%% 6. Direct visualization
plot_highmap_dashboard_local(highMap, responseSurface, opr);

fprintf('\nStep06 complete.\n');
fprintf('highMap saved to:\n  %s\n', matFile);
fprintf('Tagged highMap copy saved to:\n  %s\n', taggedMatFile);
fprintf('Window summary saved to:\n  %s\n', csvFile);
fprintf('Baseline summary saved to:\n  %s\n', baselineCsvFile);
fprintf('Constructed points: %d\n', highMap.pointCount);

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

function mode = parse_mode_env_local(name, defaultMode, allowedModes)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    mode = defaultMode;
    return;
end
if ~ismember(raw, allowedModes)
    error('%s must be one of: %s', name, strjoin(allowedModes, ', '));
end
mode = raw;
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

function eventRef = load_rotating_event_reference_local(rotDir, analysisSensors, bladeSlots, sampleRateHz, highDir)
eventRef = struct();
if exist(rotDir, 'dir') ~= 7
    error('Rotating calibration folder not found: %s', rotDir);
end
addpath(rotDir);
cleanupPath = onCleanup(@() rmpath(rotDir));
P = NewFlow_Config_20241106();
if ~strcmpi(string(P.data.highSpeedDir), string(highDir))
    warning('Rotating high-speed dir differs from current Step06 highDir. Using rotating event references anyway.');
end

require_rotating_file_local(P.files.regionSelection, 'rotating Step03 region selection');
require_rotating_file_local(P.files.highSpeedNumbering, 'rotating Step04 high-speed numbering');

S03 = load(P.files.regionSelection, 'RegionSelection');
S04 = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
RegionSelection = S03.RegionSelection;
HighSpeedNumbering = S04.HighSpeedNumbering;

eventRef.regionSelection = RegionSelection;
eventRef.highSpeedNumbering = HighSpeedNumbering;
eventRef.regionIds = collect_rotating_region_ids_local(RegionSelection);
eventRef.timeWindowSec = collect_rotating_time_window_local(HighSpeedNumbering, RegionSelection, sampleRateHz);
eventRef.analysisSensors = analysisSensors(:).';
eventRef.targetBlades = bladeSlots(:).';
end

function require_rotating_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function regionIds = collect_rotating_region_ids_local(RegionSelection)
if isfield(RegionSelection, 'region_id') && isfinite(RegionSelection.region_id)
    regionIds = RegionSelection.region_id;
else
    regionIds = 0;
end
end

function timeWindowSec = collect_rotating_time_window_local(HighSpeedNumbering, RegionSelection, sampleRateHz)
timeWindowSec = [NaN NaN];
if isfield(HighSpeedNumbering, 'time_window') && numel(HighSpeedNumbering.time_window) >= 2 && ...
        all(isfinite(HighSpeedNumbering.time_window))
    timeWindowSec = HighSpeedNumbering.time_window(1:2);
elseif isfield(RegionSelection, 'region_start_sec') && isfield(RegionSelection, 'end_time_sec') && ...
        isfinite(RegionSelection.region_start_sec) && isfinite(RegionSelection.end_time_sec)
    timeWindowSec = [RegionSelection.region_start_sec, RegionSelection.end_time_sec];
elseif isfield(HighSpeedNumbering, 'selected_revolution_ids') && ~isempty(HighSpeedNumbering.selected_revolution_ids)
    revIds = HighSpeedNumbering.selected_revolution_ids(:);
    timeWindowSec = [revIds(1), revIds(end)] ./ sampleRateHz;
end
end

function outRange = intersect_time_ranges_local(rangeA, rangeB)
if isempty(rangeA)
    outRange = rangeB;
    return;
end
if isempty(rangeB) || any(~isfinite(rangeB))
    outRange = rangeA;
    return;
end
outRange = [max(rangeA(1), rangeB(1)), min(rangeA(2), rangeB(2))];
end

function [tStart, tEnd, ok] = resolve_region_time_window_local(eventSourceMode, rid, regionTable, regionCsvFile, manualTimeWindowSec, rotEventRef)
ok = true;
switch lower(strtrim(eventSourceMode))
    case 'rotating_step03_step04'
        if isfield(rotEventRef, 'timeWindowSec') && all(isfinite(rotEventRef.timeWindowSec))
            tStart = rotEventRef.timeWindowSec(1);
            tEnd = rotEventRef.timeWindowSec(2);
        else
            warning('Rotating event reference has no valid time window. Skipping region %d.', rid);
            tStart = NaN;
            tEnd = NaN;
            ok = false;
            return;
        end
    otherwise
        rrow = regionTable(regionTable.regionId == rid, :);
        if isempty(rrow)
            warning('Region %d is not in %s. Skipped.', rid, regionCsvFile);
            tStart = NaN;
            tEnd = NaN;
            ok = false;
            return;
        end
        tStart = rrow.bttStartSec(1);
        tEnd = rrow.bttEndSec(1);
end

if ~isempty(manualTimeWindowSec)
    tStart = max(tStart, manualTimeWindowSec(1));
    tEnd = min(tEnd, manualTimeWindowSec(2));
end
if ~(isfinite(tStart) && isfinite(tEnd) && tEnd > tStart)
    warning('Resolved region %d does not overlap active time window %s. Skipped.', rid, mat2str(manualTimeWindowSec));
    ok = false;
end
end

function [eventTimes, xBtt, eventMeta] = load_event_reference_local(eventSourceMode, sid, bladeSlot, tStart, tEnd, highDisp, rotEventRef)
eventMeta = struct('bladeLabel', bladeSlot, 'regionId', 0);
if strcmpi(strtrim(eventSourceMode), 'rotating_step03_step04')
    [eventTimes, xBtt, eventMeta] = load_rotating_events_local(sid, bladeSlot, tStart, tEnd, rotEventRef);
else
    sensorIndex = find([highDisp.sensorId] == sid, 1, 'first');
    if isempty(sensorIndex) || bladeSlot < 1 || bladeSlot > numel(highDisp(sensorIndex).slot)
        eventTimes = [];
        xBtt = [];
        return;
    end
    slotData = highDisp(sensorIndex).slot(bladeSlot);
    eventTimesAll = slotData.time(:);
    xBttAll = slotData.xMm(:);
    keepEvent = eventTimesAll >= tStart & eventTimesAll <= tEnd & isfinite(xBttAll);
    eventTimes = eventTimesAll(keepEvent);
    xBtt = xBttAll(keepEvent);
    eventMeta.bladeLabel = bladeSlot;
end
end

function [eventTimes, xBtt, eventMeta] = load_rotating_events_local(sid, bladeId, tStart, tEnd, rotEventRef)
eventTimes = [];
xBtt = [];
eventMeta = struct('bladeLabel', bladeId, 'regionId', first_region_id_local(rotEventRef));
HighSpeedNumbering = rotEventRef.highSpeedNumbering;
sensorIdx = find([HighSpeedNumbering.sensor.sensor_id] == sid, 1, 'first');
if isempty(sensorIdx)
    warning('Rotating Step04 numbering does not contain sensor %d.', sid);
    return;
end
S = HighSpeedNumbering.sensor(sensorIdx);
if ~isfield(S, 'selected_rows_by_physical') || isempty(S.selected_rows_by_physical)
    warning('Rotating Step04 sensor %d has no selected_rows_by_physical.', sid);
    return;
end
if bladeId < 1 || bladeId > size(S.selected_rows_by_physical, 2)
    warning('Requested blade %d is outside rotating selected_rows_by_physical width for sensor %d.', bladeId, sid);
    return;
end
peakTimes = [];
if isfield(S, 'revolution_peak_times') && ~isempty(S.revolution_peak_times)
    peakTimes = S.revolution_peak_times(:, bladeId);
elseif isfield(S, 'numbering_table') && ~isempty(S.numbering_table)
    colName = sprintf('B%d', bladeId);
    if any(strcmpi(S.numbering_table.Properties.VariableNames, colName))
        peakTimes = S.numbering_table.(colName);
    end
end
if isempty(peakTimes)
    return;
end
keep = isfinite(peakTimes) & peakTimes >= tStart & peakTimes <= tEnd;
eventTimes = peakTimes(keep);
xBtt = NaN(size(eventTimes));
end

function regionId = first_region_id_local(rotEventRef)
regionId = 0;
if isfield(rotEventRef, 'regionIds') && ~isempty(rotEventRef.regionIds)
    regionId = rotEventRef.regionIds(1);
end
end

function opr = load_opr_local(caseDir, sampleRateHz)
f = fullfile(caseDir, 'jiluOPR.mat');
if ~isfile(f)
    error('Missing OPR file: %s', f);
end
S = load(f, 'jiluOPR');
raw = S.jiluOPR;
opr.raw = raw;
opr.tCenter = mean(raw(:, 1:2), 2) / sampleRateHz;
opr.speedTime = opr.tCenter(1:end-1);
opr.rpm = 60 ./ diff(opr.tCenter);
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

function theta = local_rot_phase_local(oprTime, t)
if numel(oprTime) < 2
    theta = NaN(size(t));
    return;
end
phaseAnchor = 2 * pi * (0:numel(oprTime)-1).';
theta = interp1(oprTime(:), phaseAnchor, t(:), 'linear', NaN);
theta = reshape(theta, size(t));
end

function idx = pick_evenly_local(n, maxCount)
maxCount = max(1, round(maxCount));
if n <= maxCount
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxCount))).';
end
end

function [tq, vq, ok] = read_raw_pulse_local(caseDir, sensorId, tCenter, halfWindowSec, ...
    pointsPerPulse, sampleRateHz, blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries)
tq = linspace(tCenter - halfWindowSec, tCenter + halfWindowSec, pointsPerPulse).';
sampleStart = floor((tCenter - halfWindowSec) * sampleRateHz);
sampleEnd = ceil((tCenter + halfWindowSec) * sampleRateHz);

blockStart = sample_to_block_label_local(sampleStart, blockStepSamples, blockLabelStep);
blockEnd = sample_to_block_label_local(sampleEnd, blockStepSamples, blockLabelStep);
blockLabels = blockStart:blockLabelStep:blockEnd;

rawAll = zeros(0, 2);
for blockLabel = blockLabels
    raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries);
    if ~isempty(raw)
        rawAll = [rawAll; raw]; %#ok<AGROW>
    end
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
vSmooth = movmedian(v(:), min(9, numel(v)));
dv = abs(gradient(vSmooth, t(:)));
if max(dv) > 0
    w = dv ./ max(dv);
else
    w = ones(size(v));
end
w = max(0.05, w);
end

function plot_highmap_dashboard_local(highMap, responseSurface, opr)
figure('Name', '20241106 Step06 dynamic highMap', 'Color', 'w', ...
    'Position', [80, 60, 1400, 850], 'NumberTitle', 'off');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot(opr.speedTime, opr.rpm, 'k-', 'LineWidth', 1.0);
yl = [min(opr.rpm, [], 'omitnan'), max(opr.rpm, [], 'omitnan')];
if diff(yl) <= 0
    yl = yl + [-1, 1];
end
for i = 1:height(highMap.regionTable)
    if ~ismember(highMap.regionTable.regionId(i), highMap.selectedRegionIds)
        continue;
    end
    x1 = highMap.regionTable.bttStartSec(i);
    x2 = highMap.regionTable.bttEndSec(i);
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1.0 0.85 0.25], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    text(mean([x1 x2]), yl(2), sprintf('R%d', highMap.regionTable.regionId(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
uistack(findobj(gca, 'Type', 'line'), 'top');
xlabel('BTT time (s)');
ylabel('Speed (rpm)');
title('OPR speed and selected resonance regions');

nexttile; hold on; grid on; box on;
colors = lines(numel(highMap.analysisSensors));
for i = 1:numel(highMap.analysisSensors)
    sid = highMap.analysisSensors(i);
    mask = highMap.S_v == sid;
    scatter(highMap.x_v(mask), highMap.V_a(mask), 5, colors(i, :), 'filled', ...
        'MarkerFaceAlpha', 0.16, 'DisplayName', sprintf('P%d', sid));
end
xlabel('Pulse-local x (mm)');
ylabel('Voltage after edge baseline (mV)');
title('All selected dynamic waveform samples');
legend('Location', 'bestoutside');

nexttile; hold on; grid on; box on;
T = highMap.windowTable;
if ~isempty(T)
    scatter(T.regionId, T.maxPeakMv, 35, T.sensorId, 'filled');
    xlabel('Region id');
    ylabel('Max pulse amplitude (mV)');
    title('Large waveform responses by region and probe');
    cb = colorbar;
    cb.Label.String = 'Probe channel';
end

nexttile; hold on; grid on; box on;
for i = 1:min(numel(highMap.examplePulse), 12)
    E = highMap.examplePulse(i);
    plot(E.xMm, E.vCorrMv, 'LineWidth', 0.9, ...
        'DisplayName', sprintf('R%d P%d B%d', E.regionId, E.sensorId, E.bladeSlot));
end
xlabel('Pulse-local x (mm)');
ylabel('Voltage after edge baseline (mV)');
title('Example extracted raw pulses');
legend('Location', 'bestoutside');

nexttile; hold on; grid on; box on;
if isfield(responseSurface, 'Yfit')
    gapCount = numel(responseSurface.gTrainMm);
    cmap = turbo(gapCount);
    for ig = 1:gapCount
        plot(responseSurface.xGrid, responseSurface.Yfit(:, ig), '-', ...
            'Color', cmap(ig, :), 'LineWidth', 0.8);
    end
end
sampleIdx = pick_evenly_local(numel(highMap.x_v), min(6000, numel(highMap.x_v)));
scatter(highMap.x_v(sampleIdx), highMap.V_a(sampleIdx), 4, [0.05 0.05 0.05], ...
    'filled', 'MarkerFaceAlpha', 0.12);
xlabel('x (mm)');
ylabel('Voltage (mV)');
title('Dynamic samples over static Step05 response curves');

nexttile; hold on; grid on; box on;
if ~isempty(T)
    groupLabels = strcat("R", string(T.regionId), "-P", string(T.sensorId), "-B", string(T.bladeSlot));
    bar(categorical(groupLabels), T.pulseCount);
    ylabel('Pulse count');
    title('Accepted pulses per region/probe/blade');
    xtickangle(75);
end
end
