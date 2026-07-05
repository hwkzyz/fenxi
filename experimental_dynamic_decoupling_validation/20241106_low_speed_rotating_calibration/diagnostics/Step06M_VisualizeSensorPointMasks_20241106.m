%% Step06M_VisualizeSensorPointMasks_20241106
% Visualize raw/query-safe/dynamic-effective/final-valid point masks for
% CH2/CH3/CH5/CH7 in the same blade window. This helps judge whether the
% edge selection is too conservative before changing thresholds.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local parameter block. Edit here directly.
S06M = struct();
S06M.analysisSensors = [2 3 5 7];
S06M.startTimeSec = 75.0;
S06M.bladeId = 4;
S06M.windowIds = [1 4];          % suspicious windows first; set [] for all
S06M.windowLaps = 3;
S06M.slidingStepLaps = 1;
S06M.pulsePadSec = 0;
S06M.dynamicWindowMode = 'legacy_row_bounds';
S06M.windowRule = 't_start=jilublade(row,1)-pulsePadSec; t_end=jilublade(row,2)+pulsePadSec';
S06M.pulseSelectionMode = 'single';
S06M.domainSelectionMode = 'hard';
S06M.domainMarginMM = 0.02;
S06M.domainSoftMarginMM = 0.10;
S06M.queryGuardMode = 'adaptive';
S06M.queryGuardMM = 0.90;
S06M.queryGuardQuantile = 95;
S06M.queryGuardSafetyMM = 0.05;
S06M.queryGuardMinMM = 0.12;
S06M.queryGuardMaxMM = 0.90;
S06M.sensorDomainExpandMMTable = [];
S06M.sensorQueryGuardScaleTable = [5 0.80; 7 0.75];
% sensorDomainExpandMMTable: positive = wider, negative = narrower, applied
% symmetrically to both low/high domains.
% sensorQueryGuardScaleTable: scale < 1 widens both sides by shrinking the
% adaptive query guard for that sensor only.
S06M.dynamicEffectiveMode = 'gradient';
S06M.dynamicTemplateGradientMinRatio = 0.08;
S06M.dynamicTimeGradientMinRatio = 0.15;
S06M.dynamicPeakQuantile = 85;
S06M.defaultSensorThreshold = 0.5;
S06M.weightFloor = 0.05;
S06M.markerSizeRaw = 4;
S06M.markerSizeMask = 7;
S06M.showLegend = true;
S06M.showImpactTable = true;

P = apply_local_options_local(P, S06M);
[WaveformLibrary, sourceInfo] = build_step06m_waveform_inputs_inline_local(P, S06M.bladeId);

B = WaveformLibrary.Blade(S06M.bladeId);
if isempty(B.Sensor)
    error('WaveformLibrary B%d is empty.', S06M.bladeId);
end

lapCount = min(arrayfun(@(s) numel(s.Lap), B.Sensor));
numWindows = floor((lapCount - P.waveform.windowLaps) / P.waveform.slidingStepLaps) + 1;
if numWindows < 1
    error('No valid windows for B%d.', S06M.bladeId);
end

if isempty(S06M.windowIds)
    windowIds = 1:numWindows;
else
    windowIds = S06M.windowIds(:).';
end

method = build_filter_cfg_local(P);

fprintf('\n=== Step06M: sensor point-mask visualization ===\n');
fprintf('Waveform source: %s\n', sourceInfo.waveform_file);
if isfield(sourceInfo, 'note') && ~isempty(sourceInfo.note)
    fprintf('Waveform source note: %s\n', sourceInfo.note);
end
fprintf('Blade: B%d, sensors: %s\n', S06M.bladeId, mat2str(P.sensors.analysis));
fprintf('Windows: %s\n', mat2str(windowIds));

for windowId = windowIds
    lapStart = 1 + (windowId - 1) * P.waveform.slidingStepLaps;
    lapRange = lapStart:(lapStart + P.waveform.windowLaps - 1);
    summaryTable = plot_window_masks_local(B, lapRange, windowId, method, S06M);
    if S06M.showImpactTable
        fprintf('\nStep06M impact summary: B%d W%02d\n', B.blade_id, windowId);
        disp(summaryTable);
    end
end

function P = apply_local_options_local(P, S06M)
P.sensors.analysis = S06M.analysisSensors(:).';
P.region.startTimeSec = S06M.startTimeSec;
P.waveform.windowLaps = S06M.windowLaps;
P.waveform.slidingStepLaps = S06M.slidingStepLaps;
P.waveform.pulsePadSec = S06M.pulsePadSec;
P.waveform.dynamicWindowMode = char(S06M.dynamicWindowMode);
P.waveform.windowRule = char(S06M.windowRule);
P.identification.pulseSelectionMode = S06M.pulseSelectionMode;
P.identification.domainSelectionMode = S06M.domainSelectionMode;
P.identification.domainMarginMM = S06M.domainMarginMM;
P.identification.domainSoftMarginMM = S06M.domainSoftMarginMM;
P.identification.queryGuardMode = S06M.queryGuardMode;
P.identification.queryGuardMM = S06M.queryGuardMM;
P.identification.queryGuardQuantile = S06M.queryGuardQuantile;
P.identification.queryGuardSafetyMM = S06M.queryGuardSafetyMM;
P.identification.queryGuardMinMM = S06M.queryGuardMinMM;
P.identification.queryGuardMaxMM = S06M.queryGuardMaxMM;
P.identification.sensorDomainExpandMMTable = S06M.sensorDomainExpandMMTable;
P.identification.sensorQueryGuardScaleTable = S06M.sensorQueryGuardScaleTable;
P.identification.dynamicEffectiveMode = S06M.dynamicEffectiveMode;
P.identification.dynamicTemplateGradientMinRatio = S06M.dynamicTemplateGradientMinRatio;
P.identification.dynamicTimeGradientMinRatio = S06M.dynamicTimeGradientMinRatio;
P.identification.dynamicPeakQuantile = S06M.dynamicPeakQuantile;
P.identification.defaultSensorThreshold = S06M.defaultSensorThreshold;
P.identification.weightFloor = S06M.weightFloor;
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function [WaveformLibrary, sourceInfo] = build_step06m_waveform_inputs_inline_local(P, bladeIds)
sourceInfo = struct('waveform_file', "", 'note', "");
require_file_local(P.files.lowSpeedFingerprint, 'Step01 low-speed fingerprint');
[HighSpeedNumbering, LowSpeedTemplateLibrary, inlineSourceInfo] = load_step06m_inline_inputs_local(P);
loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;

fileRanges = build_case_file_ranges_inline_local(P.data.highSpeedDir, P.machine.oprChannel, P.machine.sampleRateHz);
rawWindow = build_numbering_time_window_inline_local(HighSpeedNumbering, P, bladeIds);
raw = load_raw_subset_inline_local(P.data.highSpeedDir, fileRanges, P.sensors.analysis, P.machine.sampleRateHz, rawWindow);
oprTimes = load_opr_center_times_inline_local(P);
F_omega_deg = build_phase_speed_inline_local(P, oprTimes);

WaveformLibrary = struct();
WaveformLibrary.mode = 'step06m_inline_high_speed_extraction';
WaveformLibrary.analysis_sensors = P.sensors.analysis;
WaveformLibrary.blade_count = P.machine.bladeCount;
WaveformLibrary.region_selection = HighSpeedNumbering.region_selection;
WaveformLibrary.high_speed_numbering_file = inlineSourceInfo.numbering_file;
WaveformLibrary.extraction_window_mode = P.waveform.dynamicWindowMode;
WaveformLibrary.extraction_window_rule = P.waveform.windowRule;
WaveformLibrary.note = ['Step06M builds the same high-speed waveform slices inline as Step06, ' ...
    'but only for the requested blade windows.'];
WaveformLibrary.Blade = repmat(struct('blade_id', NaN, 'Sensor', []), P.machine.bladeCount, 1);

for bladeId = bladeIds(:).'
    Sensor = repmat(struct('sensor_id', NaN, 'Lap', []), numel(P.sensors.analysis), 1);
    for is = 1:numel(P.sensors.analysis)
        sid = P.sensors.analysis(is);
        sensorNumbering = HighSpeedNumbering.sensor([HighSpeedNumbering.sensor.sensor_id] == sid);
        probe = load_probe_inline_local(P, sid);
        templateSensor = load_template_sensor_inline_local(LowSpeedTemplateLibrary, bladeId, sid);
        thetaStd = read_opr_center_standard_angle_inline_local(LowSpeedReference, sid, bladeId);
        rows = sensorNumbering.selected_rows_by_physical(:, bladeId);
        Lap = repmat(struct( ...
            'lap_id', NaN, ...
            'row_id', NaN, ...
            't', [], ...
            'V', [], ...
            't_peak', NaN, ...
            't_start', NaN, ...
            't_end', NaN, ...
            'x_abs', [], ...
            'x_rel', [], ...
            'theta', []), numel(rows), 1);
        for k = 1:numel(rows)
            rowId = rows(k);
            tPeak = probe.jilublade(rowId, 3);
            [t0, t1] = build_dynamic_segment_window_inline_local(probe.jilublade, rowId, P);
            keep = raw(sid).T >= t0 & raw(sid).T <= t1;
            tSeg = raw(sid).T(keep);
            vSeg = raw(sid).V(keep);
            idxPrev = find(oprTimes < tPeak, 1, 'last');
            if isempty(idxPrev)
                xAbs = nan(size(tSeg));
                thetaRot = nan(size(tSeg));
            else
                thetaPointsDeg = map_segment_to_relative_angle_inline_local(oprTimes(idxPrev), tSeg, F_omega_deg);
                thetaDiffDeg = wrap_to_signed_period_inline_local(thetaPointsDeg - thetaStd, 360);
                xAbs = thetaDiffDeg * (pi / 180) * P.machine.tipRadiusMM;
                thetaRot = map_time_to_rotor_phase_inline_local(oprTimes, tSeg, P.machine.oprPulsesPerRev);
            end
            xRel = xAbs - templateSensor.xc;
            valid = isfinite(tSeg) & isfinite(vSeg) & isfinite(xRel);
            Lap(k).lap_id = k;
            Lap(k).row_id = rowId;
            Lap(k).t = tSeg(valid);
            Lap(k).V = vSeg(valid);
            Lap(k).t_peak = tPeak;
            Lap(k).t_start = t0;
            Lap(k).t_end = t1;
            Lap(k).x_abs = xAbs(valid);
            Lap(k).x_rel = xRel(valid);
            Lap(k).theta = thetaRot(valid);
        end
        Sensor(is).sensor_id = sid;
        Sensor(is).theta_std_deg = thetaStd;
        Sensor(is).template_x = templateSensor.x_grid(:);
        Sensor(is).template_v = templateSensor.v_grid(:);
        Sensor(is).template_dv_dx = templateSensor.dv_dx(:);
        Sensor(is).template_xc = templateSensor.xc;
        Sensor(is).template_domain = templateSensor.x_domain(:).';
        Sensor(is).Lap = Lap;
    end
    WaveformLibrary.Blade(bladeId).blade_id = bladeId;
    WaveformLibrary.Blade(bladeId).Sensor = Sensor;
end

sourceInfo.waveform_file = string(inlineSourceInfo.numbering_file);
sourceInfo.note = sprintf(['inline raw extraction from high-speed data; numbering source: %s; ' ...
    'template source: %s'], inlineSourceInfo.numbering_note, inlineSourceInfo.template_note);
end

function [HighSpeedNumbering, LowSpeedTemplateLibrary, sourceInfo] = load_step06m_inline_inputs_local(P)
sourceInfo = struct('numbering_file', "", 'template_file', "", 'numbering_note', "", 'template_note', "");
sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
exactNumberingFile = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
exactTemplateFile = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');

if exist(exactNumberingFile, 'file') == 2
    loadedNumbering = load(exactNumberingFile, 'HighSpeedNumbering');
    HighSpeedNumbering = loadedNumbering.HighSpeedNumbering;
    sourceInfo.numbering_file = string(exactNumberingFile);
    sourceInfo.numbering_note = "exact Step04 numbering artifact";
else
    [HighSpeedNumbering, sourceFile] = load_subset_high_speed_numbering_inline_local(P);
    sourceInfo.numbering_file = string(sourceFile);
    sourceInfo.numbering_note = "subset sensors cut from compatible Step04 superset artifact";
end

if exist(exactTemplateFile, 'file') == 2
    loadedTemplate = load(exactTemplateFile, 'LowSpeedTemplateLibrary');
    LowSpeedTemplateLibrary = loadedTemplate.LowSpeedTemplateLibrary;
    sourceInfo.template_file = string(exactTemplateFile);
    sourceInfo.template_note = "exact Step05 template artifact";
else
    [LowSpeedTemplateLibrary, sourceFile] = load_compatible_low_speed_template_library_inline_local(P);
    sourceInfo.template_file = string(sourceFile);
    sourceInfo.template_note = "compatible Step05 superset artifact";
end
end

function [HighSpeedNumbering, sourceFile] = load_subset_high_speed_numbering_inline_local(P)
sourceFile = '';
timeLabel = time_label_local(P.region.startTimeSec);
pattern = sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '04_high_speed_numbering', 'S*', pattern));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'HighSpeedNumbering');
    if ~isfield(loaded, 'HighSpeedNumbering') || ~isfield(loaded.HighSpeedNumbering, 'sensor')
        continue;
    end
    available = [loaded.HighSpeedNumbering.sensor.sensor_id];
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing Step04 numbering for sensors %s at %.3f s.', mat2str(requested), P.region.startTimeSec);
end
sourceFile = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(sourceFile, 'HighSpeedNumbering');
HighSpeedNumbering = loaded.HighSpeedNumbering;
keep = ismember([HighSpeedNumbering.sensor.sensor_id], requested);
HighSpeedNumbering.sensor = HighSpeedNumbering.sensor(keep);
end

function [LowSpeedTemplateLibrary, sourceFile] = load_compatible_low_speed_template_library_inline_local(P)
sourceFile = '';
candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'LowSpeedTemplateLibrary');
    if ~isfield(loaded, 'LowSpeedTemplateLibrary')
        continue;
    end
    if isfield(loaded.LowSpeedTemplateLibrary, 'sensor_ids')
        available = unique(loaded.LowSpeedTemplateLibrary.sensor_ids(:).', 'stable');
    elseif isfield(loaded.LowSpeedTemplateLibrary, 'analysis_sensors')
        available = unique(loaded.LowSpeedTemplateLibrary.analysis_sensors(:).', 'stable');
    else
        continue;
    end
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing Step05 template library for sensors %s.', mat2str(requested));
end
sourceFile = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(sourceFile, 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
end

function probe = load_probe_inline_local(P, sid)
probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
loaded = load(probeFile, 'jilublade');
probe = struct('sensor_id', sid, 'jilublade', loaded.jilublade);
end

function templateSensor = load_template_sensor_inline_local(LowSpeedTemplateLibrary, bladeId, sid)
entry = LowSpeedTemplateLibrary.entry(bladeId);
if ~entry.exists
    error('Missing low-speed template bundle for B%d.', bladeId);
end
loadedBundle = load(entry.template_file, 'TemplateBundle');
bundle = loadedBundle.TemplateBundle;
row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
if height(row) ~= 1
    error('Template bundle for B%d does not contain CH%d.', bladeId, sid);
end
loadedSensor = load(char(row.template_file(1)), 'sensor_template');
templateSensor = loadedSensor.sensor_template.Sensor;
end

function thetaStd = read_opr_center_standard_angle_inline_local(LowSpeedReference, sid, bladeId)
if isfield(LowSpeedReference, 'standard_angles_opr_center') && ...
        size(LowSpeedReference.standard_angles_opr_center, 1) >= sid && ...
        size(LowSpeedReference.standard_angles_opr_center, 2) >= bladeId
    thetaStd = LowSpeedReference.standard_angles_opr_center(sid, bladeId);
else
    error('Low-speed reference does not contain OPR-center angle for CH%d B%d.', sid, bladeId);
end
end

function oprTimes = load_opr_center_times_inline_local(P)
oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
require_file_local(oprFile, 'high-speed OPR timing file');
loaded = load(oprFile, 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function F = build_phase_speed_inline_local(P, oprTimes)
omegaFile = fullfile(P.data.highSpeedPulseDir, 'omega.mat');
if exist(omegaFile, 'file') == 2
    omegaVars = whos('-file', omegaFile);
    if ismember('omega', {omegaVars.name})
        loadedOmega = load(omegaFile, 'omega');
        omega = loadedOmega.omega;
        if ~isempty(omega)
            F = griddedInterpolant(omega(:, 1), omega(:, 2), 'linear', 'nearest');
            return;
        end
    end
end
spdT = oprTimes(1:(end - P.machine.oprPulsesPerRev));
spdV = 360 ./ max(oprTimes((P.machine.oprPulsesPerRev + 1):end) - ...
    oprTimes(1:(end - P.machine.oprPulsesPerRev)), eps);
F = griddedInterpolant(spdT, spdV, 'linear', 'nearest');
end

function thetaPoints = map_segment_to_relative_angle_inline_local(tRef, tSeg, F_omega_deg)
if isempty(tSeg)
    thetaPoints = zeros(size(tSeg));
    return;
end
dtFirst = linspace(tRef, tSeg(1), 10);
thetaBase = trapz(dtFirst, F_omega_deg(dtFirst));
wSeg = F_omega_deg(tSeg);
thetaRel = cumtrapz(tSeg, wSeg);
thetaPoints = thetaBase + thetaRel;
end

function thetaRot = map_time_to_rotor_phase_inline_local(oprTimes, sampleTimes, pulsesPerRev)
thetaRot = nan(size(sampleTimes));
if isempty(sampleTimes)
    return;
end
for i = 1:numel(sampleTimes)
    prevIdx = find(oprTimes <= sampleTimes(i), 1, 'last');
    nextIdx = prevIdx + pulsesPerRev;
    if isempty(prevIdx) || nextIdx > numel(oprTimes)
        continue;
    end
    revDt = oprTimes(nextIdx) - oprTimes(prevIdx);
    if revDt <= eps
        continue;
    end
    thetaRot(i) = 2 * pi * (sampleTimes(i) - oprTimes(prevIdx)) / revDt;
end
end

function angle = wrap_to_signed_period_inline_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function timeWindow = build_numbering_time_window_inline_local(HighSpeedNumbering, P, bladeIds)
timeWindow = [inf, -inf];
for is = 1:numel(HighSpeedNumbering.sensor)
    sid = HighSpeedNumbering.sensor(is).sensor_id;
    probe = load_probe_inline_local(P, sid);
    rowsByBlade = HighSpeedNumbering.sensor(is).selected_rows_by_physical(:, bladeIds);
    rows = rowsByBlade(isfinite(rowsByBlade));
    for rowId = rows(:).'
        [t0, t1] = build_dynamic_segment_window_inline_local(probe.jilublade, rowId, P);
        timeWindow(1) = min(timeWindow(1), t0);
        timeWindow(2) = max(timeWindow(2), t1);
    end
end
if ~all(isfinite(timeWindow))
    error('Could not build inline high-speed raw-data window.');
end
end

function [t0, t1] = build_dynamic_segment_window_inline_local(jilublade, row, P)
t0 = jilublade(row, 1) - P.waveform.pulsePadSec;
t1 = jilublade(row, 2) + P.waveform.pulsePadSec;
if ~all(isfinite([t0, t1])) || t1 <= t0
    error('Invalid jilublade row-bound window for row %d: [%.9f, %.9f].', row, t0, t1);
end
end

function fileRanges = build_case_file_ranges_inline_local(caseDir, oprChannel, sampleRateHz)
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
    fileIds(i) = str2double(token{1});
end
[fileIds, order] = sort(fileIds);
files = files(order); %#ok<NASGU>
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = 0;
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_channel_inline_local(caseDir, oprChannel, fileIds(i), sampleRateHz);
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

function raw = load_raw_subset_inline_local(caseDir, fileRanges, sensorIds, sampleRateHz, timeWindow)
raw(max(sensorIds)) = struct('T', [], 'V', []);
useFiles = fileRanges([fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2));
for k = 1:numel(useFiles)
    fileId = useFiles(k).file_id;
    offset = useFiles(k).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_channel_inline_local(caseDir, sid, fileId, sampleRateHz);
        tGlobal = tLocal(:) + offset;
        keep = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        raw(sid).T = [raw(sid).T; tGlobal(keep)]; %#ok<AGROW>
        raw(sid).V = [raw(sid).V; vLocal(keep)]; %#ok<AGROW>
    end
end
end

function [tSec, v] = load_raw_channel_inline_local(caseDir, channelId, fileId, sampleRateHz)
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

function value = safe_divide_local(numerator, denominator)
if ~isfinite(denominator) || abs(denominator) <= eps
    value = NaN;
else
    value = numerator / denominator;
end
end

function method = build_filter_cfg_local(P)
method = struct();
method.pulse_selection_mode = lower(strtrim(P.identification.pulseSelectionMode));
method.domain_selection_mode = lower(strtrim(P.identification.domainSelectionMode));
method.domain_margin_mm = P.identification.domainMarginMM;
method.domain_soft_margin_mm = P.identification.domainSoftMarginMM;
method.query_guard_mode = lower(strtrim(P.identification.queryGuardMode));
method.query_guard_mm = P.identification.queryGuardMM;
method.query_guard_quantile = P.identification.queryGuardQuantile;
method.query_guard_safety_mm = P.identification.queryGuardSafetyMM;
method.query_guard_min_mm = P.identification.queryGuardMinMM;
method.query_guard_max_mm = P.identification.queryGuardMaxMM;
method.sensor_domain_expand_mm_table = P.identification.sensorDomainExpandMMTable;
method.sensor_query_guard_scale_table = P.identification.sensorQueryGuardScaleTable;
method.dynamic_effective_mode = lower(strtrim(P.identification.dynamicEffectiveMode));
method.dynamic_template_gradient_min_ratio = P.identification.dynamicTemplateGradientMinRatio;
method.dynamic_time_gradient_min_ratio = P.identification.dynamicTimeGradientMinRatio;
method.dynamic_peak_quantile = P.identification.dynamicPeakQuantile;
method.default_sensor_threshold = P.identification.defaultSensorThreshold;
method.weight_floor = P.identification.weightFloor;
end

function summaryTable = plot_window_masks_local(B, lapRange, windowId, method, S06M)
sensorsNow = B.Sensor;
keep = ismember([sensorsNow.sensor_id], S06M.analysisSensors);
sensorsNow = sensorsNow(keep);
if isempty(sensorsNow)
    error('No requested sensors found for B%d.', B.blade_id);
end

fig = figure( ...
    'Name', sprintf('Step06M B%d W%02d masks', B.blade_id, windowId), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1, 1, 24, 20], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorsNow), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

rows = repmat(struct( ...
    'BladeID', NaN, ...
    'WindowID', NaN, ...
    'SensorID', NaN, ...
    'RawCount', NaN, ...
    'PulseCount', NaN, ...
    'QueryCount', NaN, ...
    'DynamicCount', NaN, ...
    'FinalCount', NaN, ...
    'QueryLossCount', NaN, ...
    'FinalKeepRatioVsPulse', NaN, ...
    'FinalKeepRatioVsRaw', NaN, ...
    'TemplateWidthMM', NaN, ...
    'PulseWidthMM', NaN, ...
    'DynamicWidthMM', NaN, ...
    'QueryOnlyWidthMM', NaN, ...
    'QueryCoreWidthMM', NaN, ...
    'FinalWidthMM', NaN, ...
    'FinalWidthRatioVsTemplate', NaN, ...
    'QueryGuardMM', NaN), numel(sensorsNow), 1);

for is = 1:numel(sensorsNow)
    Sdata = sensorsNow(is);
    Tpl = template_from_waveform_sensor_local(Sdata);
    [tRaw, xRaw, vRaw, thetaRaw] = concatenate_laps_local(Sdata.Lap(lapRange));
    debug = build_sensor_debug_masks_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, Sdata.sensor_id);

    ax = nexttile;
    hold(ax, 'on');
    plot(ax, xRaw, vRaw, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', S06M.markerSizeRaw, ...
        'DisplayName', sprintf('raw (%d)', nnz(debug.maskRaw)));
    plot(ax, xRaw(debug.maskQueryOnly), vRaw(debug.maskQueryOnly), '.', ...
        'Color', [0 0.45 0.74], 'MarkerSize', S06M.markerSizeMask, ...
        'DisplayName', sprintf('query-safe (%d)', nnz(debug.maskQueryOnly)));
    plot(ax, xRaw(debug.maskDynamicOnly), vRaw(debug.maskDynamicOnly), '.', ...
        'Color', [0.93 0.69 0.13], 'MarkerSize', S06M.markerSizeMask, ...
        'DisplayName', sprintf('dynamic-effective (%d)', nnz(debug.maskDynamicOnly)));
    plot(ax, xRaw(debug.maskFinal), vRaw(debug.maskFinal), '.', ...
        'Color', [0.85 0.33 0.10], 'MarkerSize', S06M.markerSizeMask, ...
        'DisplayName', sprintf('final valid (%d)', nnz(debug.maskFinal)));

    xline(ax, debug.templateDomain(1), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, ...
        'DisplayName', 'template domain');
    xline(ax, debug.templateDomain(2), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    xGrid = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    xSafe = [xGrid(1) + debug.queryGuardMM, xGrid(2) - debug.queryGuardMM];
    xline(ax, xSafe(1), ':', 'Color', [0 0.45 0.74], 'LineWidth', 1.0, ...
        'DisplayName', 'query-safe core');
    xline(ax, xSafe(2), ':', 'Color', [0 0.45 0.74], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    yline(ax, debug.threshold, '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.9, ...
        'DisplayName', sprintf('threshold %.3f V', debug.threshold));

    xlabel(ax, 'x_{rel} (mm)');
    ylabel(ax, 'Voltage (V)');
    title(ax, sprintf(['B%d W%02d CH%d | raw=%d, pulse=%d, query=%d, dynamic=%d, final=%d | ', ...
        'guard=%.3f mm, selWidth=%.3f mm'], ...
        B.blade_id, windowId, Sdata.sensor_id, ...
        nnz(debug.maskRaw), nnz(debug.maskPulseFiniteDomain), nnz(debug.maskQueryOnly), ...
        nnz(debug.maskDynamicOnly), nnz(debug.maskFinal), ...
        debug.queryGuardMM, range_selected_local(xRaw, debug.maskFinal)), ...
        'FontWeight', 'normal');
    grid(ax, 'on');
    box(ax, 'on');
    style_axes_local(ax);
    if S06M.showLegend && is == 1
        legend(ax, 'Location', 'eastoutside', 'Box', 'off');
    end

    pulseCount = nnz(debug.maskPulseFiniteDomain);
    rawCount = nnz(debug.maskRaw);
    finalCount = nnz(debug.maskFinal);
    queryCount = nnz(debug.maskQueryOnly);
    templateWidth = diff(debug.templateDomain);
    queryCoreWidth = diff(xSafe);
    finalWidth = range_selected_local(xRaw, debug.maskFinal);
    rows(is).BladeID = B.blade_id;
    rows(is).WindowID = windowId;
    rows(is).SensorID = Sdata.sensor_id;
    rows(is).RawCount = rawCount;
    rows(is).PulseCount = pulseCount;
    rows(is).QueryCount = queryCount;
    rows(is).DynamicCount = nnz(debug.maskDynamicOnly);
    rows(is).FinalCount = finalCount;
    rows(is).QueryLossCount = max(pulseCount - queryCount, 0);
    rows(is).FinalKeepRatioVsPulse = safe_divide_local(finalCount, pulseCount);
    rows(is).FinalKeepRatioVsRaw = safe_divide_local(finalCount, rawCount);
    rows(is).TemplateWidthMM = templateWidth;
    rows(is).PulseWidthMM = range_selected_local(xRaw, debug.maskPulseFiniteDomain);
    rows(is).DynamicWidthMM = range_selected_local(xRaw, debug.maskDynamicOnly);
    rows(is).QueryOnlyWidthMM = range_selected_local(xRaw, debug.maskQueryOnly);
    rows(is).QueryCoreWidthMM = queryCoreWidth;
    rows(is).FinalWidthMM = finalWidth;
    rows(is).FinalWidthRatioVsTemplate = safe_divide_local(finalWidth, templateWidth);
    rows(is).QueryGuardMM = debug.queryGuardMM;
end
summaryTable = struct2table(rows);
summaryTable = sortrows(summaryTable, {'FinalKeepRatioVsPulse', 'FinalWidthRatioVsTemplate', 'SensorID'}, ...
    {'ascend', 'ascend', 'ascend'});
end

function debug = build_sensor_debug_masks_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, sensorId)
debug = struct();
debug.maskRaw = true(size(tRaw));

f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw(:), 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw(:), 'pchip', NaN);
threshold = method.default_sensor_threshold;
if isfield(Tpl, 'threshold') && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
end
debug.threshold = threshold;

xGrid = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
if strcmpi(method.pulse_selection_mode, 'all')
    [maskPulse, pulseSegmentCount] = isolate_all_pulses_local(vRaw, threshold);
else
    [maskPulse, pulseSegmentCount] = isolate_main_pulse_local(vRaw, threshold);
end
maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, method);
xDomainUse = expand_sensor_domain_local(Tpl.x_domain, method, sensorId);
maskDomain = xRaw >= xDomainUse(1) + method.domain_margin_mm & ...
             xRaw <= xDomainUse(2) - method.domain_margin_mm;
maskFinite = isfinite(f0) & isfinite(fx0) & isfinite(vRaw) & isfinite(thetaRaw);
maskBase = maskPulse & maskEffective & maskDomain & maskFinite;
queryGuardMM = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, maskBase, method, sensorId);
maskQuerySafe = xRaw >= xGrid(1) + queryGuardMM & xRaw <= xGrid(2) - queryGuardMM;
maskFinal = maskBase;
if strcmpi(method.domain_selection_mode, 'hard')
    safeMask = maskFinal & maskQuerySafe;
    if nnz(safeMask) >= 8
        maskFinal = safeMask;
    end
end
if nnz(maskFinal) < 8
    maskFinal = maskEffective & maskDomain & maskFinite;
    if strcmpi(method.domain_selection_mode, 'hard')
        safeMask = maskFinal & maskQuerySafe;
        if nnz(safeMask) >= 8
            maskFinal = safeMask;
        end
    end
end

debug.maskPulseFiniteDomain = maskPulse & maskDomain & maskFinite;
debug.maskQueryOnly = debug.maskPulseFiniteDomain & maskQuerySafe;
debug.maskDynamicOnly = debug.maskPulseFiniteDomain & maskEffective;
debug.maskFinal = maskFinal;
debug.maskQuerySafe = maskQuerySafe;
debug.maskEffective = maskEffective;
debug.templateDomain = xDomainUse;
debug.queryGuardMM = queryGuardMM;
debug.pulseSegmentCount = pulseSegmentCount;
end

function Tpl = template_from_waveform_sensor_local(Sdata)
Tpl = struct();
Tpl.x_grid = Sdata.template_x(:);
Tpl.v_grid = Sdata.template_v(:);
Tpl.dv_dx = Sdata.template_dv_dx(:);
Tpl.x_domain = Sdata.template_domain(:).';
Tpl.threshold = NaN;
end

function [t, x, v, theta] = concatenate_laps_local(Lap)
t = [];
x = [];
v = [];
theta = [];
for k = 1:numel(Lap)
    t = [t; Lap(k).t(:)]; %#ok<AGROW>
    x = [x; Lap(k).x_rel(:)]; %#ok<AGROW>
    v = [v; Lap(k).V(:)]; %#ok<AGROW>
    theta = [theta; Lap(k).theta(:)]; %#ok<AGROW>
end
[t, order] = sort(t);
x = x(order);
v = v(order);
theta = theta(order);
end

function [mask, segmentCount] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segmentCount = numel(starts);
end

function [mask, segmentCount] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segmentCount = 1;
    return;
end
segmentGap = idx(starts(2:end)) - idx(ends(1:end-1));
largeGapThreshold = max(10, round(0.02 * numel(v)));
pulseBreaks = find(segmentGap > largeGapThreshold);
groupStarts = [1; pulseBreaks + 1];
groupEnds = [pulseBreaks; numel(starts)];
for ig = 1:numel(groupStarts)
    segIds = groupStarts(ig):groupEnds(ig);
    bestSeg = segIds(1);
    bestPeak = -inf;
    for iseg = segIds
        seg = idx(starts(iseg):ends(iseg));
        peakVal = max(v(seg));
        if peakVal > bestPeak
            bestPeak = peakVal;
            bestSeg = iseg;
        end
    end
    mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end
segmentCount = numel(groupStarts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end

gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gTplMax = max(gTpl(finite), [], 'omitnan');
if ~isfinite(gTplMax) || gTplMax <= 0
    maskTpl = finite;
else
    maskTpl = gTpl >= method.dynamic_template_gradient_min_ratio * gTplMax;
end

gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if ~isfinite(gTimeMax) || gTimeMax <= 0
    maskTime = false(size(v));
else
    maskTime = gTime >= method.dynamic_time_gradient_min_ratio * gTimeMax;
end
peakLevel = prctile(v(finite), min(max(method.dynamic_peak_quantile, 0), 100));
maskPeak = v >= max(threshold, peakLevel);
mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function queryGuardMM = resolve_query_guard_mm_local(Tpl, x, v, baseMask, method, sensorId)
queryGuardMM = method.query_guard_mm;
if ~strcmpi(method.query_guard_mode, 'adaptive') || nnz(baseMask) < 8
    queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, method, sensorId);
    return;
end
xSel = x(baseMask);
vSel = v(baseMask);
xStat = invert_template_voltage_local(Tpl, vSel, xSel);
uApp = abs(xSel(:) - xStat(:));
uApp = uApp(isfinite(uApp));
if isempty(uApp)
    return;
end
adaptiveGuard = prctile(uApp, min(max(method.query_guard_quantile, 0), 100)) + ...
    method.query_guard_safety_mm;
queryGuardMM = min(max(adaptiveGuard, method.query_guard_min_mm), method.query_guard_max_mm);
queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, method, sensorId);
end

function queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, method, sensorId)
if ~isfield(method, 'sensor_query_guard_scale_table') || isempty(method.sensor_query_guard_scale_table)
    return;
end
tableNow = method.sensor_query_guard_scale_table;
if size(tableNow, 2) < 2
    return;
end
idx = find(tableNow(:, 1) == sensorId, 1);
if isempty(idx)
    return;
end
scale = tableNow(idx, 2);
if ~isfinite(scale) || scale <= 0
    return;
end
queryGuardMM = queryGuardMM * scale;
if isfield(method, 'query_guard_min_mm') && isfinite(method.query_guard_min_mm)
    queryGuardMM = max(queryGuardMM, 0.5 * method.query_guard_min_mm);
end
if isfield(method, 'query_guard_max_mm') && isfinite(method.query_guard_max_mm)
    queryGuardMM = min(queryGuardMM, method.query_guard_max_mm);
end
end

function xDomainUse = expand_sensor_domain_local(xDomain, method, sensorId)
xDomainUse = xDomain;
if ~isfield(method, 'sensor_domain_expand_mm_table') || isempty(method.sensor_domain_expand_mm_table)
    return;
end
tableNow = method.sensor_domain_expand_mm_table;
if size(tableNow, 2) < 2
    return;
end
idx = find(tableNow(:, 1) == sensorId, 1);
if isempty(idx)
    return;
end
expandMM = tableNow(idx, 2);
if ~isfinite(expandMM) || abs(expandMM) <= eps
    return;
end
xCandidate = [xDomain(1) - expandMM, xDomain(2) + expandMM];
if xCandidate(2) <= xCandidate(1)
    return;
end
xDomainUse = xCandidate;
end

function xStat = invert_template_voltage_local(Tpl, v, xRef)
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

function width = range_selected_local(x, mask)
if nnz(mask) < 2
    width = NaN;
else
    xSel = x(mask);
    width = max(xSel) - min(xSel);
end
end

function style_axes_local(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = time_label_local(tSec)
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

