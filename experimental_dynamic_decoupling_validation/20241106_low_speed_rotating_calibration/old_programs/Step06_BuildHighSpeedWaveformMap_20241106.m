%% Step06_BuildHighSpeedWaveformMap_20241106
% Build a physical-blade waveform library from the shared region numbering.
% The library is sliced by blade after numbering, not by rerunning numbering.
%
% Semantic role in the 20241106 main chain:
%   This is the high-speed waveform-region extraction step.
%   It corresponds to the "dynamic map / waveform-window construction"
%   role that is Step02 in the 20251222 OPRCenterStd route.

clear; close all; clc;

P = NewFlow_Config_20241106();

% Step06 local parameter block.
% Edit here first so this script can be rerun directly without jumping back
% to NewFlow_Config_20241106.
S06 = struct();
S06.analysisSensors = [2 3 5 7];
S06.startTimeSec = 75.0;
S06.windowLaps = 3;
S06.slidingStepLaps = 1;
S06.pulsePadSec = 2 / 5e6;
S06.dynamicWindowMode = 'legacy_row_bounds';
S06.windowRule = 't_start=jilublade(row,1)-pulsePadSec; t_end=jilublade(row,2)+pulsePadSec';
S06.pulseSelectionMode = 'single';
S06.domainSelectionMode = 'hard';
S06.domainMarginMM = 0.02;
S06.domainSoftMarginMM = 0.10;
S06.queryGuardMode = 'adaptive';
S06.queryGuardMM = 0.90;
S06.queryGuardQuantile = 95;
S06.queryGuardSafetyMM = 0.05;
S06.queryGuardMinMM = 0.12;
S06.queryGuardMaxMM = 0.90;
S06.sensorDomainExpandMMTable = [2 0.25; 3 0.25];
S06.sensorQueryGuardScaleTable = [2 0.80; 3 0.80];
% Each row is [sensor_id, scale]. scale < 1 widens both sides by shrinking
% the adaptive query guard for that sensor only.
S06.dynamicEffectiveMode = 'gradient';
S06.dynamicTemplateGradientMinRatio = 0.08;
S06.dynamicTimeGradientMinRatio = 0.15;
S06.dynamicPeakQuantile = 85;
S06.defaultSensorThreshold = 0.5;
S06.weightFloor = 0.05;
S06.outputLabel = '';
S06.viewEnable = true;
P = apply_step06_local_options_local(P, S06);

require_file_local(P.files.lowSpeedFingerprint, 'Step01 low-speed fingerprint');
[HighSpeedNumbering, LowSpeedTemplateLibrary, step06SourceInfo] = load_step06_inputs_local(P);
loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;

outFile = P.files.waveformLibrary;
outDir = fileparts(outFile);
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

fileRanges = build_case_file_ranges_local(P.data.highSpeedDir, P.machine.oprChannel, P.machine.sampleRateHz);
rawWindow = build_numbering_time_window_local(HighSpeedNumbering, P);
raw = load_raw_subset_local(P.data.highSpeedDir, fileRanges, P.sensors.analysis, P.machine.sampleRateHz, rawWindow);
oprTimes = load_opr_center_times_local(P);
F_omega_deg = build_phase_speed_local(P, oprTimes);

WaveformLibrary = struct();
WaveformLibrary.mode = 'physical_blade_slices_from_shared_region_numbering';
WaveformLibrary.analysis_sensors = P.sensors.analysis;
WaveformLibrary.blade_count = P.machine.bladeCount;
WaveformLibrary.region_selection = HighSpeedNumbering.region_selection;
WaveformLibrary.high_speed_numbering_file = P.files.highSpeedNumbering;
WaveformLibrary.extraction_window_mode = P.waveform.dynamicWindowMode;
WaveformLibrary.extraction_window_rule = P.waveform.windowRule;
WaveformLibrary.note = ['Blade(b).Sensor(s).Lap(k) contains the same physical blade b across all sensors. ' ...
    'Each raw high-speed pulse is sliced by jilublade(row,1:2) plus pulsePadSec, not by a fixed t_peak window.'];
WaveformLibrary.Blade = repmat(struct('blade_id', NaN, 'Sensor', []), P.machine.bladeCount, 1);

for bladeId = 1:P.machine.bladeCount
    Sensor = repmat(struct('sensor_id', NaN, 'Lap', []), numel(P.sensors.analysis), 1);
    for is = 1:numel(P.sensors.analysis)
        sid = P.sensors.analysis(is);
        sensorNumbering = HighSpeedNumbering.sensor([HighSpeedNumbering.sensor.sensor_id] == sid);
        probe = load_probe_local(P, sid);
        templateSensor = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid);
        thetaStd = read_opr_center_standard_angle_local(LowSpeedReference, sid, bladeId);
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
            [t0, t1] = build_dynamic_segment_window_local(probe.jilublade, rowId, tPeak, P);
            keep = raw(sid).T >= t0 & raw(sid).T <= t1;
            tSeg = raw(sid).T(keep);
            vSeg = raw(sid).V(keep);
            tPeak = probe.jilublade(rowId, 3);
            idxPrev = find(oprTimes < tPeak, 1, 'last');
            if isempty(idxPrev)
                xAbs = nan(size(tSeg));
                thetaRot = nan(size(tSeg));
            else
                thetaPointsDeg = map_segment_to_relative_angle_local(oprTimes(idxPrev), tSeg, F_omega_deg);
                thetaDiffDeg = wrap_to_signed_period_local(thetaPointsDeg - thetaStd, 360);
                xAbs = thetaDiffDeg * (pi / 180) * P.machine.tipRadiusMM;
                thetaRot = map_time_to_rotor_phase_local(oprTimes, tSeg, P.machine.oprPulsesPerRev);
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

method = build_waveform_filter_method_local(P);
[FilteredWaveformLibrary, FilteredSummary] = build_filtered_waveform_library_local( ...
    WaveformLibrary, LowSpeedTemplateLibrary, P, method);

save(outFile, 'WaveformLibrary', '-v7.3');
save(P.files.filteredWaveformLibrary, 'FilteredWaveformLibrary', '-v7.3');
writetable(build_waveform_summary_local(WaveformLibrary), strrep(outFile, '.mat', '.csv'));
writetable(FilteredSummary, P.files.filteredWaveformSummary);
export_single_sensor_waveform_libraries_local(WaveformLibrary, FilteredWaveformLibrary, P);

fprintf('\n=== Step05: waveform library in region ===\n');
fprintf('Raw extraction window: [%.6f, %.6f] s\n', rawWindow(1), rawWindow(2));
fprintf('High-speed numbering source: %s\n', step06SourceInfo.numbering_file);
fprintf('Low-speed template source: %s\n', step06SourceInfo.template_file);
if isfield(step06SourceInfo, 'note') && ~isempty(step06SourceInfo.note)
    fprintf('Step06 source note: %s\n', step06SourceInfo.note);
end
fprintf('Saved row-bound waveform library: %s\n', outFile);
fprintf('Saved filtered identification-input waveform library: %s\n', P.files.filteredWaveformLibrary);
fprintf('Saved filtered waveform summary: %s\n', P.files.filteredWaveformSummary);
fprintf('Saved single-sensor waveform libraries under: %s\n', ...
    fullfile(P.outputDir, '05_waveform_library', 'by_sensor'));
if P.view.enable
    visualize_filtered_waveform_summary_local(FilteredSummary, P);
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function P = apply_step06_local_options_local(P, S06)
P.sensors.analysis = S06.analysisSensors(:).';
P.region.startTimeSec = S06.startTimeSec;
P.waveform.windowLaps = S06.windowLaps;
P.waveform.slidingStepLaps = S06.slidingStepLaps;
P.waveform.pulsePadSec = S06.pulsePadSec;
P.waveform.dynamicWindowMode = S06.dynamicWindowMode;
P.waveform.windowRule = char(S06.windowRule);
P.identification.pulseSelectionMode = S06.pulseSelectionMode;
P.identification.domainSelectionMode = S06.domainSelectionMode;
P.identification.domainMarginMM = S06.domainMarginMM;
P.identification.domainSoftMarginMM = S06.domainSoftMarginMM;
P.identification.queryGuardMode = S06.queryGuardMode;
P.identification.queryGuardMM = S06.queryGuardMM;
P.identification.queryGuardQuantile = S06.queryGuardQuantile;
P.identification.queryGuardSafetyMM = S06.queryGuardSafetyMM;
P.identification.queryGuardMinMM = S06.queryGuardMinMM;
P.identification.queryGuardMaxMM = S06.queryGuardMaxMM;
P.identification.sensorDomainExpandMMTable = S06.sensorDomainExpandMMTable;
P.identification.sensorQueryGuardScaleTable = S06.sensorQueryGuardScaleTable;
P.identification.dynamicEffectiveMode = S06.dynamicEffectiveMode;
P.identification.dynamicTemplateGradientMinRatio = S06.dynamicTemplateGradientMinRatio;
P.identification.dynamicTimeGradientMinRatio = S06.dynamicTimeGradientMinRatio;
P.identification.dynamicPeakQuantile = S06.dynamicPeakQuantile;
P.identification.defaultSensorThreshold = S06.defaultSensorThreshold;
P.identification.weightFloor = S06.weightFloor;
P.view.enable = logical(S06.viewEnable);
P.step06 = struct();
P.step06.outputLabel = char(S06.outputLabel);
P = refresh_step06_artifact_paths_local(P);
end

function P = refresh_step06_artifact_paths_local(P)
sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
P.files.lowSpeedNumbering = fullfile(P.outputDir, '02_low_speed_numbering', sensorTag, ...
    'LowSpeedNumbering_20241106.mat');
P.files.regionSelection = fullfile(P.outputDir, '03_region_selection', ...
    sprintf('HighSpeedRegion_%s_20241106.mat', timeLabel));
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
P.files.waveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('WaveformLibrary_%s_20241106.mat', timeLabel));
P.files.filteredWaveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('FilteredWaveformLibrary_%s_20241106.mat', timeLabel));
P.files.filteredWaveformSummary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('FilteredWaveformSummary_%s_20241106.csv', timeLabel));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');
end

function [HighSpeedNumbering, LowSpeedTemplateLibrary, sourceInfo] = load_step06_inputs_local(P)
sourceInfo = struct('numbering_file', "", 'template_file', "", 'note', "");

if exist(P.files.highSpeedNumbering, 'file') == 2
    loaded = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
    HighSpeedNumbering = loaded.HighSpeedNumbering;
    sourceInfo.numbering_file = string(P.files.highSpeedNumbering);
else
    [HighSpeedNumbering, sourceInfo.numbering_file] = ...
        load_subset_high_speed_numbering_local(P);
end

if exist(P.files.lowSpeedTemplateLibrary, 'file') == 2
    loadedTemplate = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
    LowSpeedTemplateLibrary = loadedTemplate.LowSpeedTemplateLibrary;
    sourceInfo.template_file = string(P.files.lowSpeedTemplateLibrary);
else
    [LowSpeedTemplateLibrary, sourceInfo.template_file] = ...
        load_compatible_low_speed_template_library_local(P);
end

sourceInfo.note = "exact artifacts if available; otherwise subset/compatible superset artifacts";
end

function [HighSpeedNumbering, sourceFile] = load_subset_high_speed_numbering_local(P)
timeLabel = time_label_local(P.region.startTimeSec);
pattern = sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '04_high_speed_numbering', 'S*', pattern));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    filePath = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(filePath, 'HighSpeedNumbering');
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
    error(['Missing Step04 high-speed numbering for sensors %s at %.3f s, and no ' ...
        'compatible superset numbering artifact was found.'], ...
        mat2str(requested), P.region.startTimeSec);
end
sourceFile = string(fullfile(candidates(bestIdx).folder, candidates(bestIdx).name));
loaded = load(char(sourceFile), 'HighSpeedNumbering');
HighSpeedNumbering = loaded.HighSpeedNumbering;
keep = ismember([HighSpeedNumbering.sensor.sensor_id], requested);
HighSpeedNumbering.sensor = HighSpeedNumbering.sensor(keep);
end

function [LowSpeedTemplateLibrary, sourceFile] = load_compatible_low_speed_template_library_local(P)
candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    filePath = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(filePath, 'LowSpeedTemplateLibrary');
    if ~isfield(loaded, 'LowSpeedTemplateLibrary') || ~isfield(loaded.LowSpeedTemplateLibrary, 'sensor_ids')
        continue;
    end
    available = unique(loaded.LowSpeedTemplateLibrary.sensor_ids(:).', 'stable');
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error(['Missing Step05 low-speed template library for sensors %s, and no compatible ' ...
        'superset template artifact was found.'], mat2str(requested));
end
sourceFile = string(fullfile(candidates(bestIdx).folder, candidates(bestIdx).name));
loaded = load(char(sourceFile), 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
end

function probe = load_probe_local(P, sid)
probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
loaded = load(probeFile, 'jilublade');
probe = struct('sensor_id', sid, 'jilublade', loaded.jilublade);
end

function templateSensor = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid)
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

function thetaStd = read_opr_center_standard_angle_local(LowSpeedReference, sid, bladeId)
if isfield(LowSpeedReference, 'standard_angles_opr_center') && ...
        size(LowSpeedReference.standard_angles_opr_center, 1) >= sid && ...
        size(LowSpeedReference.standard_angles_opr_center, 2) >= bladeId
    thetaStd = LowSpeedReference.standard_angles_opr_center(sid, bladeId);
else
    error('Low-speed reference does not contain OPR-center angle for CH%d B%d.', sid, bladeId);
end
end

function oprTimes = load_opr_center_times_local(P)
oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
require_file_local(oprFile, 'high-speed OPR timing file');
loaded = load(oprFile, 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function F = build_phase_speed_local(P, oprTimes)
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

function thetaPoints = map_segment_to_relative_angle_local(tRef, tSeg, F_omega_deg)
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

function thetaRot = map_time_to_rotor_phase_local(oprTimes, sampleTimes, pulsesPerRev)
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

function angle = wrap_to_signed_period_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function timeWindow = build_numbering_time_window_local(HighSpeedNumbering, P)
timeWindow = [inf, -inf];
for is = 1:numel(HighSpeedNumbering.sensor)
    sid = HighSpeedNumbering.sensor(is).sensor_id;
    probe = load_probe_local(P, sid);
    rows = HighSpeedNumbering.sensor(is).selected_rows_by_physical;
    rows = rows(isfinite(rows));
    for rowId = rows(:).'
        tPeak = probe.jilublade(rowId, 3);
        [t0, t1] = build_dynamic_segment_window_local(probe.jilublade, rowId, tPeak, P);
        timeWindow(1) = min(timeWindow(1), t0);
        timeWindow(2) = max(timeWindow(2), t1);
    end
end
if ~all(isfinite(timeWindow))
    error('Could not build waveform-library raw-data window.');
end
end

function [t0, t1] = build_dynamic_segment_window_local(jilublade, row, ~, P)
t0 = jilublade(row, 1) - P.waveform.pulsePadSec;
t1 = jilublade(row, 2) + P.waveform.pulsePadSec;
if ~all(isfinite([t0, t1])) || t1 <= t0
    error('Invalid jilublade row-bound window for row %d: [%.9f, %.9f].', row, t0, t1);
end
end

function T = build_waveform_summary_local(WaveformLibrary)
rows = [];
for b = 1:numel(WaveformLibrary.Blade)
    B = WaveformLibrary.Blade(b);
    for is = 1:numel(B.Sensor)
        S = B.Sensor(is);
        pointCounts = arrayfun(@(x) numel(x.t), S.Lap);
        rows = [rows; struct( ... %#ok<AGROW>
            'BladeID', B.blade_id, ...
            'SensorID', S.sensor_id, ...
            'LapCount', numel(S.Lap), ...
            'MinPointCount', min(pointCounts), ...
            'MedianPointCount', median(pointCounts), ...
            'MaxPointCount', max(pointCounts))];
    end
end
T = struct2table(rows);
end

function method = build_waveform_filter_method_local(P)
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
if ~ismember(method.pulse_selection_mode, {'single', 'all'})
    error('P.identification.pulseSelectionMode must be "single" or "all".');
end
if ~ismember(method.domain_selection_mode, {'hard', 'soft'})
    error('P.identification.domainSelectionMode must be "hard" or "soft".');
end
if ~ismember(method.query_guard_mode, {'fixed', 'adaptive'})
    error('P.identification.queryGuardMode must be "fixed" or "adaptive".');
end
if ~ismember(method.dynamic_effective_mode, {'legacy', 'gradient'})
    error('P.identification.dynamicEffectiveMode must be "legacy" or "gradient".');
end
end

function [FilteredWaveformLibrary, Summary] = build_filtered_waveform_library_local( ...
        WaveformLibrary, LowSpeedTemplateLibrary, P, method)
FilteredWaveformLibrary = struct();
FilteredWaveformLibrary.mode = 'step05_filtered_identification_input_waveforms';
FilteredWaveformLibrary.source_waveform_library = P.files.waveformLibrary;
FilteredWaveformLibrary.extraction_window_mode = P.waveform.dynamicWindowMode;
FilteredWaveformLibrary.extraction_window_rule = P.waveform.windowRule;
FilteredWaveformLibrary.note = ['Step06 first stores jilublade row-bound waveforms, then filters them ' ...
    'with old Step03-style pulse/domain/effective/query-guard rules. Step07 should consume this library.'];
FilteredWaveformLibrary.analysis_sensors = P.sensors.analysis;
FilteredWaveformLibrary.window_laps = P.waveform.windowLaps;
FilteredWaveformLibrary.sliding_step_laps = P.waveform.slidingStepLaps;
FilteredWaveformLibrary.method = method;
FilteredWaveformLibrary.Blade = repmat(struct('blade_id', NaN, 'Window', []), ...
    P.machine.bladeCount, 1);

rows = [];
for bladeId = 1:P.machine.bladeCount
    templateEntry = LowSpeedTemplateLibrary.entry(bladeId);
    if ~templateEntry.exists
        warning('Skipping B%d waveform filtering because its low-speed template bundle is missing.', bladeId);
        continue;
    end
    B = WaveformLibrary.Blade(bladeId);
    lapCount = min(arrayfun(@(s) numel(s.Lap), B.Sensor));
    numWindows = floor((lapCount - P.waveform.windowLaps) / P.waveform.slidingStepLaps) + 1;
    if numWindows < 1
        warning('Skipping B%d waveform filtering because only %d laps are available.', bladeId, lapCount);
        continue;
    end
    Window = repmat(empty_filtered_window_local(), numWindows, 1);
    for w = 1:numWindows
        lapStart = 1 + (w - 1) * P.waveform.slidingStepLaps;
        lapRange = lapStart:(lapStart + P.waveform.windowLaps - 1);
        [bundle, sensorMeta] = build_filtered_bundle_for_blade_window_local(B, lapRange, method);
        Window(w).window_id = w;
        Window(w).lap_range = lapRange;
        Window(w).Bundle = bundle;
        Window(w).SensorMeta = sensorMeta;

        for im = 1:numel(sensorMeta)
            M = sensorMeta(im);
            rows = [rows; struct( ... %#ok<AGROW>
                'BladeID', bladeId, ...
                'WindowID', w, ...
                'SensorID', M.sensor_id, ...
                'LapRange', mat2str(lapRange), ...
                'RawPointCount', M.raw_points, ...
                'SelectedPointCount', M.valid_points, ...
                'PulseSegmentCount', M.pulse_segment_count, ...
                'DynamicEffectivePoints', M.dynamic_effective_points, ...
                'QuerySafePoints', M.query_safe_points, ...
                'QueryGuardMM', M.query_guard_mm, ...
                'SelectedXMinMM', M.selected_x_range(1), ...
                'SelectedXMaxMM', M.selected_x_range(2), ...
                'ShapeRMSE', M.shape_rmse, ...
                'TemplateFile', templateEntry.template_file, ...
                'Status', 'step06_filtered_ready_for_step07_identification')];
        end
    end
    FilteredWaveformLibrary.Blade(bladeId).blade_id = bladeId;
    FilteredWaveformLibrary.Blade(bladeId).Window = Window;
end
Summary = struct2table(rows);
end

function export_single_sensor_waveform_libraries_local(WaveformLibrary, FilteredWaveformLibrary, P)
timeLabel = time_label_local(P.region.startTimeSec);
baseDir = fullfile(P.outputDir, '05_waveform_library', 'by_sensor');
for sid = P.sensors.analysis(:).'
    sensorDir = fullfile(baseDir, sprintf('CH%d', sid));
    if exist(sensorDir, 'dir') ~= 7
        mkdir(sensorDir);
    end
    WaveformLibrarySingle = filter_waveform_library_single_sensor_local(WaveformLibrary, sid);
    FilteredWaveformLibrarySingle = filter_filtered_waveform_library_single_sensor_local( ...
        FilteredWaveformLibrary, sid);
    save(fullfile(sensorDir, sprintf('WaveformLibrary_%s_CH%d_20241106.mat', timeLabel, sid)), ...
        'WaveformLibrarySingle', '-v7.3');
    save(fullfile(sensorDir, sprintf('FilteredWaveformLibrary_%s_CH%d_20241106.mat', timeLabel, sid)), ...
        'FilteredWaveformLibrarySingle', '-v7.3');
end
end

function Wout = filter_waveform_library_single_sensor_local(Win, sid)
Wout = Win;
Wout.analysis_sensors = sid;
for ib = 1:numel(Win.Blade)
    S = Win.Blade(ib).Sensor;
    keep = [S.sensor_id] == sid;
    Wout.Blade(ib).Sensor = S(keep);
end
end

function Fout = filter_filtered_waveform_library_single_sensor_local(Fin, sid)
Fout = Fin;
Fout.analysis_sensors = sid;
for ib = 1:numel(Fin.Blade)
    B = Fin.Blade(ib);
    for iw = 1:numel(B.Window)
        Wsrc = B.Window(iw);
        bundle = Wsrc.Bundle;
        sensorMeta = Wsrc.SensorMeta;
        keepSensor = [sensorMeta.sensor_id] == sid;
        sensorMeta = sensorMeta(keepSensor);
        keepPoint = bundle.S == sid;
        bundle.X = bundle.X(keepPoint);
        bundle.T = bundle.T(keepPoint);
        bundle.T_rel = bundle.T_rel(keepPoint);
        bundle.V = bundle.V(keepPoint);
        bundle.S = bundle.S(keepPoint);
        bundle.W = bundle.W(keepPoint);
        bundle.Theta = bundle.Theta(keepPoint);
        bundle.F0 = bundle.F0(keepPoint);
        bundle.Fx = bundle.Fx(keepPoint);
        bundle.sensor_index = ones(nnz(keepPoint), 1);
        bundle.sensor_ids = sid;
        bundle.window_meta = sensorMeta;
        bundle.valid_segment_count = sum([sensorMeta.pulse_segment_count]);
        bundle.point_count = numel(bundle.T);
        Wsrc.Bundle = bundle;
        Wsrc.SensorMeta = sensorMeta;
        B.Window(iw) = Wsrc;
    end
    Fout.Blade(ib) = B;
end
end

function W = empty_filtered_window_local()
W = struct('window_id', NaN, 'lap_range', [], 'Bundle', [], 'SensorMeta', []);
end

function [bundle, sensorMeta] = build_filtered_bundle_for_blade_window_local(B, lapRange, method)
X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
F0 = [];
Fx = [];
sensorIndex = [];
sensorMeta = repmat(empty_sensor_meta_local(), numel(B.Sensor), 1);

for is = 1:numel(B.Sensor)
    Sdata = B.Sensor(is);
    Tpl = template_from_waveform_sensor_local(Sdata);
    [tRaw, xRaw, vRaw, thetaRaw] = concatenate_laps_local(Sdata.Lap(lapRange));
    [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        filter_sensor_waveform_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, Sdata.sensor_id);
    meta.sensor_id = Sdata.sensor_id;
    sensorMeta(is) = meta;

    X = [X; xSel(:)]; %#ok<AGROW>
    T = [T; tSel(:)]; %#ok<AGROW>
    V = [V; vSel(:)]; %#ok<AGROW>
    S = [S; repmat(Sdata.sensor_id, numel(tSel), 1)]; %#ok<AGROW>
    W = [W; wSel(:)]; %#ok<AGROW>
    Theta = [Theta; thetaSel(:)]; %#ok<AGROW>
    F0 = [F0; f0Sel(:)]; %#ok<AGROW>
    Fx = [Fx; fxSel(:)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, numel(tSel), 1)]; %#ok<AGROW>
end

if isempty(T)
    error('No valid filtered waveform points were built for B%d laps %s.', ...
        B.blade_id, mat2str(lapRange));
end

bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.sensor_index = sensorIndex;
bundle.sensor_ids = [B.Sensor.sensor_id];
bundle.window_meta = sensorMeta;
bundle.valid_segment_count = sum([sensorMeta.pulse_segment_count]);
bundle.point_count = numel(T);
bundle.time_window = [min(T), max(T)];
bundle.selection_pass = 'step05_old_step03_core_filtering';
end

function meta = empty_sensor_meta_local()
meta = struct( ...
    'sensor_id', NaN, ...
    'raw_points', 0, ...
    'valid_points', 0, ...
    'pulse_segment_count', 0, ...
    'dynamic_effective_points', 0, ...
    'query_safe_points', 0, ...
    'query_guard_mm', NaN, ...
    'domain_selection_mode', '', ...
    'domain_soft_margin_mm', NaN, ...
    'template_domain', [NaN NaN], ...
    'selected_x_range', [NaN NaN], ...
    'shape_rmse', NaN);
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

function [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        filter_sensor_waveform_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, sensorId)
meta = empty_sensor_meta_local();
meta.raw_points = numel(tRaw);
meta.domain_selection_mode = method.domain_selection_mode;
meta.domain_soft_margin_mm = method.domain_soft_margin_mm;
meta.template_domain = expand_sensor_domain_local(Tpl.x_domain, method, sensorId);

f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw(:), 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw(:), 'pchip', NaN);
if isfield(Tpl, 'threshold') && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = method.default_sensor_threshold;
end
xDomainGrid = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
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
maskQuerySafe = xRaw >= xDomainGrid(1) + queryGuardMM & ...
                xRaw <= xDomainGrid(2) - queryGuardMM;
mask = maskBase;
if strcmpi(method.domain_selection_mode, 'hard')
    safeMask = mask & maskQuerySafe;
    if nnz(safeMask) >= 8
        mask = safeMask;
    end
end
if nnz(mask) < 8
    mask = maskEffective & maskDomain & maskFinite;
    if strcmpi(method.domain_selection_mode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
end

tSel = tRaw(mask);
xSel = xRaw(mask);
vSel = vRaw(mask);
thetaSel = thetaRaw(mask);
f0Sel = f0(mask);
fxSel = fx0(mask);
wSel = build_filtered_weight_local(tSel, xSel, vSel, Tpl, xDomainGrid, queryGuardMM, method);

if isempty(tSel)
    meta.valid_points = 0;
else
    vStatic = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xSel(:), 'pchip', NaN);
    meta.valid_points = numel(tSel);
    meta.selected_x_range = [min(xSel), max(xSel)];
    meta.shape_rmse = sqrt(mean((vSel(:) - vStatic(:)).^2, 'omitnan'));
end
meta.pulse_segment_count = pulseSegmentCount;
meta.dynamic_effective_points = nnz(maskEffective);
meta.query_safe_points = nnz(maskQuerySafe(mask));
meta.query_guard_mm = queryGuardMM;
end

function wTotal = build_filtered_weight_local(t, x, v, Tpl, xDomain, queryGuardMM, method)
wEdge = build_edge_weight_local(t, v, method.weight_floor);
wDomain = build_domain_soft_weight_local(x, xDomain, method.domain_soft_margin_mm, method.weight_floor);
wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMM, method.domain_soft_margin_mm, method.weight_floor);
wGradient = build_template_gradient_weight_local(Tpl, x, method.domain_selection_mode, method.weight_floor);
wTotal = max(method.weight_floor, wEdge(:) .* wDomain(:) .* wQuery(:) .* wGradient(:));
if max(wTotal) > 0
    wTotal = max(method.weight_floor, wTotal ./ max(wTotal));
end
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

function xDomainUse = expand_sensor_domain_local(xDomain, method, sensorId)
xDomainUse = xDomain(:).';
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
extraMM = tableNow(idx, 2);
if ~isfinite(extraMM) || extraMM <= 0
    return;
end
xDomainUse = [xDomainUse(1) - extraMM, xDomainUse(2) + extraMM];
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

function wEdge = build_edge_weight_local(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_local(x, xDomain, marginMM, floorW)
if marginMM <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMM, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMM, marginMM, floorW)
if queryGuardMM <= 0 || marginMM <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMM - marginMM, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMM, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_local(Tpl, x, domainSelectionMode, floorW)
if ~strcmpi(domainSelectionMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end

function visualize_filtered_waveform_summary_local(FilteredSummary, P)
figDir = fullfile(P.view.figureDir, '05_waveform_library');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end
bladeIds = unique(FilteredSummary.BladeID).';
sensorIds = unique(FilteredSummary.SensorID).';
pointMat = nan(numel(sensorIds), numel(bladeIds));
ratioMat = nan(numel(sensorIds), numel(bladeIds));
for is = 1:numel(sensorIds)
    for ib = 1:numel(bladeIds)
        rows = FilteredSummary(FilteredSummary.SensorID == sensorIds(is) & ...
            FilteredSummary.BladeID == bladeIds(ib), :);
        if height(rows) >= 1
            pointMat(is, ib) = median(rows.SelectedPointCount, 'omitnan');
            ratioMat(is, ib) = median(rows.SelectedPointCount ./ max(rows.RawPointCount, 1), 'omitnan');
        end
    end
end

fig = figure('Name', 'Step06 filtered waveform summary', 'Color', 'w');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
imagesc(bladeIds, sensorIds, pointMat);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Median selected points');
colorbar;
box on;

nexttile;
imagesc(bladeIds, sensorIds, ratioMat);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Selected / raw point ratio');
colorbar;
box on;

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step06V_FilteredWaveformSummary_20241106.png'), 'Resolution', 300);
end
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

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end
