%% Step07D_AuditEO12AcrossBlades_20241106
% Check whether EO12 appears naturally for other blades/windows without
% changing the Step06 solver.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local parameter block. Edit here directly.
A08D = struct();
A08D.analysisSensors = [5 7];
A08D.startTimeSec = 75.0;
A08D.bladeIds = 1:6;
A08D.windowIds = [];              % [] -> all windows
A08D.targetEO = 12;
A08D.freqSearchHz = [0 1000];
A08D.eoPad = 2;
A08D.viewEnable = true;

P = apply_local_options_local(P, A08D);
[FilteredWaveformLibrary, sourceInfo] = load_compatible_filtered_waveform_library_local(P);
[WaveformLibrary, sourceInfo] = load_paired_waveform_library_local(sourceInfo);

EOAuditTable = build_eo_audit_table_local(FilteredWaveformLibrary, WaveformLibrary, P, A08D);
EOAuditSummary = build_eo_audit_summary_local(EOAuditTable, A08D.targetEO);

fprintf('\n=== Step07D: EO12 audit across blades ===\n');
fprintf('Sensors: %s\n', mat2str(P.sensors.analysis));
fprintf('Filtered waveform source: %s\n', sourceInfo.filtered_file);
fprintf('Waveform library source: %s\n', sourceInfo.waveform_file);
disp(EOAuditSummary);
disp(EOAuditTable);

if A08D.viewEnable
    plot_eo_rank_heatmap_local(EOAuditTable, A08D);
    plot_eo_gap_heatmap_local(EOAuditTable, A08D);
end

function P = apply_local_options_local(P, A08D)
P.sensors.analysis = A08D.analysisSensors(:).';
P.region.startTimeSec = A08D.startTimeSec;

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
P.files.filteredWaveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('FilteredWaveformLibrary_%s_20241106.mat', timeLabel));
end

function [FilteredWaveformLibrary, sourceInfo] = load_compatible_filtered_waveform_library_local(P)
sourceInfo = struct('filtered_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');
if exist(P.files.filteredWaveformLibrary, 'file') == 2
    loaded = load(P.files.filteredWaveformLibrary, 'FilteredWaveformLibrary');
    FilteredWaveformLibrary = loaded.FilteredWaveformLibrary;
    sourceInfo.filtered_file = string(P.files.filteredWaveformLibrary);
    return;
end

timeLabel = time_label_local(P.region.startTimeSec);
pattern = sprintf('FilteredWaveformLibrary_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '05_waveform_library', 'S*', pattern));
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'FilteredWaveformLibrary');
    if ~isfield(loaded, 'FilteredWaveformLibrary') || ~isfield(loaded.FilteredWaveformLibrary, 'analysis_sensors')
        continue;
    end
    available = unique(loaded.FilteredWaveformLibrary.analysis_sensors(:).', 'stable');
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing compatible filtered waveform library for sensors %s.', mat2str(requested));
end
fileNow = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(fileNow, 'FilteredWaveformLibrary');
FilteredWaveformLibrary = loaded.FilteredWaveformLibrary;
sourceInfo.filtered_file = string(fileNow);
end

function [WaveformLibrary, sourceInfo] = load_paired_waveform_library_local(sourceInfo)
waveformFile = strrep(char(sourceInfo.filtered_file), 'FilteredWaveformLibrary_', 'WaveformLibrary_');
if exist(waveformFile, 'file') ~= 2
    error('Paired waveform library not found for filtered file:%s%s', newline, waveformFile);
end
loaded = load(waveformFile, 'WaveformLibrary');
if ~isfield(loaded, 'WaveformLibrary')
    error('Variable WaveformLibrary is missing in file:%s%s', newline, waveformFile);
end
WaveformLibrary = loaded.WaveformLibrary;
sourceInfo.waveform_file = string(waveformFile);
end

function T = build_eo_audit_table_local(FilteredWaveformLibrary, WaveformLibrary, P, A08D)
rows = [];
for bladeId = A08D.bladeIds(:).'
    B = FilteredWaveformLibrary.Blade(bladeId);
    if isempty(B.Window)
        continue;
    end
    for iw = 1:numel(B.Window)
        W = B.Window(iw);
        if ~isempty(A08D.windowIds) && ~ismember(W.window_id, A08D.windowIds)
            continue;
        end
        bundle = subset_bundle_by_sensor_local(W.Bundle, A08D.analysisSensors);
        if isempty(bundle) || bundle.point_count < 20
            continue;
        end
        bundle = attach_template_and_speed_local(bundle, WaveformLibrary, bladeId, W.lap_range, P);
        eoCandidates = build_eo_candidates_local(bundle.rot_freq_mean_hz, A08D.freqSearchHz, A08D.eoPad);
        seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
        eoList = [seedTable.EO];
        rankTarget = find(eoList == A08D.targetEO, 1, 'first');
        top3 = eoList(1:min(3, numel(eoList)));
        row = struct();
        row.BladeID = bladeId;
        row.WindowID = W.window_id;
        row.LapStart = W.lap_range(1);
        row.LapEnd = W.lap_range(end);
        row.RotFreqHz = bundle.rot_freq_mean_hz;
        row.Top1EO = seedTable(1).EO;
        row.Top1RMSE = seedTable(1).weighted_voltage_rmse;
        row.Top3EO = string(mat2str(top3));
        row.TargetEO = A08D.targetEO;
        row.TargetEORank = rank_or_nan_local(rankTarget);
        row.TargetEORMSE = rmse_or_nan_local(seedTable, rankTarget);
        row.TargetEOGapV = row.TargetEORMSE - row.Top1RMSE;
        rows = [rows; row]; %#ok<AGROW>
    end
end
T = struct2table(rows);
end

function bundle = subset_bundle_by_sensor_local(bundleIn, requestedSensors)
bundle = bundleIn;
requestedSensors = requestedSensors(:).';
available = unique(bundleIn.sensor_ids(:).', 'stable');
sensorKeep = available(ismember(available, requestedSensors));
if isempty(sensorKeep)
    bundle = [];
    return;
end

pointMask = ismember(bundleIn.S(:).', sensorKeep);
pointFields = {'X', 'T', 'T_rel', 'V', 'S', 'W', 'Theta', 'F0', 'Fx', 'sensor_index'};
for i = 1:numel(pointFields)
    name = pointFields{i};
    if isfield(bundle, name) && numel(bundle.(name)) == numel(pointMask)
        bundle.(name) = bundle.(name)(pointMask);
    end
end

bundle.sensor_ids = sensorKeep;
sensorIndex = nan(size(bundle.S));
for is = 1:numel(sensorKeep)
    sensorIndex(bundle.S == sensorKeep(is)) = is;
end
bundle.sensor_index = sensorIndex;
bundle.point_count = nnz(pointMask);

if isfield(bundle, 'window_meta') && ~isempty(bundle.window_meta)
    metaKeep = ismember([bundle.window_meta.sensor_id], sensorKeep);
    bundle.window_meta = bundle.window_meta(metaKeep);
end

if isfield(bundle, 'valid_segment_count')
    if isfield(bundle, 'window_meta') && ~isempty(bundle.window_meta) && isfield(bundle.window_meta, 'pulse_segment_count')
        bundle.valid_segment_count = sum([bundle.window_meta.pulse_segment_count], 'omitnan');
    else
        bundle.valid_segment_count = nnz(isfinite(bundle.V));
    end
end
end

function bundle = attach_template_and_speed_local(bundle, WaveformLibrary, bladeId, lapRange, P)
bundle.interp_v = cell(numel(bundle.sensor_ids), 1);
bundle.x_domain_by_sensor = nan(numel(bundle.sensor_ids), 2);
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    sensorsNow = WaveformLibrary.Blade(bladeId).Sensor;
    sensorMask = [sensorsNow.sensor_id] == sid;
    Sdata = sensorsNow(sensorMask);
    if isempty(Sdata)
        error('WaveformLibrary B%d is missing CH%d.', bladeId, sid);
    end
    x = Sdata.template_x(:);
    v = Sdata.template_v(:);
    [x, ia] = unique(x, 'stable');
    v = v(ia);
    bundle.interp_v{is} = griddedInterpolant(x, v, 'pchip', 'none');
    if isfield(Sdata, 'template_domain') && numel(Sdata.template_domain) >= 2
        bundle.x_domain_by_sensor(is, :) = Sdata.template_domain(:).';
    else
        bundle.x_domain_by_sensor(is, :) = [min(x), max(x)];
    end
end
[rotFreqHz, rotRpm] = estimate_window_speed_local(WaveformLibrary, bladeId, lapRange, bundle.sensor_ids);
bundle.rot_freq_mean_hz = rotFreqHz;
bundle.rot_rpm_mean = rotRpm;
bundle.overshoot_penalty_weight = P.identification.overshootPenaltyWeight;
bundle.sensor_eta_limit_mm = P.identification.sensorEtaLimitMM;
bundle.sensor_eta_reg_weight_v_per_mm = P.identification.sensorEtaRegWeightVPerMM;
end

function [rotFreqHz, rotRpm] = estimate_window_speed_local(WaveformLibrary, bladeId, lapRange, sensorIds)
periods = [];
for sid = sensorIds(:).'
    sensorsNow = WaveformLibrary.Blade(bladeId).Sensor;
    sensorMask = [sensorsNow.sensor_id] == sid;
    Sdata = sensorsNow(sensorMask);
    if isempty(Sdata) || numel(Sdata.Lap) < max(lapRange)
        continue;
    end
    tPeak = arrayfun(@(x) x.t_peak, Sdata.Lap(lapRange));
    dt = diff(tPeak(:));
    periods = [periods; dt(isfinite(dt) & dt > eps)]; %#ok<AGROW>
end
if isempty(periods)
    error('Could not estimate rotation speed for B%d laps %s.', bladeId, mat2str(lapRange));
end
rotFreqHz = 1 / median(periods, 'omitnan');
rotRpm = 60 * rotFreqHz;
end

function Tsum = build_eo_audit_summary_local(T, targetEO)
bladeIds = unique(T.BladeID).';
rows = repmat(struct( ...
    'BladeID', NaN, ...
    'WindowCount', NaN, ...
    'Top1Count', NaN, ...
    'Top3Count', NaN, ...
    'MedianRank', NaN, ...
    'MedianGapV', NaN, ...
    'TargetEO', targetEO), numel(bladeIds), 1);
for i = 1:numel(bladeIds)
    b = bladeIds(i);
    Tb = T(T.BladeID == b, :);
    rows(i).BladeID = b;
    rows(i).WindowCount = height(Tb);
    rows(i).Top1Count = nnz(Tb.TargetEORank == 1);
    rows(i).Top3Count = nnz(Tb.TargetEORank <= 3);
    rows(i).MedianRank = median(Tb.TargetEORank, 'omitnan');
    rows(i).MedianGapV = median(Tb.TargetEOGapV, 'omitnan');
end
Tsum = struct2table(rows);
end

function plot_eo_rank_heatmap_local(T, A08D)
bladeIds = unique(T.BladeID).';
windowIds = unique(T.WindowID).';
rankMat = nan(numel(bladeIds), numel(windowIds));
for ib = 1:numel(bladeIds)
    for iw = 1:numel(windowIds)
        row = T(T.BladeID == bladeIds(ib) & T.WindowID == windowIds(iw), :);
        if height(row) == 1
            rankMat(ib, iw) = row.TargetEORank;
        end
    end
end
fig = figure('Name', sprintf('Step07D EO%d rank across blades', A08D.targetEO), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 20, 9], 'NumberTitle', 'off');
imagesc(windowIds, bladeIds, rankMat);
set(gca, 'YDir', 'normal');
xlabel('Window ID');
ylabel('Blade ID');
title(sprintf('EO%d rank in VP seed table', A08D.targetEO), 'FontWeight', 'normal');
cb = colorbar;
cb.Label.String = 'Rank (1 = best)';
box on;
style_axes_local();
end

function plot_eo_gap_heatmap_local(T, A08D)
bladeIds = unique(T.BladeID).';
windowIds = unique(T.WindowID).';
gapMat = nan(numel(bladeIds), numel(windowIds));
for ib = 1:numel(bladeIds)
    for iw = 1:numel(windowIds)
        row = T(T.BladeID == bladeIds(ib) & T.WindowID == windowIds(iw), :);
        if height(row) == 1
            gapMat(ib, iw) = row.TargetEOGapV;
        end
    end
end
fig = figure('Name', sprintf('Step07D EO%d gap across blades', A08D.targetEO), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 20, 9], 'NumberTitle', 'off');
imagesc(windowIds, bladeIds, gapMat);
set(gca, 'YDir', 'normal');
xlabel('Window ID');
ylabel('Blade ID');
title(sprintf('EO%d RMSE gap relative to best EO', A08D.targetEO), 'FontWeight', 'normal');
cb = colorbar;
cb.Label.String = 'Target EO RMSE - best RMSE (V)';
box on;
style_axes_local();
end

function eoCandidates = build_eo_candidates_local(rotFreqHz, freqSearchHz, eoPad)
rotFreqHz = max(rotFreqHz, eps);
freqLo = min(freqSearchHz);
freqHi = max(freqSearchHz);
eoMin = max(1, ceil(freqLo / rotFreqHz) - max(0, floor(eoPad)));
eoMax = max(eoMin, floor(freqHi / rotFreqHz) + max(0, floor(eoPad)));
eoAll = eoMin:eoMax;
freqAll = eoAll .* rotFreqHz;
eoCandidates = eoAll(freqAll >= freqLo & freqAll <= freqHi);
if isempty(eoCandidates)
    eoCenter = max(1, round(mean([freqLo, freqHi]) / rotFreqHz));
    eoCandidates = max(1, eoCenter - eoPad):max(1, eoCenter + eoPad);
end
eoCandidates = unique(round(eoCandidates(:).'));
end

function seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'weighted_voltage_rmse', inf, 'plain_voltage_rmse', inf, ...
    'linear_vp_rmse', inf), numel(eoCandidates), 1);
finite = isfinite(bundle.V) & isfinite(bundle.F0) & isfinite(bundle.Fx) & ...
    isfinite(bundle.W) & isfinite(bundle.Theta);
if nnz(finite) < 3
    error('Too few finite points for VP seed scan.');
end
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo * bundle.Theta(finite));
    c1 = cos(eo * bundle.Theta(finite));
    basis = [-bundle.Fx(finite), -bundle.Fx(finite) .* s1, -bundle.Fx(finite) .* c1];
    wSqrt = sqrt(bundle.W(finite));
    y = bundle.V(finite) - bundle.F0(finite);
    coeff = (basis .* wSqrt) \ (y .* wSqrt);
    dxC = coeff(1);
    eta = zeros(1, numel(bundle.sensor_ids));
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y - basis * coeff;
    [obj, plainRmse] = template_synchronous_objective_local([A, phi, dxC, eta], eo, bundle);
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx_c = dxC;
    rows(i).sensor_eta = eta;
    rows(i).sensor_eta_max_abs = max(abs(eta), [], 'omitnan');
    rows(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    rows(i).plain_voltage_rmse = plainRmse;
    rows(i).linear_vp_rmse = sqrt(sum(bundle.W(finite) .* linRes.^2) / max(sum(bundle.W(finite)), eps));
end
score = [[rows.weighted_voltage_rmse].', [rows.linear_vp_rmse].', [rows.plain_voltage_rmse].', [rows.EO].'];
[~, order] = sortrows(score, [1, 2, 3, 4]);
seedTable = rows(order);
end

function [obj, plainRmse, vPred, uEst, coverage] = template_synchronous_objective_local(params, eo, bundle)
A = params(1);
phi = params(2);
dxC = params(3);
params = bound_sensor_eta_params_local(params, abs(bundle.sensor_eta_limit_mm), numel(bundle.sensor_ids));
eta = params(4:end).';
uEst = A .* sin(eo .* bundle.Theta + phi);
etaSample = eta(bundle.sensor_index(:));
xIn = bundle.X - dxC - etaSample - uEst;
vPred = nan(size(bundle.V));
clamped = false(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    xLo = bundle.x_domain_by_sensor(is, 1);
    xHi = bundle.x_domain_by_sensor(is, 2);
    xEvalRaw = xIn(mask);
    xEval = min(max(xEvalRaw, xLo), xHi);
    clamped(mask) = xEvalRaw < xLo | xEvalRaw > xHi;
    overshoot(mask) = max(xLo - xEvalRaw, 0) + max(xEvalRaw - xHi, 0);
    vPred(mask) = bundle.interp_v{is}(xEval);
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
res = bundle.V(valid) - vPred(valid);
overshootPenalty = bundle.overshoot_penalty_weight * sum(bundle.W(valid) .* (overshoot(valid) .^ 2));
etaRegPenalty = bundle.point_count * (bundle.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta .^ 2))) ^ 2;
obj = sum(bundle.W(valid) .* (res .^ 2)) + overshootPenalty + etaRegPenalty + 1e6 * nnz(~valid);
plainRmse = sqrt(mean(res .^ 2));
coverage = struct();
coverage.clamp_fraction = nnz(clamped) / max(numel(clamped), 1);
coverage.max_query_overshoot_mm = max(overshoot, [], 'omitnan');
coverage.overshoot_penalty = overshootPenalty;
coverage.sensor_eta_reg_penalty = etaRegPenalty;
end

function p = bound_sensor_eta_params_local(p, sensorEtaLimit, nSensor)
need = 3 + nSensor;
if numel(p) < need
    p(numel(p)+1:need) = 0;
end
p = p(1:need);
eta = p(4:end);
eta = max(min(eta, sensorEtaLimit), -sensorEtaLimit);
if ~isempty(eta)
    eta = eta - eta(1);
    eta = max(min(eta, sensorEtaLimit), -sensorEtaLimit);
end
p(4:end) = eta;
end

function y = rank_or_nan_local(rankValue)
if isempty(rankValue)
    y = NaN;
else
    y = rankValue;
end
end

function y = rmse_or_nan_local(seedTable, rankValue)
if isempty(rankValue)
    y = NaN;
else
    y = seedTable(rankValue).weighted_voltage_rmse;
end
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
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

