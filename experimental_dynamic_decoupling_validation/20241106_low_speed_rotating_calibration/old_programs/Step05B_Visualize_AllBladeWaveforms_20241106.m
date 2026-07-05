%% Step05B_Visualize_AllBladeWaveforms_20241106
% Show all blade-slot waveforms in the same high-speed window and compare
% them against the low-speed template of the target blade.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
legacyDir = fullfile(validationRoot, '20241106_low_speed_gap_prior_decoupling', 'legacy');
addpath(legacyDir);

%% Parameters
P.targetBlade = 6;
P.analysisSensors = [2 3 5 7];
P.analysisStartTime = 75.0;
P.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', make_time_label_local(P.analysisStartTime), P.templateSuffix);
P.resultSuffix = sprintf('Direct_%s', make_time_label_local(P.analysisStartTime));
P.windowID = [];
P.maxLapsPerBlade = 3;
P.saveFigures = true;

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
templateBundleFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.templateSuffix));
resultFile = fullfile(routeDir, 'output', 'identification', ...
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.resultSuffix));

legacyCfg = Get_20241106_BTT_Config();
dynamicCaseOutputDir = fullfile(legacyCfg.output_root, legacyCfg.dynamic_cases{1});
legacyCfg.dynamic_data_dir = fullfile(legacyCfg.dataset_root, legacyCfg.dynamic_cases{1});
sensorConfigFile = fullfile(legacyCfg.reference_output_dir, 'Sensor_Config_20241106.mat');
figureDir = fullfile(routeDir, 'output', 'figures');
if exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

fprintf('\n=== Step05B: all-blade waveform overview ===\n');
fprintf('Target blade template: B%d, sensors %s\n', P.targetBlade, mat2str(P.analysisSensors));

Template = load_template_from_bundle_local(templateBundleFile, P.analysisSensors);
loadedResult = load(resultFile, 'Result');
Result = loadedResult.Result;
windowID = choose_window_id_local(Result, P);
timeWindow = Result.WindowResult(windowID).time_window;
loadedDynamic = load(Result.DynamicMapFile, 'DynamicMap');
DynamicMap = loadedDynamic.DynamicMap;

loadedSensorConfig = load(sensorConfigFile, 'Sensor_Config');
Sensor_Config = loadedSensorConfig.Sensor_Config;
oprReference = build_opr_reference_from_case_local(dynamicCaseOutputDir, legacyCfg);
F_omega_deg = build_phase_speed_local_from_case_local(dynamicCaseOutputDir, legacyCfg);
rawStream = load_case_raw_streams_window_local(legacyCfg, P.analysisSensors, timeWindow);
probeCache = load_probe_cache_local(dynamicCaseOutputDir, P.analysisSensors, ...
    load_opr_center_times_local(dynamicCaseOutputDir), legacyCfg.opr_pulses_per_rev);
selectedRevIds = get_window_anchor_revolutions_local(DynamicMap, probeCache, windowID);

allBladeData = build_all_blade_waveforms_local(rawStream, dynamicCaseOutputDir, ...
    Sensor_Config, oprReference, F_omega_deg, legacyCfg, Template, P, probeCache, selectedRevIds);
scoreTable = build_blade_score_table_local(allBladeData, P.analysisSensors, P.targetBlade);
disp(scoreTable);

plot_all_blade_grid_local(allBladeData, P, windowID, figureDir);
plot_blade_score_summary_local(scoreTable, P, windowID, figureDir);
writetable(scoreTable, fullfile(figureDir, ...
    sprintf('Step05B_AllBladeScores_B%d_%s_W%02d_20241106.csv', ...
    P.targetBlade, sensorTag, windowID)));

function windowID = choose_window_id_local(Result, P)
if ~isempty(P.windowID)
    windowID = P.windowID;
    return;
end
if isfield(Result, 'BestWindow') && isfield(Result.BestWindow, 'window_id')
    windowID = Result.BestWindow.window_id;
else
    windowID = 1;
end
end

function allBladeData = build_all_blade_waveforms_local(rawStream, caseOutputDir, Sensor_Config, oprReference, F_omega_deg, cfg, Template, P, probeCache, selectedRevIds)
bladeCount = cfg.blades_num;
sensorIds = P.analysisSensors(:).';
allBladeData = repmat(struct('blade_id', [], 'Sensor', []), bladeCount, 1);

for bladeId = 1:bladeCount
    allBladeData(bladeId).blade_id = bladeId;
    sensorStruct = repmat(struct( ...
        'sensor_id', [], ...
        'theta_std_deg', nan, ...
        'laps', [], ...
        'template_corr_mean', nan, ...
        'template_corr_min', nan), numel(sensorIds), 1);

    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        probe = probeCache([probeCache.sensor_id] == sid);
        jilublade = probe.jilublade;
        thetaStd = read_opr_center_standard_angle_local(Sensor_Config, sid, bladeId, oprReference);
        bladeRows = project_rows_by_revolution_local(probe, selectedRevIds, thetaStd);
        bladeRows = bladeRows(1:min(numel(bladeRows), P.maxLapsPerBlade));

        laps = repmat(struct( ...
            'lap_id', [], ...
            'x_rel', [], ...
            'V', [], ...
            'peak_x_mm', nan, ...
            'peak_v', nan, ...
            'corr_with_target_template', nan), numel(bladeRows), 1);

        tplSensor = Template.Sensor([Template.Sensor.sensor_id] == sid);
        corrVals = nan(numel(bladeRows), 1);
        for il = 1:numel(bladeRows)
            rowId = bladeRows(il);
            tPeak = jilublade(rowId, 3);
            [tStart, tEnd] = build_dynamic_segment_window_local(cfg, jilublade, rowId, tPeak);
            keepMask = rawStream(sid).T >= tStart & rawStream(sid).T <= tEnd;
            tSeg = rawStream(sid).T(keepMask);
            vSeg = rawStream(sid).V(keepMask);
            if isempty(tSeg)
                continue;
            end
            idxPrev = find(probe.opr_times < tPeak, 1, 'last');
            if isempty(idxPrev)
                continue;
            end
            xRel = map_segment_to_relative_x_local(probe.opr_times(idxPrev), tSeg, F_omega_deg, thetaStd, cfg.r_tip_mm, tplSensor.xc);
            [xRel, order] = sort(xRel(:));
            vSeg = vSeg(order);

            vTpl = interp1(tplSensor.x_grid(:), tplSensor.v_grid(:), xRel, 'linear', 'extrap');
            valid = isfinite(vSeg) & isfinite(vTpl);
            corrVals(il) = local_corrcoef_scalar(vSeg(valid), vTpl(valid));
            [peakV, peakIdx] = max(vSeg);

            laps(il).lap_id = il;
            laps(il).x_rel = xRel;
            laps(il).V = vSeg;
            laps(il).peak_x_mm = xRel(peakIdx);
            laps(il).peak_v = peakV;
            laps(il).corr_with_target_template = corrVals(il);
        end

        sensorStruct(is).sensor_id = sid;
        sensorStruct(is).theta_std_deg = thetaStd;
        sensorStruct(is).laps = laps;
        sensorStruct(is).template_corr_mean = mean(corrVals, 'omitnan');
        sensorStruct(is).template_corr_min = min(corrVals, [], 'omitnan');
    end

    allBladeData(bladeId).Sensor = sensorStruct;
end
end

function scoreTable = build_blade_score_table_local(allBladeData, sensorIds, targetBlade)
rows = repmat(struct( ...
    'bladeId', nan, ...
    'sensorCount', nan, ...
    'meanTemplateCorr', nan, ...
    'minTemplateCorr', nan, ...
    'meanPeakXAbsMM', nan, ...
    'isTargetBlade', false), numel(allBladeData), 1);

for ib = 1:numel(allBladeData)
    corrVals = [];
    peakVals = [];
    for is = 1:numel(sensorIds)
        laps = allBladeData(ib).Sensor(is).laps;
        for il = 1:numel(laps)
            if isfinite(laps(il).corr_with_target_template)
                corrVals(end + 1, 1) = laps(il).corr_with_target_template; %#ok<AGROW>
            end
            if isfinite(laps(il).peak_x_mm)
                peakVals(end + 1, 1) = abs(laps(il).peak_x_mm); %#ok<AGROW>
            end
        end
    end
    rows(ib).bladeId = allBladeData(ib).blade_id;
    rows(ib).sensorCount = numel(sensorIds);
    rows(ib).meanTemplateCorr = mean(corrVals, 'omitnan');
    rows(ib).minTemplateCorr = min(corrVals, [], 'omitnan');
    rows(ib).meanPeakXAbsMM = mean(peakVals, 'omitnan');
    rows(ib).isTargetBlade = (allBladeData(ib).blade_id == targetBlade);
end

scoreTable = struct2table(rows);
scoreTable = sortrows(scoreTable, {'meanTemplateCorr', 'minTemplateCorr'}, {'descend', 'descend'});
end

function plot_all_blade_grid_local(allBladeData, P, windowID, figureDir)
sensorIds = P.analysisSensors(:).';
fig = figure('Name', sprintf('20241106 all-blade waveform overview B%d W%02d', P.targetBlade, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 24, 18]);
tiledlayout(fig, numel(sensorIds), numel(allBladeData), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    for ib = 1:numel(allBladeData)
        nexttile;
        hold on;
        laps = allBladeData(ib).Sensor(is).laps;
        cmap = lines(max(1, numel(laps)));
        for il = 1:numel(laps)
            plot(laps(il).x_rel, laps(il).V, '-', 'Color', cmap(il, :), 'LineWidth', 1.0);
        end
        title(sprintf('CH%d B%d', sid, allBladeData(ib).blade_id), 'FontWeight', 'normal');
        if ib == 1
            ylabel('V (V)');
        end
        if is == numel(sensorIds)
            xlabel('x rel. target-template center (mm)');
        end
        style_axes_local();
    end
end

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05B_AllBladeWaveforms_B%d_S%s_W%02d_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05B_AllBladeWaveforms_B%d_S%s_W%02d_20241106.pdf', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'ContentType', 'vector');
end
end

function plot_blade_score_summary_local(scoreTable, P, windowID, figureDir)
fig = figure('Name', sprintf('20241106 all-blade score summary B%d W%02d', P.targetBlade, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 15, 10]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(scoreTable.bladeId, scoreTable.meanTemplateCorr, 0.65, 'FaceColor', [0.2 0.45 0.75]); hold on;
targetIdx = find(scoreTable.isTargetBlade, 1, 'first');
if ~isempty(targetIdx)
    plot(scoreTable.bladeId(targetIdx), scoreTable.meanTemplateCorr(targetIdx), 'ro', 'MarkerSize', 7, 'LineWidth', 1.2);
end
xlabel('Blade ID');
ylabel('Mean corr.');
title(sprintf('Target template = B%d', P.targetBlade), 'FontWeight', 'normal');
style_axes_local();

nexttile;
bar(scoreTable.bladeId, scoreTable.minTemplateCorr, 0.65, 'FaceColor', [0.85 0.55 0.18]); hold on;
if ~isempty(targetIdx)
    plot(scoreTable.bladeId(targetIdx), scoreTable.minTemplateCorr(targetIdx), 'ro', 'MarkerSize', 7, 'LineWidth', 1.2);
end
xlabel('Blade ID');
ylabel('Min corr.');
title('Worst-pass correlation across sensors/laps', 'FontWeight', 'normal');
style_axes_local();

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05B_AllBladeScores_B%d_S%s_W%02d_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05B_AllBladeScores_B%d_S%s_W%02d_20241106.pdf', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'ContentType', 'vector');
end
end

function Template = load_template_from_bundle_local(bundleFile, sensorIds)
loaded = load(bundleFile, 'TemplateBundle');
bundle = loaded.TemplateBundle;
sensorEntries = cell(numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
    loadedSensor = load(char(row.template_file(1)), 'sensor_template');
    sensorEntries{i} = loadedSensor.sensor_template.Sensor;
end
Template = struct();
Template.Sensor = vertcat(sensorEntries{:});
end

function rawStream = load_case_raw_streams_window_local(cfg, sensorIds, timeWindow)
fileRanges = build_dynamic_file_ranges_local(cfg);
selectMask = [fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2);
fileRanges = fileRanges(selectMask);
maxSid = max(sensorIds);
rawStream(maxSid) = struct('T', [], 'V', []);
for i = 1:numel(fileRanges)
    fileId = fileRanges(i).file_id;
    offset = fileRanges(i).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_case_channel_local(cfg.dynamic_data_dir, sid, fileId, cfg.pinlv);
        tGlobal = tLocal(:) + offset;
        keepMask = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        rawStream(sid).T = [rawStream(sid).T; tGlobal(keepMask)]; %#ok<AGROW>
        rawStream(sid).V = [rawStream(sid).V; vLocal(keepMask)]; %#ok<AGROW>
    end
end
end

function fileRanges = build_dynamic_file_ranges_local(cfg)
fileIds = list_case_file_ids_local(cfg.dynamic_data_dir, cfg.opr_id);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = [];
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_case_channel_local(cfg.dynamic_data_dir, cfg.opr_id, fileIds(i), cfg.pinlv);
    if isempty(lastEnd)
        offset = 0;
    else
        offset = lastEnd + 1 / cfg.pinlv - tOpr(1);
    end
    lastEnd = tOpr(end) + offset;
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
end
end

function [tSec, v] = load_raw_case_channel_local(caseDir, sid, fileId, pinlv)
filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', sid, fileId));
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / pinlv;
v = raw(:, 2);
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

function F = build_phase_speed_local_from_case_local(caseOutputDir, cfg)
loadedOpr = load(fullfile(caseOutputDir, 'jiluOPR.mat'));
loadedOmega = load(fullfile(caseOutputDir, 'omega.mat'));
if isfield(loadedOmega, 'omega') && ~isempty(loadedOmega.omega)
    omega = loadedOmega.omega;
    F = griddedInterpolant(omega(:, 1), omega(:, 2), 'linear', 'nearest');
    return;
end
jiluOPR = loadedOpr.jiluOPR;
centerTimes = jiluOPR(:, 1);
spdT = centerTimes(1:(end - cfg.opr_pulses_per_rev));
spdV = 360 ./ max(centerTimes((cfg.opr_pulses_per_rev + 1):end) - centerTimes(1:(end - cfg.opr_pulses_per_rev)), eps);
F = griddedInterpolant(spdT, spdV, 'linear', 'nearest');
end

function oprReference = build_opr_reference_from_case_local(caseOutputDir, cfg)
loaded = load(fullfile(caseOutputDir, 'jiluOPR.mat'), 'jiluOPR');
jiluOPR = loaded.jiluOPR;
oprReference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0);
if size(jiluOPR, 2) < 2
    return;
end
centerTime = jiluOPR(:, 1);
startTime = jiluOPR(:, 2);
n = min(numel(centerTime) - cfg.opr_pulses_per_rev, numel(startTime));
if n < 1
    return;
end
dtCenter = centerTime(1:n) - startTime(1:n);
dtRev = centerTime((1:n) + cfg.opr_pulses_per_rev) - centerTime(1:n);
valid = isfinite(dtCenter) & isfinite(dtRev) & dtRev > eps;
shiftDeg = 360 * dtCenter(valid) ./ dtRev(valid);
oprReference.phase_shift_deg = median(shiftDeg, 'omitnan');
oprReference.phase_shift_mm = oprReference.phase_shift_deg * (pi / 180) * cfg.r_tip_mm;
end

function thetaStd = read_opr_center_standard_angle_local(Sensor_Config, sid, bladeId, oprReference)
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    thetaStd = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, bladeId);
    return;
end
thetaStd = Sensor_Config.Standard_Relative_Angles(sid, bladeId);
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
thetaStd = thetaStd - oprReference.phase_shift_deg;
end

function [tStart, tEnd] = build_dynamic_segment_window_local(~, jilublade, rowId, tPeak)
syncCfg = build_single_sync_experiment_config_20241106();
if strcmpi(syncCfg.dynamic_window_mode, 'legacy_row_bounds')
    tStart = jilublade(rowId, 1) - syncCfg.pulse_pad_sec;
    tEnd = jilublade(rowId, 2) + syncCfg.pulse_pad_sec;
else
    tStart = tPeak - syncCfg.pulse_window_sec;
    tEnd = tPeak + syncCfg.pulse_window_sec;
end
end

function xRel = map_segment_to_relative_x_local(tRef, tSeg, F_omega_deg, thetaStd, rTipMm, xc)
thetaPoints = map_segment_to_relative_angle_local(tRef, tSeg, F_omega_deg);
thetaDiff = wrap_to_signed_period_local(thetaPoints - thetaStd, 360);
xAbs = deg2rad(thetaDiff) * rTipMm;
xRel = xAbs - xc;
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

function value = local_corrcoef_scalar(a, b)
if numel(a) < 3 || numel(b) < 3
    value = nan;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function probeCache = load_probe_cache_local(caseOutputDir, sensorIds, oprTimes, oprPulsesPerRev)
probeCache = repmat(struct( ...
    'sensor_id', NaN, ...
    'jilublade', [], ...
    'rel_angles_deg', [], ...
    'prev_opr_index', [], ...
    'revolution_index', [], ...
    'opr_times', []), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    loadedProbe = load(fullfile(caseOutputDir, sprintf('jilublade_probe%d.mat', sid)), 'jilublade');
    jilublade = loadedProbe.jilublade;
    [relAnglesDeg, prevOprIndex, revolutionIndex] = ...
        compute_pulse_reference_geometry_local(jilublade(:, 3), oprTimes, oprPulsesPerRev);
    probeCache(is).sensor_id = sid;
    probeCache(is).jilublade = jilublade;
    probeCache(is).rel_angles_deg = relAnglesDeg;
    probeCache(is).prev_opr_index = prevOprIndex;
    probeCache(is).revolution_index = revolutionIndex;
    probeCache(is).opr_times = oprTimes;
end
end

function oprTimes = load_opr_center_times_local(caseOutputDir)
loaded = load(fullfile(caseOutputDir, 'jiluOPR.mat'), 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function selectedRevIds = get_window_anchor_revolutions_local(DynamicMap, probeCache, windowID)
if ~isfield(DynamicMap, 'SelectionInfo') || ~isfield(DynamicMap.SelectionInfo, 'anchor_sensor_id') || ...
        ~isfield(DynamicMap.SelectionInfo, 'anchor_selected_rows')
    error('DynamicMap.SelectionInfo is missing anchor-row metadata.');
end
anchorSid = DynamicMap.SelectionInfo.anchor_sensor_id;
anchorRows = DynamicMap.SelectionInfo.anchor_selected_rows(:);
lapRange = DynamicMap.Window(windowID).lap_range;
anchorProbe = probeCache([probeCache.sensor_id] == anchorSid);
selectedRevIds = anchorProbe.revolution_index(anchorRows(lapRange));
end

function selectedRows = project_rows_by_revolution_local(probe, selectedRevIds, thetaStd)
selectedRows = nan(numel(selectedRevIds), 1);
for k = 1:numel(selectedRevIds)
    revId = selectedRevIds(k);
    candidates = find(probe.revolution_index == revId & isfinite(probe.rel_angles_deg));
    if isempty(candidates)
        continue;
    end
    candidateError = wrap_to_signed_period_local(probe.rel_angles_deg(candidates) - thetaStd, 360);
    [~, bestPos] = min(abs(candidateError));
    selectedRows(k) = candidates(bestPos);
end
selectedRows = selectedRows(isfinite(selectedRows));
end

function [relAnglesDeg, prevOprIndex, revolutionIndex] = compute_pulse_reference_geometry_local(arrivalTimes, oprTimes, oprPulsesPerRev)
relAnglesDeg = nan(size(arrivalTimes));
prevOprIndex = nan(size(arrivalTimes));
revolutionIndex = nan(size(arrivalTimes));
if isempty(arrivalTimes) || numel(oprTimes) <= oprPulsesPerRev
    return;
end
for i = 1:numel(arrivalTimes)
    tMeas = arrivalTimes(i);
    idx = find(oprTimes < tMeas, 1, 'last');
    if isempty(idx) || idx + oprPulsesPerRev > numel(oprTimes)
        continue;
    end
    t0 = oprTimes(idx);
    t1 = oprTimes(idx + oprPulsesPerRev);
    dtRev = t1 - t0;
    if ~isfinite(dtRev) || dtRev <= eps
        continue;
    end
    relAnglesDeg(i) = 360 * (tMeas - t0) / dtRev;
    prevOprIndex(i) = idx;
    revolutionIndex(i) = idx;
end
end

function angle = wrap_to_signed_period_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function style_axes_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function label = make_time_label_local(tSec)
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
