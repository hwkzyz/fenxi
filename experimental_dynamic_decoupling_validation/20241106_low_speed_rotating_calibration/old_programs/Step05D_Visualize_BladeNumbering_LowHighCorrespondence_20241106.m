%% Step05D_Visualize_BladeNumbering_LowHighCorrespondence_20241106
% Visual audit for the current fingerprint-locked blade numbering route.
% Figure 1: one selected high-speed revolution with local-slot / physical-blade labels.
% Figure 2: low-speed template of the target blade overlaid with the matched
%           high-speed passes from one dynamic window.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
legacyDir = fullfile(validationRoot, '20241106_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacyDir, 'dir') ~= 7
    error('Legacy helper folder not found: %s', legacyDir);
end
addpath(legacyDir);

%% Parameters
P.targetBlade = 6;
P.analysisSensors = [2 3 5 7];
P.analysisStartTime = 75.0;
P.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', make_time_label_local(P.analysisStartTime), P.templateSuffix);
P.windowID = 1;
P.revolutionInWindow = 1;
P.maxOverlayPasses = 3;
P.saveFigures = true;

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
templateBundleFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.templateSuffix));
dynamicMapFile = fullfile(routeDir, 'output', 'dynamic_maps', ...
    sprintf('DynamicMap_B%d_%s_SlidingWindows_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.dynamicSuffix));
numberingFile = fullfile(routeDir, 'output', 'blade_numbering', ...
    sprintf('BladeNumbering_%s_%s_20241106.mat', sensorTag, P.dynamicSuffix));
figureDir = fullfile(routeDir, 'output', 'figures');
if exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

fprintf('\n=== Step05D: blade numbering + low/high correspondence ===\n');
fprintf('Target blade: B%d\n', P.targetBlade);
fprintf('Sensors: %s\n', mat2str(P.analysisSensors));
fprintf('Dynamic map: %s\n', dynamicMapFile);
fprintf('Blade numbering: %s\n', numberingFile);
fprintf('Template bundle: %s\n', templateBundleFile);

DynamicMap = load_required_variable_local(dynamicMapFile, 'DynamicMap');
BladeNumbering = load_required_variable_local(numberingFile, 'BladeNumbering');
Template = load_template_from_bundle_local(templateBundleFile, P.analysisSensors);

legacyCfg = Get_20241106_BTT_Config();
legacyCfg.dynamic_data_dir = fullfile(legacyCfg.dataset_root, legacyCfg.dynamic_cases{1});
caseOutputDir = fullfile(legacyCfg.output_root, legacyCfg.dynamic_cases{1});
syncCfg = build_single_sync_experiment_config_20241106();
oprTimes = load_opr_center_times_local(caseOutputDir);
probeCache = load_probe_cache_local(caseOutputDir, P.analysisSensors, oprTimes, legacyCfg.opr_pulses_per_rev);

selectedRevIds = get_window_anchor_revolutions_local(DynamicMap, probeCache, P.windowID);
if isempty(selectedRevIds)
    error('No anchor revolutions found for dynamic window %d.', P.windowID);
end
if P.revolutionInWindow < 1 || P.revolutionInWindow > numel(selectedRevIds)
    error('revolutionInWindow=%d exceeds the %d laps in window %d.', ...
        P.revolutionInWindow, numel(selectedRevIds), P.windowID);
end
selectedRevId = selectedRevIds(P.revolutionInWindow);
rawAuditWindow = build_revolution_raw_window_local(probeCache, P.analysisSensors, selectedRevId, syncCfg);
rawAuditStream = load_case_raw_streams_window_local(legacyCfg, P.analysisSensors, rawAuditWindow);

summaryTable = build_numbering_summary_table_local(BladeNumbering);
disp(summaryTable);
writetable(summaryTable, fullfile(figureDir, ...
    sprintf('Step05D_NumberingSummary_B%d_%s_W%02d_R%02d_20241106.csv', ...
    P.targetBlade, sensorTag, P.windowID, P.revolutionInWindow)));

plot_numbering_audit_local(rawAuditStream, probeCache, BladeNumbering, ...
    selectedRevId, rawAuditWindow, P, figureDir);
plot_low_high_correspondence_local(DynamicMap, Template, BladeNumbering, P, figureDir);

function S = load_required_variable_local(filepath, varName)
if exist(filepath, 'file') ~= 2
    error('Required file not found: %s', filepath);
end
loaded = load(filepath, varName);
if ~isfield(loaded, varName)
    error('Variable %s is missing in %s.', varName, filepath);
end
S = loaded.(varName);
end

function Template = load_template_from_bundle_local(bundleFile, sensorIds)
loaded = load_required_variable_local(bundleFile, 'TemplateBundle');
bundle = loaded;
sensorEntries = cell(numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
    if isempty(row)
        error('Bundle %s does not contain CH%d.', bundleFile, sid);
    end
    loadedSensor = load_required_variable_local(char(row.template_file(1)), 'sensor_template');
    sensorEntries{i} = loadedSensor.Sensor;
end
Template = struct();
Template.Route = bundle.Route;
Template.TargetBlade = bundle.TargetBlade;
Template.SensorIDs = sensorIds(:).';
Template.Sensor = vertcat(sensorEntries{:});
end

function summaryTable = build_numbering_summary_table_local(BladeNumbering)
sensorIds = BladeNumbering.sensor_ids(:);
physicalToLocal = BladeNumbering.physical_to_local_blade_ids;
localToPhysical = BladeNumbering.local_to_physical_blade_ids;
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'target_blade', NaN, ...
    'target_local_slot', NaN, ...
    'match_score', NaN, ...
    'match_margin', NaN, ...
    'physical_to_local', "", ...
    'local_to_physical', ""), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    rows(i).sensor_id = sensorIds(i);
    rows(i).target_blade = BladeNumbering.target_blade;
    rows(i).target_local_slot = physicalToLocal(BladeNumbering.target_blade, i);
    rows(i).match_score = BladeNumbering.match_score(i);
    rows(i).match_margin = BladeNumbering.match_margin(i);
    rows(i).physical_to_local = string(mat2str(physicalToLocal(:, i).'));
    rows(i).local_to_physical = string(mat2str(localToPhysical(:, i).'));
end
summaryTable = struct2table(rows);
end

function plot_numbering_audit_local(rawStream, probeCache, BladeNumbering, selectedRevId, rawWindow, P, figureDir)
fig = figure('Name', sprintf('20241106 numbering audit B%d W%02d R%02d', ...
    P.targetBlade, P.windowID, P.revolutionInWindow), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 17, 14]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(P.analysisSensors)
    sid = P.analysisSensors(i);
    nexttile;
    hold on;

    probe = probeCache([probeCache.sensor_id] == sid);
    if isempty(probe)
        error('Probe cache is missing CH%d.', sid);
    end
    candidateRows = find(probe.revolution_index == selectedRevId);
    [~, order] = sort(probe.jilublade(candidateRows, 3), 'ascend');
    candidateRows = candidateRows(order);
    bladeCount = numel(candidateRows);

    if bladeCount < 1
        error('CH%d does not contain any pulse in revolution %d.', sid, selectedRevId);
    end

    tRef = min(probe.jilublade(candidateRows, 1));
    keep = rawStream(sid).T >= rawWindow(1) & rawStream(sid).T <= rawWindow(2);
    tPlotMs = 1000 * (rawStream(sid).T(keep) - tRef);
    vPlot = rawStream(sid).V(keep);
    plot(tPlotMs, vPlot, 'k-', 'LineWidth', 0.9, 'DisplayName', 'high-speed raw waveform');

    localToPhysical = BladeNumbering.local_to_physical_blade_ids(:, i);
    targetLocal = BladeNumbering.physical_to_local_blade_ids(P.targetBlade, i);
    cOther = [0.25, 0.45, 0.80];
    cTarget = [0.85, 0.20, 0.20];

    for k = 1:bladeCount
        rowId = candidateRows(k);
        tPeak = probe.jilublade(rowId, 3);
        tStart = probe.jilublade(rowId, 1) - syncCfg_local().pulse_pad_sec;
        tEnd = probe.jilublade(rowId, 2) + syncCfg_local().pulse_pad_sec;
        inSeg = rawStream(sid).T >= tStart & rawStream(sid).T <= tEnd;
        segT = rawStream(sid).T(inSeg);
        segV = rawStream(sid).V(inSeg);
        if isempty(segT)
            continue;
        end

        if k <= numel(localToPhysical) && isfinite(localToPhysical(k))
            physicalBlade = localToPhysical(k);
        else
            physicalBlade = NaN;
        end
        if physicalBlade == P.targetBlade
            cNow = cTarget;
            lw = 1.4;
        else
            cNow = cOther;
            lw = 1.0;
        end

        plot(1000 * (segT - tRef), segV, '-', 'Color', cNow, 'LineWidth', lw, ...
            'HandleVisibility', 'off');
        xline(1000 * (tPeak - tRef), ':', 'Color', cNow, 'LineWidth', 0.8, 'HandleVisibility', 'off');

        [peakV, peakPos] = max(segV);
        peakTms = 1000 * (segT(peakPos) - tRef);
        scatter(peakTms, peakV, 20, cNow, 'filled', 'HandleVisibility', 'off');
        labelText = sprintf('L%d / P%d', k, physicalBlade);
        text(peakTms, peakV, ['  ' labelText], ...
            'Color', cNow, 'FontName', 'Times New Roman', 'FontSize', 7.8, ...
            'VerticalAlignment', 'bottom', 'Clipping', 'on');
    end

    xlabel('Time in selected revolution (ms)');
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('CH%d  rev %d  B%d -> L%d  score=%.3f', ...
        sid, selectedRevId, P.targetBlade, targetLocal, BladeNumbering.match_score(i)), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 8.5);
    subtitle(sprintf('P->L %s', mat2str(BladeNumbering.physical_to_local_blade_ids(:, i).')), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 7.8);
    style_axes_local();
end

sgtitle(sprintf('Fingerprint-locked numbering audit, B%d, W%02d lap %d, raw window %.6f-%.6f s', ...
    P.targetBlade, P.windowID, P.revolutionInWindow, rawWindow(1), rawWindow(2)), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.saveFigures
    pngFile = fullfile(figureDir, sprintf( ...
        'Step05D_NumberingAudit_B%d_%s_W%02d_R%02d_20241106.png', ...
        P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID, P.revolutionInWindow));
    pdfFile = fullfile(figureDir, sprintf( ...
        'Step05D_NumberingAudit_B%d_%s_W%02d_R%02d_20241106.pdf', ...
        P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID, P.revolutionInWindow));
    exportgraphics(fig, pngFile, 'Resolution', 300);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
end
end

function plot_low_high_correspondence_local(DynamicMap, Template, BladeNumbering, P, figureDir)
Wmap = DynamicMap.Window(P.windowID);
fig = figure('Name', sprintf('20241106 low/high correspondence B%d W%02d', ...
    P.targetBlade, P.windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 17, 18]);
tiledlayout(fig, numel(P.analysisSensors), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(P.analysisSensors)
    sid = P.analysisSensors(i);
    nexttile;
    hold on;

    rawSensor = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
    tplSensor = Template.Sensor([Template.Sensor.sensor_id] == sid);
    if isempty(rawSensor) || isempty(tplSensor)
        error('Window/template is missing CH%d.', sid);
    end

    segments = split_sensor_points_into_passes_local(rawSensor.t);
    maxPasses = min(P.maxOverlayPasses, numel(segments));
    plot(tplSensor.x_grid, tplSensor.v_grid, 'k-', 'LineWidth', 1.5, ...
        'DisplayName', 'low-speed template');
    cmap = lines(max(1, maxPasses));
    corrText = strings(maxPasses, 1);

    for ip = 1:maxPasses
        rows = segments{ip};
        [xRel, order] = sort(rawSensor.x_rel(rows));
        vObs = rawSensor.V(rows);
        vObs = vObs(order);
        plot(xRel, vObs, '-', 'Color', cmap(ip, :), 'LineWidth', 1.0, ...
            'DisplayName', sprintf('high-speed pass %d', ip));
        [peakV, peakPos] = max(vObs);
        scatter(xRel(peakPos), peakV, 18, cmap(ip, :), 'filled', 'HandleVisibility', 'off');

        vTpl = interp1(tplSensor.x_grid(:), tplSensor.v_grid(:), xRel(:), 'linear', 'extrap');
        valid = isfinite(vObs(:)) & isfinite(vTpl(:));
        corrVal = local_corrcoef_scalar_local(vObs(valid), vTpl(valid));
        corrText(ip) = sprintf('P%d: r=%.3f', ip, corrVal);
    end

    xline(tplSensor.x_domain(1), ':', 'Color', [0.45, 0.45, 0.45], 'HandleVisibility', 'off');
    xline(tplSensor.x_domain(2), ':', 'Color', [0.45, 0.45, 0.45], 'HandleVisibility', 'off');
    xlabel('x relative to template center (mm)');
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('CH%d  physical B%d  local slot L%d', sid, P.targetBlade, ...
        BladeNumbering.physical_to_local_blade_ids(P.targetBlade, i)), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 8.5);
    subtitle(strjoin(cellstr(corrText(corrText ~= "")), ', '), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 7.8);
    if i == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();
end

sgtitle(sprintf('Low-speed template and matched high-speed B%d waveform, window %d', ...
    P.targetBlade, P.windowID), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.saveFigures
    pngFile = fullfile(figureDir, sprintf( ...
        'Step05D_LowHighCorrespondence_B%d_%s_W%02d_20241106.png', ...
        P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID));
    pdfFile = fullfile(figureDir, sprintf( ...
        'Step05D_LowHighCorrespondence_B%d_%s_W%02d_20241106.pdf', ...
        P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID));
    exportgraphics(fig, pngFile, 'Resolution', 300);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
end
end

function loadWindow = build_revolution_raw_window_local(probeCache, sensorIds, revId, syncCfg)
loadWindow = [inf, -inf];
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    probe = probeCache([probeCache.sensor_id] == sid);
    candidateRows = find(probe.revolution_index == revId);
    if isempty(candidateRows)
        continue;
    end
    tStart = min(probe.jilublade(candidateRows, 1)) - syncCfg.pulse_pad_sec;
    tEnd = max(probe.jilublade(candidateRows, 2)) + syncCfg.pulse_pad_sec;
    loadWindow(1) = min(loadWindow(1), tStart);
    loadWindow(2) = max(loadWindow(2), tEnd);
end
if ~all(isfinite(loadWindow)) || loadWindow(2) <= loadWindow(1)
    error('Failed to build the raw window for revolution %d.', revId);
end
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
anchorSid = DynamicMap.SelectionInfo.anchor_sensor_id;
anchorRows = DynamicMap.SelectionInfo.anchor_selected_rows(:);
lapRange = DynamicMap.Window(windowID).lap_range;
anchorProbe = probeCache([probeCache.sensor_id] == anchorSid);
selectedRevIds = anchorProbe.revolution_index(anchorRows(lapRange));
selectedRevIds = selectedRevIds(isfinite(selectedRevIds));
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

function segments = split_sensor_points_into_passes_local(t)
t = t(:);
if isempty(t)
    segments = {};
    return;
end
[tSorted, order] = sort(t);
dt = diff(tSorted);
positiveDt = dt(dt > 0);
if isempty(positiveDt)
    segments = {order};
    return;
end
gapThreshold = max(5 * median(positiveDt, 'omitnan'), prctile(positiveDt, 95));
breaks = find(dt > gapThreshold);
edges = [0; breaks(:); numel(tSorted)];
segments = cell(numel(edges) - 1, 1);
for i = 1:numel(segments)
    idx = (edges(i) + 1):edges(i + 1);
    segments{i} = order(idx);
end
end

function value = local_corrcoef_scalar_local(a, b)
if numel(a) < 3 || numel(b) < 3
    value = nan;
    return;
end
C = corrcoef(a(:), b(:));
if numel(C) < 4
    value = nan;
else
    value = C(1, 2);
end
end

function style_axes_local()
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'LineWidth', 0.8);
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

function cfg = syncCfg_local()
cfg = build_single_sync_experiment_config_20241106();
end
