%% Step07A_AuditIdentificationChain_20241106
% Audit the whole identification chain before changing the Step06 solver.
% Focus:
%   1) is the low-speed template for the chosen blade/sensors trustworthy,
%   2) is the high-speed physical-blade numbering for the chosen blade stable,
%   3) do the chosen Step06 windows contain enough per-sensor information,
%   4) do the Step06 window-wise results look self-consistent.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local audit parameter block. Edit here directly.
A07 = struct();
A07.analysisSensors = [5 7];
A07.startTimeSec = 75.0;
A07.outputLabel = 'B6_only';
A07.bladeId = 6;
A07.windowIds = [];                 % [] -> audit all windows in Step06 result
A07.maxLapsToShow = 6;              % for numbering-time plots
A07.viewEnable = true;
A07.saveFigures = false;

P = apply_local_audit_options_local(P, A07);

[LowSpeedTemplateLibrary, LowSpeedTemplateAudit, sourceTpl] = ...
    load_compatible_template_artifacts_local(P);
[HighSpeedNumbering, sourceNum] = load_compatible_high_speed_numbering_local(P);
IdentificationResult = load_identification_result_local(P);

bladeId = A07.bladeId;
windowTable = build_window_audit_table_local(IdentificationResult, bladeId, P.sensors.analysis, A07.windowIds);
templateTable = build_template_audit_table_local(LowSpeedTemplateAudit, bladeId, P.sensors.analysis);
numberingTable = build_numbering_lap_table_local(HighSpeedNumbering, P, bladeId, P.sensors.analysis);

fprintf('\n=== Step07A: identification-chain audit ===\n');
fprintf('Blade: B%d\n', bladeId);
fprintf('Sensors: %s\n', mat2str(P.sensors.analysis));
fprintf('Template source: %s\n', sourceTpl.template_file);
fprintf('Template audit source: %s\n', sourceTpl.audit_file);
fprintf('High-speed numbering source: %s\n', sourceNum.numbering_file);
fprintf('Identification source: %s\n', P.files.identificationResult);
disp(templateTable);
disp(windowTable);

if P.view.enable
    plot_low_speed_template_audit_local(LowSpeedTemplateLibrary, LowSpeedTemplateAudit, bladeId, P.sensors.analysis);
    plot_high_speed_numbering_audit_local(numberingTable, bladeId, P.sensors.analysis, A07.maxLapsToShow);
    plot_identification_window_audit_local(windowTable, bladeId, P.sensors.analysis);
    plot_step06_bundle_domain_audit_local(IdentificationResult, bladeId, P.sensors.analysis, A07.windowIds);
end

function P = apply_local_audit_options_local(P, A07)
P.sensors.analysis = A07.analysisSensors(:).';
P.region.startTimeSec = A07.startTimeSec;
P.view.enable = logical(A07.viewEnable);
P.view.saveFigures = logical(A07.saveFigures);

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
suffix = output_suffix_local(A07.outputLabel);

P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');
P.files.lowSpeedTemplateAudit = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateAudit_20241106.mat');
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
end

function IdentificationResult = load_identification_result_local(P)
require_file_local(P.files.identificationResult, 'Step06 identification result');
loaded = load(P.files.identificationResult, 'IdentificationResult');
IdentificationResult = loaded.IdentificationResult;
end

function [LowSpeedTemplateLibrary, LowSpeedTemplateAudit, sourceInfo] = ...
        load_compatible_template_artifacts_local(P)
sourceInfo = struct('template_file', "", 'audit_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');

if exist(P.files.lowSpeedTemplateLibrary, 'file') == 2 && exist(P.files.lowSpeedTemplateAudit, 'file') == 2
    loadedLib = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
    loadedAudit = load(P.files.lowSpeedTemplateAudit, 'LowSpeedTemplateAudit');
    LowSpeedTemplateLibrary = loadedLib.LowSpeedTemplateLibrary;
    LowSpeedTemplateAudit = loadedAudit.LowSpeedTemplateAudit;
    sourceInfo.template_file = string(P.files.lowSpeedTemplateLibrary);
    sourceInfo.audit_file = string(P.files.lowSpeedTemplateAudit);
    return;
end

candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    libFile = fullfile(candidates(i).folder, candidates(i).name);
    auditFile = fullfile(candidates(i).folder, 'LowSpeedTemplateAudit_20241106.mat');
    if exist(auditFile, 'file') ~= 2
        continue;
    end
    loaded = load(libFile, 'LowSpeedTemplateLibrary');
    if ~isfield(loaded, 'LowSpeedTemplateLibrary')
        continue;
    end
    if isfield(loaded.LowSpeedTemplateLibrary, 'sensor_ids')
        available = unique(loaded.LowSpeedTemplateLibrary.sensor_ids(:).', 'stable');
    else
        continue;
    end
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing compatible Step05 template artifacts for sensors %s.', mat2str(requested));
end

libFile = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
auditFile = fullfile(candidates(bestIdx).folder, 'LowSpeedTemplateAudit_20241106.mat');
loadedLib = load(libFile, 'LowSpeedTemplateLibrary');
loadedAudit = load(auditFile, 'LowSpeedTemplateAudit');
LowSpeedTemplateLibrary = loadedLib.LowSpeedTemplateLibrary;
LowSpeedTemplateAudit = loadedAudit.LowSpeedTemplateAudit;
sourceInfo.template_file = string(libFile);
sourceInfo.audit_file = string(auditFile);
end

function [HighSpeedNumbering, sourceInfo] = load_compatible_high_speed_numbering_local(P)
sourceInfo = struct('numbering_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');

if exist(P.files.highSpeedNumbering, 'file') == 2
    loaded = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
    HighSpeedNumbering = loaded.HighSpeedNumbering;
    sourceInfo.numbering_file = string(P.files.highSpeedNumbering);
    return;
end

timeLabel = time_label_local(P.region.startTimeSec);
pattern = sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '04_high_speed_numbering', 'S*', pattern));
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
    error('Missing compatible Step04 numbering artifact for sensors %s.', mat2str(requested));
end

fileNow = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(fileNow, 'HighSpeedNumbering');
HighSpeedNumbering = loaded.HighSpeedNumbering;
keep = ismember([HighSpeedNumbering.sensor.sensor_id], requested);
HighSpeedNumbering.sensor = HighSpeedNumbering.sensor(keep);
sourceInfo.numbering_file = string(fileNow);
end

function T = build_template_audit_table_local(LowSpeedTemplateAudit, bladeId, sensorIds)
rows = repmat(struct( ...
    'BladeID', NaN, ...
    'SensorID', NaN, ...
    'SelectedRevolutionCount', NaN, ...
    'SelectedPulseCount', NaN, ...
    'WidePointCount', NaN, ...
    'TrustPointCount', NaN, ...
    'DomainLeftMM', NaN, ...
    'DomainRightMM', NaN, ...
    'DetectedCenterMM', NaN, ...
    'ThresholdV', NaN, ...
    'BaselineV', NaN), numel(sensorIds), 1);

bladeAudit = LowSpeedTemplateAudit.blade(bladeId);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    A = bladeAudit.sensor([bladeAudit.sensor.sensor_id] == sid);
    rows(i).BladeID = bladeId;
    rows(i).SensorID = sid;
    rows(i).SelectedRevolutionCount = numel(A.selected_revolution_ids);
    rows(i).SelectedPulseCount = numel(A.selected_pulse_indices);
    rows(i).WidePointCount = numel(A.x_wide);
    rows(i).TrustPointCount = numel(A.x_selected);
    rows(i).DomainLeftMM = A.x_domain_rel(1);
    rows(i).DomainRightMM = A.x_domain_rel(2);
    rows(i).DetectedCenterMM = A.xc_detected_rel;
    rows(i).ThresholdV = A.threshold;
    rows(i).BaselineV = A.baseline;
end
T = struct2table(rows);
end

function T = build_numbering_lap_table_local(HighSpeedNumbering, P, bladeId, sensorIds)
rows = [];
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    S = HighSpeedNumbering.sensor([HighSpeedNumbering.sensor.sensor_id] == sid);
    probe = load_probe_local(P, sid);
    rowIds = S.selected_rows_by_physical(:, bladeId);
    for k = 1:numel(rowIds)
        rowId = rowIds(k);
        rows = [rows; struct( ... %#ok<AGROW>
            'BladeID', bladeId, ...
            'SensorID', sid, ...
            'LapID', k, ...
            'RowID', rowId, ...
            'PeakTimeSec', probe.jilublade(rowId, 3), ...
            'StartTimeSec', probe.jilublade(rowId, 1), ...
            'EndTimeSec', probe.jilublade(rowId, 2), ...
            'PulseWidthUs', 1e6 * (probe.jilublade(rowId, 2) - probe.jilublade(rowId, 1)), ...
            'BestScore', S.best_score, ...
            'ScoreMargin', S.score_margin)];
    end
end
T = struct2table(rows);
end

function T = build_window_audit_table_local(IdentificationResult, bladeId, sensorIds, windowIds)
windowResults = IdentificationResult.WindowResult;
windowResults = windowResults([windowResults.blade_id] == bladeId);
if ~isempty(windowIds)
    keep = ismember([windowResults.window_id], windowIds);
    windowResults = windowResults(keep);
end

rows = [];
for iw = 1:numel(windowResults)
    WR = windowResults(iw);
    bundle = WR.Bundle;
    result = WR.Result;
    row = struct();
    row.BladeID = bladeId;
    row.WindowID = WR.window_id;
    row.LapStart = WR.lap_range(1);
    row.LapEnd = WR.lap_range(end);
    row.CenterTimeSec = median(bundle.T, 'omitnan');
    row.EO = result.EO_id;
    row.FrequencyHz = result.fn_id;
    row.AmplitudeMM = result.A_id;
    row.DxCMM = result.dx_c_id;
    row.WeightedVoltageRMSE = result.weighted_voltage_rmse;
    row.PlainVoltageRMSE = result.plain_voltage_rmse;
    row.PointCount = result.point_count;
    row.ValidSegmentCount = result.valid_segment_count;
    row.ChosenPass = string(WR.ChosenPass);
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        keep = bundle.S == sid;
        row.(sprintf('CH%d_PointCount', sid)) = nnz(keep);
        row.(sprintf('CH%d_XMinMM', sid)) = min_or_nan_local(bundle.X(keep));
        row.(sprintf('CH%d_XMaxMM', sid)) = max_or_nan_local(bundle.X(keep));
        row.(sprintf('CH%d_XMedianMM', sid)) = median_or_nan_local(bundle.X(keep));
        row.(sprintf('CH%d_TMinSec', sid)) = min_or_nan_local(bundle.T(keep));
        row.(sprintf('CH%d_TMaxSec', sid)) = max_or_nan_local(bundle.T(keep));
    end
    rows = [rows; row]; %#ok<AGROW>
end
T = struct2table(rows);
end

function plot_low_speed_template_audit_local(LowSpeedTemplateLibrary, LowSpeedTemplateAudit, bladeId, sensorIds)
bladeAudit = LowSpeedTemplateAudit.blade(bladeId);
fig = figure('Name', sprintf('Step07A low-speed template audit B%d', bladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 24, 8 * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorIds), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    A = bladeAudit.sensor([bladeAudit.sensor.sensor_id] == sid);
    Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid);
    [xWide, vWide] = downsample_cloud_local(A.x_wide, A.v_wide, 18000);
    [xSel, vSel] = downsample_cloud_local(A.x_selected, A.v_selected, 18000);

    nexttile;
    plot(xWide, vWide, '.', 'Color', [0.78 0.78 0.78], 'MarkerSize', 3, 'DisplayName', 'wide cloud'); hold on;
    plot(xSel, vSel, 'k.', 'MarkerSize', 3, 'DisplayName', 'trust cloud');
    plot(Tpl.x_grid, Tpl.v_grid, 'r-', 'LineWidth', 1.3, 'DisplayName', 'template');
    xline(A.x_domain_rel(1), 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(A.x_domain_rel(2), 'b-.', 'LineWidth', 1.0, 'DisplayName', 'domain');
    xline(0, 'k--', 'LineWidth', 0.9, 'DisplayName', 'final center');
    xlabel('x relative to final center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('B%d CH%d low-speed template cloud', bladeId, sid), 'FontWeight', 'normal');
    legend('Location', 'best', 'Box', 'off');
    grid on; box on; style_axes_local();

    nexttile;
    yyaxis left;
    plot(A.x_profile, A.v_profile, 'k-', 'LineWidth', 1.1);
    ylabel('Smoothed profile (V)');
    yyaxis right;
    plot(A.x_profile, A.gradient_abs, 'b-', 'LineWidth', 1.1);
    hold on;
    yline(A.gradient_threshold, 'r--', 'LineWidth', 1.0);
    ylabel('|dV/dx| (V/mm)');
    xlabel('x relative to final center (mm)');
    title(sprintf('B%d CH%d gradient-based trust decision', bladeId, sid), 'FontWeight', 'normal');
    grid on; box on; style_axes_local();
end
end

function plot_high_speed_numbering_audit_local(numberingTable, bladeId, sensorIds, maxLapsToShow)
fig = figure('Name', sprintf('Step07A high-speed numbering audit B%d', bladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 23, 14], 'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    T = numberingTable(numberingTable.SensorID == sid, :);
    keep = T.LapID <= maxLapsToShow;
    plot(T.LapID(keep), T.RowID(keep), 'o-', 'LineWidth', 1.1, 'DisplayName', sprintf('CH%d row id', sid));
end
xlabel('Lap ID');
ylabel('Selected row ID');
title(sprintf('B%d high-speed numbering: selected row ids', bladeId), 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on; style_axes_local();

nexttile;
hold on;
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    T = numberingTable(numberingTable.SensorID == sid, :);
    keep = T.LapID <= maxLapsToShow;
    plot(T.LapID(keep), 1e3 * T.PeakTimeSec(keep), 'o-', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('CH%d peak time', sid));
end
xlabel('Lap ID');
ylabel('Peak time (ms)');
title(sprintf('B%d high-speed numbering: selected peak times', bladeId), 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on; style_axes_local();
end

function plot_identification_window_audit_local(windowTable, bladeId, sensorIds)
fig = figure('Name', sprintf('Step07A Step06 window audit B%d', bladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 24, 16], 'NumberTitle', 'off');
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(windowTable.WindowID, windowTable.FrequencyHz, 'o-r', 'LineWidth', 1.2, ...
    'MarkerFaceColor', [0.9 0.1 0.1], 'DisplayName', 'frequency'); hold on;
yyaxis right;
stem(windowTable.WindowID, windowTable.AmplitudeMM, 'Color', [0.2 0.55 0.25], ...
    'LineWidth', 1.0, 'Marker', 'none', 'DisplayName', 'amplitude');
yyaxis left;
for k = 1:height(windowTable)
    text(windowTable.WindowID(k), windowTable.FrequencyHz(k), sprintf(' EO%d', windowTable.EO(k)), ...
        'FontName', 'Times New Roman', 'FontSize', 7.5, 'Color', [0.2 0.2 0.2]);
end
xlabel('Window ID');
ylabel('Frequency (Hz)');
yyaxis right;
ylabel('Amplitude (mm)');
title(sprintf('B%d Step06 window results', bladeId), 'FontWeight', 'normal');
grid on; box on; style_axes_local();

nexttile;
hold on;
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    plot(windowTable.WindowID, windowTable.(sprintf('CH%d_PointCount', sid)), 'o-', ...
        'LineWidth', 1.1, 'DisplayName', sprintf('CH%d points', sid));
end
xlabel('Window ID');
ylabel('Selected point count');
title('Per-sensor point contribution by window', 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on; style_axes_local();

nexttile;
plot(windowTable.WindowID, windowTable.WeightedVoltageRMSE, 'o-k', ...
    'LineWidth', 1.2, 'MarkerFaceColor', [0.2 0.2 0.2]);
xlabel('Window ID');
ylabel('Weighted voltage RMSE (V)');
title('Step06 fit error by window', 'FontWeight', 'normal');
grid on; box on; style_axes_local();
end

function plot_step06_bundle_domain_audit_local(IdentificationResult, bladeId, sensorIds, windowIds)
windowResults = IdentificationResult.WindowResult;
windowResults = windowResults([windowResults.blade_id] == bladeId);
if ~isempty(windowIds)
    keep = ismember([windowResults.window_id], windowIds);
    windowResults = windowResults(keep);
end

for iw = 1:numel(windowResults)
    WR = windowResults(iw);
    if isempty(WR.Bundle) || isempty(WR.Result)
        continue;
    end
    bundle = WR.Bundle;
    fig = figure('Name', sprintf('Step07A bundle domain audit B%d W%02d', bladeId, WR.window_id), ...
        'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 4.8 * numel(sensorIds)], ...
        'NumberTitle', 'off');
    tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(sensorIds)
        sid = sensorIds(i);
        rows = find(bundle.S == sid);
        xLo = bundle.x_domain_by_sensor(i, 1);
        xHi = bundle.x_domain_by_sensor(i, 2);
        xGrid = linspace(xLo, xHi, 1000).';
        vTpl = bundle.interp_v{i}(xGrid);
        nexttile;
        plot(xGrid, vTpl, 'k-', 'LineWidth', 1.2, 'DisplayName', 'template'); hold on;
        scatter(bundle.X(rows), bundle.V(rows), 10, bundle.T(rows), 'filled', 'DisplayName', 'Step06 selected points');
        xline(xLo, 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        xline(xHi, 'b-.', 'LineWidth', 1.0, 'DisplayName', 'domain');
        xlabel('x query (mm)');
        ylabel(sprintf('CH%d voltage (V)', sid));
        title(sprintf('B%d W%02d CH%d selected query points', bladeId, WR.window_id, sid), ...
            'FontWeight', 'normal');
        if i == 1
            legend('Location', 'best', 'Box', 'off');
        end
        grid on; box on; style_axes_local();
    end
end
end

function probe = load_probe_local(P, sid)
probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
require_file_local(probeFile, sprintf('high-speed CH%d pulse timing file', sid));
loaded = load(probeFile, 'jilublade');
probe = loaded;
end

function Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid)
entry = LowSpeedTemplateLibrary.entry(bladeId);
if ~entry.exists
    error('Missing low-speed template bundle for B%d.', bladeId);
end
loadedBundle = load(char(entry.template_file), 'TemplateBundle');
bundle = loadedBundle.TemplateBundle;
row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
if height(row) ~= 1
    error('Low-speed template bundle for B%d does not contain CH%d.', bladeId, sid);
end
loadedSensor = load(char(row.template_file(1)), 'sensor_template');
Tpl = loadedSensor.sensor_template.Sensor;
end

function [xDs, vDs] = downsample_cloud_local(x, v, maxPoints)
keep = isfinite(x) & isfinite(v);
x = x(keep);
v = v(keep);
if numel(x) > maxPoints
    idx = round(linspace(1, numel(x), maxPoints));
else
    idx = 1:numel(x);
end
xDs = x(idx);
vDs = v(idx);
end

function y = min_or_nan_local(x)
if isempty(x) || all(~isfinite(x))
    y = NaN;
else
    y = min(x, [], 'omitnan');
end
end

function y = max_or_nan_local(x)
if isempty(x) || all(~isfinite(x))
    y = NaN;
else
    y = max(x, [], 'omitnan');
end
end

function y = median_or_nan_local(x)
if isempty(x) || all(~isfinite(x))
    y = NaN;
else
    y = median(x, 'omitnan');
end
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
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

function suffix = output_suffix_local(label)
suffix = '';
if isempty(label)
    return;
end
label = regexprep(char(label), '[^\w\d-]', '_');
label = regexprep(label, '_+', '_');
label = strtrim(label);
if ~isempty(label)
    suffix = ['_', label];
end
end

