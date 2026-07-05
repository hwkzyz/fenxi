%% Step07C_VisualizeLowVsHighWaveformAlignment_20241106
% Put the low-speed template and high-speed waveform data on the same axes.
% This script is only for checking physical alignment, not for rerunning
% identification.
%
% Key idea:
%   low-speed template and high-speed waveform are both mapped into the same
%   OPR-center-based x coordinate before comparison.
%
%   thetaStd = standard_angles_opr_center(sensor, blade)
%   x_abs    = wrap(theta_point - thetaStd) * r_tip
%   x_rel    = x_abs - x_c_template
%
% Therefore the sensor installation angle is already embedded in thetaStd.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local parameter block. Edit here directly.
S08C = struct();
S08C.analysisSensors = [5 7];
S08C.startTimeSec = 75.0;
S08C.bladeId = 6;
S08C.lapIds = 1:6;                 % [] -> first up to maxLapsAuto
S08C.maxLapsAuto = 6;
S08C.showXAbs = true;              % also show pre-center-shift alignment
S08C.showXRel = true;              % show template vs high-speed overlay
S08C.saveFigures = false;

P = apply_local_options_local(P, S08C);
[LowSpeedTemplateLibrary, templateSource] = load_compatible_template_library_local(P);
[WaveformLibrary, waveformSource] = load_compatible_waveform_library_local(P);

bladeId = S08C.bladeId;
sensorIds = P.sensors.analysis(:).';

fprintf('\n=== Step07C: low-speed template vs high-speed waveform ===\n');
fprintf('Blade: B%d\n', bladeId);
fprintf('Sensors: %s\n', mat2str(sensorIds));
fprintf('Low-speed template source: %s\n', templateSource.template_file);
fprintf('High-speed waveform source: %s\n', waveformSource.waveform_file);
fprintf(['Installation angle usage:\n' ...
    '  thetaStd = standard_angles_opr_center(sensor, blade)\n' ...
    '  x_abs    = wrap(theta_point - thetaStd) * r_tip\n' ...
    '  x_rel    = x_abs - x_c_template\n']);

if S08C.showXRel
    plot_xrel_overlay_local(LowSpeedTemplateLibrary, WaveformLibrary, bladeId, sensorIds, S08C, P);
end
if S08C.showXAbs
    plot_xabs_alignment_local(LowSpeedTemplateLibrary, WaveformLibrary, bladeId, sensorIds, S08C, P);
end

function P = apply_local_options_local(P, S08C)
P.sensors.analysis = S08C.analysisSensors(:).';
P.region.startTimeSec = S08C.startTimeSec;
P.view.saveFigures = logical(S08C.saveFigures);

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');
P.files.waveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('WaveformLibrary_%s_20241106.mat', timeLabel));
end

function [LowSpeedTemplateLibrary, sourceInfo] = load_compatible_template_library_local(P)
sourceInfo = struct('template_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');
if exist(P.files.lowSpeedTemplateLibrary, 'file') == 2
    loaded = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
    LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
    sourceInfo.template_file = string(P.files.lowSpeedTemplateLibrary);
    return;
end

candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'LowSpeedTemplateLibrary');
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
    error('Missing compatible low-speed template library for sensors %s.', mat2str(requested));
end
fileNow = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(fileNow, 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
sourceInfo.template_file = string(fileNow);
end

function [WaveformLibrary, sourceInfo] = load_compatible_waveform_library_local(P)
sourceInfo = struct('waveform_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');
if exist(P.files.waveformLibrary, 'file') == 2
    loaded = load(P.files.waveformLibrary, 'WaveformLibrary');
    WaveformLibrary = loaded.WaveformLibrary;
    sourceInfo.waveform_file = string(P.files.waveformLibrary);
    return;
end

timeLabel = time_label_local(P.region.startTimeSec);
pattern = sprintf('WaveformLibrary_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '05_waveform_library', 'S*', pattern));
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'WaveformLibrary');
    if ~isfield(loaded, 'WaveformLibrary') || ~isfield(loaded.WaveformLibrary, 'analysis_sensors')
        continue;
    end
    available = unique(loaded.WaveformLibrary.analysis_sensors(:).', 'stable');
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing compatible high-speed waveform library for sensors %s.', mat2str(requested));
end
fileNow = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(fileNow, 'WaveformLibrary');
WaveformLibrary = loaded.WaveformLibrary;
sourceInfo.waveform_file = string(fileNow);
end

function plot_xrel_overlay_local(LowSpeedTemplateLibrary, WaveformLibrary, bladeId, sensorIds, S08C, P)
fig = figure('Name', sprintf('Step07C x_rel overlay B%d', bladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 4.8 * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid);
    Hs = get_waveform_sensor_local(WaveformLibrary, bladeId, sid);
    lapIds = resolve_lap_ids_local(Hs, S08C);

    nexttile;
    plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.4, 'DisplayName', 'low-speed template'); hold on;
    cmap = lines(max(1, numel(lapIds)));
    for k = 1:numel(lapIds)
        lapId = lapIds(k);
        if lapId > numel(Hs.Lap)
            continue;
        end
        Lap = Hs.Lap(lapId);
        [xRel, order] = sort(Lap.x_rel(:));
        vNow = Lap.V(order);
        plot(xRel, vNow, '-', 'Color', cmap(k, :), 'LineWidth', 0.95, ...
            'DisplayName', sprintf('high-speed lap %d', lapId));
    end
    xline(Tpl.x_domain(1), 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(Tpl.x_domain(2), 'b-.', 'LineWidth', 1.0, 'DisplayName', 'template domain');
    xline(0, 'k--', 'LineWidth', 0.9, 'DisplayName', 'template center');
    xlabel('x relative to template center (mm)');
    ylabel(sprintf('CH%d voltage (V)', sid));
    title(sprintf(['B%d CH%d: low-speed template vs high-speed waveform\n' ...
        'thetaStd = %.3f deg, x_c = %.4f mm'], ...
        bladeId, sid, Hs.theta_std_deg, Tpl.xc), 'FontWeight', 'normal');
    if i == 1
        legend('Location', 'best', 'Box', 'off');
    end
    grid on; box on; style_axes_local();
end

if P.view.saveFigures
    figDir = fullfile(P.view.figureDir, '08_alignment_audit');
    if exist(figDir, 'dir') ~= 7
        mkdir(figDir);
    end
    exportgraphics(fig, fullfile(figDir, sprintf('Step07C_XRelOverlay_B%d_%s.png', ...
        bladeId, sensor_tag_local(sensorIds))), 'Resolution', 300);
end
end

function plot_xabs_alignment_local(LowSpeedTemplateLibrary, WaveformLibrary, bladeId, sensorIds, S08C, P)
fig = figure('Name', sprintf('Step07C x_abs alignment B%d', bladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 4.8 * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid);
    Hs = get_waveform_sensor_local(WaveformLibrary, bladeId, sid);
    lapIds = resolve_lap_ids_local(Hs, S08C);

    nexttile;
    cmap = lines(max(1, numel(lapIds)));
    for k = 1:numel(lapIds)
        lapId = lapIds(k);
        if lapId > numel(Hs.Lap)
            continue;
        end
        Lap = Hs.Lap(lapId);
        [xAbs, order] = sort(Lap.x_abs(:));
        vNow = Lap.V(order);
        plot(xAbs, vNow, '-', 'Color', cmap(k, :), 'LineWidth', 0.95, ...
            'DisplayName', sprintf('lap %d', lapId)); hold on;
    end
    xline(Tpl.xc, 'k--', 'LineWidth', 1.0, 'DisplayName', 'template x_c');
    xlabel('x absolute under OPR-center reference (mm)');
    ylabel(sprintf('CH%d voltage (V)', sid));
    title(sprintf(['B%d CH%d: x_{abs} before subtracting template center\n' ...
        'x_rel = x_abs - x_c,  thetaStd = %.3f deg'], ...
        bladeId, sid, Hs.theta_std_deg), 'FontWeight', 'normal');
    if i == 1
        legend('Location', 'best', 'Box', 'off');
    end
    grid on; box on; style_axes_local();
end

if P.view.saveFigures
    figDir = fullfile(P.view.figureDir, '08_alignment_audit');
    if exist(figDir, 'dir') ~= 7
        mkdir(figDir);
    end
    exportgraphics(fig, fullfile(figDir, sprintf('Step07C_XAbsAlignment_B%d_%s.png', ...
        bladeId, sensor_tag_local(sensorIds))), 'Resolution', 300);
end
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

function Hs = get_waveform_sensor_local(WaveformLibrary, bladeId, sid)
S = WaveformLibrary.Blade(bladeId).Sensor;
idx = find([S.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('High-speed waveform library for B%d does not contain CH%d.', bladeId, sid);
end
Hs = S(idx);
end

function lapIds = resolve_lap_ids_local(Hs, S08C)
if ~isempty(S08C.lapIds)
    lapIds = S08C.lapIds(:).';
else
    lapIds = 1:min(S08C.maxLapsAuto, numel(Hs.Lap));
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

