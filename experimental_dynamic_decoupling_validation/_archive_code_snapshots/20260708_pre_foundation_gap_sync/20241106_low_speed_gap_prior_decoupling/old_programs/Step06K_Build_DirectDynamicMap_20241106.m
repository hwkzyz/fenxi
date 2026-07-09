%% Step06K: build direct DynamicMap for 20241106
% Direct-only dynamic waveform map matching the 20250527/20251222 Step02
% structure:
%   DynamicMap.Window(w).Sensor(s).t / x_abs / x_rel / V / W / theta
%
% This script keeps continuous raw waveform samples in each selected blade
% pass. It does not use the gap-library highMap and does not downsample each
% pulse to a small point cloud.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
outDir = fullfile(routeDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06k_20241106_direct_dynamic_map');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Parameters to tune
P = struct();
P.targetBlade = 1;              % blade slot relative to the 20241106 OPR pulse
P.analysisSensors = [5 7];      % capacitive probes used by the direct route
P.analysisStartTimeSec = 72.0;  % visible start of the high-speed waveform region
P.targetLaps = 20;              % number of selected blade passes per probe
P.windowLaps = 3;               % laps per sliding identification window
P.slidingStepLaps = 1;
P.pulseHalfWindowSec = 3.0e-4;  % raw waveform window around each blade pass
P.edgeFractionForBaseline = 0.18;
P.minFinitePointsPerPulse = 100;
P.blockStepSamples = 1e7;
P.blockLabelStep = 1000;
P.rawCacheMaxEntries = 10;
P.templateFile = '';

P.targetBlade = parse_scalar_env_local('STEP06K_TARGET_BLADE', P.targetBlade);
P.analysisSensors = parse_int_env_local('STEP06K_ANALYSIS_SENSORS', P.analysisSensors);
P.analysisStartTimeSec = parse_scalar_env_local('STEP06K_START_TIME_SEC', P.analysisStartTimeSec);
P.targetLaps = round(parse_scalar_env_local('STEP06K_TARGET_LAPS', P.targetLaps));
P.windowLaps = round(parse_scalar_env_local('STEP06K_WINDOW_LAPS', P.windowLaps));
P.slidingStepLaps = round(parse_scalar_env_local('STEP06K_SLIDING_STEP_LAPS', P.slidingStepLaps));
P.pulseHalfWindowSec = parse_scalar_env_local('STEP06K_PULSE_HALF_WINDOW_SEC', P.pulseHalfWindowSec);

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
if isempty(P.templateFile)
    P.templateFile = fullfile(outDir, sprintf( ...
        'Step06I_Generated_LowSpeed_Template_20241106_B%d_%s.mat', ...
        P.targetBlade, sensorTag));
end

step02File = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
if ~isfile(step02File)
    error('Run Step02 first. Missing file: %s', step02File);
end
if ~isfile(P.templateFile)
    error('Run Step06I first. Missing low-speed template: %s', P.templateFile);
end

S02 = load(step02File);
S = load(P.templateFile, 'Template');
Template = clean_template_local(S.Template, P.analysisSensors);

fprintf('\n=== Step06K 20241106 direct DynamicMap ===\n');
fprintf('Low-speed template: %s\n', P.templateFile);
fprintf('High-speed raw folder: %s\n', S02.highDir);
fprintf('Blade B%d, sensors %s\n', P.targetBlade, mat2str(P.analysisSensors));
fprintf('Start time %.6f s, target laps %d, window %d laps, step %d lap(s)\n', ...
    P.analysisStartTimeSec, P.targetLaps, P.windowLaps, P.slidingStepLaps);
fprintf('Raw pulse half-window %.1f us\n', P.pulseHalfWindowSec * 1e6);

opr = load_opr_local(S02.highDir, S02.sampleRateHz);
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
        error('Blade slot %d is not available for sensor %d.', P.targetBlade, sid);
    end

    slot = S02.highDisp(sensorIdx).slot(P.targetBlade);
    eventTimesAll = slot.time(:);
    keepEvent = eventTimesAll >= P.analysisStartTimeSec & isfinite(eventTimesAll);
    eventTimes = eventTimesAll(keepEvent);
    if numel(eventTimes) < P.targetLaps
        error('P%d has only %d blade passes after %.6f s.', ...
            sid, numel(eventTimes), P.analysisStartTimeSec);
    end
    eventTimes = eventTimes(1:P.targetLaps);
    selectionRows{is} = eventTimes(:);
    globalTime(1) = min(globalTime(1), min(eventTimes) - P.pulseHalfWindowSec);
    globalTime(2) = max(globalTime(2), max(eventTimes) + P.pulseHalfWindowSec);

    Lap = repmat(struct('lap_id', NaN, 'event_time', NaN, 't', [], ...
        'x_abs', [], 'V', [], 'W', [], 'theta', [], 'rot_freq_hz', NaN), ...
        P.targetLaps, 1);

    for lapId = 1:P.targetLaps
        tEvent = eventTimes(lapId);
        rotFreqHz = local_rot_freq_local(opr.tCenter, tEvent);
        tipSpeedMmS = 2 * pi * S02.rTipMm * rotFreqHz;

        [tRaw, vRaw, ok] = read_raw_segment_local(S02.highDir, sid, tEvent, ...
            P.pulseHalfWindowSec, S02.sampleRateHz, P.blockStepSamples, ...
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

numWindows = floor((P.targetLaps - P.windowLaps) / P.slidingStepLaps) + 1;
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
DynamicMap.method = 'direct_continuous_raw_waveform_window_map';
DynamicMap.description = ['Continuous raw high-speed waveform windows for the direct ' ...
    'low-speed-template method, matching the 20250527/20251222 DynamicMap layout.'];
DynamicMap.TemplateFile = P.templateFile;
DynamicMap.Step02File = step02File;
DynamicMap.lowDir = S02.lowDir;
DynamicMap.highDir = S02.highDir;
DynamicMap.SourceSettings = P;
DynamicMap.analysisSensors = P.analysisSensors(:).';
DynamicMap.targetBlade = P.targetBlade;
DynamicMap.globalTimeWindow = globalTime;
DynamicMap.selectionEventTimesBySensor = selectionRows;
DynamicMap.Window = Window;

matFile = fullfile(outDir, sprintf( ...
    'Step06K_DirectDynamicMap_20241106_B%d_%s.mat', P.targetBlade, sensorTag));
save(matFile, 'DynamicMap', 'Template', 'P', '-v7.3');

summary = build_window_summary_local(DynamicMap);
csvFile = fullfile(outDir, sprintf( ...
    'Step06K_DirectDynamicMap_WindowSummary_20241106_B%d_%s.csv', ...
    P.targetBlade, sensorTag));
writetable(summary, csvFile);

figFile = fullfile(figDir, sprintf( ...
    'Step06K_DirectDynamicMap_20241106_B%d_%s.png', P.targetBlade, sensorTag));
plot_dynamic_map_local(DynamicMap, Template, figFile);

fprintf('\nStep06K complete.\n');
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

function opr = load_opr_local(caseDir, sampleRateHz)
oprFile = fullfile(caseDir, 'jiluOPR.mat');
if ~isfile(oprFile)
    error('Missing OPR file: %s', oprFile);
end
S = load(oprFile, 'jiluOPR');
raw = S.jiluOPR;
raw(raw(:, 1) == 0, :) = [];
t = raw(:, 1) / sampleRateHz;
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

function [tRaw, vRaw, ok] = read_raw_segment_local(caseDir, sensorId, tCenter, ...
    halfWindowSec, sampleRateHz, blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries)
sampleStart = floor((tCenter - halfWindowSec) * sampleRateHz);
sampleEnd = ceil((tCenter + halfWindowSec) * sampleRateHz);
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
keep = rawT >= tCenter - halfWindowSec & rawT <= tCenter + halfWindowSec & ...
    isfinite(rawAll(:, 2));
tRaw = rawT(keep);
vRaw = rawAll(keep, 2);
ok = numel(tRaw) >= 10;
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
