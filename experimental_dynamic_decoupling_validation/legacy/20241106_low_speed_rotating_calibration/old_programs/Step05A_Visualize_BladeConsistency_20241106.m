%% Step05A_Visualize_BladeConsistency_20241106
% Visualize whether the low-speed template and high-speed selected waveforms
% from different sensors are consistent with the same blade.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));

%% Parameters
P.targetBlade = 6;
P.analysisSensors = [2 3 5 7];
P.analysisStartTime = 75.0;
P.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', make_time_label_local(P.analysisStartTime), P.templateSuffix);
P.resultSuffix = sprintf('Direct_%s', make_time_label_local(P.analysisStartTime));
P.windowID = [];
P.maxPassCount = 3;
P.useBestWindowFromResult = true;
P.saveFigures = true;

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
templateBundleFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.templateSuffix));
dynamicMapFile = fullfile(routeDir, 'output', 'dynamic_maps', ...
    sprintf('DynamicMap_B%d_%s_SlidingWindows_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.dynamicSuffix));
resultFile = fullfile(routeDir, 'output', 'identification', ...
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.resultSuffix));

figureDir = fullfile(routeDir, 'output', 'figures');
if exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

fprintf('\n=== Step05A: blade-consistency visualization ===\n');
fprintf('Blade: %d, sensors: %s\n', P.targetBlade, mat2str(P.analysisSensors));
fprintf('Template bundle: %s\n', templateBundleFile);
fprintf('Dynamic map: %s\n', dynamicMapFile);
fprintf('Result file: %s\n', resultFile);

Template = load_template_from_bundle_local(templateBundleFile, P.analysisSensors);
loadedMap = load(dynamicMapFile, 'DynamicMap');
DynamicMap = filter_dynamic_map_sensors_local(loadedMap.DynamicMap, P.analysisSensors);
loadedResult = load(resultFile, 'Result');
Result = loadedResult.Result;

windowID = choose_window_id_local(Result, P);
fprintf('Window used for consistency audit: %d\n', windowID);

Wmap = DynamicMap.Window(windowID);
sensorIds = P.analysisSensors(:).';
passData = build_pass_data_local(Wmap, Template, sensorIds);
pairTable = build_pair_consistency_table_local(passData, sensorIds);
disp(pairTable);

plot_low_vs_high_overlay_local(passData, Result, windowID, figureDir, P);
plot_peak_consistency_local(passData, Result, windowID, figureDir, P);
writetable(pairTable, fullfile(figureDir, ...
    sprintf('Step05A_BladeConsistency_Table_B%d_%s_W%02d_20241106.csv', ...
    P.targetBlade, sensorTag, windowID)));

function windowID = choose_window_id_local(Result, P)
if ~isempty(P.windowID)
    windowID = P.windowID;
    return;
end
if P.useBestWindowFromResult && isfield(Result, 'BestWindow') && isstruct(Result.BestWindow) ...
        && isfield(Result.BestWindow, 'window_id') && isfinite(Result.BestWindow.window_id)
    windowID = Result.BestWindow.window_id;
    return;
end
windowID = 1;
end

function passData = build_pass_data_local(Wmap, Template, sensorIds)
passData = repmat(struct( ...
    'sensor_id', [], ...
    'template_x', [], ...
    'template_v', [], ...
    'x_domain', [], ...
    'passes', []), numel(sensorIds), 1);

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rawSensor = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
    tplSensor = Template.Sensor([Template.Sensor.sensor_id] == sid);
    if isempty(rawSensor) || isempty(tplSensor)
        error('Missing CH%d in window/template.', sid);
    end

    segments = split_sensor_points_into_passes_local(rawSensor.t);
    passStruct = repmat(struct( ...
        'pass_id', [], ...
        't', [], ...
        'x_rel', [], ...
        'V', [], ...
        'W', [], ...
        'peak_x_mm', nan, ...
        'peak_v', nan, ...
        'corr_with_template', nan, ...
        'rmse_with_template', nan), numel(segments), 1);

    for ip = 1:numel(segments)
        rows = segments{ip};
        x_rel = rawSensor.x_rel(rows);
        v_obs = rawSensor.V(rows);
        t_obs = rawSensor.t(rows);
        w_obs = rawSensor.W(rows);
        [x_rel, order] = sort(x_rel(:));
        v_obs = v_obs(order);
        t_obs = t_obs(order);
        w_obs = w_obs(order);

        v_tpl = interp1(tplSensor.x_grid(:), tplSensor.v_grid(:), x_rel, 'linear', 'extrap');
        valid = isfinite(v_obs) & isfinite(v_tpl);
        corr_val = local_corrcoef_scalar(v_obs(valid), v_tpl(valid));
        rmse_val = sqrt(mean((v_obs(valid) - v_tpl(valid)).^2, 'omitnan'));
        [peak_v, peak_idx] = max(v_obs);

        passStruct(ip).pass_id = ip;
        passStruct(ip).t = t_obs;
        passStruct(ip).x_rel = x_rel;
        passStruct(ip).V = v_obs;
        passStruct(ip).W = w_obs;
        passStruct(ip).peak_x_mm = x_rel(peak_idx);
        passStruct(ip).peak_v = peak_v;
        passStruct(ip).corr_with_template = corr_val;
        passStruct(ip).rmse_with_template = rmse_val;
    end

    passData(is).sensor_id = sid;
    passData(is).template_x = tplSensor.x_grid(:);
    passData(is).template_v = tplSensor.v_grid(:);
    passData(is).x_domain = tplSensor.x_domain(:).';
    passData(is).passes = passStruct;
end
end

function pairTable = build_pair_consistency_table_local(passData, sensorIds)
rows = repmat(struct( ...
    'sensorA', nan, ...
    'sensorB', nan, ...
    'passCountCompared', nan, ...
    'meanPeakXDiffMM', nan, ...
    'maxPeakXDiffMM', nan, ...
    'meanTemplateCorr', nan, ...
    'minTemplateCorr', nan), 0, 1);

for ia = 1:(numel(sensorIds) - 1)
    for ib = (ia + 1):numel(sensorIds)
        A = passData(ia).passes;
        B = passData(ib).passes;
        n = min(numel(A), numel(B));
        peakDiff = nan(n, 1);
        corrBoth = nan(n, 1);
        for ip = 1:n
            peakDiff(ip) = abs(A(ip).peak_x_mm - B(ip).peak_x_mm);
            corrBoth(ip) = min(A(ip).corr_with_template, B(ip).corr_with_template);
        end
        row = struct();
        row.sensorA = sensorIds(ia);
        row.sensorB = sensorIds(ib);
        row.passCountCompared = n;
        row.meanPeakXDiffMM = mean(peakDiff, 'omitnan');
        row.maxPeakXDiffMM = max(peakDiff, [], 'omitnan');
        row.meanTemplateCorr = mean(corrBoth, 'omitnan');
        row.minTemplateCorr = min(corrBoth, [], 'omitnan');
        rows(end + 1, 1) = row; %#ok<AGROW>
    end
end

pairTable = struct2table(rows);
end

function plot_low_vs_high_overlay_local(passData, Result, windowID, figureDir, P)
sensorIds = [passData.sensor_id];
maxPassCount = min(P.maxPassCount, min(arrayfun(@(s) numel(s.passes), passData)));
fig = figure('Name', sprintf('20241106 blade-consistency overlay B%d W%02d', P.targetBlade, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 18]);
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    S = passData(is);
    nexttile;
    plot(S.template_x, S.template_v, 'k-', 'LineWidth', 1.4, 'DisplayName', 'low-speed template'); hold on;
    cmap = lines(maxPassCount);
    for ip = 1:maxPassCount
        Pk = S.passes(ip);
        plot(Pk.x_rel, Pk.V, '-', 'Color', cmap(ip, :), 'LineWidth', 1.0, ...
            'DisplayName', sprintf('high-speed pass %d', ip));
        scatter(Pk.peak_x_mm, Pk.peak_v, 18, cmap(ip, :), 'filled', ...
            'HandleVisibility', 'off');
    end
    xline(S.x_domain(1), ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    xline(S.x_domain(2), ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    ylabel(sprintf('CH%d V (V)', S.sensor_id));
    title(sprintf('B%d W%02d CH%d: first %d passes, template overlap', ...
        P.targetBlade, windowID, S.sensor_id, maxPassCount), 'FontWeight', 'normal');
    if is == numel(sensorIds)
        xlabel('x relative to template center (mm)');
    end
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();
end

sgtitle(sprintf('Low-speed template vs selected high-speed waveform passes, EO=%d, f=%.3f Hz, A=%.4f mm', ...
    Result.BestWindow.EO_id, Result.BestWindow.fn_id, Result.BestWindow.A_id), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05A_LowVsHighOverlay_B%d_S%s_W%02d_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05A_LowVsHighOverlay_B%d_S%s_W%02d_20241106.pdf', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'ContentType', 'vector');
end
end

function plot_peak_consistency_local(passData, Result, windowID, figureDir, P)
sensorIds = [passData.sensor_id];
maxPassCount = min(P.maxPassCount, min(arrayfun(@(s) numel(s.passes), passData)));
fig = figure('Name', sprintf('20241106 blade-consistency peak audit B%d W%02d', P.targetBlade, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

peakX = nan(maxPassCount, numel(sensorIds));
corrMat = nan(maxPassCount, numel(sensorIds));
for is = 1:numel(sensorIds)
    for ip = 1:maxPassCount
        peakX(ip, is) = passData(is).passes(ip).peak_x_mm;
        corrMat(ip, is) = passData(is).passes(ip).corr_with_template;
    end
end

nexttile;
plot(1:maxPassCount, peakX, 'LineWidth', 1.1, 'Marker', 'o');
xlabel('High-speed pass ID');
ylabel('Peak position (mm)');
title(sprintf('B%d W%02d peak position across sensors', P.targetBlade, windowID), ...
    'FontWeight', 'normal');
legend(compose('CH%d', sensorIds), 'Location', 'best', 'Box', 'off');
style_axes_local();

nexttile;
bar(corrMat, 'grouped');
xlabel('High-speed pass ID');
ylabel('Corr(template, pass)');
title(sprintf('B%d W%02d template correlation across sensors', P.targetBlade, windowID), ...
    'FontWeight', 'normal');
legend(compose('CH%d', sensorIds), 'Location', 'best', 'Box', 'off');
ylim([0, 1.05]);
style_axes_local();

sgtitle(sprintf('Blade-consistency summary, EO=%d, f=%.3f Hz, A=%.4f mm', ...
    Result.BestWindow.EO_id, Result.BestWindow.fn_id, Result.BestWindow.A_id), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05A_PeakConsistency_B%d_S%s_W%02d_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors), windowID)), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05A_PeakConsistency_B%d_S%s_W%02d_20241106.pdf', ...
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
    if isempty(row)
        error('Bundle %s does not contain CH%d.', bundleFile, sid);
    end
    loadedSensor = load(char(row.template_file(1)), 'sensor_template');
    sensorEntries{i} = loadedSensor.sensor_template.Sensor;
end

Template = struct();
Template.Route = bundle.Route;
Template.TargetBlade = bundle.TargetBlade;
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = ['S', sprintf('%d', sensorIds)];
Template.Sensor = vertcat(sensorEntries{:});
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, sensorIds)
DynamicMap.SensorIDs = sensorIds(:).';
DynamicMap.SensorTag = ['S', sprintf('%d', sensorIds)];
for iw = 1:numel(DynamicMap.Window)
    keepMask = ismember([DynamicMap.Window(iw).Sensor.sensor_id], sensorIds);
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(keepMask);
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

function value = local_corrcoef_scalar(a, b)
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
