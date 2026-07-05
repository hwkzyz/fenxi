%% Step06V_VisualizeIdentificationWaveformComparison_20241106
% Visualize Step06 identification waveform comparison.
% Focus:
%   1) observed selected high-speed points vs identified predicted waveform,
%   2) after displacement compensation, whether high-speed points fall back
%      onto the low-speed template,
%   3) for one pass (one circle), whether different sensors show consistent
%      blade-wise waveform alignment.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Step06V local parameter block.
% Edit here directly when auditing one blade / one sensor group.
S07V = struct();
S07V.analysisSensors = [5 7];
S07V.startTimeSec = 75.0;
S07V.outputLabel = 'B4_only';
S07V.blades = 4;
S07V.windowIds = [];                 % [] -> use all sliding windows for each blade
S07V.useBestWindowWhenEmpty = false;
S07V.passIds = [];                   % [] -> use center pass only
S07V.maxPassCountWhenAuto = 1;
S07V.saveFigures = true;
S07V.showTimeDomain = true;
S07V.showQueryDomain = true;
S07V.showStackedPass = true;
S07V.showGapComparison = true;
S07V.showResidual = true;

P = apply_step07v_local_options_local(P, S07V);
require_file_local(P.files.identificationResult, 'Step06 identification result');

loaded = load(P.files.identificationResult, 'IdentificationResult');
IdentificationResult = loaded.IdentificationResult;

figDir = fullfile(P.view.figureDir, '06_identification');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

windowResults = IdentificationResult.WindowResult;
bestByBlade = IdentificationResult.BestByBlade;

fprintf('\n=== Step06V: identification waveform comparison ===\n');
fprintf('Identification source: %s\n', P.files.identificationResult);
fprintf('Sensors: %s\n', mat2str(P.sensors.analysis));
fprintf('Blades: %s\n', mat2str(P.step07v.blades));

for bladeId = P.step07v.blades
    bladeWindows = windowResults([windowResults.blade_id] == bladeId);
    if isempty(bladeWindows)
        warning('IdentificationResult does not contain blade B%d.', bladeId);
        continue;
    end

    windowIds = resolve_window_ids_local(bladeWindows, bestByBlade, bladeId, P.step07v);
    fprintf('B%d windows to visualize: %s\n', bladeId, mat2str(windowIds));

    for windowId = windowIds
        WR = get_blade_window_result_local(bladeWindows, windowId);
        if isempty(WR) || isempty(WR.Bundle) || isempty(WR.Result)
            warning('B%d W%02d has incomplete result content. Skipping.', bladeId, windowId);
            continue;
        end

        if P.step07v.showTimeDomain
            plot_time_domain_comparison_local(WR, bladeId, figDir, P);
        end
        if P.step07v.showQueryDomain
            plot_query_domain_alignment_local(WR, bladeId, figDir, P);
        end
        if P.step07v.showStackedPass
            plot_stacked_pass_comparison_local(WR, bladeId, figDir, P);
        end
        if P.step07v.showGapComparison
            plot_direct_vs_gap_comparison_local(WR, bladeId, figDir, P);
        end
    end
end

fprintf('Figures: %s\n', figDir);

function P = apply_step07v_local_options_local(P, S07V)
P.step07v = struct();
P.step07v.blades = S07V.blades(:).';
P.step07v.windowIds = S07V.windowIds(:).';
P.step07v.useBestWindowWhenEmpty = logical(S07V.useBestWindowWhenEmpty);
P.step07v.passIds = S07V.passIds(:).';
P.step07v.maxPassCountWhenAuto = max(1, floor(S07V.maxPassCountWhenAuto));
P.step07v.showTimeDomain = logical(S07V.showTimeDomain);
P.step07v.showQueryDomain = logical(S07V.showQueryDomain);
P.step07v.showStackedPass = logical(S07V.showStackedPass);
P.step07v.showGapComparison = logical(S07V.showGapComparison);
P.step07v.showResidual = logical(S07V.showResidual);

P.sensors.analysis = S07V.analysisSensors(:).';
P.region.startTimeSec = S07V.startTimeSec;
P.view.saveFigures = logical(S07V.saveFigures);

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
suffix = output_suffix_local(S07V.outputLabel);
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
end

function windowIds = resolve_window_ids_local(bladeWindows, bestByBlade, bladeId, cfg)
if ~isempty(cfg.windowIds)
    windowIds = cfg.windowIds(:).';
    return;
end

if cfg.useBestWindowWhenEmpty
    idx = find([bestByBlade.blade_id] == bladeId, 1, 'first');
    if ~isempty(idx) && isfield(bestByBlade(idx), 'window_id') && isfinite(bestByBlade(idx).window_id)
        windowIds = bestByBlade(idx).window_id;
        return;
    end
end

windowIds = unique([bladeWindows.window_id], 'stable');
end

function WR = get_blade_window_result_local(bladeWindows, windowId)
WR = [];
idx = find([bladeWindows.window_id] == windowId, 1, 'first');
if ~isempty(idx)
    WR = bladeWindows(idx);
end
end

function plot_time_domain_comparison_local(WR, bladeId, figDir, P)
bundle = WR.Bundle;
R = WR.Result;
sensorIds = bundle.sensor_ids(:).';

fig = figure('Name', sprintf('Step06V time-domain B%d W%02d', bladeId, WR.window_id), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 17, 4.2 * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rows = find(bundle.sensor_index == is);
    [tMs, order] = sort(1e3 * (bundle.T(rows) - min(bundle.T(rows), [], 'omitnan')));
    rows = rows(order);
    nexttile;
    plot(tMs, bundle.V(rows), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.95, ...
        'DisplayName', 'high-speed observed'); hold on;
    plot(tMs, R.V_pred(rows), 'r-', 'LineWidth', 1.25, ...
        'DisplayName', 'identified prediction');
    ylabel(sprintf('CH%d (V)', sid), 'Interpreter', 'tex');
    title(sprintf('B%d W%02d CH%d, chosen=%s, RMSE=%.5f V, points=%d', ...
        bladeId, WR.window_id, sid, WR.ChosenPass, R.weighted_voltage_rmse, numel(rows)), ...
        'FontWeight', 'normal');
    if is == numel(sensorIds)
        xlabel('Selected-point time (ms)', 'Interpreter', 'tex');
    end
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();
end

sgtitle(sprintf('Time-domain comparison: EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm', ...
    R.EO_id, R.fn_id, R.A_id, R.dx_c_id), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, sprintf( ...
        'Step06V_TimeDomain_B%d_W%02d_%s.png', bladeId, WR.window_id, sensor_tag_local(sensorIds))), ...
        'Resolution', 300);
end
end

function plot_query_domain_alignment_local(WR, bladeId, figDir, P)
bundle = WR.Bundle;
R = WR.Result;
sensorIds = bundle.sensor_ids(:).';

fig = figure('Name', sprintf('Step06V query-domain B%d W%02d', bladeId, WR.window_id), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 17, 4.5 * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rows = find(bundle.sensor_index == is);
    xLo = bundle.x_domain_by_sensor(is, 1);
    xHi = bundle.x_domain_by_sensor(is, 2);
    xGrid = linspace(xLo, xHi, 1000).';
    vTpl = bundle.interp_v{is}(xGrid);
    eta = 0;
    if numel(R.sensor_eta_id) >= is
        eta = R.sensor_eta_id(is);
    end
    xQuery = bundle.X(rows) - R.dx_c_id - eta - R.u_est(rows);
    [xQuery, order] = sort(xQuery);
    rows = rows(order);

    nexttile;
    plot(xGrid, vTpl, 'k-', 'LineWidth', 1.4, 'DisplayName', 'low-speed template'); hold on;
    scatter(xQuery, bundle.V(rows), 10, [0.2 0.45 0.85], 'filled', ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55, ...
        'DisplayName', 'observed mapped to query x');
    plot(xQuery, R.V_pred(rows), 'r-', 'LineWidth', 1.1, ...
        'DisplayName', 'prediction on query x');
    xline(xLo, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
    xline(xHi, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
    ylabel(sprintf('CH%d (V)', sid), 'Interpreter', 'tex');
    title(sprintf('B%d W%02d CH%d, x_{query} range = [%.3f, %.3f] mm', ...
        bladeId, WR.window_id, sid, min(xQuery, [], 'omitnan'), max(xQuery, [], 'omitnan')), ...
        'FontWeight', 'normal');
    if is == numel(sensorIds)
        xlabel('Query coordinate x (mm)', 'Interpreter', 'tex');
    end
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();
end

sgtitle(sprintf('Query-domain template alignment: EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm', ...
    R.EO_id, R.fn_id, R.A_id, R.dx_c_id), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, sprintf( ...
        'Step06V_QueryDomain_B%d_W%02d_%s.png', bladeId, WR.window_id, sensor_tag_local(sensorIds))), ...
        'Resolution', 300);
end
end

function plot_stacked_pass_comparison_local(WR, bladeId, figDir, P)
bundle = WR.Bundle;
R = WR.Result;
sensorIds = bundle.sensor_ids(:).';
passIds = resolve_pass_ids_local(bundle, P.step07v);
if isempty(passIds)
    warning('B%d W%02d has no valid pass segmentation.', bladeId, WR.window_id);
    return;
end

    % One pass per figure makes cross-sensor phase/numbering easier to read.
for ip = 1:numel(passIds)
    passId = passIds(ip);
    [segmentsBySensor, validSensors] = collect_pass_rows_local(bundle, passId);
    if nnz(validSensors) < 1
        continue;
    end

    fig = figure('Name', sprintf('Step06V stacked-pass B%d W%02d P%02d', bladeId, WR.window_id, passId), ...
        'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 17, 10], ...
        'NumberTitle', 'off');
    ax = axes(fig); hold(ax, 'on');

    offsetStep = estimate_offset_step_local(bundle, R, segmentsBySensor);
    yTick = nan(1, numel(sensorIds));
    yLabel = cell(1, numel(sensorIds));
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        rows = segmentsBySensor{is};
        if isempty(rows)
            continue;
        end
        [tMs, order] = sort(1e3 * (bundle.T(rows) - min(bundle.T(rows), [], 'omitnan')));
        rows = rows(order);
        yOffset = (numel(sensorIds) - is) * offsetStep;
        yTick(is) = yOffset;
        yLabel{is} = sprintf('CH%d', sid);

        if is == 1
            rawHandleVisibility = 'on';
            predHandleVisibility = 'on';
        else
            rawHandleVisibility = 'off';
            predHandleVisibility = 'off';
        end
        plot(ax, tMs, bundle.V(rows) + yOffset, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, ...
            'DisplayName', ternary_label_local(is == 1, 'observed', ''), ...
            'HandleVisibility', rawHandleVisibility);
        plot(ax, tMs, R.V_pred(rows) + yOffset, 'r-', 'LineWidth', 1.2, ...
            'DisplayName', ternary_label_local(is == 1, 'predicted', ''), ...
            'HandleVisibility', predHandleVisibility);
        text(ax, tMs(1), yOffset + 0.04 * offsetStep, sprintf('CH%d', sid), ...
            'FontName', 'Times New Roman', 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end

    xlabel(ax, 'Pass-local time (ms)', 'Interpreter', 'tex');
    ylabel(ax, 'Stacked voltage', 'Interpreter', 'tex');
    title(ax, sprintf(['B%d W%02d pass %d, chosen=%s, EO=%d, f=%.3f Hz, ' ...
        'A=%.4f mm, dx_c=%.4f mm'], ...
        bladeId, WR.window_id, passId, WR.ChosenPass, R.EO_id, R.fn_id, R.A_id, R.dx_c_id), ...
        'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'Box', 'off');
    keepTick = find(isfinite(yTick));
    [yTickSorted, order] = sort(yTick(keepTick), 'ascend');
    yLabelSorted = yLabel(keepTick);
    yLabelSorted = yLabelSorted(order);
    set(ax, 'YTick', yTickSorted, 'YTickLabel', yLabelSorted);
    style_axes_local();

    if P.view.saveFigures
        exportgraphics(fig, fullfile(figDir, sprintf( ...
            'Step06V_StackedPass_B%d_W%02d_P%02d_%s.png', ...
            bladeId, WR.window_id, passId, sensor_tag_local(sensorIds))), ...
            'Resolution', 300);
    end
end
end

function passIds = resolve_pass_ids_local(bundle, cfg)
segmentsBySensor = cell(numel(bundle.sensor_ids), 1);
minSegmentCount = inf;
for is = 1:numel(bundle.sensor_ids)
    rows = find(bundle.sensor_index == is);
    segments = split_sensor_points_into_passes_local(bundle.T(rows));
    segmentsBySensor{is} = segments; %#ok<NASGU>
    minSegmentCount = min(minSegmentCount, numel(segments));
end

if ~isfinite(minSegmentCount) || minSegmentCount < 1
    passIds = [];
    return;
end

if ~isempty(cfg.passIds)
    passIds = cfg.passIds(cfg.passIds >= 1 & cfg.passIds <= minSegmentCount);
    return;
end

centerPass = max(1, ceil(minSegmentCount / 2));
halfSpan = floor((cfg.maxPassCountWhenAuto - 1) / 2);
firstPass = max(1, centerPass - halfSpan);
lastPass = min(minSegmentCount, firstPass + cfg.maxPassCountWhenAuto - 1);
firstPass = max(1, lastPass - cfg.maxPassCountWhenAuto + 1);
passIds = firstPass:lastPass;
end

function [segmentsBySensor, validSensors] = collect_pass_rows_local(bundle, passId)
segmentsBySensor = cell(numel(bundle.sensor_ids), 1);
validSensors = false(numel(bundle.sensor_ids), 1);
for is = 1:numel(bundle.sensor_ids)
    rows = find(bundle.sensor_index == is);
    segments = split_sensor_points_into_passes_local(bundle.T(rows));
    if numel(segments) >= passId
        segmentsBySensor{is} = rows(segments{passId});
        validSensors(is) = true;
    else
        segmentsBySensor{is} = [];
    end
end
end

function offsetStep = estimate_offset_step_local(bundle, R, segmentsBySensor)
spanAll = [];
for is = 1:numel(segmentsBySensor)
    rows = segmentsBySensor{is};
    if isempty(rows)
        continue;
    end
    spanAll(end + 1, 1) = range(bundle.V(rows)); %#ok<AGROW>
    spanAll(end + 1, 1) = range(R.V_pred(rows)); %#ok<AGROW>
end
offsetStep = 1.3 * max(spanAll, [], 'omitnan');
if ~isfinite(offsetStep) || offsetStep <= 0
    offsetStep = 1;
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
    sortedIdx = (edges(i) + 1):edges(i + 1);
    segments{i} = order(sortedIdx);
end
end

function plot_direct_vs_gap_comparison_local(WR, bladeId, figDir, P)
if ~isfield(WR, 'GapLibraryCompare') || isempty(WR.GapLibraryCompare) || ...
        ~isstruct(WR.GapLibraryCompare) || ...
        ~isfield(WR.GapLibraryCompare, 'gap_best_fit') || isempty(WR.GapLibraryCompare.gap_best_fit) || ...
        ~isfield(WR.GapLibraryCompare, 'gap_bundle') || isempty(WR.GapLibraryCompare.gap_bundle)
    warning('B%d W%02d has no gap-library comparison result to visualize.', bladeId, WR.window_id);
    return;
end

bundleDir = WR.Bundle;
bundleGap = WR.GapLibraryCompare.gap_bundle;
Rdir = WR.Result;
Rgap = WR.GapLibraryCompare.gap_best_fit;
sensorIds = bundleDir.sensor_ids(:).';

rowsPerSensor = 2 + double(P.step07v.showResidual);
fig = figure('Name', sprintf('Step06V direct-vs-gap B%d W%02d', bladeId, WR.window_id), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1, 1, 17, 3.8 * rowsPerSensor * numel(sensorIds)], ...
    'NumberTitle', 'off');
tiledlayout(fig, rowsPerSensor * numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rowsDir = find(bundleDir.sensor_index == is);
    rowsGap = find(bundleGap.sensor_index == is);
    if isempty(rowsDir) || isempty(rowsGap)
        continue;
    end

    [tDirMs, orderDir] = sort(1e3 * (bundleDir.T(rowsDir) - min(bundleDir.T(rowsDir), [], 'omitnan')));
    rowsDir = rowsDir(orderDir);
    [tGapMs, orderGap] = sort(1e3 * (bundleGap.T(rowsGap) - min(bundleGap.T(rowsGap), [], 'omitnan')));
    rowsGap = rowsGap(orderGap);

    nexttile;
    plot(tDirMs, bundleDir.V(rowsDir), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 0.9, ...
        'DisplayName', 'observed'); hold on;
    plot(tDirMs, Rdir.V_pred(rowsDir), 'r-', 'LineWidth', 1.15, ...
        'DisplayName', 'direct');
    plot(tGapMs, Rgap.Vpred(rowsGap) / 1000, '-', 'Color', [0.10 0.45 0.85], 'LineWidth', 1.10, ...
        'DisplayName', sprintf('gap-%s', char(WR.GapLibraryCompare.gap_best_model_name)));
    ylabel(sprintf('CH%d (V)', sid), 'Interpreter', 'tex');
    title(sprintf('B%d W%02d CH%d time-domain overlay', bladeId, WR.window_id, sid), ...
        'FontWeight', 'normal');
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();

    nexttile;
    xLo = bundleDir.x_domain_by_sensor(is, 1);
    xHi = bundleDir.x_domain_by_sensor(is, 2);
    xGrid = linspace(xLo, xHi, 1000).';
    vTpl = bundleDir.interp_v{is}(xGrid);
    etaDir = 0;
    if numel(Rdir.sensor_eta_id) >= is
        etaDir = Rdir.sensor_eta_id(is);
    end
    xQueryDir = bundleDir.X(rowsDir) - Rdir.dx_c_id - etaDir - Rdir.u_est(rowsDir);

    etaGap = 0;
    if numel(Rgap.sensorEtaMm) >= is
        etaGap = Rgap.sensorEtaMm(is);
    end
    xQueryGap = bundleGap.X(rowsGap) - Rgap.dxMm - etaGap - Rgap.uMm(rowsGap);

    plot(xGrid, vTpl, 'k-', 'LineWidth', 1.25, 'DisplayName', 'template'); hold on;
    scatter(xQueryDir, bundleDir.V(rowsDir), 10, [0.85 0.25 0.20], 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35, 'DisplayName', 'direct mapped');
    scatter(xQueryGap, bundleGap.V(rowsGap) / 1000, 10, [0.10 0.45 0.85], 'filled', ...
        'MarkerFaceAlpha', 0.30, 'MarkerEdgeAlpha', 0.30, 'DisplayName', 'gap mapped');
    xline(xLo, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    xline(xHi, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    ylabel(sprintf('CH%d (V)', sid), 'Interpreter', 'tex');
    title(sprintf('CH%d query-domain alignment', sid), 'FontWeight', 'normal');
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();

    if P.step07v.showResidual
        nexttile;
        resDir = bundleDir.V(rowsDir) - Rdir.V_pred(rowsDir);
        resGap = bundleGap.V(rowsGap) / 1000 - Rgap.Vpred(rowsGap) / 1000;
        plot(tDirMs, resDir, '-', 'Color', [0.85 0.25 0.20], 'LineWidth', 1.0, ...
            'DisplayName', sprintf('direct %.4f V', Rdir.weighted_voltage_rmse)); hold on;
        plot(tGapMs, resGap, '-', 'Color', [0.10 0.45 0.85], 'LineWidth', 1.0, ...
            'DisplayName', sprintf('gap %.4f V', Rgap.weightedRmseMv / 1000));
        yline(0, 'k:', 'HandleVisibility', 'off');
        ylabel(sprintf('CH%d res. (V)', sid), 'Interpreter', 'tex');
        if is == numel(sensorIds)
            xlabel('Selected-point time (ms)', 'Interpreter', 'tex');
        end
        if is == 1
            legend('Location', 'best', 'Box', 'off');
        end
        style_axes_local();
    end
end

sgtitle(sprintf(['Direct vs gap waveform comparison: B%d W%02d, direct EO=%d, gap model=%s, gap EO=%d, ' ...
    'RMSE %.5f V vs %.5f V'], ...
    bladeId, WR.window_id, Rdir.EO_id, char(WR.GapLibraryCompare.gap_best_model_name), ...
    Rgap.EO, Rdir.weighted_voltage_rmse, Rgap.weightedRmseMv / 1000), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, sprintf( ...
        'Step06V_DirectVsGap_B%d_W%02d_%s.png', ...
        bladeId, WR.window_id, sensor_tag_local(sensorIds))), 'Resolution', 300);
end
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
grid(gca, 'off');
end

function require_file_local(filePath, label)
if exist(filePath, 'file') ~= 2
    error('Missing %s:\n  %s', label, filePath);
end
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

function out = ternary_label_local(tf, valueIfTrue, valueIfFalse)
if tf
    out = valueIfTrue;
else
    out = valueIfFalse;
end
end

