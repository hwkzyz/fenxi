%% Step05_Main_Visualize_OPRCenterStd_DirectTemplate_20241106
% Visualize and summarize the 20241106 direct-template identification result.
% Tune parameters here, then run this file directly.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));

%% Parameters to tune
P.result.bladeIds = 1:6;
P.result.sensorTag = 'S57';
P.result.analysisStartTimeSec = 75.0;
P.result.eoLock = 12;
P.result.runLabel = make_time_label_local(P.result.analysisStartTimeSec);
P.result.resultSuffix = make_result_suffix_local(P.result.eoLock, P.result.runLabel);
P.compare.freeSuffix = make_result_suffix_local([], P.result.runLabel);  % optional; skipped if the file is absent.

P.view.reconstructionBladeID = 6;
P.view.reconstructionWindowID = [];  % [] uses Result.BestWindow.window_id.
P.view.reconstructionPassID = [];    % [] uses the middle pass in the selected window.
P.view.maxDetailPasses = 3;
P.view.dxCLimitMM = 0.35;
P.view.saveFigures = true;
P.view.saveSummary = true;

fprintf('\n=== Step05: 20241106 direct-template visualization ===\n');
step05_visualize_direct_template_embedded(routeDir, P);

function step05_visualize_direct_template_embedded(routeDir, P)
resultDir = fullfile(routeDir, 'output', 'identification');
figureDir = fullfile(routeDir, 'output', 'figures', 'main_direct_template_oprcenterstd_20241106');
if exist(figureDir, 'dir') ~= 7; mkdir(figureDir); end

[Summary, Results] = load_all_blade_results_local(resultDir, P);
if isempty(Summary)
    error('No direct-template results found in %s.', resultDir);
end

trendFig = plot_all_blade_trends_local(Summary, Results, P);
trendPng = fullfile(figureDir, sprintf('Step05_AllBlade_DirectTemplate_%s_Trend_20241106.png', P.result.resultSuffix));
trendPdf = fullfile(figureDir, sprintf('Step05_AllBlade_DirectTemplate_%s_Trend_20241106.pdf', P.result.resultSuffix));
if P.view.saveFigures
    exportgraphics(trendFig, trendPng, 'Resolution', 300);
    exportgraphics(trendFig, trendPdf, 'ContentType', 'vector');
end

summaryCsv = fullfile(resultDir, sprintf('Step05_DirectTemplate_OPRCenterStd_Summary_%s_20241106.csv', P.result.resultSuffix));
if P.view.saveSummary
    writetable(Summary, summaryCsv);
end

plot_reconstruction_compare_local(Results, P, figureDir);
plot_selection_audit_local(Results, P, figureDir);
plot_single_pass_reconstruction_detail_local(Results, P, figureDir);

fprintf('Loaded %d blade result(s).\n', height(Summary));
if P.view.saveFigures; fprintf('Saved trend figure: %s\n', trendPng); end
if P.view.saveSummary; fprintf('Saved summary: %s\n', summaryCsv); end
disp(Summary);
end

function [Summary, Results] = load_all_blade_results_local(resultDir, P)
rows = repmat(empty_summary_row_local(), numel(P.result.bladeIds), 1);
Results = struct('blade_id', {}, 'Result', {});
rowCount = 0;
for i = 1:numel(P.result.bladeIds)
    bladeID = P.result.bladeIds(i);
    resultFile = fullfile(resultDir, sprintf( ...
        'Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_20241106.mat', ...
        bladeID, P.result.sensorTag, P.result.resultSuffix));
    if exist(resultFile, 'file') ~= 2
        warning('Result file not found for B%d: %s', bladeID, resultFile);
        continue;
    end

    loaded = load(resultFile, 'Result');
    Result = loaded.Result;
    rowCount = rowCount + 1;
    rows(rowCount) = summarize_result_local(bladeID, resultFile, Result, P);
    Results(rowCount).blade_id = bladeID;
    Results(rowCount).Result = Result;
end

if rowCount == 0
    Summary = table();
    Results = struct('blade_id', {}, 'Result', {});
else
    Summary = struct2table(rows(1:rowCount));
end
end

function row = summarize_result_local(bladeID, resultFile, Result, P)
T = Result.Trend;
row = empty_summary_row_local();
row.BladeID = bladeID;
row.ResultFile = string(resultFile);
row.DominantEO = mode(T.EO_id);
row.EOConsistency = mean(T.EO_id == row.DominantEO, 'omitnan');
row.MeanFrequencyHz = mean(T.fn_id, 'omitnan');
row.MedianFrequencyHz = median(T.fn_id, 'omitnan');
row.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
row.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
row.MeanAmplitudeMM = mean(T.A_id, 'omitnan');
row.MedianAmplitudeMM = median(T.A_id, 'omitnan');
row.MedianDxCMM = median(T.dx_c_id, 'omitnan');
row.MinDxCMM = min(T.dx_c_id, [], 'omitnan');
row.MaxDxCMM = max(T.dx_c_id, [], 'omitnan');
row.DxCAtLimitCount = nnz(abs(T.dx_c_id) >= (P.view.dxCLimitMM - 1e-3));
row.WindowCount = height(T);

if isfield(Result, 'DynamicMapFile') && exist(Result.DynamicMapFile, 'file') == 2
    loadedMap = load(Result.DynamicMapFile, 'DynamicMap');
    DynamicMap = loadedMap.DynamicMap;
    row.DynamicStartSec = DynamicMap.GlobalTimeWindow(1);
    row.DynamicEndSec = DynamicMap.GlobalTimeWindow(2);
    if isfield(DynamicMap, 'Window')
        rot = [DynamicMap.Window.rot_freq_mean_hz];
        row.RotFreqMedianHz = median(rot, 'omitnan');
        row.RpmMedian = 60 * row.RotFreqMedianHz;
    end
end

freeFile = fullfile(fileparts(resultFile), sprintf( ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_20241106.mat', ...
    bladeID, P.result.sensorTag, P.compare.freeSuffix));
if exist(freeFile, 'file') == 2
    F = load(freeFile, 'Result');
    TF = F.Result.Trend;
    row.FreeDominantEO = F.Result.ResonanceSummary.dominant_eo;
    row.FreeMedianRMSE = median(TF.weighted_voltage_rmse, 'omitnan');
    row.FreeMedianDxCMM = median(TF.dx_c_id, 'omitnan');
    row.FreeDxCAtLimitCount = nnz(abs(TF.dx_c_id) >= (P.view.dxCLimitMM - 1e-3));
    row.LockMinusFreeMedianRMSE = row.MedianRMSE - row.FreeMedianRMSE;
end
end

function row = empty_summary_row_local()
row = struct( ...
    'BladeID', NaN, ...
    'ResultFile', string(missing), ...
    'DominantEO', NaN, ...
    'EOConsistency', NaN, ...
    'MeanFrequencyHz', NaN, ...
    'MedianFrequencyHz', NaN, ...
    'MeanRMSE', NaN, ...
    'MedianRMSE', NaN, ...
    'MeanAmplitudeMM', NaN, ...
    'MedianAmplitudeMM', NaN, ...
    'MedianDxCMM', NaN, ...
    'MinDxCMM', NaN, ...
    'MaxDxCMM', NaN, ...
    'DxCAtLimitCount', NaN, ...
    'WindowCount', NaN, ...
    'DynamicStartSec', NaN, ...
    'DynamicEndSec', NaN, ...
    'RotFreqMedianHz', NaN, ...
    'RpmMedian', NaN, ...
    'FreeDominantEO', NaN, ...
    'FreeMedianRMSE', NaN, ...
    'FreeMedianDxCMM', NaN, ...
    'FreeDxCAtLimitCount', NaN, ...
    'LockMinusFreeMedianRMSE', NaN);
end

function fig = plot_all_blade_trends_local(Summary, Results, P)
fig = figure('Name', '20241106 direct OPRCenterStd all-blade trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 16]);
tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(Summary.BladeID, Summary.DominantEO);
ylabel('EO');
title('Dominant engine order');
style_axes_for_step05_local();

nexttile;
bar(Summary.BladeID, Summary.MedianFrequencyHz);
ylabel('Frequency (Hz)');
title('Median frequency');
style_axes_for_step05_local();

nexttile;
bar(Summary.BladeID, Summary.MedianAmplitudeMM);
ylabel('A (mm)');
title('Median amplitude');
style_axes_for_step05_local();

nexttile;
bar(Summary.BladeID, Summary.MedianDxCMM);
yline(P.view.dxCLimitMM, 'r:', 'LineWidth', 1.0);
yline(-P.view.dxCLimitMM, 'r:', 'LineWidth', 1.0);
ylabel('dx_c (mm)');
title('Median static shift');
style_axes_for_step05_local();

nexttile;
bar(Summary.BladeID, Summary.DxCAtLimitCount);
ylabel('Count');
title('dx_c boundary hits');
style_axes_for_step05_local();

nexttile;
bar(Summary.BladeID, Summary.MedianRMSE);
ylabel('RMSE (V)');
title('Median weighted RMSE');
style_axes_for_step05_local();

nexttile([1 2]);
hold on;
colors = lines(numel(Results));
for i = 1:numel(Results)
    T = Results(i).Result.Trend;
    plot(T.window_id, T.A_id, 'o-', 'LineWidth', 1.0, 'MarkerSize', 3.5, ...
        'Color', colors(i, :), 'DisplayName', sprintf('B%d', Results(i).blade_id));
end
xlabel('Window');
ylabel('A (mm)');
title('Amplitude trend by window');
legend('Location', 'bestoutside', 'Box', 'off');
style_axes_for_step05_local();
end

function plot_reconstruction_compare_local(Results, P, figureDir)
[Result, bladeID] = find_result_for_blade_local(Results, P.view.reconstructionBladeID);
if isempty(Result)
    warning('No result available for reconstruction blade B%d.', P.view.reconstructionBladeID);
    return;
end
if isempty(P.view.reconstructionWindowID)
    windowID = Result.BestWindow.window_id;
else
    windowID = P.view.reconstructionWindowID;
end
if windowID < 1 || windowID > numel(Result.WindowResult)
    warning('Requested reconstruction window %d is outside available range.', windowID);
    return;
end
WR = Result.WindowResult(windowID);
if ~isfield(WR, 'bundle') || ~isfield(WR, 'Result') || ~isfield(WR.Result, 'V_pred')
    warning('Window %d does not contain bundle/Result.V_pred; reconstruction comparison skipped.', windowID);
    return;
end

bundle = WR.bundle;
R = WR.Result;
fig = figure('Name', sprintf('20241106 reconstructed waveform comparison B%d W%02d', bladeID, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 22, 14]);
tiledlayout(fig, numel(bundle.sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    mask = bundle.sensor_index == is;
    t_ms = 1e3 * (bundle.T(mask) - min(bundle.T(mask)));
    [t_ms, order] = sort(t_ms);
    v_obs = bundle.V(mask);
    v_pred = R.V_pred(mask);
    nexttile;
    plot(t_ms, v_obs(order), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 1.0); hold on;
    plot(t_ms, v_pred(order), 'r-', 'LineWidth', 1.2);
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('B%d W%02d, EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm, RMSE=%.5f V', ...
        bladeID, windowID, R.EO_id, R.fn_id, R.A_id, R.dx_c_id, R.weighted_voltage_rmse), ...
        'FontWeight', 'normal');
    if is == numel(bundle.sensor_ids)
        xlabel('Time in selected points (ms)');
    end
    if is == 1
        legend({'observed', 'reconstructed'}, 'Location', 'best', 'Box', 'off');
    end
    style_axes_for_step05_local();
end

pngFile = fullfile(figureDir, sprintf('Step05_Reconstruction_B%d_W%02d_%s_20241106.png', ...
    bladeID, windowID, P.result.resultSuffix));
if P.view.saveFigures
    exportgraphics(fig, pngFile, 'Resolution', 300);
end
end

function plot_selection_audit_local(Results, P, figureDir)
[Result, bladeID] = find_result_for_blade_local(Results, P.view.reconstructionBladeID);
if isempty(Result)
    warning('No result available for selection audit blade B%d.', P.view.reconstructionBladeID);
    return;
end
if isempty(P.view.reconstructionWindowID)
    windowID = Result.BestWindow.window_id;
else
    windowID = P.view.reconstructionWindowID;
end
if windowID < 1 || windowID > numel(Result.WindowResult)
    warning('Requested selection-audit window %d is outside available range.', windowID);
    return;
end
WR = Result.WindowResult(windowID);
if ~isfield(WR, 'bundle') || isempty(WR.bundle)
    warning('Window %d does not contain bundle; selection audit skipped.', windowID);
    return;
end
if ~isfield(Result, 'DynamicMapFile') || exist(Result.DynamicMapFile, 'file') ~= 2
    warning('DynamicMap file is unavailable; selection audit skipped.');
    return;
end
if ~isfield(Result, 'TemplateFile') || exist(Result.TemplateFile, 'file') ~= 2
    warning('Template file is unavailable; selection audit skipped.');
    return;
end

loadedMap = load(Result.DynamicMapFile, 'DynamicMap');
loadedTemplate = load(Result.TemplateFile, 'Template');
DynamicMap = loadedMap.DynamicMap;
Template = loadedTemplate.Template;
if windowID > numel(DynamicMap.Window)
    warning('Window %d is outside DynamicMap.Window; selection audit skipped.', windowID);
    return;
end

Wmap = DynamicMap.Window(windowID);
bundle = WR.bundle;
R = WR.Result;
sensorIds = bundle.sensor_ids;
fig = figure('Name', sprintf('20241106 selected dynamic waveform points B%d W%02d', bladeID, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 22, 15]);
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rawSensor = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
    tplSensor = Template.Sensor([Template.Sensor.sensor_id] == sid);
    selectedMask = bundle.sensor_index == is;

    nexttile;
    plot(tplSensor.x_grid, tplSensor.v_grid, 'k-', 'LineWidth', 1.5, ...
        'DisplayName', 'low-speed template'); hold on;
    scatter(rawSensor.x_rel, rawSensor.V, 7, [0.78 0.78 0.78], 'filled', ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', 'window samples');
    scatter(bundle.X(selectedMask), bundle.V(selectedMask), 11, bundle.W(selectedMask), 'filled', ...
        'DisplayName', 'used in fit');
    xline(tplSensor.x_domain(1), ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xline(tplSensor.x_domain(2), ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('B%d W%02d CH%d: used %d/%d points, x_{used}=[%.3f, %.3f] mm', ...
        bladeID, windowID, sid, nnz(selectedMask), numel(rawSensor.V), ...
        min(bundle.X(selectedMask), [], 'omitnan'), max(bundle.X(selectedMask), [], 'omitnan')), ...
        'FontWeight', 'normal');
    if is == numel(sensorIds)
        xlabel('x relative to template center (mm)');
    end
    if is == 1
        legend('Location', 'best', 'Box', 'off');
    end
    cb = colorbar;
    cb.Label.String = 'fit weight';
    style_axes_for_step05_local();
end

sgtitle(sprintf('Selected points for direct-template fit: EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm, RMSE=%.5f V', ...
    R.EO_id, R.fn_id, R.A_id, R.dx_c_id, R.weighted_voltage_rmse), ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'normal');

pngFile = fullfile(figureDir, sprintf('Step05_SelectionAudit_B%d_W%02d_%s_20241106.png', ...
    bladeID, windowID, P.result.resultSuffix));
pdfFile = fullfile(figureDir, sprintf('Step05_SelectionAudit_B%d_W%02d_%s_20241106.pdf', ...
    bladeID, windowID, P.result.resultSuffix));
if P.view.saveFigures
    exportgraphics(fig, pngFile, 'Resolution', 300);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
end
end

function plot_single_pass_reconstruction_detail_local(Results, P, figureDir)
[Result, bladeID] = find_result_for_blade_local(Results, P.view.reconstructionBladeID);
if isempty(Result)
    return;
end
if isempty(P.view.reconstructionWindowID)
    windowID = Result.BestWindow.window_id;
else
    windowID = P.view.reconstructionWindowID;
end
WR = Result.WindowResult(windowID);
if ~isfield(WR, 'bundle') || ~isfield(WR, 'Result') || ~isfield(WR.Result, 'V_pred')
    warning('Window %d does not contain bundle/Result.V_pred; single-pass reconstruction skipped.', windowID);
    return;
end

bundle = WR.bundle;
R = WR.Result;
segments_by_sensor = cell(numel(bundle.sensor_ids), 1);
min_segments = inf;
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    segments = split_sensor_points_into_passes_local(bundle.T(mask));
    segments_by_sensor{is} = segments;
    min_segments = min(min_segments, numel(segments));
end
if ~isfinite(min_segments) || min_segments < 1
    warning('No single-pass segments found for window %d.', windowID);
    return;
end

if isempty(P.view.reconstructionPassID)
    centerPass = max(1, ceil(min_segments / 2));
else
    centerPass = min(max(1, P.view.reconstructionPassID), min_segments);
end
maxCols = max(1, min(P.view.maxDetailPasses, min_segments));
firstPass = max(1, centerPass - floor((maxCols - 1) / 2));
lastPass = min(min_segments, firstPass + maxCols - 1);
firstPass = max(1, lastPass - maxCols + 1);
passIDs = firstPass:lastPass;

fig = figure('Name', sprintf('20241106 single-pass reconstructed waveform detail B%d W%02d', bladeID, windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 25, 14]);
tiledlayout(fig, numel(bundle.sensor_ids), numel(passIDs), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    sensorMask = bundle.sensor_index == is;
    sensorRows = find(sensorMask);
    segments = segments_by_sensor{is};
    for ip = 1:numel(passIDs)
        passID = passIDs(ip);
        rows = sensorRows(segments{passID});
        [t_ms, order] = sort(1e3 * (bundle.T(rows) - min(bundle.T(rows))));
        rows = rows(order);
        nexttile;
        plot(t_ms, bundle.V(rows), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 1.05); hold on;
        plot(t_ms, R.V_pred(rows), 'r-', 'LineWidth', 1.25);
        ylabel(sprintf('CH%d V', sid));
        if is == 1
            title(sprintf('Pass %d', passID), 'FontWeight', 'normal');
        end
        if is == numel(bundle.sensor_ids)
            xlabel('Time in pass (ms)');
        end
        if is == 1 && ip == numel(passIDs)
            legend({'observed', 'reconstructed'}, 'Location', 'best', 'Box', 'off');
        end
        style_axes_for_step05_local();
    end
end
sgtitle(sprintf('B%d W%02d, laps %s, EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm, RMSE=%.5f V', ...
    bladeID, windowID, mat2str(WR.lap_range), R.EO_id, R.fn_id, R.A_id, R.dx_c_id, R.weighted_voltage_rmse), ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'normal');

pngFile = fullfile(figureDir, sprintf('Step05_SinglePass_B%d_W%02d_%s_20241106.png', ...
    bladeID, windowID, P.result.resultSuffix));
if P.view.saveFigures
    exportgraphics(fig, pngFile, 'Resolution', 300);
end
end

function [Result, bladeID] = find_result_for_blade_local(Results, requestedBladeID)
Result = [];
bladeID = requestedBladeID;
idx = find([Results.blade_id] == requestedBladeID, 1, 'first');
if isempty(idx)
    if isempty(Results)
        return;
    end
    idx = 1;
    bladeID = Results(idx).blade_id;
end
Result = Results(idx).Result;
end

function segments = split_sensor_points_into_passes_local(t)
t = t(:);
if isempty(t)
    segments = {};
    return;
end
[t_sorted, order] = sort(t);
dt = diff(t_sorted);
positive_dt = dt(dt > 0);
if isempty(positive_dt)
    segments = {order};
    return;
end
gap_threshold = max(5 * median(positive_dt, 'omitnan'), prctile(positive_dt, 95));
breaks = find(dt > gap_threshold);
edges = [0; breaks(:); numel(t_sorted)];
segments = cell(numel(edges) - 1, 1);
for i = 1:numel(segments)
    sorted_idx = (edges(i) + 1):edges(i + 1);
    segments{i} = order(sorted_idx);
end
end

function style_axes_for_step05_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function suffix = make_result_suffix_local(eoLock, runLabel)
if isempty(eoLock)
    suffix = sprintf('Direct_%s', runLabel);
else
    suffix = sprintf('DirectEO%d_%s', eoLock, runLabel);
end
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
