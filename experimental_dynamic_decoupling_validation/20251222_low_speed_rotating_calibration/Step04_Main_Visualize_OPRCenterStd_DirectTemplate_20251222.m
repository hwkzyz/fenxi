%% Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222
% Visualize and summarize the latest direct-template identification result.
% Tune parameters here, then run this file directly.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
C = CaseConfig();

%% Parameters to tune
P.name.resultSuffix = 'Main_DirectTemplate_OPRCenterStd_Dx035_WithEta020_Reg002_PrevWinPhaseSafe';
P.view.reconstructionWindowID = [];  % [] uses Result.BestWindow.window_id.
P.view.reconstructionPassID = [];    % [] uses the middle pass in the selected window.
P.view.maxDetailPasses = 3;          % number of single-pass columns to show.
P.view.saveFigures = false;
P.view.saveSummary = true;

fprintf('\n=== Step04: summary and trend figure ===\n');
step04_visualize_direct_template_embedded(routeDir, P, C);

function step04_visualize_direct_template_embedded(routeDir, P, C)
resultDir = fullfile(routeDir, 'output', 'identification');
figureDir = fullfile(routeDir, 'output', 'figures', 'main_direct_template_oprcenterstd');
if exist(figureDir, 'dir') ~= 7; mkdir(figureDir); end

resultFile = resolve_step03_result_file_local(resultDir, C, P.name.resultSuffix);

loaded = load(resultFile, 'Result');
Result = loaded.Result;
T = Result.Trend;

fig = figure('Name', sprintf('%s direct OPRCenterStd identification trend', C.caseTag), 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(T.window_id, T.EO_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5); hold on;
yline(mode(T.EO_id), ':', 'LineWidth', 1.0);
ylabel('EO');
title('Engine order');
style_axes_for_step04_local();

nexttile;
plot(T.window_id, T.fn_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('Frequency (Hz)');
title('Identified frequency');
style_axes_for_step04_local();

nexttile;
yyaxis left;
plot(T.window_id, T.weighted_voltage_rmse, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('RMSE (V)');
yyaxis right;
plot(T.window_id, T.A_id, 's-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('A (mm)');
xlabel('Window');
title('Fit quality and amplitude');
style_axes_for_step04_local();

pngFile = fullfile(figureDir, sprintf('Step04_DirectTemplate_OPRCenterStd_Trend_%s_%s.png', C.dataset, C.caseTag));
pdfFile = fullfile(figureDir, sprintf('Step04_DirectTemplate_OPRCenterStd_Trend_%s_%s.pdf', C.dataset, C.caseTag));
if P.view.saveFigures
    exportgraphics(fig, pngFile, 'Resolution', 300);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
end

summary = table();
summary.ResultFile = string(resultFile);
summary.DominantEO = mode(T.EO_id);
summary.EOConsistency = mean(T.EO_id == summary.DominantEO, 'omitnan');
summary.MeanFrequencyHz = mean(T.fn_id, 'omitnan');
summary.MedianFrequencyHz = median(T.fn_id, 'omitnan');
summary.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
summary.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
summary.MeanAmplitudeMM = mean(T.A_id, 'omitnan');
summaryCsv = fullfile(resultDir, sprintf('Step04_DirectTemplate_OPRCenterStd_Summary_%s_%s.csv', C.dataset, C.caseTag));
if P.view.saveSummary
    writetable(summary, summaryCsv);
end

plot_reconstruction_compare_local(Result, P);
plot_single_pass_reconstruction_detail_local(Result, P);

fprintf('Loaded result: %s\n', resultFile);
if P.view.saveFigures; fprintf('Saved figure: %s\n', pngFile); end
if P.view.saveSummary; fprintf('Saved summary: %s\n', summaryCsv); end
disp(summary);
end

function plot_reconstruction_compare_local(Result, P)
if isempty(Result.WindowResult)
    warning('Result.WindowResult is empty; reconstruction comparison skipped.');
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
fig = figure('Name', sprintf('20251222 reconstructed waveform comparison W%02d', windowID), ...
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
    title(sprintf('Window %02d, EO=%d, f=%.3f Hz, A=%.4f mm, RMSE=%.5f V', ...
        windowID, R.EO_id, R.fn_id, R.A_id, R.weighted_voltage_rmse), 'FontWeight', 'normal');
    if is == numel(bundle.sensor_ids)
        xlabel('Time in selected points (ms)');
    end
    if is == 1
        legend({'observed', 'reconstructed'}, 'Location', 'best', 'Box', 'off');
    end
    style_axes_for_step04_local();
end
end

function plot_single_pass_reconstruction_detail_local(Result, P)
if isempty(Result.WindowResult)
    warning('Result.WindowResult is empty; single-pass reconstruction skipped.');
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

fig = figure('Name', sprintf('20251222 single-pass reconstructed waveform detail W%02d', windowID), ...
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
        style_axes_for_step04_local();
    end
end
sgtitle(sprintf('Window %02d, laps %s, EO=%d, f=%.3f Hz, A=%.4f mm, RMSE=%.5f V', ...
    windowID, mat2str(WR.lap_range), R.EO_id, R.fn_id, R.A_id, R.weighted_voltage_rmse), ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'normal');
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

function style_axes_for_step04_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function resultFile = resolve_step03_result_file_local(resultDir, C, resultSuffix)
candidates = {
    fullfile(resultDir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s_%s_%s.mat', C.caseTag, resultSuffix, C.dataset))
    fullfile(resultDir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_%s.mat', C.bladeId, C.sensorTag, resultSuffix, C.dataset))
    };
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        resultFile = candidates{i};
        return;
    end
end
error('Direct OPRCenterStd Step03 result not found. Tried:\n  %s\n  %s', candidates{1}, candidates{2});
end
