%% Step06_VisualizeIdentificationDiagnostics_20241106
% Visualize Step06 direct-template identification results without rerunning
% identification. The figures help diagnose unstable sliding windows.

clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));

resultFile = '';          % leave empty to use the newest S57 result file
windowIds = 1:6;           % e.g. 1:6; empty = all windows in the result
bladeIds = [];            % e.g. 4; empty = all blades in the result
showFigures = true;
saveFigures = true;

sensorTag = 'S57';
resultDir = fullfile(scriptDir, 'output', 'new_flow', '06_identification', sensorTag);
if isempty(resultFile)
    preferredFile = fullfile(resultDir, 'IdentificationResult_T075s_20241106_B4_only.mat');
    if exist(preferredFile, 'file') == 2
        resultFile = preferredFile;
    else
        files = dir(fullfile(resultDir, 'IdentificationResult_*_20241106_B4_only.mat'));
        if isempty(files)
            files = dir(fullfile(resultDir, 'IdentificationResult_*_20241106*.mat'));
        end
        if isempty(files)
            error('No Step06 IdentificationResult file found in %s.', resultDir);
        end
        [~, idx] = max([files.datenum]);
        resultFile = fullfile(files(idx).folder, files(idx).name);
    end
end

if ~exist(resultFile, 'file')
    if isempty(resultFile)
        files = dir(fullfile(resultDir, 'IdentificationResult_*_20241106*.mat'));
        if isempty(files)
            error('No Step06 IdentificationResult file found in %s.', resultDir);
        end
        [~, idx] = max([files.datenum]);
        resultFile = fullfile(files(idx).folder, files(idx).name);
    else
        error('Step06 result file does not exist: %s', resultFile);
    end
end

S = load(resultFile, 'IdentificationResult');
IdentificationResult = S.IdentificationResult;

figRoot = fullfile(scriptDir, 'output', 'new_flow', '07_step06_visualization', sensorTag);
if saveFigures && exist(figRoot, 'dir') ~= 7
    mkdir(figRoot);
end

fprintf('Loaded Step06 result: %s\n', resultFile);
plot_step06_trend_dashboard_local(IdentificationResult, figRoot, showFigures, saveFigures);
plot_step06_window_diagnostics_local(IdentificationResult, windowIds, bladeIds, figRoot, showFigures, saveFigures);


function plot_step06_trend_dashboard_local(R, figRoot, showFigures, saveFigures)
T = R.Trend;
fig = figure('Name', 'Step06 trend dashboard', 'Color', 'w', ...
    'Position', [80, 60, 1450, 900], 'Visible', visible_state_local(showFigures));
tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_trend_metric_local(T, 'EO', 'EO', 'Selected EO');
plot_trend_metric_local(T, 'FrequencyHz', 'Frequency (Hz)', 'Identified frequency');
plot_trend_metric_local(T, 'AmplitudeMM', 'A (mm)', 'Amplitude');
plot_trend_metric_local(T, 'DxCMM', 'dx_c (mm)', 'Center offset');
plot_trend_metric_local(T, 'WeightedVoltageRMSE', 'RMSE (V)', 'Weighted voltage RMSE');
plot_trend_metric_local(T, 'SensorEtaMaxAbsMM', 'max |eta_s| (mm)', 'Sensor eta usage');
plot_trend_metric_local(T, 'PointCount', 'points', 'Selected point count');
plot_trend_metric_local(T, 'ClampFraction', 'fraction', 'Query clamp fraction');

if saveFigures
    exportgraphics(fig, fullfile(figRoot, 'Step06_TrendDashboard.png'), 'Resolution', 250);
end
end


function plot_trend_metric_local(T, fieldName, yLabelText, titleText)
nexttile; hold on; grid on; box on;
bladeIds = unique(T.BladeID).';
for bladeId = bladeIds
    rows = T(T.BladeID == bladeId, :);
    y = rows.(fieldName);
    if iscell(y)
        continue;
    end
    plot(rows.WindowID, y, 'o-', 'LineWidth', 1.2, 'DisplayName', sprintf('B%d', bladeId));
end
xlabel('Window ID');
ylabel(yLabelText);
title(titleText);
legend('Location', 'best');
end


function plot_step06_window_diagnostics_local(R, windowIds, bladeIds, figRoot, showFigures, saveFigures)
windows = R.WindowResult(:);
keep = true(size(windows));
if ~isempty(windowIds)
    keep = keep & ismember([windows.window_id], windowIds);
end
if ~isempty(bladeIds)
    keep = keep & ismember([windows.blade_id], bladeIds);
end
windows = windows(keep);
if isempty(windows)
    warning('No Step06 windows match the requested filters.');
    return;
end

for i = 1:numel(windows)
    W = windows(i);
    plot_one_step06_window_local(W, figRoot, showFigures, saveFigures);
end
end


function plot_one_step06_window_local(W, figRoot, showFigures, saveFigures)
bundle = W.Bundle;
result = W.Result;
sensorIds = bundle.sensor_ids(:).';
nSensor = numel(sensorIds);

fig = figure('Name', sprintf('Step06 B%d W%02d diagnostics', W.blade_id, W.window_id), ...
    'Color', 'w', 'Position', [60, 40, 1550, 920], 'Visible', visible_state_local(showFigures));
tiledlayout(fig, max(2, nSensor + 1), 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:nSensor
    sid = sensorIds(is);
    mask = bundle.sensor_index == is;
    x = bundle.X(mask);
    t = bundle.T(mask);
    v = bundle.V(mask);
    vp = result.V_pred(mask);
    u = result.u_est(mask);
    w = bundle.W(mask);
    Tpl = template_from_window_sensor_local(bundle.WaveformSensor(is));
    eta = sensor_eta_at_local(result, is);
    dx = result.dx_c_id;

    nexttile; hold on; grid on; box on;
    plot(Tpl.x, Tpl.v, 'k-', 'LineWidth', 1.4, 'DisplayName', 'low-speed template');
    scatter(x, v, 8, w, 'filled', 'DisplayName', 'selected high-speed points');
    scatter(x - dx - eta - u, vp, 8, 'r', 'filled', 'DisplayName', 'predicted query');
    xlabel('x / query x (mm)');
    ylabel('V (V)');
    title(sprintf('CH%d waveform match', sid));
    legend('Location', 'best');
    colorbar;

    nexttile; hold on; grid on; box on;
    [tOrder, order] = sort(t);
    plot(tOrder, v(order), '.', 'Color', [0.2 0.45 0.8], 'DisplayName', 'measured');
    plot(tOrder, vp(order), '-', 'Color', [0.85 0.2 0.15], 'LineWidth', 1.1, 'DisplayName', 'predicted');
    xlabel('Time (s)');
    ylabel('V (V)');
    title(sprintf('CH%d time waveform', sid));
    legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    residual = v - vp;
    plot(tOrder, residual(order), '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0);
    yline(0, 'k:');
    xlabel('Time (s)');
    ylabel('V_{meas}-V_{pred} (V)');
    title(sprintf('CH%d residual, RMSE %.4g V', sid, sqrt(mean(residual.^2, 'omitnan'))));
end

nexttile([1 3]); hold on; grid on; box on;
[tAll, orderAll] = sort(bundle.T);
uAll = result.u_est(orderAll);
plot(tAll, uAll, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('u(t) (mm)');
title(sprintf(['B%d W%02d laps %s | EO=%d, f=%.3f Hz, A=%.4f mm, dx_c=%.4f mm, ' ...
    'weighted RMSE=%.5f V, chosen=%s'], ...
    W.blade_id, W.window_id, mat2str(W.lap_range), result.EO_id, result.fn_id, ...
    result.A_id, result.dx_c_id, result.weighted_voltage_rmse, string(W.ChosenPass)));

if saveFigures
    bladeDir = fullfile(figRoot, sprintf('B%d', W.blade_id));
    if exist(bladeDir, 'dir') ~= 7
        mkdir(bladeDir);
    end
    fileName = sprintf('Step06_B%d_W%02d_WaveformDiagnostics.png', W.blade_id, W.window_id);
    exportgraphics(fig, fullfile(bladeDir, fileName), 'Resolution', 250);
end
end


function Tpl = template_from_window_sensor_local(S)
Tpl = struct();
Tpl.x = S.template_x(:);
Tpl.v = S.template_v(:);
end


function eta = sensor_eta_at_local(result, is)
eta = 0;
if isfield(result, 'sensor_eta_id') && numel(result.sensor_eta_id) >= is && isfinite(result.sensor_eta_id(is))
    eta = result.sensor_eta_id(is);
end
end


function state = visible_state_local(tf)
if tf
    state = 'on';
else
    state = 'off';
end
end
