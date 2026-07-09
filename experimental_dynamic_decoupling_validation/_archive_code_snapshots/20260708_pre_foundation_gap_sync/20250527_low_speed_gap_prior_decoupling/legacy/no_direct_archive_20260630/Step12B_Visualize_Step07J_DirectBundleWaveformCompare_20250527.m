%% Visualize Step07J waveform comparison: direct bundle vs no-direct bundle
% Reads saved Step07J main gap_tilt results and compares the reconstructed
% waveform in the same window for each sensor.
%
% Optional environment variables:
%   STEP07J_BUNDLE_WAVE_WINDOW      window index, default direct result best window
%   STEP07J_BUNDLE_WAVE_SENSORS     sensor ids separated by spaces, e.g. "1 3 6"
%   STEP07J_BUNDLE_WAVE_MAX_POINTS  plotted high-speed points per sensor/case

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_bundle_waveform_compare');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

directFile = fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', ...
    C0.dataset, C0.caseTag));
nodirectFile = fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt_nodirect.mat', ...
    C0.dataset, C0.caseTag));
if exist(directFile, 'file') ~= 2
    error('Missing direct-assisted Step07J result: %s', directFile);
end
if exist(nodirectFile, 'file') ~= 2
    error('Missing no-direct Step07J result. Run the no-direct ablation first: %s', nodirectFile);
end

Sd = load(directFile, 'Result');
Sn = load(nodirectFile, 'Result');
Direct = Sd.Result;
NoDirect = Sn.Result;

windowId = parse_positive_integer_env_local('STEP07J_BUNDLE_WAVE_WINDOW', 0);
if windowId <= 0
    windowId = Direct.BestWindowIndex;
end
if windowId < 1 || windowId > numel(Direct.WindowResult) || windowId > numel(NoDirect.WindowResult)
    error('Window %d is outside the available Step07J result range.', windowId);
end

sensorOverride = parse_integer_list_env_local('STEP07J_BUNDLE_WAVE_SENSORS');
maxPointsPerCase = parse_positive_integer_env_local('STEP07J_BUNDLE_WAVE_MAX_POINTS', 900);

wrD = Direct.WindowResult(windowId);
wrN = NoDirect.WindowResult(windowId);
fitD = wrD.modelFits.gap_tilt;
fitN = wrN.modelFits.gap_tilt;
bundleD = wrD.bundle;
bundleN = wrN.bundle;
if isempty(sensorOverride)
    sensorIds = intersect(bundleD.sensorIds(:).', bundleN.sensorIds(:).', 'stable');
else
    sensorIds = sensorOverride(:).';
end

style = paper_style_local();
fig = figure('Name', sprintf('Step07J direct bundle waveform compare W%d', windowId), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18.5, max(7.5, 3.65 * numel(sensorIds))]);
tiledlayout(numel(sensorIds), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

curveRows = {};
caseInfo = struct( ...
    'label', {'Direct bundle', 'No direct bundle'}, ...
    'shortLabel', {'direct_bundle', 'no_direct_bundle'}, ...
    'bundle', {bundleD, bundleN}, ...
    'fit', {fitD, fitN}, ...
    'color', {style.blue, style.orange});

fprintf('\n=== Step07J direct-bundle waveform comparison ===\n');
fprintf('Window: %d, laps direct %s, no-direct %s\n', ...
    windowId, mat2str(wrD.lapRange), mat2str(wrN.lapRange));
fprintf('Direct: EO%d, %.3f Hz, A %.4f mm, RMSE %.2f mV, points %d\n', ...
    fitD.EO, fitD.freqHz, fitD.amplitudeMm, fitD.weightedRmseMv, bundleD.pointCount);
fprintf('No-direct: EO%d, %.3f Hz, A %.4f mm, RMSE %.2f mV, points %d\n', ...
    fitN.EO, fitN.freqHz, fitN.amplitudeMm, fitN.weightedRmseMv, bundleN.pointCount);

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    yLimits = [inf, -inf];
    panelData = cell(1, numel(caseInfo));
    for ic = 1:numel(caseInfo)
        [panelData{ic}, yLimits] = prepare_case_sensor_waveform_local( ...
            caseInfo(ic).bundle, caseInfo(ic).fit, sid, maxPointsPerCase, yLimits);
    end

    for ic = 1:numel(caseInfo)
        nexttile; hold on; grid on; box on;
        P = panelData{ic};
        plot(P.xObserved, P.vObserved, '.', 'Color', style.gray, ...
            'MarkerSize', 3.3, 'DisplayName', sprintf('measured (%d pts)', P.rawPointCount));
        plot(P.xPred, P.vPred, '-', 'Color', caseInfo(ic).color, ...
            'LineWidth', 1.35, 'DisplayName', sprintf('gap\\_tilt prediction, EO%d', caseInfo(ic).fit.EO));
        xlabel('Equivalent x after vibration compensation (mm)', 'Interpreter', 'none');
        ylabel('Voltage (mV)');
        title(sprintf('CH%d | %s | RMSE %.2f mV', ...
            sid, caseInfo(ic).label, caseInfo(ic).fit.weightedRmseMv), ...
            'Interpreter', 'none');
        ylim(yLimits + [-1, 1] * 0.05 * diff(yLimits));
        if is == 1
            legend('Location', 'best', 'Box', 'off', 'Interpreter', 'none');
        end
        format_axes_local(gca, style);

        curveRows = append_curve_local(curveRows, windowId, sid, ...
            caseInfo(ic).shortLabel, "observed", P.xObserved, P.vObserved);
        curveRows = append_curve_local(curveRows, windowId, sid, ...
            caseInfo(ic).shortLabel, "prediction", P.xPred, P.vPred);
    end
end

sgtitle(sprintf(['Step07J gap\\_tilt waveform comparison, %s %s W%02d | ' ...
    'direct %.2f mV / no-direct %.2f mV'], ...
    C0.dataset, C0.caseTag, windowId, fitD.weightedRmseMv, fitN.weightedRmseMv), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

figBase = sprintf('Step07J_BundleWaveformCompare_%s_%s_W%02d', ...
    C0.dataset, C0.caseTag, windowId);
pngFile = fullfile(figDir, [figBase, '.png']);
pdfFile = fullfile(figDir, [figBase, '.pdf']);
svgFile = fullfile(figDir, [figBase, '.svg']);
csvFile = fullfile(figDir, [figBase, '_Curves.csv']);
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
print(fig, svgFile, '-dsvg');
if ~isempty(curveRows)
    CurveTable = struct2table(vertcat(curveRows{:}));
    writetable(CurveTable, csvFile);
end

fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', pngFile, pdfFile, svgFile, csvFile);

function [P, yLimits] = prepare_case_sensor_waveform_local(bundle, fit, sid, maxPoints, yLimits)
idxSensor = find(bundle.sensorIds == sid, 1, 'first');
if isempty(idxSensor)
    error('CH%d not found in bundle.', sid);
end
idxAll = find(bundle.sensorIndex == idxSensor);
idxPlot = decimate_index_local(idxAll, maxPoints);
xComp = bundle.X(idxPlot) - fit.dxMm - fit.uMm(idxPlot);
vObs = bundle.V(idxPlot);
vPred = fit.VPred(idxPlot);
[xSorted, order] = sort(xComp(:));
P = struct();
P.xObserved = xSorted;
P.vObserved = vObs(order);
P.xPred = xSorted;
P.vPred = vPred(order);
P.rawPointCount = numel(idxAll);
allY = [P.vObserved(:); P.vPred(:)];
allY = allY(isfinite(allY));
if ~isempty(allY)
    yLimits(1) = min(yLimits(1), min(allY));
    yLimits(2) = max(yLimits(2), max(allY));
end
if ~all(isfinite(yLimits)) || diff(yLimits) <= 0
    yLimits = [-1, 1];
end
end

function idx = decimate_index_local(idx, maxCount)
if numel(idx) <= maxCount
    return;
end
idx = idx(unique(round(linspace(1, numel(idx), maxCount))));
end

function rows = append_curve_local(rows, windowId, sensorId, caseName, curveName, x, y)
n = numel(x);
newRows = repmat(struct('window_id', windowId, 'sensor_id', sensorId, ...
    'case_name', string(caseName), 'curve_name', string(curveName), ...
    'point_index', NaN, 'x_mm', NaN, 'voltage_mV', NaN), n, 1);
for i = 1:n
    newRows(i).point_index = i;
    newRows(i).x_mm = x(i);
    newRows(i).voltage_mV = y(i);
end
rows{end+1, 1} = newRows;
end

function style = paper_style_local()
style = struct();
style.blue = [0.000, 0.447, 0.698];
style.orange = [0.835, 0.369, 0.000];
style.gray = [0.45, 0.45, 0.45];
style.fontName = 'Arial';
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'TickDir', 'out');
end

function value = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 0
    error('%s must be a nonnegative number.', name);
end
value = floor(value);
end

function values = parse_integer_list_env_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    values = [];
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    error('%s must contain integer IDs separated by spaces.', name);
end
end
