%% Visualize Step07J method-specific waveform reconstruction, 20250527
% Make one figure per method: direct/fixed, gap_only, and gap_tilt.
% Each figure shows calibrated low-speed waveforms, vibration-compensated
% high-speed samples, and the method-specific reconstructed waveform.
%
% Optional environment variables:
%   STEP07J_METHOD_RESULT_FILE   explicit Step07J .mat file for one run
%   STEP07J_METHOD_WINDOW        window index, default Result.BestWindowIndex
%   STEP07J_METHOD_SENSORS       sensor ids separated by spaces, e.g. "1 3 6"
%   STEP07J_METHOD_MAX_POINTS    plotted high-speed points per sensor

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C = current_case_local(thisDir);
explicitResultFile = strtrim(getenv('STEP07J_METHOD_RESULT_FILE'));
windowOverride = parse_positive_integer_env_local('STEP07J_METHOD_WINDOW', 0);
sensorOverride = parse_integer_list_env_local('STEP07J_METHOD_SENSORS');
maxPointsPerSensor = parse_positive_integer_env_local('STEP07J_METHOD_MAX_POINTS', 900);
if ~isempty(explicitResultFile)
    C.resultFile = explicitResultFile;
end

make_method_figures_local(C, windowOverride, sensorOverride, maxPointsPerSensor);

function make_method_figures_local(C, windowOverride, sensorOverride, maxPointsPerSensor)
S = load(C.resultFile, 'Result');
Result = S.Result;
if windowOverride > 0
    windowId = windowOverride;
elseif isfield(Result, 'BestWindowIndex') && isfinite(Result.BestWindowIndex)
    windowId = Result.BestWindowIndex;
else
    windowId = 1;
end
if windowId < 1 || windowId > numel(Result.WindowResult)
    error('Window %d is outside result range 1-%d.', windowId, numel(Result.WindowResult));
end

wr = Result.WindowResult(windowId);
bundle = wr.bundle;
if isempty(sensorOverride)
    sensorIds = bundle.sensorIds(:).';
else
    sensorIds = sensorOverride(:).';
end

figDir = fullfile(C.caseDir, 'outputs', 'figures_step07j_method_separate_pipeline');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

methods = struct( ...
    'field', {'fixed','gap_only','gap_tilt'}, ...
    'label', {'direct/fixed','gap only','gap + tilt'});

fprintf('\n=== Step07J method-specific waveform pipeline ===\n');
fprintf('Dataset: %s, blade/sensors: %s\n', C.dataset, C.bladeTag);
fprintf('Result: %s\n', C.resultFile);
fprintf('Window: %d\n', windowId);

for im = 1:numel(methods)
    methodField = methods(im).field;
    if ~isfield(wr.modelFits, methodField)
        warning('Window %d does not contain model fit "%s". Skipping.', windowId, methodField);
        continue;
    end
    fit = wr.modelFits.(methodField);
    [fig, curveRows] = draw_one_method_local(C, wr, bundle, fit, methods(im), ...
        windowId, sensorIds, maxPointsPerSensor);

    [~, resultBase] = fileparts(C.resultFile);
    figBase = sprintf('Step07J_MethodSeparate_%s_W%02d_%s', resultBase, windowId, methodField);
    pngFile = fullfile(figDir, [figBase, '.png']);
    pdfFile = fullfile(figDir, [figBase, '.pdf']);
    csvFile = fullfile(figDir, [figBase, '_Curves.csv']);
    exportgraphics(fig, pngFile, 'Resolution', 300);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
    if ~isempty(curveRows)
        T = struct2table(vertcat(curveRows{:}));
        writetable(T, csvFile);
    end
    fprintf('Saved %s:\n  %s\n  %s\n  %s\n', methodField, pngFile, pdfFile, csvFile);
end
end

function [fig, curveRows] = draw_one_method_local(C, ~, bundle, fit, methodInfo, ...
    windowId, sensorIds, maxPointsPerSensor)
style = paper_style_local();
fig = figure('Name', sprintf('%s %s W%d', C.dataset, methodInfo.label, windowId), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 17.5, max(7.0, 4.2 * numel(sensorIds))]);
tiledlayout(numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
curveRows = {};

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    idxSensor = find(bundle.sensorIds == sid, 1, 'first');
    if isempty(idxSensor)
        warning('CH%d not found in bundle. Skipping.', sid);
        continue;
    end
    idxAll = find(bundle.sensorIndex == idxSensor);
    idxPlot = decimate_index_local(idxAll, maxPointsPerSensor);
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);

    xLow = Tpl.x_grid(:);
    vLow = (Tpl.v_grid(:) - Tpl.baseline) * 1000;
    xGapLow = corr.x(:);
    vGapLow = corr.vFitMv(:);

    xRaw = bundle.X(idxPlot) - fit.dxMm;
    xComp = xRaw - fit.uMm(idxPlot);
    vObs = bundle.V(idxPlot);
    vPred = fit.VPred(idxPlot);
    [xCompSort, orderComp] = sort(xComp(:));
    vObsSort = vObs(orderComp);
    vPredSort = vPred(orderComp);

    dg = fit.deltaGapMm(idxSensor);
    dmu = fit.deltaMuGapPerXMm(idxSensor);
    dtau = fit.deltaTauMm(idxSensor);
    xDense = build_dense_x_local(xLow, xCompSort);
    [vFinalLine, ~] = eval_static_warp_template_local(Tpl, bundle.responseSurface, corr, dg, dmu, dtau, xDense);

    nexttile; hold on; grid on; box on;
    plot(xCompSort, vObsSort, '.', 'Color', style.gray, 'MarkerSize', 3.4, ...
        'DisplayName', 'high-speed data');
    plot(xLow, vLow, '--', 'Color', style.black, 'LineWidth', 1.0, ...
        'DisplayName', 'low-speed template');
    plot(xGapLow, vGapLow, '-', 'Color', style.blue, 'LineWidth', 1.15, ...
        'DisplayName', 'static gap library');
    plot(xCompSort, vPredSort, '.', 'Color', style.greenDark, 'MarkerSize', 2.8, ...
        'HandleVisibility', 'off');
    plot(xDense, vFinalLine, '-', 'Color', style.green, 'LineWidth', 1.55, ...
        'DisplayName', sprintf('%s final', methodInfo.label));

    title(sprintf('%s CH%d W%d | %s, EO%d, f %.2f Hz, A %.4f mm, RMSE %.2f mV, dg %.4f, dmu %.5f', ...
        C.dataset, sid, windowId, methodInfo.label, fit.EO, fit.freqHz, ...
        fit.amplitudeMm, fit.weightedRmseMv, dg, dmu), 'Interpreter', 'none');
    xlabel('Equivalent x after vibration compensation (mm)', 'Interpreter', 'none');
    ylabel('Voltage (mV)');
    if is == 1
        legend('Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
    end
    format_axes_local(gca, style);

    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, methodInfo.field, "low_template", xLow, vLow);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, methodInfo.field, "low_gap_library_fit", xGapLow, vGapLow);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, methodInfo.field, "high_observed_vib_comp", xCompSort, vObsSort);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, methodInfo.field, "method_reconstruction_line", xDense, vFinalLine);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, methodInfo.field, "method_prediction_samples", xCompSort, vPredSort);

    fprintf('  %-8s CH%d: points %d -> %d, EO%d, A %.4f, RMSE %.2f, dg %.4f, dmu %.5f\n', ...
        methodInfo.field, sid, numel(idxAll), numel(idxPlot), fit.EO, fit.amplitudeMm, ...
        fit.weightedRmseMv, dg, dmu);
end
end

function C = current_case_local(thisDir)
C = CaseConfig();
C.bladeTag = C.caseTag;
C.caseDir = thisDir;
C.resultFile = fullfile(C.caseDir, 'outputs', ...
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', C.dataset, C.caseTag));
if ~isfile(C.resultFile)
    legacyFile = fullfile(C.caseDir, 'outputs', ...
        sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_method_compare.mat', ...
        C.dataset, C.bladeId, C.sensorTag));
    if isfile(legacyFile)
        C.resultFile = legacyFile;
    else
        error(['Missing Step07J comparison result file: %s\n' ...
            'Legacy fallback also missing: %s\n' ...
            'Run Step09_Compare_Step07J_Methods_20250527 first, or set STEP07J_METHOD_RESULT_FILE.'], ...
            C.resultFile, legacyFile);
    end
end
end

function xDense = build_dense_x_local(xLow, xHigh)
lo = max(min(xLow(:)), prctile(xHigh(isfinite(xHigh)), 1));
hi = min(max(xLow(:)), prctile(xHigh(isfinite(xHigh)), 99));
if ~isfinite(lo) || ~isfinite(hi) || lo >= hi
    lo = min(xLow(:));
    hi = max(xLow(:));
end
xDense = linspace(lo, hi, 700).';
end

function rows = append_curve_local(rows, dataset, windowId, sensorId, method, curveName, x, v)
x = x(:);
v = v(:);
mask = isfinite(x) & isfinite(v);
x = x(mask);
v = v(mask);
n = numel(x);
if n == 0
    return;
end
newRows = repmat(struct('dataset', string(dataset), 'windowId', windowId, ...
    'sensorId', sensorId, 'method', string(method), 'curve', string(curveName), ...
    'xMm', NaN, 'vMv', NaN), n, 1);
for i = 1:n
    newRows(i).xMm = x(i);
    newRows(i).vMv = v(i);
end
rows{end+1, 1} = newRows;
end

function idx = decimate_index_local(idx, maxCount)
if isinf(maxCount) || numel(idx) <= maxCount
    return;
end
idx = idx(unique(round(linspace(1, numel(idx), maxCount))));
end

function [model, overshoot] = eval_static_warp_template_local(Tpl, responseSurface, corr, dg, dmu, dtau, xOpr)
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr(:) - tau);
gEffRaw = g0 + mu .* (xOpr(:) - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
xLib = min(max(xLibRaw, xMin), xMax);
gEff = min(max(gEffRaw, gMin), gMax);
overshoot = hypot(max(xMin - xLibRaw, 0) + max(xLibRaw - xMax, 0), ...
    max(gMin - gEffRaw, 0) + max(gEffRaw - gMax, 0));
model = eval_response_surface_local(responseSurface, gEff, xLib);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1);
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function n = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    n = defaultValue;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp < 0
    error('%s must be a nonnegative integer.', name);
end
n = floor(tmp);
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

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 8.5;
style.black = [0.00, 0.00, 0.00];
style.blue = [0.05, 0.26, 0.62];
style.green = [0.00, 0.48, 0.30];
style.greenDark = [0.00, 0.30, 0.18];
style.gray = [0.62, 0.62, 0.62];
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', 0.7, 'TickDir', 'out', 'Layer', 'top');
end
