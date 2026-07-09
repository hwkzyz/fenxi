%% Visualize Step07J waveform reconstruction pipeline, 20251222
% This script reads saved Step07J results and overlays four waveform levels:
%   1) low-speed rotating template;
%   2) low-speed template reconstructed by the calibrated gap library;
%   3) high-speed samples after vibration compensation;
%   4) final Step07J reconstructed waveform from the main model.
%
% Optional environment variables:
%   STEP07J_PIPELINE_RESULT_FILE   explicit Step07J .mat file for one run
%   STEP07J_PIPELINE_WINDOW        window index, default Result.BestWindowIndex
%   STEP07J_PIPELINE_SENSORS       sensor ids separated by spaces, e.g. "1 3 6"
%   STEP07J_PIPELINE_MAX_POINTS    plotted high-speed points per sensor

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C = current_case_local(thisDir);
explicitResultFile = strtrim(getenv('STEP07J_PIPELINE_RESULT_FILE'));
windowOverride = parse_positive_integer_env_local('STEP07J_PIPELINE_WINDOW', 0);
sensorOverride = parse_integer_list_env_local('STEP07J_PIPELINE_SENSORS');
maxPointsPerSensor = parse_positive_integer_env_local('STEP07J_PIPELINE_MAX_POINTS', 700);
if ~isempty(explicitResultFile)
    C.resultFile = explicitResultFile;
end

make_case_figure_local(C, windowOverride, sensorOverride, maxPointsPerSensor);

function make_case_figure_local(C, windowOverride, sensorOverride, maxPointsPerSensor)
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
mainModel = resolve_main_model_local(Result);
if ~isfield(wr.modelFits, mainModel)
    error('Window %d does not contain model fit "%s".', windowId, mainModel);
end
fit = wr.modelFits.(mainModel);
if isempty(sensorOverride)
    sensorIds = bundle.sensorIds(:).';
else
    sensorIds = sensorOverride(:).';
end

figDir = fullfile(C.caseDir, 'outputs', 'figures_step07j_waveform_reconstruction_pipeline');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end
style = paper_style_local();
fig = figure('Name', sprintf('Step07J reconstruction pipeline %s W%d', C.dataset, windowId), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 17.5, max(7.0, 4.2 * numel(sensorIds))]);
tiledlayout(numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

curveRows = {};
fprintf('\n=== Step07J waveform reconstruction pipeline ===\n');
fprintf('Dataset: %s, blade/sensors: %s\n', C.dataset, C.bladeTag);
fprintf('Result: %s\n', C.resultFile);
fprintf('Window: %d, main model: %s\n', windowId, mainModel);

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    idxSensor = find(bundle.sensorIds == sid, 1, 'first');
    if isempty(idxSensor)
        warning('CH%d not found in bundle. Skipping.', sid);
        continue;
    end
    maskAll = bundle.sensorIndex == idxSensor;
    idxAll = find(maskAll);
    idxPlot = decimate_index_local(idxAll, maxPointsPerSensor);
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);

    xLow = Tpl.x_grid(:);
    vLow = (Tpl.v_grid(:) - Tpl.baseline) * 1000;
    xGapLow = corr.x(:);
    vGapLow = corr.vFitMv(:);

    xRaw = bundle.X(idxPlot) - fit.dxMm;
    uPlot = fit.uMm(idxPlot);
    xComp = xRaw - uPlot;
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
    plot(xLow, vLow, '-', 'Color', style.black, 'LineWidth', 1.1, ...
        'DisplayName', 'low-speed template');
    plot(xGapLow, vGapLow, '-', 'Color', style.blue, 'LineWidth', 1.2, ...
        'DisplayName', 'low + gap library');
    plot(xCompSort, vObsSort, '.', 'Color', style.gray, 'MarkerSize', 3.7, ...
        'DisplayName', 'high observed, vib-comp');
    plot(xDense, vFinalLine, '-', 'Color', style.green, 'LineWidth', 1.45, ...
        'DisplayName', sprintf('final %s', mainModel));
    plot(xCompSort, vPredSort, '.', 'Color', style.greenDark, 'MarkerSize', 3.2, ...
        'DisplayName', 'final prediction samples');

    title(sprintf('%s CH%d W%d | EO%d, f %.2f Hz, A %.4f mm, RMSE %.2f mV', ...
        C.dataset, sid, windowId, fit.EO, fit.freqHz, fit.amplitudeMm, fit.weightedRmseMv), ...
        'Interpreter', 'none');
    xlabel('Equivalent x after vibration compensation (mm)', 'Interpreter', 'none');
    ylabel('Voltage (mV)');
    if is == 1
        legend('Location', 'best', 'Box', 'off', 'Interpreter', 'none');
    end
    format_axes_local(gca, style);

    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, "low_template", xLow, vLow);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, "low_gap_library_fit", xGapLow, vGapLow);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, "high_observed_vib_comp", xCompSort, vObsSort);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, "final_reconstruction_line", xDense, vFinalLine);
    curveRows = append_curve_local(curveRows, C.dataset, windowId, sid, "final_prediction_samples", xCompSort, vPredSort);

    fprintf('  CH%d: high points %d -> plotted %d, dg %.4f, dmu %.5f, dtau %.4f\n', ...
        sid, numel(idxAll), numel(idxPlot), dg, dmu, dtau);
end

[~, resultBase] = fileparts(C.resultFile);
figBase = sprintf('Step07J_WaveformReconstructionPipeline_%s_W%02d', resultBase, windowId);
pngFile = fullfile(figDir, [figBase, '.png']);
pdfFile = fullfile(figDir, [figBase, '.pdf']);
csvFile = fullfile(figDir, [figBase, '_Curves.csv']);
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
if ~isempty(curveRows)
    T = struct2table(vertcat(curveRows{:}));
    writetable(T, csvFile);
end
fprintf('Saved:\n  %s\n  %s\n  %s\n', pngFile, pdfFile, csvFile);
end

function C = current_case_local(thisDir)
C = CaseConfig();
C.bladeTag = C.caseTag;
C.caseDir = thisDir;
C.resultFile = resolve_step07j_result_local(fullfile(C.caseDir, 'outputs'), C);
end

function resultFile = resolve_step07j_result_local(outDir, C)
regionTag = '';
if isfield(C, 'flowConfig') && isfield(C.flowConfig, 'identification') && ...
        isfield(C.flowConfig.identification, 'resonanceRegionShortTag')
    regionTag = C.flowConfig.identification.resonanceRegionShortTag;
end

candidates = {
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_%s_main_gaptilt.mat', C.dataset, C.caseTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', C.dataset, C.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', C.dataset, C.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_%s_main_gaptilt.mat', C.dataset, C.bladeId, C.sensorTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaptilt.mat', C.dataset, C.bladeId, C.sensorTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_method_compare.mat', C.dataset, C.bladeId, C.sensorTag)
    };
candidates = candidates(~cellfun(@(s) contains(s, '__'), candidates));
for i = 1:numel(candidates)
    candidate = fullfile(outDir, candidates{i});
    if isfile(candidate)
        resultFile = candidate;
        return;
    end
end

files = dir(fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s*.mat', C.dataset, C.caseTag)));
if isempty(files)
    files = dir(fullfile(outDir, sprintf( ...
        'Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s*.mat', ...
        C.dataset, C.bladeId, C.sensorTag)));
end
if ~isempty(files)
    [~, idx] = max([files.datenum]);
    resultFile = fullfile(files(idx).folder, files(idx).name);
    warning('Using latest matching Step07J result because canonical files were not found:\n  %s', resultFile);
    return;
end

error(['Missing Step07J result for %s. Expected a main_gaptilt or method_compare file under:\n' ...
    '  %s\nRun Step07J/Step00 first.'], C.caseTag, outDir);
end

function mainModel = resolve_main_model_local(Result)
if isfield(Result, 'MainModel') && ~isempty(Result.MainModel)
    mainModel = char(Result.MainModel);
elseif isfield(Result, 'cfg') && isfield(Result.cfg, 'mainModel') && ~isempty(Result.cfg.mainModel)
    mainModel = char(Result.cfg.mainModel);
else
    mainModel = 'gap_tilt';
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

function rows = append_curve_local(rows, dataset, windowId, sensorId, curveName, x, v)
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
    'sensorId', sensorId, 'curve', string(curveName), 'xMm', NaN, 'vMv', NaN), n, 1);
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
n = str2double(raw);
if ~isfinite(n) || n <= 0
    error('%s must be a positive number.', name);
end
n = round(n);
end

function values = parse_integer_list_env_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    values = [];
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    error('%s must contain integer sensor ids.', name);
end
end

function style = paper_style_local()
style = struct();
style.fontName = 'Arial';
style.fontSize = 8.5;
style.black = [0.05 0.05 0.05];
style.gray = [0.50 0.50 0.50];
style.blue = [0.10 0.32 0.68];
style.green = [0.00 0.55 0.28];
style.greenDark = [0.00 0.34 0.17];
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.fontSize, ...
    'LineWidth', 0.8, 'TickDir', 'out');
end

