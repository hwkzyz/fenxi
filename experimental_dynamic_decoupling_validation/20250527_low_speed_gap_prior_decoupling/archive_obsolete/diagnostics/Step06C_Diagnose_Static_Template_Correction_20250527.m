%% Step06C: Diagnose low-speed-template plus gap-library static correction
% This step does not identify vibration.  It checks whether the static gap
% library can improve the high-speed no-vibration template built from the
% measured low-speed template:
%
%   M_s(g,x) = T_low,s(x) + F(g, x - dx_lib,s) - F(g0_s, x - dx_lib,s)
%
% No affine gain or voltage offset is fitted in this diagnostic.  The goal
% is to verify the template construction before adding VP or waveform
% vibration identification.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06c_static_template_correction');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.anchorGapGridCount = 41;
cfg.anchorDxGridMm = -0.80:0.02:0.80;
cfg.staticDxGridMm = -0.25:0.01:0.25;
cfg.deltaGapGridMm = -0.35:0.025:0.35;
cfg.lowResidualWeightGrid = 0:0.10:1.00;
cfg.highFitMaxPointsPerSensor = 260;
cfg.minFitPoints = 30;

sensorOverride = strtrim(getenv('STEP06C_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06C_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
templateDir = fullfile(rootDir, '20250527_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(highMapFile)
    error('Run Step06 first. Missing file: %s', highMapFile);
end
if ~isfile(templateFile)
    error('Missing low-speed template file: %s', templateFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
H = load(highMapFile, 'highMap');
highMap = H.highMap;
St = load(templateFile, 'Template');
Template = St.Template;

if isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(highMap.t_v)
    highMap.x_model_v = highMap.x_fit_v(:);
else
    highMap.x_model_v = highMap.x_v(:);
end

fprintf('\n=== Step06C: static template correction diagnostic ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('HighMap: %s\n', highMapFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Sensors: %s\n', mat2str(cfg.analysisSensors));

%% 2. Anchor each low-speed template to the gap library without affine terms
anchorRows = cell(numel(cfg.analysisSensors), 1);
anchor = repmat(struct('sensorId', NaN, 'g0Mm', NaN, 'dxLibMm', NaN, ...
    'rmseMv', NaN, 'pointCount', NaN, 'lowX', [], 'lowV', [], 'fitV', []), ...
    numel(cfg.analysisSensors), 1);
gGrid = linspace(min(responseSurface.gTrainMm), max(responseSurface.gTrainMm), ...
    cfg.anchorGapGridCount);

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    [g0, dxLib, rmse, cache] = fit_low_template_to_gap_library_no_affine_local( ...
        Tsen, responseSurface, gGrid, cfg.anchorDxGridMm, cfg.minFitPoints);
    anchor(is).sensorId = sid;
    anchor(is).g0Mm = g0;
    anchor(is).dxLibMm = dxLib;
    anchor(is).rmseMv = rmse;
    anchor(is).pointCount = numel(cache.x);
    anchor(is).lowX = cache.x;
    anchor(is).lowV = cache.vLow;
    anchor(is).fitV = cache.vFit;
    anchorRows{is} = table(sid, g0, dxLib, rmse, numel(cache.x), ...
        'VariableNames', {'sensorId', 'g0Mm', 'dxLibMm', 'lowAnchorRmseMv', 'pointCount'});
    fprintf('  CH%d low-template anchor: g0 %.4f mm, dx_lib %.4f mm, RMSE %.3f mV\n', ...
        sid, g0, dxLib, rmse);
end
anchorTable = vertcat(anchorRows{:});

%% 3. Static high-speed comparison: low template vs gap-only vs corrected template
windowSpecs = build_sliding_window_specs_local(highMap);
compareRows = {};
bestCache = struct([]);
rowIdx = 0;
for iw = 1:numel(windowSpecs)
    W = subset_highmap_by_laps_local(highMap, windowSpecs(iw).lapRange);
    for is = 1:numel(cfg.analysisSensors)
        sid = cfg.analysisSensors(is);
        Tsen = get_template_sensor_local(Template, sid);
        A = anchor([anchor.sensorId] == sid);
        sensorData = extract_sensor_static_points_local(W, sid, cfg.highFitMaxPointsPerSensor);
        if numel(sensorData.x) < cfg.minFitPoints
            continue;
        end

        lowFit = fit_high_static_low_template_local(sensorData, Tsen, cfg.staticDxGridMm);
        gapFit = fit_high_static_gap_only_local(sensorData, responseSurface, ...
            gGrid, cfg.staticDxGridMm);
        corrFit = fit_high_static_corrected_template_local(sensorData, Tsen, ...
            responseSurface, A.g0Mm, A.dxLibMm, cfg.deltaGapGridMm, cfg.staticDxGridMm);
        blendFit = fit_high_static_gap_plus_low_residual_local(sensorData, A, ...
            responseSurface, gGrid, cfg.staticDxGridMm, cfg.lowResidualWeightGrid);

        rowIdx = rowIdx + 1;
        compareRows{rowIdx, 1} = make_compare_row_local(iw, windowSpecs(iw).lapRange, ...
            sid, 'low_template_only', NaN, lowFit.dxMm, NaN, lowFit.rmseMv, lowFit.pointCount);
        rowIdx = rowIdx + 1;
        compareRows{rowIdx, 1} = make_compare_row_local(iw, windowSpecs(iw).lapRange, ...
            sid, 'gap_library_only', gapFit.gMm, gapFit.dxMm, NaN, gapFit.rmseMv, gapFit.pointCount);
        rowIdx = rowIdx + 1;
        compareRows{rowIdx, 1} = make_compare_row_local(iw, windowSpecs(iw).lapRange, ...
            sid, 'low_template_plus_gap_delta', corrFit.gMm, corrFit.dxMm, NaN, corrFit.rmseMv, corrFit.pointCount);
        rowIdx = rowIdx + 1;
        compareRows{rowIdx, 1} = make_compare_row_local(iw, windowSpecs(iw).lapRange, ...
            sid, 'gap_library_plus_low_residual', blendFit.gMm, blendFit.dxMm, ...
            blendFit.lowResidualWeight, blendFit.rmseMv, blendFit.pointCount);

        if iw == 1
            bestCache(end+1).sensorId = sid; %#ok<SAGROW>
            bestCache(end).sensorData = sensorData;
            bestCache(end).lowFit = lowFit;
            bestCache(end).gapFit = gapFit;
            bestCache(end).corrFit = corrFit;
            bestCache(end).blendFit = blendFit;
            bestCache(end).anchor = A;
        end
    end
end
compareTable = vertcat(compareRows{:});
summaryTable = summarize_static_compare_local(compareTable);

%% 4. Save diagnostics
matFile = fullfile(outDir, sprintf('Step06C_Static_Template_Correction_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
anchorCsv = fullfile(outDir, sprintf('Step06C_LowTemplate_GapLibrary_Anchor_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
compareCsv = fullfile(outDir, sprintf('Step06C_Static_Template_Correction_Comparison_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
summaryCsv = fullfile(outDir, sprintf('Step06C_Static_Template_Correction_Summary_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figAnchor = fullfile(figDir, sprintf('Step06C_LowTemplate_Anchor_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figHigh = fullfile(figDir, sprintf('Step06C_HighStatic_Window01_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'cfg', 'responseFile', 'highMapFile', 'templateFile', ...
    'anchorTable', 'compareTable', 'summaryTable', 'anchor', 'bestCache', '-v7.3');
writetable(anchorTable, anchorCsv);
writetable(compareTable, compareCsv);
writetable(summaryTable, summaryCsv);
plot_anchor_fit_local(anchor, figAnchor);
plot_high_static_window_local(bestCache, figHigh);

fprintf('\nStep06C complete.\n');
disp(anchorTable);
disp(summaryTable);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', matFile, anchorCsv, compareCsv, summaryCsv);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', targetBlade, exactTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20250527.mat', targetBlade)));
if isempty(files)
    files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_GradientXRange030_20250527.mat', targetBlade)));
end
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)(?:_GradientXRange030)?_20250527\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    available = arrayfun(@(c) str2double(c), token{1});
    if all(ismember(sensorIds, available))
        templateFile = fullfile(files(i).folder, files(i).name);
        return;
    end
end
error('No low-speed template file found under %s for sensors %s.', templateDir, mat2str(sensorIds));
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function [gBest, dxBest, rmseBest, cache] = fit_low_template_to_gap_library_no_affine_local( ...
    Tsen, responseSurface, gGrid, dxGrid, minFitPoints)
x = Tsen.x_grid(:);
vLow = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
mask = isfinite(x) & isfinite(vLow);
mask = mask & x >= min(responseSurface.xGrid) + max(abs(dxGrid)) & ...
    x <= max(responseSurface.xGrid) - max(abs(dxGrid));
x = x(mask);
vLow = vLow(mask);
if numel(x) < minFitPoints
    error('CH%d has too few low-template points for anchoring.', Tsen.sensor_id);
end
rmseBest = inf;
gBest = NaN;
dxBest = NaN;
vBest = nan(size(vLow));
for g = gGrid(:).'
    for dx = dxGrid(:).'
        vFit = eval_response_surface_local(responseSurface, g, x - dx);
        rmse = rmse_no_affine_local(vLow, vFit, ones(size(vLow)));
        if rmse < rmseBest
            rmseBest = rmse;
            gBest = g;
            dxBest = dx;
            vBest = vFit;
        end
    end
end
cache = struct('x', x, 'vLow', vLow, 'vFit', vBest);
end

function sensorData = extract_sensor_static_points_local(W, sid, maxPoints)
mask = W.S_v == sid & isfinite(W.x_model_v) & isfinite(W.V_a) & isfinite(W.W_v);
idx = find(mask);
if isempty(idx)
    sensorData = struct('x', [], 'v', [], 'w', []);
    return;
end
if numel(idx) > maxPoints
    pick = round(linspace(1, numel(idx), maxPoints));
    idx = idx(pick);
end
sensorData = struct('x', W.x_model_v(idx), 'v', W.V_a(idx), 'w', max(W.W_v(idx), 0.05));
end

function fit = fit_high_static_low_template_local(sensorData, Tsen, dxGrid)
vTpl = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
fit = init_fit_result_local();
for dx = dxGrid(:).'
    pred = interp1(Tsen.x_grid(:), vTpl, sensorData.x - dx, 'pchip', NaN);
    rmse = rmse_no_affine_local(sensorData.v, pred, sensorData.w);
    if rmse < fit.rmseMv
        fit = struct('gMm', NaN, 'dxMm', dx, 'deltaGapMm', NaN, ...
            'rmseMv', rmse, 'pointCount', nnz(isfinite(pred)), 'pred', pred);
    end
end
end

function fit = fit_high_static_gap_only_local(sensorData, responseSurface, gGrid, dxGrid)
fit = init_fit_result_local();
for g = gGrid(:).'
    for dx = dxGrid(:).'
        pred = eval_response_surface_local(responseSurface, g, sensorData.x - dx);
        rmse = rmse_no_affine_local(sensorData.v, pred, sensorData.w);
        if rmse < fit.rmseMv
            fit = struct('gMm', g, 'dxMm', dx, 'deltaGapMm', NaN, ...
                'rmseMv', rmse, 'pointCount', nnz(isfinite(pred)), 'pred', pred);
        end
    end
end
end

function fit = fit_high_static_corrected_template_local(sensorData, Tsen, responseSurface, g0, dxLib, deltaGrid, dxGrid)
vTpl = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
fit = init_fit_result_local();
for dg = deltaGrid(:).'
    g = min(max(g0 + dg, min(responseSurface.gTrainMm)), max(responseSurface.gTrainMm));
    for dx = dxGrid(:).'
        xLow = sensorData.x - dx;
        low = interp1(Tsen.x_grid(:), vTpl, xLow, 'pchip', NaN);
        corr = eval_response_surface_local(responseSurface, g, xLow - dxLib) - ...
            eval_response_surface_local(responseSurface, g0, xLow - dxLib);
        pred = low + corr;
        rmse = rmse_no_affine_local(sensorData.v, pred, sensorData.w);
        if rmse < fit.rmseMv
            fit = struct('gMm', g, 'dxMm', dx, 'deltaGapMm', dg, ...
                'rmseMv', rmse, 'pointCount', nnz(isfinite(pred)), 'pred', pred);
        end
    end
end
end

function fit = fit_high_static_gap_plus_low_residual_local(sensorData, anchor, responseSurface, gGrid, dxGrid, lambdaGrid)
residualLow = anchor.lowV(:) - anchor.fitV(:);
fit = init_fit_result_local();
for g = gGrid(:).'
    for dx = dxGrid(:).'
        xLow = sensorData.x - dx;
        base = eval_response_surface_local(responseSurface, g, xLow - anchor.dxLibMm);
        residual = interp1(anchor.lowX(:), residualLow, xLow, 'pchip', NaN);
        for lambda = lambdaGrid(:).'
            pred = base + lambda .* residual;
            rmse = rmse_no_affine_local(sensorData.v, pred, sensorData.w);
            if rmse < fit.rmseMv
                fit = struct('gMm', g, 'dxMm', dx, 'deltaGapMm', NaN, ...
                    'lowResidualWeight', lambda, 'rmseMv', rmse, ...
                    'pointCount', nnz(isfinite(pred)), 'pred', pred);
            end
        end
    end
end
end

function fit = init_fit_result_local()
fit = struct('gMm', NaN, 'dxMm', NaN, 'deltaGapMm', NaN, ...
    'lowResidualWeight', NaN, 'rmseMv', inf, 'pointCount', 0, 'pred', []);
end

function rmse = rmse_no_affine_local(v, pred, w)
valid = isfinite(v) & isfinite(pred) & isfinite(w);
if nnz(valid) < 3
    rmse = inf;
    return;
end
rmse = sqrt(sum(w(valid) .* (v(valid) - pred(valid)).^2) / max(sum(w(valid)), eps));
end

function row = make_compare_row_local(windowId, lapRange, sid, methodName, gMm, dxMm, lowResidualWeight, rmseMv, pointCount)
row = table(windowId, lapRange(1), lapRange(end), sid, string(methodName), ...
    gMm, dxMm, lowResidualWeight, rmseMv, pointCount, ...
    'VariableNames', {'window_id', 'lap_start', 'lap_end', 'sensorId', ...
    'method', 'gMm', 'dxMm', 'lowResidualWeight', 'rmseMv', 'pointCount'});
end

function summaryTable = summarize_static_compare_local(compareTable)
methods = unique(compareTable.method, 'stable');
rows = cell(numel(methods), 1);
for i = 1:numel(methods)
    mask = compareTable.method == methods(i);
    rows{i} = table(methods(i), nnz(mask), ...
        mean(compareTable.rmseMv(mask), 'omitnan'), ...
        median(compareTable.rmseMv(mask), 'omitnan'), ...
        min(compareTable.rmseMv(mask), [], 'omitnan'), ...
        max(compareTable.rmseMv(mask), [], 'omitnan'), ...
        'VariableNames', {'method', 'rowCount', 'meanRmseMv', ...
        'medianRmseMv', 'minRmseMv', 'maxRmseMv'});
end
summaryTable = vertcat(rows{:});
end

function windowSpecs = build_sliding_window_specs_local(highMap)
laps = unique(highMap.rev_v(:).');
laps = sort(laps(isfinite(laps)));
winSize = highMap.analysisWinSize;
step = highMap.slidingStep;
nWin = floor((numel(laps) - winSize) / step) + 1;
windowSpecs = repmat(struct('window_id', NaN, 'lapRange', []), nWin, 1);
for iw = 1:nWin
    startIdx = 1 + (iw - 1) * step;
    windowSpecs(iw).window_id = iw;
    windowSpecs(iw).lapRange = laps(startIdx:startIdx + winSize - 1);
end
end

function W = subset_highmap_by_laps_local(highMap, lapRange)
mask = ismember(highMap.rev_v, lapRange);
W = highMap;
fields = {'t_v','x_v','x_fit_v','x_model_v','x_timebased_ref_v','x_mm_v', ...
    'x_peak_mm_v','V_raw_mV','V_a','baseline_mV_v','S_v','rev_v','W_v', ...
    'theta_v','pulsePeakTime_v'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(W, f) && numel(W.(f)) == numel(mask)
        W.(f) = W.(f)(mask);
    end
end
W.pointCount = nnz(mask);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function plot_anchor_fit_local(anchor, figFile)
fig = figure('Name', 'Step06C low-template gap-library anchor', 'Color', 'w', ...
    'Position', [80, 80, 1300, 310 * numel(anchor)]);
tiledlayout(numel(anchor), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(anchor)
    nexttile; hold on; grid on; box on;
    plot(anchor(i).lowX, anchor(i).lowV, '.', 'MarkerSize', 5, ...
        'DisplayName', 'T low');
    plot(anchor(i).lowX, anchor(i).fitV, '-', 'LineWidth', 1.5, ...
        'DisplayName', 'F(g0,x-dx)');
    title(sprintf('CH%d | g0 %.4f mm | dx %.4f mm | RMSE %.2f mV', ...
        anchor(i).sensorId, anchor(i).g0Mm, anchor(i).dxLibMm, anchor(i).rmseMv));
    xlabel('x (mm)');
    ylabel('mV');
    legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_high_static_window_local(bestCache, figFile)
if isempty(bestCache)
    return;
end
fig = figure('Name', 'Step06C high-speed static template comparison', 'Color', 'w', ...
    'Position', [80, 80, 1350, 330 * numel(bestCache)]);
tiledlayout(numel(bestCache), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(bestCache)
    D = bestCache(i).sensorData;
    [xSort, order] = sort(D.x);
    nexttile; hold on; grid on; box on;
    plot(D.x, D.v, '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 5, ...
        'DisplayName', 'high-speed samples');
    plot(xSort, bestCache(i).lowFit.pred(order), '-', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('T low %.1f mV', bestCache(i).lowFit.rmseMv));
    plot(xSort, bestCache(i).gapFit.pred(order), '-', 'LineWidth', 1.1, ...
        'DisplayName', sprintf('F(g) %.1f mV', bestCache(i).gapFit.rmseMv));
    plot(xSort, bestCache(i).corrFit.pred(order), '-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T low + deltaF %.1f mV', bestCache(i).corrFit.rmseMv));
    plot(xSort, bestCache(i).blendFit.pred(order), '-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('F(g)+lambda R %.1f mV', bestCache(i).blendFit.rmseMv));
    title(sprintf('Window 1 CH%d | corrected g %.4f mm', ...
        bestCache(i).sensorId, bestCache(i).corrFit.gMm));
    xlabel('x (mm)');
    ylabel('mV');
    legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end
