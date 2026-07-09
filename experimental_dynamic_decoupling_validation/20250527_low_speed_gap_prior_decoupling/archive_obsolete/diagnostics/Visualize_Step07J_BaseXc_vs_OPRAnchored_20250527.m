%% Visualize Step07J BaseXc_GradientTpl vs OPRAnchored diagnostics
% This script opens diagnostic figures only. It does not save image files.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
rotDir = fullfile(rootDir, '20250527_low_speed_rotating_calibration');
templateDir = fullfile(rotDir, 'output', 'templates');
dynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');

targetBlade = 1;
sensorIds = [1, 3, 6];
sensorTag = ['S', sprintf('%d', sensorIds)];
windowId = parse_positive_integer_env_local('STEP07J_VIS_WINDOW', 1);

Route = struct([]);
Route(1).name = 'BaseXc + GradientTpl';
Route(1).shortName = 'BaseXc';
Route(1).templateFile = fullfile(templateDir, sprintf( ...
    'Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', targetBlade, sensorTag));
Route(1).dynamicFile = fullfile(dynamicDir, sprintf( ...
    'DynamicMap_B%d_%s_SlidingWindows_Diag80L_W3S1_BaseXc_GradientTpl_20250527.mat', targetBlade, sensorTag));
Route(1).resultFile = choose_existing_file_local({
    fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20250527_B%d_%s_Default10After.mat', targetBlade, sensorTag))
    fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20250527_B%d_%s_BaseXcSingleHardSmoke.mat', targetBlade, sensorTag))
    fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20250527_B%d_%s_Diag80L_W3S1_BaseXc_GradientTpl.mat', targetBlade, sensorTag))
    });

Route(2).name = 'OPRAnchored';
Route(2).shortName = 'OPRA';
Route(2).templateFile = fullfile(templateDir, sprintf( ...
    'Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20250527.mat', targetBlade, sensorTag));
Route(2).dynamicFile = fullfile(dynamicDir, sprintf( ...
    'DynamicMap_B%d_%s_SlidingWindows_Diag80L_W3S1_OPRAnchored_20250527.mat', targetBlade, sensorTag));
Route(2).resultFile = choose_existing_file_local({
    fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20250527_B%d_%s_OPRASmokeSingleHard.mat', targetBlade, sensorTag))
    fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_20250527_B%d_%s_OPRA80.mat', targetBlade, sensorTag))
    });

for ir = 1:numel(Route)
    require_file_local(Route(ir).templateFile);
    require_file_local(Route(ir).dynamicFile);
    S = load(Route(ir).templateFile, 'Template');
    D = load(Route(ir).dynamicFile, 'DynamicMap');
    Route(ir).Template = S.Template;
    Route(ir).DynamicMap = D.DynamicMap;
    if ~isempty(Route(ir).resultFile)
        R = load(Route(ir).resultFile, 'Result');
        Route(ir).Result = R.Result;
    else
        Route(ir).Result = [];
    end
end

fprintf('\n=== Step07J BaseXc vs OPRAnchored visual diagnostic ===\n');
fprintf('Window id: %d, blade: B%d, sensors: %s\n', windowId, targetBlade, mat2str(sensorIds));
for ir = 1:numel(Route)
    fprintf('\n[%s]\nTemplate: %s\nDynamic:  %s\nResult:   %s\n', ...
        Route(ir).name, Route(ir).templateFile, Route(ir).dynamicFile, Route(ir).resultFile);
end

plot_coordinate_centers_local(Route, sensorIds);
plot_template_dynamic_overlay_local(Route, sensorIds, windowId);
plot_gap_only_candidate_ranking_local(Route, windowId);
plot_gap_only_collapse_local(Route, sensorIds, windowId);

fprintf('\nInterpretation guide:\n');
fprintf('1) Coordinate center plot shows whether template x=0 and dynamic x=0 use the same origin.\n');
fprintf('2) Template/dynamic overlays show whether the raw selected points sit on the static template before fitting.\n');
fprintf('3) Candidate ranking shows the actual gap_only EO decision: lower RMSE/score wins.\n');
fprintf('4) Collapse plots show the selected gap_only dynamic compensation for each route.\n');
fprintf('\nNo figures were saved. Close the figure windows when finished.\n');

%% Local functions
function plot_coordinate_centers_local(Route, sensorIds)
rows = [];
for ir = 1:numel(Route)
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        Tpl = get_template_sensor_local(Route(ir).Template, sid);
        dynXc = infer_dynamic_xcenter_local(Route(ir).DynamicMap, sid);
        rows = [rows; {Route(ir).shortName, sid, Tpl.xc, dynXc, Tpl.xc - dynXc}]; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'route','sensor_id','template_xc_mm','dynamic_xc_mm','delta_mm'});
disp(T);

fig = figure('Name', 'Step07J coordinate centers: BaseXc vs OPRAnchored', 'Color', 'w');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for ir = 1:numel(Route)
    idx = strcmp(T.route, Route(ir).shortName);
    plot(T.sensor_id(idx), T.template_xc_mm(idx), '-o', 'LineWidth', 1.2, ...
        'DisplayName', [Route(ir).shortName ' template xc']);
    plot(T.sensor_id(idx), T.dynamic_xc_mm(idx), '--s', 'LineWidth', 1.2, ...
        'DisplayName', [Route(ir).shortName ' dynamic xc']);
end
xlabel('Sensor'); ylabel('x center (mm)');
title('Template center vs dynamic-map center');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on; grid on; box on;
barData = nan(numel(sensorIds), numel(Route));
for ir = 1:numel(Route)
    idx = strcmp(T.route, Route(ir).shortName);
    barData(:, ir) = T.delta_mm(idx);
end
bar(sensorIds, barData);
yline(0, '--k');
xlabel('Sensor'); ylabel('template xc - dynamic xc (mm)');
title('Origin mismatch');
legend({Route.shortName}, 'Location', 'best', 'Box', 'off');
drawnow;
end

function plot_template_dynamic_overlay_local(Route, sensorIds, windowId)
fig = figure('Name', 'Step07J raw dynamic points on static templates', 'Color', 'w');
tiledlayout(fig, numel(sensorIds), numel(Route), 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    for ir = 1:numel(Route)
        nexttile; hold on; grid on; box on;
        Tpl = get_template_sensor_local(Route(ir).Template, sid);
        W = Route(ir).DynamicMap.Window(min(windowId, numel(Route(ir).DynamicMap.Window)));
        D = get_dynamic_sensor_local(W, sid);
        vDynMv = (D.V(:) - Tpl.baseline) * 1000;
        plot(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, '-', ...
            'LineWidth', 1.4, 'DisplayName', 'static template');
        plot(D.x_rel(:), vDynMv, '.', 'MarkerSize', 3.5, 'DisplayName', 'dynamic samples');
        xlabel('x_{rel} (mm)', 'Interpreter', 'tex'); ylabel('mV');
        title(sprintf('%s CH%d W%d', Route(ir).shortName, sid, windowId));
        legend('Location', 'best', 'Box', 'off');
    end
end
drawnow;
end

function plot_gap_only_candidate_ranking_local(Route, windowId)
fig = figure('Name', 'Step07J gap_only EO candidate ranking', 'Color', 'w');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ir = 1:numel(Route)
    nexttile; hold on; grid on; box on;
    if isempty(Route(ir).Result)
        title([Route(ir).shortName ' no Result file']);
        continue;
    end
    wr = Route(ir).Result.WindowResult(min(windowId, numel(Route(ir).Result.WindowResult)));
    CT = wr.modelFits.gap_only.CandidateTable;
    [~, order] = sort(CT.EO);
    CT = CT(order, :);
    bar(categorical(string(CT.EO)), CT.weightedRmseMv);
    ylabel('gap\_only weighted RMSE (mV)', 'Interpreter', 'tex');
    xlabel('EO candidate');
    title(sprintf('%s W%d: selected EO%d', Route(ir).shortName, windowId, wr.modelFits.gap_only.EO));
    text(1:height(CT), CT.weightedRmseMv, compose('%.1f', CT.weightedRmseMv), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
drawnow;
end

function plot_gap_only_collapse_local(Route, sensorIds, windowId)
fig = figure('Name', 'Step07J gap_only collapse: selected route fit', 'Color', 'w');
tiledlayout(fig, numel(sensorIds), 2 * numel(Route), 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    for ir = 1:numel(Route)
        if isempty(Route(ir).Result)
            nexttile; axis off; nexttile; axis off;
            continue;
        end
        wr = Route(ir).Result.WindowResult(min(windowId, numel(Route(ir).Result.WindowResult)));
        bundle = decimate_bundle_like_step07j_local(wr.bundle, Route(ir).Result.cfg.maxPointsPerWindow);
        fit = wr.modelFits.gap_only;
        nFit = numel(fit.uMm);
        if numel(bundle.X) ~= nFit
            bundle = truncate_bundle_local(bundle, nFit);
        end
        if isfield(bundle, 'S') && numel(bundle.S) == numel(bundle.X)
            mask = bundle.S == sid;
        else
            mask = bundle.sensorIds(bundle.sensorIndex).' == sid;
        end
        if ~any(mask)
            nexttile; axis off; nexttile; axis off;
            continue;
        end
        xBefore = bundle.X(mask) - fit.dxMm;
        xAfter = xBefore - fit.uMm(mask);
        vObs = bundle.V(mask);
        vPred = fit.VPred(mask);
        [xSort, order] = sort(xAfter);

        nexttile; hold on; grid on; box on;
        plot(xBefore, vObs, '.', 'MarkerSize', 3.5, 'DisplayName', 'samples');
        xlabel('x-dx (mm)'); ylabel('mV');
        title(sprintf('%s CH%d before, EO%d', Route(ir).shortName, sid, fit.EO));

        nexttile; hold on; grid on; box on;
        plot(xAfter, vObs, '.', 'MarkerSize', 3.5, 'DisplayName', 'samples');
        plot(xSort, vPred(order), '-', 'LineWidth', 1.2, 'DisplayName', 'fit');
        xlabel('x-dx-u (mm)'); ylabel('mV');
        title(sprintf('%s CH%d after, RMSE %.1f', Route(ir).shortName, sid, fit.weightedRmseMv));
        legend('Location', 'best', 'Box', 'off');
    end
end
drawnow;
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function bundle = decimate_bundle_like_step07j_local(bundle, maxPoints)
if ~isfinite(maxPoints) || numel(bundle.T) <= maxPoints
    return;
end
idx = unique(round(linspace(1, numel(bundle.T), maxPoints)));
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(bundle, f) && numel(bundle.(f)) >= max(idx)
        bundle.(f) = bundle.(f)(idx);
    end
end
if isfield(bundle, 'pointCount')
    bundle.pointCount = numel(idx);
end
end

function bundle = truncate_bundle_local(bundle, n)
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(bundle, f) && numel(bundle.(f)) >= n
        bundle.(f) = bundle.(f)(1:n);
    end
end
if isfield(bundle, 'pointCount')
    bundle.pointCount = n;
end
end

function D = get_dynamic_sensor_local(W, sid)
idx = find([W.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Dynamic window does not contain CH%d.', sid);
end
D = W.Sensor(idx);
end

function xc = infer_dynamic_xcenter_local(DynamicMap, sid)
xc = NaN;
for iw = 1:numel(DynamicMap.Window)
    D = get_dynamic_sensor_local(DynamicMap.Window(iw), sid);
    if isfield(D, 'x_abs') && isfield(D, 'x_rel') && ~isempty(D.x_abs) && ~isempty(D.x_rel)
        v = D.x_abs(:) - D.x_rel(:);
        xc = median(v(isfinite(v)), 'omitnan');
        return;
    end
end
end

function file = choose_existing_file_local(files)
file = '';
for i = 1:numel(files)
    if isfile(files{i})
        file = files{i};
        return;
    end
end
end

function require_file_local(file)
if ~isfile(file)
    error('Missing required file: %s', file);
end
end

function n = parse_positive_integer_env_local(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    n = defaultValue;
    return;
end
v = str2double(txt);
if ~isfinite(v) || v < 1
    error('%s must be a positive integer.', name);
end
n = floor(v);
end
