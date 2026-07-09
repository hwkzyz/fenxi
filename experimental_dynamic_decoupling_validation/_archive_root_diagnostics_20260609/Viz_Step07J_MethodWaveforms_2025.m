%% Viz_Step07J_MethodWaveforms_2025
% Compare the final best-window waveform fits from the direct-template
% route and the nested Step07J static-warp models.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));

cases = [
    struct('dataset', '20250527', ...
    'caseDir', fullfile(rootDir, '20250527_low_speed_gap_prior_decoupling'), ...
    'resultFile', fullfile(rootDir, '20250527_low_speed_gap_prior_decoupling', ...
    'outputs', 'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136.mat'))
    struct('dataset', '20251222', ...
    'caseDir', fullfile(rootDir, '20251222_low_speed_gap_prior_decoupling'), ...
    'resultFile', fullfile(rootDir, '20251222_low_speed_gap_prior_decoupling', ...
    'outputs', 'Step07J_NestedStaticWarp_VPFullWave_20251222_B1_S123.mat'))
    ];

for ic = 1:numel(cases)
    C = cases(ic);
    if ~isfile(C.resultFile)
        warning('Missing Step07J result for %s: %s', C.dataset, C.resultFile);
        continue;
    end
    outDir = fullfile(C.caseDir, 'outputs', 'figures_step07j_method_waveform_compare');
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    S = load(C.resultFile, 'Result', 'Summary');
    Result = S.Result;
    if isfield(Result, 'BestWindowIndex')
        bestIdx = Result.BestWindowIndex;
    else
        bestIdx = 1;
    end
    wr = Result.WindowResult(bestIdx);
    sensorIds = read_sensor_ids_local(wr, Result);

    methods = build_method_series_local(Result, bestIdx);
    summary = build_method_summary_table_local(C.dataset, bestIdx, methods);
    csvFile = fullfile(outDir, sprintf('Step07J_MethodWaveformCompare_%s_BestWindow%d.csv', ...
        C.dataset, bestIdx));
    writetable(summary, csvFile);

    figFile = fullfile(outDir, sprintf('Step07J_MethodWaveformCompare_Compensated_%s_BestWindow%d.png', ...
        C.dataset, bestIdx));
    plot_method_waveform_compare_local(C.dataset, bestIdx, sensorIds, methods, figFile);

    nestedRawFigFile = fullfile(outDir, sprintf('Step07J_NestedRawX_CommonObserved_%s_BestWindow%d.png', ...
        C.dataset, bestIdx));
    plot_nested_rawx_common_observed_local(C.dataset, bestIdx, sensorIds, methods, nestedRawFigFile);

    fprintf('Saved %s method waveform comparison:\n  %s\n  %s\n  %s\n', ...
        C.dataset, figFile, nestedRawFigFile, csvFile);
end

%% Local functions
function sensorIds = read_sensor_ids_local(wr, Result)
sensorIds = [];
if isfield(wr, 'bundle') && isfield(wr.bundle, 'sensorIds')
    sensorIds = wr.bundle.sensorIds(:).';
elseif isfield(Result, 'cfg') && isfield(Result.cfg, 'analysisSensors')
    sensorIds = Result.cfg.analysisSensors(:).';
end
if isempty(sensorIds)
    error('Cannot infer sensor IDs from Step07J result.');
end
end

function methods = build_method_series_local(Result, bestIdx)
methods = struct('name', {}, 'label', {}, 'bundle', {}, 'fit', {}, ...
    'x', {}, 't', {}, 'v', {}, 'sensorIndex', {}, 'sensorIds', {}, ...
    'vPred', {}, 'uMm', {}, 'xComp', {}, 'rmseMv', {}, 'eo', {}, 'ampMm', {}, 'dxMm', {});

if isfield(Result, 'directMainFile') && isfile(Result.directMainFile)
    direct = load_direct_method_local(Result.directMainFile, bestIdx);
    if ~isempty(direct)
        methods(end+1) = direct; %#ok<AGROW>
    end
end

wr = Result.WindowResult(bestIdx);
nestedNames = {'fixed', 'gap_only', 'gap_tilt', 'gap_tilt_shift'};
nestedLabels = {'fixed', 'gap only', 'gap + tilt', 'gap + tilt + shift'};
for i = 1:numel(nestedNames)
    name = nestedNames{i};
    if ~isfield(wr.modelFits, name)
        continue;
    end
    fit = wr.modelFits.(name);
    B = align_bundle_to_prediction_local(wr.bundle, numel(fit.VPred));
    M = struct();
    M.name = name;
    M.label = nestedLabels{i};
    M.bundle = B;
    M.fit = fit;
    M.x = B.X(:);
    M.t = B.TRel(:);
    M.v = B.V(:);
    M.sensorIndex = B.sensorIndex(:);
    M.sensorIds = B.sensorIds(:).';
    M.vPred = fit.VPred(:);
    M.uMm = fit.uMm(:);
    M.rmseMv = fit.weightedRmseMv;
    M.eo = fit.EO;
    M.ampMm = fit.amplitudeMm;
    M.dxMm = fit.dxMm;
    M.xComp = M.x - M.dxMm - M.uMm;
    methods(end+1) = M; %#ok<AGROW>
end
end

function M = load_direct_method_local(directFile, bestIdx)
M = [];
try
    S = load(directFile, 'Result');
catch
    return;
end
if ~isfield(S, 'Result') || ~isfield(S.Result, 'WindowResult') || ...
        numel(S.Result.WindowResult) < bestIdx
    return;
end
wr = S.Result.WindowResult(bestIdx);
if ~isfield(wr, 'bundle') || ~isfield(wr, 'Result')
    return;
end
B0 = wr.bundle;
R = wr.Result;
if ~isfield(R, 'V_pred') || isempty(R.V_pred)
    return;
end
B = convert_direct_bundle_local(B0, numel(R.V_pred));
voltageScale = infer_direct_voltage_scale_local(B.V, R.V_pred);
M = struct();
M.name = 'direct_low_template_main';
M.label = 'direct';
M.bundle = B;
M.fit = R;
M.x = B.X(:);
M.t = B.TRel(:);
M.v = voltageScale .* B.V(:);
M.sensorIndex = B.sensorIndex(:);
M.sensorIds = B.sensorIds(:).';
M.vPred = voltageScale .* R.V_pred(:);
if isfield(R, 'u_est')
    M.uMm = R.u_est(:);
else
    M.uMm = zeros(size(M.x));
end
if isfield(R, 'weighted_voltage_rmse')
    M.rmseMv = voltageScale .* R.weighted_voltage_rmse;
else
    M.rmseMv = sqrt(mean((M.v - M.vPred).^2, 'omitnan'));
end
if isfield(R, 'EO_id'), M.eo = R.EO_id; else, M.eo = NaN; end
if isfield(R, 'A_id'), M.ampMm = R.A_id; else, M.ampMm = NaN; end
if isfield(R, 'dx_c_id'), M.dxMm = R.dx_c_id; else, M.dxMm = NaN; end
if isfinite(M.dxMm)
    M.xComp = M.x - M.dxMm - M.uMm;
else
    M.xComp = M.x - M.uMm;
end
end

function voltageScale = infer_direct_voltage_scale_local(vObs, vPred)
vAll = [vObs(:); vPred(:)];
vRange = max(abs(vAll), [], 'omitnan');
if isfinite(vRange) && vRange < 20
    voltageScale = 1000;
else
    voltageScale = 1;
end
end

function B = convert_direct_bundle_local(B0, nPred)
B = struct();
B.X = B0.X(:);
if isfield(B0, 'T_rel')
    B.TRel = B0.T_rel(:);
elseif isfield(B0, 'T')
    B.TRel = B0.T(:) - min(B0.T(:));
else
    B.TRel = (1:numel(B.X)).';
end
if isfield(B0, 'V')
    B.V = B0.V(:);
else
    B.V = B0.Y_obs(:);
end
if isfield(B0, 'sensor_index')
    B.sensorIndex = B0.sensor_index(:);
else
    B.sensorIndex = B0.sensorIndex(:);
end
if isfield(B0, 'sensor_ids')
    B.sensorIds = B0.sensor_ids(:).';
else
    B.sensorIds = B0.sensorIds(:).';
end
B = decimate_bundle_struct_local(B, nPred);
end

function B = align_bundle_to_prediction_local(B0, nPred)
B = struct();
B.X = B0.X(:);
B.TRel = B0.TRel(:);
B.V = B0.V(:);
B.sensorIndex = B0.sensorIndex(:);
B.sensorIds = B0.sensorIds(:).';
B = decimate_bundle_struct_local(B, nPred);
end

function B = decimate_bundle_struct_local(B, nTarget)
n = numel(B.V);
if nTarget == n
    return;
end
if nTarget <= 0 || nTarget > n
    error('Cannot align bundle with %d points to prediction with %d points.', n, nTarget);
end
idx = unique(round(linspace(1, n, nTarget)));
if numel(idx) ~= nTarget
    idx = round(linspace(1, n, nTarget));
end
B.X = B.X(idx);
B.TRel = B.TRel(idx);
B.V = B.V(idx);
B.sensorIndex = B.sensorIndex(idx);
end

function T = build_method_summary_table_local(dataset, bestIdx, methods)
rows = cell(numel(methods), 1);
for i = 1:numel(methods)
    M = methods(i);
    rows{i} = table(string(dataset), bestIdx, string(M.name), string(M.label), ...
        M.eo, M.ampMm, M.dxMm, M.rmseMv, numel(M.v), ...
        'VariableNames', {'dataset', 'bestWindowIndex', 'method', 'label', ...
        'EO', 'amplitudeMm', 'dxMm', 'weightedRmseMv', 'pointCount'});
end
T = vertcat(rows{:});
end

function plot_method_waveform_compare_local(dataset, bestIdx, sensorIds, methods, figFile)
style = paper_style_local();
colors = method_colors_local();
fig = figure('Name', sprintf('%s Step07J method waveform comparison', dataset), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18.0, max(12.0, 4.6 * numel(sensorIds))]);
tiledlayout(numel(sensorIds), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    nexttile; hold on; box on;
    for im = 1:numel(methods)
        M = methods(im);
        mask = sensor_mask_local(M, sid, is);
        if ~any(mask)
            continue;
        end
        [xSort, order] = sort(M.xComp(mask));
        vObs = M.v(mask);
        vPred = M.vPred(mask);
        c = colors.(M.name);
        plot(xSort, vObs(order), '.', 'Color', lighten_color_local(c, 0.68), ...
            'MarkerSize', 2.6, 'HandleVisibility', 'off');
        plot(xSort, vPred(order), '-', 'Color', colors.(M.name), ...
            'LineWidth', 1.0, 'DisplayName', method_legend_local(M));
    end
    xlabel('x_{OPR}-dx-u(t) (mm)', 'Interpreter', 'tex');
    ylabel('Voltage (mV)');
    title(sprintf('%s W%d CH%d compensated waveform', dataset, bestIdx, sid), 'Interpreter', 'none');
    legend('Location', 'best', 'Box', 'off');
    format_axes_local(gca, style);

    nexttile; hold on; box on;
    for im = 1:numel(methods)
        M = methods(im);
        mask = sensor_mask_local(M, sid, is);
        if ~any(mask)
            continue;
        end
        [xSort, order] = sort(M.xComp(mask));
        res = M.v(mask) - M.vPred(mask);
        plot(xSort, res(order), '-', 'Color', colors.(M.name), ...
            'LineWidth', 0.9, 'DisplayName', method_legend_local(M));
    end
    yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
    xlabel('x_{OPR}-dx-u(t) (mm)', 'Interpreter', 'tex');
    ylabel('Residual (mV)');
    title(sprintf('%s W%d CH%d residual', dataset, bestIdx, sid), 'Interpreter', 'none');
    legend('Location', 'best', 'Box', 'off');
    format_axes_local(gca, style);
end

export_paper_figure_local(fig, figFile);
end

function plot_nested_rawx_common_observed_local(dataset, bestIdx, sensorIds, methods, figFile)
nestedMask = ~strcmp({methods.name}, 'direct_low_template_main');
nestedMethods = methods(nestedMask);
if isempty(nestedMethods)
    return;
end

style = paper_style_local();
colors = method_colors_local();
ref = nestedMethods(1);
fig = figure('Name', sprintf('%s Step07J nested raw-x waveform comparison', dataset), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18.0, max(12.0, 4.6 * numel(sensorIds))]);
tiledlayout(numel(sensorIds), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    refMask = sensor_mask_local(ref, sid, is);

    nexttile; hold on; box on;
    if any(refMask)
        [xObs, obsOrder] = sort(ref.x(refMask));
        vObs = ref.v(refMask);
        plot(xObs, vObs(obsOrder), '.', 'Color', [0.55, 0.55, 0.55], ...
            'MarkerSize', 3.0, 'DisplayName', 'observed');
    end
    for im = 1:numel(nestedMethods)
        M = nestedMethods(im);
        mask = sensor_mask_local(M, sid, is);
        if ~any(mask)
            continue;
        end
        [xSort, order] = sort(M.x(mask));
        vPred = M.vPred(mask);
        plot(xSort, vPred(order), '-', 'Color', colors.(M.name), ...
            'LineWidth', 1.0, 'DisplayName', method_legend_local(M));
    end
    xlabel('x_{OPR} (mm)', 'Interpreter', 'tex');
    ylabel('Voltage (mV)');
    title(sprintf('%s W%d CH%d nested raw-x waveform', dataset, bestIdx, sid), 'Interpreter', 'none');
    legend('Location', 'best', 'Box', 'off');
    format_axes_local(gca, style);

    nexttile; hold on; box on;
    for im = 1:numel(nestedMethods)
        M = nestedMethods(im);
        mask = sensor_mask_local(M, sid, is);
        if ~any(mask)
            continue;
        end
        [xSort, order] = sort(M.x(mask));
        res = M.v(mask) - M.vPred(mask);
        plot(xSort, res(order), '-', 'Color', colors.(M.name), ...
            'LineWidth', 0.9, 'DisplayName', method_legend_local(M));
    end
    yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
    xlabel('x_{OPR} (mm)', 'Interpreter', 'tex');
    ylabel('Residual (mV)');
    title(sprintf('%s W%d CH%d nested raw-x residual', dataset, bestIdx, sid), 'Interpreter', 'none');
    legend('Location', 'best', 'Box', 'off');
    format_axes_local(gca, style);
end

export_paper_figure_local(fig, figFile);
end

function mask = sensor_mask_local(M, sid, sensorOrdinal)
if isfield(M, 'sensorIds') && numel(M.sensorIds) >= sensorOrdinal && M.sensorIds(sensorOrdinal) == sid
    mask = M.sensorIndex == sensorOrdinal;
else
    hit = find(M.sensorIds == sid, 1, 'first');
    if isempty(hit)
        mask = false(size(M.sensorIndex));
    else
        mask = M.sensorIndex == hit;
    end
end
end

function label = method_legend_local(M)
label = sprintf('%s EO%d RMSE %.1f', M.label, round(M.eo), M.rmseMv);
end

function colors = method_colors_local()
colors = struct();
colors.direct_low_template_main = [0.05, 0.05, 0.05];
colors.fixed = [0.00, 0.28, 0.70];
colors.gap_only = [0.82, 0.10, 0.10];
colors.gap_tilt = [0.00, 0.45, 0.28];
colors.gap_tilt_shift = [0.85, 0.37, 0.05];
end

function c2 = lighten_color_local(c, amount)
c2 = c + amount .* (1 - c);
c2 = min(max(c2, 0), 1);
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.tickFontSize = 8.2;
style.gray = [0.38, 0.38, 0.38];
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on');
grid(ax, 'off');
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
xlabel(ax, get(get(ax, 'XLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
ylabel(ax, get(get(ax, 'YLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
lgd = legend(ax);
if ~isempty(lgd) && isvalid(lgd)
    set(lgd, 'FontName', style.fontName, 'FontSize', 7.5);
end
end

function export_paper_figure_local(fig, pngFile)
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
try
    exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
catch
end
try
    print(fig, fullfile(folder, [name, '.emf']), '-dmeta', '-r300');
catch
end
end
