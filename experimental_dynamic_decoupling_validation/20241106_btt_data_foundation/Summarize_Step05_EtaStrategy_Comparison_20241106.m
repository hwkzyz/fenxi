clc; close all;

%SUMMARIZE_STEP05_ETASTRATEGY_COMPARISON_20241106
% Collect Step05 eta-strategy comparison results into one table.

cfg = BTTProjectConfig_20241106();

S = struct();
S.case_name = cfg.dynamic_cases{1};
S.blade_id = cfg.step05_target_blades(1);
S.analysis_sensors = cfg.step05_analysis_sensors;
S.sensor_tag = ['S', sprintf('%d', S.analysis_sensors)];
S.output_dir = fullfile(cfg.output_root, 'step05_eta_strategy_summary');

if exist(S.output_dir, 'dir') ~= 7
    mkdir(S.output_dir);
end

Route = [
    make_route_local('MainFixedJointStaticEta', ...
        fullfile(cfg.output_root, 'step05_single_sync_direct_template', S.case_name, ...
        sprintf('Result_Step05_FoundationMainPulseAdaptiveFixedJointEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'main_fixed_joint_static_eta')
    make_route_local('NoEta', ...
        fullfile(cfg.output_root, 'step05_eta_controlled_compare', 'NoEta', S.case_name, ...
        sprintf('Result_Step05_EtaCompareNoEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'no_eta')
    make_route_local('FreeEta', ...
        fullfile(cfg.output_root, 'step05_eta_controlled_compare', 'WithEta', S.case_name, ...
        sprintf('Result_Step05_EtaCompareWithEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'per_window_free_eta')
    make_route_local('FixedStep04XcEta', ...
        fullfile(cfg.output_root, 'step05_static_eta_compare', 'FixedStep04XcEta', S.case_name, ...
        sprintf('Result_Step05_FixedStep04XcEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'fixed_step04_xc_eta')
    make_route_local('BoundedStep04XcEta_pm002', ...
        fullfile(cfg.output_root, 'step05_static_eta_compare', 'BoundedStep04XcEta', S.case_name, ...
        sprintf('Result_Step05_BoundedStep04XcEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'bounded_step04_xc_eta')
    make_route_local('BoundedStep04XcEta_pm003', ...
        fullfile(cfg.output_root, 'step05_static_eta_compare_bound003', 'BoundedStep04XcEta', S.case_name, ...
        sprintf('Result_Step05_BoundedStep04XcEta_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'bounded_step04_xc_eta')
    make_route_local('FixedJointStaticEtaPreview', ...
        fullfile(cfg.output_root, 'step05_static_eta_compare', 'FixedJointStaticEtaPreview', S.case_name, ...
        sprintf('Result_Step05_FixedJointStaticEtaPreview_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'fixed_joint_static_eta_preview')
    make_route_local('BoundedJointStaticEtaPreview_pm020', ...
        fullfile(cfg.output_root, 'step05_static_eta_compare', 'BoundedJointStaticEtaPreview', S.case_name, ...
        sprintf('Result_Step05_BoundedJointStaticEtaPreview_B%d_%s_%s.mat', ...
        S.blade_id, S.sensor_tag, cfg.dataset)), ...
        'bounded_joint_static_eta_preview')
    ];

summary_rows = repmat(make_empty_summary_row_local(S), 0, 1);
window_rows = repmat(make_empty_window_row_local(S), 0, 1);

for i = 1:numel(Route)
    if ~isfile(Route(i).result_file)
        warning('Missing route result: %s', Route(i).result_file);
        continue;
    end

    loaded = load(Route(i).result_file, 'Result');
    Result = loaded.Result;

    summary_rows(end+1, 1) = build_route_summary_row_local(Route(i), Result, S); %#ok<SAGROW>
    window_rows = [window_rows; build_route_window_rows_local(Route(i), Result, S)]; %#ok<AGROW>
end

SummaryTable = struct2table(summary_rows);
WindowTable = struct2table(window_rows);

summary_csv = fullfile(S.output_dir, ...
    sprintf('EtaStrategySummary_B%d_%s_%s.csv', ...
    S.blade_id, S.sensor_tag, cfg.dataset));
window_csv = fullfile(S.output_dir, ...
    sprintf('EtaStrategyWindow_B%d_%s_%s.csv', ...
    S.blade_id, S.sensor_tag, cfg.dataset));
mat_file = fullfile(S.output_dir, ...
    sprintf('EtaStrategySummary_B%d_%s_%s.mat', ...
    S.blade_id, S.sensor_tag, cfg.dataset));

writetable(SummaryTable, summary_csv);
writetable(WindowTable, window_csv);
save(mat_file, 'SummaryTable', 'WindowTable', 'Route', 'S', 'cfg', '-v7.3');

plot_eta_strategy_summary_local(WindowTable, SummaryTable, S, cfg);

fprintf('\n=== Step05 eta strategy comparison summary ===\n');
fprintf('Saved summary CSV:\n  %s\n', summary_csv);
fprintf('Saved window CSV:\n  %s\n', window_csv);
disp(SummaryTable);


function route = make_route_local(label, result_file, model_class)
route = struct();
route.label = string(label);
route.result_file = result_file;
route.model_class = string(model_class);
end


function row = build_route_summary_row_local(route, Result, S)
T = Result.Trend;
ok = string(T.status) == "ok" & isfinite(T.weighted_voltage_rmse);
eta = collect_eta_matrix_local(Result, S);

row = make_empty_summary_row_local(S);
row.route = route.label;
row.model_class = route.model_class;
row.result_file = string(route.result_file);
row.ok_window_count = nnz(ok);
row.window_count = height(T);
row.eo12_window_count = nnz(ok & T.EO_id == 12);
row.eo4_window_count = nnz(ok & T.EO_id == 4);
row.eo11_window_count = nnz(ok & T.EO_id == 11);
row.other_eo_window_count = nnz(ok & ~ismember(T.EO_id, [4 11 12]));
row.median_rmse_V = median(T.weighted_voltage_rmse(ok), 'omitnan');
row.mean_rmse_V = mean(T.weighted_voltage_rmse(ok), 'omitnan');
row.median_A_mm = median(T.A_id(ok), 'omitnan');
row.median_dx_c_mm = median(T.dx_c_id(ok), 'omitnan');
row.reported_eo = get_struct_numeric_scalar_local(Result.ResonanceSummary, 'reported_eo');
row.reported_freq_hz = get_struct_numeric_scalar_local(Result.ResonanceSummary, 'reported_freq_hz');
row.joint_best_eo = get_struct_numeric_scalar_local(Result.ResonanceSummary, 'joint_best_eo');
row.joint_status = get_struct_string_scalar_local(Result.ResonanceSummary, 'joint_status');
row.joint_gap_ratio = get_struct_numeric_scalar_local(Result.ResonanceSummary, 'joint_eo_gap_ratio');

for i = 1:numel(S.analysis_sensors)
    sid = S.analysis_sensors(i);
    row.(sprintf('eta_CH%d_median_mm', sid)) = median(eta(ok, i), 'omitnan');
    row.(sprintf('eta_CH%d_min_mm', sid)) = min(eta(ok, i));
    row.(sprintf('eta_CH%d_max_mm', sid)) = max(eta(ok, i));
end
end


function rows = build_route_window_rows_local(route, Result, S)
T = Result.Trend;
eta = collect_eta_matrix_local(Result, S);
rows = repmat(make_empty_window_row_local(S), height(T), 1);

for i = 1:height(T)
    rows(i).route = route.label;
    rows(i).model_class = route.model_class;
    rows(i).window_id = T.window_id(i);
    rows(i).lap_start = T.lap_start(i);
    rows(i).lap_end = T.lap_end(i);
    rows(i).status = string(T.status{i});
    rows(i).EO_id = T.EO_id(i);
    rows(i).A_mm = T.A_id(i);
    rows(i).dx_c_mm = T.dx_c_id(i);
    rows(i).weighted_rmse_V = T.weighted_voltage_rmse(i);
    rows(i).plain_rmse_V = T.plain_voltage_rmse(i);
    rows(i).point_count = T.point_count(i);

    for j = 1:numel(S.analysis_sensors)
        sid = S.analysis_sensors(j);
        rows(i).(sprintf('eta_CH%d_mm', sid)) = eta(i, j);
    end
end
end


function eta = collect_eta_matrix_local(Result, S)
eta = nan(height(Result.Trend), numel(S.analysis_sensors));
for i = 1:height(Result.Trend)
    if i > numel(Result.WindowResult) || ...
            isempty(Result.WindowResult(i).Result) || ...
            ~isfield(Result.WindowResult(i).Result, 'sensor_eta_id') || ...
            isempty(Result.WindowResult(i).Result.sensor_eta_id)
        continue;
    end

    e = Result.WindowResult(i).Result.sensor_eta_id(:).';
    eta(i, 1:min(numel(e), size(eta, 2))) = e(1:min(numel(e), size(eta, 2)));
end
end


function row = make_empty_summary_row_local(S)
row = struct( ...
    'route', "", ...
    'model_class', "", ...
    'result_file', "", ...
    'window_count', NaN, ...
    'ok_window_count', NaN, ...
    'eo12_window_count', NaN, ...
    'eo4_window_count', NaN, ...
    'eo11_window_count', NaN, ...
    'other_eo_window_count', NaN, ...
    'reported_eo', NaN, ...
    'reported_freq_hz', NaN, ...
    'joint_best_eo', NaN, ...
    'joint_status', "", ...
    'joint_gap_ratio', NaN, ...
    'median_rmse_V', NaN, ...
    'mean_rmse_V', NaN, ...
    'median_A_mm', NaN, ...
    'median_dx_c_mm', NaN);

for sid = S.analysis_sensors(:).'
    row.(sprintf('eta_CH%d_median_mm', sid)) = NaN;
    row.(sprintf('eta_CH%d_min_mm', sid)) = NaN;
    row.(sprintf('eta_CH%d_max_mm', sid)) = NaN;
end
end


function row = make_empty_window_row_local(S)
row = struct( ...
    'route', "", ...
    'model_class', "", ...
    'window_id', NaN, ...
    'lap_start', NaN, ...
    'lap_end', NaN, ...
    'status', "", ...
    'EO_id', NaN, ...
    'A_mm', NaN, ...
    'dx_c_mm', NaN, ...
    'weighted_rmse_V', NaN, ...
    'plain_rmse_V', NaN, ...
    'point_count', NaN);

for sid = S.analysis_sensors(:).'
    row.(sprintf('eta_CH%d_mm', sid)) = NaN;
end
end


function v = get_struct_numeric_scalar_local(S, name)
v = NaN;
if isstruct(S) && isfield(S, name) && isnumeric(S.(name)) && ...
        isscalar(S.(name)) && isfinite(S.(name))
    v = S.(name);
end
end


function s = get_struct_string_scalar_local(S, name)
s = "";
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    s = string(S.(name));
end
end


function plot_eta_strategy_summary_local(WindowTable, SummaryTable, S, cfg)
if isempty(WindowTable)
    return;
end

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Name', sprintf('Step05 eta strategy comparison B%d %s', ...
    S.blade_id, S.sensor_tag));
fig.Position(3:4) = [1100 760];
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

routes = unique(string(WindowTable.route), 'stable');
colors = lines(numel(routes));

ax1 = nexttile;
hold(ax1, 'on');
grid(ax1, 'on');
for i = 1:numel(routes)
    idx = string(WindowTable.route) == routes(i);
    plot(ax1, WindowTable.window_id(idx), WindowTable.EO_id(idx), ...
        'o-', 'LineWidth', 1.1, 'Color', colors(i, :));
end
ylabel(ax1, 'EO');
title(ax1, 'Identified EO');
legend(ax1, cellstr(routes), 'Location', 'bestoutside', 'Interpreter', 'none');

ax2 = nexttile;
hold(ax2, 'on');
grid(ax2, 'on');
for i = 1:numel(routes)
    idx = string(WindowTable.route) == routes(i);
    plot(ax2, WindowTable.window_id(idx), WindowTable.weighted_rmse_V(idx), ...
        'o-', 'LineWidth', 1.1, 'Color', colors(i, :));
end
ylabel(ax2, 'weighted RMSE (V)');
title(ax2, 'Residual');

ax3 = nexttile;
hold(ax3, 'on');
grid(ax3, 'on');
for i = 1:numel(routes)
    idx = string(WindowTable.route) == routes(i);
    if ismember('eta_CH7_mm', WindowTable.Properties.VariableNames)
        plot(ax3, WindowTable.window_id(idx), WindowTable.eta_CH7_mm(idx), ...
            'o-', 'LineWidth', 1.1, 'Color', colors(i, :));
    end
end
xlabel(ax3, 'window id');
ylabel(ax3, '\eta_{CH7} (mm)');
title(ax3, 'CH7 eta trend');

sgtitle(fig, sprintf('Step05 eta strategy comparison | B%d S%s %s', ...
    S.blade_id, S.sensor_tag, cfg.dataset), 'Interpreter', 'none');

fig_file = fullfile(S.output_dir, ...
    sprintf('EtaStrategyComparison_B%d_%s_%s.png', ...
    S.blade_id, S.sensor_tag, cfg.dataset));
try
    exportgraphics(fig, fig_file, 'Resolution', 220);
catch
    saveas(fig, fig_file);
end
close(fig);
fprintf('Saved comparison figure:\n  %s\n', fig_file);
end
