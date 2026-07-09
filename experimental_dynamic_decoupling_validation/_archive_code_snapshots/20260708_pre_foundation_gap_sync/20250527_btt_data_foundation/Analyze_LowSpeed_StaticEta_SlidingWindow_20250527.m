clc; close all;

%ANALYZE_LOWSPEED_STATICETA_SLIDINGWINDOW_20250527
% Estimate static inter-sensor x-frame offsets from low-speed waveforms.
%
% Low-speed assumption:
%   V ~= T_s,b(x - d_l - eta_s), with no vibration term.
%
% This script uses the Step04 cached low-speed OPRCenterStd point cloud and
% computes eta_s in 3-lap sliding windows:
%   eta_s = x_center(reference sensor) - x_center(sensor s)
%
% The primary center is an aggregate 3-lap super-Gaussian center.  A
% threshold/centroid pulse-median center is also reported as a diagnostic.

cfg = BTTDataConfig_20250527();

if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end
if ~isfield(cfg, 'output_root') || isempty(cfg.output_root)
    cfg.output_root = fullfile(pwd, 'output');
end

S = struct();
S.blade_id = 1;
S.sensor_ids = [1 3 6];
S.reference_sensor_id = 1;
S.window_laps = 3;
S.sliding_step_laps = 1;
S.min_points_per_window = 150;
S.min_pulses_per_window = 2;
S.output_dir = fullfile(cfg.output_root, 'step05_static_eta_low_speed_sliding');
S.source_file = fullfile(cfg.step04_output_dir, ...
    sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));
S.template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_S%s_%s.mat', ...
    sprintf('%d', cfg.sensor_ids), cfg.dataset));

if exist(S.output_dir, 'dir') ~= 7
    mkdir(S.output_dir);
end

if ~isfile(S.source_file)
    error('Missing Step04 source point cloud: %s', S.source_file);
end

loaded = load(S.source_file, 'point_cloud', 'speed_diagnostic');
point_cloud = loaded.point_cloud;

Template = [];
if isfile(S.template_file)
    loaded_template = load(S.template_file, 'Template');
    Template = loaded_template.Template;
end

PulseCenterTable = build_low_speed_pulse_center_table_local(point_cloud, cfg, S);
WindowEtaTable = build_low_speed_eta_window_table_local(point_cloud, PulseCenterTable, cfg, S);
SummaryTable = build_low_speed_eta_summary_table_local(WindowEtaTable, S);

csv_window = fullfile(S.output_dir, ...
    sprintf('LowSpeed_StaticEta_3LapSliding_B%d_S%s_%s.csv', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));
csv_pulse = fullfile(S.output_dir, ...
    sprintf('LowSpeed_PulseCenters_B%d_S%s_%s.csv', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));
csv_summary = fullfile(S.output_dir, ...
    sprintf('LowSpeed_StaticEta_3LapSliding_Summary_B%d_S%s_%s.csv', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));
mat_file = fullfile(S.output_dir, ...
    sprintf('LowSpeed_StaticEta_3LapSliding_B%d_S%s_%s.mat', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));

writetable(WindowEtaTable, csv_window);
writetable(PulseCenterTable, csv_pulse);
writetable(SummaryTable, csv_summary);
save(mat_file, 'WindowEtaTable', 'PulseCenterTable', 'SummaryTable', ...
    'S', 'cfg', 'Template', '-v7.3');

fprintf('\n=== Low-speed static eta sliding-window analysis ===\n');
fprintf('Dataset: %s, blade B%d, sensors S%s, reference CH%d\n', ...
    cfg.dataset, S.blade_id, sprintf('%d', S.sensor_ids), S.reference_sensor_id);
fprintf('Window: %d laps, step: %d lap\n', S.window_laps, S.sliding_step_laps);
fprintf('Saved window table:\n  %s\n', csv_window);
fprintf('Saved summary table:\n  %s\n', csv_summary);
disp(SummaryTable);

print_low_speed_eta_window_preview_local(WindowEtaTable, S);
plot_low_speed_eta_sliding_local(WindowEtaTable, SummaryTable, S, cfg);


function PulseCenterTable = build_low_speed_pulse_center_table_local(point_cloud, cfg, S)
rows = repmat(make_empty_pulse_center_row_local(), 0, 1);

for i = 1:numel(point_cloud)
    pc = point_cloud(i);
    if pc.blade_id ~= S.blade_id || ~ismember(pc.sensor_id, S.sensor_ids)
        continue;
    end

    pulse_ids = unique(pc.pulse_index(:), 'stable');
    threshold = get_sensor_threshold_local(cfg, pc.sensor_id);

    for ip = 1:numel(pulse_ids)
        pid = pulse_ids(ip);
        mask = pc.pulse_index == pid;
        x = pc.x_mm(mask);
        v = pc.v(mask);
        lap_values = pc.lap_index(mask);

        finite = isfinite(x) & isfinite(v);
        x = x(finite);
        v = v(finite);

        row = make_empty_pulse_center_row_local();
        row.sensor_id = pc.sensor_id;
        row.blade_id = pc.blade_id;
        row.pulse_id = pid;
        row.lap_index = median(lap_values(isfinite(lap_values)), 'omitnan');
        row.point_count = numel(x);
        row.threshold_v = threshold;

        if numel(x) < 10
            row.status = "too_few_points";
            rows(end+1, 1) = row; %#ok<SAGROW>
            continue;
        end

        row.baseline_v = estimate_pulse_baseline_local(v);
        vz = max(v(:) - row.baseline_v, 0);
        row.feature_v = max(v(:), [], 'omitnan') - row.baseline_v;

        high = v(:) >= threshold;
        if nnz(high) < 5 && max(vz) > 0
            high = vz >= prctile(vz, 75);
        end

        if nnz(high) >= 3 && sum(vz(high), 'omitnan') > eps
            row.center_threshold_mm = sum(x(high) .* vz(high), 'omitnan') ./ ...
                sum(vz(high), 'omitnan');
        end

        high85 = vz >= prctile(vz, 85) & vz > 0;
        if nnz(high85) >= 3 && sum(vz(high85), 'omitnan') > eps
            row.center_top85_mm = sum(x(high85) .* vz(high85), 'omitnan') ./ ...
                sum(vz(high85), 'omitnan');
        end

        if isfinite(row.center_threshold_mm)
            row.status = "ok";
        else
            row.status = "no_valid_center";
        end

        rows(end+1, 1) = row; %#ok<SAGROW>
    end
end

PulseCenterTable = struct2table(rows);
end


function WindowEtaTable = build_low_speed_eta_window_table_local(point_cloud, PulseCenterTable, cfg, S)
lap_min = inf;
lap_max = -inf;
for i = 1:numel(point_cloud)
    pc = point_cloud(i);
    if pc.blade_id == S.blade_id && ismember(pc.sensor_id, S.sensor_ids)
        lap_min = min(lap_min, min(pc.lap_index, [], 'omitnan'));
        lap_max = max(lap_max, max(pc.lap_index, [], 'omitnan'));
    end
end

if ~isfinite(lap_min) || ~isfinite(lap_max)
    error('No low-speed point cloud found for blade %d.', S.blade_id);
end

lap_min = ceil(lap_min);
lap_max = floor(lap_max);
starts = lap_min:S.sliding_step_laps:(lap_max - S.window_laps + 1);
rows = repmat(make_empty_eta_window_row_local(S), 0, 1);

for iw = 1:numel(starts)
    lap_start = starts(iw);
    lap_end = lap_start + S.window_laps - 1;

    row = make_empty_eta_window_row_local(S);
    row.window_id = iw;
    row.lap_start = lap_start;
    row.lap_end = lap_end;
    row.lap_count = S.window_laps;

    centers_agg = nan(1, numel(S.sensor_ids));
    centers_pulse = nan(1, numel(S.sensor_ids));
    centers_top85 = nan(1, numel(S.sensor_ids));
    point_counts = nan(1, numel(S.sensor_ids));
    pulse_counts = nan(1, numel(S.sensor_ids));
    pulse_center_mad = nan(1, numel(S.sensor_ids));

    for is = 1:numel(S.sensor_ids)
        sid = S.sensor_ids(is);
        pc = get_point_cloud_entry_local(point_cloud, sid, S.blade_id);
        if isempty(pc)
            continue;
        end

        mask = pc.lap_index >= lap_start & pc.lap_index <= lap_end;
        x = pc.x_mm(mask);
        v = pc.v(mask);
        pulse = pc.pulse_index(mask);
        threshold = get_sensor_threshold_local(cfg, sid);

        [centers_agg(is), agg_info] = estimate_aggregate_window_center_local( ...
            x, v, threshold, S);

        point_counts(is) = agg_info.point_count;
        pulse_counts(is) = numel(unique(pulse(isfinite(pulse))));

        pt = PulseCenterTable( ...
            PulseCenterTable.sensor_id == sid & ...
            PulseCenterTable.blade_id == S.blade_id & ...
            PulseCenterTable.lap_index >= lap_start & ...
            PulseCenterTable.lap_index <= lap_end & ...
            string(PulseCenterTable.status) == "ok", :);

        centers_pulse(is) = median(pt.center_threshold_mm, 'omitnan');
        centers_top85(is) = median(pt.center_top85_mm, 'omitnan');
        pulse_center_mad(is) = mad(pt.center_threshold_mm, 1);
    end

    ref_idx = find(S.sensor_ids == S.reference_sensor_id, 1, 'first');
    ref_center = centers_agg(ref_idx);
    ref_pulse = centers_pulse(ref_idx);
    ref_top85 = centers_top85(ref_idx);

    row.reference_center_agg_mm = ref_center;
    row.reference_center_pulse_mm = ref_pulse;

    for is = 1:numel(S.sensor_ids)
        sid = S.sensor_ids(is);
        row.(sprintf('center_CH%d_agg_mm', sid)) = centers_agg(is);
        row.(sprintf('center_CH%d_pulse_mm', sid)) = centers_pulse(is);
        row.(sprintf('center_CH%d_top85_mm', sid)) = centers_top85(is);
        row.(sprintf('point_count_CH%d', sid)) = point_counts(is);
        row.(sprintf('pulse_count_CH%d', sid)) = pulse_counts(is);
        row.(sprintf('pulse_center_mad_CH%d_mm', sid)) = pulse_center_mad(is);
        row.(sprintf('eta_CH%d_agg_mm', sid)) = ref_center - centers_agg(is);
        row.(sprintf('eta_CH%d_pulse_mm', sid)) = ref_pulse - centers_pulse(is);
        row.(sprintf('eta_CH%d_top85_mm', sid)) = ref_top85 - centers_top85(is);
    end

    ok_points = all(point_counts >= S.min_points_per_window);
    ok_pulses = all(pulse_counts >= S.min_pulses_per_window);
    ok_centers = all(isfinite(centers_agg));
    if ok_points && ok_pulses && ok_centers
        row.status = "ok";
    elseif ~ok_centers
        row.status = "missing_center";
    elseif ~ok_points
        row.status = "too_few_points";
    else
        row.status = "too_few_pulses";
    end

    rows(end+1, 1) = row; %#ok<AGROW>
end

WindowEtaTable = struct2table(rows);
end


function SummaryTable = build_low_speed_eta_summary_table_local(WindowEtaTable, S)
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'reference_sensor_id', S.reference_sensor_id, ...
    'window_count', height(WindowEtaTable), ...
    'ok_window_count', NaN, ...
    'eta_agg_median_mm', NaN, ...
    'eta_agg_mean_mm', NaN, ...
    'eta_agg_std_mm', NaN, ...
    'eta_agg_mad_mm', NaN, ...
    'eta_agg_min_mm', NaN, ...
    'eta_agg_max_mm', NaN, ...
    'eta_pulse_median_mm', NaN, ...
    'eta_top85_median_mm', NaN), 0, 1);

ok = string(WindowEtaTable.status) == "ok";

for sid = S.sensor_ids(:).'
    row = struct();
    row.sensor_id = sid;
    row.reference_sensor_id = S.reference_sensor_id;
    row.window_count = height(WindowEtaTable);
    row.ok_window_count = nnz(ok);

    agg = WindowEtaTable.(sprintf('eta_CH%d_agg_mm', sid));
    pulse = WindowEtaTable.(sprintf('eta_CH%d_pulse_mm', sid));
    top85 = WindowEtaTable.(sprintf('eta_CH%d_top85_mm', sid));

    agg = agg(ok & isfinite(agg));
    pulse = pulse(ok & isfinite(pulse));
    top85 = top85(ok & isfinite(top85));

    row.eta_agg_median_mm = median(agg, 'omitnan');
    row.eta_agg_mean_mm = mean(agg, 'omitnan');
    row.eta_agg_std_mm = std(agg, 'omitnan');
    row.eta_agg_mad_mm = mad(agg, 1);
    row.eta_agg_min_mm = min(agg);
    row.eta_agg_max_mm = max(agg);
    row.eta_pulse_median_mm = median(pulse, 'omitnan');
    row.eta_top85_median_mm = median(top85, 'omitnan');

    rows(end+1, 1) = row; %#ok<AGROW>
end

SummaryTable = struct2table(rows);
end


function row = make_empty_pulse_center_row_local()
row = struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'pulse_id', NaN, ...
    'lap_index', NaN, ...
    'point_count', NaN, ...
    'threshold_v', NaN, ...
    'baseline_v', NaN, ...
    'feature_v', NaN, ...
    'center_threshold_mm', NaN, ...
    'center_top85_mm', NaN, ...
    'status', "");
end


function row = make_empty_eta_window_row_local(S)
row = struct( ...
    'window_id', NaN, ...
    'lap_start', NaN, ...
    'lap_end', NaN, ...
    'lap_count', NaN, ...
    'status', "", ...
    'reference_center_agg_mm', NaN, ...
    'reference_center_pulse_mm', NaN);

for sid = S.sensor_ids(:).'
    row.(sprintf('center_CH%d_agg_mm', sid)) = NaN;
    row.(sprintf('center_CH%d_pulse_mm', sid)) = NaN;
    row.(sprintf('center_CH%d_top85_mm', sid)) = NaN;
    row.(sprintf('point_count_CH%d', sid)) = NaN;
    row.(sprintf('pulse_count_CH%d', sid)) = NaN;
    row.(sprintf('pulse_center_mad_CH%d_mm', sid)) = NaN;
    row.(sprintf('eta_CH%d_agg_mm', sid)) = NaN;
    row.(sprintf('eta_CH%d_pulse_mm', sid)) = NaN;
    row.(sprintf('eta_CH%d_top85_mm', sid)) = NaN;
end
end


function pc = get_point_cloud_entry_local(point_cloud, sensor_id, blade_id)
pc = [];
for i = 1:numel(point_cloud)
    if point_cloud(i).sensor_id == sensor_id && point_cloud(i).blade_id == blade_id
        pc = point_cloud(i);
        return;
    end
end
end


function [xc, info] = estimate_aggregate_window_center_local(x, v, threshold, S)
info = struct('point_count', 0, 'high_point_count', 0, 'mode', '');
xc = NaN;

x = x(:);
v = v(:);
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
info.point_count = numel(x);

if numel(x) < S.min_points_per_window
    info.mode = 'too_few_points';
    return;
end

baseline = estimate_template_selection_baseline_local(v);
vz = max(v - baseline, 0);

high = vz >= prctile(vz, 85) & vz > 0;
info.high_point_count = nnz(high);
if nnz(high) >= 5 && sum(vz(high), 'omitnan') > eps
    xc0 = sum(x(high) .* vz(high), 'omitnan') ./ sum(vz(high), 'omitnan');
else
    [~, imax] = max(vz);
    xc0 = x(imax);
end

strict = v >= threshold;
if nnz(strict) < 30
    strict = vz >= prctile(vz, 60) & vz > 0;
end
if nnz(strict) < 30
    strict = true(size(x));
end

xc = refine_center_by_sg_fit_local(x(strict), v(strict), baseline, xc0);
if ~isfinite(xc)
    xc = xc0;
end
info.mode = 'aggregate_sgfit';
end


function baseline = estimate_pulse_baseline_local(v)
v = v(:);
if isempty(v)
    baseline = NaN;
    return;
end
n = numel(v);
n_edge = max(3, min(floor(0.15 * n), 50));
edge_values = [v(1:n_edge); v((n-n_edge+1):n)];
baseline = median(edge_values, 'omitnan');
if ~isfinite(baseline)
    baseline = prctile(v, 10);
end
end


function baseline = estimate_template_selection_baseline_local(v)
v = v(:);
if isempty(v)
    baseline = NaN;
    return;
end
low_mask = v <= prctile(v, 30);
if nnz(low_mask) >= 10
    baseline = median(v(low_mask), 'omitnan');
else
    baseline = prctile(v, 10);
end
if ~isfinite(baseline)
    baseline = median(v, 'omitnan');
end
end


function xc = refine_center_by_sg_fit_local(x, v, baseline, xc_initial)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
if numel(x) < 50
    xc = xc_initial;
    return;
end
[peak_val, idx_peak] = max(v);
B0 = max(peak_val - baseline, 0.1);
w0 = max(std(x), 0.2);
n0 = 3.0;
xc0 = xc_initial;
if ~isfinite(xc0)
    xc0 = x(idx_peak);
end

sg_model = @(p, xx) p(1) .* exp(-abs((xx - p(4)) ./ max(p(2), 1e-6)).^max(p(3), 1e-6)) + baseline;
obj_fun = @(p) sum((v - sg_model(p, x)).^2);
p0 = [B0, w0, n0, xc0];
lb = [0.05, 0.05, 1.2, min(x) - 0.5];
ub = [max(10, 2 * B0 + 0.5), 6.0, 10.0, max(x) + 0.5];

try
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        p_opt = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], opts);
    else
        bounded_obj = @(p) obj_fun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
        p_opt = fminsearch(bounded_obj, p0, opts);
        p_opt = min(max(p_opt, lb), ub);
    end
    xc = p_opt(4);
catch
    xc = xc_initial;
end
end


function threshold = get_sensor_threshold_local(cfg, sid)
threshold = cfg.sensor_threshold_default;
if isfield(cfg, 'sensor_thresholds') && isa(cfg.sensor_thresholds, 'containers.Map')
    try
        if isKey(cfg.sensor_thresholds, sid)
            threshold = cfg.sensor_thresholds(sid);
        end
    catch
        threshold = cfg.sensor_threshold_default;
    end
elseif isfield(cfg, 'sensor_threshold_ids') && isfield(cfg, 'sensor_threshold_values')
    idx = find(cfg.sensor_threshold_ids == sid, 1, 'first');
    if ~isempty(idx) && idx <= numel(cfg.sensor_threshold_values)
        threshold = cfg.sensor_threshold_values(idx);
    end
elseif isfield(cfg, 'sensor_thresholds') && isnumeric(cfg.sensor_thresholds)
    if sid >= 1 && sid <= numel(cfg.sensor_thresholds)
        threshold = cfg.sensor_thresholds(sid);
    end
end
if ~isfinite(threshold)
    threshold = cfg.sensor_threshold_default;
end
end


function print_low_speed_eta_window_preview_local(T, S)
ok = string(T.status) == "ok";
fprintf('\nWindow eta preview (aggregate sgfit centers):\n');
cols = {'window_id', 'lap_start', 'lap_end', 'status'};
for sid = S.sensor_ids(:).'
    if sid == S.reference_sensor_id
        continue;
    end
    cols{end+1} = sprintf('eta_CH%d_agg_mm', sid); %#ok<AGROW>
end
disp(T(:, cols));

fprintf('OK windows: %d/%d\n', nnz(ok), height(T));
end


function plot_low_speed_eta_sliding_local(T, SummaryTable, S, cfg)
fig = figure('Visible', 'off', ...
    'Name', sprintf('Low-speed static eta B%d S%s', ...
    S.blade_id, sprintf('%d', S.sensor_ids)), ...
    'Color', 'w');
fig.Position(3:4) = [920 520];
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

colors = lines(max(3, numel(S.sensor_ids)));
legend_text = {};

for is = 1:numel(S.sensor_ids)
    sid = S.sensor_ids(is);
    if sid == S.reference_sensor_id
        continue;
    end

    y = T.(sprintf('eta_CH%d_agg_mm', sid));
    plot(ax, T.window_id, y, 'o-', ...
        'LineWidth', 1.2, 'MarkerSize', 4.0, ...
        'Color', colors(is, :));
    legend_text{end+1} = sprintf('CH%d aggregate sgfit', sid); %#ok<AGROW>

    idx = SummaryTable.sensor_id == sid;
    if any(idx)
        med = SummaryTable.eta_agg_median_mm(find(idx, 1, 'first'));
        yline(ax, med, '--', ...
            sprintf('CH%d median %.4f', sid, med), ...
            'Color', colors(is, :), ...
            'LabelHorizontalAlignment', 'left');
    end
end

xlabel(ax, '3-lap sliding window index');
ylabel(ax, sprintf('\\eta_s relative to CH%d (mm)', S.reference_sensor_id));
title(ax, sprintf('Low-speed static eta, B%d, S%s, %s', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset), ...
    'Interpreter', 'none');
legend(ax, legend_text, 'Location', 'best');

fig_file = fullfile(S.output_dir, ...
    sprintf('LowSpeed_StaticEta_3LapSliding_B%d_S%s_%s.png', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));

try
    exportgraphics(fig, fig_file, 'Resolution', 220);
catch
    saveas(fig, fig_file);
end
close(fig);
fprintf('Saved eta trend figure:\n  %s\n', fig_file);
end
