%% Visualize_DynamicSelection_SoftVsHard_20251222
% Compare dynamic data selection between:
%   1) single-pulse + hard query-safe selection
%   2) all-pulse + soft weighted selection

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
template_file = fullfile(route_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
dynamic_file = fullfile(route_dir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_20251222.mat');
soft_result_file = fullfile(route_dir, 'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_GradientXRange030_20251222.mat');
hard_result_file = fullfile(route_dir, 'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_HardSingleUnifiedTest_GradientXRange030_20251222.mat');

if exist(template_file, 'file') ~= 2
    error('Template file not found: %s', template_file);
end
if exist(dynamic_file, 'file') ~= 2
    error('Dynamic map file not found: %s', dynamic_file);
end
if exist(soft_result_file, 'file') ~= 2
    error('Soft main result not found: %s', soft_result_file);
end
if exist(hard_result_file, 'file') ~= 2
    warning('Hard single test result not found. Run the hard/single test first to plot EO comparison.');
end

loaded_template = load(template_file, 'Template');
loaded_dynamic = load(dynamic_file, 'DynamicMap');
Template = loaded_template.Template;
DynamicMap = loaded_dynamic.DynamicMap;

figure_dir = fullfile(route_dir, 'output', 'figures', 'dynamic_selection_compare');
table_dir = fullfile(route_dir, 'output', 'identification');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

query_guard_mm = 0.75;
hard_method = make_method_local('single', 'hard', 0, query_guard_mm);
soft_method = make_method_local('all', 'soft', 0.2, query_guard_mm);

summary_rows = build_selection_summary_local(DynamicMap, Template, hard_method, soft_method);
summary_table = struct2table(summary_rows);
summary_csv = fullfile(table_dir, 'Step03_DynamicSelection_SoftVsHard_20251222.csv');
writetable(summary_table, summary_csv);

plot_window_id = 1;
plot_window_selection_local(DynamicMap.Window(plot_window_id), Template, hard_method, soft_method, figure_dir);
plot_summary_local(summary_table, soft_result_file, hard_result_file, figure_dir);

fprintf('Saved selection summary: %s\n', summary_csv);
fprintf('Saved window comparison figure: %s\n', fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_W01_20251222.png'));
fprintf('Saved trend comparison figure: %s\n', fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_Summary_20251222.png'));

function method = make_method_local(pulse_mode, domain_mode, soft_margin_mm, query_guard_mm)
method = struct();
method.pulse_selection_mode = pulse_mode;
method.domain_selection_mode = domain_mode;
method.domain_soft_margin_mm = soft_margin_mm;
method.query_guard_mm = query_guard_mm;
method.dynamic_effective_mode = 'gradient';
method.dynamic_template_gradient_min_ratio = 0.08;
method.dynamic_time_gradient_min_ratio = 0.15;
method.dynamic_peak_quantile = 85;
method.domain_margin_mm = 0.02;
method.weight_floor = 0.05;
end

function rows = build_selection_summary_local(DynamicMap, Template, hard_method, soft_method)
rows = repmat(struct( ...
    'window_id', NaN, 'sensor_id', NaN, ...
    'hard_points', NaN, 'soft_points', NaN, ...
    'soft_weight_mean', NaN, 'soft_weight_min', NaN, 'soft_weight_max', NaN, ...
    'hard_x_min', NaN, 'hard_x_max', NaN, 'soft_x_min', NaN, 'soft_x_max', NaN), ...
    numel(DynamicMap.Window) * numel(Template.Sensor), 1);
idx = 0;
for iw = 1:numel(DynamicMap.Window)
    W = DynamicMap.Window(iw);
    for is = 1:numel(Template.Sensor)
        Tpl = Template.Sensor(is);
        D = W.Sensor([W.Sensor.sensor_id] == Tpl.sensor_id);
        hard = select_points_local(D, Tpl, hard_method);
        soft = select_points_local(D, Tpl, soft_method);
        idx = idx + 1;
        rows(idx).window_id = W.window_id;
        rows(idx).sensor_id = Tpl.sensor_id;
        rows(idx).hard_points = nnz(hard.mask);
        rows(idx).soft_points = nnz(soft.mask);
        rows(idx).soft_weight_mean = mean(soft.weight(soft.mask), 'omitnan');
        rows(idx).soft_weight_min = min(soft.weight(soft.mask), [], 'omitnan');
        rows(idx).soft_weight_max = max(soft.weight(soft.mask), [], 'omitnan');
        rows(idx).hard_x_min = min(D.x_rel(hard.mask), [], 'omitnan');
        rows(idx).hard_x_max = max(D.x_rel(hard.mask), [], 'omitnan');
        rows(idx).soft_x_min = min(D.x_rel(soft.mask), [], 'omitnan');
        rows(idx).soft_x_max = max(D.x_rel(soft.mask), [], 'omitnan');
    end
end
rows = rows(1:idx);
end

function plot_window_selection_local(W, Template, hard_method, soft_method, figure_dir)
fig = figure('Color', 'w', 'Position', [80 80 1500 900]);
tl = tiledlayout(numel(Template.Sensor), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('20251222 dynamic selection comparison, window %02d, laps %s', ...
    W.window_id, mat2str(W.lap_range)));

for is = 1:numel(Template.Sensor)
    Tpl = Template.Sensor(is);
    D = W.Sensor([W.Sensor.sensor_id] == Tpl.sensor_id);
    hard = select_points_local(D, Tpl, hard_method);
    soft = select_points_local(D, Tpl, soft_method);
    x_query_safe = [min(Tpl.x_grid(:)) + hard_method.query_guard_mm, ...
        max(Tpl.x_grid(:)) - hard_method.query_guard_mm];

    nexttile;
    plot(D.t, D.V, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, ...
        'DisplayName', 'raw waveform'); hold on;
    scatter(D.t(hard.mask), D.V(hard.mask), 10, [0.00 0.32 0.74], 'filled', ...
        'MarkerFaceAlpha', 0.60, 'DisplayName', 'single+hard selected');
    scatter(D.t(soft.mask), D.V(soft.mask), 8, soft.weight(soft.mask), 'filled', ...
        'MarkerFaceAlpha', 0.70, 'DisplayName', 'all+soft selected');
    yline(Tpl.threshold, 'k--', 'LineWidth', 0.8, 'DisplayName', 'pulse threshold');
    colormap(gca, turbo);
    cb = colorbar;
    cb.Label.String = 'soft weight';
    grid on;
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title(sprintf('S%d time-domain points', Tpl.sensor_id));
    if is == 1
        legend('Location', 'best');
    end

    nexttile;
    plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.2, 'DisplayName', 'static template'); hold on;
    scatter(D.x_rel, D.V, 6, [0.82 0.82 0.82], 'filled', ...
        'MarkerFaceAlpha', 0.25, 'DisplayName', 'all dynamic points');
    scatter(D.x_rel(hard.mask), D.V(hard.mask), 13, [0.00 0.32 0.74], 'filled', ...
        'MarkerFaceAlpha', 0.70, 'DisplayName', 'single+hard');
    scatter(D.x_rel(soft.mask), D.V(soft.mask), 11, soft.weight(soft.mask), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'DisplayName', 'all+soft');
    xline(Tpl.x_domain(1), 'k-.', 'LineWidth', 1.0, 'DisplayName', 'x domain');
    xline(Tpl.x_domain(2), 'k-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(x_query_safe(1), 'm--', 'LineWidth', 1.0, 'DisplayName', 'hard query-safe');
    xline(x_query_safe(2), 'm--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    colormap(gca, turbo);
    cb = colorbar;
    cb.Label.String = 'soft weight';
    grid on;
    xlabel('x_{dynamic} (mm)');
    ylabel('Voltage (V)');
    title(sprintf('S%d x-domain selection: hard=%d, soft=%d', ...
        Tpl.sensor_id, nnz(hard.mask), nnz(soft.mask)));
    if is == 1
        legend('Location', 'best');
    end
end

exportgraphics(fig, fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_W01_20251222.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_W01_20251222.pdf'), 'ContentType', 'image', 'Resolution', 300);
end

function plot_summary_local(summary_table, soft_result_file, hard_result_file, figure_dir)
soft_loaded = load(soft_result_file, 'Result');
soft_trend = soft_loaded.Result.Trend;
has_hard = exist(hard_result_file, 'file') == 2;
if has_hard
    hard_loaded = load(hard_result_file, 'Result');
    hard_trend = hard_loaded.Result.Trend;
else
    hard_trend = table();
end

window_ids = unique(summary_table.window_id);
hard_counts = zeros(size(window_ids));
soft_counts = zeros(size(window_ids));
soft_wmean = zeros(size(window_ids));
for i = 1:numel(window_ids)
    rows = summary_table(summary_table.window_id == window_ids(i), :);
    hard_counts(i) = sum(rows.hard_points, 'omitnan');
    soft_counts(i) = sum(rows.soft_points, 'omitnan');
    soft_wmean(i) = mean(rows.soft_weight_mean, 'omitnan');
end

fig = figure('Color', 'w', 'Position', [120 120 1200 820]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, '20251222 dynamic selection: single+hard vs all+soft');

nexttile;
if has_hard
    plot(hard_trend.window_id, hard_trend.EO_id, 'o-', 'Color', [0.00 0.32 0.74], ...
        'LineWidth', 1.2, 'DisplayName', 'single+hard EO'); hold on;
end
plot(soft_trend.window_id, soft_trend.EO_id, 's-', 'Color', [0.85 0.20 0.10], ...
    'LineWidth', 1.2, 'DisplayName', 'all+soft EO');
yline(14, 'k:', 'DisplayName', 'EO14');
grid on;
ylabel('EO');
legend('Location', 'best');

nexttile;
bar(window_ids - 0.18, hard_counts ./ 1000, 0.36, 'FaceColor', [0.00 0.32 0.74], ...
    'DisplayName', 'single+hard selected points'); hold on;
bar(window_ids + 0.18, soft_counts ./ 1000, 0.36, 'FaceColor', [0.85 0.20 0.10], ...
    'DisplayName', 'all+soft selected points');
grid on;
ylabel('Selected points (x10^3)');
legend('Location', 'best');

nexttile;
plot(window_ids, soft_wmean, 'o-', 'Color', [0.25 0.45 0.20], 'LineWidth', 1.2);
grid on;
xlabel('Window');
ylabel('Mean soft weight');
ylim([0 1.05]);

exportgraphics(fig, fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_Summary_20251222.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, 'Step03_DynamicSelection_SoftVsHard_Summary_20251222.pdf'), 'ContentType', 'image', 'Resolution', 300);
end

function sel = select_points_local(D, Tpl, method)
v_raw = D.V(:);
x_raw = D.x_rel(:);
t_raw = D.t(:);
theta_raw = D.theta(:);
f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_raw, 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x_raw, 'pchip', NaN);
threshold = Tpl.threshold;
x_query_domain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];

if strcmpi(method.pulse_selection_mode, 'all')
    mask_pulse = isolate_all_pulses_local(v_raw, threshold);
else
    mask_pulse = isolate_main_pulse_local(v_raw, threshold);
end
mask_effective = build_dynamic_effective_mask_local(x_raw, v_raw, t_raw, Tpl, threshold, method);
mask_domain = x_raw >= Tpl.x_domain(1) + method.domain_margin_mm & ...
              x_raw <= Tpl.x_domain(2) - method.domain_margin_mm;
mask_query_safe = x_raw >= x_query_domain(1) + method.query_guard_mm & ...
                  x_raw <= x_query_domain(2) - method.query_guard_mm;

mask = mask_pulse & mask_effective & mask_domain & ...
    isfinite(f0) & isfinite(fx0) & isfinite(v_raw) & isfinite(theta_raw);
if strcmpi(method.domain_selection_mode, 'hard')
    safe_mask = mask & mask_query_safe;
    if nnz(safe_mask) >= 8
        mask = safe_mask;
    end
end
if nnz(mask) < 8
    mask = mask_effective & mask_domain & isfinite(f0) & isfinite(fx0) & isfinite(v_raw) & isfinite(theta_raw);
    if strcmpi(method.domain_selection_mode, 'hard')
        safe_mask = mask & mask_query_safe;
        if nnz(safe_mask) >= 8
            mask = safe_mask;
        end
    end
end

weight = zeros(size(v_raw));
if any(mask)
    w_edge = build_edge_weight_local(t_raw(mask), v_raw(mask), method.weight_floor);
    if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
        w_template = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x_raw(mask), 'linear', method.weight_floor);
    else
        w_template = ones(nnz(mask), 1);
    end
    w_domain = build_domain_soft_weight_local(x_raw(mask), x_query_domain, method.domain_soft_margin_mm, method.weight_floor);
    w_query = build_query_guard_soft_weight_local( ...
        x_raw(mask), x_query_domain, method.query_guard_mm, method.domain_soft_margin_mm, method.weight_floor);
    w_gradient = build_template_gradient_weight_local(Tpl, x_raw(mask), method.domain_selection_mode, method.weight_floor);
    w_total = max(method.weight_floor, w_edge(:) .* w_template(:) .* w_domain(:) .* w_query(:) .* w_gradient(:));
    if max(w_total) > 0
        w_total = max(method.weight_floor, w_total ./ max(w_total));
    end
    weight(mask) = w_total;
end
sel = struct('mask', mask, 'weight', weight);
end

function mask = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
end

function mask = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    return;
end
segment_gap = idx(starts(2:end)) - idx(ends(1:end-1));
large_gap_threshold = max(10, round(0.02 * numel(v)));
pulse_breaks = find(segment_gap > large_gap_threshold);
group_starts = [1; pulse_breaks + 1];
group_ends = [pulse_breaks; numel(starts)];
for ig = 1:numel(group_starts)
    seg_ids = group_starts(ig):group_ends(ig);
    best_seg = seg_ids(1);
    best_peak = -inf;
    for iseg = seg_ids
        seg = idx(starts(iseg):ends(iseg));
        peak_val = max(v(seg));
        if peak_val > best_peak
            best_peak = peak_val;
            best_seg = iseg;
        end
    end
    mask(idx(starts(best_seg)):idx(ends(best_seg))) = true;
end
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end
g_tpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    g_tpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
g_tpl_max = max(g_tpl(finite), [], 'omitnan');
if ~isfinite(g_tpl_max) || g_tpl_max <= 0
    mask_tpl = finite;
else
    mask_tpl = g_tpl >= method.dynamic_template_gradient_min_ratio * g_tpl_max;
end
g_time = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    g_time(finite) = abs(gradient(v(finite), t(finite)));
end
g_time_max = max(g_time(finite), [], 'omitnan');
if ~isfinite(g_time_max) || g_time_max <= 0
    mask_time = false(size(v));
else
    mask_time = g_time >= method.dynamic_time_gradient_min_ratio * g_time_max;
end
peak_q = min(max(method.dynamic_peak_quantile, 0), 100);
peak_level = prctile(v(finite), peak_q);
mask_peak = v >= max(threshold, peak_level);
mask = finite & (mask_tpl | mask_time | mask_peak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function w_edge = build_edge_weight_local(t, v, floor_w)
if numel(v) < 3 || range(t) <= 0
    w_edge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    w_edge = dv ./ max(dv);
else
    w_edge = ones(size(dv));
end
w_edge = max(floor_w, w_edge);
end

function w_domain = build_domain_soft_weight_local(x, x_domain, margin_mm, floor_w)
if margin_mm <= 0
    w_domain = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
ratio = min(max(dist_to_edge ./ margin_mm, 0), 1);
w_domain = floor_w + (1 - floor_w) .* ratio;
w_domain(~isfinite(w_domain)) = floor_w;
end

function w_query = build_query_guard_soft_weight_local(x, x_domain, query_guard_mm, margin_mm, floor_w)
if query_guard_mm <= 0 || margin_mm <= 0
    w_query = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
soft_start = max(query_guard_mm - margin_mm, 0);
ratio = min(max((dist_to_edge - soft_start) ./ max(margin_mm, eps), 0), 1);
w_query = floor_w + (1 - floor_w) .* ratio;
w_query(~isfinite(w_query)) = floor_w;
end

function w_gradient = build_template_gradient_weight_local(Tpl, x, domain_selection_mode, floor_w)
if ~strcmpi(domain_selection_mode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    w_gradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
g_max = max(g, [], 'omitnan');
if ~isfinite(g_max) || g_max <= 0
    w_gradient = ones(size(x));
    return;
end
w_gradient = floor_w + (1 - floor_w) .* g ./ g_max;
w_gradient(~isfinite(w_gradient)) = floor_w;
end
