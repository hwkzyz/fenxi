%% Compare_LowSpeedTemplates_AcrossRoutes_20251222
% Compare legacy low-speed OPRCenterStd templates against the current
% 20251222 foundation Step04 templates. This script is read-only
% diagnostics: it does not feed legacy templates into the foundation route.

clear; clc; close all;

cfg = BTTDataConfig_20251222();
foundation_standard_angle_source = 'oprcenter';
foundation_dir = fileparts(mfilename('fullpath'));
validation_root = fileparts(foundation_dir);

blade_ids = 1:cfg.blades_num;
sensor_ids = cfg.sensor_ids(:).';
sensor_tag = sprintf('S%s', sprintf('%d', sensor_ids));

legacy_template_dir = fullfile(validation_root, ...
    '20251222_low_speed_rotating_calibration', 'output', 'templates');
foundation_template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', sensor_tag, cfg.dataset));

out_dir = fullfile(cfg.output_root, 'diagnostics', 'low_speed_template_compare');
fig_dir = fullfile(cfg.figure_root, 'diagnostics', 'low_speed_template_compare');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

fprintf('\n=== Low-speed template route comparison: %s ===\n', cfg.dataset);
fprintf('Legacy template dir:   %s\n', legacy_template_dir);
fprintf('Foundation template:   %s\n', foundation_template_file);
fprintf('Output diagnostics dir: %s\n', out_dir);

if ~isfile(foundation_template_file)
    error('Foundation Step04 template not found: %s', foundation_template_file);
end

F = load(foundation_template_file, 'Template', 'metadata');
FoundationTemplate = F.Template;

rows = repmat(empty_compare_row_local(), 0, 1);
CurveCache = repmat(struct('blade_id', NaN, 'sensor_id', NaN, 'curve', []), 0, 1);
for ib = 1:numel(blade_ids)
    blade_id = blade_ids(ib);
    legacy_template_file = fullfile(legacy_template_dir, ...
        sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_%s.mat', ...
        blade_id, sensor_tag, cfg.dataset));
    if ~isfile(legacy_template_file)
        warning('Legacy template missing for B%d: %s', blade_id, legacy_template_file);
        continue;
    end

    L = load(legacy_template_file, 'Template');
    LegacyTemplate = L.Template;

    for sid = sensor_ids
        legacy_entry = get_legacy_entry_local(LegacyTemplate, sid);
        foundation_entry = get_foundation_entry_local(FoundationTemplate, sid, blade_id);
        [row, curve] = compare_one_template_pair_local( ...
            cfg, blade_id, sid, legacy_template_file, foundation_template_file, ...
            legacy_entry, foundation_entry, foundation_standard_angle_source);
        rows(end + 1, 1) = row; %#ok<SAGROW>
        CurveCache(end + 1).blade_id = blade_id; %#ok<SAGROW>
        CurveCache(end).sensor_id = sid;
        CurveCache(end).curve = curve;
    end
end

CompareTable = struct2table(rows);
summary_file = fullfile(out_dir, sprintf('LowSpeedTemplate_RouteCompare_%s.csv', cfg.dataset));
writetable(CompareTable, summary_file);

role_file = fullfile(out_dir, sprintf('LowSpeedTemplate_RouteRoles_%s.csv', cfg.dataset));
RoleTable = build_route_role_table_local(cfg, legacy_template_dir, foundation_template_file);
writetable(RoleTable, role_file);

fprintf('\nSaved tables:\n  %s\n  %s\n', summary_file, role_file);
disp(CompareTable(:, {'blade_id','sensor_id','raw_corr','raw_rmse_v', ...
    'best_shift_mm','best_shift_corr','best_shift_rmse_v','xc_delta_mm', ...
    'legacy_domain_width_mm','foundation_domain_width_mm'}));

plot_summary_local(CompareTable, blade_ids, sensor_ids, fig_dir, cfg.dataset);
for blade_id = blade_ids
    plot_blade_overlay_local(CompareTable, CurveCache, blade_id, sensor_ids, fig_dir, cfg.dataset);
end

fprintf('\nSaved figures under:\n  %s\n', fig_dir);


function row = empty_compare_row_local()
row = struct( ...
    'dataset', "", ...
    'blade_id', NaN, ...
    'sensor_id', NaN, ...
    'legacy_template_role', "legacy_generated_template", ...
    'foundation_template_role', "foundation_rebuilt_template", ...
    'legacy_template_file', "", ...
    'foundation_template_file', "", ...
    'legacy_standard_angle_source', "GradientXRange030_OPRCenterStd", ...
    'foundation_standard_angle_source', "", ...
    'legacy_xc_mm', NaN, ...
    'foundation_xc_mm', NaN, ...
    'xc_delta_mm', NaN, ...
    'legacy_domain_left_mm', NaN, ...
    'legacy_domain_right_mm', NaN, ...
    'foundation_domain_left_mm', NaN, ...
    'foundation_domain_right_mm', NaN, ...
    'legacy_query_safe_left_mm', NaN, ...
    'legacy_query_safe_right_mm', NaN, ...
    'foundation_query_safe_left_mm', NaN, ...
    'foundation_query_safe_right_mm', NaN, ...
    'legacy_domain_width_mm', NaN, ...
    'foundation_domain_width_mm', NaN, ...
    'domain_overlap_left_mm', NaN, ...
    'domain_overlap_right_mm', NaN, ...
    'domain_overlap_width_mm', NaN, ...
    'legacy_grid_left_mm', NaN, ...
    'legacy_grid_right_mm', NaN, ...
    'foundation_grid_left_mm', NaN, ...
    'foundation_grid_right_mm', NaN, ...
    'legacy_point_count', NaN, ...
    'foundation_point_count', NaN, ...
    'foundation_final_point_count', NaN, ...
    'legacy_wide_point_count', NaN, ...
    'foundation_wide_point_count', NaN, ...
    'legacy_lap_count', NaN, ...
    'foundation_stable_window_lap_count', NaN, ...
    'legacy_v_p2p_v', NaN, ...
    'foundation_v_p2p_v', NaN, ...
    'legacy_dv_dx_peak_v_per_mm', NaN, ...
    'foundation_dv_dx_peak_v_per_mm', NaN, ...
    'raw_corr', NaN, ...
    'raw_rmse_v', NaN, ...
    'demean_rmse_v', NaN, ...
    'best_shift_mm', NaN, ...
    'best_shift_corr', NaN, ...
    'best_shift_rmse_v', NaN, ...
    'best_shift_demean_rmse_v', NaN, ...
    'best_shift_overlap_width_mm', NaN, ...
    'source_point_policy', "", ...
    'quality_status', "");
end


function [row, curve] = compare_one_template_pair_local(cfg, blade_id, sid, ...
    legacy_template_file, foundation_template_file, legacy_entry, foundation_entry, ...
    foundation_standard_angle_source)
row = empty_compare_row_local();
row.dataset = string(cfg.dataset);
row.blade_id = blade_id;
row.sensor_id = sid;
row.legacy_template_file = string(legacy_template_file);
row.foundation_template_file = string(foundation_template_file);
row.foundation_standard_angle_source = string(get_nested_string_local(foundation_entry, ...
    {'standard_angle_source', 'xc_reference_source'}, foundation_standard_angle_source));

[x_old, v_old] = get_curve_local(legacy_entry);
[x_new, v_new] = get_curve_local(foundation_entry);
dv_old = get_numeric_vector_field_local(legacy_entry, {'dv_dx'});
dv_new = get_numeric_vector_field_local(foundation_entry, {'dv_dx'});

old_domain = get_domain_local(legacy_entry, 'x_domain', [min_safe_local(x_old), max_safe_local(x_old)]);
new_domain = get_domain_local(foundation_entry, 'x_domain', [min_safe_local(x_new), max_safe_local(x_new)]);
old_query = get_domain_local(legacy_entry, 'x_query_safe_domain', [NaN NaN]);
new_query = get_domain_local(foundation_entry, 'x_query_safe_domain', [NaN NaN]);

row.legacy_xc_mm = get_numeric_scalar_field_local(legacy_entry, {'xc_mm', 'x_center_mm', 'xc'});
row.foundation_xc_mm = get_numeric_scalar_field_local(foundation_entry, {'xc_mm', 'x_center_mm', 'xc'});
row.xc_delta_mm = row.foundation_xc_mm - row.legacy_xc_mm;

row.legacy_domain_left_mm = old_domain(1);
row.legacy_domain_right_mm = old_domain(2);
row.foundation_domain_left_mm = new_domain(1);
row.foundation_domain_right_mm = new_domain(2);
row.legacy_query_safe_left_mm = old_query(1);
row.legacy_query_safe_right_mm = old_query(2);
row.foundation_query_safe_left_mm = new_query(1);
row.foundation_query_safe_right_mm = new_query(2);
row.legacy_domain_width_mm = diff(old_domain);
row.foundation_domain_width_mm = diff(new_domain);

common_domain = intersect_domains_local(old_domain, new_domain);
row.domain_overlap_left_mm = common_domain(1);
row.domain_overlap_right_mm = common_domain(2);
row.domain_overlap_width_mm = diff(common_domain);

row.legacy_grid_left_mm = min_safe_local(x_old);
row.legacy_grid_right_mm = max_safe_local(x_old);
row.foundation_grid_left_mm = min_safe_local(x_new);
row.foundation_grid_right_mm = max_safe_local(x_new);

row.legacy_point_count = get_numeric_scalar_field_local(legacy_entry, {'point_count'});
row.foundation_point_count = get_numeric_scalar_field_local(foundation_entry, {'point_count'});
row.foundation_final_point_count = get_numeric_scalar_field_local(foundation_entry, {'final_point_count'});
row.legacy_wide_point_count = get_numeric_scalar_field_local(legacy_entry, {'wide_point_count'});
row.foundation_wide_point_count = get_numeric_scalar_field_local(foundation_entry, {'wide_point_count'});
row.legacy_lap_count = get_numeric_scalar_field_local(legacy_entry, {'lap_count'});
row.foundation_stable_window_lap_count = get_numeric_scalar_field_local(foundation_entry, {'stable_window_lap_count'});

row.legacy_v_p2p_v = range_safe_local(v_old);
row.foundation_v_p2p_v = range_safe_local(v_new);
row.legacy_dv_dx_peak_v_per_mm = max(abs(dv_old), [], 'omitnan');
row.foundation_dv_dx_peak_v_per_mm = max(abs(dv_new), [], 'omitnan');
row.source_point_policy = string(get_nested_string_local(foundation_entry, {'source_point_policy'}, ""));
row.quality_status = string(get_nested_string_local(foundation_entry, {'quality_status'}, ""));

[raw_metrics, raw_curve] = compare_curves_on_domain_local(x_old, v_old, x_new, v_new, common_domain, 0);
shift_grid = -1.50:0.01:1.50;
[best_metrics, best_curve] = find_best_shift_local(x_old, v_old, old_domain, x_new, v_new, new_domain, shift_grid);

row.raw_corr = raw_metrics.corr;
row.raw_rmse_v = raw_metrics.rmse_v;
row.demean_rmse_v = raw_metrics.demean_rmse_v;
row.best_shift_mm = best_metrics.shift_mm;
row.best_shift_corr = best_metrics.corr;
row.best_shift_rmse_v = best_metrics.rmse_v;
row.best_shift_demean_rmse_v = best_metrics.demean_rmse_v;
row.best_shift_overlap_width_mm = best_metrics.overlap_width_mm;

curve = struct();
curve.x_old = x_old;
curve.v_old = v_old;
curve.x_new = x_new;
curve.v_new = v_new;
curve.old_domain = old_domain;
curve.new_domain = new_domain;
curve.old_query = old_query;
curve.new_query = new_query;
curve.raw_curve = raw_curve;
curve.best_curve = best_curve;
end


function RoleTable = build_route_role_table_local(cfg, legacy_template_dir, foundation_template_file)
rows = [
    struct('dataset', string(cfg.dataset), 'route_name', "legacy_low_speed_rotating_calibration", ...
    'template_role', "legacy_generated_template", ...
    'main_builder', "Step01_Main_Build_OPRCenterStd_Template_20251222.m", ...
    'publisher_or_wrapper', "Step01L_Build_AllBladeSensor_TemplateLibrary_20251222.m / analysis_newflow Step05 wraps accepted template", ...
    'standard_angle_source', "GradientXRange030_OPRCenterStd", ...
    'template_location', string(legacy_template_dir))
    struct('dataset', string(cfg.dataset), 'route_name', "foundation_step04", ...
    'template_role', "foundation_rebuilt_template", ...
    'main_builder', "Step04_Build_LowSpeed_OPRCenterStd_Template_20251222.m", ...
    'publisher_or_wrapper', "", ...
    'standard_angle_source', "oprcenter", ...
    'template_location', string(foundation_template_file))
    ];
RoleTable = struct2table(rows);
end


function entry = get_legacy_entry_local(Template, sid)
entry = struct();
if isfield(Template, 'Sensor')
    S = Template.Sensor;
    for i = 1:numel(S)
        if isfield(S(i), 'sensor_id') && S(i).sensor_id == sid
            entry = S(i);
            return;
        end
    end
    if numel(S) >= sid
        entry = S(sid);
    end
end
end


function entry = get_foundation_entry_local(Template, sid, blade_id)
entry = struct();
if ~isfield(Template, 'SensorBlade')
    return;
end
SB = Template.SensorBlade;
for i = 1:numel(SB)
    if isfield(SB(i), 'sensor_id') && isfield(SB(i), 'blade_id') && ...
            SB(i).sensor_id == sid && SB(i).blade_id == blade_id
        entry = SB(i);
        return;
    end
end
end


function [x, v] = get_curve_local(entry)
x = get_numeric_vector_field_local(entry, {'x_grid', 'X_grid', 'x'});
v = get_numeric_vector_field_local(entry, {'v_grid', 'V_grid', 'v_template', 'v'});
x = x(:);
v = v(:);
n = min(numel(x), numel(v));
x = x(1:n);
v = v(1:n);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
[x, ia] = unique(x, 'stable');
v = v(ia);
end


function vec = get_numeric_vector_field_local(S, names)
vec = [];
if isempty(S) || ~isstruct(S)
    return;
end
for i = 1:numel(names)
    if isfield(S, names{i}) && isnumeric(S.(names{i}))
        vec = S.(names{i})(:);
        return;
    end
end
end


function val = get_numeric_scalar_field_local(S, names)
val = NaN;
if isempty(S) || ~isstruct(S)
    return;
end
for i = 1:numel(names)
    if isfield(S, names{i}) && isnumeric(S.(names{i})) && ~isempty(S.(names{i}))
        tmp = S.(names{i});
        val = tmp(1);
        return;
    end
end
end


function txt = get_nested_string_local(S, names, default_txt)
txt = default_txt;
if isempty(S) || ~isstruct(S)
    return;
end
for i = 1:numel(names)
    if isfield(S, names{i}) && ~isempty(S.(names{i}))
        txt = string(S.(names{i}));
        return;
    end
end
end


function domain = get_domain_local(S, name, fallback)
domain = fallback;
if isempty(S) || ~isstruct(S) || ~isfield(S, name)
    domain = sanitize_domain_local(domain);
    return;
end
d = S.(name);
if isnumeric(d) && numel(d) >= 2
    domain = d(1:2);
end
domain = sanitize_domain_local(domain);
end


function domain = sanitize_domain_local(domain)
domain = double(domain(:).');
if numel(domain) < 2
    domain = [NaN NaN];
else
    domain = domain(1:2);
end
if ~all(isfinite(domain)) || domain(2) <= domain(1)
    domain = [NaN NaN];
end
end


function common = intersect_domains_local(a, b)
a = sanitize_domain_local(a);
b = sanitize_domain_local(b);
common = [max(a(1), b(1)), min(a(2), b(2))];
if ~all(isfinite(common)) || common(2) <= common(1)
    common = [NaN NaN];
end
end


function [metrics, curve] = compare_curves_on_domain_local(x_old, v_old, x_new, v_new, domain, shift_mm)
metrics = empty_metrics_local();
metrics.shift_mm = shift_mm;
curve = struct('xq', [], 'legacy_v', [], 'foundation_v', [], 'residual_v', []);
domain = sanitize_domain_local(domain);
if ~all(isfinite(domain)) || numel(x_old) < 5 || numel(x_new) < 5
    return;
end
nq = max(200, min(1200, round(diff(domain) / 0.005)));
xq = linspace(domain(1), domain(2), nq).';
vo = interp1(x_old, v_old, xq, 'pchip', NaN);
vn = interp1(x_new, v_new, xq - shift_mm, 'pchip', NaN);
valid = isfinite(vo) & isfinite(vn);
if nnz(valid) < 40
    return;
end
xq = xq(valid);
vo = vo(valid);
vn = vn(valid);
res = vo - vn;
metrics.corr = corr_local(vo, vn);
metrics.rmse_v = sqrt(mean(res .^ 2, 'omitnan'));
metrics.demean_rmse_v = sqrt(mean(((vo - mean(vo, 'omitnan')) - ...
    (vn - mean(vn, 'omitnan'))) .^ 2, 'omitnan'));
metrics.overlap_width_mm = max(xq) - min(xq);
curve.xq = xq;
curve.legacy_v = vo;
curve.foundation_v = vn;
curve.residual_v = res;
end


function [best_metrics, best_curve] = find_best_shift_local(x_old, v_old, old_domain, x_new, v_new, new_domain, shift_grid)
best_metrics = empty_metrics_local();
best_curve = struct('xq', [], 'legacy_v', [], 'foundation_v', [], 'residual_v', []);
for shift_mm = shift_grid
    shifted_new_domain = new_domain + shift_mm;
    domain = intersect_domains_local(old_domain, shifted_new_domain);
    [m, c] = compare_curves_on_domain_local(x_old, v_old, x_new, v_new, domain, shift_mm);
    if ~isfinite(m.corr)
        continue;
    end
    if ~isfinite(best_metrics.corr) || ...
            m.corr > best_metrics.corr + 1e-6 || ...
            (abs(m.corr - best_metrics.corr) <= 1e-6 && m.demean_rmse_v < best_metrics.demean_rmse_v)
        best_metrics = m;
        best_curve = c;
    end
end
end


function metrics = empty_metrics_local()
metrics = struct('shift_mm', NaN, 'corr', NaN, 'rmse_v', NaN, ...
    'demean_rmse_v', NaN, 'overlap_width_mm', NaN);
end


function c = corr_local(a, b)
a = a(:);
b = b(:);
valid = isfinite(a) & isfinite(b);
a = a(valid);
b = b(valid);
if numel(a) < 3
    c = NaN;
    return;
end
a = a - mean(a, 'omitnan');
b = b - mean(b, 'omitnan');
den = sqrt(sum(a .^ 2, 'omitnan') * sum(b .^ 2, 'omitnan'));
if den <= 0
    c = NaN;
else
    c = sum(a .* b, 'omitnan') / den;
end
end


function val = min_safe_local(x)
if isempty(x)
    val = NaN;
else
    val = min(x, [], 'omitnan');
end
end


function val = max_safe_local(x)
if isempty(x)
    val = NaN;
else
    val = max(x, [], 'omitnan');
end
end


function val = range_safe_local(x)
if isempty(x) || all(~isfinite(x))
    val = NaN;
else
    val = max(x, [], 'omitnan') - min(x, [], 'omitnan');
end
end


function plot_summary_local(T, blade_ids, sensor_ids, fig_dir, dataset)
fig = figure('Color', 'w', 'Position', [60, 60, 1280, 840], 'Visible', 'off');
tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

plot_heatmap_metric_local(T, blade_ids, sensor_ids, 'raw_corr', 'raw corr');
plot_heatmap_metric_local(T, blade_ids, sensor_ids, 'best_shift_corr', 'best-shift corr');
plot_heatmap_metric_local(T, blade_ids, sensor_ids, 'best_shift_mm', 'best x-shift (mm)');
plot_heatmap_metric_local(T, blade_ids, sensor_ids, 'raw_rmse_v', 'raw RMSE (V)');
plot_heatmap_metric_local(T, blade_ids, sensor_ids, 'best_shift_rmse_v', 'best-shift RMSE (V)');

nexttile;
hold on;
for sid = sensor_ids
    idx = T.sensor_id == sid;
    plot(T.blade_id(idx), T.foundation_domain_width_mm(idx) - T.legacy_domain_width_mm(idx), ...
        '-o', 'LineWidth', 1.2, 'DisplayName', sprintf('S%d', sid));
end
grid on;
xlabel('Blade');
ylabel('foundation - legacy width (mm)');
title('trusted-domain width difference');
legend('Location', 'best');

sgtitle(sprintf('%s low-speed template route comparison', dataset), 'Interpreter', 'none');
png_file = fullfile(fig_dir, sprintf('LowSpeedTemplate_RouteCompare_Summary_%s.png', dataset));
pdf_file = fullfile(fig_dir, sprintf('LowSpeedTemplate_RouteCompare_Summary_%s.pdf', dataset));
exportgraphics(fig, png_file, 'Resolution', 220);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');
close(fig);
end


function plot_heatmap_metric_local(T, blade_ids, sensor_ids, field_name, title_text)
M = NaN(numel(sensor_ids), numel(blade_ids));
for is = 1:numel(sensor_ids)
    for ib = 1:numel(blade_ids)
        idx = T.sensor_id == sensor_ids(is) & T.blade_id == blade_ids(ib);
        if any(idx)
            M(is, ib) = T.(field_name)(find(idx, 1, 'first'));
        end
    end
end
nexttile;
imagesc(blade_ids, sensor_ids, M);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Blade');
ylabel('Sensor');
title(title_text, 'Interpreter', 'none');
for is = 1:numel(sensor_ids)
    for ib = 1:numel(blade_ids)
        if isfinite(M(is, ib))
            text(blade_ids(ib), sensor_ids(is), sprintf('%.3g', M(is, ib)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 8, ...
                'FontWeight', 'bold');
        end
    end
end
end


function plot_blade_overlay_local(T, CurveCache, blade_id, sensor_ids, fig_dir, dataset)
fig = figure('Color', 'w', 'Position', [80, 80, 1320, 880], 'Visible', 'off');
tiledlayout(fig, numel(sensor_ids), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for sid = sensor_ids
    idx = find([CurveCache.blade_id] == blade_id & [CurveCache.sensor_id] == sid, 1, 'first');
    row_idx = find(T.blade_id == blade_id & T.sensor_id == sid, 1, 'first');
    if isempty(idx) || isempty(row_idx)
        nexttile; title(sprintf('B%d S%d missing', blade_id, sid)); nexttile;
        continue;
    end
    C = CurveCache(idx).curve;
    R = T(row_idx, :);

    nexttile;
    plot(C.x_old, C.v_old, 'k-', 'LineWidth', 1.25, 'DisplayName', 'legacy');
    hold on;
    plot(C.x_new, C.v_new, 'r-', 'LineWidth', 1.10, 'DisplayName', 'foundation');
    plot_domain_lines_local(C.old_domain, [0.20 0.20 0.20], 'legacy domain');
    plot_domain_lines_local(C.new_domain, [0.85 0.10 0.10], 'foundation domain');
    grid on;
    xlabel('x (mm)');
    ylabel('V');
    title(sprintf('B%d S%d raw: corr %.3f, RMSE %.3g V', ...
        blade_id, sid, R.raw_corr, R.raw_rmse_v), 'Interpreter', 'none');
    legend('Location', 'best');

    nexttile;
    if ~isempty(C.best_curve.xq)
        plot(C.best_curve.xq, C.best_curve.legacy_v, 'k-', 'LineWidth', 1.25, ...
            'DisplayName', 'legacy');
        hold on;
        plot(C.best_curve.xq, C.best_curve.foundation_v, 'r-', 'LineWidth', 1.10, ...
            'DisplayName', sprintf('foundation shifted %.2f mm', R.best_shift_mm));
        yyaxis right;
        plot(C.best_curve.xq, C.best_curve.residual_v, '-', 'Color', [0.20 0.45 0.85], ...
            'LineWidth', 0.85, 'DisplayName', 'legacy - foundation');
        ylabel('residual V');
        yyaxis left;
    end
    grid on;
    xlabel('x in legacy coordinate (mm)');
    ylabel('V');
    title(sprintf('best shift: corr %.3f, RMSE %.3g V, shift %.2f mm', ...
        R.best_shift_corr, R.best_shift_rmse_v, R.best_shift_mm), 'Interpreter', 'none');
    legend('Location', 'best');
end
sgtitle(sprintf('%s low-speed template overlay: B%d / S%s', ...
    dataset, blade_id, sprintf('%d', sensor_ids)), 'Interpreter', 'none');
png_file = fullfile(fig_dir, sprintf('LowSpeedTemplate_RouteCompare_B%d_S%s_%s.png', ...
    blade_id, sprintf('%d', sensor_ids), dataset));
pdf_file = fullfile(fig_dir, sprintf('LowSpeedTemplate_RouteCompare_B%d_S%s_%s.pdf', ...
    blade_id, sprintf('%d', sensor_ids), dataset));
exportgraphics(fig, png_file, 'Resolution', 220);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');
close(fig);
end


function plot_domain_lines_local(domain, color, label_text)
domain = sanitize_domain_local(domain);
if ~all(isfinite(domain))
    return;
end
xline(domain(1), '--', 'Color', color, 'LineWidth', 1.0, 'DisplayName', label_text);
xline(domain(2), '--', 'Color', color, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
