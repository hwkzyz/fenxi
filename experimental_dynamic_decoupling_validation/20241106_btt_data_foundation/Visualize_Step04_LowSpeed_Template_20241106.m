%% Visualize_Step04_LowSpeed_Template_20241106
% Read-only visualization for Step04 low-speed OPRCenterStd templates.
% This script loads the saved Template and optional source cache. It does not
% rebuild low-speed point clouds or recalibrate templates.

clear; clc; close all;

cfg = BTTDataConfig_20241106();
cfg = apply_visualization_defaults_local(cfg);
show_plots = parse_bool_env_local('STEP04_VIS_SHOW_PLOTS', true);
save_figures = parse_bool_env_local('STEP04_VIS_SAVE_FIGURES', true);

sensor_tag = sprintf('S%s', sprintf('%d', cfg.sensor_ids));
template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', sensor_tag, cfg.dataset));
source_file = fullfile(cfg.step04_output_dir, ...
    sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));

fig_dir = fullfile(cfg.figure_root, 'step04_low_speed_template_visualization');
if save_figures && exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

if ~isfile(template_file)
    error('Missing Step04 template. Run Step04 first:\n  %s', template_file);
end
loaded = load(template_file, 'Template');
Template = loaded.Template;
validate_visualization_template_local(Template);

target_blades = get_template_target_blades_local(Template, cfg);

point_cloud = [];
if isfile(source_file)
    loaded_source = load(source_file, 'point_cloud');
    if isfield(loaded_source, 'point_cloud')
        point_cloud = loaded_source.point_cloud;
    end
else
    warning('Step04Vis:MissingSourceCache', ...
        'Source cache is missing; source point-cloud plots will be skipped.');
end

fprintf('\n=== Step04 read-only visualization: %s ===\n', cfg.dataset);
fprintf('Template: %s\n', template_file);
fprintf('Source cache: %s\n', source_file);
fprintf('Target blades: %s\n', mat2str(target_blades));
fprintf('Figure dir: %s\n', fig_dir);

plot_clean_templates_local(Template, target_blades, fig_dir, show_plots, save_figures);
plot_pointcloud_template_fit_local(Template, target_blades, fig_dir, show_plots, save_figures);
plot_gradient_weight_local(Template, target_blades, fig_dir, show_plots, save_figures);
plot_quality_heatmap_local(Template, target_blades, fig_dir, show_plots, save_figures);
plot_summary_table_local(Template, fig_dir, show_plots, save_figures);
if ~isempty(point_cloud)
    plot_source_pointcloud_overview_local(Template, point_cloud, target_blades, fig_dir, show_plots, save_figures);
    plot_representative_target_pulses_local(Template, point_cloud, target_blades, cfg, fig_dir, show_plots, save_figures);
end

fprintf('>>> Saved Step04 visualization figures under:\n  %s\n', fig_dir);


function cfg = apply_visualization_defaults_local(cfg)
if ~isfield(cfg, 'output_root') || isempty(cfg.output_root)
    cfg.output_root = fullfile(cfg.route_dir, 'output');
end
if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end
if ~isfield(cfg, 'figure_root') || isempty(cfg.figure_root)
    cfg.figure_root = fullfile(cfg.output_root, 'figures');
end
end


function validate_visualization_template_local(Template)
if ~isfield(Template, 'SensorBlade') || isempty(Template.SensorBlade)
    error('Template missing SensorBlade entries.');
end
if ~isfield(Template, 'Summary_Table') || ~istable(Template.Summary_Table)
    error('Template missing Summary_Table.');
end
if ~isfield(Template, 'SchemaVersion')
    warning('Step04Vis:MissingSchemaVersion', ...
        'Template has no SchemaVersion; visualization will continue.');
end
end


function target_blades = get_template_target_blades_local(Template, cfg)
if isfield(Template, 'Target_Blades') && ~isempty(Template.Target_Blades)
    target_blades = Template.Target_Blades(:).';
elseif isfield(cfg, 'step04_target_blades') && ~isempty(cfg.step04_target_blades)
    target_blades = cfg.step04_target_blades(:).';
else
    target_blades = 1:Template.Blades_Num;
end
target_blades = unique(target_blades(isfinite(target_blades) & target_blades >= 1), 'stable');
if isempty(target_blades)
    error('Step04 visualization cannot determine target blades from Template or cfg.');
end
end


function plot_clean_templates_local(Template, target_blades, fig_dir, show_plots, save_figures)
sensor_ids = Template.Sensor_IDs;
colors = blade_colors_local(Template, target_blades);
for sid = sensor_ids(:).'
    fig = make_fig_local(sprintf('Step04 CH%d clean templates', sid), [28 17], show_plots);
    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile;
    hold(ax1, 'on'); box(ax1, 'on');
    for blade_id = target_blades
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            continue;
        end
        plot(ax1, entry.x_grid, entry.v_grid, '-', 'LineWidth', 1.45, ...
            'Color', colors(blade_id, :), 'DisplayName', sprintf('B%d', blade_id));
    end
    draw_all_domains_local(ax1, Template, sid, target_blades, colors);
    xlabel(ax1, 'x (mm)');
    ylabel(ax1, 'V');
    title(ax1, sprintf('CH%d final templates', sid));
    legend(ax1, 'Location', 'eastoutside', 'FontSize', 8);
    style_axis_local(ax1);

    ax2 = nexttile;
    hold(ax2, 'on'); box(ax2, 'on');
    for blade_id = target_blades
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            continue;
        end
        plot(ax2, entry.x_grid, entry.dv_dx, '-', 'LineWidth', 1.15, ...
            'Color', colors(blade_id, :), 'DisplayName', sprintf('B%d', blade_id));
    end
    yline(ax2, 0, 'k:', 'HandleVisibility', 'off');
    draw_all_domains_local(ax2, Template, sid, target_blades, colors);
    xlabel(ax2, 'x (mm)');
    ylabel(ax2, 'dV/dx (V/mm)');
    title(ax2, sprintf('CH%d gradients', sid));
    style_axis_local(ax2);

    save_figure_local(fig, fig_dir, sprintf('Step04Vis_Template_Curves_CH%d', sid), save_figures, true);
    close_if_needed_local(fig, show_plots);
end
end


function plot_pointcloud_template_fit_local(Template, target_blades, fig_dir, show_plots, save_figures)
sensor_ids = Template.Sensor_IDs;
colors = blade_colors_local(Template, target_blades);
n_col = min(3, numel(target_blades));
n_row = ceil(numel(target_blades) / n_col);
for sid = sensor_ids(:).'
    fig = make_fig_local(sprintf('Step04 CH%d point cloud fit', sid), [30 18.5], show_plots);
    tiledlayout(fig, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
    for blade_id = target_blades
        ax = nexttile;
        hold(ax, 'on'); box(ax, 'on');
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry)
            title(ax, sprintf('CH%d B%d missing', sid, blade_id));
            style_axis_local(ax);
            continue;
        end
        if has_entry_points_local(entry, 'x_points_preview', 'v_points_preview')
            plot(ax, entry.x_points_preview, entry.v_points_preview, '.', ...
                'Color', [0.82 0.82 0.82], 'MarkerSize', 2, ...
                'DisplayName', 'template input preview');
        end
        if all(isfinite(entry.x_domain))
            shade_domain_local(ax, entry.x_domain, colors(blade_id, :));
        end
        if has_entry_points_local(entry, 'x_used_points_preview', 'v_used_points_preview')
            plot(ax, entry.x_used_points_preview, entry.v_used_points_preview, '.', ...
                'Color', [0.18 0.18 0.18], 'MarkerSize', 2, ...
                'DisplayName', 'final-template used preview');
        end
        if ~isempty(entry.x_grid)
            plot(ax, entry.x_grid, entry.v_grid_raw, '-', 'LineWidth', 0.75, ...
                'Color', [0.35 0.35 0.35], 'DisplayName', 'binned median');
            plot(ax, entry.x_grid, entry.v_grid, '-', 'LineWidth', 1.5, ...
                'Color', colors(blade_id, :), 'DisplayName', 'final template');
        end
        draw_domain_boundaries_local(ax, entry.x_domain, colors(blade_id, :));
        draw_query_safe_boundaries_local(ax, entry.x_query_safe_domain);
        xlabel(ax, 'x (mm)');
        ylabel(ax, 'V');
        title(ax, sprintf('CH%d B%d | %s | domain %.3g mm | final N=%d', ...
            sid, blade_id, char(string(entry.quality_status)), diff(entry.x_domain), ...
            entry.final_point_count), 'Interpreter', 'none');
        style_axis_local(ax);
        if blade_id == target_blades(1)
            legend(ax, 'Location', 'best', 'FontSize', 7);
        end
    end
    save_figure_local(fig, fig_dir, sprintf('Step04Vis_PointCloud_TemplateFit_CH%d', sid), save_figures, 'png');
    close_if_needed_local(fig, show_plots);
end
end


function plot_gradient_weight_local(Template, target_blades, fig_dir, show_plots, save_figures)
sensor_ids = Template.Sensor_IDs;
colors = blade_colors_local(Template, target_blades);
n_col = min(3, numel(target_blades));
n_row = ceil(numel(target_blades) / n_col);
for sid = sensor_ids(:).'
    fig = make_fig_local(sprintf('Step04 CH%d gradient and weight', sid), [30 18.5], show_plots);
    tiledlayout(fig, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
    for blade_id = target_blades
        ax = nexttile;
        hold(ax, 'on'); box(ax, 'on');
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            title(ax, sprintf('CH%d B%d missing', sid, blade_id));
            style_axis_local(ax);
            continue;
        end
        yyaxis(ax, 'left');
        plot(ax, entry.x_grid, entry.dv_dx, '-', 'LineWidth', 1.1, ...
            'Color', colors(blade_id, :), 'DisplayName', 'dV/dx');
        ylabel(ax, 'dV/dx');
        yyaxis(ax, 'right');
        plot(ax, entry.x_grid, entry.weight_grid, '-', 'LineWidth', 0.9, ...
            'Color', [0.20 0.20 0.20], 'DisplayName', 'weight');
        ylabel(ax, 'weight');
        draw_domain_boundaries_local(ax, entry.x_domain, colors(blade_id, :));
        draw_query_safe_boundaries_local(ax, entry.x_query_safe_domain);
        xlabel(ax, 'x (mm)');
        title(ax, sprintf('CH%d B%d gradient/weight', sid, blade_id));
        style_axis_local(ax);
    end
    save_figure_local(fig, fig_dir, sprintf('Step04Vis_Gradient_Weight_CH%d', sid), save_figures, true);
    close_if_needed_local(fig, show_plots);
end
end


function plot_quality_heatmap_local(Template, target_blades, fig_dir, show_plots, save_figures)
S = Template.Summary_Table;
sensor_ids = Template.Sensor_IDs;
[pulse_mat, width_mat, final_mat, amp_mat] = summary_to_matrices_local(S, sensor_ids, target_blades);

fig = make_fig_local('Step04 quality heatmap', [28 12], show_plots);
tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_heatmap_local(pulse_mat, sensor_ids, target_blades, 'Pulse count');
plot_heatmap_local(width_mat, sensor_ids, target_blades, 'x-domain width (mm)');
plot_heatmap_local(final_mat, sensor_ids, target_blades, 'Final points');
plot_heatmap_local(amp_mat, sensor_ids, target_blades, 'Amplitude (V)');
save_figure_local(fig, fig_dir, 'Step04Vis_Quality_Heatmap', save_figures, true);
close_if_needed_local(fig, show_plots);
end


function plot_summary_table_local(Template, fig_dir, show_plots, save_figures)
S = Template.Summary_Table;
show_cols = {'sensor_id', 'blade_id', 'pulse_count', 'point_count', ...
    'wide_point_count', 'final_point_count', 'domain_width_mm', ...
    'query_safe_left_mm', 'query_safe_right_mm', ...
    'amplitude_v', 'source_point_policy', 'quality_status'};
show_cols = show_cols(ismember(show_cols, S.Properties.VariableNames));
fig = make_fig_local('Step04 summary table', [30 11], show_plots);
uitable('Data', sanitize_table_cells_local(table2cell(S(:, show_cols))), ...
    'ColumnName', show_cols, 'Units', 'normalized', 'Position', [0 0 1 1]);
save_figure_local(fig, fig_dir, 'Step04Vis_Summary_Table', save_figures, false);
close_if_needed_local(fig, show_plots);
end


function plot_source_pointcloud_overview_local(Template, point_cloud, target_blades, fig_dir, show_plots, save_figures)
sensor_ids = Template.Sensor_IDs;
n_col = min(4, numel(sensor_ids));
n_row = ceil(numel(sensor_ids) / n_col);
fig = make_fig_local('Step04 source point cloud overview', [30 13], show_plots);
tiledlayout(fig, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
for sid = sensor_ids(:).'
    ax = nexttile;
    hold(ax, 'on'); box(ax, 'on');
    for blade_id = target_blades
        pc = get_point_cloud_entry_local(point_cloud, sid, blade_id);
        if isempty(pc)
            continue;
        end
        x = pc.x_mm(:);
        v = pc.v(:);
        if numel(x) > 3500
            idx = unique(round(linspace(1, numel(x), 3500)));
            x = x(idx);
            v = v(idx);
        end
        plot(ax, x, v, '.', 'MarkerSize', 1.5, 'DisplayName', sprintf('B%d', blade_id));
    end
    xlabel(ax, 'x (mm)');
    ylabel(ax, 'V');
    title(ax, sprintf('CH%d source point clouds, not final calibration points', sid));
    style_axis_local(ax);
    legend(ax, 'Location', 'best', 'FontSize', 7);
end
save_figure_local(fig, fig_dir, 'Step04Vis_Source_PointCloud_Overview', save_figures, 'png');
close_if_needed_local(fig, show_plots);
end


function plot_representative_target_pulses_local(Template, point_cloud, target_blades, cfg, fig_dir, show_plots, save_figures)
target_blade = choose_target_blade_local(cfg, target_blades);
target_sensors = choose_target_sensors_local(cfg, Template.Sensor_IDs);
if isempty(target_sensors) || ~isfinite(target_blade)
    return;
end
fig = make_fig_local(sprintf('Step04 target representative pulses B%d', target_blade), [28 16], show_plots);
tiledlayout(fig, numel(target_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for sid = target_sensors(:).'
    pc = get_point_cloud_entry_local(point_cloud, sid, target_blade);
    entry = get_template_entry_local(Template, sid, target_blade);
    pulse_id = choose_representative_pulse_local(pc);
    [t, x, v] = extract_one_pulse_local(pc, entry, pulse_id);

    ax_time = nexttile;
    hold(ax_time, 'on'); box(ax_time, 'on');
    if ~isempty(t)
        plot(ax_time, 1000 * (t - median(t, 'omitnan')), v, '-', ...
            'Color', [0.00 0.20 0.65], 'LineWidth', 1.1);
    end
    xlabel(ax_time, 'Time from pulse center (ms)');
    ylabel(ax_time, sprintf('CH%d V', sid));
    title(ax_time, sprintf('CH%d B%d representative pulse', sid, target_blade));
    style_axis_local(ax_time);

    ax_x = nexttile;
    hold(ax_x, 'on'); box(ax_x, 'on');
    if ~isempty(entry) && all(isfinite(entry.x_domain))
        shade_domain_local(ax_x, entry.x_domain, [0.95 0.75 0.20]);
    end
    if ~isempty(x)
        plot(ax_x, x, v, '-', 'Color', [0.30 0.30 0.30], 'LineWidth', 0.9, ...
            'DisplayName', 'single source pulse');
    end
    if ~isempty(entry) && ~isempty(entry.x_grid)
        plot(ax_x, entry.x_grid, entry.v_grid, '-', 'Color', [0.80 0.05 0.05], ...
            'LineWidth', 1.5, 'DisplayName', 'final template');
        draw_domain_boundaries_local(ax_x, entry.x_domain, [0.95 0.75 0.20]);
        draw_query_safe_boundaries_local(ax_x, entry.x_query_safe_domain);
    end
    xlabel(ax_x, 'x (mm)');
    ylabel(ax_x, sprintf('CH%d V', sid));
    title(ax_x, sprintf('Pulse in x-domain | %s', pulse_label_local(pulse_id)));
    style_axis_local(ax_x);
    legend(ax_x, 'Location', 'best', 'FontSize', 7);
end
save_figure_local(fig, fig_dir, sprintf('Step04Vis_Target_RepresentativePulses_B%d', target_blade), save_figures, 'png');
close_if_needed_local(fig, show_plots);
end


function target_blade = choose_target_blade_local(cfg, target_blades)
requested = get_cfg_value_local(cfg, 'step05_target_blades', []);
requested = requested(:).';
requested = requested(ismember(requested, target_blades));
if isempty(requested)
    target_blade = target_blades(1);
else
    target_blade = requested(1);
end
end


function target_sensors = choose_target_sensors_local(cfg, sensor_ids)
requested = get_cfg_value_local(cfg, 'step05_analysis_sensors', []);
if isempty(requested)
    target_sensors = sensor_ids(:).';
else
    target_sensors = intersect(requested(:).', sensor_ids(:).', 'stable');
end
end


function value = get_cfg_value_local(cfg, field_name, default_value)
if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
    value = cfg.(field_name);
else
    value = default_value;
end
end


function [pulse_mat, width_mat, final_mat, amp_mat] = summary_to_matrices_local(S, sensor_ids, target_blades)
pulse_mat = nan(numel(sensor_ids), numel(target_blades));
width_mat = nan(numel(sensor_ids), numel(target_blades));
final_mat = nan(numel(sensor_ids), numel(target_blades));
amp_mat = nan(numel(sensor_ids), numel(target_blades));
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    for iBlade = 1:numel(target_blades)
        blade_id = target_blades(iBlade);
        idx = find(S.sensor_id == sid & S.blade_id == blade_id, 1, 'first');
        if isempty(idx)
            continue;
        end
        pulse_mat(iSensor, iBlade) = get_table_value_local(S, idx, 'pulse_count');
        width_mat(iSensor, iBlade) = get_table_value_local(S, idx, 'domain_width_mm');
        final_mat(iSensor, iBlade) = get_table_value_local(S, idx, 'final_point_count');
        amp_mat(iSensor, iBlade) = get_table_value_local(S, idx, 'amplitude_v');
    end
end
end


function value = get_table_value_local(T, idx, name)
value = NaN;
if any(strcmp(T.Properties.VariableNames, name))
    col = T.(name);
    if idx <= numel(col)
        value = col(idx);
    end
end
end


function plot_heatmap_local(M, sensor_ids, target_blades, title_text)
ax = nexttile;
imagesc(ax, 1:numel(target_blades), 1:numel(sensor_ids), M);
set(ax, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), ...
    'XTick', 1:numel(target_blades), 'XTickLabel', compose('B%d', target_blades));
xlabel(ax, 'Blade');
ylabel(ax, 'Sensor');
title(ax, title_text);
colorbar(ax);
style_axis_local(ax);
end


function entry = get_template_entry_local(Template, sid, blade_id)
entry = [];
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && Template.SensorBlade(i).blade_id == blade_id
        entry = Template.SensorBlade(i);
        return;
    end
end
end


function pc = get_point_cloud_entry_local(point_cloud, sid, blade_id)
pc = [];
for i = 1:numel(point_cloud)
    if point_cloud(i).sensor_id == sid && point_cloud(i).blade_id == blade_id
        pc = point_cloud(i);
        return;
    end
end
end


function draw_all_domains_local(ax, Template, sid, target_blades, colors)
yl = ylim(ax);
for blade_id = target_blades
    entry = get_template_entry_local(Template, sid, blade_id);
    if isempty(entry)
        continue;
    end
    draw_domain_boundaries_local(ax, entry.x_domain, colors(blade_id, :), yl);
    draw_query_safe_boundaries_local(ax, entry.x_query_safe_domain, yl);
end
end


function draw_domain_boundaries_local(ax, domain, color, yl)
if nargin < 4
    yl = ylim(ax);
end
if numel(domain) < 2 || any(~isfinite(domain)) || domain(2) <= domain(1)
    return;
end
plot(ax, [domain(1), domain(1)], yl, ':', 'Color', color, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax, [domain(2), domain(2)], yl, ':', 'Color', color, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
end


function draw_query_safe_boundaries_local(ax, domain, yl)
if nargin < 3
    yl = ylim(ax);
end
if numel(domain) < 2 || any(~isfinite(domain)) || domain(2) <= domain(1)
    return;
end
plot(ax, [domain(1), domain(1)], yl, '--', 'Color', [0.05 0.55 0.15], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax, [domain(2), domain(2)], yl, '--', 'Color', [0.05 0.55 0.15], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
end


function shade_domain_local(ax, domain, color)
if numel(domain) < 2 || any(~isfinite(domain)) || domain(2) <= domain(1)
    return;
end
yl = ylim(ax);
patch(ax, [domain(1) domain(2) domain(2) domain(1)], [yl(1) yl(1) yl(2) yl(2)], ...
    color, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end


function colors = blade_colors_local(Template, target_blades)
n = max([Template.Blades_Num, target_blades(:).']);
colors = lines(n);
end


function tf = has_entry_points_local(entry, x_field, v_field)
tf = isfield(entry, x_field) && isfield(entry, v_field) && ...
    ~isempty(entry.(x_field)) && ~isempty(entry.(v_field));
end


function pulse_id = choose_representative_pulse_local(pc)
if isempty(pc) || isempty(pc.pulse_index)
    pulse_id = NaN;
    return;
end
pulse_ids = unique(pc.pulse_index(:));
pulse_ids = pulse_ids(isfinite(pulse_ids));
if isempty(pulse_ids)
    pulse_id = NaN;
    return;
end
pulse_id = pulse_ids(round((numel(pulse_ids) + 1) / 2));
end


function [t, x, v] = extract_one_pulse_local(pc, entry, pulse_id)
t = [];
x = [];
v = [];
if isempty(pc) || ~isfinite(pulse_id)
    return;
end
pulse_mask = pc.pulse_index(:) == pulse_id;
t = pc.t_s(pulse_mask);
xc = 0;
if ~isempty(entry) && isfield(entry, 'xc') && isscalar(entry.xc) && isfinite(entry.xc)
    xc = entry.xc;
elseif isfield(pc, 'xc') && isscalar(pc.xc) && isfinite(pc.xc)
    xc = pc.xc;
end
x = pc.x_mm(pulse_mask) - xc;
v = pc.v(pulse_mask);
[t, order] = sort(t);
x = x(order);
v = v(order);
end


function label = pulse_label_local(pulse_id)
if isfinite(pulse_id)
    label = sprintf('pulse %d', pulse_id);
else
    label = 'pulse unavailable';
end
end


function fig = make_fig_local(name, size_cm, show_plots)
fig = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 size_cm(1) size_cm(2)], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
end


function style_axis_local(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 8, 'TickDir', 'in', 'Box', 'on');
end


function save_figure_local(fig, fig_dir, base_name, save_figures, export_mode)
if ~save_figures
    return;
end
if nargin < 5 || isempty(export_mode)
    export_mode = 'vector';
end
if islogical(export_mode)
    if export_mode
        export_mode = 'vector';
    else
        export_mode = 'image';
    end
end
export_mode = lower(char(string(export_mode)));
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end
png_file = fullfile(fig_dir, [base_name '.png']);
pdf_file = fullfile(fig_dir, [base_name '.pdf']);
emf_file = fullfile(fig_dir, [base_name '.emf']);
exportgraphics(fig, png_file, 'Resolution', 220);
switch export_mode
    case 'vector'
        exportgraphics(fig, pdf_file, 'ContentType', 'vector');
        try
            exportgraphics(fig, emf_file, 'ContentType', 'vector');
        catch
            saveas(fig, emf_file);
        end
    case 'image'
        exportgraphics(fig, pdf_file, 'ContentType', 'image', 'Resolution', 220);
    case 'png'
        % Point-cloud figures remain PNG-only: vector exports are large.
    otherwise
        error('Step04Vis:InvalidExportMode', 'Unknown export mode: %s', export_mode);
end
end


function close_if_needed_local(fig, show_plots)
if ~show_plots
    close(fig);
end
end


function state = visibility_state_local(show_plots)
if show_plots
    state = 'on';
else
    state = 'off';
end
end


function value = parse_bool_env_local(name, default_value)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    value = default_value;
elseif ismember(txt, {'1', 'true', 'yes', 'on'})
    value = true;
elseif ismember(txt, {'0', 'false', 'no', 'off'})
    value = false;
else
    warning('Step04Vis:InvalidBooleanEnv', '%s=%s is not logical; using default.', name, txt);
    value = default_value;
end
end


function cells_out = sanitize_table_cells_local(cells_in)
cells_out = cells_in;
for i = 1:numel(cells_out)
    value = cells_out{i};
    if isstring(value)
        if isscalar(value)
            cells_out{i} = char(value);
        else
            cells_out{i} = char(join(value, ', '));
        end
    elseif ismissing(value)
        cells_out{i} = '';
    elseif isnumeric(value) && isscalar(value)
        if isfinite(value)
            cells_out{i} = sprintf('%.6g', value);
        else
            cells_out{i} = '';
        end
    end
end
end
