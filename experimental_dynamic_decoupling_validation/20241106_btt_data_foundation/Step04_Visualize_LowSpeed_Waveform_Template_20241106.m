clc; clear; close all;

%STEP04_VISUALIZE_LOWSPEED_WAVEFORM_TEMPLATE_20241106
% Visualize the low-speed OPRCenterStd waveform calibration products.
%
% This script is intentionally read-only: it loads the Step04 Template and
% point_cloud products, then draws diagnostic figures for the low-speed
% waveform calibration and trusted-domain selection.

cfg = BTTDataConfig_20241106();

target_blade = get_cfg_value_local(cfg, 'step05_target_blades', 4);
target_blade = target_blade(1);
target_sensors = get_cfg_value_local(cfg, 'step05_analysis_sensors', [2 5 7]);
max_scatter_points = 8000;
save_figures = true;
show_figures = true;

template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_S%s_%s.mat', ...
    sprintf('%d', cfg.sensor_ids), cfg.dataset));
source_file = fullfile(cfg.step04_output_dir, ...
    sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));

if ~isfile(template_file)
    error('Missing Step04 template file. Run Step04 first:\n  %s', template_file);
end
if ~isfile(source_file)
    error('Missing Step04 source-data file. Run Step04 first:\n  %s', source_file);
end

loaded_template = load(template_file, 'Template');
Template = loaded_template.Template;
loaded_source = load(source_file, 'point_cloud');
point_cloud = loaded_source.point_cloud;

fig_dir = fullfile(cfg.figure_root, 'step04_low_speed_waveform_visualization');
if save_figures && exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

fprintf('\n=== Step04 low-speed waveform visualization ===\n');
fprintf('Template file: %s\n', template_file);
fprintf('Target blade: B%d\n', target_blade);
fprintf('Target sensors: %s\n', mat2str(target_sensors));

style = make_figure_style_local();

fig1 = plot_target_pointcloud_template_local( ...
    Template, point_cloud, target_sensors, target_blade, max_scatter_points, style);
save_figure_local(fig1, fig_dir, 'Step04_LowSpeed_TargetB4_S257_PointCloud_Template', save_figures);

fig2 = plot_target_gradient_weight_local(Template, target_sensors, target_blade, style);
save_figure_local(fig2, fig_dir, 'Step04_LowSpeed_TargetB4_S257_Gradient_Weight', save_figures);

fig3 = plot_representative_pulses_local( ...
    Template, point_cloud, target_sensors, target_blade, style);
save_figure_local(fig3, fig_dir, 'Step04_LowSpeed_TargetB4_S257_RepresentativePulses', save_figures);

fig4 = plot_domain_overview_local(Template, style);
save_figure_local(fig4, fig_dir, 'Step04_LowSpeed_AllSensorBlade_DomainOverview', save_figures);

if ~show_figures
    close([fig1 fig2 fig3 fig4]);
end

fprintf('Saved low-speed visualization figures to:\n  %s\n', fig_dir);


function value = get_cfg_value_local(cfg, field_name, default_value)
if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
    value = cfg.(field_name);
else
    value = default_value;
end
end


function style = make_figure_style_local()
style.font = 'Times New Roman';
style.label_fs = 9;
style.tick_fs = 8;
style.line_w = 1.25;
style.template_w = 1.8;
style.colors = lines(8);
style.domain_color = [0.95 0.75 0.20];
style.grid_color = [0.80 0.80 0.80];
end


function fig = plot_target_pointcloud_template_local( ...
    Template, point_cloud, target_sensors, target_blade, max_scatter_points, style)

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 12]);
tiledlayout(numel(target_sensors), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(target_sensors)
    sid = target_sensors(i);
    entry = get_template_entry_local(Template, sid, target_blade);
    pc = get_point_cloud_entry_local(point_cloud, sid, target_blade);

    ax = nexttile;
    hold(ax, 'on');
    if ~isempty(entry) && all(isfinite(entry.x_domain))
        yl_guess = estimate_ylim_local(pc.v, entry.v_grid);
        shade_domain_local(ax, entry.x_domain, yl_guess, style.domain_color);
    end

    if ~isempty(pc) && ~isempty(pc.x_mm)
        xc = 0;
        if ~isempty(entry) && isfield(entry, 'xc') && isscalar(entry.xc) && isfinite(entry.xc)
            xc = entry.xc;
        elseif isfield(pc, 'xc') && isscalar(pc.xc) && isfinite(pc.xc)
            xc = pc.xc;
        end
        x = pc.x_mm(:) - xc;
        v = pc.v(:);
        valid = isfinite(x) & isfinite(v);
        used = false(size(x));
        if isfield(pc, 'is_used_for_final_template') && ...
                numel(pc.is_used_for_final_template) == numel(x)
            used = pc.is_used_for_final_template(:);
        elseif ~isempty(entry) && all(isfinite(entry.x_domain))
            used = x >= entry.x_domain(1) & x <= entry.x_domain(2);
        end

        x_wide = x(valid);
        v_wide = v(valid);
        if numel(x_wide) > max_scatter_points
            idx = unique(round(linspace(1, numel(x_wide), max_scatter_points)));
            x_wide = x_wide(idx);
            v_wide = v_wide(idx);
        end
        scatter(ax, x_wide, v_wide, 4, 'filled', ...
            'MarkerFaceColor', [0.55 0.55 0.55], ...
            'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.12, ...
            'DisplayName', 'wide waveform points');

        x_used = x(valid & used);
        v_used = v(valid & used);
        if numel(x_used) > max_scatter_points
            idx = unique(round(linspace(1, numel(x_used), max_scatter_points)));
            x_used = x_used(idx);
            v_used = v_used(idx);
        end
        scatter(ax, x_used, v_used, 5, 'filled', ...
            'MarkerFaceColor', [0.18 0.18 0.18], ...
            'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22, ...
            'DisplayName', 'used calibration points');
    end

    if ~isempty(entry) && ~isempty(entry.x_grid)
        if isfield(entry, 'wide_x_grid') && ~isempty(entry.wide_x_grid)
            plot(ax, entry.wide_x_grid, entry.wide_v_grid, '--', ...
                'Color', [0.35 0.55 0.95], 'LineWidth', 0.9, ...
                'DisplayName', 'preliminary wide template');
        end
        plot(ax, entry.x_grid, entry.v_grid_raw, '-', ...
            'Color', [0.40 0.60 0.95], 'LineWidth', 0.8, ...
            'DisplayName', 'binned median, trusted');
        plot(ax, entry.x_grid, entry.v_grid, '-', ...
            'Color', [0.00 0.20 0.65], 'LineWidth', style.template_w, ...
            'DisplayName', 'final smoothed template');
        draw_domain_boundaries_local(ax, entry.x_domain, style.domain_color);
    end

    format_axes_local(ax, style);
    ylabel(ax, sprintf('CH%d V (V)', sid), 'FontName', style.font, 'FontSize', style.label_fs);
    title(ax, sprintf('B%d low-speed point cloud and template | CH%d', target_blade, sid), ...
        'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
    if i == numel(target_sensors)
        xlabel(ax, 'x - xc in OPRCenterStd frame (mm)', 'FontName', style.font, 'FontSize', style.label_fs);
    end
    if i == 1
        legend(ax, 'Location', 'best', 'Box', 'off', ...
            'FontName', style.font, 'FontSize', style.tick_fs);
    end
end
end


function fig = plot_target_gradient_weight_local(Template, target_sensors, target_blade, style)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3 3 18 10]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
for i = 1:numel(target_sensors)
    sid = target_sensors(i);
    entry = get_template_entry_local(Template, sid, target_blade);
    if isempty(entry) || isempty(entry.x_grid)
        continue;
    end
    color = style.colors(i, :);
    plot(ax1, entry.x_grid, entry.dv_dx, '-', 'Color', color, 'LineWidth', style.line_w, ...
        'DisplayName', sprintf('CH%d', sid));
    draw_domain_boundaries_local(ax1, entry.x_domain, color);
end
format_axes_local(ax1, style);
ylabel(ax1, 'dV/dx (V/mm)', 'FontName', style.font, 'FontSize', style.label_fs);
title(ax1, sprintf('B%d template gradient and trusted-domain boundaries', target_blade), ...
    'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
legend(ax1, 'Location', 'best', 'Box', 'off', 'FontName', style.font, 'FontSize', style.tick_fs);

ax2 = nexttile;
hold(ax2, 'on');
for i = 1:numel(target_sensors)
    sid = target_sensors(i);
    entry = get_template_entry_local(Template, sid, target_blade);
    if isempty(entry) || isempty(entry.x_grid)
        continue;
    end
    color = style.colors(i, :);
    plot(ax2, entry.x_grid, entry.weight_grid, '-', 'Color', color, 'LineWidth', style.line_w, ...
        'DisplayName', sprintf('CH%d', sid));
    draw_domain_boundaries_local(ax2, entry.x_domain, color);
end
format_axes_local(ax2, style);
xlabel(ax2, 'x - xc in OPRCenterStd frame (mm)', 'FontName', style.font, 'FontSize', style.label_fs);
ylabel(ax2, 'Weight', 'FontName', style.font, 'FontSize', style.label_fs);
ylim(ax2, [0 1.05]);
title(ax2, 'Template weights used by later waveform matching', ...
    'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
end


function fig = plot_representative_pulses_local(Template, point_cloud, target_sensors, target_blade, style)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 18 11]);
tiledlayout(numel(target_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(target_sensors)
    sid = target_sensors(i);
    pc = get_point_cloud_entry_local(point_cloud, sid, target_blade);
    entry = get_template_entry_local(Template, sid, target_blade);
    pulse_id = choose_representative_pulse_local(pc);
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

    ax_time = nexttile;
    hold(ax_time, 'on');
    if ~isempty(t)
        t_ms = 1000 * (t - median(t, 'omitnan'));
        plot(ax_time, t_ms, v, '-', 'Color', [0.00 0.20 0.65], 'LineWidth', style.line_w);
    end
    format_axes_local(ax_time, style);
    ylabel(ax_time, sprintf('CH%d V (V)', sid), 'FontName', style.font, 'FontSize', style.label_fs);
    title(ax_time, sprintf('Representative low-speed pulse | CH%d B%d', sid, target_blade), ...
        'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
    if i == numel(target_sensors)
        xlabel(ax_time, 'Time from pulse center (ms)', 'FontName', style.font, 'FontSize', style.label_fs);
    end

    ax_x = nexttile;
    hold(ax_x, 'on');
    if ~isempty(entry) && all(isfinite(entry.x_domain))
        yl_guess = estimate_ylim_local(v, entry.v_grid);
        shade_domain_local(ax_x, entry.x_domain, yl_guess, style.domain_color);
    end
    if ~isempty(x)
        plot(ax_x, x, v, '-', 'Color', [0.30 0.30 0.30], 'LineWidth', 0.9, ...
            'DisplayName', 'single pulse');
    end
    if ~isempty(entry) && ~isempty(entry.x_grid)
        plot(ax_x, entry.x_grid, entry.v_grid, '-', 'Color', [0.80 0.05 0.05], ...
            'LineWidth', style.template_w, 'DisplayName', 'template');
        draw_domain_boundaries_local(ax_x, entry.x_domain, style.domain_color);
    end
    format_axes_local(ax_x, style);
    title(ax_x, sprintf('Pulse in x-domain | pulse %d', pulse_id), ...
        'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
    if i == numel(target_sensors)
        xlabel(ax_x, 'x - xc in OPRCenterStd frame (mm)', 'FontName', style.font, 'FontSize', style.label_fs);
    end
end
end


function fig = plot_domain_overview_local(Template, style)
S = Template.Summary_Table;
sensor_ids = Template.Sensor_IDs(:).';
blades_num = Template.Blades_Num;
width_mat = nan(numel(sensor_ids), blades_num);
amp_mat = nan(numel(sensor_ids), blades_num);
point_mat = nan(numel(sensor_ids), blades_num);

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    for blade_id = 1:blades_num
        idx = S.sensor_id == sid & S.blade_id == blade_id;
        if any(idx)
            width_mat(i, blade_id) = S.domain_width_mm(idx);
            amp_mat(i, blade_id) = S.amplitude_v(idx);
            point_mat(i, blade_id) = S.point_count(idx);
        end
    end
end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 18 7]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_heatmap_local(nexttile, width_mat, sensor_ids, 'Domain width (mm)', style);
plot_heatmap_local(nexttile, amp_mat, sensor_ids, 'Amplitude (V)', style);
plot_heatmap_local(nexttile, point_mat ./ 1000, sensor_ids, 'Points (x10^3)', style);
end


function plot_heatmap_local(ax, mat, sensor_ids, title_text, style)
imagesc(ax, mat);
axis(ax, 'tight');
set(ax, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), ...
    'XTick', 1:size(mat, 2), 'TickDir', 'in', 'Box', 'on', ...
    'FontName', style.font, 'FontSize', style.tick_fs);
xlabel(ax, 'Blade ID', 'FontName', style.font, 'FontSize', style.label_fs);
ylabel(ax, 'Sensor', 'FontName', style.font, 'FontSize', style.label_fs);
title(ax, title_text, 'FontName', style.font, 'FontSize', style.label_fs, 'FontWeight', 'normal');
cb = colorbar(ax);
cb.FontName = style.font;
cb.FontSize = style.tick_fs;
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


function yl = estimate_ylim_local(varargin)
vals = [];
for i = 1:nargin
    v = varargin{i};
    vals = [vals; v(:)]; %#ok<AGROW>
end
vals = vals(isfinite(vals));
if isempty(vals)
    yl = [0 1];
    return;
end
lo = prctile(vals, 0.5);
hi = prctile(vals, 99.5);
if hi <= lo
    lo = min(vals);
    hi = max(vals);
end
pad = 0.08 * max(hi - lo, eps);
yl = [lo - pad, hi + pad];
end


function shade_domain_local(ax, x_domain, yl, color)
if any(~isfinite(x_domain)) || any(~isfinite(yl)) || diff(x_domain) <= 0 || diff(yl) <= 0
    return;
end
patch(ax, [x_domain(1) x_domain(2) x_domain(2) x_domain(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], color, ...
    'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
ylim(ax, yl);
end


function draw_domain_boundaries_local(ax, x_domain, color)
if any(~isfinite(x_domain)) || diff(x_domain) <= 0
    return;
end
yl = ylim(ax);
plot(ax, [x_domain(1) x_domain(1)], yl, ':', 'Color', color, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax, [x_domain(2) x_domain(2)], yl, ':', 'Color', color, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
ylim(ax, yl);
end


function format_axes_local(ax, style)
set(ax, 'FontName', style.font, 'FontSize', style.tick_fs, ...
    'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.8);
grid(ax, 'off');
end


function save_figure_local(fig, fig_dir, base_name, save_figures)
if ~save_figures
    return;
end
png_file = fullfile(fig_dir, [base_name '.png']);
pdf_file = fullfile(fig_dir, [base_name '.pdf']);
exportgraphics(fig, png_file, 'Resolution', 300);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');
fprintf('Saved figure: %s\n', png_file);
end
