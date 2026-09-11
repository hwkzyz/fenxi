%% Compare Low-Speed And High-Speed Reconstructed Waveforms (20250527)
clc
clear
close all

cfg = build_single_sync_experiment_config_20250527();

%% User Parameters
cfg.analysis_start_time = 1.5;
cfg.target_laps = 8;
cfg.dynamic_window_mode = 'peak_centered_fixed';

compare_blades = 1:cfg.num_blades;
compare_sensors = cfg.base_cfg.sensor_ids;
show_individual_laps = true;

cfg.analysis_sensors = unique(compare_sensors, 'stable');
cfg = resolve_single_sync_experiment_paths_20250527(cfg);

%% Prepare Low-Speed Case
if ~isfile(cfg.sensor_config_file)
    Step1_Build_Sensor_Config_20250527(false);
end
loaded_cfg = load(cfg.sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;

extract_opts = struct('sensor_ids', compare_sensors);
low_case_data = Extract_BTT_Features_20250527(cfg.reference_case_name, cfg.base_cfg, extract_opts);
opr_reference_low = read_sensor_config_opr_reference_local(Sensor_Config);
opr_times_low = low_case_data.opr_times(:) + opr_reference_low.median_center_minus_start_s;
F_omega_deg_low = build_phase_speed_local_local(opr_times_low, cfg.opr_pulses_per_rev);

comparison_rows = repmat(struct( ...
    'BladeID', NaN, ...
    'SensorID', NaN, ...
    'LowPeakMeanV', NaN, ...
    'HighPeakMeanV', NaN, ...
    'PeakDeltaV', NaN, ...
    'LowPeakStdV', NaN, ...
    'HighPeakStdV', NaN), numel(compare_blades) * numel(compare_sensors), 1);

row_idx = 0;
panel_cache = cell(numel(compare_blades), numel(compare_sensors));

for ib = 1:numel(compare_blades)
    blade_id = compare_blades(ib);
    cfg.target_blade = blade_id;
    cfg = resolve_single_sync_experiment_paths_20250527(cfg);
    dynamic_case = prepare_identification_case_single_sync_20250527(cfg, false, false);

    for is = 1:numel(compare_sensors)
        sid = compare_sensors(is);
        low_laps = build_low_speed_laps_local( ...
            cfg, Sensor_Config, low_case_data, opr_times_low, F_omega_deg_low, ...
            opr_reference_low, sid, blade_id);
        dyn_idx = find([dynamic_case.Extracted_Data.sensor_id] == sid, 1, 'first');
        high_laps = dynamic_case.Extracted_Data(dyn_idx).Laps;

        [xq, low_mean, high_mean, low_peak_stats, high_peak_stats] = ...
            summarize_low_high_local(low_laps, high_laps);

        row_idx = row_idx + 1;
        comparison_rows(row_idx).BladeID = blade_id;
        comparison_rows(row_idx).SensorID = sid;
        comparison_rows(row_idx).LowPeakMeanV = low_peak_stats.mean_peak_v;
        comparison_rows(row_idx).HighPeakMeanV = high_peak_stats.mean_peak_v;
        comparison_rows(row_idx).PeakDeltaV = high_peak_stats.mean_peak_v - low_peak_stats.mean_peak_v;
        comparison_rows(row_idx).LowPeakStdV = low_peak_stats.std_peak_v;
        comparison_rows(row_idx).HighPeakStdV = high_peak_stats.std_peak_v;

        panel_cache{ib, is} = struct( ...
            'xq', xq, ...
            'low_mean', low_mean, ...
            'high_mean', high_mean, ...
            'low_peak_stats', low_peak_stats, ...
            'high_peak_stats', high_peak_stats, ...
            'low_laps', low_laps, ...
            'high_laps', high_laps);
    end
end

comparison_rows = comparison_rows(1:row_idx);
comparison_table = struct2table(comparison_rows);
disp(comparison_table);

delta_matrix = nan(numel(compare_blades), numel(compare_sensors));
for ib = 1:numel(compare_blades)
    for is = 1:numel(compare_sensors)
        mask = comparison_table.BladeID == compare_blades(ib) & comparison_table.SensorID == compare_sensors(is);
        if any(mask)
            delta_matrix(ib, is) = comparison_table.PeakDeltaV(find(mask, 1, 'first'));
        end
    end
end

fig_summary = figure('Name', '20250527 Low vs High Peak Delta Summary', ...
    'Color', 'w', 'Position', [80, 60, 1320, 860], 'NumberTitle', 'off');
tiledlayout(fig_summary, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(delta_matrix);
axis tight;
colormap(turbo);
colorbar;
set(gca, 'XTick', 1:numel(compare_sensors), 'XTickLabel', compose('CH%d', compare_sensors), ...
    'YTick', 1:numel(compare_blades), 'YTickLabel', compose('B%d', compare_blades));
xlabel('Sensor');
ylabel('Blade');
title('High-speed peak mean minus low-speed peak mean (V)');
for ib = 1:numel(compare_blades)
    for is = 1:numel(compare_sensors)
        if isfinite(delta_matrix(ib, is))
            text(is, ib, sprintf('%.3f', delta_matrix(ib, is)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end
end

nexttile;
uitable('Parent', fig_summary, ...
    'Data', table2cell(comparison_table), ...
    'ColumnName', comparison_table.Properties.VariableNames, ...
    'Units', 'normalized', ...
    'Position', [0.03 0.05 0.94 0.35]);

for ib = 1:numel(compare_blades)
    blade_id = compare_blades(ib);
    fig_blade = figure('Name', sprintf('20250527 Low vs High Waveform Compare - Blade %d', blade_id), ...
        'Color', 'w', 'Position', [120, 80, 1500, 900], 'NumberTitle', 'off');
    tiledlayout(fig_blade, numel(compare_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for is = 1:numel(compare_sensors)
        sid = compare_sensors(is);
        cache = panel_cache{ib, is};
        row_mask = comparison_table.BladeID == blade_id & comparison_table.SensorID == sid;
        row = comparison_table(find(row_mask, 1, 'first'), :);

        nexttile;
        hold on;
        grid on;
        box on;
        if show_individual_laps
            for k = 1:numel(cache.low_laps)
                plot(cache.low_laps(k).x_points, cache.low_laps(k).v_points, '-', ...
                    'Color', [0.72 0.82 0.96], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            for k = 1:numel(cache.high_laps)
                plot(cache.high_laps(k).x_points, cache.high_laps(k).v_points, '-', ...
                    'Color', [0.96 0.78 0.62], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
        plot(cache.xq, cache.low_mean, '-', 'Color', [0.10 0.45 0.85], 'LineWidth', 2.2, 'DisplayName', 'Low-speed mean');
        plot(cache.xq, cache.high_mean, '--', 'Color', [0.85 0.30 0.10], 'LineWidth', 2.2, 'DisplayName', 'High-speed mean');
        xlabel('Restored position (mm)');
        ylabel('Voltage (V)');
        title(sprintf('B%d CH%d waveform compare | peak delta = %.4f V', blade_id, sid, row.PeakDeltaV));
        legend('Location', 'best');

        nexttile;
        hold on;
        grid on;
        box on;
        plot(1:numel(cache.low_peak_stats.peak_vals), cache.low_peak_stats.peak_vals, '-o', ...
            'Color', [0.10 0.45 0.85], 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Low-speed peaks');
        plot(1:numel(cache.high_peak_stats.peak_vals), cache.high_peak_stats.peak_vals, '-s', ...
            'Color', [0.85 0.30 0.10], 'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'High-speed peaks');
        yline(cache.low_peak_stats.mean_peak_v, ':', 'Color', [0.10 0.45 0.85], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(cache.high_peak_stats.mean_peak_v, ':', 'Color', [0.85 0.30 0.10], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        xlabel('Lap index');
        ylabel('Peak voltage (V)');
        title(sprintf('B%d CH%d peak compare | low %.4f V | high %.4f V', blade_id, sid, ...
            cache.low_peak_stats.mean_peak_v, cache.high_peak_stats.mean_peak_v));
        legend('Location', 'best');
    end
end


function laps = build_low_speed_laps_local(cfg, Sensor_Config, case_data, opr_times, F_omega_deg, opr_reference, sid, blade_id)
selected_rows = select_static_rows_local(cfg, Sensor_Config, case_data, sid, blade_id);
arrival_times = case_data.channels(sid).arrival_times(:);
theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference);

window_start = min(arrival_times(selected_rows) - cfg.pulse_window_sec);
window_end = max(arrival_times(selected_rows) + cfg.pulse_window_sec);
file_ranges = case_data.file_ranges;
overlap_mask = [file_ranges.t_end] >= window_start & [file_ranges.t_start] <= window_end;
file_ranges = file_ranges(overlap_mask);

raw_t = [];
raw_v = [];
for iFile = 1:numel(file_ranges)
    file_id = file_ranges(iFile).file_id;
    offset = file_ranges(iFile).offset;
    filepath = fullfile(cfg.static_data_dir, sprintf('4-%d-%d.mat', sid, file_id));
    if ~isfile(filepath)
        continue;
    end
    [t_local, v_local] = load_raw_case_channel_local(cfg.static_data_dir, sid, file_id, cfg.pinlv);
    t_global = t_local(:) + offset;
    keep_mask = t_global >= window_start & t_global <= window_end;
    if ~any(keep_mask)
        continue;
    end
    raw_t = [raw_t; t_global(keep_mask)]; %#ok<AGROW>
    raw_v = [raw_v; v_local(keep_mask)]; %#ok<AGROW>
end

if isempty(raw_t)
    error('No low-speed raw files were found for CH%d in %s.', sid, cfg.static_data_dir);
end

laps = struct('lap_id', {}, 'row_id', {}, 'x_points', {}, 't_points', {}, 'v_points', {}, 'peak_val', {}, 'peak_time', {});
for k = 1:numel(selected_rows)
    row_id = selected_rows(k);
    t_peak = arrival_times(row_id);
    t_start = t_peak - cfg.pulse_window_sec;
    t_end = t_peak + cfg.pulse_window_sec;
    mask = raw_t >= t_start & raw_t <= t_end;
    if nnz(mask) < 5
        continue;
    end

    t_seg = raw_t(mask);
    v_seg = raw_v(mask);
    idx_prev = find(opr_times < t_peak, 1, 'last');
    if isempty(idx_prev)
        continue;
    end
    t_ref = opr_times(idx_prev);
    theta_points = map_segment_to_relative_angle_local_local(t_ref, t_seg, F_omega_deg);
    theta_diff = mod(theta_points - theta_std + 180, 360) - 180;
    x_points_mm = theta_diff * (pi / 180) * cfg.r_tip_mm;

    [peak_val, idx_pk] = max(v_seg);
    laps(end + 1).lap_id = k; %#ok<AGROW>
    laps(end).row_id = row_id;
    laps(end).x_points = x_points_mm(:);
    laps(end).t_points = t_seg(:);
    laps(end).v_points = v_seg(:);
    laps(end).peak_val = peak_val;
    laps(end).peak_time = t_seg(idx_pk);
end
end


function selected_rows = select_static_rows_local(cfg, Sensor_Config, case_data, sid, blade_id)
start_idx = Sensor_Config.Target_Indices(sid);
arrival_times = case_data.channels(sid).arrival_times(:);
row_ids = start_idx + (blade_id - 1):cfg.num_blades:numel(arrival_times);

if ~isempty(cfg.static_calibration_time_range)
    t_lo = cfg.static_calibration_time_range(1);
    t_hi = cfg.static_calibration_time_range(min(2, end));
    time_mask = arrival_times(row_ids) >= t_lo & arrival_times(row_ids) <= t_hi;
else
    t_lo = cfg.static_calibration_start_time;
    time_mask = arrival_times(row_ids) >= t_lo;
end

candidate_rows = row_ids(time_mask);
selected_rows = candidate_rows(1:min(cfg.static_target_laps, numel(candidate_rows)));
end


function [xq, low_mean, high_mean, low_peak_stats, high_peak_stats] = summarize_low_high_local(low_laps, high_laps)
low_peak_stats.peak_vals = [low_laps.peak_val];
low_peak_stats.mean_peak_v = mean(low_peak_stats.peak_vals, 'omitnan');
low_peak_stats.std_peak_v = std(low_peak_stats.peak_vals, 'omitnan');

high_peak_stats.peak_vals = [high_laps.peak_val];
high_peak_stats.mean_peak_v = mean(high_peak_stats.peak_vals, 'omitnan');
high_peak_stats.std_peak_v = std(high_peak_stats.peak_vals, 'omitnan');

low_x = vertcat(low_laps.x_points);
high_x = vertcat(high_laps.x_points);
low_lo = quantile(low_x, 0.02);
low_hi = quantile(low_x, 0.98);
high_lo = quantile(high_x, 0.02);
high_hi = quantile(high_x, 0.98);
span_lo = max(low_lo, high_lo);
span_hi = min(low_hi, high_hi);
if span_hi <= span_lo
    span_lo = max(min(low_x), min(high_x));
    span_hi = min(max(low_x), max(high_x));
end
if span_hi <= span_lo
    span_lo = min([low_x; high_x]);
    span_hi = max([low_x; high_x]);
end
xq = linspace(span_lo, span_hi, 500).';

low_mat = interpolate_laps_local(low_laps, xq);
high_mat = interpolate_laps_local(high_laps, xq);
low_mean = mean(low_mat, 2, 'omitnan');
high_mean = mean(high_mat, 2, 'omitnan');
end


function mat = interpolate_laps_local(laps, xq)
mat = nan(numel(xq), numel(laps));
for k = 1:numel(laps)
    x = laps(k).x_points(:);
    v = laps(k).v_points(:);
    [x, idx] = sort(x);
    v = v(idx);
    [x_u, ia] = unique(x, 'stable');
    v_u = v(ia);
    mat(:, k) = interp1(x_u, v_u, xq, 'linear', NaN);
end
end


function [t_sec, v] = load_raw_case_channel_local(case_dir, sid, file_id, pinlv)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id));
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / pinlv;
v = raw(:, 2);
end


function F_omega_deg = build_phase_speed_local_local(opr_times, blades_num)
spd_t = opr_times(1:(end - blades_num));
spd_v = 360 ./ max(opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num)), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end


function opr_reference = read_sensor_config_opr_reference_local(Sensor_Config)
opr_reference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if isfield(Sensor_Config, 'OPRReference')
    ref = Sensor_Config.OPRReference;
    if isfield(ref, 'phase_shift_deg')
        opr_reference.phase_shift_deg = ref.phase_shift_deg;
    end
    if isfield(ref, 'phase_shift_mm')
        opr_reference.phase_shift_mm = ref.phase_shift_mm;
    end
    if isfield(ref, 'median_center_minus_start_s')
        opr_reference.median_center_minus_start_s = ref.median_center_minus_start_s;
    elseif isfield(ref, 'center_minus_start_s')
        opr_reference.median_center_minus_start_s = ref.center_minus_start_s;
    end
end
end


function theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference)
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    theta_std = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, blade_id);
    return;
end
theta_std = Sensor_Config.Standard_Relative_Angles(sid, blade_id);
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
if isfield(opr_reference, 'phase_shift_deg') && isfinite(opr_reference.phase_shift_deg)
    theta_std = theta_std - opr_reference.phase_shift_deg;
end
end


function theta_points = map_segment_to_relative_angle_local_local(t_ref, t_seg, F_omega_deg)
if isempty(t_seg)
    theta_points = zeros(size(t_seg));
    return;
end
dt_first = linspace(t_ref, t_seg(1), 10);
theta_base = trapz(dt_first, F_omega_deg(dt_first));
w_seg = F_omega_deg(t_seg);
theta_rel = cumtrapz(t_seg, w_seg);
theta_points = theta_base + theta_rel;
end
