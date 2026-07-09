function experiment_case = prepare_identification_case_single_sync_20251222(cfg, show_plots, save_output)
%PREPARE_IDENTIFICATION_CASE_SINGLE_SYNC_20251222
% Build a compact Step5 case around the requested identification time window.

if nargin < 1 || isempty(cfg)
    cfg = build_single_sync_experiment_config_20251222();
end
if nargin < 2
    show_plots = cfg.show_figures;
end
if nargin < 3
    save_output = isfield(cfg, 'save_experiment_case') && cfg.save_experiment_case;
end

if ~isfile(cfg.sensor_config_file)
    Step1_Build_Sensor_Config_20251222(false);
end
if ~isfile(fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', cfg.analysis_sensors(1))))
    Step2_Extract_JiluBlade_20251222(cfg.dynamic_case_name, false);
end

loaded_cfg = load(cfg.sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;

opr_loaded = load(fullfile(cfg.case_output_dir, 'jiluOPR.mat'), 'jiluOPR');
jiluOPR = opr_loaded.jiluOPR;
opr_times = jiluOPR(:, 1);
opr_reference = build_opr_reference_from_jilu_local(jiluOPR, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
omega_loaded = load(fullfile(cfg.case_output_dir, 'omega.mat'));
f_omega_deg = build_phase_speed_local(opr_times, cfg.opr_pulses_per_rev);

[selection, global_window] = select_identification_rows_local(cfg, Sensor_Config);
file_ranges = build_dynamic_file_ranges_local(cfg);
[raw_stream, selected_file_ids] = load_case_raw_streams_window_local(cfg, file_ranges, global_window);

experiment_case = struct();
experiment_case.Config = cfg;
experiment_case.Sensor_Config = Sensor_Config;
experiment_case.Raw_Stream = raw_stream;
experiment_case.Global_OPR_Times = opr_times(:);
experiment_case.jiluOPR = jiluOPR;
experiment_case.omega = omega_loaded;
experiment_case.Extracted_Data = struct([]);
experiment_case.Diagnostics = struct();

for i = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(i);
    probe_path = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_path, 'jilublade');
    jilublade = loaded_probe.jilublade;

    blade_rows = selection(i).selected_rows;
    laps = struct('lap_id', {}, 'row_id', {}, 'x_points', {}, 't_points', {}, ...
        'v_points', {}, 'peak_val', {}, 'peak_time', {});
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, cfg.target_blade, opr_reference);
    peak_features = nan(numel(blade_rows), 1);

    for k = 1:numel(blade_rows)
        row_id = blade_rows(k);
        t_peak = jilublade(row_id, 3);
        [t_start, t_end] = build_dynamic_segment_window_local(cfg, jilublade, row_id, t_peak);
        mask = raw_stream(sid).T >= t_start & raw_stream(sid).T <= t_end;
        if nnz(mask) < 5
            continue;
        end

        t_seg = raw_stream(sid).T(mask);
        v_seg = raw_stream(sid).V(mask);
        idx_prev = find(opr_times < t_peak, 1, 'last');
        if isempty(idx_prev)
            continue;
        end

        t_ref = opr_times(idx_prev);
        theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, f_omega_deg);
        theta_diff = mod(theta_points - theta_std + 180, 360) - 180;
        x_points_mm = theta_diff * (pi / 180) * cfg.r_tip_mm;

        laps(end + 1).lap_id = k; %#ok<AGROW>
        laps(end).row_id = row_id;
        laps(end).x_points = x_points_mm(:);
        laps(end).t_points = t_seg(:);
        laps(end).v_points = v_seg(:);
        [peak_features(k), idx_pk] = max(v_seg);
        laps(end).peak_val = peak_features(k);
        if strcmpi(cfg.dynamic_window_mode, 'peak_centered_fixed')
            laps(end).peak_time = t_peak;
        else
            laps(end).peak_time = t_seg(idx_pk);
        end
    end

    experiment_case.Extracted_Data(end + 1).sensor_id = sid; %#ok<AGROW>
    experiment_case.Extracted_Data(end).selected_rows = blade_rows(:);
    experiment_case.Extracted_Data(end).Laps = laps;
    experiment_case.Extracted_Data(end).fingerprint_corr = NaN;
    experiment_case.Extracted_Data(end).peak_features = peak_features(:);
end

experiment_case.Diagnostics.opr_pulses_per_rev = cfg.opr_pulses_per_rev;
experiment_case.Diagnostics.target_laps = cfg.target_laps;
experiment_case.Diagnostics.analysis_start_time = cfg.analysis_start_time;
experiment_case.Diagnostics.time_window = global_window;
experiment_case.Diagnostics.file_ids_loaded = selected_file_ids(:).';
experiment_case.Diagnostics.OPRReference = opr_reference;
experiment_case.Diagnostics.standard_angle_reference = 'opr_pulse_center';
experiment_case.Diagnostics.standard_angle_phase_shift_deg = opr_reference.phase_shift_deg;
experiment_case.Diagnostics.standard_angle_phase_shift_mm = opr_reference.phase_shift_mm;
experiment_case.Diagnostics.speed_diagnostics = compute_opr_speed_diagnostics_20251222( ...
    opr_times, cfg.opr_pulses_per_rev, global_window);
if isfield(omega_loaded, 'omega_rpm')
    experiment_case.Diagnostics.rot_freq_global_hz = median(omega_loaded.omega_rpm, 'omitnan') / 60;
else
    experiment_case.Diagnostics.rot_freq_global_hz = ...
        experiment_case.Diagnostics.speed_diagnostics.global_rot_freq_mean_hz;
end
experiment_case.Diagnostics.rot_freq_mean_hz = ...
    experiment_case.Diagnostics.speed_diagnostics.local_rot_freq_mean_hz;
if ~isfinite(experiment_case.Diagnostics.rot_freq_mean_hz)
    experiment_case.Diagnostics.rot_freq_mean_hz = experiment_case.Diagnostics.rot_freq_global_hz;
end

if save_output
    save(cfg.case_file, 'experiment_case', '-v7.3');
    fprintf('Saved compact Step5 case to %s\n', cfg.case_file);
else
    fprintf('Prepared compact Step5 case in memory for Blade %d.\n', cfg.target_blade);
end

if show_plots
    plot_identification_case_local(experiment_case, cfg);
end
end


function [selection, global_window] = select_identification_rows_local(cfg, Sensor_Config)
selection = repmat(struct( ...
    'sensor_id', NaN, ...
    'selected_rows', [], ...
    'window_start', NaN, ...
    'window_end', NaN), numel(cfg.analysis_sensors), 1);
global_window = [inf, -inf];

for i = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(i);
    probe_path = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    if ~isfile(probe_path)
        error('Missing Step2 output: %s', probe_path);
    end

    loaded_probe = load(probe_path, 'jilublade');
    jilublade = loaded_probe.jilublade;

    blade_mask = jilublade(:, 4) == cfg.target_blade & jilublade(:, 3) >= cfg.analysis_start_time;
    blade_rows = find(blade_mask);
    if isempty(blade_rows)
        error('No target-blade rows found for CH%d after %.3f s.', sid, cfg.analysis_start_time);
    end
    blade_rows = blade_rows(1:min(cfg.target_laps, numel(blade_rows)));

    selection(i).sensor_id = sid;
    selection(i).selected_rows = blade_rows(:);
    if isfield(cfg, 'dynamic_window_mode') && strcmpi(cfg.dynamic_window_mode, 'peak_centered_fixed')
        selection(i).window_start = min(jilublade(blade_rows, 3)) - cfg.pulse_window_sec;
        selection(i).window_end = max(jilublade(blade_rows, 3)) + cfg.pulse_window_sec;
    else
        selection(i).window_start = min(jilublade(blade_rows, 1)) - cfg.pulse_pad_sec;
        selection(i).window_end = max(jilublade(blade_rows, 2)) + cfg.pulse_pad_sec;
    end

    global_window(1) = min(global_window(1), selection(i).window_start);
    global_window(2) = max(global_window(2), selection(i).window_end);

    if ~isfinite(Sensor_Config.Standard_Relative_Angles(sid, cfg.target_blade))
        error('Standard relative angle is missing for CH%d blade %d.', sid, cfg.target_blade);
    end
end
end


function file_ranges = build_dynamic_file_ranges_local(cfg)
file_ids = list_case_file_ids_local(cfg.dynamic_data_dir, cfg.opr_channel);
file_ranges = repmat(struct( ...
    'file_id', NaN, ...
    'offset', NaN, ...
    't_start', NaN, ...
    't_end', NaN, ...
    'is_selected', false), numel(file_ids), 1);

last_end = [];
for i = 1:numel(file_ids)
    file_id = file_ids(i);
    [t_opr, ~] = load_raw_case_channel_local(cfg.dynamic_data_dir, cfg.opr_channel, file_id, cfg.pinlv);
    if isempty(t_opr)
        continue;
    end

    if isempty(last_end)
        offset = 0;
    else
        offset = last_end + 1 / cfg.pinlv - t_opr(1);
    end
    last_end = t_opr(end) + offset;

    file_ranges(i).file_id = file_id;
    file_ranges(i).offset = offset;
    file_ranges(i).t_start = t_opr(1) + offset;
    file_ranges(i).t_end = t_opr(end) + offset;
end
end


function [raw_stream, selected_file_ids] = load_case_raw_streams_window_local(cfg, file_ranges, time_window)
select_mask = [file_ranges.t_end] >= time_window(1) & [file_ranges.t_start] <= time_window(2);
file_ranges = file_ranges(select_mask);
if isempty(file_ranges)
    error('No dynamic files overlap the requested identification window.');
end

selected_file_ids = [file_ranges.file_id];

max_sid = max(cfg.analysis_sensors);
raw_stream(max_sid) = struct('T', [], 'V', []);

for i = 1:numel(file_ranges)
    file_id = file_ranges(i).file_id;
    offset = file_ranges(i).offset;
    for sid = cfg.analysis_sensors
        filepath = fullfile(cfg.dynamic_data_dir, sprintf('4-%d-%d.mat', sid, file_id));
        if ~isfile(filepath)
            continue;
        end
        [t_local, v_local] = load_raw_case_channel_local(cfg.dynamic_data_dir, sid, file_id, cfg.pinlv);
        t_global = t_local(:) + offset;
        keep_mask = t_global >= time_window(1) & t_global <= time_window(2);
        if ~any(keep_mask)
            continue;
        end
        raw_stream(sid).T = [raw_stream(sid).T; t_global(keep_mask)]; %#ok<AGROW>
        raw_stream(sid).V = [raw_stream(sid).V; v_local(keep_mask)]; %#ok<AGROW>
    end
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


function file_ids = list_case_file_ids_local(case_dir, channel_id)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        file_ids(i) = str2double(tok{1});
    end
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));
end


function F_omega_deg = build_phase_speed_local(opr_times, blades_num)
spd_t = opr_times(1:(end - blades_num));
spd_v = 360 ./ max(opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num)), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end


function theta_std_center = convert_standard_angle_to_opr_center_local(theta_std_start, opr_reference)
theta_std_center = theta_std_start;
if isstruct(opr_reference) && isfield(opr_reference, 'phase_shift_deg') && ...
        isfinite(opr_reference.phase_shift_deg)
    theta_std_center = theta_std_start - opr_reference.phase_shift_deg;
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
theta_std = convert_standard_angle_to_opr_center_local(theta_std, opr_reference);
end


function opr_reference = build_opr_reference_from_jilu_local(jiluOPR, pulses_per_rev, r_tip_mm)
opr_reference = struct('mode', 'multi_threshold_center', ...
    'standard_angle_reference', 'opr_pulse_center', ...
    'phase_shift_deg', 0, ...
    'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if size(jiluOPR, 2) < 2 || pulses_per_rev < 1
    return;
end
center_time = jiluOPR(:, 1);
start_time = jiluOPR(:, 2);
if size(jiluOPR, 2) == 2 && median(center_time - start_time, 'omitnan') < 0
    warning(['Two-column jiluOPR appears to store [start,end], not [center,start]. ' ...
        'OPR center phase shift is left at zero.']);
    return;
end
n = min(numel(center_time) - pulses_per_rev, numel(start_time));
if n < 1
    return;
end
dt_center = center_time(1:n) - start_time(1:n);
dt_rev = center_time((1:n) + pulses_per_rev) - center_time(1:n);
valid = isfinite(dt_center) & isfinite(dt_rev) & dt_rev > eps;
if ~any(valid)
    return;
end
shift_deg = 360 * dt_center(valid) ./ dt_rev(valid);
opr_reference.phase_shift_deg = median(shift_deg, 'omitnan');
opr_reference.phase_shift_mm = opr_reference.phase_shift_deg * (pi / 180) * r_tip_mm;
opr_reference.median_center_minus_start_s = median(dt_center(valid), 'omitnan');
end


function theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg)
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


function [t_start, t_end] = build_dynamic_segment_window_local(cfg, jilublade, row_id, t_peak)
if isfield(cfg, 'dynamic_window_mode') && strcmpi(cfg.dynamic_window_mode, 'legacy_row_bounds')
    t_start = jilublade(row_id, 1) - cfg.pulse_pad_sec;
    t_end = jilublade(row_id, 2) + cfg.pulse_pad_sec;
else
    t_start = t_peak - cfg.pulse_window_sec;
    t_end = t_peak + cfg.pulse_window_sec;
end
end


function plot_identification_case_local(experiment_case, cfg)
fig = figure('Name', sprintf('Step5 Compact Case - Blade %d', cfg.target_blade), ...
    'Color', 'w', 'Position', [120, 90, 1450, 860], 'NumberTitle', 'off');
tiledlayout(fig, numel(cfg.analysis_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(i);
    data_idx = find([experiment_case.Extracted_Data.sensor_id] == sid, 1, 'first');

    nexttile;
    hold on; grid on; box on;
    plot(experiment_case.Raw_Stream(sid).T, experiment_case.Raw_Stream(sid).V, '-', ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d compact raw window', sid));

    nexttile;
    hold on; grid on; box on;
    if ~isempty(data_idx)
        laps = experiment_case.Extracted_Data(data_idx).Laps;
        cmap = lines(max(1, numel(laps)));
        for k = 1:numel(laps)
            plot(laps(k).x_points, laps(k).v_points, '-', 'Color', cmap(k, :), 'LineWidth', 1.2);
        end
    end
    xlabel('Restored position (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d mapped compact laps', sid));
end
end

