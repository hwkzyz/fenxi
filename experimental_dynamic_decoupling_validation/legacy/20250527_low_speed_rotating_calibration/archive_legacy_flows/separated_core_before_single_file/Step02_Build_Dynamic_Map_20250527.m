%% Step02_Build_Dynamic_Map_20250527
% Build sliding-window dynamic waveform maps for the proposed method.
%
% The vibration-region settings are kept the same as the old SG Step5 case:
%   cfg.target_blades        = 1
%   cfg.analysis_sensors     = [1 3 6]
%   cfg.analysis_start_time  = 1.5 s
%   cfg.target_laps          = 20
%   cfg.analysis_win_size    = 3 laps
%   cfg.sliding_step         = 1 lap
%
% This script only prepares waveform windows. The identification model is
% changed in Step03.

clear; clc; close all;

%% Settings
core_dir = fileparts(mfilename('fullpath'));
route_dir = fileparts(core_dir);
old_result_file = ...
    'E:\0小论文+程序\0博士期间小论文+程序\7超高斯模型-权重-瞬态\程序\直叶片验证\实验验证\超高斯\20250527适配\Result_single_sync_20250527_B1_S136_20250526_2500_3500_t400_Start1p5s.mat';

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1, 3, 6];
cfg.analysis_start_time = 1.5;
cfg.target_laps = parse_positive_integer_env_local('STEP02_TARGET_LAPS', 20);
cfg.analysis_win_size = parse_positive_integer_env_local('STEP02_ANALYSIS_WIN_SIZE', 3);
cfg.sliding_step = parse_positive_integer_env_local('STEP02_SLIDING_STEP', 1);
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
source_result_sensor_tag = 'S136';
dynamic_suffix = sanitize_dynamic_suffix_local(strtrim(getenv('STEP02_DYNAMIC_SUFFIX')));

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
figure_dir = fullfile(output_dir, 'figures');
if exist(dynamic_dir, 'dir') ~= 7; mkdir(dynamic_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_suffix = strtrim(getenv('STEP02_TEMPLATE_SUFFIX'));
if isempty(template_suffix)
    template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', ...
        cfg.target_blades, sensor_tag));
else
    template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s_%s_20250527.mat', ...
        cfg.target_blades, sensor_tag, template_suffix));
end
dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows%s_20250527.mat', ...
    cfg.target_blades, sensor_tag, dynamic_suffix));

fprintf('\n=== Step02: sliding-window dynamic waveform maps ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('Old SG result: %s\n', old_result_file);

if ~isfile(template_file)
    error('Template file not found. Run Step01 first: %s', template_file);
end
if ~isfile(old_result_file)
    error('Old SG result file not found: %s', old_result_file);
end

loaded_template = load(template_file, 'Template');
Template = loaded_template.Template;
loaded_result = load(old_result_file, 'Result_Struct');
OldResult = loaded_result.Result_Struct;
current_cfg = build_single_sync_experiment_config_20250527();

% Use the current workspace config for OPR timing and raw-data paths.  The
% old SG result is retained only as a trend/window reference.
cfg.dynamic_data_dir = current_cfg.dynamic_data_dir;
cfg.case_output_dir = current_cfg.case_output_dir;
cfg.sensor_config_file = current_cfg.sensor_config_file;
cfg.opr_channel = current_cfg.opr_channel;
cfg.opr_pulses_per_rev = current_cfg.opr_pulses_per_rev;
cfg.pinlv = current_cfg.pinlv;
cfg.r_tip_mm = current_cfg.r_tip_mm;
cfg.pulse_window_sec = current_cfg.pulse_window_sec;
cfg.pulse_pad_sec = current_cfg.pulse_pad_sec;
cfg.dynamic_window_mode = current_cfg.dynamic_window_mode;

loaded_sensor_config = load(cfg.sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_sensor_config.Sensor_Config;
loaded_opr = load(fullfile(cfg.case_output_dir, 'jiluOPR.mat'), 'jiluOPR');
opr_times = loaded_opr.jiluOPR(:, 1);
opr_reference = build_opr_reference_from_jilu_local(loaded_opr.jiluOPR, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
F_omega_deg = build_phase_speed_local(opr_times, cfg.opr_pulses_per_rev);
fprintf('Dynamic OPR center shift: %.4f deg, %.4f mm from rising edge.\n', ...
    opr_reference.phase_shift_deg, opr_reference.phase_shift_mm);

%% Select the same 20 target-blade laps after the same start time
selection = repmat(struct('sensor_id', NaN, 'selected_rows', []), numel(cfg.analysis_sensors), 1);
global_window = [inf, -inf];
for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_file = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;

    blade_mask = jilublade(:, 4) == cfg.target_blades & jilublade(:, 3) >= cfg.analysis_start_time;
    blade_rows = find(blade_mask);
    if numel(blade_rows) < cfg.target_laps
        error('CH%d has only %d target-blade laps after %.3f s.', sid, numel(blade_rows), cfg.analysis_start_time);
    end
    blade_rows = blade_rows(1:cfg.target_laps);

    selection(is).sensor_id = sid;
    selection(is).selected_rows = blade_rows(:);
    global_window(1) = min(global_window(1), min(jilublade(blade_rows, 1)) - cfg.pulse_window_sec);
    global_window(2) = max(global_window(2), max(jilublade(blade_rows, 2)) + cfg.pulse_window_sec);
end
fprintf('Global raw-data window: %.6f-%.6f s\n', global_window(1), global_window(2));

%% Load raw dynamic voltage streams only once
file_ranges = build_dynamic_file_ranges_local(cfg.dynamic_data_dir, cfg.opr_channel, cfg.pinlv);
selected_file_mask = [file_ranges.t_end] >= global_window(1) & [file_ranges.t_start] <= global_window(2);
selected_file_ranges = file_ranges(selected_file_mask);
if isempty(selected_file_ranges)
    error('No raw dynamic files overlap the requested window.');
end

raw_stream(max(cfg.analysis_sensors)) = struct('T', [], 'V', []);
for ir = 1:numel(selected_file_ranges)
    file_id = selected_file_ranges(ir).file_id;
    offset = selected_file_ranges(ir).offset;
    for sid = cfg.analysis_sensors
        [t_local, v_local] = load_raw_case_channel_local(cfg.dynamic_data_dir, sid, file_id, cfg.pinlv);
        if isempty(t_local)
            continue;
        end
        t_global = t_local(:) + offset;
        keep = t_global >= global_window(1) & t_global <= global_window(2);
        raw_stream(sid).T = [raw_stream(sid).T; t_global(keep)];
        raw_stream(sid).V = [raw_stream(sid).V; v_local(keep)];
    end
end

%% Extract each selected lap waveform once
LapData = repmat(struct('sensor_id', NaN, 'Lap', []), numel(cfg.analysis_sensors), 1);
for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_file = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;
    theta_std = Sensor_Config.Standard_Relative_Angles(sid, cfg.target_blades);
    theta_std = convert_standard_angle_to_opr_center_local(theta_std, opr_reference);

    Lap = repmat(struct('lap_id', NaN, 't', [], 'x_abs', [], 'V', [], 'theta', []), cfg.target_laps, 1);
    for lap_id = 1:cfg.target_laps
        row_id = selection(is).selected_rows(lap_id);
        t_peak = jilublade(row_id, 3);
        if strcmpi(cfg.dynamic_window_mode, 'legacy_row_bounds')
            t_start = jilublade(row_id, 1) - cfg.pulse_pad_sec;
            t_end = jilublade(row_id, 2) + cfg.pulse_pad_sec;
        else
            t_start = t_peak - cfg.pulse_window_sec;
            t_end = t_peak + cfg.pulse_window_sec;
        end

        mask = raw_stream(sid).T >= t_start & raw_stream(sid).T <= t_end;
        t_seg = raw_stream(sid).T(mask);
        v_seg = raw_stream(sid).V(mask);
        idx_prev = find(opr_times < t_peak, 1, 'last');
        if isempty(idx_prev) || numel(t_seg) < 5
            continue;
        end

        theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
        theta_diff_deg = mod(theta_points_deg - theta_std + 180, 360) - 180;
        x_abs = theta_diff_deg * (pi / 180) * cfg.r_tip_mm;
        theta_rot = map_time_to_rotor_phase_local(opr_times, t_seg, cfg.opr_pulses_per_rev);
        valid = isfinite(theta_rot);

        Lap(lap_id).lap_id = lap_id;
        Lap(lap_id).t = t_seg(valid);
        Lap(lap_id).x_abs = x_abs(valid);
        Lap(lap_id).V = v_seg(valid);
        Lap(lap_id).theta = theta_rot(valid);
    end

    LapData(is).sensor_id = sid;
    LapData(is).Lap = Lap;
end

%% Build all sliding-window maps
num_windows = floor((cfg.target_laps - cfg.analysis_win_size) / cfg.sliding_step) + 1;
Window = repmat(struct( ...
    'window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, ...
    'Sensor', []), num_windows, 1);

for w_idx = 1:num_windows
    lap_start = 1 + (w_idx - 1) * cfg.sliding_step;
    lap_end = lap_start + cfg.analysis_win_size - 1;
    lap_range = lap_start:lap_end;
    Sensor = repmat(struct( ...
        'sensor_id', NaN, 't', [], 'x_abs', [], 'x_rel', [], ...
        'V', [], 'W', [], 'theta', [], 'point_count', NaN), numel(cfg.analysis_sensors), 1);

    all_t_window = [];
    for is = 1:numel(cfg.analysis_sensors)
        sid = cfg.analysis_sensors(is);
        Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
        t = [];
        x_abs = [];
        V = [];
        theta = [];
        for lap_id = lap_range
            D = LapData(is).Lap(lap_id);
            t = [t; D.t(:)]; %#ok<AGROW>
            x_abs = [x_abs; D.x_abs(:)]; %#ok<AGROW>
            V = [V; D.V(:)]; %#ok<AGROW>
            theta = [theta; D.theta(:)]; %#ok<AGROW>
        end
        x_rel = x_abs - Tpl.xc;
        W = build_simple_waveform_weight_local(V);

        Sensor(is).sensor_id = sid;
        Sensor(is).t = t(:);
        Sensor(is).x_abs = x_abs(:);
        Sensor(is).x_rel = x_rel(:);
        Sensor(is).V = V(:);
        Sensor(is).W = W(:);
        Sensor(is).theta = theta(:);
        Sensor(is).point_count = numel(t);
        all_t_window = [all_t_window; t(:)]; %#ok<AGROW>
    end

    Window(w_idx).window_id = w_idx;
    Window(w_idx).lap_range = lap_range;
    Window(w_idx).time_window = [min(all_t_window), max(all_t_window)];
    Window(w_idx).rot_freq_mean_hz = compute_local_rot_freq_local(opr_times, cfg.opr_pulses_per_rev, Window(w_idx).time_window);
    Window(w_idx).rot_rpm_mean = 60 * Window(w_idx).rot_freq_mean_hz;
    Window(w_idx).Sensor = Sensor;
end

DynamicMap = struct();
DynamicMap.Route = 'low_speed_template_sliding_window_dynamic_map_no_gap';
DynamicMap.TargetBlade = cfg.target_blades;
DynamicMap.SensorIDs = cfg.analysis_sensors;
DynamicMap.SensorTag = sensor_tag;
DynamicMap.SourceSettings = cfg;
DynamicMap.SourceResultFile = old_result_file;
DynamicMap.SourceResultSensorTag = source_result_sensor_tag;
DynamicMap.TemplateFileForXRel = template_file;
DynamicMap.TemplateSuffixForXRel = template_suffix;
DynamicMap.XCenterBySensor = build_xcenter_table_local(Template, cfg.analysis_sensors);
DynamicMap.OPRReference = opr_reference;
if strcmp(sensor_tag, source_result_sensor_tag)
    DynamicMap.OldSGTrends = OldResult.Trends;
    DynamicMap.OldSGBestWindowID = OldResult.BestWindow.window_id;
    DynamicMap.OldSGBestLapRange = OldResult.BestWindow.lap_range;
end
DynamicMap.GlobalTimeWindow = global_window;
DynamicMap.SelectedRawFileIDs = [selected_file_ranges.file_id];
DynamicMap.Window = Window;

save(dynamic_file, 'DynamicMap', '-v7.3');
fprintf('Saved sliding-window dynamic map: %s\n', dynamic_file);
fprintf('Built %d windows.\n', num_windows);
if isfield(DynamicMap, 'OldSGBestWindowID')
    fprintf('Old SG best window: %d, laps %s\n', ...
        DynamicMap.OldSGBestWindowID, mat2str(DynamicMap.OldSGBestLapRange));
else
    fprintf('Old SG comparison skipped: source result is %s, current map is %s.\n', ...
        source_result_sensor_tag, sensor_tag);
end

%% Visualization of one sliding window for checking
if isfield(DynamicMap, 'OldSGBestWindowID')
    best_plot_idx = DynamicMap.OldSGBestWindowID;
else
    best_plot_idx = 1;
end
fig = figure('Name', 'Step02 sliding-window dynamic map check', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
tiledlayout(fig, numel(cfg.analysis_sensors), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = DynamicMap.Window(best_plot_idx).Sensor(is);

    nexttile;
    plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.8, 'DisplayName', 'low-speed spline template'); hold on;
    scatter(D.x_rel, D.V, 8, D.t, 'filled', 'DisplayName', 'dynamic samples');
    xlabel('x relative to template center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d dynamic samples, window %d, laps %s', ...
        sid, best_plot_idx, mat2str(DynamicMap.Window(best_plot_idx).lap_range)));
    cb = colorbar;
    cb.Label.String = 'Time (s)';
    box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');
    legend('Location', 'best');
end

exportgraphics(fig, fullfile(figure_dir, sprintf('Step02_DynamicMap_B%d_%s_SlidingWindows%s.png', ...
    cfg.target_blades, sensor_tag, dynamic_suffix)), 'Resolution', 300);

%% Local helpers kept at the end so the main workflow above stays readable
function v = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    error('%s must be a positive integer.', name);
end
v = max(1, floor(tmp));
end

function suffix = sanitize_dynamic_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function T = build_xcenter_table_local(Template, sensor_ids)
sensor_id = sensor_ids(:);
xc_mm = nan(numel(sensor_id), 1);
for i = 1:numel(sensor_id)
    idx = find([Template.Sensor.sensor_id] == sensor_id(i), 1, 'first');
    if ~isempty(idx) && isfield(Template.Sensor(idx), 'xc')
        xc_mm(i) = Template.Sensor(idx).xc;
    end
end
T = table(sensor_id, xc_mm);
end

function theta_std_center = convert_standard_angle_to_opr_center_local(theta_std_start, opr_reference)
theta_std_center = theta_std_start;
if isstruct(opr_reference) && isfield(opr_reference, 'phase_shift_deg') && ...
        isfinite(opr_reference.phase_shift_deg)
    theta_std_center = theta_std_start - opr_reference.phase_shift_deg;
end
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
if size(jiluOPR, 2) >= 3
    start_time = jiluOPR(:, 2);
else
    start_time = jiluOPR(:, 2);
    if median(center_time - start_time, 'omitnan') < 0
        warning(['Two-column jiluOPR appears to store [start,end], not [center,start]. ' ...
            'OPR center phase shift is left at zero.']);
        return;
    end
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

function file_ranges = build_dynamic_file_ranges_local(case_dir, opr_channel, pinlv)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', opr_channel)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(opr_channel) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        file_ids(i) = str2double(tok{1});
    end
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));

file_ranges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(file_ids), 1);
last_end = [];
for i = 1:numel(file_ids)
    file_id = file_ids(i);
    [t_opr, ~] = load_raw_case_channel_local(case_dir, opr_channel, file_id, pinlv);
    if isempty(t_opr)
        continue;
    end
    if isempty(last_end)
        offset = 0;
    else
        offset = last_end + 1 / pinlv - t_opr(1);
    end
    last_end = t_opr(end) + offset;
    file_ranges(i).file_id = file_id;
    file_ranges(i).offset = offset;
    file_ranges(i).t_start = t_opr(1) + offset;
    file_ranges(i).t_end = t_opr(end) + offset;
end
end

function [t_sec, v] = load_raw_case_channel_local(case_dir, sid, file_id, pinlv)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id));
if ~isfile(filepath)
    t_sec = [];
    v = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function F_omega_deg = build_phase_speed_local(opr_times, blades_num)
spd_t = opr_times(1:end-blades_num);
spd_v = 360 ./ max(opr_times(blades_num+1:end) - opr_times(1:end-blades_num), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end

function theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg)
dt_first = linspace(t_ref, t_seg(1), 10);
theta_base = trapz(dt_first, F_omega_deg(dt_first));
w_seg = F_omega_deg(t_seg);
theta_rel = cumtrapz(t_seg, w_seg);
theta_points = theta_base + theta_rel;
end

function theta_rot = map_time_to_rotor_phase_local(opr_times, sample_times, num_blades)
theta_rot = nan(size(sample_times));
for i = 1:numel(opr_times)-1
    t0 = opr_times(i);
    t1 = opr_times(i+1);
    mask = sample_times >= t0 & sample_times <= t1;
    if ~any(mask)
        continue;
    end
    frac = (sample_times(mask) - t0) ./ max(t1 - t0, eps);
    theta_rot(mask) = (2*pi/num_blades) * ((i - 1) + frac);
end
end

function W = build_simple_waveform_weight_local(V)
if isempty(V)
    W = [];
    return;
end
V = V(:);
v_floor = prctile(V, 5);
v_peak = prctile(V, 99);
span = max(v_peak - v_floor, eps);
W = (V - v_floor) ./ span;
W = min(max(W, 0.05), 1.0);
end

function rot_freq_hz = compute_local_rot_freq_local(opr_times, pulses_per_rev, time_window)
mask = opr_times >= time_window(1) & opr_times <= time_window(2);
t = opr_times(mask);
if numel(t) > pulses_per_rev
    rot_freq_hz = median(1 ./ max(t(1+pulses_per_rev:end) - t(1:end-pulses_per_rev), eps), 'omitnan');
else
    rot_freq_hz = NaN;
end
end
