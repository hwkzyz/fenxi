%% Step02_Main_Build_OPRCenterStd_DynamicMap_20250527
% Build the latest dynamic sliding-window map in the OPRCenterStd coordinate.
% Tune parameters here, then run this file directly.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
C = CaseConfig();
validationRoot = fileparts(routeDir);

%% Parameters to tune
P.name.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.name.dynamicSuffix = 'Main20L_W3S1_GradientXRange030_OPRCenterStd';
templateSuffixOverride = strtrim(getenv('STEP02_TEMPLATE_SUFFIX'));
dynamicSuffixOverride = strtrim(getenv('STEP02_DYNAMIC_SUFFIX'));
if ~isempty(templateSuffixOverride)
    P.name.templateSuffix = templateSuffixOverride;
end
if ~isempty(dynamicSuffixOverride)
    P.name.dynamicSuffix = dynamicSuffixOverride;
end

P.path.dynamicDataDir = fullfile('E:\试验数据\20250527\试验20250527', ...
    '20250526_2500-3500_t400');
P.path.caseOutputDir = fullfile(validationRoot, '20250527_low_speed_gap_prior_decoupling', 'legacy', ...
    'output', '20250526_2500-3500_t400');
P.path.sensorConfigFile = fullfile(validationRoot, '20250527_low_speed_gap_prior_decoupling', 'legacy', ...
    'output', '20250526_910_reference', 'Sensor_Config_20250527.mat');

P.step02.analysisSensors = C.sensorIds;
P.step02.analysisStartTime = 1.5;
P.step02.targetLaps = 20;
P.step02.windowLaps = 3;
P.step02.slidingStepLaps = 1;
P.step02.checkPlotWindowID = 1;
P.step02.oldSgResultFile = '';   % Optional, only used for historical comparison.
P.step02.oprChannel = 4;
P.step02.oprPulsesPerRev = 6;
P.step02.sampleRateHz = 5e6;
P.step02.tipRadiusMM = 62.0;
P.step02.pulseWindowSec = 6e-4;
P.step02.pulsePadSec = 2 / P.step02.sampleRateHz;
P.step02.dynamicWindowMode = 'peak_centered_fixed';
P.view.saveFigures = false;

%% Files
[templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20250527(routeDir, C, P.name.templateSuffix);
fprintf('Template source: %s [%s]\n', templateFile, templateSourceMode);

fprintf('\n=== Step02: dynamic sliding-window map ===\n');
step02_build_dynamic_map_embedded(routeDir, P, C, templateFile);

function step02_build_dynamic_map_embedded(rootDir, P, C, template_file)
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

%% Settings
route_dir = rootDir;
old_result_file = strtrim(P.step02.oldSgResultFile);

cfg = struct();
cfg.target_blades = C.bladeId;
cfg.analysis_sensors = P.step02.analysisSensors;
cfg.analysis_start_time = P.step02.analysisStartTime;
cfg.target_laps = P.step02.targetLaps;
cfg.analysis_win_size = P.step02.windowLaps;
cfg.sliding_step = P.step02.slidingStepLaps;
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
dynamic_suffix = sanitize_dynamic_suffix_local(P.name.dynamicSuffix);
template_suffix = P.name.templateSuffix;

output_dir = fullfile(route_dir, 'output');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
figure_dir = fullfile(output_dir, 'figures');
if exist(dynamic_dir, 'dir') ~= 7; mkdir(dynamic_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_%s_SlidingWindows%s_%s.mat', ...
    C.caseTag, dynamic_suffix, C.dataset));

fprintf('\n=== Step02: sliding-window dynamic waveform maps ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
if isempty(old_result_file)
    fprintf('Old SG comparison: skipped\n');
else
    fprintf('Old SG comparison file: %s\n', old_result_file);
end

if ~isfile(template_file)
    error('Template file not found. Run Step01 first: %s', template_file);
end

loaded_template = load(template_file, 'Template');
Template = loaded_template.Template;
hasOldResult = ~isempty(old_result_file) && isfile(old_result_file);
source_result_sensor_tag = sensor_tag;
if hasOldResult
    loaded_result = load(old_result_file, 'Result_Struct');
    OldResult = loaded_result.Result_Struct;
elseif ~isempty(old_result_file)
    warning('Old SG result file not found. Continuing without historical comparison: %s', old_result_file);
    old_result_file = '';
    hasOldResult = false;
end
% Current workspace paths and timing constants are explicit Step02 parameters.
cfg.dynamic_data_dir = P.path.dynamicDataDir;
cfg.case_output_dir = P.path.caseOutputDir;
cfg.sensor_config_file = P.path.sensorConfigFile;
cfg.opr_channel = P.step02.oprChannel;
cfg.opr_pulses_per_rev = P.step02.oprPulsesPerRev;
cfg.pinlv = P.step02.sampleRateHz;
cfg.r_tip_mm = P.step02.tipRadiusMM;
cfg.pulse_window_sec = P.step02.pulseWindowSec;
cfg.pulse_pad_sec = P.step02.pulsePadSec;
cfg.dynamic_window_mode = P.step02.dynamicWindowMode;

if exist(cfg.dynamic_data_dir, 'dir') ~= 7
    error('Dynamic raw-data folder not found: %s', cfg.dynamic_data_dir);
end
if exist(cfg.case_output_dir, 'dir') ~= 7
    error('Step01 legacy output folder not found: %s', cfg.case_output_dir);
end
if exist(cfg.sensor_config_file, 'file') ~= 2
    error('Sensor_Config file not found: %s', cfg.sensor_config_file);
end

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
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, cfg.target_blades, opr_reference);

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
if hasOldResult && strcmp(sensor_tag, source_result_sensor_tag)
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
    fprintf('Old SG comparison skipped.\n');
end

%% Visualization of one sliding window for checking
if isfield(DynamicMap, 'OldSGBestWindowID')
    best_plot_idx = DynamicMap.OldSGBestWindowID;
else
    best_plot_idx = min(max(1, round(P.step02.checkPlotWindowID)), num_windows);
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

if P.view.saveFigures
    exportgraphics(fig, fullfile(figure_dir, sprintf('Step02_DynamicMap_%s_SlidingWindows%s_%s.png', ...
        C.caseTag, dynamic_suffix, C.dataset)), 'Resolution', 300);
end

%% Local helpers kept at the end so the main workflow above stays readable
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
if numel(opr_times) <= num_blades
    theta_rot = nan(size(sample_times));
    return;
end
rev_anchor_times = opr_times(1:num_blades:end);
rev_anchor_times = rev_anchor_times(:);
rev_phase = 2*pi*(0:numel(rev_anchor_times)-1).';
theta_vec = interp1(rev_anchor_times, rev_phase, sample_times(:), 'linear', 'extrap');
theta_rot = reshape(theta_vec, size(sample_times));
theta_rot(sample_times < opr_times(1) | sample_times > opr_times(end)) = nan;
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
end
