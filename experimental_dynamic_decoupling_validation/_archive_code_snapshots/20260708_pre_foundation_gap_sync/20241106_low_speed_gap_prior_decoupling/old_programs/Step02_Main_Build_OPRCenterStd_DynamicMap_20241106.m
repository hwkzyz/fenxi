%% Step02_Main_Build_OPRCenterStd_DynamicMap_20241106
% Build the latest dynamic sliding-window map in the OPRCenterStd coordinate.
% Tune parameters here, then run this file directly.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
legacyDir = fullfile(validationRoot, '20241106_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacyDir, 'dir') ~= 7
    error('Legacy 20241106 helper folder not found: %s', legacyDir);
end
addpath(legacyDir);

baseCfg = Get_20241106_BTT_Config();
syncCfg = build_single_sync_experiment_config_20241106();
dynamicCase = baseCfg.dynamic_cases{1};

%% Parameters to tune
P.name.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.name.dynamicSuffix = 'Main20L_W3S1_GradientXRange030_OPRCenterStd';

P.path.dynamicDataDir = fullfile(baseCfg.dataset_root, dynamicCase);
P.path.caseOutputDir = fullfile(baseCfg.output_root, dynamicCase);
P.path.sensorConfigFile = fullfile(baseCfg.reference_output_dir, 'Sensor_Config_20241106.mat');
P.path.regionPlanFile = fullfile(routeDir, 'outputs', 'Step05_BTT_WindowPlan_20241106.csv');

P.step02.targetBlade = 2;
P.step02.analysisSensors = [5 7];
P.step02.analysisStartTime = 55.3999459833259;
P.step02.targetLaps = 20;
P.step02.windowLaps = 3;
P.step02.slidingStepLaps = 1;
P.step02.checkPlotWindowID = 1;
P.step02.useRegionPlan = true;
P.step02.useRegionBlade = false;
P.step02.regionId = 2;
P.step02.regionStartMode = 'peak';  % 'peak' uses peakBttTimeSec; 'start' uses bttStartSec.
P.step02.oldSgResultFile = '';   % Optional, only used for historical comparison.
P.step02.oprChannel = baseCfg.opr_id;
P.step02.oprPulsesPerRev = baseCfg.opr_pulses_per_rev;
P.step02.sampleRateHz = baseCfg.pinlv;
P.step02.tipRadiusMM = baseCfg.r_tip_mm;
P.step02.pulseWindowSec = syncCfg.pulse_window_sec;
P.step02.pulsePadSec = syncCfg.pulse_pad_sec;
P.step02.dynamicWindowMode = syncCfg.dynamic_window_mode;
P.view.saveFigures = true;

targetBladeEnv = str2double(strtrim(getenv('STEP20241106_TARGET_BLADE')));
if isfinite(targetBladeEnv) && targetBladeEnv >= 1
    P.step02.targetBlade = round(targetBladeEnv);
end
analysisSensorsEnv = sscanf(strtrim(getenv('STEP20241106_ANALYSIS_SENSORS')), '%d').';
if ~isempty(analysisSensorsEnv)
    P.step02.analysisSensors = unique(analysisSensorsEnv, 'stable');
end
regionIdEnv = str2double(strtrim(getenv('STEP20241106_REGION_ID')));
if isfinite(regionIdEnv) && regionIdEnv >= 1
    P.step02.regionId = round(regionIdEnv);
end
regionStartModeEnv = strtrim(getenv('STEP20241106_REGION_START_MODE'));
if ~isempty(regionStartModeEnv)
    P.step02.regionStartMode = regionStartModeEnv;
end
useRegionPlanEnv = strtrim(getenv('STEP20241106_USE_REGION_PLAN'));
if ~isempty(useRegionPlanEnv)
    useRegionPlanText = lower(strtrim(useRegionPlanEnv));
    P.step02.useRegionPlan = ismember(useRegionPlanText, {'1', 'true', 'yes', 'y', 'on'});
end
analysisStartEnv = str2double(strtrim(getenv('STEP20241106_ANALYSIS_START_TIME')));
if isfinite(analysisStartEnv)
    P.step02.analysisStartTime = analysisStartEnv;
end
dynamicSuffixEnv = strtrim(getenv('STEP20241106_DYNAMIC_SUFFIX'));
if ~isempty(dynamicSuffixEnv)
    P.name.dynamicSuffix = dynamicSuffixEnv;
end

%% Files
sensorTag = ['S', sprintf('%d', P.step02.analysisSensors)];
templateFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('Template_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    P.step02.targetBlade, sensorTag, P.name.templateSuffix));
if exist(templateFile, 'file') ~= 2
    explain_missing_template_local(templateFile, routeDir, P.name.templateSuffix);
end

fprintf('\n=== Step02: dynamic sliding-window map ===\n');
step02_build_dynamic_map_embedded(routeDir, P);

function step02_build_dynamic_map_embedded(rootDir, P)
%% Step02_Build_Dynamic_Map_20241106
% Build sliding-window dynamic waveform maps for the proposed method.
%
% The vibration-region settings are kept the same as the old SG Step5 case:
%   cfg.target_blades        = 1
%   cfg.analysis_sensors     = [1 2 3]
%   cfg.analysis_start_time  = 50.2 s
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
cfg.target_blades = P.step02.targetBlade;
cfg.analysis_sensors = P.step02.analysisSensors;
cfg.analysis_start_time = P.step02.analysisStartTime;
cfg.target_laps = P.step02.targetLaps;
cfg.analysis_win_size = P.step02.windowLaps;
cfg.sliding_step = P.step02.slidingStepLaps;
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
source_result_sensor_tag = sensor_tag;
dynamic_suffix = sanitize_dynamic_suffix_local(P.name.dynamicSuffix);
selection_info = load_region_selection_local(P.path.regionPlanFile, P.step02);
if selection_info.use_region_plan
    cfg.target_blades = selection_info.target_blade;
    cfg.analysis_start_time = selection_info.start_time_sec;
end

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
figure_dir = fullfile(output_dir, 'figures');
if exist(dynamic_dir, 'dir') ~= 7; mkdir(dynamic_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_suffix = P.name.templateSuffix;
if isempty(template_suffix)
    template_file = fullfile(template_dir, sprintf( ...
        'Template_LowSpeedRotating_B%d_%s_20241106.mat', cfg.target_blades, sensor_tag));
else
    template_file = fullfile(template_dir, sprintf( ...
        'Template_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
        cfg.target_blades, sensor_tag, template_suffix));
end
dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows%s_20241106.mat', ...
    cfg.target_blades, sensor_tag, dynamic_suffix));

fprintf('\n=== Step02: sliding-window dynamic waveform maps ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('High-speed raw-data folder: %s\n', P.path.dynamicDataDir);
if selection_info.use_region_plan
    fprintf(['Region plan row %d (%s mode): blade %d, selected start %.6f s, ' ...
        'region [%.6f, %.6f] s, peak %.6f s, EO %.0f, dominant freq %.3f Hz.\n'], ...
        selection_info.region_id, selection_info.start_mode, ...
        selection_info.target_blade, ...
        selection_info.start_time_sec, selection_info.region_start_sec, selection_info.end_time_sec, ...
        selection_info.peak_time_sec, selection_info.dominant_order, ...
        selection_info.dominant_freq_hz);
else
    fprintf('Region plan not used. Using explicit target blade/time settings at the top of this file.\n');
end
if isempty(old_result_file)
    fprintf('Old SG comparison: skipped\n');
else
    fprintf('Old SG comparison file: %s\n', old_result_file);
end

if ~isfile(template_file)
    explain_missing_template_local(template_file, route_dir, template_suffix);
end

loaded_template = load(template_file, 'Template');
Template = loaded_template.Template;
hasOldResult = ~isempty(old_result_file) && isfile(old_result_file);
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
DynamicMap.SelectionInfo = selection_info;
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
tiledlayout(fig, numel(cfg.analysis_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = DynamicMap.Window(best_plot_idx).Sensor(is);

    nexttile;
    plot(D.t, D.V, '-', 'LineWidth', 1.0, 'Color', [0 0.45 0.74]);
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d raw waveform, window %d, time [%.6f, %.6f] s', ...
        sid, best_plot_idx, min(D.t), max(D.t)));
    box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');

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
    exportgraphics(fig, fullfile(figure_dir, sprintf('Step02_DynamicMap_B%d_%s_SlidingWindows%s.png', ...
        cfg.target_blades, sensor_tag, dynamic_suffix)), 'Resolution', 300);
end

%% Local helpers kept at the end so the main workflow above stays readable
function info = load_region_selection_local(region_plan_file, step02_cfg)
info = struct( ...
    'use_region_plan', false, ...
    'region_id', NaN, ...
    'start_mode', 'manual', ...
    'target_blade', step02_cfg.targetBlade, ...
    'region_selected_blade', NaN, ...
    'start_time_sec', step02_cfg.analysisStartTime, ...
    'region_start_sec', NaN, ...
    'end_time_sec', NaN, ...
    'peak_time_sec', NaN, ...
    'dominant_order', NaN, ...
    'dominant_freq_hz', NaN, ...
    'region_plan_file', region_plan_file);
if ~isfield(step02_cfg, 'useRegionPlan') || ~step02_cfg.useRegionPlan
    return;
end
if exist(region_plan_file, 'file') ~= 2
    warning('Region plan file not found. Falling back to explicit Step02 settings: %s', region_plan_file);
    return;
end

T = readtable(region_plan_file);
row = T(T.regionId == step02_cfg.regionId, :);
if height(row) ~= 1
    warning('Region plan file must contain exactly one row for regionId=%d. Falling back to explicit Step02 settings.', step02_cfg.regionId);
    return;
end

info.use_region_plan = true;
info.region_id = step02_cfg.regionId;
info.start_mode = lower(strtrim(step02_cfg.regionStartMode));
info.region_selected_blade = row.selectedBladeSlot(1);
if isfield(step02_cfg, 'useRegionBlade') && step02_cfg.useRegionBlade
    info.target_blade = info.region_selected_blade;
else
    info.target_blade = step02_cfg.targetBlade;
end
info.region_start_sec = row.bttStartSec(1);
info.end_time_sec = row.bttEndSec(1);
info.peak_time_sec = row.peakBttTimeSec(1);
info.dominant_order = row.dominantOrder(1);
info.dominant_freq_hz = row.dominantFreqHz(1);
if strcmpi(info.start_mode, 'peak')
    info.start_time_sec = info.peak_time_sec;
elseif strcmpi(info.start_mode, 'start')
    info.start_time_sec = info.region_start_sec;
else
    error('P.step02.regionStartMode must be "peak" or "start".');
end
end

function explain_missing_template_local(template_file, route_dir, template_suffix)
template_dir = fullfile(route_dir, 'output', 'templates');
nearby = dir(fullfile(template_dir, 'Template_LowSpeedRotating_B*_S*_20241106.mat'));
nearby_names = {nearby.name};
if isempty(nearby_names)
    nearby_text = '  (no 20241106 low-speed templates found)';
else
    nearby_text = sprintf('  %s\n', nearby_names{:});
end
error(['OPRCenterStd template file not found:\n  %s\n\n' ...
    'This Step02 route must use the template generated by:\n' ...
    '  Step01_Main_Build_OPRCenterStd_Template_20241106\n\n' ...
    'Requested suffix:\n  %s\n\n' ...
    'Existing 20241106 templates in output/templates:\n%s'], ...
    template_file, template_suffix, nearby_text);
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

