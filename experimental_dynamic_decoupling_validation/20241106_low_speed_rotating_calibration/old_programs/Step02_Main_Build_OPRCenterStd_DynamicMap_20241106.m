%% Step02_Main_Build_OPRCenterStd_DynamicMap_20241106
% Build the latest dynamic sliding-window map in the OPRCenterStd coordinate.
% Tune parameters here, then run this file directly.

clear; close all; clc;

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

P.path.dynamicDataDir = fullfile(baseCfg.dataset_root, dynamicCase);
P.path.caseOutputDir = fullfile(baseCfg.output_root, dynamicCase);
P.path.sensorConfigFile = fullfile(baseCfg.reference_output_dir, 'Sensor_Config_20241106.mat');
P.path.regionPlanFile = fullfile(routeDir, 'outputs', 'Step05_BTT_WindowPlan_20241106.csv');

P.step02.targetBlade = 6;
P.step02.analysisSensors = [2 3 5 7];
P.step02.analysisStartTime = 75.0;
P.step02.targetLaps = 20;
P.step02.windowLaps = 3;
P.step02.slidingStepLaps = 1;
P.step02.checkPlotWindowID = 1;
P.step02.crossSensorMatchMode = 'poly_peak_vector';
P.step02.matchPolyDegree = 6;
P.step02.matchFeatureWeights = struct('peak_v', 0.7, 'peak_x', 0.3);
P.step02.numberingSeedLaps = 1;
P.step02.buildAllBladeDynamicMaps = true;
P.step02.useRegionPlan = false;
P.step02.useRegionBlade = false;
P.step02.regionId = 2;
P.step02.regionStartMode = 'manual';  % 'manual' uses analysisStartTime. Region plan supports 'peak' or 'start'.
P.step02.oldSgResultFile = '';   % Optional, only used for historical comparison.
P.step02.oprChannel = baseCfg.opr_id;
P.step02.oprPulsesPerRev = baseCfg.opr_pulses_per_rev;
P.step02.sampleRateHz = baseCfg.pinlv;
P.step02.tipRadiusMM = baseCfg.r_tip_mm;
P.step02.pulseWindowSec = syncCfg.pulse_window_sec;
P.step02.pulsePadSec = syncCfg.pulse_pad_sec;
P.step02.dynamicWindowMode = syncCfg.dynamic_window_mode;
P.view.saveFigures = true;
runLabelTimeMs = round(1000 * P.step02.analysisStartTime);
if mod(runLabelTimeMs, 1000) == 0
    P.name.runLabel = sprintf('T%03ds', runLabelTimeMs / 1000);
else
    P.name.runLabel = sprintf('T%03dp%03ds', floor(runLabelTimeMs / 1000), mod(runLabelTimeMs, 1000));
end
P.name.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', P.name.runLabel, P.name.templateSuffix);

envOverrideText = lower(strtrim(getenv('STEP20241106_USE_ENV_OVERRIDES')));
if ismember(envOverrideText, {'1', 'true', 'yes', 'y', 'on'})
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
        runLabelTimeMs = round(1000 * P.step02.analysisStartTime);
        if mod(runLabelTimeMs, 1000) == 0
            P.name.runLabel = sprintf('T%03ds', runLabelTimeMs / 1000);
        else
            P.name.runLabel = sprintf('T%03dp%03ds', floor(runLabelTimeMs / 1000), mod(runLabelTimeMs, 1000));
        end
        P.name.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', P.name.runLabel, P.name.templateSuffix);
    end
    dynamicSuffixEnv = strtrim(getenv('STEP20241106_DYNAMIC_SUFFIX'));
    if ~isempty(dynamicSuffixEnv)
        P.name.dynamicSuffix = dynamicSuffixEnv;
    end
    numberingSeedLapsEnv = str2double(strtrim(getenv('STEP20241106_NUMBERING_SEED_LAPS')));
    if isfinite(numberingSeedLapsEnv) && numberingSeedLapsEnv >= 1
        P.step02.numberingSeedLaps = round(numberingSeedLapsEnv);
    end
    buildAllBladeMapsEnv = strtrim(getenv('STEP20241106_BUILD_ALL_BLADE_MAPS'));
    if ~isempty(buildAllBladeMapsEnv)
        P.step02.buildAllBladeDynamicMaps = ismember(lower(strtrim(buildAllBladeMapsEnv)), ...
            {'1', 'true', 'yes', 'y', 'on'});
    end
end

%% Files
sensorTag = ['S', sprintf('%d', P.step02.analysisSensors)];
defaultRegionSelectionFile = fullfile(routeDir, 'output', 'region_selection', ...
    sprintf('HighSpeedRegionSelection_%s_20241106.mat', P.name.runLabel));
regionSelectionFileEnv = strtrim(getenv('STEP20241106_REGION_SELECTION_FILE'));
if isempty(regionSelectionFileEnv)
    P.path.regionSelectionFile = defaultRegionSelectionFile;
else
    P.path.regionSelectionFile = regionSelectionFileEnv;
end
templateFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
    P.step02.targetBlade, sensorTag, P.name.templateSuffix));

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
dynamic_suffix = sanitize_dynamic_suffix_local(P.name.dynamicSuffix);
selection_info = load_region_selection_local(P.path.regionSelectionFile, P.path.regionPlanFile, P.step02);
if selection_info.use_saved_selection || selection_info.use_region_plan
    cfg.target_blades = selection_info.target_blade;
    cfg.analysis_start_time = selection_info.start_time_sec;
    if ~isempty(selection_info.analysis_sensors)
        cfg.analysis_sensors = unique(selection_info.analysis_sensors, 'stable');
    end
    if isfinite(selection_info.target_laps) && selection_info.target_laps >= 1
        cfg.target_laps = round(selection_info.target_laps);
    end
end
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
source_result_sensor_tag = sensor_tag;

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
numbering_dir = fullfile(output_dir, 'blade_numbering');
figure_dir = fullfile(output_dir, 'figures');
if exist(dynamic_dir, 'dir') ~= 7; mkdir(dynamic_dir); end
if exist(numbering_dir, 'dir') ~= 7; mkdir(numbering_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_suffix = P.name.templateSuffix;
if isempty(template_suffix)
    template_file = fullfile(template_dir, sprintf( ...
        'TemplateBundle_LowSpeedRotating_B%d_%s_20241106.mat', cfg.target_blades, sensor_tag));
else
    template_file = fullfile(template_dir, sprintf( ...
        'TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
        cfg.target_blades, sensor_tag, template_suffix));
end
dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows%s_20241106.mat', ...
    cfg.target_blades, sensor_tag, dynamic_suffix));
numbering_file = fullfile(numbering_dir, sprintf('BladeNumbering_%s%s_20241106.mat', ...
    sensor_tag, dynamic_suffix));

fprintf('\n=== Step02: sliding-window dynamic waveform maps ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('High-speed raw-data folder: %s\n', P.path.dynamicDataDir);
if selection_info.use_saved_selection
    fprintf(['Saved region selection (%s): blade %d, start %.6f s, ' ...
        'region [%.6f, %.6f] s, peak %.6f s, EO %.0f, dominant freq %.3f Hz.\n'], ...
        selection_info.selection_source, selection_info.target_blade, ...
        selection_info.start_time_sec, selection_info.region_start_sec, selection_info.end_time_sec, ...
        selection_info.peak_time_sec, selection_info.dominant_order, ...
        selection_info.dominant_freq_hz);
    fprintf('Region selection file: %s\n', selection_info.region_selection_file);
elseif selection_info.use_region_plan
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
fprintf('Dynamic map output: %s\n', dynamic_file);
fprintf('Blade numbering output: %s\n', numbering_file);
if isempty(old_result_file)
    fprintf('Old SG comparison: skipped\n');
else
    fprintf('Old SG comparison file: %s\n', old_result_file);
end

if ~isfile(template_file)
    explain_missing_template_local(template_file, route_dir, template_suffix);
end

Template = load_template_for_sensors_local(template_file, cfg.analysis_sensors);
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
cfg.cross_sensor_match_mode = lower(strtrim(P.step02.crossSensorMatchMode));
cfg.match_poly_degree = P.step02.matchPolyDegree;
cfg.match_feature_weights = P.step02.matchFeatureWeights;
cfg.numbering_seed_laps = P.step02.numberingSeedLaps;
cfg.build_all_blade_dynamic_maps = logical(P.step02.buildAllBladeDynamicMaps);

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
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    cfg.blades_num = size(Sensor_Config.Standard_Relative_Angles_OPRCenter, 2);
else
    cfg.blades_num = size(Sensor_Config.Standard_Relative_Angles, 2);
end
loaded_opr = load(fullfile(cfg.case_output_dir, 'jiluOPR.mat'), 'jiluOPR');
opr_times = loaded_opr.jiluOPR(:, 1);
opr_reference = build_opr_reference_from_jilu_local(loaded_opr.jiluOPR, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
F_omega_deg = build_phase_speed_local(opr_times, cfg.opr_pulses_per_rev);
fprintf('Dynamic OPR center shift: %.4f deg, %.4f mm from rising edge.\n', ...
    opr_reference.phase_shift_deg, opr_reference.phase_shift_mm);

%% Select the same 20 target-blade laps once on an anchor sensor, then project to all sensors
probe_cache = load_probe_cache_local(cfg.case_output_dir, cfg.analysis_sensors, opr_times, cfg.opr_pulses_per_rev);
[selection, selection_diag, global_window] = select_target_rows_by_anchor_projection_local( ...
    probe_cache, cfg, Sensor_Config, opr_reference);
selection_info.row_selection_mode = selection_diag.mode;
selection_info.anchor_sensor_id = selection_diag.anchor_sensor_id;
selection_info.anchor_selected_rows = selection_diag.anchor_selected_rows;
selection_info.sensor_projection = selection_diag.sensor_projection;
fprintf('Row selection mode: %s, anchor sensor: CH%d\n', ...
    selection_diag.mode, selection_diag.anchor_sensor_id);
for is = 1:numel(selection_diag.sensor_projection)
    D = selection_diag.sensor_projection(is);
    fprintf('  CH%d projected rows: %d laps, mean|max angle error = %.3f|%.3f deg\n', ...
        D.sensor_id, numel(D.selected_rows), D.mean_abs_angle_error_deg, D.max_abs_angle_error_deg);
end
fprintf('Global raw-data window: %.6f-%.6f s\n', global_window(1), global_window(2));
load_window = global_window;
if strcmpi(cfg.cross_sensor_match_mode, 'poly_peak_vector')
    load_window = build_poly_match_candidate_window_local(probe_cache, selection_diag.anchor_sensor_id, ...
        selection_diag.anchor_selected_rows, cfg);
    fprintf('Poly-peak matching candidate raw-data window: %.6f-%.6f s\n', load_window(1), load_window(2));
end

%% Load raw dynamic voltage streams only once
file_ranges = build_dynamic_file_ranges_local(cfg.dynamic_data_dir, cfg.opr_channel, cfg.pinlv);
selected_file_mask = [file_ranges.t_end] >= load_window(1) & [file_ranges.t_start] <= load_window(2);
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
        keep = t_global >= load_window(1) & t_global <= load_window(2);
        raw_stream(sid).T = [raw_stream(sid).T; t_global(keep)];
        raw_stream(sid).V = [raw_stream(sid).V; v_local(keep)];
    end
end

if strcmpi(cfg.cross_sensor_match_mode, 'poly_peak_vector')
    [selection, selection_diag, global_window] = refine_selection_by_poly_peak_vector_local( ...
        selection, selection_diag, probe_cache, raw_stream, cfg, Sensor_Config, opr_reference, F_omega_deg);
    selection_info.row_selection_mode = selection_diag.mode;
    selection_info.sensor_projection = selection_diag.sensor_projection;
    selection_info.blade_numbering = build_blade_numbering_info_local(selection_diag, cfg);
    selection_info.blade_numbering_file = numbering_file;
    fprintf('Row selection refined by six-peak fingerprint locking (seed laps = %d).\n', cfg.numbering_seed_laps);
    for is = 1:numel(selection_diag.sensor_projection)
        D = selection_diag.sensor_projection(is);
        fprintf(['  CH%d physical->local %s, target B%d->local B%d, group score %.3f, ' ...
            'mean|max angle error = %.3f|%.3f deg\n'], ...
            D.sensor_id, mat2str(D.physical_to_local_blade_ids), ...
            cfg.target_blades, D.matched_blade_id, D.match_score, ...
            D.mean_abs_angle_error_deg, D.max_abs_angle_error_deg);
    end
    fprintf('Refined global raw-data window: %.6f-%.6f s\n', global_window(1), global_window(2));
end

%% Extract each selected lap waveform once
LapData = repmat(struct('sensor_id', NaN, 'Lap', []), numel(cfg.analysis_sensors), 1);
for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_file = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;
    theta_blade_id = selection(is).theta_blade_id;
    if ~isfinite(theta_blade_id)
        theta_blade_id = cfg.target_blades;
    end
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, theta_blade_id, opr_reference);

    Lap = repmat(struct('lap_id', NaN, 't', [], 'x_abs', [], 'V', [], 'theta', []), cfg.target_laps, 1);
    for lap_id = 1:cfg.target_laps
        row_id = selection(is).selected_rows(lap_id);
        t_peak = jilublade(row_id, 3);
        [t_start, t_end] = build_dynamic_segment_window_local(jilublade, row_id, t_peak, cfg);

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
if isfield(selection_info, 'blade_numbering')
    BladeNumbering = selection_info.blade_numbering; %#ok<NASGU>
    save(numbering_file, 'BladeNumbering', '-v7.3');
end
fprintf('Saved sliding-window dynamic map: %s\n', dynamic_file);
if isfield(selection_info, 'blade_numbering')
    fprintf('Saved blade numbering map: %s\n', numbering_file);
end
fprintf('Built %d windows.\n', num_windows);
if cfg.build_all_blade_dynamic_maps && isfield(selection_info, 'blade_numbering')
    [all_map_files, skipped_blades] = save_all_physical_blade_dynamic_maps_local( ...
        DynamicMap, selection_info, selection_diag, raw_stream, probe_cache, opr_times, ...
        F_omega_deg, Sensor_Config, opr_reference, cfg, template_dir, template_suffix, ...
        dynamic_dir, sensor_tag, dynamic_suffix, selected_file_ranges);
    fprintf('All-physical-blade DynamicMaps saved: %d file(s).\n', numel(all_map_files));
    if ~isempty(skipped_blades)
        fprintf('All-physical-blade DynamicMaps skipped for blades: %s\n', mat2str(skipped_blades));
    end
end
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

function [saved_files, skipped_blades] = save_all_physical_blade_dynamic_maps_local( ...
        base_map, selection_info, selection_diag, raw_stream, probe_cache, opr_times, ...
        F_omega_deg, Sensor_Config, opr_reference, cfg, template_dir, template_suffix, ...
        dynamic_dir, sensor_tag, dynamic_suffix, selected_file_ranges)
saved_files = {};
skipped_blades = [];

for blade_id = 1:cfg.blades_num
    template_file = build_template_file_name_local(template_dir, blade_id, sensor_tag, template_suffix);
    if exist(template_file, 'file') ~= 2
        warning('Skipping B%d DynamicMap export because low-speed template is missing: %s', ...
            blade_id, template_file);
        skipped_blades(end + 1) = blade_id; %#ok<AGROW>
        continue;
    end

    try
        TemplateBlade = load_template_for_sensors_local(template_file, cfg.analysis_sensors);
        selection_blade = build_physical_blade_selection_local(selection_diag, cfg, blade_id);
        [WindowBlade, global_window_blade] = build_dynamic_windows_for_physical_blade_local( ...
            selection_blade, raw_stream, probe_cache, opr_times, F_omega_deg, ...
            Sensor_Config, opr_reference, TemplateBlade, cfg, blade_id);
    catch ME
        warning('Skipping B%d DynamicMap export: %s', blade_id, ME.message);
        skipped_blades(end + 1) = blade_id; %#ok<AGROW>
        continue;
    end

    DynamicMap = base_map; %#ok<NASGU>
    DynamicMap.TargetBlade = blade_id;
    DynamicMap.TemplateFileForXRel = template_file;
    DynamicMap.XCenterBySensor = build_xcenter_table_local(TemplateBlade, cfg.analysis_sensors);
    DynamicMap.GlobalTimeWindow = global_window_blade;
    DynamicMap.Window = WindowBlade;
    DynamicMap.SelectedRawFileIDs = [selected_file_ranges.file_id];
    DynamicMap.SourceSettings.target_blades = blade_id;
    DynamicMap.SelectionInfo = selection_info;
    DynamicMap.SelectionInfo.target_blade = blade_id;
    DynamicMap.SelectionInfo.physical_blade_id = blade_id;
    DynamicMap.SelectionInfo.selected_rows_for_physical_blade = selection_blade;
    DynamicMap.SelectionInfo.workflow_note = ['High-speed numbering was locked once for the sensor subset ' ...
        'and region; this DynamicMap is a physical-blade slice from that shared numbering.'];

    dynamic_file_blade = fullfile(dynamic_dir, sprintf( ...
        'DynamicMap_B%d_%s_SlidingWindows%s_20241106.mat', ...
        blade_id, sensor_tag, dynamic_suffix));
    save(dynamic_file_blade, 'DynamicMap', '-v7.3');
    saved_files{end + 1} = dynamic_file_blade; %#ok<AGROW>
end
end

function template_file = build_template_file_name_local(template_dir, blade_id, sensor_tag, template_suffix)
if isempty(template_suffix)
    template_file = fullfile(template_dir, sprintf( ...
        'TemplateBundle_LowSpeedRotating_B%d_%s_20241106.mat', blade_id, sensor_tag));
else
    template_file = fullfile(template_dir, sprintf( ...
        'TemplateBundle_LowSpeedRotating_B%d_%s_%s_20241106.mat', ...
        blade_id, sensor_tag, template_suffix));
end
end

function selection_blade = build_physical_blade_selection_local(selection_diag, cfg, blade_id)
selection_blade = repmat(struct( ...
    'sensor_id', NaN, ...
    'selected_rows', [], ...
    'matched_blade_id', NaN, ...
    'theta_blade_id', blade_id), numel(cfg.analysis_sensors), 1);

for is = 1:numel(cfg.analysis_sensors)
    D = selection_diag.sensor_projection(is);
    if isfield(D, 'selected_rows_by_physical') && ...
            size(D.selected_rows_by_physical, 2) >= blade_id
        selected_rows = D.selected_rows_by_physical(:, blade_id);
    elseif blade_id == cfg.target_blades
        selected_rows = D.selected_rows(:);
    else
        selected_rows = nan(cfg.target_laps, 1);
    end
    selected_rows = enforce_strictly_increasing_rows_local(selected_rows, D.sensor_id);
    selection_blade(is).sensor_id = D.sensor_id;
    selection_blade(is).selected_rows = selected_rows(:);
    if isfield(D, 'physical_to_local_blade_ids') && numel(D.physical_to_local_blade_ids) >= blade_id
        selection_blade(is).matched_blade_id = D.physical_to_local_blade_ids(blade_id);
    else
        selection_blade(is).matched_blade_id = blade_id;
    end
    selection_blade(is).theta_blade_id = blade_id;
end
end

function [Window, global_window] = build_dynamic_windows_for_physical_blade_local( ...
        selection_blade, raw_stream, probe_cache, opr_times, F_omega_deg, ...
        Sensor_Config, opr_reference, Template, cfg, blade_id)
LapData = repmat(struct('sensor_id', NaN, 'Lap', []), numel(cfg.analysis_sensors), 1);
global_window = [inf, -inf];

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_idx = find([probe_cache.sensor_id] == sid, 1, 'first');
    if isempty(probe_idx)
        error('Probe cache missing CH%d.', sid);
    end
    probe = probe_cache(probe_idx);
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference);
    if ~isfinite(theta_std)
        error('Standard relative angle is missing for CH%d blade %d.', sid, blade_id);
    end

    Lap = repmat(struct('lap_id', NaN, 't', [], 'x_abs', [], 'V', [], 'theta', []), cfg.target_laps, 1);
    selected_rows = selection_blade(is).selected_rows(:);
    for lap_id = 1:cfg.target_laps
        row_id = selected_rows(lap_id);
        t_peak = probe.jilublade(row_id, 3);
        [t_start, t_end] = build_dynamic_segment_window_local(probe.jilublade, row_id, t_peak, cfg);

        mask = raw_stream(sid).T >= t_start & raw_stream(sid).T <= t_end;
        t_seg = raw_stream(sid).T(mask);
        v_seg = raw_stream(sid).V(mask);
        idx_prev = find(opr_times < t_peak, 1, 'last');
        if isempty(idx_prev) || numel(t_seg) < 5
            continue;
        end

        theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
        theta_diff_deg = wrap_to_signed_period_local(theta_points_deg - theta_std, 360);
        x_abs = theta_diff_deg * (pi / 180) * cfg.r_tip_mm;
        theta_rot = map_time_to_rotor_phase_local(opr_times, t_seg, cfg.opr_pulses_per_rev);
        valid = isfinite(theta_rot);

        Lap(lap_id).lap_id = lap_id;
        Lap(lap_id).t = t_seg(valid);
        Lap(lap_id).x_abs = x_abs(valid);
        Lap(lap_id).V = v_seg(valid);
        Lap(lap_id).theta = theta_rot(valid);
        if any(valid)
            global_window(1) = min(global_window(1), min(t_seg(valid)));
            global_window(2) = max(global_window(2), max(t_seg(valid)));
        end
    end

    LapData(is).sensor_id = sid;
    LapData(is).Lap = Lap;
end

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

    if isempty(all_t_window)
        error('No dynamic samples were extracted for B%d window %d.', blade_id, w_idx);
    end
    Window(w_idx).window_id = w_idx;
    Window(w_idx).lap_range = lap_range;
    Window(w_idx).time_window = [min(all_t_window), max(all_t_window)];
    Window(w_idx).rot_freq_mean_hz = compute_local_rot_freq_local(opr_times, cfg.opr_pulses_per_rev, Window(w_idx).time_window);
    Window(w_idx).rot_rpm_mean = 60 * Window(w_idx).rot_freq_mean_hz;
    Window(w_idx).Sensor = Sensor;
end

if ~all(isfinite(global_window))
    error('No finite global time window was built for B%d.', blade_id);
end
end

%% Local helpers kept at the end so the main workflow above stays readable
function info = load_region_selection_local(region_selection_file, region_plan_file, step02_cfg)
info = struct( ...
    'use_saved_selection', false, ...
    'use_region_plan', false, ...
    'selection_source', 'explicit_step02', ...
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
    'target_laps', step02_cfg.targetLaps, ...
    'analysis_sensors', step02_cfg.analysisSensors, ...
    'region_selection_file', region_selection_file, ...
    'region_plan_file', region_plan_file);
if exist(region_selection_file, 'file') == 2
    loaded = load(region_selection_file);
    if isfield(loaded, 'RegionSelection') && isstruct(loaded.RegionSelection)
        info = merge_region_selection_local(info, loaded.RegionSelection);
        info.use_saved_selection = true;
        return;
    end
    warning('Saved region selection file exists but does not contain struct RegionSelection: %s', region_selection_file);
end
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

function info = merge_region_selection_local(info, source)
fields = fieldnames(info);
for i = 1:numel(fields)
    name = fields{i};
    if isfield(source, name) && ~isempty(source.(name))
        info.(name) = source.(name);
    end
end
if isfield(source, 'selection_source')
    info.selection_source = source.selection_source;
end
if isfield(source, 'target_laps') && isfinite(source.target_laps)
    info.target_laps = source.target_laps;
end
if isfield(source, 'analysis_sensors') && ~isempty(source.analysis_sensors)
    info.analysis_sensors = source.analysis_sensors;
end
if ~isfinite(info.start_time_sec) && isfinite(info.region_start_sec)
    info.start_time_sec = info.region_start_sec;
end
end

function probe_cache = load_probe_cache_local(case_output_dir, sensor_ids, opr_times, opr_pulses_per_rev)
probe_cache = repmat(struct( ...
    'sensor_id', NaN, ...
    'jilublade', [], ...
    'rel_angles_deg', [], ...
    'prev_opr_index', [], ...
    'revolution_index', [], ...
    'opr_times', []), numel(sensor_ids), 1);
for is = 1:numel(sensor_ids)
    sid = sensor_ids(is);
    probe_file = fullfile(case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    if exist(probe_file, 'file') ~= 2
        error('Missing Step2 pulse file for CH%d: %s', sid, probe_file);
    end
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;
    [rel_angles_deg, prev_opr_index, revolution_index] = ...
        compute_pulse_reference_geometry_local(jilublade(:, 3), opr_times, opr_pulses_per_rev);
    probe_cache(is).sensor_id = sid;
    probe_cache(is).jilublade = jilublade;
    probe_cache(is).rel_angles_deg = rel_angles_deg;
    probe_cache(is).prev_opr_index = prev_opr_index;
    probe_cache(is).revolution_index = revolution_index;
    probe_cache(is).opr_times = opr_times;
end
end

function [selection, diag, global_window] = select_target_rows_by_anchor_projection_local( ...
        probe_cache, cfg, Sensor_Config, opr_reference)
selection = repmat(struct( ...
    'sensor_id', NaN, ...
    'selected_rows', [], ...
    'selected_rows_by_physical', [], ...
    'angle_error_by_physical', [], ...
    'matched_blade_id', NaN, ...
    'theta_blade_id', NaN), numel(cfg.analysis_sensors), 1);
sensor_projection = repmat(struct( ...
    'sensor_id', NaN, ...
    'selected_rows', [], ...
    'matched_blade_id', NaN, ...
    'theta_blade_id', NaN, ...
    'physical_to_local_blade_ids', [], ...
    'local_to_physical_blade_ids', [], ...
    'blade_index_offset', NaN, ...
    'blade_index_offset_scores', [], ...
    'mean_abs_angle_error_deg', NaN, ...
    'max_abs_angle_error_deg', NaN, ...
    'match_score', NaN, ...
    'corr_peak_v', NaN, ...
    'corr_peak_x', NaN, ...
    'match_margin', NaN, ...
    'fit_valid_count', NaN), numel(cfg.analysis_sensors), 1);
global_window = [inf, -inf];

anchor_sid = choose_anchor_sensor_local(probe_cache, cfg, Sensor_Config);
anchor_idx = find([probe_cache.sensor_id] == anchor_sid, 1, 'first');
anchor_probe = probe_cache(anchor_idx);
anchor_mask = anchor_probe.jilublade(:, 4) == cfg.target_blades & ...
    anchor_probe.jilublade(:, 3) >= cfg.analysis_start_time;
anchor_rows = find(anchor_mask);
if numel(anchor_rows) < cfg.target_laps
    error('Anchor CH%d has only %d target-blade laps after %.3f s.', ...
        anchor_sid, numel(anchor_rows), cfg.analysis_start_time);
end
anchor_rows = anchor_rows(1:cfg.target_laps);

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    cache_idx = find([probe_cache.sensor_id] == sid, 1, 'first');
    probe = probe_cache(cache_idx);
    selected_rows = nan(cfg.target_laps, 1);
    angle_error_deg = nan(cfg.target_laps, 1);

    if sid == anchor_sid
        selected_rows = anchor_rows(:);
        angle_error_deg(:) = 0;
        matched_blade_id = cfg.target_blades;
        match_score = 1.0;
    else
        theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, cfg.target_blades, opr_reference);
        if ~isfinite(theta_std)
            error('Standard relative angle is missing for CH%d blade %d.', sid, cfg.target_blades);
        end
        for k = 1:cfg.target_laps
            rev_id = anchor_probe.revolution_index(anchor_rows(k));
            if ~isfinite(rev_id)
                error('Anchor CH%d row %d does not have a valid revolution index.', anchor_sid, anchor_rows(k));
            end
            candidates = find(probe.revolution_index == rev_id & isfinite(probe.rel_angles_deg));
            if isempty(candidates)
                error('CH%d does not have any pulse candidate in revolution %d.', sid, rev_id);
            end
            candidate_error = wrap_to_signed_period_local(probe.rel_angles_deg(candidates) - theta_std, 360);
            [~, best_pos] = min(abs(candidate_error));
            selected_rows(k) = candidates(best_pos);
            angle_error_deg(k) = candidate_error(best_pos);
        end
        matched_blade_id = cfg.target_blades;
        match_score = nan;
    end

    selected_rows = enforce_strictly_increasing_rows_local(selected_rows, sid);
    selection(is).sensor_id = sid;
    selection(is).selected_rows = selected_rows(:);
    selection(is).matched_blade_id = matched_blade_id;
    selection(is).theta_blade_id = matched_blade_id;
    sensor_projection(is).sensor_id = sid;
    sensor_projection(is).selected_rows = selected_rows(:);
    rows_by_physical = nan(numel(selected_rows), cfg.blades_num);
    errors_by_physical = nan(numel(selected_rows), cfg.blades_num);
    rows_by_physical(:, cfg.target_blades) = selected_rows(:);
    errors_by_physical(:, cfg.target_blades) = angle_error_deg(:);
    sensor_projection(is).selected_rows_by_physical = rows_by_physical;
    sensor_projection(is).angle_error_by_physical = errors_by_physical;
    sensor_projection(is).matched_blade_id = matched_blade_id;
    sensor_projection(is).theta_blade_id = matched_blade_id;
    sensor_projection(is).mean_abs_angle_error_deg = mean(abs(angle_error_deg), 'omitnan');
    sensor_projection(is).max_abs_angle_error_deg = max(abs(angle_error_deg), [], 'omitnan');
    sensor_projection(is).match_score = match_score;

    jilublade = probe.jilublade;
    for k = 1:numel(selected_rows)
        [t_start, t_end] = build_dynamic_segment_window_local(jilublade, selected_rows(k), jilublade(selected_rows(k), 3), cfg);
        global_window(1) = min(global_window(1), t_start);
        global_window(2) = max(global_window(2), t_end);
    end
end

diag = struct();
diag.mode = 'anchor_projection_opr_angle';
diag.anchor_sensor_id = anchor_sid;
diag.anchor_selected_rows = anchor_rows(:);
diag.sensor_projection = sensor_projection;
end

function load_window = build_poly_match_candidate_window_local(probe_cache, anchor_sensor_id, anchor_selected_rows, cfg)
anchor_idx = find([probe_cache.sensor_id] == anchor_sensor_id, 1, 'first');
anchor_probe = probe_cache(anchor_idx);
selected_rev_ids = anchor_probe.revolution_index(anchor_selected_rows(:));
selected_rev_ids = selected_rev_ids(isfinite(selected_rev_ids));
load_window = [inf, -inf];

for is = 1:numel(probe_cache)
    probe = probe_cache(is);
    for k = 1:numel(selected_rev_ids)
        row_ids = find(probe.revolution_index == selected_rev_ids(k));
        for ir = 1:numel(row_ids)
            row_id = row_ids(ir);
            [t_start, t_end] = build_dynamic_segment_window_local(probe.jilublade, row_id, probe.jilublade(row_id, 3), cfg);
            load_window(1) = min(load_window(1), t_start);
            load_window(2) = max(load_window(2), t_end);
        end
    end
end

if ~all(isfinite(load_window)) || load_window(2) <= load_window(1)
    error('Failed to build the candidate raw-data window for polynomial peak-vector matching.');
end
end

function [selection, diag, global_window] = refine_selection_by_poly_peak_vector_local( ...
        selection, diag, probe_cache, raw_stream, cfg, Sensor_Config, opr_reference, F_omega_deg)
anchor_sid = diag.anchor_sensor_id;
anchor_rows = diag.anchor_selected_rows(:);
anchor_idx = find([probe_cache.sensor_id] == anchor_sid, 1, 'first');
anchor_probe = probe_cache(anchor_idx);
selected_rev_ids = anchor_probe.revolution_index(anchor_rows);
if any(~isfinite(selected_rev_ids))
    error('Anchor rows contain invalid revolution indices for polynomial peak-vector matching.');
end
global_window = [inf, -inf];

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_idx = find([probe_cache.sensor_id] == sid, 1, 'first');
    probe = probe_cache(probe_idx);
    match_result = lock_blade_numbering_by_fingerprint_local( ...
        raw_stream(sid), probe, selected_rev_ids, cfg, Sensor_Config);
    if ~match_result.locked
        warning('CH%d fingerprint locking failed. Falling back to initial target-blade projection.', sid);
        theta_std_fallback = read_opr_center_standard_angle_local(Sensor_Config, sid, cfg.target_blades, opr_reference);
        [selected_rows, angle_error_deg] = project_rows_by_revolution_and_angle_local(probe, selected_rev_ids, theta_std_fallback);
        matched_blade_id = cfg.target_blades;
        theta_blade_id = cfg.target_blades;
        best_score = nan;
        best_corr_peak_v = nan;
        best_corr_peak_x = nan;
        score_margin = nan;
        fit_valid_count = 0;
        seed_lap_count = 0;
        physical_to_local_blade_ids = nan(1, cfg.blades_num);
        local_to_physical_blade_ids = nan(1, cfg.blades_num);
        blade_index_offset = nan;
        blade_index_offset_scores = nan(1, cfg.blades_num);
        selected_rows_by_physical_matrix = nan(numel(selected_rows), cfg.blades_num);
        angle_error_by_physical_matrix = nan(numel(selected_rows), cfg.blades_num);
        selected_rows_by_physical_matrix(:, cfg.target_blades) = selected_rows(:);
        angle_error_by_physical_matrix(:, cfg.target_blades) = angle_error_deg(:);
    else
        matched_blade_id = match_result.physical_to_local_blade_ids(cfg.target_blades);
        theta_blade_id = cfg.target_blades;
        selected_rows = match_result.selected_rows_by_physical{cfg.target_blades};
        angle_error_deg = match_result.angle_error_by_physical{cfg.target_blades};
        best_score = match_result.best_score;
        best_corr_peak_v = match_result.target_corr_peak_v;
        best_corr_peak_x = nan;
        score_margin = match_result.score_margin;
        fit_valid_count = match_result.valid_lap_count;
        seed_lap_count = match_result.seed_lap_count;
        physical_to_local_blade_ids = match_result.physical_to_local_blade_ids;
        local_to_physical_blade_ids = match_result.local_to_physical_blade_ids;
        blade_index_offset = match_result.best_shift;
        blade_index_offset_scores = match_result.shift_scores;
        selected_rows_by_physical_matrix = cell_rows_to_matrix_local( ...
            match_result.selected_rows_by_physical, numel(selected_rev_ids), cfg.blades_num);
        angle_error_by_physical_matrix = cell_rows_to_matrix_local( ...
            match_result.angle_error_by_physical, numel(selected_rev_ids), cfg.blades_num);
    end

    selected_rows = enforce_strictly_increasing_rows_local(selected_rows, sid);
    selection(is).selected_rows = selected_rows(:);
    selection(is).matched_blade_id = matched_blade_id;
    selection(is).theta_blade_id = theta_blade_id;
    diag.sensor_projection(is).selected_rows = selected_rows(:);
    diag.sensor_projection(is).selected_rows_by_physical = selected_rows_by_physical_matrix;
    diag.sensor_projection(is).angle_error_by_physical = angle_error_by_physical_matrix;
    diag.sensor_projection(is).matched_blade_id = matched_blade_id;
    diag.sensor_projection(is).theta_blade_id = theta_blade_id;
    diag.sensor_projection(is).physical_to_local_blade_ids = physical_to_local_blade_ids(:).';
    diag.sensor_projection(is).local_to_physical_blade_ids = local_to_physical_blade_ids(:).';
    diag.sensor_projection(is).blade_index_offset = blade_index_offset;
    diag.sensor_projection(is).blade_index_offset_scores = blade_index_offset_scores(:).';
    diag.sensor_projection(is).mean_abs_angle_error_deg = mean(abs(angle_error_deg), 'omitnan');
    diag.sensor_projection(is).max_abs_angle_error_deg = max(abs(angle_error_deg), [], 'omitnan');
    diag.sensor_projection(is).match_score = best_score;
    diag.sensor_projection(is).corr_peak_v = best_corr_peak_v;
    diag.sensor_projection(is).corr_peak_x = best_corr_peak_x;
    diag.sensor_projection(is).match_margin = score_margin;
    diag.sensor_projection(is).fit_valid_count = fit_valid_count;
    diag.sensor_projection(is).seed_lap_count = seed_lap_count;

    for k = 1:numel(selected_rows)
        row_id = selected_rows(k);
        [t_start, t_end] = build_dynamic_segment_window_local(probe.jilublade, row_id, probe.jilublade(row_id, 3), cfg);
        global_window(1) = min(global_window(1), t_start);
        global_window(2) = max(global_window(2), t_end);
    end
end

diag.mode = 'anchor_projection_then_fingerprint_lock';
end

function M = cell_rows_to_matrix_local(C, row_count, blade_count)
M = nan(row_count, blade_count);
for blade_id = 1:min(numel(C), blade_count)
    values = C{blade_id};
    n = min(numel(values), row_count);
    if n > 0
        M(1:n, blade_id) = values(1:n);
    end
end
end

function [selected_rows, angle_error_deg] = project_rows_by_revolution_and_angle_local(probe, selected_rev_ids, theta_std)
selected_rows = nan(numel(selected_rev_ids), 1);
angle_error_deg = nan(numel(selected_rev_ids), 1);
for k = 1:numel(selected_rev_ids)
    candidates = find(probe.revolution_index == selected_rev_ids(k) & isfinite(probe.rel_angles_deg));
    if isempty(candidates)
        continue;
    end
    candidate_error = wrap_to_signed_period_local(probe.rel_angles_deg(candidates) - theta_std, 360);
    [~, best_pos] = min(abs(candidate_error));
    selected_rows(k) = candidates(best_pos);
    angle_error_deg(k) = candidate_error(best_pos);
end
end

function feature_bank = build_poly_peak_feature_bank_local( ...
        raw_sensor, probe, selected_rev_ids, cfg, Sensor_Config, opr_reference, F_omega_deg)
feature_bank = repmat(struct( ...
    'blade_id', NaN, ...
    'theta_std', NaN, ...
    'selected_rows', [], ...
    'angle_error_deg', [], ...
    'peak_v', [], ...
    'peak_x', [], ...
    'fit_ok', [], ...
    'peak_v_norm', [], ...
    'peak_x_norm', []), cfg.blades_num, 1);
for blade_id = 1:cfg.blades_num
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, probe.sensor_id, blade_id, opr_reference);
    [selected_rows, angle_error_deg] = project_rows_by_revolution_and_angle_local(probe, selected_rev_ids, theta_std);
    feature = build_poly_peak_feature_for_rows_local(raw_sensor, probe, selected_rows, theta_std, cfg, F_omega_deg);
    feature_bank(blade_id).blade_id = blade_id;
    feature_bank(blade_id).theta_std = theta_std;
    feature_bank(blade_id).selected_rows = selected_rows(:);
    feature_bank(blade_id).angle_error_deg = angle_error_deg(:);
    feature_bank(blade_id).peak_v = feature.peak_v;
    feature_bank(blade_id).peak_x = feature.peak_x;
    feature_bank(blade_id).fit_ok = feature.fit_ok;
    feature_bank(blade_id).peak_v_norm = feature.peak_v_norm;
    feature_bank(blade_id).peak_x_norm = feature.peak_x_norm;
end
end

function match_result = lock_blade_numbering_by_fingerprint_local( ...
        raw_sensor, probe, selected_rev_ids, cfg, Sensor_Config)
match_result = struct();
match_result.locked = false;
match_result.best_shift = NaN;
match_result.shift_scores = nan(1, cfg.blades_num);
match_result.best_score = NaN;
match_result.score_margin = NaN;
match_result.valid_lap_count = 0;
match_result.seed_lap_count = 0;
match_result.local_to_physical_blade_ids = nan(1, cfg.blades_num);
match_result.physical_to_local_blade_ids = nan(1, cfg.blades_num);
match_result.target_corr_peak_v = NaN;
match_result.selected_rows_by_physical = repmat({nan(numel(selected_rev_ids), 1)}, cfg.blades_num, 1);
match_result.angle_error_by_physical = repmat({nan(numel(selected_rev_ids), 1)}, cfg.blades_num, 1);

sid = probe.sensor_id;
if ~isfield(Sensor_Config, 'Fingerprints') || ~isKey(Sensor_Config.Fingerprints, sid)
    return;
end
ref_fingerprint = normalize_fingerprint_max_local(Sensor_Config.Fingerprints(sid));
if numel(ref_fingerprint) ~= cfg.blades_num || any(~isfinite(ref_fingerprint))
    return;
end

shift_corr = nan(numel(selected_rev_ids), cfg.blades_num);
valid_lap_mask = false(numel(selected_rev_ids), 1);
candidate_rows_cache = cell(numel(selected_rev_ids), 1);
seed_lap_indices = [];
seed_lap_target = max(1, min(cfg.numbering_seed_laps, numel(selected_rev_ids)));

for k = 1:numel(selected_rev_ids)
    [candidate_rows, peak_vec_norm] = build_revolution_fingerprint_local(raw_sensor, probe, selected_rev_ids(k), cfg);
    if numel(candidate_rows) ~= cfg.blades_num || numel(peak_vec_norm) ~= cfg.blades_num
        continue;
    end
    candidate_rows_cache{k} = candidate_rows;
    valid_lap_mask(k) = true;
    if numel(seed_lap_indices) < seed_lap_target
        seed_lap_indices(end + 1) = k; %#ok<AGROW>
        for shift = 0:(cfg.blades_num - 1)
            curr_shifted = circshift(peak_vec_norm(:), -shift);
            corr_val = local_corrcoef_scalar_local(curr_shifted, ref_fingerprint(:));
            shift_corr(k, shift + 1) = corr_val;
        end
    end
end

if isempty(seed_lap_indices)
    return;
end

shift_scores = mean(shift_corr(seed_lap_indices, :), 1, 'omitnan');
[best_score, best_shift_pos] = max(shift_scores);
if isempty(best_shift_pos) || ~isfinite(best_score)
    return;
end
best_shift = best_shift_pos - 1;
local_to_physical = mod((1:cfg.blades_num) - 1 - best_shift, cfg.blades_num) + 1;
physical_to_local = invert_blade_mapping_local(local_to_physical);

selected_rows_by_physical = repmat({nan(numel(selected_rev_ids), 1)}, cfg.blades_num, 1);
angle_error_by_physical = repmat({nan(numel(selected_rev_ids), 1)}, cfg.blades_num, 1);
for k = 1:numel(selected_rev_ids)
    candidate_rows = candidate_rows_cache{k};
    if numel(candidate_rows) ~= cfg.blades_num
        continue;
    end
    for blade_id = 1:cfg.blades_num
        local_slot = physical_to_local(blade_id);
        if ~isfinite(local_slot) || local_slot < 1 || local_slot > numel(candidate_rows)
            continue;
        end
        row_id = candidate_rows(local_slot);
        selected_rows_by_physical{blade_id}(k) = row_id;
        angle_error_by_physical{blade_id}(k) = 0;
    end
end

sorted_scores = sort(shift_scores, 'descend', 'MissingPlacement', 'last');
if numel(sorted_scores) >= 2 && isfinite(sorted_scores(1)) && isfinite(sorted_scores(2))
    score_margin = sorted_scores(1) - sorted_scores(2);
else
    score_margin = NaN;
end

target_physical = cfg.target_blades;
target_local = physical_to_local(target_physical);
target_corr_peak_v = shift_scores(best_shift_pos);
if isfinite(target_local)
    target_corr_peak_v = shift_scores(best_shift_pos);
end

match_result.locked = true;
match_result.best_shift = best_shift;
match_result.shift_scores = shift_scores;
match_result.best_score = best_score;
match_result.score_margin = score_margin;
match_result.valid_lap_count = nnz(valid_lap_mask);
match_result.seed_lap_count = numel(seed_lap_indices);
match_result.local_to_physical_blade_ids = local_to_physical;
match_result.physical_to_local_blade_ids = physical_to_local;
match_result.target_corr_peak_v = target_corr_peak_v;
match_result.selected_rows_by_physical = selected_rows_by_physical;
match_result.angle_error_by_physical = angle_error_by_physical;
end

function [candidate_rows, peak_vec_norm] = build_revolution_fingerprint_local(raw_sensor, probe, rev_id, cfg)
candidate_rows = find(probe.revolution_index == rev_id);
if isempty(candidate_rows)
    peak_vec_norm = [];
    return;
end
[~, ord] = sort(probe.jilublade(candidate_rows, 3), 'ascend');
candidate_rows = candidate_rows(ord);
if numel(candidate_rows) ~= cfg.blades_num
    peak_vec_norm = [];
    return;
end

peak_vec = nan(cfg.blades_num, 1);
for i = 1:cfg.blades_num
    row_id = candidate_rows(i);
    t_peak = probe.jilublade(row_id, 3);
    [t_start, t_end] = build_dynamic_segment_window_local(probe.jilublade, row_id, t_peak, cfg);
    keep = raw_sensor.T >= t_start & raw_sensor.T <= t_end;
    t_seg = raw_sensor.T(keep);
    v_seg = raw_sensor.V(keep);
    if numel(t_seg) < max(8, cfg.match_poly_degree + 2)
        peak_vec_norm = [];
        return;
    end
    peak_vec(i) = fit_polynomial_peak_value_only_local(t_seg, v_seg, cfg.match_poly_degree);
    if ~isfinite(peak_vec(i))
        peak_vec_norm = [];
        return;
    end
end

peak_vec_norm = normalize_fingerprint_max_local(peak_vec);
if any(~isfinite(peak_vec_norm))
    candidate_rows = [];
    peak_vec_norm = [];
end
end

function peak_v = fit_polynomial_peak_value_only_local(t_seg, v_seg, degree)
peak_v = nan;
t_seg = t_seg(:);
v_seg = v_seg(:);
valid = isfinite(t_seg) & isfinite(v_seg);
t_seg = t_seg(valid);
v_seg = v_seg(valid);
if numel(t_seg) < max(3, degree + 1) || range(t_seg) <= eps
    return;
end
[t_seg, ord] = sort(t_seg);
v_seg = v_seg(ord);
real_degree = min(degree, numel(t_seg) - 1);
[p, ~, mu] = polyfit(t_seg, v_seg, real_degree);
t_fine = linspace(min(t_seg), max(t_seg), 200).';
v_fine = polyval(p, t_fine, [], mu);
peak_v = max(v_fine);
end

function vec_norm = normalize_fingerprint_max_local(vec)
vec = vec(:);
vec_norm = nan(size(vec));
valid = isfinite(vec);
if nnz(valid) ~= numel(vec)
    return;
end
vmax = max(vec);
if ~isfinite(vmax) || vmax <= eps
    return;
end
vec_norm = vec / vmax;
end

function feature = build_poly_peak_feature_for_rows_local(raw_sensor, probe, row_ids, theta_std, cfg, F_omega_deg)
n = numel(row_ids);
peak_v = nan(n, 1);
peak_x = nan(n, 1);
fit_ok = false(n, 1);

for k = 1:n
    row_id = row_ids(k);
    if ~isfinite(row_id) || row_id < 1 || row_id > size(probe.jilublade, 1)
        continue;
    end
    t_peak = probe.jilublade(row_id, 3);
    [t_start, t_end] = build_dynamic_segment_window_local(probe.jilublade, row_id, t_peak, cfg);
    keep = raw_sensor.T >= t_start & raw_sensor.T <= t_end;
    t_seg = raw_sensor.T(keep);
    v_seg = raw_sensor.V(keep);
    if numel(t_seg) < max(8, cfg.match_poly_degree + 2)
        continue;
    end
    idx_prev = probe.prev_opr_index(row_id);
    if ~isfinite(idx_prev) || idx_prev < 1
        continue;
    end
    t_ref = probe.opr_times(idx_prev);
    x_abs = map_segment_to_relative_xabs_local(t_ref, t_seg, F_omega_deg, theta_std, cfg.r_tip_mm);
    [peak_v(k), peak_x(k), fit_ok(k)] = fit_polynomial_peak_local(x_abs, v_seg, cfg.match_poly_degree);
end

feature = struct();
feature.peak_v = peak_v;
feature.peak_x = peak_x;
feature.fit_ok = fit_ok;
feature.peak_v_norm = normalize_vector_local(peak_v, fit_ok);
feature.peak_x_norm = normalize_vector_local(peak_x, fit_ok);
end

function match_result = solve_cyclic_blade_match_local(anchor_bank, test_bank, target_blade, weights)
blade_count = numel(anchor_bank);
Corr = nan(blade_count, blade_count);
CorrPeakV = nan(blade_count, blade_count);
CorrPeakX = nan(blade_count, blade_count);
FitCount = nan(blade_count, blade_count);
use_peak_x = use_peak_x_local(weights);

for ba = 1:blade_count
    for bb = 1:blade_count
        anchor_feature = anchor_bank(ba);
        test_feature = test_bank(bb);
        CorrPeakV(ba, bb) = compute_single_feature_corr_local(anchor_feature, test_feature, 'peak_v_norm');
        CorrPeakX(ba, bb) = compute_single_feature_corr_local(anchor_feature, test_feature, 'peak_x_norm');
        Corr(ba, bb) = combine_feature_score_local(CorrPeakV(ba, bb), CorrPeakX(ba, bb), weights, use_peak_x);
        n = min(numel(anchor_feature.fit_ok), numel(test_feature.fit_ok));
        FitCount(ba, bb) = nnz(anchor_feature.fit_ok(1:n) & test_feature.fit_ok(1:n));
    end
end

shift_scores = nan(blade_count, 1);
for shift = 0:(blade_count - 1)
    mapping = apply_cyclic_shift_local(blade_count, shift);
    diag_scores = nan(blade_count, 1);
    for ba = 1:blade_count
        diag_scores(ba) = Corr(ba, mapping(ba));
    end
    shift_scores(shift + 1) = mean(diag_scores, 'omitnan');
end

[best_score, best_shift_pos] = max(shift_scores);
if isempty(best_shift_pos) || ~isfinite(best_score)
    best_shift = NaN;
    best_mapping = 1:blade_count;
    matched_blade_id = NaN;
    target_corr_peak_v = NaN;
    target_corr_peak_x = NaN;
    target_fit_valid_count = NaN;
    score_margin = NaN;
else
    best_shift = best_shift_pos - 1;
    best_mapping = apply_cyclic_shift_local(blade_count, best_shift);
    matched_blade_id = best_mapping(target_blade);
    target_corr_peak_v = CorrPeakV(target_blade, matched_blade_id);
    target_corr_peak_x = CorrPeakX(target_blade, matched_blade_id);
    target_fit_valid_count = FitCount(target_blade, matched_blade_id);
    sorted_scores = sort(shift_scores, 'descend', 'MissingPlacement', 'last');
    if numel(sorted_scores) >= 2 && isfinite(sorted_scores(1)) && isfinite(sorted_scores(2))
        score_margin = sorted_scores(1) - sorted_scores(2);
    else
        score_margin = NaN;
    end
end

match_result = struct();
match_result.Corr = Corr;
match_result.CorrPeakV = CorrPeakV;
match_result.CorrPeakX = CorrPeakX;
match_result.FitCount = FitCount;
match_result.shift_scores = shift_scores;
match_result.best_shift = best_shift;
match_result.best_mapping = best_mapping;
match_result.best_score = best_score;
match_result.matched_blade_id = matched_blade_id;
match_result.target_corr_peak_v = target_corr_peak_v;
match_result.target_corr_peak_x = target_corr_peak_x;
match_result.target_fit_valid_count = target_fit_valid_count;
match_result.score_margin = score_margin;
end

function mapping = apply_cyclic_shift_local(blade_count, shift)
mapping = mod((1:blade_count) - 1 + shift, blade_count) + 1;
end

function inverse_mapping = invert_blade_mapping_local(mapping)
inverse_mapping = nan(size(mapping));
for i = 1:numel(mapping)
    local_blade_id = mapping(i);
    if isfinite(local_blade_id) && local_blade_id >= 1 && local_blade_id <= numel(mapping)
        inverse_mapping(local_blade_id) = i;
    end
end
end

function blade_numbering = build_blade_numbering_info_local(diag, cfg)
sensor_projection = diag.sensor_projection;
sensor_ids = reshape([sensor_projection.sensor_id], [], 1);
physical_to_local = nan(cfg.blades_num, numel(sensor_projection));
local_to_physical = nan(cfg.blades_num, numel(sensor_projection));
blade_index_offset = nan(numel(sensor_projection), 1);
match_score = nan(numel(sensor_projection), 1);
match_margin = nan(numel(sensor_projection), 1);
selected_rows_by_physical = cell(numel(sensor_projection), 1);

for i = 1:numel(sensor_projection)
    physical_to_local(:, i) = sensor_projection(i).physical_to_local_blade_ids(:);
    local_to_physical(:, i) = sensor_projection(i).local_to_physical_blade_ids(:);
    blade_index_offset(i) = sensor_projection(i).blade_index_offset;
    match_score(i) = sensor_projection(i).match_score;
    match_margin(i) = sensor_projection(i).match_margin;
    selected_rows_by_physical{i} = sensor_projection(i).selected_rows_by_physical;
end

blade_numbering = struct();
blade_numbering.mode = 'six_peak_fingerprint_lock';
blade_numbering.numbering_seed_laps = cfg.numbering_seed_laps;
blade_numbering.anchor_sensor_id = diag.anchor_sensor_id;
blade_numbering.sensor_ids = sensor_ids;
blade_numbering.physical_blade_ids = (1:cfg.blades_num).';
blade_numbering.target_blade = cfg.target_blades;
blade_numbering.physical_to_local_blade_ids = physical_to_local;
blade_numbering.local_to_physical_blade_ids = local_to_physical;
blade_numbering.blade_index_offset = blade_index_offset;
blade_numbering.match_score = match_score;
blade_numbering.match_margin = match_margin;
blade_numbering.selected_rows_by_physical = selected_rows_by_physical;
blade_numbering.description = ['Columns correspond to sensor_ids. physical_to_local_blade_ids(p,s) gives the ' ...
    'local pulse slot on sensor s that corresponds to physical blade p after six-peak fingerprint locking.'];
end

function [t_start, t_end] = build_dynamic_segment_window_local(jilublade, row_id, t_peak, cfg)
if strcmpi(cfg.dynamic_window_mode, 'legacy_row_bounds')
    t_start = jilublade(row_id, 1) - cfg.pulse_pad_sec;
    t_end = jilublade(row_id, 2) + cfg.pulse_pad_sec;
else
    t_start = t_peak - cfg.pulse_window_sec;
    t_end = t_peak + cfg.pulse_window_sec;
end
end

function x_abs = map_segment_to_relative_xabs_local(t_ref, t_seg, F_omega_deg, theta_std, r_tip_mm)
theta_points_deg = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg);
theta_diff_deg = wrap_to_signed_period_local(theta_points_deg - theta_std, 360);
x_abs = theta_diff_deg * (pi / 180) * r_tip_mm;
end

function [peak_v, peak_x, ok] = fit_polynomial_peak_local(x, v, degree)
peak_v = nan;
peak_x = nan;
ok = false;
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
if numel(x) < degree + 2 || range(x) <= eps
    return;
end
[x, order] = sort(x);
v = v(order);
x_center = mean(x);
x_scale = max(std(x), eps);
xn = (x - x_center) / x_scale;
coef = polyfit(xn, v, degree);
x_dense = linspace(min(x), max(x), 500).';
v_fit = polyval(coef, (x_dense - x_center) / x_scale);
[peak_v, idx] = max(v_fit);
peak_x = x_dense(idx);
range_x = max(x_dense) - min(x_dense);
left_guard = min(x_dense) + 0.10 * range_x;
right_guard = max(x_dense) - 0.10 * range_x;
ok = isfinite(peak_v) && isfinite(peak_x) && peak_x >= left_guard && peak_x <= right_guard;
if ~ok
    peak_v = nan;
    peak_x = nan;
end
end

function corr_value = compute_single_feature_corr_local(A, B, field_name)
va = A.(field_name);
vb = B.(field_name);
n = min(numel(va), numel(vb));
if n < 3
    corr_value = nan;
    return;
end
valid = A.fit_ok(1:n) & B.fit_ok(1:n) & isfinite(va(1:n)) & isfinite(vb(1:n));
corr_value = local_corrcoef_scalar_local(va(valid), vb(valid));
end

function value = local_corrcoef_scalar_local(a, b)
if numel(a) < 3 || numel(b) < 3 || std(a, 'omitnan') <= eps || std(b, 'omitnan') <= eps
    value = nan;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function score = combine_feature_score_local(corr_peak_v, corr_peak_x, weights, use_peak_x)
parts = [];
part_weights = [];
if isfinite(corr_peak_v)
    parts(end + 1, 1) = corr_peak_v; %#ok<AGROW>
    part_weights(end + 1, 1) = weights.peak_v; %#ok<AGROW>
end
if use_peak_x && isfinite(corr_peak_x)
    parts(end + 1, 1) = corr_peak_x; %#ok<AGROW>
    part_weights(end + 1, 1) = weights.peak_x; %#ok<AGROW>
end
if isempty(parts)
    score = nan;
    return;
end
score = sum(parts .* part_weights) / sum(part_weights);
end

function tf = use_peak_x_local(weights)
tf = isstruct(weights) && isfield(weights, 'peak_x') && isfinite(weights.peak_x) && weights.peak_x > 0;
end

function v_norm = normalize_vector_local(v, fit_ok)
v = v(:);
fit_ok = fit_ok(:);
v_norm = nan(size(v));
valid = fit_ok & isfinite(v);
if nnz(valid) < 3
    return;
end
mu = mean(v(valid), 'omitnan');
sigma = std(v(valid), 'omitnan');
if ~isfinite(sigma) || sigma <= eps
    v_norm(valid) = 0;
    return;
end
v_norm(valid) = (v(valid) - mu) / sigma;
end

function anchor_sid = choose_anchor_sensor_local(probe_cache, cfg, Sensor_Config)
anchor_candidates = cfg.analysis_sensors(:).';
if isfield(Sensor_Config, 'Capacitance_IDs') && ~isempty(Sensor_Config.Capacitance_IDs)
    preferred = Sensor_Config.Capacitance_IDs(:).';
    preferred = preferred(ismember(preferred, anchor_candidates));
    anchor_candidates = unique([preferred, anchor_candidates], 'stable');
end
for sid = anchor_candidates
    cache_idx = find([probe_cache.sensor_id] == sid, 1, 'first');
    if isempty(cache_idx)
        continue;
    end
    jilublade = probe_cache(cache_idx).jilublade;
    mask = jilublade(:, 4) == cfg.target_blades & jilublade(:, 3) >= cfg.analysis_start_time;
    if nnz(mask) >= cfg.target_laps
        anchor_sid = sid;
        return;
    end
end
error('No anchor sensor has at least %d identified B%d laps after %.3f s.', ...
    cfg.target_laps, cfg.target_blades, cfg.analysis_start_time);
end

function selected_rows = enforce_strictly_increasing_rows_local(selected_rows, sid)
selected_rows = round(selected_rows(:));
if any(~isfinite(selected_rows))
    error('Projected row selection contains NaN for CH%d.', sid);
end
if any(diff(selected_rows) <= 0)
    error('Projected rows are not strictly increasing for CH%d.', sid);
end
end

function angle = wrap_to_signed_period_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function explain_missing_template_local(template_file, route_dir, template_suffix)
template_dir = fullfile(route_dir, 'output', 'templates');
nearby = [ ...
    dir(fullfile(template_dir, 'TemplateBundle_LowSpeedRotating_B*_S*_20241106.mat')); ...
    dir(fullfile(template_dir, 'Template_LowSpeedRotating_B*_CH*_20241106.mat')); ...
    dir(fullfile(template_dir, 'Template_LowSpeedRotating_B*_S*_20241106.mat'))];
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

function Template = load_template_for_sensors_local(template_file, analysis_sensors)
loaded = load(template_file);
if isfield(loaded, 'TemplateBundle')
    Template = load_template_from_bundle_local(loaded.TemplateBundle, analysis_sensors);
    return;
end
if isfield(loaded, 'sensor_template')
    Template = normalize_sensor_template_local(loaded.sensor_template, analysis_sensors, template_file);
    return;
end
if isfield(loaded, 'Template')
    Template = normalize_sensor_template_local(loaded.Template, analysis_sensors, template_file);
    return;
end
error('Unsupported template file format: %s', template_file);
end

function Template = load_template_from_bundle_local(bundle, analysis_sensors)
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
    if height(row) ~= 1
        error('Template bundle does not contain CH%d.', sid);
    end
    sensor_loaded = load(char(row.template_file(1)));
    sensor_template = extract_single_sensor_template_local(sensor_loaded, sid, char(row.template_file(1)));
    if is == 1
        Template = sensor_template;
        sensor_entries = repmat(sensor_template.Sensor, numel(analysis_sensors), 1);
    end
    sensor_entries(is) = sensor_template.Sensor;
end
Template.Sensor = sensor_entries;
Template.SensorIDs = analysis_sensors;
Template.SensorTag = ['S', sprintf('%d', analysis_sensors)];
end

function Template = normalize_sensor_template_local(template_in, analysis_sensors, source_file)
if ~isfield(template_in, 'Sensor') || isempty(template_in.Sensor)
    error('Template file missing Sensor field: %s', source_file);
end
sensor_struct = template_in.Sensor;
if numel(sensor_struct) == 1
    sid = sensor_struct.sensor_id;
    if numel(analysis_sensors) ~= 1 || sid ~= analysis_sensors(1)
        error('Single-sensor template %s does not match requested sensors %s.', ...
            source_file, mat2str(analysis_sensors));
    end
    Template = template_in;
    Template.SensorIDs = analysis_sensors;
    Template.SensorTag = ['S', sprintf('%d', analysis_sensors)];
    return;
end
available = [sensor_struct.sensor_id];
sensor_idx = zeros(size(analysis_sensors));
for is = 1:numel(analysis_sensors)
    idx = find(available == analysis_sensors(is), 1, 'first');
    if isempty(idx)
        error('Template file %s does not contain CH%d.', source_file, analysis_sensors(is));
    end
    sensor_idx(is) = idx;
end
Template = template_in;
Template.Sensor = sensor_struct(sensor_idx);
Template.SensorIDs = analysis_sensors;
Template.SensorTag = ['S', sprintf('%d', analysis_sensors)];
end

function sensor_template = extract_single_sensor_template_local(sensor_loaded, sid, source_file)
if isfield(sensor_loaded, 'sensor_template')
    sensor_template = normalize_sensor_template_local(sensor_loaded.sensor_template, sid, source_file);
    return;
end
if isfield(sensor_loaded, 'Template')
    sensor_template = normalize_sensor_template_local(sensor_loaded.Template, sid, source_file);
    return;
end
error('Sensor template file is unsupported: %s', source_file);
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

function [rel_angles_deg, prev_opr_index, revolution_index] = compute_pulse_reference_geometry_local(arrival_times, opr_times, opr_pulses_per_rev)
rel_angles_deg = nan(size(arrival_times));
prev_opr_index = nan(size(arrival_times));
revolution_index = nan(size(arrival_times));
if numel(opr_times) <= opr_pulses_per_rev
    return;
end
F_omega_deg = build_phase_speed_local(opr_times, opr_pulses_per_rev);
for i = 1:numel(arrival_times)
    t_meas = arrival_times(i);
    idx_prev_opr = find(opr_times < t_meas, 1, 'last');
    if isempty(idx_prev_opr)
        continue;
    end
    t_ref = opr_times(idx_prev_opr);
    t_grid = linspace(t_ref, t_meas, 10);
    rel_angles_deg(i) = trapz(t_grid, F_omega_deg(t_grid));
    prev_opr_index(i) = idx_prev_opr;
    revolution_index(i) = floor((idx_prev_opr - 1) / opr_pulses_per_rev) + 1;
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

