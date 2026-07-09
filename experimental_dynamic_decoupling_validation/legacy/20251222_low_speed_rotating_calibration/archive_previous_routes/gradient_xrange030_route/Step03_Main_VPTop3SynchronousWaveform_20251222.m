%% Step03_Main_VPTop3SynchronousWaveform_20251222
% Main low-speed-template-only identification route.
% First-order template VP screens EO candidates; the final result still
% comes from the full synchronous waveform objective.
%
% Every integer EO candidate is scored by the first-order model
% V - T(x) ~= -T'(x) * [dx_c + a*sin(EO*theta) + b*cos(EO*theta)].
% By default, the top-3 EO candidates enter full waveform refinement in
% every window. For synchronous vibration, frequency is constrained by
% f = EO * rot_freq_mean; the final refinement optimizes amplitude/phase
% and the spatial alignment term dx_c.

clear; clc; close all;

%% Settings matching the super-Gaussian Step5 route
route_dir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1, 2, 3];
sensor_override = strtrim(getenv('STEP03_ANALYSIS_SENSORS'));
if ~isempty(sensor_override)
    parsed_sensors = sscanf(sensor_override, '%d').';
    if isempty(parsed_sensors)
        error('STEP03_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 2 3".');
    end
    cfg.analysis_sensors = parsed_sensors;
end
cfg.analysis_start_time = 1.5;
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
cfg.freq_search_hz = [100, 1000];
cfg.eo_pad = 2;
cfg.amplitude_limit_mm = 0.50;
cfg.dx_c_limit_mm = 0.20;
cfg.debug_max_windows = parse_positive_integer_env_local('STEP03_04_MAX_WINDOWS', inf);
cfg.vp_top_k_eo = parse_positive_integer_env_local('STEP03_MAIN_TOP_K_EO', 3);
cfg.vp_all_eo_warmup_windows = parse_nonnegative_integer_env_local('STEP03_04_ALL_EO_WARMUP_WINDOWS', 0);
cfg.vp_gap_ratio_fallback = parse_numeric_env_local('STEP03_04_GAP_RATIO_FALLBACK', -inf);
cfg.vp_linear_gap_ratio_fallback = parse_numeric_env_local('STEP03_04_LINEAR_GAP_RATIO_FALLBACK', -inf);
cfg.diagnostic_eo = parse_optional_integer_env_local('STEP03_DIAGNOSTIC_EO');
cfg.sensor_eta_limit_mm = parse_nonnegative_numeric_env_local('STEP03D_SENSOR_ETA_LIMIT_MM', 0.08);
cfg.sensor_eta_reg_weight_v_per_mm = parse_nonnegative_numeric_env_local('STEP03D_SENSOR_ETA_REG_WEIGHT_V_PER_MM', 1.00);
cfg.pulse_selection_mode = lower(strtrim(getenv('STEP03D_PULSE_MODE')));
if isempty(cfg.pulse_selection_mode)
    cfg.pulse_selection_mode = 'single';
end
if ~ismember(cfg.pulse_selection_mode, {'single', 'all'})
    error('STEP03D_PULSE_MODE must be "single" or "all".');
end

method = struct();
method.name = 'template_only_main_vp_top3_synchronous_waveform';
method.selection_rule = 'first_order_template_vp_top3_then_synchronous_waveform_rmse';
method.use_reference_eo_constraint = false;
method.use_vp = true;
method.vp_top_k_eo = cfg.vp_top_k_eo;
method.refine_params = {'A', 'phi', 'dx_c', 'sensor_eta'};
method.template_forward = 'interp1_low_speed_template';
method.vp_seed_model = 'V_minus_Tx_equals_minus_Tprime_times_u';
method.final_waveform_objective = 'fixed_eo_direct_template_voltage_residual_without_sensor_affine_projection';
method.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
method.sensor_eta_reg_weight_v_per_mm = cfg.sensor_eta_reg_weight_v_per_mm;
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.query_guard_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm + method.coverage_safety_margin_mm;
method.query_guard_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MM', method.query_guard_mm);
method.query_guard_mode = lower(strtrim(getenv('STEP03D_QUERY_GUARD_MODE')));
if isempty(method.query_guard_mode)
    method.query_guard_mode = 'adaptive';
end
if ~ismember(method.query_guard_mode, {'fixed', 'adaptive'})
    error('STEP03D_QUERY_GUARD_MODE must be "fixed" or "adaptive".');
end
method.query_guard_quantile = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_QUANTILE', 95);
method.query_guard_safety_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_SAFETY_MM', method.coverage_safety_margin_mm);
method.query_guard_min_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MIN_MM', 0.12);
method.query_guard_max_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MAX_MM', method.query_guard_mm);
method.domain_selection_mode = lower(strtrim(getenv('STEP03D_DOMAIN_SELECTION_MODE')));
if isempty(method.domain_selection_mode)
    method.domain_selection_mode = 'hard';
end
if ~ismember(method.domain_selection_mode, {'hard', 'soft'})
    error('STEP03D_DOMAIN_SELECTION_MODE must be "hard" or "soft".');
end
method.domain_soft_margin_mm = parse_nonnegative_numeric_env_local('STEP03D_DOMAIN_SOFT_MARGIN_MM', 0);
method.overshoot_penalty_weight = parse_nonnegative_numeric_env_local('STEP03D_OVERSHOOT_PENALTY_WEIGHT', 100);
method.coverage_mode = 'prefer_points_safe_for_all_bounded_u';
method.extrapolation_mode = 'clamp_to_template_edge';
method.default_sensor_threshold = 0.5;
method.pulse_selection_mode = cfg.pulse_selection_mode;
method.dynamic_effective_mode = lower(strtrim(getenv('STEP03D_DYNAMIC_EFFECTIVE_MODE')));
if isempty(method.dynamic_effective_mode)
    method.dynamic_effective_mode = 'gradient';
end
if ~ismember(method.dynamic_effective_mode, {'legacy', 'gradient'})
    error('STEP03D_DYNAMIC_EFFECTIVE_MODE must be "legacy" or "gradient".');
end
method.dynamic_template_gradient_min_ratio = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_TEMPLATE_GRADIENT_MIN_RATIO', 0.08);
method.dynamic_time_gradient_min_ratio = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_TIME_GRADIENT_MIN_RATIO', 0.15);
method.dynamic_peak_quantile = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_PEAK_QUANTILE', 85);

sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
result_dir = fullfile(output_dir, 'identification');
if exist(result_dir, 'dir') ~= 7; mkdir(result_dir); end

template_file_override = strtrim(getenv('STEP03_TEMPLATE_FILE'));
if ~isempty(template_file_override)
    if exist(template_file_override, 'file') ~= 2
        error('STEP03_TEMPLATE_FILE does not exist: %s', template_file_override);
    end
    template_file = template_file_override;
else
    template_file = find_template_file_for_sensors_local(template_dir, cfg.target_blades, cfg.analysis_sensors);
end
dynamic_file_override = strtrim(getenv('STEP03_DYNAMIC_MAP_FILE'));
if ~isempty(dynamic_file_override)
    if exist(dynamic_file_override, 'file') ~= 2
        error('STEP03_DYNAMIC_MAP_FILE does not exist: %s', dynamic_file_override);
    end
    dynamic_file = dynamic_file_override;
else
    dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, cfg.target_blades, cfg.analysis_sensors);
end
result_suffix = sanitize_result_suffix_local(strtrim(getenv('STEP03_RESULT_SUFFIX')));
result_file = fullfile(result_dir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s%s_20251222.mat', ...
    cfg.target_blades, sensor_tag, result_suffix));

loaded_template = load(template_file, 'Template');
loaded_dynamic = load(dynamic_file, 'DynamicMap');
Template = filter_template_sensors_local(loaded_template.Template, cfg.analysis_sensors, sensor_tag);
DynamicMap = filter_dynamic_map_sensors_local(loaded_dynamic.DynamicMap, cfg.analysis_sensors, sensor_tag);
coordinate_check = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysis_sensors, 1e-6);
source_cfg = DynamicMap.SourceSettings;
cfg.analysis_start_time = source_cfg.analysis_start_time;
cfg.target_laps = source_cfg.target_laps;
cfg.analysis_win_size = source_cfg.analysis_win_size;
cfg.sliding_step = source_cfg.sliding_step;

fprintf('\n=== Step03 main: VP top-3 synchronous waveform identification ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Template source: %s\n', template_file);
fprintf('Dynamic map source: %s\n', dynamic_file);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinate_check.max_abs_delta_mm, coordinate_check.status);
if ~coordinate_check.is_consistent
    warning('Template/DynamicMap x-center mismatch. This run is not a clean GradientXRange030 coordinate chain.');
end
fprintf('Frequency search: %.1f-%.1f Hz; no reference EO; first-order VP only screens candidates.\n', ...
    cfg.freq_search_hz(1), cfg.freq_search_hz(2));
fprintf('Model: V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi)); synchronous f = EO*rot_freq_mean.\n');
fprintf(['Template coverage guard: %s, fixed %.3f mm, min %.3f mm, max %.3f mm, ' ...
    'q%.1f + %.3f mm.\n'], method.query_guard_mode, method.query_guard_mm, ...
    method.query_guard_min_mm, method.query_guard_max_mm, method.query_guard_quantile, ...
    method.query_guard_safety_mm);
fprintf('Pulse selection mode: %s.\n', method.pulse_selection_mode);
fprintf('Dynamic effective selection: %s.\n', method.dynamic_effective_mode);
fprintf('Domain selection: %s; soft margin %.3f mm; overshoot penalty %.3g.\n', ...
    method.domain_selection_mode, method.domain_soft_margin_mm, method.overshoot_penalty_weight);
fprintf('First-order VP/adaptive: keep top %d EO candidates in every window by default.\n', ...
    method.vp_top_k_eo);
fprintf('Sensor eta limit: %.4f mm (no EO prior; per-sensor x-zero correction).\n', ...
    method.sensor_eta_limit_mm);
fprintf('Sensor eta regularization: %.4f V/mm.\n', cfg.sensor_eta_reg_weight_v_per_mm);
if cfg.vp_all_eo_warmup_windows > 0 || isfinite(cfg.vp_gap_ratio_fallback) || isfinite(cfg.vp_linear_gap_ratio_fallback)
    fprintf('Optional fallback enabled: warmup=%d, weighted gap=%.6g, linear gap=%.6g.\n', ...
        cfg.vp_all_eo_warmup_windows, cfg.vp_gap_ratio_fallback, cfg.vp_linear_gap_ratio_fallback);
else
    fprintf('Optional fallback disabled: no warmup and no all-EO ambiguity fallback.\n');
end
if isfinite(cfg.diagnostic_eo)
    fprintf('Diagnostic only: report EO%d rank after VP screening and full waveform refinement.\n', cfg.diagnostic_eo);
end

num_windows = numel(DynamicMap.Window);
if isfinite(cfg.debug_max_windows)
    num_windows = min(num_windows, cfg.debug_max_windows);
end

trend_rows = repmat(struct( ...
    'window_id', NaN, 'lap_start', NaN, 'lap_end', NaN, 'window_center_time', NaN, ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, 'A_id', NaN, 'EO_id', NaN, ...
    'fn_id', NaN, 'phi_id_wrapped', NaN, 'dx_c_id', NaN, 'd0_id', NaN, ...
    'sensor_eta_max_abs_mm', NaN, ...
    'weighted_voltage_rmse', NaN, 'plain_voltage_rmse', NaN, ...
    'valid_segment_count', NaN, 'point_count', NaN), num_windows, 1);
window_results = repmat(struct(), num_windows, 1);
best_window = struct('weighted_voltage_rmse', inf);

for w_idx = 1:num_windows
    Wmap = DynamicMap.Window(w_idx);
    bundle = build_template_observation_bundle_local(Wmap, Template, cfg.analysis_sensors, method);
    eo_candidates = build_eo_candidates_local(bundle.rot_freq_mean_hz, cfg.freq_search_hz, cfg.eo_pad);
    seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
    [final_eo_candidates, selection_info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
    result = refine_template_waveform_fit_local(bundle, seed_table, final_eo_candidates, cfg);
    result.VPSeedTable = struct2table(seed_table);
    result.VPSelectedEO = final_eo_candidates(:).';
    result.VPSelectionInfo = selection_info;

    trend_rows(w_idx).window_id = w_idx;
    trend_rows(w_idx).lap_start = Wmap.lap_range(1);
    trend_rows(w_idx).lap_end = Wmap.lap_range(end);
    trend_rows(w_idx).window_center_time = median(bundle.T, 'omitnan');
    trend_rows(w_idx).rot_freq_mean_hz = bundle.rot_freq_mean_hz;
    trend_rows(w_idx).rot_rpm_mean = Wmap.rot_rpm_mean;
    trend_rows(w_idx).A_id = result.A_id;
    trend_rows(w_idx).EO_id = result.EO_id;
    trend_rows(w_idx).fn_id = result.fn_id;
    trend_rows(w_idx).phi_id_wrapped = result.phi_id_wrapped;
    trend_rows(w_idx).dx_c_id = result.dx_c_id;
    trend_rows(w_idx).d0_id = result.dx_c_id;
    trend_rows(w_idx).sensor_eta_max_abs_mm = max(abs(result.sensor_eta_id), [], 'omitnan');
    trend_rows(w_idx).weighted_voltage_rmse = result.weighted_voltage_rmse;
    trend_rows(w_idx).plain_voltage_rmse = result.plain_voltage_rmse;
    trend_rows(w_idx).valid_segment_count = result.valid_segment_count;
    trend_rows(w_idx).point_count = result.point_count;

    window_results(w_idx).window_id = w_idx;
    window_results(w_idx).lap_range = Wmap.lap_range;
    window_results(w_idx).time_window = Wmap.time_window;
    window_results(w_idx).bundle = bundle;
    window_results(w_idx).seed_table = seed_table;
    window_results(w_idx).CandidateTable = result.CandidateTable;
    window_results(w_idx).Result = result;

    if result.weighted_voltage_rmse < best_window.weighted_voltage_rmse
        best_window = result;
        best_window.window_id = w_idx;
        best_window.lap_range = Wmap.lap_range;
        best_window.seed_table = seed_table;
    end

    diagnostic_text = '';
    if isfinite(cfg.diagnostic_eo)
        diag_seed_rank = find([seed_table.EO] == cfg.diagnostic_eo, 1);
        diag_final_rank = find(result.CandidateTable.EO == cfg.diagnostic_eo, 1);
        diagnostic_text = sprintf(', EO%d seed/final-rank=%s/%s', ...
            cfg.diagnostic_eo, rank_to_string_local(diag_seed_rank), rank_to_string_local(diag_final_rank));
    end
    fprintf('Window %02d/%02d laps %s: EO=%d, f=%.3f Hz, A=%.4f mm, RMSE=%.5f V, clamp=%.2f%%%s\n', ...
        w_idx, num_windows, mat2str(Wmap.lap_range), result.EO_id, result.fn_id, ...
        result.A_id, result.weighted_voltage_rmse, 100 * result.Coverage.clamp_fraction, diagnostic_text);
end

Trend = struct2table(trend_rows);
Result = struct();
Result.Method = method.name;
Result.AnalysisSettings = cfg;
Result.MethodSettings = method;
Result.TemplateFile = template_file;
Result.DynamicMapFile = dynamic_file;
Result.CoordinateCheck = coordinate_check;
Result.TargetBlade = cfg.target_blades;
Result.SensorIDs = cfg.analysis_sensors;
Result.SensorTag = sensor_tag;
Result.Trend = Trend;
Result.WindowResult = window_results;
Result.BestWindow = best_window;
Result.ResonanceSummary = build_resonance_summary_local(trend_rows);

save(result_file, 'Result', '-v7.3');

fprintf('\n=== Step03 main trend ===\n');
disp(Trend(:, {'window_id','lap_start','lap_end','EO_id','fn_id','A_id','weighted_voltage_rmse'}));
fprintf('Dominant EO: %d; mean f = %.6f Hz; median f = %.6f Hz\n', ...
    Result.ResonanceSummary.dominant_eo, Result.ResonanceSummary.mean_freq_hz, ...
    Result.ResonanceSummary.median_freq_hz);
fprintf('Saved result: %s\n', result_file);

%% Local functions
function v = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = default_value;
else
    v = max(1, floor(tmp));
end
end

function v = parse_nonnegative_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp < 0
    v = default_value;
else
    v = floor(tmp);
end
end

function v = parse_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp)
    v = default_value;
else
    v = tmp;
end
end

function v = parse_nonnegative_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
parsed = str2double(raw);
if ~isfinite(parsed) || parsed < 0
    error('%s must be a nonnegative numeric value.', name);
end
v = parsed;
end

function v = parse_optional_integer_env_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    v = NaN;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = NaN;
else
    v = floor(tmp);
end
end

function suffix = sanitize_result_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function template_file = find_template_file_for_sensors_local(template_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20251222.mat', target_blade, sensor_tag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', target_blade, sensor_tag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', target_blade, sensor_tag)
    };
for ip = 1:numel(patterns)
    candidate = fullfile(template_dir, patterns{ip});
    if exist(candidate, 'file') == 2
        template_file = candidate;
        return;
    end
end
files = dir(fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_*.mat', target_blade)));
if isempty(files)
    error('No low-speed template file found under %s.', template_dir);
end
error('Template for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
patterns = {
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20251222.mat', target_blade, sensor_tag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20251222.mat', target_blade, sensor_tag)
    };
for ip = 1:numel(patterns)
    candidate = fullfile(dynamic_dir, patterns{ip});
    if exist(candidate, 'file') == 2
        dynamic_file = candidate;
        return;
    end
end
files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_B%d_*.mat', target_blade)));
if isempty(files)
    error('No dynamic map file found under %s.', dynamic_dir);
end
error('DynamicMap for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function Template = filter_template_sensors_local(Template, analysis_sensors, sensor_tag)
available = [Template.Sensor.sensor_id];
sensor_idx = zeros(size(analysis_sensors));
for is = 1:numel(analysis_sensors)
    idx = find(available == analysis_sensors(is), 1);
    if isempty(idx)
        error('Template does not contain sensor %d.', analysis_sensors(is));
    end
    sensor_idx(is) = idx;
end
Template.Sensor = Template.Sensor(sensor_idx);
Template.SensorIDs = analysis_sensors;
Template.SensorTag = sensor_tag;
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, analysis_sensors, sensor_tag)
for iw = 1:numel(DynamicMap.Window)
    available = [DynamicMap.Window(iw).Sensor.sensor_id];
    sensor_idx = zeros(size(analysis_sensors));
    for is = 1:numel(analysis_sensors)
        idx = find(available == analysis_sensors(is), 1);
        if isempty(idx)
            error('DynamicMap window %d does not contain sensor %d.', iw, analysis_sensors(is));
        end
        sensor_idx(is) = idx;
    end
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(sensor_idx);
end
DynamicMap.SensorIDs = analysis_sensors;
DynamicMap.SensorTag = sensor_tag;
end

function check = check_template_dynamic_xcenter_local(Template, DynamicMap, analysis_sensors, tolerance_mm)
rows = repmat(struct('sensor_id', NaN, 'template_xc_mm', NaN, ...
    'dynamic_xc_mm', NaN, 'delta_mm', NaN), numel(analysis_sensors), 1);
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    rows(is).sensor_id = sid;
    if ~isempty(Tpl) && isfield(Tpl, 'xc')
        rows(is).template_xc_mm = Tpl.xc;
    end
    rows(is).dynamic_xc_mm = infer_dynamic_xcenter_local(DynamicMap, sid);
    rows(is).delta_mm = rows(is).template_xc_mm - rows(is).dynamic_xc_mm;
end
delta = [rows.delta_mm];
check = struct();
check.table = struct2table(rows);
check.tolerance_mm = tolerance_mm;
check.max_abs_delta_mm = max(abs(delta), [], 'omitnan');
check.is_consistent = isfinite(check.max_abs_delta_mm) && check.max_abs_delta_mm <= tolerance_mm;
if check.is_consistent
    check.status = 'consistent';
else
    check.status = 'mismatch';
end
if isfield(DynamicMap, 'TemplateFileForXRel')
    check.dynamic_template_file_for_xrel = DynamicMap.TemplateFileForXRel;
end
end

function xc = infer_dynamic_xcenter_local(DynamicMap, sid)
xc = NaN;
for iw = 1:numel(DynamicMap.Window)
    S = DynamicMap.Window(iw).Sensor;
    idx = find([S.sensor_id] == sid, 1, 'first');
    if isempty(idx) || isempty(S(idx).x_abs) || isempty(S(idx).x_rel)
        continue;
    end
    v = S(idx).x_abs(:) - S(idx).x_rel(:);
    xc = median(v(isfinite(v)), 'omitnan');
    return;
end
end

function bundle = build_template_observation_bundle_local(Wmap, Template, analysis_sensors, method)
X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
Y_obs = [];
F0_all = [];
Fx_all = [];
sensor_index_pt = [];
window_meta = struct('sensor_id', {}, 'valid_points', {}, 'shape_rmse', {});

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    v_raw = D.V(:);
    x_raw = D.x_rel(:);
    t_raw = D.t(:);
    theta_raw = D.theta(:);
    f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_raw, 'pchip', NaN);
    fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x_raw, 'pchip', NaN);

    if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
        threshold = Tpl.threshold;
    else
        threshold = method.default_sensor_threshold;
    end
    x_domain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    if strcmpi(method.pulse_selection_mode, 'all')
        [mask_pulse, pulse_segment_count] = isolate_all_pulses_local(v_raw, threshold);
    else
        [mask_pulse, pulse_segment_count] = isolate_main_pulse_local(v_raw, threshold);
    end
    mask_effective = build_dynamic_effective_mask_local(x_raw, v_raw, t_raw, Tpl, threshold, method);
    mask_domain = x_raw >= Tpl.x_domain(1) + method.domain_margin_mm & ...
                  x_raw <= Tpl.x_domain(2) - method.domain_margin_mm;
    mask_finite = isfinite(f0) & isfinite(fx0) & isfinite(v_raw) & isfinite(theta_raw);
    mask_base = mask_pulse & mask_effective & mask_domain & mask_finite;
    query_guard_mm = resolve_query_guard_mm_local(Tpl, x_raw, v_raw, mask_base, method);
    mask_query_safe = x_raw >= x_domain(1) + query_guard_mm & ...
                      x_raw <= x_domain(2) - query_guard_mm;
    mask = mask_base;
    if strcmpi(method.domain_selection_mode, 'hard')
        safe_mask = mask & mask_query_safe;
        if nnz(safe_mask) >= 8
            mask = safe_mask;
        end
    end
    if nnz(mask) < 8
        mask = mask_effective & mask_domain & mask_finite;
        if strcmpi(method.domain_selection_mode, 'hard')
            safe_mask = mask & mask_query_safe;
            if nnz(safe_mask) >= 8
                mask = safe_mask;
            end
        end
    end
    if nnz(mask) < 8
        continue;
    end

    t_sel = t_raw(mask);
    x_sel = x_raw(mask);
    v_sel = v_raw(mask);
    theta_sel = theta_raw(mask);
    f0_sel = f0(mask);
    fx_sel = fx0(mask);
    y_obs_sel = x_sel - invert_template_voltage_local(Tpl, v_sel, x_sel);

    w_edge = build_edge_weight_local(t_sel, v_sel, method.weight_floor);
    if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
        w_template = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x_sel, 'linear', method.weight_floor);
    else
        w_template = ones(size(x_sel));
    end
    w_domain = build_domain_soft_weight_local(x_sel, x_domain, method.domain_soft_margin_mm, method.weight_floor);
    w_query = build_query_guard_soft_weight_local( ...
        x_sel, x_domain, query_guard_mm, method.domain_soft_margin_mm, method.weight_floor);
    w_gradient = build_template_gradient_weight_local(Tpl, x_sel, method.domain_selection_mode, method.weight_floor);
    w_total = max(method.weight_floor, w_edge .* w_template .* w_domain .* w_query .* w_gradient);
    if max(w_total) > 0
        w_total = max(method.weight_floor, w_total ./ max(w_total));
    end

    v_static = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_sel, 'pchip', NaN);
    shape_rmse = sqrt(mean((v_sel - v_static).^2, 'omitnan'));

    X = [X; x_sel]; %#ok<AGROW>
    T = [T; t_sel]; %#ok<AGROW>
    V = [V; v_sel]; %#ok<AGROW>
    S = [S; repmat(sid, numel(t_sel), 1)]; %#ok<AGROW>
    W = [W; w_total]; %#ok<AGROW>
    Theta = [Theta; theta_sel]; %#ok<AGROW>
    Y_obs = [Y_obs; y_obs_sel]; %#ok<AGROW>
    F0_all = [F0_all; f0_sel]; %#ok<AGROW>
    Fx_all = [Fx_all; fx_sel]; %#ok<AGROW>
    sensor_index_pt = [sensor_index_pt; repmat(is, numel(t_sel), 1)]; %#ok<AGROW>

    meta_idx = numel(window_meta) + 1;
    window_meta(meta_idx).sensor_id = sid;
    window_meta(meta_idx).valid_points = numel(t_sel);
    window_meta(meta_idx).pulse_segment_count = pulse_segment_count;
    window_meta(meta_idx).dynamic_effective_points = nnz(mask_effective);
    window_meta(meta_idx).query_safe_points = nnz(mask_query_safe(mask));
    window_meta(meta_idx).query_guard_mm = query_guard_mm;
    window_meta(meta_idx).domain_selection_mode = method.domain_selection_mode;
    window_meta(meta_idx).domain_soft_margin_mm = method.domain_soft_margin_mm;
    window_meta(meta_idx).template_domain = x_domain;
    window_meta(meta_idx).selected_x_range = [min(x_sel), max(x_sel)];
    window_meta(meta_idx).shape_rmse = shape_rmse;
end

if isempty(T)
    error('No valid template observation points were constructed for window %d.', Wmap.window_id);
end

interp_v = cell(numel(analysis_sensors), 1);
x_domain_by_sensor = NaN(numel(analysis_sensors), 2);
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    interp_v{is} = griddedInterpolant(Tpl.x_grid(:), Tpl.v_grid(:), 'pchip', 'none');
    x_domain_by_sensor(is, :) = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
end

bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.Y_obs = Y_obs;
bundle.F0 = F0_all;
bundle.Fx = Fx_all;
bundle.sensor_index = sensor_index_pt;
bundle.sensor_ids = analysis_sensors;
bundle.interp_v = interp_v;
bundle.x_domain_by_sensor = x_domain_by_sensor;
bundle.query_guard_mm = max([window_meta.query_guard_mm], [], 'omitnan');
bundle.overshoot_penalty_weight = method.overshoot_penalty_weight;
bundle.sensor_eta_limit_mm = method.sensor_eta_limit_mm;
bundle.sensor_eta_reg_weight_v_per_mm = method.sensor_eta_reg_weight_v_per_mm;
bundle.window_meta = window_meta;
bundle.valid_segment_count = sum([window_meta.pulse_segment_count]);
bundle.point_count = numel(T);
bundle.time_window = [min(T), max(T)];
bundle.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
bundle.rot_rpm_mean = Wmap.rot_rpm_mean;
end

function x_stat = invert_template_voltage_local(Tpl, v, x_ref)
x_grid = Tpl.x_grid(:);
v_grid = Tpl.v_grid(:);
x_stat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diff_v = v_grid - vv;
    crossing_x = [];
    exact_idx = find(abs(diff_v) <= 1e-10);
    if ~isempty(exact_idx)
        crossing_x = x_grid(exact_idx);
    end
    for k = 1:numel(diff_v)-1
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k+1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k+1) > 0
            continue;
        end
        denom = v_grid(k+1) - v_grid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - v_grid(k)) / denom;
        crossing_x(end+1, 1) = x_grid(k) + alpha * (x_grid(k+1) - x_grid(k)); %#ok<AGROW>
    end
    if isempty(crossing_x)
        [~, idx] = min(abs(diff_v));
        x_stat(i) = x_grid(idx);
    else
        [~, idx] = min(abs(crossing_x - x_ref(i)));
        x_stat(i) = crossing_x(idx);
    end
end
end

function query_guard_mm = resolve_query_guard_mm_local(Tpl, x, v, base_mask, method)
query_guard_mm = method.query_guard_mm;
if ~strcmpi(method.query_guard_mode, 'adaptive') || nnz(base_mask) < 8
    return;
end
x_sel = x(base_mask);
v_sel = v(base_mask);
x_stat = invert_template_voltage_local(Tpl, v_sel, x_sel);
u_app = abs(x_sel(:) - x_stat(:));
u_app = u_app(isfinite(u_app));
if isempty(u_app)
    return;
end
q = min(max(method.query_guard_quantile, 0), 100);
adaptive_guard = prctile(u_app, q) + method.query_guard_safety_mm;
query_guard_mm = min(max(adaptive_guard, method.query_guard_min_mm), method.query_guard_max_mm);
end

function eo_candidates = build_eo_candidates_local(rot_freq_hz, freq_search_hz, eo_pad)
if nargin < 3 || isempty(eo_pad)
    eo_pad = 0;
end
rot_freq_hz = max(rot_freq_hz, eps);
freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min_strict = max(1, ceil(freq_lo / rot_freq_hz));
eo_max_strict = max(eo_min_strict, floor(freq_hi / rot_freq_hz));
eo_candidates = eo_min_strict:eo_max_strict;
if isempty(eo_candidates)
    eo_center = max(1, round(mean([freq_lo, freq_hi]) / rot_freq_hz));
    eo_candidates = max(1, eo_center - eo_pad):max(1, eo_center + eo_pad);
end
freq_candidates = eo_candidates .* rot_freq_hz;
eo_candidates = eo_candidates(freq_candidates >= freq_lo & freq_candidates <= freq_hi);
if isempty(eo_candidates)
    error('No EO candidates remain inside [%.3f, %.3f] Hz at rot_freq %.6f Hz.', ...
        freq_lo, freq_hi, rot_freq_hz);
end
end

function seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates)
eo_candidates = unique(round(eo_candidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'weighted_voltage_rmse', inf, 'plain_voltage_rmse', inf, ...
    'linear_vp_rmse', inf), numel(eo_candidates), 1);
n_sensor = numel(bundle.sensor_ids);
for i = 1:numel(eo_candidates)
    eo = eo_candidates(i);
    s1 = sin(eo * bundle.Theta);
    c1 = cos(eo * bundle.Theta);
    basis = [-bundle.Fx(:), -bundle.Fx(:) .* s1, -bundle.Fx(:) .* c1];
    w_sqrt = sqrt(bundle.W(:));
    y = bundle.V(:) - bundle.F0(:);
    coeff = (basis .* w_sqrt) \ (y .* w_sqrt);
    dx_c = coeff(1);
    eta = zeros(1, n_sensor);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    lin_res = y - basis * coeff;
    [obj, plain_rmse] = template_synchronous_objective_local([A, phi, dx_c, eta], eo, bundle);
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx_c = dx_c;
    rows(i).sensor_eta = eta;
    rows(i).sensor_eta_max_abs = max(abs(eta), [], 'omitnan');
    rows(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    rows(i).plain_voltage_rmse = plain_rmse;
    rows(i).linear_vp_rmse = sqrt(sum(bundle.W(:) .* lin_res.^2) / max(sum(bundle.W(:)), eps));
end
score = [[rows.weighted_voltage_rmse].', [rows.linear_vp_rmse].', [rows.plain_voltage_rmse].', [rows.EO].'];
[~, order] = sortrows(score, [1, 2, 3]);
seed_table = rows(order);
end

function eo_keep = select_vp_topk_eo_local(seed_table, top_k)
if nargin < 2 || isempty(top_k) || ~isfinite(top_k)
    eo_keep = unique(round([seed_table.EO]), 'stable');
    return;
end
top_k = min(max(1, floor(top_k)), numel(seed_table));
eo_keep = unique(round([seed_table(1:top_k).EO]), 'stable');
end

function [eo_keep, info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, window_id)
all_eo = unique(round(eo_candidates(:).'));
top_k_eo = select_vp_topk_eo_local(seed_table, method.vp_top_k_eo);
eo_keep = intersect(all_eo, top_k_eo, 'stable');

best_rmse = seed_table(1).weighted_voltage_rmse;
if numel(seed_table) >= 2
    second_rmse = seed_table(2).weighted_voltage_rmse;
else
    second_rmse = inf;
end
vp_gap_ratio = (second_rmse - best_rmse) / max(best_rmse, eps);
linear_gap_ratio = inf;
if numel(seed_table) >= 2
    linear_gap_ratio = (seed_table(2).linear_vp_rmse - seed_table(1).linear_vp_rmse) / ...
        max(seed_table(1).linear_vp_rmse, eps);
end

fallback_reason = '';
if window_id <= cfg.vp_all_eo_warmup_windows
    fallback_reason = 'warmup_all_eo';
elseif vp_gap_ratio < cfg.vp_gap_ratio_fallback
    fallback_reason = 'ambiguous_weighted_vp_gap';
elseif linear_gap_ratio < cfg.vp_linear_gap_ratio_fallback
    fallback_reason = 'ambiguous_linear_vp_gap';
end

if ~isempty(fallback_reason)
    eo_keep = all_eo;
end

info = struct();
info.window_id = window_id;
info.all_eo = all_eo;
info.top_k_eo = top_k_eo;
info.selected_eo = eo_keep;
info.best_seed_eo = seed_table(1).EO;
info.best_seed_weighted_rmse = best_rmse;
info.second_seed_weighted_rmse = second_rmse;
info.vp_gap_ratio = vp_gap_ratio;
info.linear_gap_ratio = linear_gap_ratio;
info.fallback_reason = fallback_reason;
end

function result = refine_template_waveform_fit_local(bundle, seed_table, eo_candidates, cfg)
amp_limit_mm = abs(cfg.amplitude_limit_mm);
dx_c_limit_mm = abs(cfg.dx_c_limit_mm);
sensor_eta_limit_mm = abs(cfg.sensor_eta_limit_mm);
eo_candidates = unique(round(eo_candidates(:).'));
n_sensor = numel(bundle.sensor_ids);
candidate = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'fn_hz', NaN, 'objective', inf, 'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, 'V_pred', [], 'u_est', []), numel(eo_candidates), 1);
best = struct('objective', inf);

for i = 1:numel(eo_candidates)
    eo = eo_candidates(i);
    seed_idx = find([seed_table.EO] == eo, 1, 'first');
    if isempty(seed_idx)
        seed = seed_table(1);
    else
        seed = seed_table(seed_idx);
    end
    eta0 = zeros(1, n_sensor);
    if isfield(seed, 'sensor_eta') && ~isempty(seed.sensor_eta)
        eta0(1:min(n_sensor, numel(seed.sensor_eta))) = seed.sensor_eta(1:min(n_sensor, numel(seed.sensor_eta)));
    end
    seed_params = [seed.A, seed.phi, seed.dx_c, eta0];
    fun = @(p) template_synchronous_objective_local(p, eo, bundle);
    opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
    p_opt = fminsearch(@(p) bounded_synchronous_objective_local(p, fun, amp_limit_mm, dx_c_limit_mm, sensor_eta_limit_mm, n_sensor), seed_params, opts);
    p_opt(1) = min(abs(p_opt(1)), amp_limit_mm);
    p_opt(2) = wrap_to_pi_local(p_opt(2));
    p_opt(3) = max(min(p_opt(3), dx_c_limit_mm), -dx_c_limit_mm);
    p_opt = bound_sensor_eta_params_local(p_opt, sensor_eta_limit_mm, n_sensor);
    [obj, plain_rmse, v_pred, u_est, coverage] = template_synchronous_objective_local(p_opt, eo, bundle);

    candidate(i).EO = eo;
    candidate(i).A = p_opt(1);
    candidate(i).phi = p_opt(2);
    candidate(i).dx_c = p_opt(3);
    candidate(i).sensor_eta = p_opt(4:3+n_sensor);
    candidate(i).sensor_eta_max_abs = max(abs(candidate(i).sensor_eta), [], 'omitnan');
    candidate(i).fn_hz = eo * bundle.rot_freq_mean_hz;
    candidate(i).objective = obj;
    candidate(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    candidate(i).plain_voltage_rmse = plain_rmse;
    candidate(i).clamp_fraction = coverage.clamp_fraction;
    candidate(i).max_query_overshoot_mm = coverage.max_query_overshoot_mm;
    candidate(i).overshoot_penalty = coverage.overshoot_penalty;
    candidate(i).V_pred = v_pred;
    candidate(i).u_est = u_est;

    if obj < best.objective
        best = candidate(i);
    end
end

candidate_table = struct2table(rmfield(candidate, {'V_pred','u_est','sensor_eta'}));
candidate_table = sortrows(candidate_table, {'weighted_voltage_rmse','plain_voltage_rmse','EO'}, ...
    {'ascend','ascend','ascend'});

result = struct();
result.A_id = best.A;
result.EO_id = best.EO;
result.phi_id_wrapped = best.phi;
result.dx_c_id = best.dx_c;
result.d0_id = best.dx_c;
result.sensor_eta_id = best.sensor_eta;
result.sensor_eta_max_abs_mm = best.sensor_eta_max_abs;
result.fn_id = best.fn_hz;
result.weighted_voltage_rmse = best.weighted_voltage_rmse;
result.plain_voltage_rmse = best.plain_voltage_rmse;
result.V_pred = best.V_pred;
result.u_est = best.u_est;
result.objective_score = best.objective;
result.valid_segment_count = bundle.valid_segment_count;
result.point_count = bundle.point_count;
result.CandidateTable = candidate_table;
result.Coverage = struct( ...
    'clamp_fraction', best.clamp_fraction, ...
    'max_query_overshoot_mm', best.max_query_overshoot_mm, ...
    'query_guard_mm', bundle.query_guard_mm, ...
    'overshoot_penalty', best.overshoot_penalty, ...
    'overshoot_penalty_weight', bundle.overshoot_penalty_weight);
result.metric_note = 'Template-only VP top-K synchronous waveform fit with shared dx_c plus per-sensor x-zero eta; no EO prior.';
end

function obj = bounded_synchronous_objective_local(p, fun, amp_limit_mm, dx_c_limit_mm, sensor_eta_limit_mm, n_sensor)
p_use = p;
p_use(1) = min(abs(p_use(1)), amp_limit_mm);
p_use(2) = wrap_to_pi_local(p_use(2));
p_use(3) = max(min(p_use(3), dx_c_limit_mm), -dx_c_limit_mm);
p_use = bound_sensor_eta_params_local(p_use, sensor_eta_limit_mm, n_sensor);
obj = fun(p_use);
end

function p = bound_sensor_eta_params_local(p, sensor_eta_limit_mm, n_sensor)
need = 3 + n_sensor;
if numel(p) < need
    p(numel(p)+1:need) = 0;
end
p = p(1:need);
eta = p(4:end);
eta = max(min(eta, sensor_eta_limit_mm), -sensor_eta_limit_mm);
if ~isempty(eta)
    eta = eta - eta(1);
    eta = max(min(eta, sensor_eta_limit_mm), -sensor_eta_limit_mm);
end
p(4:end) = eta;
end

function [obj, plain_rmse, v_pred, u_est, coverage] = template_synchronous_objective_local(params, eo, bundle)
A = params(1);
phi = params(2);
dx_c = params(3);
params = bound_sensor_eta_params_local(params, abs(bundle.sensor_eta_limit_mm), numel(bundle.sensor_ids));
eta = params(4:end).';
u_est = A .* sin(eo .* bundle.Theta + phi);
eta_sample = eta(bundle.sensor_index(:));
x_in = bundle.X - dx_c - eta_sample - u_est;
v_pred = NaN(size(bundle.V));
clamped = false(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    x_lo = bundle.x_domain_by_sensor(is, 1);
    x_hi = bundle.x_domain_by_sensor(is, 2);
    x_eval_raw = x_in(mask);
    x_eval = min(max(x_eval_raw, x_lo), x_hi);
    clamped(mask) = x_eval_raw < x_lo | x_eval_raw > x_hi;
    overshoot(mask) = max(x_lo - x_eval_raw, 0) + max(x_eval_raw - x_hi, 0);
    v_pred(mask) = bundle.interp_v{is}(x_eval);
end
valid = isfinite(v_pred);
res = bundle.V(valid) - v_pred(valid);
overshoot_penalty = bundle.overshoot_penalty_weight * ...
    sum(bundle.W(valid) .* (overshoot(valid) .^ 2));
eta_reg_penalty = bundle.point_count * (bundle.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta .^ 2))) ^ 2;
obj = sum(bundle.W(valid) .* (res .^ 2)) + overshoot_penalty + eta_reg_penalty + 1e6 * nnz(~valid);
plain_rmse = sqrt(mean(res .^2));
coverage = struct();
coverage.clamp_fraction = nnz(clamped) / max(numel(clamped), 1);
coverage.max_query_overshoot_mm = max(overshoot, [], 'omitnan');
coverage.overshoot_penalty = overshoot_penalty;
coverage.sensor_eta_reg_penalty = eta_reg_penalty;
end

function summary = build_resonance_summary_local(trend_rows)
eo_vals = [trend_rows.EO_id];
freq_vals = [trend_rows.fn_id];
amp_vals = [trend_rows.A_id];
wrmse_vals = [trend_rows.weighted_voltage_rmse];
summary = struct();
summary.dominant_eo = mode(eo_vals(isfinite(eo_vals)));
summary.mean_freq_hz = mean(freq_vals, 'omitnan');
summary.median_freq_hz = median(freq_vals, 'omitnan');
summary.mean_amp_mm = mean(amp_vals, 'omitnan');
summary.median_amp_mm = median(amp_vals, 'omitnan');
summary.mean_weighted_rmse_v = mean(wrmse_vals, 'omitnan');
summary.median_weighted_rmse_v = median(wrmse_vals, 'omitnan');
end

function [mask, segment_count] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segment_count = numel(starts);
end

function [mask, segment_count] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segment_count = 1;
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
segment_count = numel(group_starts);
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

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function s = rank_to_string_local(rank_value)
if isfinite(rank_value)
    s = sprintf('%d', rank_value);
else
    s = 'NA';
end
end

