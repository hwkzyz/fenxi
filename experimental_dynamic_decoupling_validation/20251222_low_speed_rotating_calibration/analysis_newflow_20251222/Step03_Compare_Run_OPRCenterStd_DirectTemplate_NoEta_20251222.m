%% Step03_Compare_Run_OPRCenterStd_DirectTemplate_NoEta_20251222
% Baseline direct-template identification without eta_s.
% Tune parameters here, then run this file directly.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
C = CaseConfig();

%% Parameters to tune
P.name.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.name.dynamicSuffix = 'Main20L_W3S1_GradientXRange030_OPRCenterStd';
P.name.resultSuffix = 'Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe';

P.step03.analysisSensors = C.sensorIds;
P.step03.freqSearchHz = [100 1000];
P.step03.eoPad = 2;
P.step03.topKEO = 3;
P.step03.amplitudeLimitMM = 0.50;
P.step03.dxCLimitMM = 0.35;
P.step03.dynamicEffectiveMode = 'gradient';
P.step03.dynamicTemplateGradientMinRatio = 0.08;
P.step03.dynamicTimeGradientMinRatio = 0.15;
P.step03.dynamicPeakQuantile = 85;
P.step03.pulseMode = 'single';
P.step03.domainSelectionMode = 'hard';
P.step03.domainSoftMarginMM = 0;
P.step03.queryGuardMM = [];
P.step03.queryGuardMode = 'adaptive';
P.step03.queryGuardQuantile = 95;
P.step03.queryGuardSafetyMM = 0.05;
P.step03.queryGuardMinMM = 0.12;
P.step03.queryGuardMaxMM = [];
P.step03.phaseSafeExpansion = true;  % Re-add points whose fitted query position stays inside the template domain.
P.step03.phaseSafeMarginMM = 0.03;
P.step03.phaseSafeReferenceMode = 'prev_window';  % Use the previous accepted window for expansion; final fit still uses full EO top-K.
P.step03.phaseSafeFallbackMode = 'linear_vp';
P.step03.phaseSafeRefreshEvery = 5;
P.step03.phaseSafePrevMaxClampFraction = 1e-6;
P.step03.phaseSafePrevMinFinalGapRatio = 0.005;
P.step03.sensorEtaLimitMM = 0;       % Keep per-sensor eta fixed at 0 in the main result.
P.step03.sensorEtaAdaptiveLimitEnable = true;
P.step03.sensorEtaAdaptiveQuantile = 85;
P.step03.sensorEtaAdaptiveSafetyFactor = 1.20;
P.step03.sensorEtaAdaptiveMinMM = 0.02;
P.step03.sensorEtaAdaptiveMaxMM = 0.20;
P.step03.sensorEtaRegWeightVPerMM = 0;
P.step03.overshootPenaltyWeight = 100;
P.step03.debugMaxWindows = inf;
P.step03.allEoWarmupWindows = 0;
P.step03.vpGapRatioFallback = -inf;
P.step03.vpLinearGapRatioFallback = -inf;
P.step03.diagnosticEO = [];

%% Files
dynamicDir = fullfile(routeDir, 'output', 'dynamic_maps');
[templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20251222(routeDir, C, P.name.templateSuffix);
fprintf('Template source: %s [%s]\n', templateFile, templateSourceMode);

dynamicCandidates = {
    fullfile(dynamicDir, sprintf('DynamicMap_%s_SlidingWindows_%s_%s.mat', C.caseTag, P.name.dynamicSuffix, C.dataset))
    fullfile(dynamicDir, sprintf('DynamicMap_B%d_%s_SlidingWindows_%s_%s.mat', C.bladeId, C.sensorTag, P.name.dynamicSuffix, C.dataset))
    };
dynamicFile = '';
for i = 1:numel(dynamicCandidates)
    if exist(dynamicCandidates{i}, 'file') == 2
        dynamicFile = dynamicCandidates{i};
        break;
    end
end
if isempty(dynamicFile)
    dynamicFile = find_dynamic_file_for_sensors_local(dynamicDir, C.bladeId, C.sensorIds);
    fprintf('Dynamic map fallback: %s\n', dynamicFile);
end

fprintf('\n=== Step03 compare: direct-template identification without eta_s ===\n');
step03_identify_direct_template_embedded(routeDir, P, templateFile, dynamicFile, C);

function step03_identify_direct_template_embedded(rootDir, P, template_file, dynamic_file, C)
%% Step03_Compare_VPTop3SynchronousWaveform_NoEta_20251222
% Baseline low-speed-template-only identification route without eta_s.
% First-order template VP screens EO candidates; the final result still
% comes from the full synchronous waveform objective.
%
% Every integer EO candidate is scored by the first-order model
% V - T(x) ~= -T'(x) * [dx_c + a*sin(EO*theta) + b*cos(EO*theta)].
% By default, the top-3 EO candidates enter full waveform refinement in
% every window. For synchronous vibration, frequency is constrained by
% f = EO * rot_freq_mean; the final refinement optimizes amplitude/phase
% and the spatial alignment term dx_c.

%% Settings matching the super-Gaussian Step5 route
route_dir = rootDir;

cfg = struct();
cfg.target_blades = C.bladeId;
cfg.analysis_sensors = P.step03.analysisSensors;
cfg.analysis_start_time = 1.5;
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
cfg.freq_search_hz = P.step03.freqSearchHz;
cfg.eo_pad = P.step03.eoPad;
cfg.amplitude_limit_mm = P.step03.amplitudeLimitMM;
cfg.dx_c_limit_mm = P.step03.dxCLimitMM;
cfg.debug_max_windows = P.step03.debugMaxWindows;
cfg.vp_top_k_eo = P.step03.topKEO;
cfg.vp_all_eo_warmup_windows = P.step03.allEoWarmupWindows;
cfg.vp_gap_ratio_fallback = P.step03.vpGapRatioFallback;
cfg.vp_linear_gap_ratio_fallback = P.step03.vpLinearGapRatioFallback;
cfg.diagnostic_eo = P.step03.diagnosticEO;
cfg.sensor_eta_limit_mm = P.step03.sensorEtaLimitMM;
cfg.sensor_eta_adaptive_limit_enable = P.step03.sensorEtaAdaptiveLimitEnable;
cfg.sensor_eta_adaptive_quantile = P.step03.sensorEtaAdaptiveQuantile;
cfg.sensor_eta_adaptive_safety_factor = P.step03.sensorEtaAdaptiveSafetyFactor;
cfg.sensor_eta_adaptive_min_mm = P.step03.sensorEtaAdaptiveMinMM;
cfg.sensor_eta_adaptive_max_mm = P.step03.sensorEtaAdaptiveMaxMM;
cfg.sensor_eta_reg_weight_v_per_mm = P.step03.sensorEtaRegWeightVPerMM;
cfg.phase_safe_expansion = P.step03.phaseSafeExpansion;
cfg.phase_safe_margin_mm = P.step03.phaseSafeMarginMM;
cfg.phase_safe_reference_mode = lower(strtrim(P.step03.phaseSafeReferenceMode));
cfg.phase_safe_fallback_mode = lower(strtrim(P.step03.phaseSafeFallbackMode));
cfg.phase_safe_refresh_every = P.step03.phaseSafeRefreshEvery;
cfg.pulse_selection_mode = lower(strtrim(P.step03.pulseMode));
if ~ismember(cfg.pulse_selection_mode, {'single', 'all'})
    error('P.step03.pulseMode must be "single" or "all".');
end
if ~ismember(cfg.phase_safe_reference_mode, {'template_inverse', 'linear_vp', 'final_core', 'prev_window'})
    error('P.step03.phaseSafeReferenceMode must be "template_inverse", "linear_vp", "final_core", or "prev_window".');
end
if ~ismember(cfg.phase_safe_fallback_mode, {'template_inverse', 'linear_vp', 'final_core'})
    error('P.step03.phaseSafeFallbackMode must be "template_inverse", "linear_vp", or "final_core".');
end

method = struct();
method.name = 'template_only_main_vp_top3_synchronous_waveform';
method.selection_rule = 'first_order_template_vp_top3_then_synchronous_waveform_rmse';
method.use_reference_eo_constraint = false;
method.use_vp = true;
method.vp_top_k_eo = cfg.vp_top_k_eo;
if cfg.sensor_eta_limit_mm > 0
    method.refine_params = {'A', 'phi', 'dx_c', 'sensor_eta'};
else
    method.refine_params = {'A', 'phi', 'dx_c'};
end
method.template_forward = 'interp1_low_speed_template';
method.vp_seed_model = 'V_minus_Tx_equals_minus_Tprime_times_u';
method.final_waveform_objective = 'fixed_eo_direct_template_voltage_residual_without_sensor_affine_projection';
method.store_y_obs = false;  % Current first-order VP route does not use template-inverse displacement observations.
method.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
method.sensor_eta_adaptive_limit_enable = cfg.sensor_eta_adaptive_limit_enable;
method.sensor_eta_reg_weight_v_per_mm = cfg.sensor_eta_reg_weight_v_per_mm;
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.query_guard_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm + method.coverage_safety_margin_mm;
if ~isempty(P.step03.queryGuardMM)
    method.query_guard_mm = P.step03.queryGuardMM;
end
method.query_guard_mode = lower(strtrim(P.step03.queryGuardMode));
if ~ismember(method.query_guard_mode, {'fixed', 'adaptive'})
    error('P.step03.queryGuardMode must be "fixed" or "adaptive".');
end
method.query_guard_quantile = P.step03.queryGuardQuantile;
method.query_guard_safety_mm = P.step03.queryGuardSafetyMM;
method.query_guard_min_mm = P.step03.queryGuardMinMM;
method.query_guard_max_mm = method.query_guard_mm;
if ~isempty(P.step03.queryGuardMaxMM)
    method.query_guard_max_mm = P.step03.queryGuardMaxMM;
end
method.phase_safe_expansion_enabled = logical(cfg.phase_safe_expansion);
method.phase_safe_margin_mm = max(0, cfg.phase_safe_margin_mm);
method.phase_safe_reference_mode = cfg.phase_safe_reference_mode;
method.phase_safe_fallback_mode = cfg.phase_safe_fallback_mode;
method.phase_safe_refresh_every = max(0, floor(cfg.phase_safe_refresh_every));
method.phase_safe_prev_max_clamp_fraction = P.step03.phaseSafePrevMaxClampFraction;
method.phase_safe_prev_min_final_gap_ratio = P.step03.phaseSafePrevMinFinalGapRatio;
method.phase_safe_u_limit_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm;
method.domain_selection_mode = lower(strtrim(P.step03.domainSelectionMode));
if ~ismember(method.domain_selection_mode, {'hard', 'soft'})
    error('P.step03.domainSelectionMode must be "hard" or "soft".');
end
method.domain_soft_margin_mm = P.step03.domainSoftMarginMM;
method.overshoot_penalty_weight = P.step03.overshootPenaltyWeight;
method.coverage_mode = 'prefer_points_safe_for_all_bounded_u';
method.extrapolation_mode = 'clamp_to_template_edge';
method.default_sensor_threshold = 0.5;
method.pulse_selection_mode = cfg.pulse_selection_mode;
method.dynamic_effective_mode = lower(strtrim(P.step03.dynamicEffectiveMode));
if ~ismember(method.dynamic_effective_mode, {'legacy', 'gradient'})
    error('P.step03.dynamicEffectiveMode must be "legacy" or "gradient".');
end
method.dynamic_template_gradient_min_ratio = P.step03.dynamicTemplateGradientMinRatio;
method.dynamic_time_gradient_min_ratio = P.step03.dynamicTimeGradientMinRatio;
method.dynamic_peak_quantile = P.step03.dynamicPeakQuantile;

sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
result_dir = fullfile(output_dir, 'identification');
if exist(result_dir, 'dir') ~= 7; mkdir(result_dir); end

if exist(template_file, 'file') ~= 2
    template_file = find_template_file_for_sensors_local(template_dir, cfg.target_blades, cfg.analysis_sensors);
end
if exist(dynamic_file, 'file') ~= 2
    dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, cfg.target_blades, cfg.analysis_sensors);
end
result_suffix = sanitize_result_suffix_local(P.name.resultSuffix);
result_file = fullfile(result_dir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_%s%s_%s.mat', ...
    C.caseTag, result_suffix, C.dataset));

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

fprintf('\n=== Step03 compare: VP top-3 synchronous waveform identification without eta_s ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Template source: %s\n', template_file);
fprintf('Dynamic map source: %s\n', dynamic_file);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinate_check.max_abs_delta_mm, coordinate_check.status);
legacy_base_xc_gradient_tpl = contains(template_file, 'GradientXRange030') && ...
    ~contains(dynamic_file, 'GradientXRange030');
if ~coordinate_check.is_consistent && legacy_base_xc_gradient_tpl
    fprintf('Coordinate note: legacy 20251222 route uses base-xc DynamicMap with GradientXRange030 template.\n');
elseif ~coordinate_check.is_consistent
    warning('Template/DynamicMap x-center mismatch. This run is not a clean GradientXRange030 coordinate chain.');
end
fprintf('Frequency search: %.1f-%.1f Hz; no reference EO; first-order VP only screens candidates.\n', ...
    cfg.freq_search_hz(1), cfg.freq_search_hz(2));
fprintf('Model: V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi)); synchronous f = EO*rot_freq_mean.\n');
fprintf(['Template coverage guard: %s, fixed %.3f mm, min %.3f mm, max %.3f mm, ' ...
    'q%.1f + %.3f mm.\n'], method.query_guard_mode, method.query_guard_mm, ...
    method.query_guard_min_mm, method.query_guard_max_mm, method.query_guard_quantile, ...
    method.query_guard_safety_mm);
if method.phase_safe_expansion_enabled
    fprintf('Phase-safe expansion enabled: margin %.3f mm, reference %s; final expanded fit keeps EO top-K.\n', ...
        method.phase_safe_margin_mm, method.phase_safe_reference_mode);
    if strcmpi(method.phase_safe_reference_mode, 'prev_window')
        fprintf('Prev-window phase-safe: fallback %s, refresh every %d window(s), max clamp %.3g, min final gap %.3g.\n', ...
            method.phase_safe_fallback_mode, method.phase_safe_refresh_every, ...
            method.phase_safe_prev_max_clamp_fraction, method.phase_safe_prev_min_final_gap_ratio);
    end
else
    fprintf('Phase-safe expansion disabled.\n');
end
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
    'valid_segment_count', NaN, 'point_count', NaN, ...
    'core_point_count', NaN, 'added_point_count', NaN, ...
    'phase_safe_reference_source', '', 'phase_safe_reference_reason', ''), num_windows, 1);
window_results = repmat(empty_window_result_local(), num_windows, 1);
best_window = struct('weighted_voltage_rmse', inf);
fprintf('Window execution: serial full top-K flow; phase-safe does not reduce EO screening.\n');

prev_phase_result = [];
for w_idx = 1:num_windows
    [trend_rows(w_idx), window_results(w_idx)] = process_step03_window_local( ...
        DynamicMap.Window(w_idx), Template, cfg, method, w_idx, prev_phase_result);
    result = window_results(w_idx).Result;
    seed_table = window_results(w_idx).seed_table;
    Wmap = DynamicMap.Window(w_idx);
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
    prev_phase_result = result;
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

tmp_result_file = sprintf('%s.tmp_%s.mat', result_file(1:end-4), datestr(now, 'yyyymmddTHHMMSSFFF'));
save(tmp_result_file, 'Result', '-v7');
if exist(result_file, 'file') == 2
    delete(result_file);
end
movefile(tmp_result_file, result_file, 'f');

fprintf('\n=== Step03 main trend ===\n');
disp(Trend(:, {'window_id','lap_start','lap_end','EO_id','fn_id','A_id','weighted_voltage_rmse'}));
fprintf('Dominant EO: %d; mean f = %.6f Hz; median f = %.6f Hz\n', ...
    Result.ResonanceSummary.dominant_eo, Result.ResonanceSummary.mean_freq_hz, ...
    Result.ResonanceSummary.median_freq_hz);
fprintf('Saved result: %s\n', result_file);

%% Local functions
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

function template_file = resolve_template_input_file_local(template_dir, C, template_suffix)
preferred = fullfile(template_dir, ...
    sprintf('Template_LowSpeedRotating_%s_%s_%s.mat', C.caseTag, template_suffix, C.dataset));
if exist(preferred, 'file') == 2
    template_file = preferred;
    return;
end
template_file = find_template_file_for_sensors_local(template_dir, C.bladeId, C.sensorIds);
fprintf('Template file fallback: %s\n', template_file);
end

function dynamic_file = resolve_dynamic_input_file_local(dynamic_dir, C, dynamic_suffix)
preferred = fullfile(dynamic_dir, ...
    sprintf('DynamicMap_%s_SlidingWindows_%s_%s.mat', C.caseTag, dynamic_suffix, C.dataset));
if exist(preferred, 'file') == 2
    dynamic_file = preferred;
    return;
end
dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, C.bladeId, C.sensorIds);
fprintf('Dynamic map fallback: %s\n', dynamic_file);
end

function template_file = find_template_file_for_sensors_local(template_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
case_tag = sprintf('B%d_%s', target_blade, sensor_tag);
patterns = {
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_OPRCenterStd_20251222.mat', case_tag)
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_20251222.mat', case_tag)
    sprintf('Template_LowSpeedRotating_%s_20251222.mat', case_tag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_20251222.mat', target_blade, sensor_tag)
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
files = dir(fullfile(template_dir, sprintf('Template_LowSpeedRotating_%s*.mat', case_tag)));
if isempty(files)
    files = dir(fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_*.mat', target_blade)));
end
if isempty(files)
    error('No low-speed template file found under %s.', template_dir);
end
error('Template for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
case_tag = sprintf('B%d_%s', target_blade, sensor_tag);
patterns = {
    sprintf('DynamicMap_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', case_tag)
    sprintf('DynamicMap_%s_SlidingWindows_20251222.mat', case_tag)
    sprintf('DynamicMap_%s_SlidingWindows_GradientXRange030_20251222.mat', case_tag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat', target_blade, sensor_tag)
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
files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_%s*.mat', case_tag)));
if isempty(files)
    files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_B%d_*.mat', target_blade)));
end
if isempty(files)
    error('No dynamic map file found under %s.', dynamic_dir);
end
error('DynamicMap for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function [trend_row, window_result] = process_step03_window_local(Wmap, Template, cfg, method, w_idx, prev_phase_result)
if nargin < 6
    prev_phase_result = [];
end
bundle = build_template_observation_bundle_local(Wmap, Template, cfg.analysis_sensors, method);
eo_candidates = build_eo_candidates_local(bundle.rot_freq_mean_hz, cfg.freq_search_hz, cfg.eo_pad);
core_bundle = bundle;
phase_ref_info = empty_phase_safe_reference_info_local();

if method.phase_safe_expansion_enabled
    if strcmpi(method.phase_safe_reference_mode, 'template_inverse')
        phase_ref = build_template_inverse_phase_safe_reference_local();
        phase_ref_info = make_phase_safe_reference_info_local( ...
            phase_ref.reference_source, 'configured_template_inverse', NaN, false);
        expanded_bundle = build_template_observation_bundle_local( ...
            Wmap, Template, cfg.analysis_sensors, method, phase_ref);
        if expanded_bundle.point_count > core_bundle.point_count
            bundle = expanded_bundle;
        end
        seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
        [final_eo_candidates, selection_info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
        result = refine_template_waveform_fit_local(bundle, seed_table, final_eo_candidates, cfg);
        result.CorePass = summarize_phase_safe_reference_local(phase_ref, core_bundle, method);
        result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
            core_bundle, bundle, phase_ref, method);
    else
        if strcmpi(method.phase_safe_reference_mode, 'prev_window')
            [phase_ref, phase_ref_info] = build_prev_window_phase_safe_reference_local( ...
                prev_phase_result, w_idx, method, cfg);
            if isempty(phase_ref)
                [phase_ref, phase_ref_info] = build_current_window_phase_safe_reference_local( ...
                    bundle, eo_candidates, cfg, method, w_idx, phase_ref_info.reason);
            end
        else
            [phase_ref, phase_ref_info] = build_current_window_phase_safe_reference_local( ...
                bundle, eo_candidates, cfg, method, w_idx, 'configured_current_window_reference');
        end
        expanded_bundle = build_template_observation_bundle_local( ...
            Wmap, Template, cfg.analysis_sensors, method, phase_ref);
        if expanded_bundle.point_count > core_bundle.point_count
            expanded_seed_table = solve_template_seed_eo_scan_local(expanded_bundle, eo_candidates);
            [expanded_eo_candidates, expanded_selection_info] = select_vp_adaptive_eo_local( ...
                expanded_seed_table, eo_candidates, method, cfg, w_idx);
            result = refine_template_waveform_fit_local( ...
                expanded_bundle, expanded_seed_table, expanded_eo_candidates, cfg);
            result.CorePass = summarize_phase_safe_reference_local(phase_ref, core_bundle, method);
            result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
                core_bundle, expanded_bundle, phase_ref, method);
            seed_table = expanded_seed_table;
            final_eo_candidates = expanded_eo_candidates;
            selection_info = expanded_selection_info;
            bundle = expanded_bundle;
        else
            seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
            [final_eo_candidates, selection_info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
            result = refine_template_waveform_fit_local(bundle, seed_table, final_eo_candidates, cfg);
            result.CorePass = summarize_phase_safe_reference_local(phase_ref, core_bundle, method);
            result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
                core_bundle, expanded_bundle, phase_ref, method);
        end
    end
else
    seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
    [final_eo_candidates, selection_info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
    result = refine_template_waveform_fit_local(bundle, seed_table, final_eo_candidates, cfg);
end
result.window_id = w_idx;
result.lap_range = Wmap.lap_range;
result.VPSeedTable = struct2table(seed_table);
result.VPSelectedEO = final_eo_candidates(:).';
result.VPSelectionInfo = selection_info;
result.PhaseSafeReferenceInfo = phase_ref_info;

trend_row = struct( ...
    'window_id', w_idx, ...
    'lap_start', Wmap.lap_range(1), ...
    'lap_end', Wmap.lap_range(end), ...
    'window_center_time', median(bundle.T, 'omitnan'), ...
    'rot_freq_mean_hz', bundle.rot_freq_mean_hz, ...
    'rot_rpm_mean', Wmap.rot_rpm_mean, ...
    'A_id', result.A_id, ...
    'EO_id', result.EO_id, ...
    'fn_id', result.fn_id, ...
    'phi_id_wrapped', result.phi_id_wrapped, ...
    'dx_c_id', result.dx_c_id, ...
    'd0_id', result.dx_c_id, ...
    'sensor_eta_max_abs_mm', max(abs(result.sensor_eta_id), [], 'omitnan'), ...
    'weighted_voltage_rmse', result.weighted_voltage_rmse, ...
    'plain_voltage_rmse', result.plain_voltage_rmse, ...
    'valid_segment_count', result.valid_segment_count, ...
    'point_count', result.point_count, ...
    'core_point_count', core_bundle.point_count, ...
    'added_point_count', result.point_count - core_bundle.point_count, ...
    'phase_safe_reference_source', phase_ref_info.source, ...
    'phase_safe_reference_reason', phase_ref_info.reason);

window_result = struct();
window_result.window_id = w_idx;
window_result.lap_range = Wmap.lap_range;
window_result.time_window = Wmap.time_window;
window_result.bundle = bundle;
window_result.seed_table = seed_table;
window_result.CandidateTable = result.CandidateTable;
window_result.Result = result;
end

function phase_ref = build_template_inverse_phase_safe_reference_local()
phase_ref = struct();
phase_ref.EO_id = NaN;
phase_ref.A_id = NaN;
phase_ref.phi_id_wrapped = NaN;
phase_ref.dx_c_id = NaN;
phase_ref.sensor_eta_id = [];
phase_ref.weighted_voltage_rmse = NaN;
phase_ref.plain_voltage_rmse = NaN;
phase_ref.reference_source = 'template_inverse_apparent_query';
end

function phase_ref = build_phase_safe_reference_from_seed_local(seed_table, cfg)
seed = seed_table(1);
phase_ref = struct();
phase_ref.EO_id = seed.EO;
phase_ref.A_id = min(abs(seed.A), abs(cfg.amplitude_limit_mm));
phase_ref.phi_id_wrapped = wrap_to_pi_local(seed.phi);
phase_ref.dx_c_id = max(min(seed.dx_c, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phase_ref.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
phase_ref.weighted_voltage_rmse = seed.weighted_voltage_rmse;
phase_ref.plain_voltage_rmse = seed.plain_voltage_rmse;
phase_ref.reference_source = 'linear_vp_seed';
end

function info = empty_phase_safe_reference_info_local()
info = make_phase_safe_reference_info_local('', '', NaN, false);
end

function info = make_phase_safe_reference_info_local(source, reason, source_window_id, used_previous_window)
info = struct();
info.source = source;
info.reason = reason;
info.source_window_id = source_window_id;
info.used_previous_window = used_previous_window;
end

function [phase_ref, info] = build_current_window_phase_safe_reference_local( ...
    bundle, eo_candidates, cfg, method, w_idx, fallback_reason)
mode = method.phase_safe_reference_mode;
if strcmpi(method.phase_safe_reference_mode, 'prev_window')
    mode = method.phase_safe_fallback_mode;
end

if strcmpi(mode, 'template_inverse')
    phase_ref = build_template_inverse_phase_safe_reference_local();
    reason = sprintf('%s_then_%s', fallback_reason, mode);
elseif strcmpi(mode, 'linear_vp')
    seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
    phase_ref = build_phase_safe_reference_from_seed_local(seed_table, cfg);
    reason = sprintf('%s_then_%s', fallback_reason, mode);
elseif strcmpi(mode, 'final_core')
    seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
    [core_eo_candidates, ~] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
    phase_ref = refine_template_waveform_fit_local(bundle, seed_table, core_eo_candidates, cfg);
    phase_ref.reference_source = 'current_core_final_waveform';
    reason = sprintf('%s_then_%s', fallback_reason, mode);
else
    error('Unsupported phase-safe current-window reference mode: %s.', mode);
end
info = make_phase_safe_reference_info_local(phase_ref.reference_source, reason, w_idx, false);
end

function [phase_ref, info] = build_prev_window_phase_safe_reference_local(prev_result, w_idx, method, cfg)
phase_ref = [];
[ok, reason] = is_prev_window_phase_reference_usable_local(prev_result, w_idx, method);
if ~ok
    info = make_phase_safe_reference_info_local('', reason, NaN, false);
    return;
end

phase_ref = struct();
phase_ref.EO_id = prev_result.EO_id;
phase_ref.A_id = min(abs(prev_result.A_id), abs(cfg.amplitude_limit_mm));
phase_ref.phi_id_wrapped = wrap_to_pi_local(prev_result.phi_id_wrapped);
phase_ref.dx_c_id = max(min(prev_result.dx_c_id, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phase_ref.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
if isfield(prev_result, 'sensor_eta_id') && ~isempty(prev_result.sensor_eta_id)
    n = min(numel(cfg.analysis_sensors), numel(prev_result.sensor_eta_id));
    phase_ref.sensor_eta_id(1:n) = prev_result.sensor_eta_id(1:n);
end
phase_ref.weighted_voltage_rmse = prev_result.weighted_voltage_rmse;
phase_ref.plain_voltage_rmse = prev_result.plain_voltage_rmse;
phase_ref.reference_source = 'prev_window_final';
phase_ref.source_window_id = prev_result.window_id;
info = make_phase_safe_reference_info_local(phase_ref.reference_source, 'prev_window_quality_pass', ...
    prev_result.window_id, true);
end

function [ok, reason] = is_prev_window_phase_reference_usable_local(prev_result, w_idx, method)
ok = false;
if isempty(prev_result) || ~isstruct(prev_result)
    reason = 'no_previous_window';
    return;
end
if method.phase_safe_refresh_every > 0 && mod(w_idx - 1, method.phase_safe_refresh_every) == 0
    reason = 'scheduled_refresh';
    return;
end
required = {'window_id', 'EO_id', 'A_id', 'phi_id_wrapped', 'dx_c_id', 'weighted_voltage_rmse', 'Coverage'};
for i = 1:numel(required)
    if ~isfield(prev_result, required{i})
        reason = ['previous_missing_', required{i}];
        return;
    end
end
vals = [prev_result.window_id, prev_result.EO_id, prev_result.A_id, ...
    prev_result.phi_id_wrapped, prev_result.dx_c_id, prev_result.weighted_voltage_rmse];
if any(~isfinite(vals))
    reason = 'previous_nonfinite_reference';
    return;
end
if isfield(prev_result.Coverage, 'clamp_fraction') && ...
        prev_result.Coverage.clamp_fraction > method.phase_safe_prev_max_clamp_fraction
    reason = 'previous_clamped';
    return;
end
if isfield(prev_result, 'VPSelectionInfo') && isfield(prev_result.VPSelectionInfo, 'fallback_reason') && ...
        ~isempty(prev_result.VPSelectionInfo.fallback_reason)
    reason = 'previous_vp_fallback';
    return;
end
gap_ratio = final_candidate_gap_ratio_local(prev_result);
if isfinite(method.phase_safe_prev_min_final_gap_ratio) && ...
        gap_ratio < method.phase_safe_prev_min_final_gap_ratio
    reason = 'previous_final_gap_too_small';
    return;
end
ok = true;
reason = 'prev_window_quality_pass';
end

function gap_ratio = final_candidate_gap_ratio_local(result)
gap_ratio = inf;
if ~isfield(result, 'CandidateTable') || height(result.CandidateTable) < 2 || ...
        ~ismember('weighted_voltage_rmse', result.CandidateTable.Properties.VariableNames)
    return;
end
rmse = result.CandidateTable.weighted_voltage_rmse;
rmse = rmse(isfinite(rmse));
if numel(rmse) < 2
    return;
end
rmse = sort(rmse(:), 'ascend');
gap_ratio = (rmse(2) - rmse(1)) / max(rmse(1), eps);
end

function window_result = empty_window_result_local()
window_result = struct( ...
    'window_id', NaN, ...
    'lap_range', [], ...
    'time_window', [], ...
    'bundle', [], ...
    'seed_table', [], ...
    'CandidateTable', table(), ...
    'Result', []);
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

function bundle = build_template_observation_bundle_local(Wmap, Template, analysis_sensors, method, phase_ref)
if nargin < 5
    phase_ref = [];
end
use_phase_safe = ~isempty(phase_ref) && method.phase_safe_expansion_enabled;
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
    mask_phase_safe = false(size(mask_query_safe));
    if use_phase_safe
        mask_phase_safe = build_phase_safe_query_mask_local( ...
            Tpl, x_raw, v_raw, theta_raw, x_domain, is, phase_ref, method);
        if nnz(mask_base & mask_phase_safe) >= 8
            mask_query_safe = mask_phase_safe;
        end
    end
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
    if method.store_y_obs
        y_obs_sel = x_sel - invert_template_voltage_local(Tpl, v_sel, x_sel);
    else
        y_obs_sel = [];
    end

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
    [eta_initial_guess_mm, eta_initial_iqr_mm] = estimate_sensor_eta_scale_local(Tpl, x_sel, v_sel);

    X = [X; x_sel]; %#ok<AGROW>
    T = [T; t_sel]; %#ok<AGROW>
    V = [V; v_sel]; %#ok<AGROW>
    S = [S; repmat(sid, numel(t_sel), 1)]; %#ok<AGROW>
    W = [W; w_total]; %#ok<AGROW>
    Theta = [Theta; theta_sel]; %#ok<AGROW>
    if method.store_y_obs
        Y_obs = [Y_obs; y_obs_sel]; %#ok<AGROW>
    end
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
    window_meta(meta_idx).phase_safe_expansion = use_phase_safe;
    window_meta(meta_idx).phase_safe_points = nnz(mask_base & mask_phase_safe);
    window_meta(meta_idx).template_domain = x_domain;
    window_meta(meta_idx).selected_x_range = [min(x_sel), max(x_sel)];
    window_meta(meta_idx).shape_rmse = shape_rmse;
    window_meta(meta_idx).eta_initial_guess_mm = eta_initial_guess_mm;
    window_meta(meta_idx).eta_initial_iqr_mm = eta_initial_iqr_mm;
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
bundle.sensor_eta_limit_by_sensor = build_sensor_eta_limit_by_sensor_local(window_meta, cfg);
bundle.sensor_eta_reg_weight_v_per_mm = method.sensor_eta_reg_weight_v_per_mm;
bundle.window_meta = window_meta;
bundle.valid_segment_count = sum([window_meta.pulse_segment_count]);
bundle.point_count = numel(T);
if use_phase_safe
    bundle.selection_pass = 'phase_safe_expanded';
else
    bundle.selection_pass = 'core_hard_adaptive';
end
bundle.phase_safe_margin_mm = method.phase_safe_margin_mm;
bundle.time_window = [min(T), max(T)];
bundle.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
bundle.rot_rpm_mean = Wmap.rot_rpm_mean;
end

function mask_safe = build_phase_safe_query_mask_local(Tpl, x, v, theta, x_domain, sensor_local_index, phase_ref, method)
if isfield(phase_ref, 'reference_source') && strcmpi(phase_ref.reference_source, 'template_inverse_apparent_query')
    mask_safe = build_template_inverse_safe_mask_local(Tpl, x, v, x_domain, method);
    return;
end
A = phase_ref.A_id;
eo = phase_ref.EO_id;
phi = phase_ref.phi_id_wrapped;
dx_c = phase_ref.dx_c_id;
eta = 0;
if isfield(phase_ref, 'sensor_eta_id') && numel(phase_ref.sensor_eta_id) >= sensor_local_index
    eta = phase_ref.sensor_eta_id(sensor_local_index);
end
u_est = A .* sin(eo .* theta(:) + phi);
x_query = x(:) - dx_c - eta - u_est;
margin = method.phase_safe_margin_mm;
mask_safe = x_query >= x_domain(1) + margin & x_query <= x_domain(2) - margin;
mask_safe = reshape(mask_safe, size(x));
end

function mask_safe = build_template_inverse_safe_mask_local(Tpl, x, v, x_domain, method)
x_stat = invert_template_voltage_local(Tpl, v(:), x(:));
margin = method.phase_safe_margin_mm;
u_app = abs(x(:) - x_stat(:));
mask_safe = isfinite(x_stat) & ...
    x_stat >= x_domain(1) + margin & x_stat <= x_domain(2) - margin & ...
    u_app <= method.phase_safe_u_limit_mm + margin;
mask_safe = reshape(mask_safe, size(x));
end

function summary = summarize_phase_safe_reference_local(phase_ref, bundle, method)
summary = struct();
summary.EO_id = phase_ref.EO_id;
summary.A_id = phase_ref.A_id;
summary.dx_c_id = phase_ref.dx_c_id;
summary.phi_id_wrapped = phase_ref.phi_id_wrapped;
summary.weighted_voltage_rmse = phase_ref.weighted_voltage_rmse;
summary.plain_voltage_rmse = phase_ref.plain_voltage_rmse;
summary.point_count = bundle.point_count;
summary.selection_pass = bundle.selection_pass;
summary.reference_mode = method.phase_safe_reference_mode;
summary.reference_source = get_optional_field_local(phase_ref, 'reference_source', method.phase_safe_reference_mode);
end

function summary = summarize_phase_safe_expansion_local(core_bundle, expanded_bundle, phase_ref, method)
summary = struct();
summary.enabled = method.phase_safe_expansion_enabled;
summary.reference_EO_id = phase_ref.EO_id;
summary.reference_A_id = phase_ref.A_id;
summary.reference_dx_c_id = phase_ref.dx_c_id;
summary.reference_phi_id_wrapped = phase_ref.phi_id_wrapped;
summary.reference_mode = method.phase_safe_reference_mode;
summary.margin_mm = method.phase_safe_margin_mm;
summary.core_point_count = core_bundle.point_count;
summary.expanded_point_count = expanded_bundle.point_count;
summary.added_point_count = expanded_bundle.point_count - core_bundle.point_count;
summary.added_fraction_of_core = summary.added_point_count / max(core_bundle.point_count, 1);
end

function value = get_optional_field_local(S, field_name, default_value)
if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
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

function [eta_initial_guess_mm, eta_initial_iqr_mm] = estimate_sensor_eta_scale_local(Tpl, x_sel, v_sel)
eta_initial_guess_mm = NaN;
eta_initial_iqr_mm = NaN;
if isempty(x_sel) || isempty(v_sel)
    return;
end
keep = isfinite(x_sel(:)) & isfinite(v_sel(:));
if nnz(keep) < 5
    return;
end
x_use = x_sel(keep);
v_use = v_sel(keep);
x_stat = invert_template_voltage_local(Tpl, v_use, x_use);
eta_app = x_use - x_stat;
eta_app = eta_app(isfinite(eta_app));
if isempty(eta_app)
    return;
end
eta_initial_guess_mm = median(eta_app, 'omitnan');
eta_initial_iqr_mm = iqr(eta_app);
end

function eta_limit_by_sensor = build_sensor_eta_limit_by_sensor_local(window_meta, cfg)
n_sensor = numel(window_meta);
eta_limit_by_sensor = abs(cfg.sensor_eta_limit_mm) * ones(1, n_sensor);
if ~isfield(cfg, 'sensor_eta_limit_mm') || abs(cfg.sensor_eta_limit_mm) <= 0
    return;
end
if ~isfield(cfg, 'sensor_eta_adaptive_limit_enable') || ~cfg.sensor_eta_adaptive_limit_enable
    return;
end
q = min(max(cfg.sensor_eta_adaptive_quantile, 0), 100);
scale = max(cfg.sensor_eta_adaptive_safety_factor, 1);
min_mm = max(cfg.sensor_eta_adaptive_min_mm, 0);
max_mm = max(abs(cfg.sensor_eta_limit_mm), abs(cfg.sensor_eta_adaptive_max_mm));
for is = 1:n_sensor
    eta_iqr = NaN;
    if isfield(window_meta(is), 'eta_initial_iqr_mm')
        eta_iqr = window_meta(is).eta_initial_iqr_mm;
    end
    if isfinite(eta_iqr) && eta_iqr > 0
        eta_lim = scale * abs(eta_iqr);
    else
        eta_lim = abs(cfg.sensor_eta_limit_mm);
    end
    eta_limit_by_sensor(is) = min(max(eta_lim, min_mm), max_mm);
end
if ~isempty(eta_limit_by_sensor)
    eta_limit_by_sensor(1) = min(max(eta_limit_by_sensor(1), min_mm), max_mm);
end
end

function sensor_eta_limit_mm = resolve_sensor_eta_limit_vector_local(bundle, cfg)
if isfield(bundle, 'sensor_eta_limit_by_sensor') && ~isempty(bundle.sensor_eta_limit_by_sensor)
    sensor_eta_limit_mm = abs(bundle.sensor_eta_limit_by_sensor(:).');
else
    sensor_eta_limit_mm = abs(cfg.sensor_eta_limit_mm) * ones(1, numel(bundle.sensor_ids));
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
sensor_eta_limit_mm = resolve_sensor_eta_limit_vector_local(bundle, cfg);
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
sensor_eta_limit_mm = abs(sensor_eta_limit_mm(:).');
if isempty(sensor_eta_limit_mm)
    sensor_eta_limit_mm = inf(1, n_sensor);
elseif numel(sensor_eta_limit_mm) == 1
    sensor_eta_limit_mm = repmat(sensor_eta_limit_mm, 1, n_sensor);
elseif numel(sensor_eta_limit_mm) < n_sensor
    sensor_eta_limit_mm(numel(sensor_eta_limit_mm)+1:n_sensor) = sensor_eta_limit_mm(end);
else
    sensor_eta_limit_mm = sensor_eta_limit_mm(1:n_sensor);
end
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
params = bound_sensor_eta_params_local(params, resolve_sensor_eta_limit_vector_local(bundle, struct('sensor_eta_limit_mm', bundle.sensor_eta_limit_mm)), numel(bundle.sensor_ids));
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
end
