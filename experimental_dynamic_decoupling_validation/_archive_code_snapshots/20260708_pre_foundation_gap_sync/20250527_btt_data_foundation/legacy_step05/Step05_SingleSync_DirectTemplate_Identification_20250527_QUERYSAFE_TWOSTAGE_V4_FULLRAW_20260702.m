clc; clear; close all;

%STEP05_SINGLESYNC_DIRECTTEMPLATE_IDENTIFICATION_20250527
% VERSION: QUERYSAFE_TWOSTAGE_V4_FULLRAW_20260702
%
% Key logic:
%   1. High-speed waveform extraction follows the original OPRCenterStd
%      DynamicMap route: fixed window around target blade passage time.
%   2. No within-pulse uniform downsampling is applied; all raw samples in
%      the fixed time window are mapped before x-domain masks.
%   3. No dynamic-gradient / time-gradient / peak-quantile mask is used as
%      the main fitting-point selector.
%   4. Core bundle uses a query-safe shrunken region:
%          core_mask = finite & inside_domain & inside_query_guard
%   5. Optional expanded bundle uses phase-safe x_query:
%          x_query = x - dx_c - eta_s - A*sin(EO*theta + phi)
%      and keeps points whose x_query stays inside the template domain.
%   6. dynamic_effective_mask and main_pulse_mask are saved only for
%      diagnostics unless S.core_mask_mode='legacy_base_query_safe'.
%
% Forward model:
%   V_i ~= T_{sensor,blade}(x_i - dx_c - eta_s - A*sin(EO*theta_i + phi))
%
% Usage:
%   Run this script directly after Step02 and Step04.

cfg = BTTProjectConfig_20250527();

%% Step05 local settings
S = struct();
S.version = 'QUERYSAFE_TWOSTAGE_V4_FULLRAW_20260702';

S.show_plots = true;
S.save_figures = cfg.save_figures;
S.target_cases = {'20250526_2500-3500_t400'};
S.target_blades = 1;
S.analysis_sensors = [1 3 6];
S.reference_sensor_id = S.analysis_sensors(1);

S.analysis_start_time_s = 1.5;
S.target_laps = 20;
S.window_laps = 3;
S.sliding_step_laps = 1;
S.debug_max_windows = inf;

S.output_dir = fullfile(cfg.output_root, 'step05_single_sync_direct_template');
S.figure_dir = fullfile(cfg.figure_root, 'step05_single_sync_direct_template');

run_tag_env = strtrim(getenv('STEP05_OUTPUT_TAG'));
if ~isempty(run_tag_env)
    run_tag_env = regexprep(run_tag_env, '[^\w\-]', '_');
    S.output_dir = fullfile(cfg.output_root, ...
        ['step05_single_sync_direct_template_', run_tag_env]);
    S.figure_dir = fullfile(cfg.figure_root, ...
        ['step05_single_sync_direct_template_', run_tag_env]);
end

show_plots_env = strtrim(getenv('STEP05_SHOW_PLOTS'));
if ~isempty(show_plots_env)
    S.show_plots = any(strcmpi(show_plots_env, {'1', 'true', 'yes', 'on'}));
end

save_figures_env = strtrim(getenv('STEP05_SAVE_FIGURES'));
if ~isempty(save_figures_env)
    S.save_figures = any(strcmpi(save_figures_env, {'1', 'true', 'yes', 'on'}));
end

debug_max_windows_env = strtrim(getenv('STEP05_DEBUG_MAX_WINDOWS'));
if ~isempty(debug_max_windows_env)
    parsed_debug_windows = str2double(debug_max_windows_env);
    if isfinite(parsed_debug_windows) && parsed_debug_windows > 0
        S.debug_max_windows = parsed_debug_windows;
    end
end

%% Step05 identification settings
S.freq_search_hz = [100 1000];
S.eo_pad = 2;
S.top_k_eo = 3;

S.amplitude_limit_mm = 0.50;
S.dx_c_limit_mm = 0.35;
S.sensor_eta_limit_mm = 0.20;
S.sensor_eta_reg_weight_v_per_mm = 0.02;

S.dynamic_window_mode = 'peak_centered_fixed';
S.pulse_window_sec = 6e-4;
S.pulse_pad_fraction = 0.30;
S.pulse_min_pad_points = 50;

S.core_mask_mode = 'query_safe_full_wave';
S.domain_margin_mm = 0.02;

S.query_guard_mode = 'adaptive';
S.query_guard_mm = [];
S.query_guard_min_mm = 0.12;
S.query_guard_max_mm = S.amplitude_limit_mm + S.dx_c_limit_mm + S.sensor_eta_limit_mm + 0.10;
S.query_guard_safety_mm = 0.05;
S.query_guard_quantile = 95;

S.phase_safe_expansion = true;
S.phase_safe_margin_mm = 0.03;
S.phase_safe_reference_mode = 'core_current';
S.choose_better_pass = true;
S.better_pass_rmse_tolerance_v = 0;

%% Step05 diagnostic, gate, and optimizer settings
S.pulse_mode = 'single';
S.dynamic_template_gradient_min_ratio = 0.08;
S.dynamic_time_gradient_min_ratio = 0.15;
S.dynamic_peak_quantile = 85;

S.min_fit_points = 80;
S.min_valid_pass_count = 2;
S.weight_floor = 0.05;

S.fminsearch_max_iter = 800;
S.fminsearch_max_fun = 2000;
S.overshoot_penalty_weight = 100;
S.objective_mode = 'rmse';
S.method_preset = 'vp_top3_bounded_rmse';
S.store_bundle_preview_points = 6000;

method_preset_env = strtrim(getenv('STEP05_METHOD_PRESET'));
if ~isempty(method_preset_env)
    S.method_preset = lower(method_preset_env);
    switch S.method_preset
        case 'vp_top3_bounded_rmse'
            S.top_k_eo = 3;
            S.amplitude_limit_mm = 0.50;
            S.dx_c_limit_mm = 0.35;
            S.sensor_eta_limit_mm = 0.20;
            S.sensor_eta_reg_weight_v_per_mm = 0.02;
            S.overshoot_penalty_weight = 100;
            S.objective_mode = 'rmse';
        case 'all_eo_unconstrained_rmse'
            S.top_k_eo = inf;
            S.amplitude_limit_mm = inf;
            S.dx_c_limit_mm = inf;
            S.sensor_eta_limit_mm = inf;
            S.sensor_eta_reg_weight_v_per_mm = 0;
            S.overshoot_penalty_weight = 0;
            S.objective_mode = 'rmse';
        case 'all_eo_unconstrained_corr'
            S.top_k_eo = inf;
            S.amplitude_limit_mm = inf;
            S.dx_c_limit_mm = inf;
            S.sensor_eta_limit_mm = inf;
            S.sensor_eta_reg_weight_v_per_mm = 0;
            S.overshoot_penalty_weight = 0;
            S.objective_mode = 'corr';
        otherwise
            error('Unsupported STEP05_METHOD_PRESET: %s', method_preset_env);
    end
    if isfinite(S.amplitude_limit_mm) && isfinite(S.dx_c_limit_mm) && ...
            isfinite(S.sensor_eta_limit_mm)
        S.query_guard_max_mm = S.amplitude_limit_mm + ...
            S.dx_c_limit_mm + S.sensor_eta_limit_mm + 0.10;
    else
        S.query_guard_max_mm = inf;
    end
end

if ischar(S.target_cases) || isstring(S.target_cases)
    S.target_cases = cellstr(S.target_cases);
end
target_cases = S.target_cases;

fprintf('\n=== Step05: single-sync direct-template identification ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Target cases: %s\n', strjoin(cellstr(target_cases), ', '));
fprintf('Target blades: %s\n', mat2str(S.target_blades));
fprintf('Analysis sensors: %s\n', mat2str(S.analysis_sensors));
fprintf('Model: V = T_{s,b}(x - dx_c - eta_s - A*sin(EO*theta + phi)).\n');
fprintf('Eta gauge: eta of reference sensor CH%d is fixed to zero.\n', S.analysis_sensors(1));
fprintf('Core mask mode: %s\n', S.core_mask_mode);
fprintf('Phase-safe expansion: %d\n', S.phase_safe_expansion);

Template = load_step04_template_local(cfg);

template_override_file = strtrim(getenv('STEP05_TEMPLATE_OVERRIDE_FILE'));
if ~isempty(template_override_file)
    if ~isfile(template_override_file)
        error('STEP05_TEMPLATE_OVERRIDE_FILE does not exist: %s', ...
            template_override_file);
    end
    override_loaded = load(template_override_file, 'Template');
    if ~isfield(override_loaded, 'Template')
        error('Template override file must contain variable Template: %s', ...
            template_override_file);
    end
    Template = override_loaded.Template;
    fprintf('Step05 template override: %s\n', template_override_file);
end

ResultSet = struct();
ResultSet.Dataset = cfg.dataset;
ResultSet.CreatedBy = mfilename;
ResultSet.CreatedOn = datestr(now, 31);
ResultSet.RunInfo = build_step05_run_info_local(cfg, S);
ResultSet.Cases = repmat(struct('case_name', '', 'BladeResult', []), numel(target_cases), 1);

for iCase = 1:numel(target_cases)
    case_name = char(target_cases{iCase});
    fprintf('\n--- Step05 case: %s ---\n', case_name);

    case_result = process_step05_case_local(case_name, cfg, S, Template);

    ResultSet.Cases(iCase).case_name = case_name;
    ResultSet.Cases(iCase).BladeResult = case_result;
end


function Template = load_step04_template_local(cfg)
if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end

sensor_tag = ['S', sprintf('%d', cfg.sensor_ids)];

preferred = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', ...
    sensor_tag, cfg.dataset));

if isfile(preferred)
    loaded = load(preferred, 'Template');
    Template = loaded.Template;
    fprintf('Step04 template: %s\n', preferred);
    return;
end

files = dir(fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_*_%s.mat', cfg.dataset)));

if isempty(files)
    error('No Step04 template file found under %s. Please run Step04_Build_LowSpeed_OPRCenterStd_Template_20250527.m first.', ...
        cfg.step04_output_dir);
end

if isempty(files)
    error('No Step04 template file found under %s.', cfg.step04_output_dir);
end

[~, idx] = max([files.datenum]);
file = fullfile(files(idx).folder, files(idx).name);

loaded = load(file, 'Template');
Template = loaded.Template;

fprintf('Step04 template fallback: %s\n', file);
end


function blade_results = process_step05_case_local(case_name, cfg, S, Template)
step02_case_dir = fullfile(cfg.step02_output_dir, case_name);

if ~isfolder(step02_case_dir) || ...
        ~isfile(fullfile(step02_case_dir, 'Step02_Dynamic_BTT_Extraction_20250527.mat'))

    error('Step02 products missing for %s. Please run Step02_Extract_Dynamic_BTT_20250527.m first.', case_name);
end

loaded = load(fullfile(step02_case_dir, ...
    'Step02_Dynamic_BTT_Extraction_20250527.mat'), ...
    'case_data', 'metadata');

case_data = loaded.case_data;

if isfield(loaded, 'metadata')
    step02_metadata = loaded.metadata;
else
    step02_metadata = struct();
end

case_dir = fullfile(cfg.dataset_root, case_name);

if ~isfolder(case_dir)
    error('Dynamic raw data folder not found: %s', case_dir);
end

opr_times = load_dynamic_opr_times_local(step02_case_dir);

[F_omega_deg_s, speed_diag] = build_dynamic_speed_interpolant_local( ...
    opr_times, cfg.blades_num);

blade_results = repmat(struct('blade_id', NaN, 'Result', []), ...
    numel(S.target_blades), 1);

for ib = 1:numel(S.target_blades)
    blade_id = S.target_blades(ib);

    fprintf('\n>>> [Step05][%s] Target blade %d\n', case_name, blade_id);

    check_step05_template_coverage_local(Template, S, blade_id);

    passes_by_sensor = build_analysis_passes_by_sensor_local( ...
        case_data, cfg, S, blade_id);

    ref_key = sprintf('CH%d', S.reference_sensor_id);

    if ~isfield(passes_by_sensor, ref_key)
        error('Reference sensor CH%d has no target blade passes.', ...
            S.reference_sensor_id);
    end

    ref_passes = passes_by_sensor.(ref_key);

    n_laps = min(S.target_laps, height(ref_passes));

    if n_laps < S.window_laps
        error('Not enough target passes: n_laps=%d, window_laps=%d.', ...
            n_laps, S.window_laps);
    end

    ref_passes = ref_passes(1:n_laps, :);

    WindowPlan = build_window_plan_local(ref_passes, S);

    if isfinite(S.debug_max_windows)
        WindowPlan = WindowPlan( ...
            1:min(height(WindowPlan), S.debug_max_windows), :);
    end

    fprintf('Windows: %d; lap window=%d, step=%d; start time %.3f s.\n', ...
        height(WindowPlan), S.window_laps, ...
        S.sliding_step_laps, S.analysis_start_time_s);

    Result = run_single_blade_step05_local( ...
        case_name, case_dir, case_data, step02_metadata, ...
        opr_times, F_omega_deg_s, speed_diag, ...
        Template, passes_by_sensor, WindowPlan, cfg, S, blade_id);

    out_dir = fullfile(S.output_dir, case_name);
    fig_case_dir = fullfile(S.figure_dir, case_name);

    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end

    if exist(fig_case_dir, 'dir') ~= 7
        mkdir(fig_case_dir);
    end

    sensor_tag = ['S', sprintf('%d', S.analysis_sensors)];

    result_file = fullfile(out_dir, sprintf( ...
        'Result_Step05_SingleSyncDirectTemplate_B%d_%s_%s.mat', ...
        blade_id, sensor_tag, cfg.dataset));

    trend_file = fullfile(out_dir, sprintf( ...
        'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
        blade_id, sensor_tag, cfg.dataset));

    save(result_file, 'Result', '-v7.3');
    writetable(Result.Trend, trend_file);

    fprintf('Saved Step05 result:\n  %s\n', result_file);
    fprintf('Saved Step05 trend:\n  %s\n', trend_file);

    if S.show_plots || S.save_figures
        plot_step05_result_local(Result, S, fig_case_dir, S.show_plots);
    end

    blade_results(ib).blade_id = blade_id;
    blade_results(ib).Result = Result;
end
end


function Result = run_single_blade_step05_local( ...
    case_name, case_dir, case_data, step02_metadata, ...
    opr_times, F_omega_deg_s, speed_diag, ...
    Template, passes_by_sensor, WindowPlan, cfg, S, blade_id)

method = build_step05_method_metadata_local(S);

trend_rows = repmat(make_empty_trend_row_local(), height(WindowPlan), 1);
WindowResult = repmat(make_empty_window_result_local(), height(WindowPlan), 1);

best = struct('weighted_voltage_rmse', inf);

for iw = 1:height(WindowPlan)
    W = WindowPlan(iw, :);

    % ---------------------------------------------------------------------
    % Candidate waveform map.
    % Full high-speed waveform extracted from fixed windows and mapped into
    % OPRCenterStd coordinates. No within-pulse downsampling and no
    % gradient/peak/main-pulse trimming here.
    % ---------------------------------------------------------------------
    bundle_all = build_direct_bundle_for_window_local( ...
        case_dir, case_data, opr_times, F_omega_deg_s, ...
        Template, passes_by_sensor, cfg, S, blade_id, W);

    bundle_all = attach_template_mini_records_local( ...
        bundle_all, Template, S.analysis_sensors, blade_id);

    WindowResult(iw).CandidateBundlePreview = ...
        downsample_bundle_for_storage_local(bundle_all, ...
        max(S.store_bundle_preview_points, 9000));

    if bundle_all.point_count < S.min_fit_points
        warning('Window %d has too few candidate points (%d).', ...
            iw, bundle_all.point_count);

        [trend_rows(iw), WindowResult(iw)] = make_failed_window_result_local( ...
            W, iw, bundle_all, 'too_few_candidate_points');
        continue;
    end

    eo_candidates = build_eo_candidates_local( ...
        bundle_all.rot_freq_mean_hz, ...
        S.freq_search_hz, ...
        S.eo_pad);

    % ---------------------------------------------------------------------
    % First pass: query-safe core bundle.
    % This is the conservative region smaller than Step04 x_domain.
    % ---------------------------------------------------------------------
    core_bundle = select_bundle_points_local(bundle_all, bundle_all.core_mask, S);

    if core_bundle.point_count < S.min_fit_points
        warning('Window %d has too few query-safe core points (%d).', ...
            iw, core_bundle.point_count);

        [trend_rows(iw), WindowResult(iw)] = make_failed_window_result_local( ...
            W, iw, bundle_all, 'too_few_query_safe_core_points');
        continue;
    end

    if core_bundle.valid_pass_count < S.min_valid_pass_count
        warning('Window %d has too few valid core passes (%d).', ...
            iw, core_bundle.valid_pass_count);

        [trend_rows(iw), WindowResult(iw)] = make_failed_window_result_local( ...
            W, iw, bundle_all, 'too_few_valid_core_passes');
        continue;
    end

    seed_table_core = solve_vp_seed_eo_scan_local(core_bundle, eo_candidates, S);
    final_eos_core = select_top_eos_local(seed_table_core, S.top_k_eo);

    fit_core = refine_direct_template_fit_local( ...
        core_bundle, seed_table_core, final_eos_core, S);

    fit_core.fit_stage = 'core_query_safe';
    fit_core.core_point_count = core_bundle.point_count;
    fit_core.expanded_point_count = core_bundle.point_count;

    % ---------------------------------------------------------------------
    % Second pass: deterministic phase-safe expansion.
    % Expanded points are selected by fitted x_query, not by fallback.
    % ---------------------------------------------------------------------
    expanded_bundle = [];
    seed_table_expanded = [];
    fit_expanded = [];

    expand_info = struct();
    expand_info.applied = false;
    expand_info.reason = 'disabled_or_not_applicable';
    expand_info.core_point_count = core_bundle.point_count;
    expand_info.expanded_point_count = NaN;
    expand_info.chosen_pass = 'core_query_safe';

    if S.phase_safe_expansion && strcmpi(fit_core.status, 'ok')
        phase_ref = fit_core;

        expanded_mask = build_phase_safe_expansion_mask_local( ...
            bundle_all, phase_ref, S);

        expanded_bundle = select_bundle_points_local(bundle_all, expanded_mask, S);

        expand_info.applied = true;
        expand_info.expanded_point_count = expanded_bundle.point_count;
        expand_info.reason = 'phase_safe_from_core_current';

        if expanded_bundle.point_count >= S.min_fit_points && ...
                expanded_bundle.valid_pass_count >= S.min_valid_pass_count

            seed_table_expanded = solve_vp_seed_eo_scan_local( ...
                expanded_bundle, eo_candidates, S);

            final_eos_expanded = select_top_eos_local( ...
                seed_table_expanded, S.top_k_eo);

            fit_expanded = refine_direct_template_fit_local( ...
                expanded_bundle, seed_table_expanded, final_eos_expanded, S);

            fit_expanded.fit_stage = 'phase_safe_expanded';
            fit_expanded.core_point_count = core_bundle.point_count;
            fit_expanded.expanded_point_count = expanded_bundle.point_count;
        else
            expand_info.reason = 'phase_safe_too_few_points';
        end
    end

    % ---------------------------------------------------------------------
    % Explicit pass choice.
    % No hidden fallback. Chosen pass is recorded.
    % ---------------------------------------------------------------------
    [final_bundle, result, final_seed_table, expand_info] = choose_step05_pass_local( ...
        core_bundle, fit_core, seed_table_core, ...
        expanded_bundle, fit_expanded, seed_table_expanded, ...
        expand_info, S);

    result.window_id = iw;
    result.lap_range = [W.lap_start, W.lap_end];
    result.time_window_s = [W.time_start_s, W.time_end_s];
    result.VPSeedTable = struct2table(final_seed_table);
    result.BundleSummary = summarize_bundle_local(final_bundle, bundle_all);
    result.ReconstructionPreview = build_reconstruction_preview_local( ...
        final_bundle, result, S.store_bundle_preview_points);

    trend_rows(iw) = build_trend_row_local(iw, W, result, final_bundle);

    WindowResult(iw).window_id = iw;
    WindowResult(iw).lap_range = [W.lap_start, W.lap_end];
    WindowResult(iw).time_window_s = [W.time_start_s, W.time_end_s];
    WindowResult(iw).Result = result;
    WindowResult(iw).CoreResult = fit_core;
    WindowResult(iw).ExpandedResult = fit_expanded;
    WindowResult(iw).ExpansionInfo = expand_info;

    WindowResult(iw).BundlePreview = ...
        downsample_bundle_for_storage_local(final_bundle, ...
        S.store_bundle_preview_points);

    WindowResult(iw).CoreBundlePreview = ...
        downsample_bundle_for_storage_local(core_bundle, ...
        S.store_bundle_preview_points);

    if ~isempty(expanded_bundle)
        WindowResult(iw).ExpandedBundlePreview = ...
            downsample_bundle_for_storage_local(expanded_bundle, ...
            S.store_bundle_preview_points);
    end

    WindowResult(iw).CorePointCount = core_bundle.point_count;
    WindowResult(iw).CandidatePointCount = bundle_all.point_count;
    WindowResult(iw).SeedTable = final_seed_table;

    if strcmpi(result.status, 'ok') && ...
            result.weighted_voltage_rmse < best.weighted_voltage_rmse

        best = result;
        best.window_id = iw;
        best.lap_range = [W.lap_start, W.lap_end];
        best.time_window_s = [W.time_start_s, W.time_end_s];
    end

    fprintf(['Window %02d/%02d laps [%d %d]: chosen=%s, EO=%d, f=%.3f Hz, ' ...
        'A=%.4f mm, dx=%.4f mm, RMSE=%.5f V, points core=%d, final=%d, candidate=%d\n'], ...
        iw, height(WindowPlan), W.lap_start, W.lap_end, result.fit_stage, ...
        result.EO_id, result.fn_id, result.A_id, result.dx_c_id, ...
        result.weighted_voltage_rmse, core_bundle.point_count, ...
        final_bundle.point_count, bundle_all.point_count);
end

Trend = struct2table(trend_rows);

Result = struct();
Result.Dataset = cfg.dataset;
Result.Case_Name = case_name;
Result.CreatedBy = mfilename;
Result.CreatedOn = datestr(now, 31);
Result.TargetBlade = blade_id;
Result.SensorIDs = S.analysis_sensors;
Result.ReferenceSensorID = S.reference_sensor_id;
Result.Method = method;
Result.AnalysisSettings = S;
Result.RunInfo = build_step05_run_info_local(cfg, S);
Result.Step02_Metadata = step02_metadata;
Result.Step04_Template_Metadata = get_optional_struct_local(Template, 'Metadata');
Result.SpeedDiagnostic = speed_diag;
Result.WindowPlan = WindowPlan;
Result.Trend = Trend;
Result.WindowResult = WindowResult;
Result.BestWindow = best;
Result.ResonanceSummary = build_resonance_summary_local(Trend);
end


function RunInfo = build_step05_run_info_local(cfg, S)
RunInfo = struct();
RunInfo.step_name = mfilename;
RunInfo.version = S.version;
RunInfo.created_on = datestr(now, 31);
RunInfo.project_config = cfg;
RunInfo.step_settings = S;
end


function method = build_step05_method_metadata_local(S)
method = struct();
method.name = 'Step05_single_sync_direct_template_OPRCenterStd_querysafe_twostage_fullraw';
method.version = 'QUERYSAFE_TWOSTAGE_V4_FULLRAW_20260702';
method.outer_workflow = 'Step5_SingleSync style sliding-window single-synchronous identification';
method.solver = 'OPRCenterStd low-speed non-parametric direct-template waveform inversion';
method.forward_model = 'V = T_{s,b}(x - dx_c - eta_s - A*sin(EO*theta + phi))';
method.vp_seed_model = 'V - T(x) = -Tprime(x) * [dx_c + a*sin(EO*theta) + b*cos(EO*theta) + eta_s]';
method.eo_selection = sprintf('top-%d VP candidates then full waveform RMSE', S.top_k_eo);
method.eta_gauge = sprintf('eta_s fixed to zero for CH%d', S.analysis_sensors(1));
method.no_template_center_subtraction = true;
method.no_sensor_affine_voltage_projection = true;
method.no_within_pulse_uniform_downsampling = true;
method.core_mask_mode = S.core_mask_mode;
method.phase_safe_expansion = S.phase_safe_expansion;
method.two_stage_region_policy = ['full raw candidate waveform -> query-safe core region -> ' ...
    'optional phase-safe expanded region'];
end


function check_step05_template_coverage_local(Template, S, blade_id)
for sid = S.analysis_sensors
    tpl = get_template_entry_local(Template, sid, blade_id);

    if isempty(tpl) || isempty(tpl.x_grid)
        error('Step04 template missing for CH%d Blade%d.', sid, blade_id);
    end

    if any(~isfinite(tpl.x_domain)) || tpl.x_domain(2) <= tpl.x_domain(1)
        warning('Step04 template CH%d Blade%d has invalid x_domain.', sid, blade_id);
    end
end
end


function [F_omega_deg_s, diagnostic] = build_dynamic_speed_interpolant_local(opr_times, blades_num)
opr_times = opr_times(:);

rev_period = opr_times((blades_num + 1):end) - ...
             opr_times(1:(end - blades_num));

speed_time = opr_times(1:(end - blades_num));
speed_deg_s = 360 ./ max(rev_period, eps);

valid = isfinite(speed_time) & isfinite(speed_deg_s) & speed_deg_s > 0;

if nnz(valid) < 2
    error('Not enough valid OPR periods for dynamic speed interpolation.');
end

F_omega_deg_s = griddedInterpolant( ...
    speed_time(valid), speed_deg_s(valid), ...
    'linear', 'nearest');

diagnostic = struct();
diagnostic.speed_time_s = speed_time(valid);
diagnostic.speed_deg_s = speed_deg_s(valid);
diagnostic.rpm = speed_deg_s(valid) / 360 * 60;
end


function passes_by_sensor = build_analysis_passes_by_sensor_local(case_data, cfg, S, blade_id)
passes_by_sensor = struct();

for sid = S.analysis_sensors
    if sid > numel(case_data.channels) || isempty(case_data.channels(sid).jilublade)
        warning('No jilublade data for CH%d.', sid);
        continue;
    end

    rows = case_data.channels(sid).jilublade;

    mask = isfinite(rows(:,4)) & rows(:,4) == blade_id;
    rows = rows(mask, :);
    rows = rows(isfinite(rows(:,3)), :);

    [~, order] = sort(rows(:,3));
    rows = rows(order, :);

    % Use the same target-blade passes after the same analysis start time for
    % every sensor. Lap indexing must be reset after this filter.
    rows = rows(rows(:,3) >= S.analysis_start_time_s, :);

    if isempty(rows)
        warning('No Blade%d passes for CH%d after %.6f s.', ...
            blade_id, sid, S.analysis_start_time_s);
        continue;
    end

    if size(rows, 1) > S.target_laps
        rows = rows(1:S.target_laps, :);
    end

    analysis_lap = (1:size(rows,1)).';

    T = table( ...
        analysis_lap, ...
        rows(:,1), ...
        rows(:,2), ...
        rows(:,3), ...
        rows(:,4), ...
        'VariableNames', {'analysis_lap','start_s','end_s','arrival_s','blade_id'});

    passes_by_sensor.(sprintf('CH%d', sid)) = T;
end
end


function WindowPlan = build_window_plan_local(ref_passes, S)
n_laps = height(ref_passes);

wins = [];
window_id = 0;

for lap_start = 1:S.sliding_step_laps:(n_laps - S.window_laps + 1)
    lap_end = lap_start + S.window_laps - 1;
    window_id = window_id + 1;

    time_start = ref_passes.start_s(lap_start);
    time_end = ref_passes.end_s(lap_end);
    time_center = median(ref_passes.arrival_s(lap_start:lap_end), 'omitnan');

    wins = [wins; window_id, lap_start, lap_end, time_start, time_end, time_center]; %#ok<AGROW>
end

WindowPlan = array2table(wins, ...
    'VariableNames', {'window_id','lap_start','lap_end','time_start_s','time_end_s','time_center_s'});
end


function bundle = build_direct_bundle_for_window_local( ...
    case_dir, case_data, opr_times, F_omega_deg_s, ...
    Template, passes_by_sensor, cfg, S, blade_id, W)

bundle = make_empty_bundle_local();

bundle.case_name = case_data.case_name;
bundle.blade_id = blade_id;
bundle.lap_range = [W.lap_start, W.lap_end];
bundle.time_window_s = [W.time_start_s, W.time_end_s];

rot_freq_samples = [];
valid_pass_count = 0;
pass_id_global = 0;

for sid = S.analysis_sensors
    key = sprintf('CH%d', sid);

    if ~isfield(passes_by_sensor, key)
        continue;
    end

    Tpass = passes_by_sensor.(key);

    pass_rows = Tpass( ...
        Tpass.analysis_lap >= W.lap_start & ...
        Tpass.analysis_lap <= W.lap_end, :);

    if isempty(pass_rows)
        continue;
    end

    tpl = get_template_entry_local(Template, sid, blade_id);

    if isempty(tpl)
        continue;
    end

    for ip = 1:height(pass_rows)
        pass_id_global = pass_id_global + 1;
        p = pass_rows(ip, :);

        % -------------------------------------------------------------
        % 1. Candidate waveform extraction.
        %    Follow original OPRCenterStd DynamicMap route:
        %    fixed window around target-blade passage time.
        %    All raw samples in this fixed time window are retained.
        %    No within-pulse downsampling, gradient / peak / query-safe
        %    trimming is performed before OPRCenterStd mapping.
        % -------------------------------------------------------------
        [t_seg, v_seg] = read_dynamic_pulse_segment_local( ...
            case_dir, case_data, sid, ...
            p.start_s, p.end_s, p.arrival_s, cfg, S);

        if numel(t_seg) < 5
            continue;
        end

        % -------------------------------------------------------------
        % 2. OPRCenterStd mapping.
        %    Every sample in one pulse segment uses the OPR immediately
        %    before the blade passage time as angular reference.
        % -------------------------------------------------------------
        theta_std_deg = Template.Standard_Relative_Angles(sid, blade_id);

        [x_mm, theta_rad, prev_opr_idx, keep] = ...
            map_dynamic_times_to_oprcenterstd_local( ...
            t_seg, opr_times, F_omega_deg_s, ...
            theta_std_deg, cfg, p.arrival_s);

        t_seg = t_seg(keep);
        v_seg = v_seg(keep);
        x_mm = x_mm(keep);
        theta_rad = theta_rad(keep);
        prev_opr_idx = prev_opr_idx(keep);

        if numel(x_mm) < 5
            continue;
        end

        % -------------------------------------------------------------
        % 3. Template evaluation at raw x.
        % -------------------------------------------------------------
        [template_v, template_dv_dx, template_weight, inside_domain_raw] = ...
            evaluate_template_at_points_local(tpl, x_mm, S);

        inside_domain = x_mm >= tpl.x_domain(1) + S.domain_margin_mm & ...
                        x_mm <= tpl.x_domain(2) - S.domain_margin_mm;

        inside_domain = inside_domain & inside_domain_raw;

        finite_mask = isfinite(x_mm) & ...
                      isfinite(v_seg) & ...
                      isfinite(theta_rad) & ...
                      isfinite(template_v) & ...
                      isfinite(template_dv_dx);

        % -------------------------------------------------------------
        % 4. Diagnostic masks only.
        %    These are retained for checking, not for the default core mask.
        % -------------------------------------------------------------
        dVdt = dynamic_time_gradient_local(t_seg, v_seg);

        dynamic_effective = build_dynamic_effective_mask_local( ...
            v_seg, dVdt, template_dv_dx, tpl, S);

        main_pulse_mask = build_main_pulse_mask_local( ...
            t_seg, v_seg, S);

        % -------------------------------------------------------------
        % 5. Query-safe shrink.
        %    Since final query coordinate is:
        %       x_query = x - dx_c - eta_s - A*sin(...)
        %    the first/core pass must use a smaller region than x_domain.
        % -------------------------------------------------------------
        base_for_guard = finite_mask & inside_domain;

        query_guard = compute_query_guard_local( ...
            x_mm, v_seg, tpl, S, base_for_guard);

        inside_guard = x_mm >= tpl.x_domain(1) + query_guard & ...
                       x_mm <= tpl.x_domain(2) - query_guard;

        switch lower(strtrim(S.core_mask_mode))
            case 'legacy_base_query_safe'
                % Historical base mask from the original direct-template
                % route, plus query-safe shrink.
                base_mask = finite_mask & ...
                            inside_domain & ...
                            dynamic_effective & ...
                            main_pulse_mask;

                core_mask = base_mask & inside_guard;

            otherwise
                % Main deterministic route for the paper:
                % complete mapped high-speed waveform inside the query-safe
                % shrunken trusted domain.
                core_mask = finite_mask & inside_domain & inside_guard;
        end

        n = numel(x_mm);

        bundle.x = [bundle.x; x_mm(:)]; %#ok<AGROW>
        bundle.v = [bundle.v; v_seg(:)]; %#ok<AGROW>
        bundle.t = [bundle.t; t_seg(:)]; %#ok<AGROW>
        bundle.theta = [bundle.theta; theta_rad(:)]; %#ok<AGROW>
        bundle.sensor_id = [bundle.sensor_id; repmat(sid, n, 1)]; %#ok<AGROW>
        bundle.blade_id_vec = [bundle.blade_id_vec; repmat(blade_id, n, 1)]; %#ok<AGROW>
        bundle.analysis_lap = [bundle.analysis_lap; repmat(p.analysis_lap, n, 1)]; %#ok<AGROW>
        bundle.pass_id = [bundle.pass_id; repmat(pass_id_global, n, 1)]; %#ok<AGROW>
        bundle.prev_opr_index = [bundle.prev_opr_index; prev_opr_idx(:)]; %#ok<AGROW>

        bundle.template_v = [bundle.template_v; template_v(:)]; %#ok<AGROW>
        bundle.template_dv_dx = [bundle.template_dv_dx; template_dv_dx(:)]; %#ok<AGROW>
        bundle.template_weight = [bundle.template_weight; template_weight(:)]; %#ok<AGROW>

        bundle.inside_domain_mask = [bundle.inside_domain_mask; inside_domain(:)]; %#ok<AGROW>
        bundle.inside_guard_mask = [bundle.inside_guard_mask; inside_guard(:)]; %#ok<AGROW>
        bundle.dynamic_effective_mask = [bundle.dynamic_effective_mask; dynamic_effective(:)]; %#ok<AGROW>
        bundle.main_pulse_mask = [bundle.main_pulse_mask; main_pulse_mask(:)]; %#ok<AGROW>
        bundle.finite_mask = [bundle.finite_mask; finite_mask(:)]; %#ok<AGROW>
        bundle.core_mask = [bundle.core_mask; core_mask(:)]; %#ok<AGROW>
        bundle.query_guard_mm = [bundle.query_guard_mm; repmat(query_guard, n, 1)]; %#ok<AGROW>

        valid_pass_count = valid_pass_count + any(core_mask);

        rot_freq_samples(end + 1, 1) = F_omega_deg_s(p.arrival_s) / 360; %#ok<AGROW>
    end
end

bundle.point_count = numel(bundle.x);
bundle.core_point_count = nnz(bundle.core_mask);
bundle.valid_pass_count = valid_pass_count;
bundle.rot_freq_mean_hz = mean(rot_freq_samples, 'omitnan');
bundle.rot_rpm_mean = 60 * bundle.rot_freq_mean_hz;

if ~isfinite(bundle.rot_freq_mean_hz)
    speed_now = F_omega_deg_s(W.time_center_s) / 360;
    bundle.rot_freq_mean_hz = speed_now;
    bundle.rot_rpm_mean = 60 * speed_now;
end
end


function bundle = make_empty_bundle_local()
bundle = struct();

bundle.case_name = '';
bundle.blade_id = NaN;
bundle.lap_range = [NaN NaN];
bundle.time_window_s = [NaN NaN];

bundle.x = [];
bundle.v = [];
bundle.t = [];
bundle.theta = [];
bundle.sensor_id = [];
bundle.blade_id_vec = [];
bundle.analysis_lap = [];
bundle.pass_id = [];
bundle.prev_opr_index = [];

bundle.template_v = [];
bundle.template_dv_dx = [];
bundle.template_weight = [];

bundle.inside_domain_mask = [];
bundle.inside_guard_mask = [];
bundle.dynamic_effective_mask = [];
bundle.main_pulse_mask = [];
bundle.finite_mask = [];
bundle.core_mask = [];
bundle.query_guard_mm = [];

bundle.point_count = 0;
bundle.core_point_count = 0;
bundle.valid_pass_count = 0;
bundle.rot_freq_mean_hz = NaN;
bundle.rot_rpm_mean = NaN;
end


function [t_seg, v_seg] = read_dynamic_pulse_segment_local( ...
    case_dir, case_data, sid, start_s, end_s, arrival_s, cfg, S)

if strcmpi(S.dynamic_window_mode, 'legacy_row_bounds')
    pulse_width = max(end_s - start_s, 2 / cfg.sample_rate_hz);
    pad_s = max(S.pulse_min_pad_points / cfg.sample_rate_hz, ...
                S.pulse_pad_fraction * pulse_width);

    t0 = start_s - pad_s;
    t1 = end_s + pad_s;
else
    % Original OPRCenterStd DynamicMap style:
    % fixed window centered at the target-blade passage time.
    t0 = arrival_s - S.pulse_window_sec;
    t1 = arrival_s + S.pulse_window_sec;
end

[t_seg, v_seg] = read_raw_time_window_local( ...
    case_dir, case_data, sid, t0, t1, cfg);
end


function [t, v] = read_raw_time_window_local(case_dir, case_data, sid, t0, t1, cfg)
t = [];
v = [];

if ~isfield(case_data, 'file_ids') || ~isfield(case_data, 'file_index_map')
    error('Step02 case_data must contain file_ids and file_index_map. Re-run improved Step02.');
end

for iFile = 1:numel(case_data.file_ids)
    info = case_data.file_index_map(iFile);

    if ~isfinite(info.global_first_sample) || ~isfinite(info.global_last_sample)
        continue;
    end

    file_t0 = info.global_first_sample / cfg.sample_rate_hz;
    file_t1 = info.global_last_sample / cfg.sample_rate_hz;

    if file_t1 < t0 || file_t0 > t1
        continue;
    end

    file_id = info.file_id;
    file = fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id));

    if ~isfile(file)
        continue;
    end

    raw = load_raw_mat_file_local(file);

    if isempty(raw)
        continue;
    end

    sample_index = raw(:,1) + info.offset_samples;
    tt = sample_index ./ cfg.sample_rate_hz;

    mask = tt >= t0 & tt <= t1;

    if any(mask)
        t = [t; tt(mask)]; %#ok<AGROW>
        v = [v; raw(mask,2)]; %#ok<AGROW>
    end
end

if ~isempty(t)
    [t, order] = sort(t);
    v = v(order);

    [t, unique_idx] = unique(t, 'stable');
    v = v(unique_idx);
end
end


function raw = load_raw_mat_file_local(file)
loaded = load(file);
fn = fieldnames(loaded);

if isempty(fn)
    raw = [];
    return;
end

raw = loaded.(fn{1});

if isempty(raw) || size(raw,2) < 2
    raw = [];
    return;
end

raw(raw(:,1) == 0, :) = [];
end


function [x_mm, theta_rad, prev_opr_idx, keep] = map_dynamic_times_to_oprcenterstd_local( ...
    t, opr_times, F_omega_deg_s, theta_std_deg, cfg, t_peak)

% Match the original OPRCenterStd DynamicMap route:
% every sample in one pulse segment uses the OPR center immediately before
% the blade passage time as its local angular reference.
t = t(:);

x_mm = nan(size(t));
theta_rad = nan(size(t));
prev_opr_idx = nan(size(t));

if isempty(t) || ~isfinite(t_peak)
    keep = false(size(t));
    return;
end

idx_ref = find(opr_times < t_peak, 1, 'last');

if isempty(idx_ref)
    keep = false(size(t));
    return;
end

t_ref = opr_times(idx_ref);

theta_actual_deg = map_segment_to_relative_angle_local( ...
    t_ref, t, F_omega_deg_s);

theta_rel_deg = wrap_to_180_local(theta_actual_deg - theta_std_deg);

x_mm = theta_rel_deg .* (pi / 180) .* cfg.r_tip_mm;

theta_rad = map_time_to_rotor_phase_local( ...
    opr_times, t, cfg.blades_num);

prev_opr_idx(:) = idx_ref;

keep = isfinite(x_mm) & isfinite(theta_rad);
end


function theta_points_deg = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg_s)
t_seg = t_seg(:);

theta_points_deg = nan(size(t_seg));

if isempty(t_seg) || ~isfinite(t_ref)
    return;
end

dt_first = linspace(t_ref, t_seg(1), 10);
theta_base = trapz(dt_first, F_omega_deg_s(dt_first));

w_seg = F_omega_deg_s(t_seg);
theta_rel = cumtrapz(t_seg, w_seg);

theta_points_deg = theta_base + theta_rel;
end


function theta_rot = map_time_to_rotor_phase_local(opr_times, sample_times, num_blades)
sample_times = sample_times(:);
theta_rot = nan(size(sample_times));

if numel(opr_times) <= num_blades
    return;
end

rev_anchor_times = opr_times(1:num_blades:end);
rev_anchor_times = rev_anchor_times(:);

rev_phase = 2*pi*(0:numel(rev_anchor_times)-1).';

theta_vec = interp1( ...
    rev_anchor_times, ...
    rev_phase, ...
    sample_times, ...
    'linear', ...
    'extrap');

theta_rot(:) = theta_vec;

theta_rot(sample_times < opr_times(1) | sample_times > opr_times(end)) = NaN;
end


function [template_v, template_dv_dx, template_weight, inside_domain] = ...
    evaluate_template_at_points_local(tpl, x, S)

x = x(:);

template_v = interp1( ...
    tpl.x_grid(:), ...
    tpl.v_grid(:), ...
    x, ...
    'linear', ...
    NaN);

template_dv_dx = interp1( ...
    tpl.x_grid(:), ...
    tpl.dv_dx(:), ...
    x, ...
    'linear', ...
    NaN);

if isfield(tpl, 'weight_grid') && ~isempty(tpl.weight_grid)
    template_weight = interp1( ...
        tpl.x_grid(:), ...
        tpl.weight_grid(:), ...
        x, ...
        'linear', ...
        S.weight_floor);
else
    template_weight = ones(size(x));
end

template_weight(~isfinite(template_weight)) = S.weight_floor;
template_weight = max(template_weight, S.weight_floor);

inside_domain = x >= tpl.x_domain(1) & x <= tpl.x_domain(2);
end


function dVdt = dynamic_time_gradient_local(t, v)
t = t(:);
v = v(:);

if numel(t) < 3
    dVdt = zeros(size(v));
    return;
end

dVdt = gradient(v) ./ max(gradient(t), eps);
dVdt(~isfinite(dVdt)) = 0;
end


function mask = build_dynamic_effective_mask_local(v, dVdt, template_dv_dx, tpl, S)
% Diagnostic mask only in the main query-safe route.
v = v(:);
dVdt = dVdt(:);
template_dv_dx = template_dv_dx(:);

max_tpl_grad = max(abs(tpl.dv_dx(:)), [], 'omitnan');

if ~isfinite(max_tpl_grad) || max_tpl_grad <= 0
    tpl_grad_mask = true(size(v));
else
    tpl_grad_mask = abs(template_dv_dx) >= ...
        S.dynamic_template_gradient_min_ratio * max_tpl_grad;
end

max_time_grad = max(abs(dVdt), [], 'omitnan');

if ~isfinite(max_time_grad) || max_time_grad <= 0
    time_grad_mask = true(size(v));
else
    time_grad_mask = abs(dVdt) >= ...
        S.dynamic_time_gradient_min_ratio * max_time_grad;
end

v_finite = v(isfinite(v));

if isempty(v_finite)
    peak_mask = true(size(v));
else
    v_thr = prctile(v_finite, S.dynamic_peak_quantile);

    if ~isfinite(v_thr)
        peak_mask = true(size(v));
    else
        peak_mask = v >= v_thr;
    end
end

mask = (tpl_grad_mask & time_grad_mask) | ...
       (tpl_grad_mask & peak_mask);

mask = mask & isfinite(v) & isfinite(template_dv_dx);
end


function mask = build_main_pulse_mask_local(t, v, S)
% Diagnostic mask only in the main query-safe route.
t = t(:); %#ok<NASGU>
v = v(:);

mask = true(size(v));

if ~strcmpi(S.pulse_mode, 'single') || numel(v) < 5
    return;
end

n = numel(v);

edge_n = max(1, floor(0.15 * n));

baseline = median([ ...
    v(1:edge_n); ...
    v((n-edge_n+1):n)], ...
    'omitnan');

amp = max(v, [], 'omitnan') - baseline;

if ~isfinite(amp) || amp <= 0
    return;
end

thr = baseline + 0.20 * amp;
above = v >= thr;

cc = contiguous_regions_local(above);

if isempty(cc)
    return;
end

best_score = -inf;
best_idx = 1;

for i = 1:size(cc,1)
    seg = cc(i,1):cc(i,2);
    score = max(v(seg), [], 'omitnan') - baseline;

    if score > best_score
        best_score = score;
        best_idx = i;
    end
end

mask = false(size(v));

span = cc(best_idx,1):cc(best_idx,2);

pad = max(2, round(0.15 * numel(span)));

a = max(1, cc(best_idx,1) - pad);
b = min(numel(v), cc(best_idx,2) + pad);

mask(a:b) = true;
end


function cc = contiguous_regions_local(mask)
mask = mask(:);

if isempty(mask)
    cc = [];
    return;
end

d = diff([false; mask; false]);

starts = find(d == 1);
ends = find(d == -1) - 1;

cc = [starts, ends];
end


function qg = compute_query_guard_local(x, v, tpl, S, base_mask)
% Resolve query guard for the first/core pass.
%
% The guard shrinks the template trusted domain because final fitting queries:
%   x_query = x - dx_c - eta_s - A*sin(EO*theta + phi)
%
% Adaptive mode estimates apparent displacement by locally inverting the
% low-speed template voltage and taking a high quantile.

if nargin < 5 || isempty(base_mask)
    base_mask = isfinite(x) & isfinite(v);
end

if strcmpi(S.query_guard_mode, 'fixed')
    if isempty(S.query_guard_mm)
        qg = S.amplitude_limit_mm + ...
             S.dx_c_limit_mm + ...
             S.sensor_eta_limit_mm + ...
             S.query_guard_safety_mm;
    else
        qg = S.query_guard_mm;
    end
else
    qg = S.query_guard_min_mm;

    if nnz(base_mask) >= 8
        x_sel = x(base_mask);
        v_sel = v(base_mask);

        x_static = invert_template_voltage_for_guard_local( ...
            tpl, v_sel, x_sel);

        u_app = abs(x_sel(:) - x_static(:));
        u_app = u_app(isfinite(u_app));

        if ~isempty(u_app)
            q = min(max(S.query_guard_quantile, 0), 100);
            qg = prctile(u_app, q) + S.query_guard_safety_mm;
        end
    end
end

if isempty(qg) || ~isfinite(qg)
    qg = S.query_guard_min_mm;
end

qg = max(qg, S.query_guard_min_mm);
qg = min(qg, S.query_guard_max_mm);
qg = max(0, qg);
end


function x_static = invert_template_voltage_for_guard_local(tpl, v, x_ref)
% Local inverse of bell-shaped non-parametric template.
% Because the template is not globally one-to-one, choose the voltage-matching
% grid point closest to the current x_ref branch.

xg = tpl.x_grid(:);
vg = tpl.v_grid(:);

valid_grid = isfinite(xg) & isfinite(vg);

xg = xg(valid_grid);
vg = vg(valid_grid);

x_static = nan(size(v(:)));

if isempty(xg)
    x_static = reshape(x_static, size(v));
    return;
end

for i = 1:numel(v)
    if ~isfinite(v(i)) || ~isfinite(x_ref(i))
        continue;
    end

    score = abs(vg - v(i));

    % Tie-breaker: prefer nearby x branch, avoid jumping across pulse peak.
    score = score + 1e-6 * abs(xg - x_ref(i));

    [~, idx] = min(score);
    x_static(i) = xg(idx);
end

x_static = reshape(x_static, size(v));
end


function eo_candidates = build_eo_candidates_local(rot_freq_mean_hz, freq_search_hz, eo_pad)
if ~isfinite(rot_freq_mean_hz) || rot_freq_mean_hz <= 0
    error('Invalid rot_freq_mean_hz for EO candidate generation.');
end

eo_min = floor(freq_search_hz(1) / rot_freq_mean_hz) - eo_pad;
eo_max = ceil(freq_search_hz(2) / rot_freq_mean_hz) + eo_pad;

eo_min = max(1, eo_min);
eo_max = max(eo_min, eo_max);

eo_candidates = eo_min:eo_max;
end


function sub = select_bundle_points_local(bundle, mask, S)
mask = mask(:) & bundle.finite_mask(:);

sub = bundle;

fields = { ...
    'x', ...
    'v', ...
    't', ...
    'theta', ...
    'sensor_id', ...
    'blade_id_vec', ...
    'analysis_lap', ...
    'pass_id', ...
    'prev_opr_index', ...
    'template_v', ...
    'template_dv_dx', ...
    'template_weight', ...
    'inside_domain_mask', ...
    'inside_guard_mask', ...
    'dynamic_effective_mask', ...
    'main_pulse_mask', ...
    'finite_mask', ...
    'core_mask', ...
    'query_guard_mm'};

for i = 1:numel(fields)
    f = fields{i};

    if isfield(sub, f) && numel(sub.(f)) == numel(mask)
        sub.(f) = sub.(f)(mask);
    end
end

sub.point_count = nnz(mask);
sub.core_point_count = nnz(mask);
sub.valid_pass_count = numel(unique(sub.pass_id));

sub.fit_weight = sub.template_weight(:);
sub.fit_weight(~isfinite(sub.fit_weight)) = S.weight_floor;
sub.fit_weight = max(sub.fit_weight, S.weight_floor);
end
function seed_table = solve_vp_seed_eo_scan_local(bundle, eo_candidates, S)
% First-order VP seed scan using non-parametric template gradient.
%
% Approximation:
%   V - T(x) ~= -T'(x) * [dx_c + eta_s + a*sin(EO*theta) + b*cos(EO*theta)]
%
% The solved sinusoid coefficients are:
%   y = a*sin(EO*theta) + b*cos(EO*theta)
%   A = sqrt(a^2 + b^2)
%   phi = atan2(b, a)
%
% because:
%   A*sin(EO*theta + phi)
%   = A*cos(phi)*sin(EO*theta) + A*sin(phi)*cos(EO*theta)

sensors = S.analysis_sensors(:).';
n_sensors = numel(sensors);

template = struct( ...
    'EO', NaN, ...
    'A', NaN, ...
    'phi', NaN, ...
    'dx_c', NaN, ...
    'sensor_eta', zeros(1, n_sensors), ...
    'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, ...
    'point_count', bundle.point_count, ...
    'rank_score', inf);

seed_table = repmat(template, numel(eo_candidates), 1);

x = bundle.x(:); %#ok<NASGU>
v = bundle.v(:);
theta = bundle.theta(:);
T0 = bundle.template_v(:);
Tp = bundle.template_dv_dx(:);
w = bundle.fit_weight(:);
sensor_id = bundle.sensor_id(:);

valid = isfinite(v) & ...
        isfinite(theta) & ...
        isfinite(T0) & ...
        isfinite(Tp) & ...
        isfinite(w) & ...
        abs(Tp) > eps;

if nnz(valid) < S.min_fit_points
    return;
end

v = v(valid);
theta = theta(valid);
T0 = T0(valid);
Tp = Tp(valid);
w = w(valid);
sensor_id = sensor_id(valid);

% Linearized displacement-like residual.
% V - T0 = -Tp * q  ->  q = -(V - T0) / Tp
q_obs = -(v - T0) ./ Tp;

% Weight should also account for template gradient reliability.
wq = w .* min(abs(Tp) ./ max(abs(Tp), [], 'omitnan'), 1);
wq(~isfinite(wq)) = S.weight_floor;
wq = max(wq, S.weight_floor);

for iEO = 1:numel(eo_candidates)
    EO = eo_candidates(iEO);

    s = sin(EO .* theta);
    c = cos(EO .* theta);

    % Columns:
    %   common dx_c
    %   sin coefficient a
    %   cos coefficient b
    %   eta for sensors 2..end, sensor 1 is gauge-fixed to zero
    X = [ones(size(q_obs)), s, c];

    for is = 2:n_sensors
        X = [X, double(sensor_id == sensors(is))]; %#ok<AGROW>
    end

    good = all(isfinite(X), 2) & isfinite(q_obs) & isfinite(wq);

    if nnz(good) < max(5, size(X,2) + 1)
        continue;
    end

    Xg = X(good, :);
    yg = q_obs(good);
    wg = sqrt(wq(good));

    Xw = Xg .* wg;
    yw = yg .* wg;

    try
        beta = Xw \ yw;
    catch
        beta = pinv(Xw) * yw;
    end

    if isempty(beta) || numel(beta) < 3 || any(~isfinite(beta))
        continue;
    end

    dx_c = beta(1);
    a = beta(2);
    b = beta(3);

    eta = zeros(1, n_sensors);
    if n_sensors > 1 && numel(beta) >= 3 + n_sensors - 1
        eta(2:end) = beta(4:end).';
    end
    eta(1) = 0;

    A = hypot(a, b);
    phi = atan2(b, a);

    A = min(max(abs(A), 0), S.amplitude_limit_mm);
    dx_c = max(min(dx_c, S.dx_c_limit_mm), -S.dx_c_limit_mm);
    eta = max(min(eta, S.sensor_eta_limit_mm), -S.sensor_eta_limit_mm);

    [wrmse, prmse] = evaluate_full_template_model_local( ...
        bundle, EO, A, phi, dx_c, eta, S);

    seed_table(iEO).EO = EO;
    seed_table(iEO).A = A;
    seed_table(iEO).phi = wrap_to_pi_local(phi);
    seed_table(iEO).dx_c = dx_c;
    seed_table(iEO).sensor_eta = eta;
    seed_table(iEO).weighted_voltage_rmse = wrmse;
    seed_table(iEO).plain_voltage_rmse = prmse;
    seed_table(iEO).point_count = bundle.point_count;
    seed_table(iEO).rank_score = wrmse;
end
end


function final_eos = select_top_eos_local(seed_table, top_k)
scores = [seed_table.rank_score];
eos = [seed_table.EO];

valid = isfinite(scores) & isfinite(eos);

if ~any(valid)
    final_eos = [];
    return;
end

scores_valid = scores(valid);
eos_valid = eos(valid);

[~, order] = sort(scores_valid, 'ascend');

if isfinite(top_k)
    order = order(1:min(top_k, numel(order)));
end

final_eos = eos_valid(order);
final_eos = unique(final_eos, 'stable');
end


function result = refine_direct_template_fit_local(bundle, seed_table, final_eos, S)
result = make_failed_result_local('no_valid_eo_candidate');

if isempty(final_eos)
    return;
end

candidate_template = struct( ...
    'EO', NaN, ...
    'A', NaN, ...
    'phi', NaN, ...
    'dx_c', NaN, ...
    'sensor_eta', [], ...
    'objective_score', inf, ...
    'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, ...
    'point_count', bundle.point_count, ...
    'exitflag', NaN);

Candidate = repmat(candidate_template, numel(final_eos), 1);

for i = 1:numel(final_eos)
    EO = final_eos(i);

    seed = get_seed_for_eo_local(seed_table, EO, S);

    p0 = pack_params_local( ...
        seed.A, ...
        seed.phi, ...
        seed.dx_c, ...
        seed.sensor_eta, ...
        S);

    obj = @(p) direct_template_objective_local(p, EO, bundle, S);

    opts = optimset( ...
        'Display', 'off', ...
        'MaxIter', S.fminsearch_max_iter, ...
        'MaxFunEvals', S.fminsearch_max_fun, ...
        'TolX', 1e-7, ...
        'TolFun', 1e-9);

    try
        [p_best, ~, exitflag] = fminsearch(obj, p0, opts);
    catch
        p_best = p0;
        exitflag = -1;
    end

    [A, phi, dx_c, eta] = unpack_params_local(p_best, S);

    [wrmse, prmse, ~, ~, objective_score] = evaluate_full_template_model_local( ...
        bundle, EO, A, phi, dx_c, eta, S);

    Candidate(i).EO = EO;
    Candidate(i).A = A;
    Candidate(i).phi = wrap_to_pi_local(phi);
    Candidate(i).dx_c = dx_c;
    Candidate(i).sensor_eta = eta;
    Candidate(i).objective_score = objective_score;
    Candidate(i).weighted_voltage_rmse = wrmse;
    Candidate(i).plain_voltage_rmse = prmse;
    Candidate(i).point_count = bundle.point_count;
    Candidate(i).exitflag = exitflag;
end

[~, best_idx] = min([Candidate.objective_score]);

best = Candidate(best_idx);

if ~isfinite(best.weighted_voltage_rmse)
    result = make_failed_result_local('all_candidates_failed');
    result.CandidateTable = struct2table(Candidate);
    return;
end

[wrmse, prmse, Vpred, coverage, objective_score] = evaluate_full_template_model_local( ...
    bundle, ...
    best.EO, ...
    best.A, ...
    best.phi, ...
    best.dx_c, ...
    best.sensor_eta, ...
    S);

result = struct();
result.status = 'ok';
result.failure_reason = '';
result.EO_id = best.EO;
result.A_id = best.A;
result.phi_id_wrapped = wrap_to_pi_local(best.phi);
result.dx_c_id = best.dx_c;
result.d0_id = best.dx_c;
result.sensor_ids = S.analysis_sensors;
result.sensor_eta_id = best.sensor_eta;
result.sensor_eta_max_abs_mm = max(abs(best.sensor_eta), [], 'omitnan');
result.fn_id = best.EO * bundle.rot_freq_mean_hz;
result.rot_freq_mean_hz = bundle.rot_freq_mean_hz;
result.rot_rpm_mean = bundle.rot_rpm_mean;
result.objective_score = objective_score;
result.weighted_voltage_rmse = wrmse;
result.plain_voltage_rmse = prmse;
result.point_count = bundle.point_count;
result.valid_segment_count = numel(unique(bundle.pass_id));
result.Coverage = coverage;
result.CandidateTable = struct2table(Candidate);
result.V_pred_preview = Vpred(1:min(end, 200));
end


function seed = get_seed_for_eo_local(seed_table, EO, S)
idx = find([seed_table.EO] == EO, 1, 'first');

if isempty(idx)
    seed = struct( ...
        'A', 0.05, ...
        'phi', 0, ...
        'dx_c', 0, ...
        'sensor_eta', zeros(1, numel(S.analysis_sensors)));
else
    seed = seed_table(idx);
end

seed.A = min(max(abs(seed.A), 0), S.amplitude_limit_mm);

if ~isfinite(seed.A)
    seed.A = min(0.05, S.amplitude_limit_mm);
end

seed.phi = wrap_to_pi_local(seed.phi);

if ~isfinite(seed.phi)
    seed.phi = 0;
end

seed.dx_c = max(min(seed.dx_c, S.dx_c_limit_mm), -S.dx_c_limit_mm);

if ~isfinite(seed.dx_c)
    seed.dx_c = 0;
end

if isempty(seed.sensor_eta) || numel(seed.sensor_eta) ~= numel(S.analysis_sensors)
    seed.sensor_eta = zeros(1, numel(S.analysis_sensors));
end

seed.sensor_eta(1) = 0;
seed.sensor_eta = max(min(seed.sensor_eta, S.sensor_eta_limit_mm), ...
                      -S.sensor_eta_limit_mm);
end


function p = pack_params_local(A, phi, dx_c, eta, S)
eta = eta(:).';

if numel(eta) ~= numel(S.analysis_sensors)
    eta = zeros(1, numel(S.analysis_sensors));
end

eta(1) = 0;

% Parameter vector:
%   [A, phi, dx_c, eta_2, eta_3, ...]
p = [A, phi, dx_c, eta(2:end)];
p(~isfinite(p)) = 0;
end


function [A, phi, dx_c, eta] = unpack_params_local(p, S)
p = p(:).';

A = p(1);
phi = p(2);
dx_c = p(3);

eta = zeros(1, numel(S.analysis_sensors));

if numel(p) > 3
    eta(2:end) = p(4:end);
end

% Soft projection for reported parameters.
% Objective also uses penalties on the unprojected values.
A = min(max(abs(A), 0), S.amplitude_limit_mm);
phi = wrap_to_pi_local(phi);
dx_c = max(min(dx_c, S.dx_c_limit_mm), -S.dx_c_limit_mm);

eta(1) = 0;
eta = max(min(eta, S.sensor_eta_limit_mm), ...
          -S.sensor_eta_limit_mm);
end


function value = direct_template_objective_local(p, EO, bundle, S)
[A, phi, dx_c, eta] = unpack_params_local(p, S);

[wrmse, ~, ~, coverage, objective_score] = evaluate_full_template_model_local( ...
    bundle, EO, A, phi, dx_c, eta, S);

if ~isfinite(objective_score)
    objective_score = 1e12;
end

penalty = 0;

% Bound penalties on original unprojected variables to guide fminsearch.
if numel(p) >= 1
    penalty = penalty + 100 * max(0, abs(p(1)) - S.amplitude_limit_mm).^2;
end

if numel(p) >= 3
    penalty = penalty + 100 * max(0, abs(p(3)) - S.dx_c_limit_mm).^2;
end

if numel(p) > 3
    penalty = penalty + 100 * sum(max(0, abs(p(4:end)) - S.sensor_eta_limit_mm).^2);
end

eta_penalty = S.sensor_eta_reg_weight_v_per_mm * ...
    sqrt(mean(eta(:).^2));

overshoot_penalty = S.overshoot_penalty_weight * ...
    coverage.overshoot_rms_mm;

value = objective_score + eta_penalty + overshoot_penalty + penalty;
end


function [wrmse, prmse, Vpred, coverage, objective_score] = evaluate_full_template_model_local( ...
    bundle, EO, A, phi, dx_c, eta, S)

sensors = S.analysis_sensors(:).';

y = A .* sin(EO .* bundle.theta(:) + phi);

eta_vec = sensor_eta_vector_local( ...
    bundle.sensor_id(:), sensors, eta);

xq = bundle.x(:) - dx_c - eta_vec(:) - y(:);

Vpred = nan(size(xq));
overshoot = zeros(size(xq));

for sid = sensors
    idx_sensor = bundle.sensor_id(:) == sid;

    if ~any(idx_sensor)
        continue;
    end

    unique_blades = unique(bundle.blade_id_vec(idx_sensor));

    for ib = 1:numel(unique_blades)
        bid = unique_blades(ib);
        idx = idx_sensor & bundle.blade_id_vec(:) == bid;

        tpl = bundle_template_from_point_local(bundle, sid, bid);

        if isempty(tpl)
            continue;
        end

        left = tpl.x_domain(1);
        right = tpl.x_domain(2);

        overshoot(idx) = max(left - xq(idx), 0) + ...
                         max(xq(idx) - right, 0);

        % Interpolate using clipped x for numerical stability.
        % Overshoot is penalized separately.
        xq_clip = min(max(xq(idx), min(tpl.x_grid)), max(tpl.x_grid));

        Vpred(idx) = interp1( ...
            tpl.x_grid(:), ...
            tpl.v_grid(:), ...
            xq_clip, ...
            'linear', ...
            NaN);
    end
end

valid = isfinite(Vpred) & ...
        isfinite(bundle.v(:)) & ...
        isfinite(bundle.fit_weight(:));

if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
    objective_score = inf;
    corr_score = NaN;
else
    res = bundle.v(:) - Vpred(:);
    w = bundle.fit_weight(:);

    residual_sse = nansum(w(valid) .* res(valid).^2);
    wrmse = sqrt(residual_sse / max(nansum(w(valid)), eps));

    prmse = sqrt(nanmean(res(valid).^2));

    rmse_objective = residual_sse / max(nansum(w(valid)), eps);
    if isfield(S, 'objective_mode') && strcmpi(S.objective_mode, 'corr')
        corr_score = compute_weighted_sensor_correlation_score_local( ...
            bundle, Vpred, valid, w);
        objective_score = max(0, 1 - corr_score);
    else
        corr_score = NaN;
        objective_score = rmse_objective;
    end
end

coverage = struct();
coverage.clamp_fraction = mean(overshoot > 0, 'omitnan');
coverage.overshoot_rms_mm = sqrt(nanmean(overshoot(:).^2));
coverage.max_overshoot_mm = max(overshoot(:), [], 'omitnan');
coverage.correlation_score = corr_score;
coverage.objective_score = objective_score;
end


function score = compute_weighted_sensor_correlation_score_local(bundle, Vpred, valid, w)
sensors = unique(bundle.sensor_id(valid));
scores = [];
weights = [];

for i = 1:numel(sensors)
    sid = sensors(i);
    idx = valid & bundle.sensor_id(:) == sid;
    if nnz(idx) < 10
        continue;
    end

    y = bundle.v(idx);
    yp = Vpred(idx);
    ww = w(idx);
    ww(~isfinite(ww)) = 0;
    if sum(ww) <= 0
        continue;
    end

    y0 = y - sum(ww .* y) / sum(ww);
    yp0 = yp - sum(ww .* yp) / sum(ww);
    denom = sqrt(sum(ww .* y0.^2) * sum(ww .* yp0.^2));
    if denom <= eps
        continue;
    end
    scores(end+1, 1) = sum(ww .* y0 .* yp0) / denom; %#ok<AGROW>
    weights(end+1, 1) = sum(ww); %#ok<AGROW>
end

if isempty(scores)
    score = -1;
else
    score = sum(weights .* scores) / sum(weights);
end
end


function tpl = bundle_template_from_point_local(bundle, sid, blade_id)
tpl = [];

if ~isfield(bundle, 'TemplateMini') || isempty(bundle.TemplateMini)
    return;
end

for i = 1:numel(bundle.TemplateMini)
    if bundle.TemplateMini(i).sensor_id == sid && ...
            bundle.TemplateMini(i).blade_id == blade_id

        tpl = bundle.TemplateMini(i);
        return;
    end
end
end


function eta_vec = sensor_eta_vector_local(sensor_id, sensors, eta)
sensor_id = sensor_id(:);
eta_vec = zeros(size(sensor_id));

for i = 1:numel(sensors)
    eta_vec(sensor_id == sensors(i)) = eta(i);
end
end


function expanded_mask = build_phase_safe_expansion_mask_local(bundle, result, S)
% Second-stage deterministic phase-safe expansion.
%
% Use fitted vibration/offset parameters to compute:
%   x_query = x - dx_c - eta_s - A*sin(EO*theta + phi)
%
% A candidate point is expanded into the fitting set only if x_query stays
% inside the Step04 trusted template domain.

sensors = S.analysis_sensors(:).';

eta = result.sensor_eta_id(:).';

if numel(eta) ~= numel(sensors)
    eta = zeros(1, numel(sensors));
end

y = result.A_id .* sin(result.EO_id .* bundle.theta(:) + result.phi_id_wrapped);

eta_vec = sensor_eta_vector_local( ...
    bundle.sensor_id(:), sensors, eta);

x_query = bundle.x(:) - result.dx_c_id - eta_vec(:) - y(:);

expanded_mask = false(size(x_query));

for sid = sensors
    idx_sensor = bundle.sensor_id(:) == sid;

    if ~any(idx_sensor)
        continue;
    end

    unique_blades = unique(bundle.blade_id_vec(idx_sensor));

    for ib = 1:numel(unique_blades)
        bid = unique_blades(ib);
        idx = idx_sensor & bundle.blade_id_vec(:) == bid;

        tpl = bundle_template_from_point_local(bundle, sid, bid);

        if isempty(tpl)
            continue;
        end

        expanded_mask(idx) = ...
            x_query(idx) >= tpl.x_domain(1) + S.phase_safe_margin_mm & ...
            x_query(idx) <= tpl.x_domain(2) - S.phase_safe_margin_mm;
    end
end

expanded_mask = expanded_mask & bundle.finite_mask(:);
end


function [chosen_bundle, chosen_result, chosen_seed_table, expand_info] = choose_step05_pass_local( ...
    core_bundle, core_result, core_seed_table, ...
    expanded_bundle, expanded_result, expanded_seed_table, ...
    expand_info, S)

chosen_bundle = core_bundle;
chosen_result = core_result;
chosen_seed_table = core_seed_table;

expand_info.chosen_pass = 'core_query_safe';

if S.choose_better_pass && ...
        ~isempty(expanded_bundle) && ...
        ~isempty(expanded_result) && ...
        isstruct(expanded_result) && ...
        isfield(expanded_result, 'status') && ...
        strcmpi(expanded_result.status, 'ok')

    improve = core_result.weighted_voltage_rmse - ...
              expanded_result.weighted_voltage_rmse;

    if isfinite(improve) && improve >= S.better_pass_rmse_tolerance_v
        chosen_bundle = expanded_bundle;
        chosen_result = expanded_result;
        chosen_seed_table = expanded_seed_table;
        expand_info.chosen_pass = 'phase_safe_expanded';
    end
end

chosen_result.fit_stage = expand_info.chosen_pass;
end


function bundle = attach_template_mini_records_local(bundle, Template, sensors, blade_id)
mini = repmat( ...
    struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'x_grid', [], ...
    'v_grid', [], ...
    'dv_dx', [], ...
    'x_domain', [NaN NaN]), ...
    numel(sensors), 1);

for i = 1:numel(sensors)
    tpl = get_template_entry_local(Template, sensors(i), blade_id);

    mini(i).sensor_id = sensors(i);
    mini(i).blade_id = blade_id;
    mini(i).x_grid = tpl.x_grid(:);
    mini(i).v_grid = tpl.v_grid(:);
    mini(i).dv_dx = tpl.dv_dx(:);
    mini(i).x_domain = tpl.x_domain;
end

bundle.TemplateMini = mini;
end


function tpl = get_template_entry_local(Template, sid, blade_id)
tpl = [];

if ~isfield(Template, 'SensorBlade')
    return;
end

for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id

        tpl = Template.SensorBlade(i);
        return;
    end
end
end
function row = make_empty_trend_row_local()
row = struct( ...
    'window_id', NaN, ...
    'lap_start', NaN, ...
    'lap_end', NaN, ...
    'window_center_time_s', NaN, ...
    'rot_freq_mean_hz', NaN, ...
    'rot_rpm_mean', NaN, ...
    'A_id', NaN, ...
    'EO_id', NaN, ...
    'fn_id', NaN, ...
    'phi_id_wrapped', NaN, ...
    'dx_c_id', NaN, ...
    'd0_id', NaN, ...
    'sensor_eta_max_abs_mm', NaN, ...
    'weighted_voltage_rmse', NaN, ...
    'plain_voltage_rmse', NaN, ...
    'valid_segment_count', NaN, ...
    'point_count', NaN, ...
    'core_point_count', NaN, ...
    'expanded_point_count', NaN, ...
    'fit_stage', '', ...
    'status', '', ...
    'failure_reason', '');
end


function Wres = make_empty_window_result_local()
Wres = struct( ...
    'window_id', NaN, ...
    'lap_range', [NaN NaN], ...
    'time_window_s', [NaN NaN], ...
    'Result', [], ...
    'CoreResult', [], ...
    'ExpandedResult', [], ...
    'ExpansionInfo', [], ...
    'CandidateBundlePreview', [], ...
    'BundlePreview', [], ...
    'CoreBundlePreview', [], ...
    'ExpandedBundlePreview', [], ...
    'CorePointCount', NaN, ...
    'CandidatePointCount', NaN, ...
    'SeedTable', []);
end


function [trend_row, window_result] = make_failed_window_result_local(W, iw, bundle, reason)
trend_row = make_empty_trend_row_local();

trend_row.window_id = iw;
trend_row.lap_start = W.lap_start;
trend_row.lap_end = W.lap_end;
trend_row.window_center_time_s = W.time_center_s;

if isstruct(bundle)
    if isfield(bundle, 'rot_freq_mean_hz')
        trend_row.rot_freq_mean_hz = bundle.rot_freq_mean_hz;
        trend_row.rot_rpm_mean = 60 * bundle.rot_freq_mean_hz;
    end

    if isfield(bundle, 'point_count')
        trend_row.point_count = bundle.point_count;
    end

    if isfield(bundle, 'core_point_count')
        trend_row.core_point_count = bundle.core_point_count;
    end

    if isfield(bundle, 'valid_pass_count')
        trend_row.valid_segment_count = bundle.valid_pass_count;
    end
end

trend_row.status = 'failed';
trend_row.failure_reason = reason;
trend_row.fit_stage = 'failed';

window_result = make_empty_window_result_local();
window_result.window_id = iw;
window_result.lap_range = [W.lap_start, W.lap_end];
window_result.time_window_s = [W.time_start_s, W.time_end_s];
window_result.Result = make_failed_result_local(reason);
window_result.CandidateBundlePreview = bundle;
window_result.BundlePreview = [];
window_result.CoreBundlePreview = [];
window_result.ExpandedBundlePreview = [];
window_result.CandidatePointCount = get_bundle_point_count_local(bundle);
window_result.CorePointCount = get_bundle_core_point_count_local(bundle);
end


function n = get_bundle_point_count_local(bundle)
if isstruct(bundle) && isfield(bundle, 'point_count')
    n = bundle.point_count;
else
    n = NaN;
end
end


function n = get_bundle_core_point_count_local(bundle)
if isstruct(bundle) && isfield(bundle, 'core_point_count')
    n = bundle.core_point_count;
else
    n = NaN;
end
end


function result = make_failed_result_local(reason)
result = struct();

result.status = 'failed';
result.failure_reason = reason;

result.EO_id = NaN;
result.A_id = NaN;
result.phi_id_wrapped = NaN;
result.dx_c_id = NaN;
result.d0_id = NaN;
result.sensor_ids = [];
result.sensor_eta_id = [];
result.sensor_eta_max_abs_mm = NaN;
result.fn_id = NaN;
result.rot_freq_mean_hz = NaN;
result.rot_rpm_mean = NaN;
result.weighted_voltage_rmse = inf;
result.plain_voltage_rmse = inf;
result.point_count = 0;
result.valid_segment_count = 0;
result.fit_stage = 'failed';

result.Coverage = struct( ...
    'clamp_fraction', NaN, ...
    'overshoot_rms_mm', NaN, ...
    'max_overshoot_mm', NaN);

result.CandidateTable = table();
result.VPSeedTable = table();
result.BundleSummary = struct();
result.ReconstructionPreview = [];
end


function row = build_trend_row_local(iw, W, result, bundle)
row = make_empty_trend_row_local();

row.window_id = iw;
row.lap_start = W.lap_start;
row.lap_end = W.lap_end;
row.window_center_time_s = W.time_center_s;

if isstruct(bundle)
    if isfield(bundle, 'rot_freq_mean_hz')
        row.rot_freq_mean_hz = bundle.rot_freq_mean_hz;
        row.rot_rpm_mean = 60 * bundle.rot_freq_mean_hz;
    end

    if isfield(bundle, 'point_count')
        row.point_count = bundle.point_count;
    end

    if isfield(bundle, 'core_point_count')
        row.core_point_count = bundle.core_point_count;
    end

    if isfield(bundle, 'valid_pass_count')
        row.valid_segment_count = bundle.valid_pass_count;
    end
end

if isempty(result) || ~isstruct(result)
    row.status = 'failed';
    row.failure_reason = 'empty_result';
    return;
end

row.status = result.status;
row.failure_reason = result.failure_reason;

if isfield(result, 'fit_stage')
    row.fit_stage = result.fit_stage;
end

if strcmpi(result.status, 'ok')
    row.A_id = result.A_id;
    row.EO_id = result.EO_id;
    row.fn_id = result.fn_id;
    row.phi_id_wrapped = result.phi_id_wrapped;
    row.dx_c_id = result.dx_c_id;
    row.d0_id = result.d0_id;
    row.sensor_eta_max_abs_mm = result.sensor_eta_max_abs_mm;
    row.weighted_voltage_rmse = result.weighted_voltage_rmse;
    row.plain_voltage_rmse = result.plain_voltage_rmse;
    row.point_count = result.point_count;
    row.valid_segment_count = result.valid_segment_count;

    if isfield(result, 'core_point_count')
        row.core_point_count = result.core_point_count;
    end

    if isfield(result, 'expanded_point_count')
        row.expanded_point_count = result.expanded_point_count;
    end
end
end


function summary = summarize_bundle_local(final_bundle, all_bundle)
summary = struct();

summary.final_point_count = get_bundle_point_count_local(final_bundle);
summary.core_point_count = get_bundle_core_point_count_local(final_bundle);
summary.candidate_point_count = get_bundle_point_count_local(all_bundle);
summary.valid_pass_count = get_bundle_valid_pass_count_local(final_bundle);

if isstruct(all_bundle) && isfield(all_bundle, 'inside_domain_mask')
    summary.inside_domain_count = nnz(all_bundle.inside_domain_mask);
else
    summary.inside_domain_count = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'inside_guard_mask')
    summary.inside_guard_count = nnz(all_bundle.inside_guard_mask);
else
    summary.inside_guard_count = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'dynamic_effective_mask')
    summary.dynamic_effective_count = nnz(all_bundle.dynamic_effective_mask);
else
    summary.dynamic_effective_count = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'main_pulse_mask')
    summary.main_pulse_count = nnz(all_bundle.main_pulse_mask);
else
    summary.main_pulse_count = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'core_mask')
    summary.core_mask_count = nnz(all_bundle.core_mask);
else
    summary.core_mask_count = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'sensor_id') && ~isempty(all_bundle.sensor_id)
    summary.sensor_count_present = numel(unique(all_bundle.sensor_id));
else
    summary.sensor_count_present = NaN;
end

if isstruct(all_bundle) && isfield(all_bundle, 'query_guard_mm') && ~isempty(all_bundle.query_guard_mm)
    summary.query_guard_mean_mm = mean(all_bundle.query_guard_mm, 'omitnan');
    summary.query_guard_max_mm = max(all_bundle.query_guard_mm, [], 'omitnan');
else
    summary.query_guard_mean_mm = NaN;
    summary.query_guard_max_mm = NaN;
end
end


function n = get_bundle_valid_pass_count_local(bundle)
if isstruct(bundle) && isfield(bundle, 'valid_pass_count')
    n = bundle.valid_pass_count;
elseif isstruct(bundle) && isfield(bundle, 'pass_id') && ~isempty(bundle.pass_id)
    n = numel(unique(bundle.pass_id));
else
    n = NaN;
end
end


function preview = downsample_bundle_for_storage_local(bundle, max_points)
preview = bundle;

if isempty(bundle) || ~isstruct(bundle) || ~isfield(bundle, 'x')
    return;
end

n = numel(bundle.x);

if n <= max_points
    preview = bundle;
    return;
end

idx = unique(round(linspace(1, n, max_points)));

fields = fieldnames(preview);

for i = 1:numel(fields)
    f = fields{i};

    value = preview.(f);

    if isnumeric(value) || islogical(value)
        if numel(value) == n
            preview.(f) = value(idx);
        end
    end
end

preview.point_count = numel(idx);

if isfield(preview, 'core_mask')
    preview.core_point_count = nnz(preview.core_mask);
end

if isfield(preview, 'pass_id')
    preview.valid_pass_count = numel(unique(preview.pass_id));
end
end


function rec = build_reconstruction_preview_local(bundle, result, max_points)
rec = struct();

if isempty(bundle) || ~isstruct(bundle) || ~strcmpi(result.status, 'ok')
    return;
end

idx = 1:numel(bundle.x);

if numel(idx) > max_points
    idx = unique(round(linspace(1, numel(bundle.x), max_points)));
end

[~, ~, Vpred] = evaluate_full_template_model_local( ...
    bundle, ...
    result.EO_id, ...
    result.A_id, ...
    result.phi_id_wrapped, ...
    result.dx_c_id, ...
    result.sensor_eta_id, ...
    make_S_like_from_result_local(result));

rec.x = bundle.x(idx);
rec.v_obs = bundle.v(idx);
rec.v_pred = Vpred(idx);
rec.t = bundle.t(idx);
rec.theta = bundle.theta(idx);
rec.sensor_id = bundle.sensor_id(idx);
rec.blade_id_vec = bundle.blade_id_vec(idx);
rec.pass_id = bundle.pass_id(idx);
rec.analysis_lap = bundle.analysis_lap(idx);
end


function S = make_S_like_from_result_local(result)
S = struct();
S.analysis_sensors = result.sensor_ids;
S.weight_floor = 0.05;
S.overshoot_penalty_weight = 100;
end


function Summary = build_resonance_summary_local(Trend)
Summary = struct();

if isempty(Trend) || height(Trend) == 0
    Summary.valid_window_count = 0;
    Summary.best_window_id = NaN;
    Summary.dominant_eo = NaN;
    Summary.mean_freq_hz = NaN;
    Summary.median_freq_hz = NaN;
    return;
end

valid = strcmpi(string(Trend.status), 'ok') & ...
        isfinite(Trend.weighted_voltage_rmse);

Summary.valid_window_count = nnz(valid);

if ~any(valid)
    Summary.best_window_id = NaN;
    Summary.dominant_eo = NaN;
    Summary.mean_freq_hz = NaN;
    Summary.median_freq_hz = NaN;
    return;
end

valid_rows = find(valid);

[~, local_best] = min(Trend.weighted_voltage_rmse(valid));

Summary.best_window_id = Trend.window_id(valid_rows(local_best));

eos = Trend.EO_id(valid);

Summary.dominant_eo = mode(eos);
Summary.mean_freq_hz = mean(Trend.fn_id(valid), 'omitnan');
Summary.median_freq_hz = median(Trend.fn_id(valid), 'omitnan');
end


function plot_step05_result_local(Result, S, fig_case_dir, show_plots)
if exist(fig_case_dir, 'dir') ~= 7
    mkdir(fig_case_dir);
end

T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];
visible = visibility_state_local(show_plots);

plot_step05_trend_overview_local(Result, S, fig_case_dir, visible);
plot_step05_quality_overview_local(Result, S, fig_case_dir, visible);
plot_step05_mask_overview_local(Result, S, fig_case_dir, visible);
plot_step05_eo_landscape_local(Result, S, fig_case_dir, visible);
plot_step05_sensor_eta_local(Result, S, fig_case_dir, visible);

ok_mask = strcmpi(string(T.status), 'ok') & isfinite(T.weighted_voltage_rmse);

if any(ok_mask)
    rmse_ok = T.weighted_voltage_rmse;
    rmse_ok(~ok_mask) = inf;
    [~, best_iw] = min(rmse_ok);

    rmse_worst = T.weighted_voltage_rmse;
    rmse_worst(~ok_mask) = -inf;
    [~, worst_iw] = max(rmse_worst);

    plot_step05_window_reconstruction_local( ...
        Result, S, fig_case_dir, visible, best_iw, 'BestWindow');

    plot_step05_window_xdomain_local( ...
        Result, S, fig_case_dir, visible, best_iw, 'BestWindow');

    plot_step05_window_mask_debug_local( ...
        Result, S, fig_case_dir, visible, best_iw, 'BestWindow');

    if isfinite(worst_iw) && worst_iw ~= best_iw
        plot_step05_window_reconstruction_local( ...
            Result, S, fig_case_dir, visible, worst_iw, 'WorstWindow');

        plot_step05_window_xdomain_local( ...
            Result, S, fig_case_dir, visible, worst_iw, 'WorstWindow');

        plot_step05_window_mask_debug_local( ...
            Result, S, fig_case_dir, visible, worst_iw, 'WorstWindow');
    end
end

if show_plots
    fprintf('Step05 plots saved for case=%s, B%d, %s\n', ...
        Result.Case_Name, blade_id, sensor_tag);
end
end


function plot_step05_trend_overview_local(Result, S, fig_case_dir, visible)
T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 trend overview B%d %s', blade_id, sensor_tag), ...
    'Color', 'w', ...
    'Position', [60 50 1500 900], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_window_series_local(nexttile, T.window_id, T.A_id, 'A (mm)', 'Amplitude');
plot_window_series_local(nexttile, T.window_id, T.EO_id, 'EO', 'Identified EO');
plot_window_series_local(nexttile, T.window_id, T.fn_id, 'f (Hz)', 'Synchronous frequency');
plot_window_series_local(nexttile, T.window_id, T.dx_c_id, 'dx_c (mm)', 'Window-level offset');
plot_window_series_local(nexttile, T.window_id, T.sensor_eta_max_abs_mm, 'max |eta_s| (mm)', 'Sensor eta max abs');
plot_window_series_local(nexttile, T.window_id, T.weighted_voltage_rmse, 'weighted RMSE (V)', 'Weighted voltage residual');
plot_window_series_local(nexttile, T.window_id, T.point_count, 'points', 'Final fit point count');
plot_window_series_local(nexttile, T.window_id, T.valid_segment_count, 'segments', 'Valid segment count');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

ok = strcmpi(string(T.status), 'ok');

stem(ax, T.window_id(ok), ones(nnz(ok), 1), ...
    'filled', 'DisplayName', 'ok');

stem(ax, T.window_id(~ok), zeros(nnz(~ok), 1), ...
    'filled', 'DisplayName', 'failed');

yticks(ax, [0 1]);
yticklabels(ax, {'failed', 'ok'});
xlabel(ax, 'Window');
ylabel(ax, 'status');
title(ax, 'Window status');
legend(ax, 'Location', 'best');

sgtitle(sprintf('Step05 trend overview | case=%s, B%d, %s', ...
    Result.Case_Name, blade_id, sensor_tag), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_TrendOverview_B%d_%s', blade_id, sensor_tag));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_quality_overview_local(Result, S, fig_case_dir, visible)
T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 quality overview B%d %s', blade_id, sensor_tag), ...
    'Color', 'w', ...
    'Position', [80 70 1500 900], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ok = strcmpi(string(T.status), 'ok') & ...
     isfinite(T.point_count) & ...
     isfinite(T.weighted_voltage_rmse);

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

scatter(ax, T.point_count(ok), T.weighted_voltage_rmse(ok), ...
    36, T.EO_id(ok), 'filled');

if any(ok)
    cb = colorbar(ax);
    cb.Label.String = 'EO';
end

xlabel(ax, 'fit point count');
ylabel(ax, 'weighted RMSE (V)');
title(ax, 'RMSE vs fit point count');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

scatter(ax, T.dx_c_id(ok), T.sensor_eta_max_abs_mm(ok), ...
    36, T.weighted_voltage_rmse(ok), 'filled');

if any(ok)
    cb = colorbar(ax);
    cb.Label.String = 'weighted RMSE (V)';
end

xlabel(ax, 'dx_c (mm)');
ylabel(ax, 'max |eta_s| (mm)');
title(ax, 'Offset coupling check');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

plot(ax, T.window_id, T.plain_voltage_rmse, ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'plain');

plot(ax, T.window_id, T.weighted_voltage_rmse, ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'weighted');

xlabel(ax, 'Window');
ylabel(ax, 'RMSE (V)');
title(ax, 'Weighted vs plain RMSE');
legend(ax, 'Location', 'best');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

plot(ax, T.window_id, T.core_point_count, ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'core');

plot(ax, T.window_id, T.point_count, ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'final');

plot(ax, T.window_id, T.expanded_point_count, ...
    'd-', 'LineWidth', 1.0, 'DisplayName', 'expanded');

xlabel(ax, 'Window');
ylabel(ax, 'point count');
title(ax, 'Core/final/expanded point count');
legend(ax, 'Location', 'best');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

coverage = collect_step05_coverage_local(Result);

plot(ax, T.window_id, coverage.clamp_fraction, ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'clamp fraction');

plot(ax, T.window_id, coverage.max_overshoot_mm, ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'max overshoot (mm)');

xlabel(ax, 'Window');
ylabel(ax, 'coverage metric');
title(ax, 'Template-domain coverage diagnostics');
legend(ax, 'Location', 'best');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

status_codes = nan(height(T), 1);
status_codes(strcmpi(string(T.status), 'ok')) = 1;
status_codes(~strcmpi(string(T.status), 'ok')) = 0;

plot(ax, T.window_id, status_codes, 'o-', 'LineWidth', 1.0);

yticks(ax, [0 1]);
yticklabels(ax, {'failed', 'ok'});
xlabel(ax, 'Window');
ylabel(ax, 'status');

for i = 1:height(T)
    if ~strcmpi(string(T.status(i)), 'ok') && ...
            strlength(string(T.failure_reason(i))) > 0

        text(ax, T.window_id(i), 0.05, ...
            char(string(T.failure_reason(i))), ...
            'Rotation', 30, ...
            'FontSize', 8, ...
            'Interpreter', 'none');
    end
end

title(ax, 'Failure-reason quick scan');

sgtitle(sprintf('Step05 quality diagnostics | case=%s, B%d, %s', ...
    Result.Case_Name, blade_id, sensor_tag), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_QualityDiagnostics_B%d_%s', blade_id, sensor_tag));

if strcmpi(visible, 'off')
    close(fig);
end
end
function plot_step05_mask_overview_local(Result, S, fig_case_dir, visible)
T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];
M = collect_step05_mask_counts_local(Result);

fig = figure( ...
    'Name', sprintf('Step05 mask overview B%d %s', blade_id, sensor_tag), ...
    'Color', 'w', ...
    'Position', [100 90 1500 900], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

plot(ax, T.window_id, M.candidate_count, ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'candidate');

plot(ax, T.window_id, M.inside_domain_count, ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'inside domain');

plot(ax, T.window_id, M.inside_guard_count, ...
    '^-', 'LineWidth', 1.0, 'DisplayName', 'query-safe guard');

plot(ax, T.window_id, M.dynamic_effective_count, ...
    'd-', 'LineWidth', 1.0, 'DisplayName', 'dynamic effective diag');

plot(ax, T.window_id, M.main_pulse_count, ...
    'v-', 'LineWidth', 1.0, 'DisplayName', 'main pulse diag');

plot(ax, T.window_id, M.core_count, ...
    'p-', 'LineWidth', 1.0, 'DisplayName', 'core');

xlabel(ax, 'Window');
ylabel(ax, 'count');
title(ax, 'Mask counts');
legend(ax, 'Location', 'best');


ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

plot(ax, T.window_id, safe_divide_local(M.inside_domain_count, M.candidate_count), ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'inside/candidate');

plot(ax, T.window_id, safe_divide_local(M.inside_guard_count, M.candidate_count), ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'guard/candidate');

plot(ax, T.window_id, safe_divide_local(M.core_count, M.candidate_count), ...
    '^-', 'LineWidth', 1.0, 'DisplayName', 'core/candidate');

plot(ax, T.window_id, safe_divide_local(M.dynamic_effective_count, M.candidate_count), ...
    'd-', 'LineWidth', 1.0, 'DisplayName', 'dynamic diag/candidate');

plot(ax, T.window_id, safe_divide_local(M.main_pulse_count, M.candidate_count), ...
    'v-', 'LineWidth', 1.0, 'DisplayName', 'main pulse diag/candidate');

xlabel(ax, 'Window');
ylabel(ax, 'ratio');
ylim(ax, [0 1.05]);
title(ax, 'Mask ratios');
legend(ax, 'Location', 'best');


ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

mask_matrix = [
    M.candidate_count(:).';
    M.inside_domain_count(:).';
    M.inside_guard_count(:).';
    M.dynamic_effective_count(:).';
    M.main_pulse_count(:).';
    M.core_count(:).'
    ];

imagesc(ax, [min(T.window_id)-0.5, max(T.window_id)+0.5], ...
    [1 6], mask_matrix);

yticks(ax, 1:6);
yticklabels(ax, {'candidate', 'inside', 'guard', 'dyn diag', 'main diag', 'core'});
xlabel(ax, 'Window');
title(ax, 'Mask-count heatmap');
colorbar(ax);


ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

plot(ax, T.window_id, M.valid_pass_count, ...
    'o-', 'LineWidth', 1.0, 'DisplayName', 'valid pass count');

plot(ax, T.window_id, M.sensor_count_present, ...
    's-', 'LineWidth', 1.0, 'DisplayName', 'sensors present');

xlabel(ax, 'Window');
ylabel(ax, 'count');
title(ax, 'Pass and sensor coverage');
legend(ax, 'Location', 'best');

sgtitle(sprintf('Step05 mask diagnostics | case=%s, B%d, %s', ...
    Result.Case_Name, blade_id, sensor_tag), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_MaskDiagnostics_B%d_%s', blade_id, sensor_tag));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_eo_landscape_local(Result, S, fig_case_dir, visible)
T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

[L1, eo_list1] = collect_step05_eo_landscape_local(Result, 'VPSeedTable', 'weighted_voltage_rmse');
[L2, eo_list2] = collect_step05_eo_landscape_local(Result, 'CandidateTable', 'weighted_voltage_rmse');

if isempty(L1) && isempty(L2)
    return;
end

fig = figure( ...
    'Name', sprintf('Step05 EO landscape B%d %s', blade_id, sensor_tag), ...
    'Color', 'w', ...
    'Position', [120 110 1450 850], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
if ~isempty(L1)
    imagesc(ax, T.window_id, eo_list1, L1);
    axis(ax, 'xy');
    colorbar(ax);
    hold(ax, 'on');
    plot(ax, T.window_id, T.EO_id, ...
        'w.-', 'LineWidth', 1.5, 'MarkerSize', 16, ...
        'DisplayName', 'selected EO');
    xlabel(ax, 'Window');
    ylabel(ax, 'EO');
    title(ax, 'VP-seed weighted RMSE landscape');
    legend(ax, 'Location', 'best');
else
    axis(ax, 'off');
    title(ax, 'No VP-seed landscape available');
end

ax = nexttile;
if ~isempty(L2)
    imagesc(ax, T.window_id, eo_list2, L2);
    axis(ax, 'xy');
    colorbar(ax);
    hold(ax, 'on');
    plot(ax, T.window_id, T.EO_id, ...
        'w.-', 'LineWidth', 1.5, 'MarkerSize', 16, ...
        'DisplayName', 'selected EO');
    xlabel(ax, 'Window');
    ylabel(ax, 'EO');
    title(ax, 'Full-wave candidate weighted RMSE landscape');
    legend(ax, 'Location', 'best');
else
    axis(ax, 'off');
    title(ax, 'No full-wave EO landscape available');
end

sgtitle(sprintf('Step05 EO-candidate landscape | case=%s, B%d, %s', ...
    Result.Case_Name, blade_id, sensor_tag), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_EOLandscape_B%d_%s', blade_id, sensor_tag));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_sensor_eta_local(Result, S, fig_case_dir, visible)
eta_mat = collect_step05_sensor_eta_matrix_local(Result);

if isempty(eta_mat)
    return;
end

T = Result.Trend;
blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 sensor eta B%d %s', blade_id, sensor_tag), ...
    'Color', 'w', ...
    'Position', [140 120 1450 820], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

for i = 1:numel(Result.SensorIDs)
    plot(ax, T.window_id, eta_mat(:,i), ...
        'o-', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('CH%d', Result.SensorIDs(i)));
end

xlabel(ax, 'Window');
ylabel(ax, 'eta_s (mm)');
title(ax, 'Per-sensor eta trend');
legend(ax, 'Location', 'best');

ax = nexttile;
imagesc(ax, T.window_id, 1:numel(Result.SensorIDs), eta_mat.');
axis(ax, 'xy');
colorbar(ax);
xlabel(ax, 'Window');
ylabel(ax, 'Sensor');
yticks(ax, 1:numel(Result.SensorIDs));
yticklabels(ax, compose('CH%d', Result.SensorIDs(:)));
title(ax, 'Per-sensor eta heatmap');

sgtitle(sprintf('Step05 sensor eta diagnostics | case=%s, B%d, %s', ...
    Result.Case_Name, blade_id, sensor_tag), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_SensorEta_B%d_%s', blade_id, sensor_tag));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_window_reconstruction_local(Result, S, fig_case_dir, visible, iw, label_tag)
if iw < 1 || iw > numel(Result.WindowResult)
    return;
end

WR = Result.WindowResult(iw);

if isempty(WR.Result) || ...
        ~isfield(WR.Result, 'ReconstructionPreview') || ...
        isempty(WR.Result.ReconstructionPreview)
    return;
end

R = WR.Result.ReconstructionPreview;
B = WR.CandidateBundlePreview;

blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 %s reconstruction B%d %s W%d', ...
    label_tag, blade_id, sensor_tag, iw), ...
    'Color', 'w', ...
    'Position', [160 130 1500 880], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, numel(Result.SensorIDs), 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(Result.SensorIDs)
    sid = Result.SensorIDs(i);

    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    if ~isempty(B) && isfield(B, 'sensor_id')
        mb = B.sensor_id == sid;

        if any(mb)
            plot(ax, B.t(mb), B.v(mb), ...
                '.', 'MarkerSize', 4, ...
                'DisplayName', 'candidate raw');
        end
    end

    mr = R.sensor_id == sid;

    if any(mr)
        plot(ax, R.t(mr), R.v_obs(mr), ...
            '.', 'MarkerSize', 6, ...
            'DisplayName', 'fit observed');

        plot(ax, R.t(mr), R.v_pred(mr), ...
            '.', 'MarkerSize', 6, ...
            'DisplayName', 'fit predicted');
    end

    xlabel(ax, 'time (s)');
    ylabel(ax, 'V');
    title(ax, sprintf('CH%d time-domain reconstruction', sid));
    legend(ax, 'Location', 'best');
end

r = WR.Result;

sgtitle(sprintf('%s W%d | stage=%s, EO=%d, A=%.4f mm, dx_c=%.4f mm, RMSE=%.5f V', ...
    label_tag, iw, r.fit_stage, r.EO_id, r.A_id, r.dx_c_id, r.weighted_voltage_rmse), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_%s_TimeReconstruction_B%d_%s_W%d', ...
    label_tag, blade_id, sensor_tag, iw));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_window_xdomain_local(Result, S, fig_case_dir, visible, iw, label_tag)
if iw < 1 || iw > numel(Result.WindowResult)
    return;
end

WR = Result.WindowResult(iw);

if isempty(WR.Result) || ~strcmpi(WR.Result.status, 'ok')
    return;
end

B = WR.CandidateBundlePreview;
C = WR.BundlePreview;

if isempty(B) || isempty(C)
    return;
end

r = WR.Result;

[~, ~, Vpred_core] = evaluate_full_template_model_local( ...
    C, ...
    r.EO_id, ...
    r.A_id, ...
    r.phi_id_wrapped, ...
    r.dx_c_id, ...
    r.sensor_eta_id, ...
    make_S_like_from_result_local(r));

blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 %s x-domain B%d %s W%d', ...
    label_tag, blade_id, sensor_tag, iw), ...
    'Color', 'w', ...
    'Position', [180 140 1500 900], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, numel(Result.SensorIDs), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(Result.SensorIDs)
    sid = Result.SensorIDs(i);

    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    mb = B.sensor_id == sid;
    mc = C.sensor_id == sid;

    tpl = bundle_template_from_point_local(C, sid, Result.TargetBlade);

    if any(mb)
        plot(ax, B.x(mb), B.v(mb), ...
            '.', 'MarkerSize', 4, ...
            'DisplayName', 'candidate raw');
    end

    if ~isempty(tpl)
        plot(ax, tpl.x_grid, tpl.v_grid, ...
            '-', 'LineWidth', 1.3, ...
            'DisplayName', 'low-speed template');

        xline(ax, tpl.x_domain(1), '--', ...
            'DisplayName', 'x domain');

        xline(ax, tpl.x_domain(2), '--', ...
            'HandleVisibility', 'off');
    end

    if any(mc)
        plot(ax, C.x(mc), C.v(mc), ...
            '.', 'MarkerSize', 6, ...
            'DisplayName', 'fit observed');
    end

    xlabel(ax, 'x in OPRCenterStd (mm)');
    ylabel(ax, 'V');
    title(ax, sprintf('CH%d: candidate/template/fit points', sid));
    legend(ax, 'Location', 'best');


    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    if any(mc)
        plot(ax, C.x(mc), C.v(mc), ...
            '.', 'MarkerSize', 5, ...
            'DisplayName', 'observed at x');

        plot(ax, C.x(mc), Vpred_core(mc), ...
            '.', 'MarkerSize', 5, ...
            'DisplayName', 'predicted');

        yv = r.A_id .* sin(r.EO_id .* C.theta(mc) + r.phi_id_wrapped);

        eta_vec = sensor_eta_vector_local( ...
            C.sensor_id(mc), ...
            Result.SensorIDs(:).', ...
            r.sensor_eta_id);

        xq = C.x(mc) - r.dx_c_id - eta_vec(:) - yv(:);

        yyaxis(ax, 'right');

        plot(ax, C.x(mc), xq, ...
            '.', 'MarkerSize', 4, ...
            'DisplayName', 'query x_q');

        ylabel(ax, 'x_q (mm)');

        yyaxis(ax, 'left');
    end

    if ~isempty(tpl)
        xline(ax, tpl.x_domain(1), '--', ...
            'DisplayName', 'x domain');

        xline(ax, tpl.x_domain(2), '--', ...
            'HandleVisibility', 'off');
    end

    xlabel(ax, 'x in OPRCenterStd (mm)');
    ylabel(ax, 'V');
    title(ax, sprintf('CH%d: prediction and query-coordinate check', sid));
    legend(ax, 'Location', 'best');
end

sgtitle(sprintf('%s W%d x-domain diagnostics | stage=%s, EO=%d, A=%.4f mm, RMSE=%.5f V', ...
    label_tag, iw, r.fit_stage, r.EO_id, r.A_id, r.weighted_voltage_rmse), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_%s_XDomain_B%d_%s_W%d', ...
    label_tag, blade_id, sensor_tag, iw));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_step05_window_mask_debug_local(Result, S, fig_case_dir, visible, iw, label_tag)
if iw < 1 || iw > numel(Result.WindowResult)
    return;
end

WR = Result.WindowResult(iw);
B = WR.CandidateBundlePreview;

if isempty(B)
    return;
end

blade_id = Result.TargetBlade;
sensor_tag = ['S', sprintf('%d', Result.SensorIDs)];

fig = figure( ...
    'Name', sprintf('Step05 %s mask debug B%d %s W%d', ...
    label_tag, blade_id, sensor_tag, iw), ...
    'Color', 'w', ...
    'Position', [200 150 1500 900], ...
    'Visible', visible, ...
    'NumberTitle', 'off');

tiledlayout(fig, numel(Result.SensorIDs), 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

mask_names = { ...
    'finite', ...
    'inside domain', ...
    'query guard', ...
    'dynamic diag', ...
    'main pulse diag', ...
    'core'};

for i = 1:numel(Result.SensorIDs)
    sid = Result.SensorIDs(i);

    ms = B.sensor_id == sid;

    if ~any(ms)
        nexttile; axis off;
        nexttile; axis off;
        continue;
    end

    x_local = B.x(ms);
    v_local = B.v(ms);

    [x_sort, idx] = sort(x_local);
    v_sort = v_local(idx);

    masks = [ ...
        double(B.finite_mask(ms)), ...
        double(B.inside_domain_mask(ms)), ...
        double(B.inside_guard_mask(ms)), ...
        double(B.dynamic_effective_mask(ms)), ...
        double(B.main_pulse_mask(ms)), ...
        double(B.core_mask(ms))];

    masks = masks(idx, :);

    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    imagesc(ax, x_sort, 1:numel(mask_names), masks.');
    axis(ax, 'xy');

    yticks(ax, 1:numel(mask_names));
    yticklabels(ax, mask_names);
    xlabel(ax, 'x in OPRCenterStd (mm)');
    title(ax, sprintf('CH%d mask-by-x map', sid));
    colorbar(ax);


    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    plot(ax, x_sort, v_sort, ...
        '.', 'MarkerSize', 5, ...
        'DisplayName', 'candidate raw');

    core_idx = ms;
    core_idx(ms) = B.core_mask(ms);

    if any(core_idx)
        plot(ax, B.x(core_idx), B.v(core_idx), ...
            '.', 'MarkerSize', 7, ...
            'DisplayName', 'core fit points');
    end

    if isfield(B, 'query_guard_mm') && any(ms)
        qg = median(B.query_guard_mm(ms), 'omitnan');
    else
        qg = NaN;
    end

    if isfinite(qg)
        ttl = sprintf('CH%d voltage vs x, query guard = %.3f mm', sid, qg);
    else
        ttl = sprintf('CH%d voltage vs x', sid);
    end

    xlabel(ax, 'x in OPRCenterStd (mm)');
    ylabel(ax, 'V');
    title(ax, ttl);
    legend(ax, 'Location', 'best');
end

sgtitle(sprintf('%s W%d mask diagnostics', label_tag, iw), ...
    'Interpreter', 'none');

save_step05_figure_local(fig, fig_case_dir, S, ...
    sprintf('Step05_%s_MaskDebug_B%d_%s_W%d', ...
    label_tag, blade_id, sensor_tag, iw));

if strcmpi(visible, 'off')
    close(fig);
end
end


function plot_window_series_local(ax, x, y, ylab, ttl)
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
plot(ax, x, y, 'o-', 'LineWidth', 1.1);
xlabel(ax, 'Window');
ylabel(ax, ylab);
title(ax, ttl);
end


function coverage = collect_step05_coverage_local(Result)
n = numel(Result.WindowResult);

coverage.clamp_fraction = nan(n, 1);
coverage.max_overshoot_mm = nan(n, 1);

for i = 1:n
    if isempty(Result.WindowResult(i).Result) || ...
            ~isfield(Result.WindowResult(i).Result, 'Coverage')
        continue;
    end

    C = Result.WindowResult(i).Result.Coverage;

    if isfield(C, 'clamp_fraction')
        coverage.clamp_fraction(i) = C.clamp_fraction;
    end

    if isfield(C, 'max_overshoot_mm')
        coverage.max_overshoot_mm(i) = C.max_overshoot_mm;
    end
end
end


function M = collect_step05_mask_counts_local(Result)
n = numel(Result.WindowResult);

fields = { ...
    'candidate_count', ...
    'inside_domain_count', ...
    'inside_guard_count', ...
    'dynamic_effective_count', ...
    'main_pulse_count', ...
    'core_count', ...
    'valid_pass_count', ...
    'sensor_count_present'};

for k = 1:numel(fields)
    M.(fields{k}) = nan(n, 1);
end

for i = 1:n
    WR = Result.WindowResult(i);

    if ~isempty(WR.Result) && ...
            isfield(WR.Result, 'BundleSummary') && ...
            ~isempty(WR.Result.BundleSummary)

        S = WR.Result.BundleSummary;

        M.candidate_count(i) = getfield_default_local(S, ...
            'candidate_point_count', NaN);

        M.inside_domain_count(i) = getfield_default_local(S, ...
            'inside_domain_count', NaN);

        M.inside_guard_count(i) = getfield_default_local(S, ...
            'inside_guard_count', NaN);

        M.dynamic_effective_count(i) = getfield_default_local(S, ...
            'dynamic_effective_count', NaN);

        M.main_pulse_count(i) = getfield_default_local(S, ...
            'main_pulse_count', NaN);

        M.core_count(i) = getfield_default_local(S, ...
            'core_mask_count', NaN);

        M.valid_pass_count(i) = getfield_default_local(S, ...
            'valid_pass_count', NaN);
    end

    if isfield(WR, 'CandidateBundlePreview') && ...
            ~isempty(WR.CandidateBundlePreview) && ...
            isfield(WR.CandidateBundlePreview, 'sensor_id')

        B = WR.CandidateBundlePreview;

        M.sensor_count_present(i) = numel(unique(B.sensor_id(:)));

        if ~isfinite(M.candidate_count(i))
            M.candidate_count(i) = numel(B.x);
        end

        if ~isfinite(M.inside_domain_count(i)) && isfield(B, 'inside_domain_mask')
            M.inside_domain_count(i) = nnz(B.inside_domain_mask);
        end

        if ~isfinite(M.inside_guard_count(i)) && isfield(B, 'inside_guard_mask')
            M.inside_guard_count(i) = nnz(B.inside_guard_mask);
        end

        if ~isfinite(M.dynamic_effective_count(i)) && isfield(B, 'dynamic_effective_mask')
            M.dynamic_effective_count(i) = nnz(B.dynamic_effective_mask);
        end

        if ~isfinite(M.main_pulse_count(i)) && isfield(B, 'main_pulse_mask')
            M.main_pulse_count(i) = nnz(B.main_pulse_mask);
        end

        if ~isfinite(M.core_count(i)) && isfield(B, 'core_mask')
            M.core_count(i) = nnz(B.core_mask);
        end

        if ~isfinite(M.valid_pass_count(i)) && isfield(B, 'pass_id')
            M.valid_pass_count(i) = numel(unique(B.pass_id));
        end
    end
end
end


function [L, eo_list] = collect_step05_eo_landscape_local(Result, table_field, metric_name)
all_eo = [];

for i = 1:numel(Result.WindowResult)
    WR = Result.WindowResult(i);
    TT = [];

    if ~isempty(WR.Result) && isfield(WR.Result, table_field)
        TT = WR.Result.(table_field);
    end

    if istable(TT) && any(strcmpi(TT.Properties.VariableNames, 'EO'))
        all_eo = [all_eo; TT.EO(:)]; %#ok<AGROW>
    end
end

eo_list = unique(all_eo(isfinite(all_eo))).';

if isempty(eo_list)
    L = [];
    return;
end

L = nan(numel(eo_list), numel(Result.WindowResult));

for i = 1:numel(Result.WindowResult)
    WR = Result.WindowResult(i);
    TT = [];

    if ~isempty(WR.Result) && isfield(WR.Result, table_field)
        TT = WR.Result.(table_field);
    end

    if ~istable(TT) || ...
            ~all(ismember({'EO', metric_name}, TT.Properties.VariableNames))
        continue;
    end

    for j = 1:height(TT)
        idx = find(eo_list == TT.EO(j), 1, 'first');

        if ~isempty(idx)
            L(idx, i) = TT.(metric_name)(j);
        end
    end
end
end


function eta_mat = collect_step05_sensor_eta_matrix_local(Result)
nw = numel(Result.WindowResult);
ns = numel(Result.SensorIDs);

eta_mat = nan(nw, ns);

for i = 1:nw
    WR = Result.WindowResult(i);

    if isempty(WR.Result) || ...
            ~isfield(WR.Result, 'sensor_eta_id') || ...
            isempty(WR.Result.sensor_eta_id)
        continue;
    end

    eta = WR.Result.sensor_eta_id(:).';
    eta_mat(i, 1:min(ns, numel(eta))) = eta(1:min(ns, numel(eta)));
end

if all(all(~isfinite(eta_mat)))
    eta_mat = [];
end
end


function y = safe_divide_local(a, b)
y = nan(size(a));

mask = isfinite(a) & isfinite(b) & b ~= 0;

y(mask) = a(mask) ./ b(mask);
end


function value = getfield_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function save_step05_figure_local(fig, figure_dir, S, tag)
if ~S.save_figures
    return;
end

if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end

png_file = fullfile(figure_dir, [tag, '.png']);
pdf_file = fullfile(figure_dir, [tag, '.pdf']);

try
    exportgraphics(fig, png_file, 'Resolution', 300);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    saveas(fig, png_file);
    saveas(fig, pdf_file);
end
end


function state = visibility_state_local(show_plots)
if show_plots
    state = 'on';
else
    state = 'off';
end
end


function value = get_optional_struct_local(S, name)
if isstruct(S) && isfield(S, name)
    value = S.(name);
else
    value = struct();
end
end


function opr_times = load_dynamic_opr_times_local(step02_case_dir)
jiluopr_file = fullfile(step02_case_dir, 'jiluOPR.mat');

if ~isfile(jiluopr_file)
    error('Missing jiluOPR.mat under %s.', step02_case_dir);
end

loaded = load(jiluopr_file);

if ~isfield(loaded, 'jiluOPR')
    error('jiluOPR.mat does not contain variable jiluOPR.');
end

jiluOPR = loaded.jiluOPR;

if isempty(jiluOPR) || size(jiluOPR, 2) < 1
    error('Invalid jiluOPR format.');
end

% Improved Step02 stores:
%   column 1 = OPR arrival/center time
%   column 2 = OPR start time
%   column 3 = OPR end time
%
% Older Step02 may store [start, end]. In that case this still loads column 1,
% but the metadata consistency warning should be checked in Step03.
opr_times = jiluOPR(:,1);
opr_times = opr_times(:);

opr_times = opr_times(isfinite(opr_times));

if numel(opr_times) < 10
    error('Too few valid OPR times in jiluOPR.mat.');
end

if any(diff(opr_times) <= 0)
    warning('jiluOPR times are not strictly increasing. Sorting them.');
    opr_times = sort(opr_times);
end
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end


function ang = wrap_to_pi_local(ang)
ang = mod(ang + pi, 2*pi) - pi;
end


function y = nanmean(x)
x = x(:);
x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end


function y = nansum(x)
x = x(:);
x = x(isfinite(x));

if isempty(x)
    y = 0;
else
    y = sum(x);
end
end
