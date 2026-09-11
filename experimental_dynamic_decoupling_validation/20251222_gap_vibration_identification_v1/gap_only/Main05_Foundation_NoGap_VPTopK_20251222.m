clc; close all;

% Resolve every Foundation helper from this frozen package before executing
% the original Step05 implementation. No parent-directory code is needed.
thisDir = fileparts(mfilename('fullpath'));
Setup_Paths_20251222();
addpath(fullfile(thisDir,'functions','foundation'), ...
    fullfile(thisDir,'functions','preparation'), ...
    fullfile(thisDir,'functions','calibration'), ...
    fullfile(thisDir,'functions','gap_aware'), ...
    fullfile(thisDir,'functions','utilities'));

%COMPARE_NOGAP_VPTOPK_DIRECTTEMPLATE_20251222
% VERSION: FOUNDATION_NOETA_VPTOPK_DIRECTTEMPLATE_OPR6SLOTCAL_20260717
%
% Key logic:
%   1. High-speed waveform extraction follows the original OPRCenterStd
%      DynamicMap route: fixed window around target blade passage time.
%   2. No within-pulse uniform downsampling is applied; all raw samples in
%      the fixed time window are mapped before x-domain masks.
%   3. Core fitting uses the complete per-pass main pulse waveform:
%          core_mask = finite & main_pulse & inside_domain & inside_query_guard
%   4. VP scans all EO values without sensor eta; an adaptive EO set enters
%      full nonlinear refinement.
%   5. eta_s is not part of this method. The low-speed x_c values are fixed,
%      and each high-speed window independently identifies A, phi and dx_c.
%   6. Optional expanded bundle uses active-model phase-safe x_query:
%          x_query = x - dx_c - A*sin(EO*theta + phi)
%      and only expands per-pass main pulse points whose x_query stays
%      inside the template domain.
%
% Forward model:
%   V_i ~= T_{sensor,blade}(x_i - dx_c - A*sin(EO*theta_i + phi))
%
% Usage:
%   Run this script directly after Step02, Step03 and Step04.
%   During method cleanup and paper/validation analysis, run only the
%   current target blade first:
%       STEP05_TARGET_BLADES=1
%       STEP05_DEBUG_MAX_WINDOWS=2     % chain smoke test
%       STEP05_DEBUG_MAX_WINDOWS=inf   % full target-window trend
%
% This is the no-gap comparison program. It is intentionally separate from
% the gap-aware main program and never optimizes a clearance increment.
% The feature-level baseline has been moved to
% optional_diagnostics/Diag01_TemplateGuided_FeatureBaseline_20251222.m and
% is diagnostic only. This file does not use pulse_feature_value as the
% identification observation. Legacy DynamicMap replay belongs under
% legacy_step05/ or optional_diagnostics/, not this official entry point.

cfg = BTTProjectConfig_20251222();
slotCal = load_opr_slot_calibration_step05_local(cfg);
cfg.opr_slot_angle_deg = double(slotCal.IntervalAngleDeg(:).');
cfg.opr_slot_anchor_deg = double(slotCal.AnchorAngleDeg(:).');
cfg.opr_slot_pulse_width_angle_deg = double(slotCal.PulseWidthAngleDeg(:).');
% Use the same selected resonance region as the identification program.
% BLADE_RESONANCE_REGION_ID is therefore the only formal selector for the
% target blade, start time and window plan (R01/B1, R04/B5, ...).
Pflow = ProjectionFlow_Config_20251222();
regionSelection = struct( ...
    'representativeBladeId', Pflow.identification.targetBlade, ...
    'analysisStartTimeSec', Pflow.identification.analysisStartTimeSec, ...
    'targetBladePasses', Pflow.identification.targetBladePasses, ...
    'windowBladePasses', Pflow.identification.windowBladePasses, ...
    'slidingStepBladePasses', Pflow.identification.slidingStepBladePasses, ...
    'shortTag', Pflow.identification.resonanceRegionShortTag);

%% Step05 local settings
S = struct();
S.method_name = 'foundation_noeta_vptopk_direct_template';
S.version = 'FOUNDATION_NOETA_VPTOPK_DIRECTTEMPLATE_OPR6SLOTCAL_20260717';
S.method_file_tag = 'NoEtaVPTopKDirectTemplate';

S.show_plots = true;
S.save_figures = cfg.save_figures;
S.target_cases = cfg.dynamic_cases;
S.target_blades = regionSelection.representativeBladeId;
S.run_scope = 'target_blade_regression';
S.validation_sequence = ['single blade + 2 windows -> single blade + all target windows -> ' ...
    'additional comparison blades only as needed'];
S.analysis_sensors = cfg.sensor_ids;
analysisSensorsEnv = strtrim(getenv('STEP05_ANALYSIS_SENSORS'));
if ~isempty(analysisSensorsEnv)
    parsedAnalysisSensors = sscanf(analysisSensorsEnv, '%d').';
    if isempty(parsedAnalysisSensors) || any(~ismember(parsedAnalysisSensors, cfg.sensor_ids))
        error('STEP05_ANALYSIS_SENSORS must be a subset of %s.', mat2str(cfg.sensor_ids));
    end
    S.analysis_sensors = unique(parsedAnalysisSensors, 'stable');
end
S.reference_sensor_id = S.analysis_sensors(1);

% 20251222 is validated against the successful SG Step5 chain:
% start at 50.2 s, use 20 target laps, and slide 3-lap windows by one lap.
% Later start times are useful diagnostics, but they are not the reference
% route for judging whether Step04/Step05 is reproducing the experiment.
S.analysis_start_time_s = regionSelection.analysisStartTimeSec;
S.target_laps = regionSelection.targetBladePasses;
S.window_laps = regionSelection.windowBladePasses;
S.sliding_step_laps = regionSelection.slidingStepBladePasses;
S.debug_max_windows = inf;

S.output_dir = fullfile(cfg.output_root, 'step05_noeta_vptopk_direct_template');
S.figure_dir = fullfile(cfg.figure_root, 'step05_noeta_vptopk_direct_template');
S.dynamic_map_override_file = '';
S.angle_map_mode = cfg.opr_angle_map_mode;
S.opr_slot_angle_deg = cfg.opr_slot_angle_deg;
S.opr_slot_anchor_deg = cfg.opr_slot_anchor_deg;
S.opr_slot_angle_file = cfg.opr_slot_calibration_file;

dynamic_map_override_env = strtrim(getenv('STEP05_DYNAMIC_MAP_OVERRIDE_FILE'));
if ~isempty(dynamic_map_override_env)
    error(['STEP05_DYNAMIC_MAP_OVERRIDE_FILE is a legacy replay route. ' ...
        'Use legacy_step05 or optional_diagnostics instead of the official foundation Step05.']);
end

run_tag_env = strtrim(getenv('STEP05_OUTPUT_TAG'));
if isempty(run_tag_env) && isfield(regionSelection, 'shortTag')
    run_tag_env = lower(regionSelection.shortTag);
end
if ~isempty(run_tag_env)
    run_tag_env = regexprep(run_tag_env, '[^\w\-]', '_');
    if strlength(run_tag_env) > 24
        error(['STEP05_OUTPUT_TAG is too long (%d chars). Use a short tag ' ...
            '(<=24 chars) to avoid MATLAB v7.3/HDF5 save failures on long Windows paths.'], ...
            strlength(run_tag_env));
    end
    % Keep the formal production path short enough for MATLAB v7.3/HDF5 on
    % Windows.  The result filename itself retains the full provenance.
    S.output_dir = fullfile(cfg.output_root, ['step05_', lower(run_tag_env)]);
    S.figure_dir = fullfile(cfg.figure_root, ['step05_', lower(run_tag_env)]);
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
    if any(strcmpi(debug_max_windows_env, {'all', 'inf', 'infinite'}))
        S.debug_max_windows = inf;
    else
        parsed_debug_windows = str2double(debug_max_windows_env);
        if isfinite(parsed_debug_windows) && parsed_debug_windows > 0
            S.debug_max_windows = parsed_debug_windows;
        else
            error('STEP05_DEBUG_MAX_WINDOWS must be a positive number, all, or inf.');
        end
    end
end

analysis_start_env = strtrim(getenv('STEP05_ANALYSIS_START_TIME'));
if ~isempty(analysis_start_env)
    parsed_start_time = str2double(analysis_start_env);
    if isfinite(parsed_start_time) && parsed_start_time >= 0
        S.analysis_start_time_s = parsed_start_time;
    end
end

target_laps_env = strtrim(getenv('STEP05_TARGET_LAPS'));
if ~isempty(target_laps_env)
    if strcmpi(target_laps_env, 'inf')
        S.target_laps = inf;
    else
        parsed_target_laps = str2double(target_laps_env);
        if isfinite(parsed_target_laps) && parsed_target_laps > 0
            S.target_laps = parsed_target_laps;
        end
    end
end

window_laps_env = strtrim(getenv('STEP05_WINDOW_LAPS'));
if ~isempty(window_laps_env)
    parsed_window_laps = str2double(window_laps_env);
    if isfinite(parsed_window_laps) && parsed_window_laps > 0
        S.window_laps = parsed_window_laps;
    end
end

sliding_step_laps_env = strtrim(getenv('STEP05_SLIDING_STEP_LAPS'));
if ~isempty(sliding_step_laps_env)
    parsed_sliding_step_laps = str2double(sliding_step_laps_env);
    if isfinite(parsed_sliding_step_laps) && parsed_sliding_step_laps > 0
        S.sliding_step_laps = parsed_sliding_step_laps;
    end
end

% Do not allow a stale shell override to silently build a Foundation bundle
% for a different resonance region.  Comparison blades remain allowed, but
% the time/window contract must stay tied to BLADE_RESONANCE_REGION_ID.
if ~isempty(analysis_start_env) && ...
        abs(S.analysis_start_time_s - Pflow.identification.analysisStartTimeSec) > 1e-9
    error('Step05:RegionStartOverride', ...
        ['STEP05_ANALYSIS_START_TIME=%.9g conflicts with selected region %s ' ...
         '(%.9g s). Select BLADE_RESONANCE_REGION_ID instead.'], ...
        S.analysis_start_time_s, Pflow.identification.resonanceRegionShortTag, ...
        Pflow.identification.analysisStartTimeSec);
end
if ~isempty(target_laps_env) && ...
        (~isfinite(S.target_laps) || S.target_laps ~= Pflow.identification.targetBladePasses)
    error('Step05:RegionTargetLapsOverride', ...
        'STEP05_TARGET_LAPS conflicts with selected region %s (%d laps).', ...
        Pflow.identification.resonanceRegionShortTag, Pflow.identification.targetBladePasses);
end
if ~isempty(window_laps_env) && S.window_laps ~= Pflow.identification.windowBladePasses
    error('Step05:RegionWindowLapsOverride', ...
        'STEP05_WINDOW_LAPS conflicts with selected region %s (%d laps).', ...
        Pflow.identification.resonanceRegionShortTag, Pflow.identification.windowBladePasses);
end
if ~isempty(sliding_step_laps_env) && ...
        S.sliding_step_laps ~= Pflow.identification.slidingStepBladePasses
    error('Step05:RegionSlidingOverride', ...
        'STEP05_SLIDING_STEP_LAPS conflicts with selected region %s (%d laps).', ...
        Pflow.identification.resonanceRegionShortTag, Pflow.identification.slidingStepBladePasses);
end

target_blades_env = strtrim(getenv('STEP05_TARGET_BLADES'));
if ~isempty(target_blades_env)
    blade_tokens = regexp(target_blades_env, '[,;\s]+', 'split');
    blade_tokens = blade_tokens(~cellfun('isempty', blade_tokens));
    parsed_blades = nan(size(blade_tokens));
    for iBladeToken = 1:numel(blade_tokens)
        parsed_blades(iBladeToken) = str2double(blade_tokens{iBladeToken});
    end
    parsed_blades = parsed_blades(isfinite(parsed_blades));
    parsed_blades = parsed_blades(parsed_blades >= 1 & parsed_blades <= cfg.blades_num);
    if ~isempty(parsed_blades)
        S.target_blades = unique(parsed_blades, 'stable');
    end
end

%% Step05 identification settings
S.freq_search_hz = [100 1000];
S.eo_pad = 2;

S.eo_candidate_policy = 'adaptive_full_fit';
S.vp_top_k_min = 3;
S.vp_top_k_max = 7;
S.vp_top_rank_keep = 7;
S.vp_gap_ratio_keep = 0.05;
S.include_neighbor_eo = true;
S.neighbor_eo_radius = 1;
S.vp_include_sensor_eta = false;
S.vp_bundle_mode = 'fit_bundle';
S.vp_seed_observation_mode = 'gradient_displacement';
S.vp_template_inverse_mode = 'branch_monotone';
S.vp_rank_score_mode = 'full_seed';
S.eo_decision_score_mode = 'weighted_voltage_rmse';
S.vp_gradient_min_ratio = 0.10;
S.vp_gradient_reference_quantile = 95;
S.vp_gradient_weight_power = 2;
S.fit_sensor_eta = false;
S.eta_model_policy = 'single';
S.active_eta_model = 'none';
S.eta_models = {S.active_eta_model};
S.eta_reference_sensor = S.analysis_sensors(1);
S.sensor_eta_mode = 'disabled';
S.static_eta_source = 'none';
S.static_eta_file = '';
S.sensor_eta_prior_mm = zeros(1, numel(S.analysis_sensors));
S.static_eta_enable = false;
S.static_eta_max_abs_report_mm = inf;
S.eo_low_gap_ratio_threshold = 0.03;
S.reject_ambiguous_frequency = false;
S.unconstrained_unidentifiable_gap_ratio = 0.08;
S.unconstrained_unidentifiable_min_candidates = 4;
S.joint_eo_gap_ratio_threshold = 0.03;
S.joint_missing_window_penalty = 0.50;
S.joint_min_window_fraction = 0.50;
% Real-time route: each window reports its own causal/local identification.
% Joint EO summaries are offline diagnostics only and must not overwrite
% the window-level EO_id/A/phase/dx results.
S.joint_reported_refit_enable = false;
S.joint_reported_refit_policy = 'refit_each_window_to_joint_eo';
S.joint_reported_refit_skip_status = {'empty', 'no_valid_candidate_table', ...
    'low_window_coverage_ambiguous'};

top_k_eo_env = strtrim(getenv('STEP05_TOP_K_EO'));
if ~isempty(top_k_eo_env)
    error(['STEP05_TOP_K_EO belongs to the archived top-k Step05 route. ' ...
        'The official foundation Step05 uses adaptive_full_fit.']);
end

vp_seed_observation_env = strtrim(getenv('STEP05_VP_SEED_OBSERVATION_MODE'));
if ~isempty(vp_seed_observation_env)
    vp_seed_observation_env = lower(strtrim(vp_seed_observation_env));
    allowed_vp_observation_modes = {'gradient_displacement', ...
        'static_residual_template_inversion', 'template_inversion', ...
        'static_residual', 'linear_voltage_residual', 'voltage_gradient'};
    if any(strcmpi(vp_seed_observation_env, allowed_vp_observation_modes))
        S.vp_seed_observation_mode = vp_seed_observation_env;
    else
        error('Unsupported STEP05_VP_SEED_OBSERVATION_MODE: %s', ...
            vp_seed_observation_env);
    end
end

S.amplitude_limit_mm = 0.50;
S.dx_c_limit_mm = 0.35;
S.sensor_eta_limit_mm = 0;
S.sensor_eta_reg_weight_v_per_mm = 0.02;

S.dynamic_window_mode = 'peak_centered_fixed';
S.pulse_window_sec = 6e-4;
S.pulse_pad_fraction = 0.30;
S.pulse_min_pad_points = 50;

S.core_mask_mode = 'sg_threshold_mainpulse_trust';
S.domain_margin_mm = 0.02;

S.query_guard_mode = 'step04';
S.query_guard_mm = [];
S.query_guard_min_mm = 0.12;
S.query_guard_max_mm = S.amplitude_limit_mm + S.dx_c_limit_mm + ...
    S.sensor_eta_limit_mm + 0.05;
S.query_guard_safety_mm = 0.05;
S.query_guard_quantile = 95;

S.phase_safe_expansion = false;
S.phase_safe_margin_mm = 0.03;
S.phase_safe_reference_mode = 'prev_window';
S.phase_safe_fallback_mode = 'linear_vp';
S.phase_safe_refresh_every = 5;
S.phase_safe_prev_max_clamp_fraction = 1e-6;
S.phase_safe_prev_min_final_gap_ratio = 0.005;
S.choose_better_pass = false;
S.better_pass_rmse_tolerance_v = 0;
S.main_pulse_mode = 'legacy_threshold';
S.default_sensor_threshold = 0.5;
if isfield(cfg, 'sensor_threshold_default') && ...
        isfinite(cfg.sensor_threshold_default)
    S.default_sensor_threshold = cfg.sensor_threshold_default;
end
S.enforce_parameter_bounds = true;
S.enforce_legacy_parameter_bounds = false;
S.use_template_eta_limit = false;

%% Step05 diagnostic, gate, and optimizer settings
S.pulse_mode = 'single';
S.dynamic_template_gradient_min_ratio = 0.08;
S.dynamic_time_gradient_min_ratio = 0.15;
S.dynamic_peak_quantile = 85;

S.min_fit_points = 80;
S.min_valid_pass_fraction = 0.50;
S.min_valid_pass_count = ceil(S.min_valid_pass_fraction * ...
    S.window_laps * numel(S.analysis_sensors));
S.min_valid_sensor_count = min(2, numel(S.analysis_sensors));
S.min_valid_lap_count = min(2, S.window_laps);
S.weight_floor = 0.05;
S.fit_weight_mode = 'sg_edge_template';
S.domain_selection_mode = 'hard';
S.domain_soft_margin_mm = 0;

S.fminsearch_max_iter = 800;
S.fminsearch_max_fun = 2000;
S.direct_multistart_enable = true;
S.direct_multistart_refine_top_n = 2;
S.overshoot_penalty_weight = 100;
S.invalid_query_penalty_scale = 5;
S.extrapolation_mode = 'clip_to_grid';
S.template_interp_method = 'pchip';
S.objective_mode = 'rmse';
% The gap-aware program consumes this neutral point bundle. Keep every
% selected raw point so its objective is not silently fitted on a preview.
S.store_bundle_preview_points = inf;

method_preset_env = strtrim(getenv('STEP05_METHOD_PRESET'));
if ~isempty(method_preset_env)
    error(['STEP05_METHOD_PRESET is no longer supported by the official ' ...
        'foundation Step05. Use legacy_step05 or optional_diagnostics for ' ...
        'legacy/corr/all-EO routes.']);
end

max_iter_env = strtrim(getenv('STEP05_FMINSEARCH_MAX_ITER'));
if ~isempty(max_iter_env)
    parsed_max_iter = str2double(max_iter_env);
    if isfinite(parsed_max_iter) && parsed_max_iter > 0
        S.fminsearch_max_iter = parsed_max_iter;
    else
        error('STEP05_FMINSEARCH_MAX_ITER must be a positive number.');
    end
end

max_fun_env = strtrim(getenv('STEP05_FMINSEARCH_MAX_FUN'));
if ~isempty(max_fun_env)
    parsed_max_fun = str2double(max_fun_env);
    if isfinite(parsed_max_fun) && parsed_max_fun > 0
        S.fminsearch_max_fun = parsed_max_fun;
    else
        error('STEP05_FMINSEARCH_MAX_FUN must be a positive number.');
    end
end

store_preview_env = strtrim(getenv('STEP05_STORE_BUNDLE_PREVIEW_POINTS'));
if ~isempty(store_preview_env)
    if any(strcmpi(store_preview_env, {'inf', 'infinite', 'all'}))
        S.store_bundle_preview_points = inf;
    else
        parsed_store_preview = str2double(store_preview_env);
        if isfinite(parsed_store_preview) && parsed_store_preview > 0
            S.store_bundle_preview_points = parsed_store_preview;
        else
            error('STEP05_STORE_BUNDLE_PREVIEW_POINTS must be positive, inf, infinite, or all.');
        end
    end
end

extrapolation_mode_env = strtrim(getenv('STEP05_EXTRAPOLATION_MODE'));
if ~isempty(extrapolation_mode_env)
    error('STEP05_EXTRAPOLATION_MODE changes the method route; use diagnostics/legacy scripts for alternatives.');
end

objective_mode_env = strtrim(getenv('STEP05_OBJECTIVE_MODE'));
if ~isempty(objective_mode_env)
    error('STEP05_OBJECTIVE_MODE changes the method route; official foundation Step05 uses RMSE.');
end

core_mask_mode_env = strtrim(getenv('STEP05_CORE_MASK_MODE'));
if ~isempty(core_mask_mode_env)
    error('STEP05_CORE_MASK_MODE changes the method route; official foundation Step05 uses per-pass main-pulse query-safe fitting.');
end

query_guard_mode_env = strtrim(getenv('STEP05_QUERY_GUARD_MODE'));
if ~isempty(query_guard_mode_env)
    error('STEP05_QUERY_GUARD_MODE changes the method route; official foundation Step05 uses Step04 query-safe domains.');
end

domain_selection_mode_env = strtrim(getenv('STEP05_DOMAIN_SELECTION_MODE'));
if ~isempty(domain_selection_mode_env)
    error('STEP05_DOMAIN_SELECTION_MODE changes the method route; official foundation Step05 uses hard domain selection.');
end

domain_soft_margin_env = strtrim(getenv('STEP05_DOMAIN_SOFT_MARGIN_MM'));
if ~isempty(domain_soft_margin_env)
    error('STEP05_DOMAIN_SOFT_MARGIN_MM is only meaningful for soft-domain diagnostic routes.');
end

query_guard_mm_env = strtrim(getenv('STEP05_QUERY_GUARD_MM'));
if ~isempty(query_guard_mm_env)
    error('STEP05_QUERY_GUARD_MM changes the method route; use Step04 x_query_safe_domain in the official Step05.');
end

query_guard_min_env = strtrim(getenv('STEP05_QUERY_GUARD_MIN_MM'));
if ~isempty(query_guard_min_env)
    error('STEP05_QUERY_GUARD_MIN_MM changes fallback guard behavior; use diagnostics for guard sweeps.');
end

query_guard_max_env = strtrim(getenv('STEP05_QUERY_GUARD_MAX_MM'));
if ~isempty(query_guard_max_env)
    error('STEP05_QUERY_GUARD_MAX_MM changes fallback guard behavior; use diagnostics for guard sweeps.');
end

vp_include_eta_env = strtrim(getenv('STEP05_VP_INCLUDE_SENSOR_ETA'));
if ~isempty(vp_include_eta_env)
    error('STEP05_VP_INCLUDE_SENSOR_ETA changes the VP seed model; official foundation Step05 keeps VP eta disabled.');
end

if ~isempty(strtrim(getenv('STEP05_FIT_SENSOR_ETA')))
    error(['Free-eta is not a mainline Step05 option. ' ...
        'Use Compare_Step05_FreeEta_20251222.m for the comparison route.']);
end

static_eta_file_env = strtrim(getenv('STEP05_STATIC_ETA_FILE'));
if ~isempty(static_eta_file_env)
    error('STEP05_STATIC_ETA_FILE is incompatible with the NoEta main route.');
end

phase_safe_env = strtrim(getenv('STEP05_PHASE_SAFE_EXPANSION'));
if ~isempty(phase_safe_env)
    error('STEP05_PHASE_SAFE_EXPANSION changes the method route; official foundation Step05 keeps it disabled.');
end

if ischar(S.target_cases) || isstring(S.target_cases)
    S.target_cases = cellstr(S.target_cases);
end

assert(strcmpi(S.active_eta_model, 'none') && ~S.fit_sensor_eta && ...
    ~S.vp_include_sensor_eta && ~S.static_eta_enable && ...
    all(S.sensor_eta_prior_mm == 0), ...
    'Step05:NoEtaContractViolation', ...
    'NoEta route requires eta_s=0 in VP, full-wave fitting and reporting.');
target_cases = S.target_cases;

fprintf('\n=== Step05: single-sync direct-template identification ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Run scope: %s\n', S.run_scope);
fprintf('Validation sequence: %s\n', S.validation_sequence);
fprintf('Target cases: %s\n', strjoin(cellstr(target_cases), ', '));
fprintf('Target blades: %s\n', mat2str(S.target_blades));
if isequal(S.target_blades(:).', 1:cfg.blades_num)
    fprintf('Target blade policy: all blades requested explicitly; expect mixed blade/template quality effects.\n');
else
    fprintf('Target blade policy: target blades only; use STEP05_TARGET_BLADES to add comparison blades.\n');
end
if isfinite(S.debug_max_windows)
    fprintf('Window limit: first %d windows for chain/regression check.\n', S.debug_max_windows);
else
    fprintf('Window limit: all target windows.\n');
end
fprintf('Analysis sensors: %s\n', mat2str(S.analysis_sensors));
if isfield(S, 'active_eta_model') && strcmpi(S.active_eta_model, 'none')
    fprintf('Model: V = T_{s,b}(x - dx_c - A*sin(EO*theta + phi)).\n');
else
    fprintf('Model: V = T_{s,b}(x - dx_c - eta_s - A*sin(EO*theta + phi)).\n');
end
fprintf('Eta policy: %s, active model=%s, reference CH%d.\n', ...
    S.eta_model_policy, S.active_eta_model, S.eta_reference_sensor);
fprintf('EO candidate policy: %s (VP all EO, adaptive full-fit set).\n', ...
    S.eo_candidate_policy);
fprintf('Core mask mode: %s\n', S.core_mask_mode);
fprintf('Query-safe source: %s\n', S.query_guard_mode);
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

validate_step05_template_contract_local(Template, cfg, S);
template_opr_method = getfield_default_local(Template, 'OPR_Timing_Method', '');
cfg_opr_method = getfield_default_local(cfg, 'opr_timing_method', '');
if ~isempty(cfg_opr_method) && ~isempty(template_opr_method) && ...
        ~strcmpi(strtrim(char(string(template_opr_method))), strtrim(char(string(cfg_opr_method))))
    error('Step05:OPRReferenceMismatch', ...
        'Template OPR_Timing_Method=%s, but cfg.opr_timing_method=%s. Rebuild Step01/Step04 consistently.', ...
        char(string(template_opr_method)), char(string(cfg_opr_method)));
end

if isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta')
    [S.sensor_eta_prior_mm, S.static_eta_prior_info] = ...
        load_joint_static_eta_prior_local(cfg, S);
    fprintf('Static eta source: %s\n', S.static_eta_source);
    fprintf('Static eta prior (sensors %s, relative to CH%d): %s mm\n', ...
        mat2str(S.analysis_sensors), S.reference_sensor_id, ...
        mat2str(S.sensor_eta_prior_mm, 8));
else
    S.static_eta_prior_info = table(S.analysis_sensors(:), ...
        zeros(numel(S.analysis_sensors), 1), ...
        'VariableNames', {'sensor_id', 'eta_prior_mm'});
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

Step05_TargetBlade_Summary = collect_step05_all_blade_summary_local(ResultSet);
Step05_AllBlade_Summary = Step05_TargetBlade_Summary; %#ok<NASGU>
if exist(S.output_dir, 'dir') ~= 7
    mkdir(S.output_dir);
end
writetable(Step05_TargetBlade_Summary, fullfile(S.output_dir, ...
    sprintf('Step05_TargetBlade_Summary_%s.csv', cfg.dataset)));
RunInfo = ResultSet.RunInfo; %#ok<NASGU>
save(fullfile(S.output_dir, sprintf('Step05_RunInfo_%s.mat', cfg.dataset)), ...
    'RunInfo', 'Step05_TargetBlade_Summary', 'Step05_AllBlade_Summary');


function Summary = collect_step05_all_blade_summary_local(ResultSet)
rows = struct( ...
    'case_name', "", ...
    'blade_id', NaN, ...
    'sensor_tag', "", ...
    'valid_window_count', NaN, ...
    'total_window_count', NaN, ...
    'dominant_eo', NaN, ...
    'unique_eo_count', NaN, ...
    'eo_sequence', "", ...
    'low_gap_window_count', NaN, ...
    'eo_low_gap_ratio_threshold', NaN, ...
    'median_eo_gap_ratio', NaN, ...
    'min_eo_gap_ratio', NaN, ...
    'joint_eo', NaN, ...
    'joint_freq_hz', NaN, ...
    'joint_eo_gap_ratio', NaN, ...
    'joint_eo_status', "", ...
    'mean_freq_hz', NaN, ...
    'median_freq_hz', NaN, ...
    'median_A_mm', NaN, ...
    'median_rmse_v', NaN, ...
    'best_window_id', NaN, ...
    'best_window_rmse_v', NaN);
rows = repmat(rows, 0, 1);

for iCase = 1:numel(ResultSet.Cases)
    case_name = string(ResultSet.Cases(iCase).case_name);
    BladeResult = ResultSet.Cases(iCase).BladeResult;
    for iBlade = 1:numel(BladeResult)
        if ~isfield(BladeResult(iBlade), 'Result') || isempty(BladeResult(iBlade).Result)
            continue;
        end
        R = BladeResult(iBlade).Result;
        T = R.Trend;
        ok = strcmpi(string(T.status), "ok");
        row = struct();
        row.case_name = case_name;
        row.blade_id = R.TargetBlade;
        row.sensor_tag = string(['S', sprintf('%d', R.SensorIDs)]);
        row.valid_window_count = nnz(ok);
        row.total_window_count = height(T);
        row.dominant_eo = NaN;
        row.unique_eo_count = NaN;
        row.eo_sequence = "";
        row.low_gap_window_count = NaN;
        row.eo_low_gap_ratio_threshold = ResultSet.RunInfo.step_settings.eo_low_gap_ratio_threshold;
        row.median_eo_gap_ratio = NaN;
        row.min_eo_gap_ratio = NaN;
        row.joint_eo = NaN;
        row.joint_freq_hz = NaN;
        row.joint_eo_gap_ratio = NaN;
        row.joint_eo_status = "";
        row.mean_freq_hz = NaN;
        row.median_freq_hz = NaN;
        row.median_A_mm = NaN;
        row.median_rmse_v = NaN;
        row.best_window_id = NaN;
        row.best_window_rmse_v = NaN;
        if isfield(R, 'JointEOSummary') && isstruct(R.JointEOSummary)
            row.joint_eo = getfield_default_local(R.JointEOSummary, ...
                'best_EO', getfield_default_local(R.JointEOSummary, 'dominant_eo', NaN));
            row.joint_eo_gap_ratio = getfield_default_local(R.JointEOSummary, 'gap_ratio', NaN);
            row.joint_eo_status = string(getfield_default_local(R.JointEOSummary, 'status', ''));
            if isfield(R, 'ResonanceSummary') && isstruct(R.ResonanceSummary) && ...
                    isfield(R.ResonanceSummary, 'reported_freq_hz') && ...
                    isfinite(R.ResonanceSummary.reported_freq_hz)
                row.joint_freq_hz = R.ResonanceSummary.reported_freq_hz;
            elseif isfinite(row.joint_eo) && any(isfinite(T.rot_freq_mean_hz))
                row.joint_freq_hz = row.joint_eo .* mean(T.rot_freq_mean_hz(isfinite(T.rot_freq_mean_hz)), 'omitnan');
            end
        end
        if any(ok)
            row.dominant_eo = mode(T.EO_id(ok));
            eo_ok = T.EO_id(ok);
            row.unique_eo_count = numel(unique(eo_ok(isfinite(eo_ok))));
            row.eo_sequence = strjoin(compose('%d', eo_ok(:).'), '|');
            if any(strcmpi(T.Properties.VariableNames, 'eo_rmse_gap_ratio'))
                gap = T.eo_rmse_gap_ratio(ok);
                row.low_gap_window_count = nnz(isfinite(gap) & ...
                    gap < row.eo_low_gap_ratio_threshold);
                row.median_eo_gap_ratio = median(gap, 'omitnan');
                row.min_eo_gap_ratio = min(gap, [], 'omitnan');
            end
            row.mean_freq_hz = mean(T.fn_id(ok), 'omitnan');
            row.median_freq_hz = median(T.fn_id(ok), 'omitnan');
            row.median_A_mm = median(T.A_id(ok), 'omitnan');
            row.median_rmse_v = median(T.weighted_voltage_rmse(ok), 'omitnan');
            [row.best_window_rmse_v, best_local] = min(T.weighted_voltage_rmse(ok));
            ok_rows = find(ok);
            row.best_window_id = T.window_id(ok_rows(best_local));
        elseif isfinite(row.joint_eo) && strcmpi(row.joint_eo_status, "ok")
            row.dominant_eo = row.joint_eo;
            row.unique_eo_count = 1;
            row.eo_sequence = sprintf('joint:%d', row.joint_eo);
            row.mean_freq_hz = row.joint_freq_hz;
            row.median_freq_hz = row.joint_freq_hz;
        end
        rows(end + 1, 1) = row; %#ok<AGROW>
    end
end

if isempty(rows)
    Summary = table();
else
    Summary = struct2table(rows);
end
end


function value = get_cfg_or_default_local(cfg, field_name, default_value)
if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
    value = cfg.(field_name);
else
    value = default_value;
end
end


function [eta, info] = load_joint_static_eta_prior_local(cfg, S)
sensors = S.analysis_sensors(:).';
blade_id = S.target_blades(1);
sensor_tag = ['S', sprintf('%d', sensors)];

eta_file = '';
if isfield(S, 'static_eta_file') && ~isempty(S.static_eta_file)
    eta_file = char(S.static_eta_file);
end
if isempty(eta_file)
    eta_file = fullfile(cfg.output_root, ...
        'step05_joint_static_eta_preview', ...
        sprintf('JointStaticEtaPreview_B%d_%s.mat', blade_id, sensor_tag));
end
if ~isfile(eta_file)
    error(['Formal Step05 requires a joint static eta preview. ' ...
        'Run Calibrate_Step05_JointStaticEta_FromPreview_20251222 first, ' ...
        'or set STEP05_STATIC_ETA_FILE. Missing file: %s'], eta_file);
end

loaded = load(eta_file, 'Summary');
if ~isfield(loaded, 'Summary')
    error('Joint static eta file must contain Summary: %s', eta_file);
end
Summary = loaded.Summary;

if isfield(Summary, 'TargetBlade') && double(Summary.TargetBlade) ~= blade_id
    error('Joint static eta target blade mismatch: file B%d, requested B%d.', ...
        double(Summary.TargetBlade), blade_id);
end
if isfield(Summary, 'SensorIDs') && ...
        ~isequal(double(Summary.SensorIDs(:).'), sensors)
    error('Joint static eta sensor mismatch: file %s, requested %s.', ...
        mat2str(double(Summary.SensorIDs(:).')), mat2str(sensors));
end
if ~isfield(Summary, 'Consensus') || ...
        ~isfield(Summary.Consensus, 'eta_median_mm')
    error('Joint static eta file missing Summary.Consensus.eta_median_mm: %s', ...
        eta_file);
end

eta = double(Summary.Consensus.eta_median_mm(:).');
if numel(eta) ~= numel(sensors) || any(~isfinite(eta))
    error('Invalid joint static eta vector in %s.', eta_file);
end
eta(1) = 0;

plan_iqr = nan(size(eta));
if isfield(Summary.Consensus, 'eta_plan_iqr_mm')
    plan_iqr = double(Summary.Consensus.eta_plan_iqr_mm(:).');
end
if numel(plan_iqr) ~= numel(eta)
    plan_iqr = nan(size(eta));
end

info = table(sensors(:), eta(:), plan_iqr(:), ...
    repmat(string(eta_file), numel(sensors), 1), ...
    'VariableNames', {'sensor_id', 'eta_prior_mm', ...
    'eta_plan_iqr_mm', 'joint_static_eta_file'});
end


function Template = load_step04_template_local(cfg)
if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end

template_override_file = strtrim(getenv('STEP05_TEMPLATE_OVERRIDE_FILE'));
if ~isempty(template_override_file)
    if ~isfile(template_override_file)
        error('STEP05_TEMPLATE_OVERRIDE_FILE does not exist: %s', template_override_file);
    end
    loaded = load(template_override_file, 'Template');
    Template = loaded.Template;
    fprintf('Step04 template override: %s\n', template_override_file);
    return;
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
    error('No Step04 template file found under %s. Please run Step04_Build_LowSpeed_OPRCenterStd_Template_20251222.m first.', ...
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
step03_case_dir = fullfile(cfg.step03_output_dir, case_name);

step02_file = fullfile(step02_case_dir, 'DynamicBTTFeature_20251222.mat');
step03_file = fullfile(step03_case_dir, 'BTT_Observation_Bundle_20251222.mat');

if ~isfile(step02_file)
    error('Step02 foundation product missing for %s. Please run Step02_Extract_Dynamic_BTT_20251222.m first:\n  %s', ...
        case_name, step02_file);
end

if ~isfile(step03_file)
    error('Step03 observation bundle missing for %s. Please run Step03_Build_Dynamic_Observation_Bundle_20251222.m first:\n  %s', ...
        case_name, step03_file);
end

loaded_step02 = load(step02_file, ...
    'DynamicBTTFeature', 'DynamicPulseTable', 'OPRTable', 'Step02Meta');
loaded_step03 = load(step03_file, 'BTT_Observation_Bundle', 'metadata');

DynamicPulseTable = loaded_step02.DynamicPulseTable;
OPRTable = loaded_step02.OPRTable;
BTT_Observation_Bundle = loaded_step03.BTT_Observation_Bundle;
ObservationTable = BTT_Observation_Bundle.Observation_Table;

if isfield(loaded_step02, 'Step02Meta')
    step02_metadata = loaded_step02.Step02Meta;
elseif isfield(loaded_step02, 'DynamicBTTFeature') && ...
        isfield(loaded_step02.DynamicBTTFeature, 'Metadata')
    step02_metadata = loaded_step02.DynamicBTTFeature.Metadata;
else
    step02_metadata = struct();
end

case_dir = fullfile(cfg.dataset_root, case_name);

if ~isfolder(case_dir)
    error('Dynamic raw data folder not found: %s', case_dir);
end

case_data = build_foundation_case_data_local( ...
    case_name, case_dir, DynamicPulseTable, ObservationTable, OPRTable, cfg, S);

opr_times = resolve_step05_opr_reference_times_local(OPRTable, cfg, Template);
opr_times = opr_times(isfinite(opr_times));
opr_events_per_revolution = resolve_step05_opr_events_per_revolution_local( ...
    cfg, Template, OPRTable, step02_metadata);
if numel(opr_times) <= opr_events_per_revolution
    error('Step02 OPRTable has too few valid OPR centers for %s with epr=%d.', ...
        case_name, opr_events_per_revolution);
end
[slotIndexOffset, slotAlignmentAudit] = align_dynamic_opr_slots_local( ...
    OPRTable, cfg.opr_slot_pulse_width_angle_deg, opr_events_per_revolution);
cfg.opr_slot_index_offset = slotIndexOffset;
S.opr_slot_index_offset = slotIndexOffset;
S.opr_slot_alignment_audit = slotAlignmentAudit;
fprintf('Unequal-slot OPR alignment: offset %+d, score %.6f deg, margin %.6f deg.\n', ...
    slotIndexOffset, slotAlignmentAudit.score_deg(1), ...
    slotAlignmentAudit.score_margin_deg(1));

[F_omega_deg_s, speed_diag] = build_dynamic_speed_interpolant_local( ...
    opr_times, opr_events_per_revolution);

blade_results = repmat(struct('blade_id', NaN, 'Result', []), ...
    numel(S.target_blades), 1);

for ib = 1:numel(S.target_blades)
    blade_id = S.target_blades(ib);

    fprintf('\n>>> [Step05][%s] Target blade %d\n', case_name, blade_id);

    check_step05_template_coverage_local(Template, S, blade_id);

    passes_by_sensor = build_analysis_passes_by_sensor_local( ...
        case_data, S, blade_id);

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
        Template, passes_by_sensor, WindowPlan, cfg, S, blade_id, ...
        opr_events_per_revolution);

    out_dir = fullfile(S.output_dir, case_name);
    fig_case_dir = fullfile(S.figure_dir, case_name);

    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end

    if exist(fig_case_dir, 'dir') ~= 7
        mkdir(fig_case_dir);
    end

    sensor_tag = ['S', sprintf('%d', S.analysis_sensors)];
    method_file_tag = getfield_default_local(S, ...
        'method_file_tag', 'FoundationMainPulseAdaptiveNoEtaVoltageDecision');

    result_file = fullfile(out_dir, sprintf( ...
        'Result_Step05_%s_B%d_%s_%s.mat', ...
        method_file_tag, blade_id, sensor_tag, cfg.dataset));

    trend_file = fullfile(out_dir, sprintf( ...
        'Trend_Step05_%s_B%d_%s_%s.csv', ...
        method_file_tag, blade_id, sensor_tag, cfg.dataset));

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


function case_data = build_foundation_case_data_local( ...
    case_name, case_dir, DynamicPulseTable, ObservationTable, OPRTable, cfg, S)
required_step02 = {'sensor_id','source_row','source_file_id','arrival_time_s', ...
    'pulse_start_time_s','pulse_end_time_s','assigned_blade_id','rev_id','is_valid'};
required_step03 = {'sensor_id','pulse_index','blade_id','arrival_time_s', ...
    'pulse_start_time_s','pulse_end_time_s','rev_id','is_valid'};
required_opr = {'opr_center_time_s'};

assert_required_columns_local(DynamicPulseTable, required_step02, 'Step02 DynamicPulseTable');
assert_required_columns_local(ObservationTable, required_step03, 'Step03 Observation_Table');
assert_required_columns_local(OPRTable, required_opr, 'Step02 OPRTable');

ObservationTable.source_file_id = nan(height(ObservationTable), 1);
ObservationTable.step02_source_row = nan(height(ObservationTable), 1);

for sid = S.analysis_sensors(:).'
    idx_obs = find(ObservationTable.sensor_id == sid);
    idx_dyn = find(DynamicPulseTable.sensor_id == sid);
    if isempty(idx_obs) || isempty(idx_dyn)
        continue;
    end

    [tf, loc] = ismember(ObservationTable.pulse_index(idx_obs), ...
        DynamicPulseTable.source_row(idx_dyn));

    matched_obs = idx_obs(tf);
    matched_dyn = idx_dyn(loc(tf));

    ObservationTable.source_file_id(matched_obs) = ...
        DynamicPulseTable.source_file_id(matched_dyn);
    ObservationTable.step02_source_row(matched_obs) = ...
        DynamicPulseTable.source_row(matched_dyn);

    dt = abs(ObservationTable.arrival_time_s(matched_obs) - ...
        DynamicPulseTable.arrival_time_s(matched_dyn));
    if any(dt > 5 / cfg.sample_rate_hz)
        error('Step02/Step03 arrival mismatch for CH%d. Max dt = %.9g s.', ...
            sid, max(dt));
    end
end

unmatched_valid = ObservationTable.is_valid & ...
    ismember(ObservationTable.sensor_id, S.analysis_sensors) & ...
    ~isfinite(ObservationTable.source_file_id);
if any(unmatched_valid)
    error('Step03 has %d valid observations without Step02 source_file_id.', ...
        nnz(unmatched_valid));
end

file_ids = unique(DynamicPulseTable.source_file_id( ...
    isfinite(DynamicPulseTable.source_file_id)));
file_ids = sort(file_ids(:).');
file_index_map = build_raw_file_index_map_local(file_ids, DynamicPulseTable, cfg);

case_data = struct();
case_data.case_name = case_name;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids;
case_data.file_index_map = file_index_map;
case_data.DynamicPulseTable = DynamicPulseTable;
case_data.ObservationTable = ObservationTable;
case_data.OPRTable = OPRTable;
case_data.source_mode = 'foundation_step02_step03_step04_fullraw';
end


function assert_required_columns_local(T, required, label)
for i = 1:numel(required)
    if ~ismember(required{i}, T.Properties.VariableNames)
        error('%s missing required column: %s', label, required{i});
    end
end
end


function file_index_map = build_raw_file_index_map_local(file_ids, DynamicPulseTable, cfg)
template = struct( ...
    'file_id', NaN, ...
    'global_first_sample', NaN, ...
    'global_last_sample', NaN, ...
    'offset_samples', 0);
file_index_map = repmat(template, numel(file_ids), 1);

for i = 1:numel(file_ids)
    fid = file_ids(i);
    idx = DynamicPulseTable.source_file_id == fid & ...
          isfinite(DynamicPulseTable.pulse_start_time_s) & ...
          isfinite(DynamicPulseTable.pulse_end_time_s);

    file_index_map(i).file_id = fid;
    file_index_map(i).offset_samples = 0;

    if ~any(idx)
        warning('No Step02 pulse rows available for source_file_id=%g.', fid);
        continue;
    end

    pad_samples = max(10 * cfg.sample_rate_hz, 2 * cfg.sample_rate_hz);
    first_sample = floor(min(DynamicPulseTable.pulse_start_time_s(idx)) * cfg.sample_rate_hz) - pad_samples;
    last_sample = ceil(max(DynamicPulseTable.pulse_end_time_s(idx)) * cfg.sample_rate_hz) + pad_samples;

    file_index_map(i).global_first_sample = max(0, first_sample);
    file_index_map(i).global_last_sample = last_sample;
end
end


function Result = run_single_blade_step05_local( ...
    case_name, case_dir, case_data, step02_metadata, ...
    opr_times, F_omega_deg_s, speed_diag, ...
    Template, passes_by_sensor, WindowPlan, cfg, S, blade_id, ...
    opr_events_per_revolution)

method = build_step05_method_metadata_local(S);

trend_rows = repmat(make_empty_trend_row_local(), height(WindowPlan), 1);
WindowResult = repmat(make_empty_window_result_local(), height(WindowPlan), 1);
FinalBundleStore = cell(height(WindowPlan), 1);
FinalSeedTableStore = cell(height(WindowPlan), 1);

best = struct('weighted_voltage_rmse', inf, 'decision_score', inf);
prev_phase_result = [];

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
        Template, passes_by_sensor, cfg, S, blade_id, W, ...
        opr_events_per_revolution);

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
    core_coverage = build_bundle_coverage_local(core_bundle, S);
    core_bundle.Coverage = core_coverage;

    if core_bundle.point_count < S.min_fit_points
        warning('Window %d has too few query-safe core points (%d).', ...
            iw, core_bundle.point_count);

        [trend_rows(iw), WindowResult(iw)] = make_failed_window_result_local( ...
            W, iw, bundle_all, 'too_few_query_safe_core_points');
        trend_rows(iw) = attach_coverage_to_trend_row_local( ...
            trend_rows(iw), core_coverage);
        WindowResult(iw).CoreCoverage = core_coverage;
        WindowResult(iw).CoreBundlePreview = ...
            downsample_bundle_for_storage_local(core_bundle, ...
            S.store_bundle_preview_points);
        continue;
    end

    if ~core_coverage.ok
        warning(['Window %d has insufficient valid core waveform coverage: %s ' ...
            '(sensor-pass %d/%d, sensor %d/%d, lap %d/%d).'], ...
            iw, core_coverage.failure_reason, ...
            core_coverage.valid_pass_count, core_coverage.min_valid_pass_count, ...
            core_coverage.valid_sensor_count, core_coverage.min_valid_sensor_count, ...
            core_coverage.valid_lap_count, core_coverage.min_valid_lap_count);

        [trend_rows(iw), WindowResult(iw)] = make_failed_window_result_local( ...
            W, iw, bundle_all, core_coverage.failure_reason);
        trend_rows(iw) = attach_coverage_to_trend_row_local( ...
            trend_rows(iw), core_coverage);
        WindowResult(iw).CoreCoverage = core_coverage;
        WindowResult(iw).CoreBundlePreview = ...
            downsample_bundle_for_storage_local(core_bundle, ...
            S.store_bundle_preview_points);
        continue;
    end

    static_eta = get_static_sensor_eta_local(core_bundle, S);
    core_bundle.static_sensor_eta = static_eta;
    core_bundle.static_eta_source = S.static_eta_source;

    [core_vp_bundle, core_vp_info] = build_vp_screening_bundle_local( ...
        bundle_all, bundle_all.core_mask, core_bundle, S, 'core_query_safe');

    seed_table_core = solve_vp_seed_eo_scan_local(core_vp_bundle, eo_candidates, S);
    [final_eos_core, eo_selection_core] = ...
        select_eo_candidates_adaptive_local(seed_table_core, S);

    fit_core = refine_direct_template_fit_local( ...
        core_bundle, seed_table_core, final_eos_core, S);

    fit_core.fit_stage = 'core_query_safe';
    fit_core.core_point_count = core_bundle.point_count;
    fit_core.expanded_point_count = core_bundle.point_count;
    fit_core.vp_point_count = core_vp_bundle.point_count;
    fit_core.vp_screening_policy = core_vp_info.policy;
    fit_core.vp_screening_fallback = core_vp_info.fallback_reason;
    fit_core.EOSelectionInfo = eo_selection_core;

    % ---------------------------------------------------------------------
    % Second pass: deterministic phase-safe expansion.
    % Expanded points are selected by fitted x_query, not by fallback.
    % ---------------------------------------------------------------------
    expanded_bundle = [];
    expanded_vp_bundle = [];
    seed_table_expanded = [];
    fit_expanded = [];
    phase_ref = [];
    expanded_vp_info = make_empty_vp_screening_info_local('phase_safe_expanded');

    expand_info = struct();
    expand_info.applied = false;
    expand_info.reason = 'disabled_or_not_applicable';
    expand_info.core_point_count = core_bundle.point_count;
    expand_info.core_vp_point_count = core_vp_bundle.point_count;
    expand_info.expanded_point_count = NaN;
    expand_info.expanded_vp_point_count = NaN;
    expand_info.chosen_pass = 'core_query_safe';

    if S.phase_safe_expansion && strcmpi(fit_core.status, 'ok')
        [phase_ref, phase_ref_info] = resolve_phase_safe_reference_local( ...
            prev_phase_result, fit_core, seed_table_core, iw, S);

        expanded_mask = build_phase_safe_expansion_mask_local( ...
            bundle_all, phase_ref, S);

        expanded_bundle = select_bundle_points_local(bundle_all, expanded_mask, S);
        expanded_bundle.static_sensor_eta = static_eta;
        expanded_bundle.static_eta_source = S.static_eta_source;
        expanded_coverage = build_bundle_coverage_local(expanded_bundle, S);
        expanded_bundle.Coverage = expanded_coverage;

        expand_info.applied = true;
        expand_info.expanded_point_count = expanded_bundle.point_count;
        expand_info.expanded_valid_pass_count = expanded_coverage.valid_pass_count;
        expand_info.expanded_valid_sensor_count = expanded_coverage.valid_sensor_count;
        expand_info.expanded_valid_lap_count = expanded_coverage.valid_lap_count;
        expand_info.reason = phase_ref_info.reason;
        expand_info.reference_source = phase_ref_info.source;
        expand_info.reference_window_id = phase_ref_info.window_id;

        if expanded_bundle.point_count >= S.min_fit_points && ...
                expanded_coverage.ok

            [expanded_vp_bundle, expanded_vp_info] = build_vp_screening_bundle_local( ...
                bundle_all, expanded_mask, expanded_bundle, S, 'phase_safe_expanded');

            seed_table_expanded = solve_vp_seed_eo_scan_local( ...
                expanded_vp_bundle, eo_candidates, S);

            [final_eos_expanded, eo_selection_expanded] = ...
                select_eo_candidates_adaptive_local(seed_table_expanded, S);
            expand_info.expanded_eo_policy = 'phase_reference_selects_points_only_adaptive_full_fit';
            expand_info.expanded_eo_candidates = final_eos_expanded(:).';
            expand_info.ExpandedEOSelectionInfo = eo_selection_expanded;
            expand_info.expanded_vp_point_count = expanded_vp_bundle.point_count;
            expand_info.expanded_vp_policy = expanded_vp_info.policy;
            expand_info.expanded_vp_fallback = expanded_vp_info.fallback_reason;

            fit_expanded = refine_direct_template_fit_local( ...
                expanded_bundle, seed_table_expanded, final_eos_expanded, S);

            fit_expanded.fit_stage = 'phase_safe_expanded';
            fit_expanded.core_point_count = core_bundle.point_count;
            fit_expanded.expanded_point_count = expanded_bundle.point_count;
            fit_expanded.vp_point_count = expanded_vp_bundle.point_count;
            fit_expanded.vp_screening_policy = expanded_vp_info.policy;
            fit_expanded.vp_screening_fallback = expanded_vp_info.fallback_reason;
            fit_expanded.EOSelectionInfo = eo_selection_expanded;
        else
            expand_info.reason = 'phase_safe_insufficient_coverage';
        end
    else
        expanded_coverage = [];
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
    result.CoreEO = getfield_default_local(fit_core, 'EO_id', NaN);
    result.ExpandedEO = getfield_default_local(fit_expanded, 'EO_id', NaN);
    result.PhaseReferenceEO = getfield_default_local(phase_ref, 'EO_id', NaN);
    result.FinalEOChangedAfterExpansion = ...
        isfinite(result.CoreEO) && isfinite(result.ExpandedEO) && ...
        result.ExpandedEO ~= result.CoreEO;
    if strcmpi(result.fit_stage, 'phase_safe_expanded') && ...
            isfield(expand_info, 'ExpandedEOSelectionInfo')
        result.EOSelectionInfo = expand_info.ExpandedEOSelectionInfo;
    else
        result.EOSelectionInfo = eo_selection_core;
    end
    result.CoreVPPointCount = core_vp_bundle.point_count;
    result.CoreVPPolicy = string(core_vp_info.policy);
    result.CoreVPFallback = string(core_vp_info.fallback_reason);
    result.ExpandedVPPointCount = get_bundle_point_count_local(expanded_vp_bundle);
    result.ExpandedVPPolicy = string(getfield_default_local( ...
        expanded_vp_info, 'policy', ''));
    result.ExpandedVPFallback = string(getfield_default_local( ...
        expanded_vp_info, 'fallback_reason', ''));
    if strcmpi(result.fit_stage, 'phase_safe_expanded') && ~isempty(expanded_vp_bundle)
        result.VPPointCount = expanded_vp_bundle.point_count;
        result.VPPolicy = string(expanded_vp_info.policy);
        result.VPFallback = string(expanded_vp_info.fallback_reason);
    else
        result.VPPointCount = core_vp_bundle.point_count;
        result.VPPolicy = string(core_vp_info.policy);
        result.VPFallback = string(core_vp_info.fallback_reason);
    end
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
    WindowResult(iw).EOSelectionInfo = eo_selection_core;
    WindowResult(iw).CoreCoverage = core_coverage;
    WindowResult(iw).ExpandedCoverage = expanded_coverage;

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
    FinalBundleStore{iw} = final_bundle;
    FinalSeedTableStore{iw} = final_seed_table;

    if isfield(result, 'status') && strcmpi(result.status, 'ok')
        prev_phase_result = result;
    else
        prev_phase_result = [];
    end

    result_best_score = getfield_default_local(result, ...
        'decision_score', result.weighted_voltage_rmse);
    current_best_score = getfield_default_local(best, ...
        'decision_score', best.weighted_voltage_rmse);

    if strcmpi(result.status, 'ok') && isfinite(result_best_score) && ...
            (result_best_score < current_best_score || ...
             (abs(result_best_score - current_best_score) <= eps(max(abs(current_best_score), 1)) && ...
              result.weighted_voltage_rmse < best.weighted_voltage_rmse))

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

sensor_tag = ['S', sprintf('%d', S.analysis_sensors)];
for iw = 1:numel(trend_rows)
    trend_rows(iw).case_name = string(case_name);
    trend_rows(iw).blade_id = blade_id;
    trend_rows(iw).sensor_tag = string(sensor_tag);
end

JointEOSummary = build_joint_eo_summary_local(WindowResult, S);
JointReportedRefit = make_empty_joint_reported_refit_summary_local();

if S.joint_reported_refit_enable
    [trend_rows, WindowResult, best, JointReportedRefit] = ...
        apply_joint_reported_eo_refit_local( ...
        trend_rows, WindowResult, FinalBundleStore, FinalSeedTableStore, ...
        WindowPlan, JointEOSummary, best, S);

    for iw = 1:numel(trend_rows)
        trend_rows(iw).case_name = string(case_name);
        trend_rows(iw).blade_id = blade_id;
        trend_rows(iw).sensor_tag = string(sensor_tag);
    end
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
Result.JointEOSummary = JointEOSummary;
Result.JointReportedRefit = JointReportedRefit;
Result.ResonanceSummary = build_resonance_summary_local(Trend, Result.JointEOSummary);
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
method.name = S.method_name;
method.version = S.version;
method.outer_workflow = 'Step5_SingleSync style sliding-window single-synchronous identification';
method.solver = 'OPRCenterStd low-speed non-parametric direct-template waveform inversion';
if isfield(S, 'active_eta_model') && strcmpi(S.active_eta_model, 'none')
    method.forward_model = 'V = T_{s,b}(x - dx_c - A*sin(EO*theta + phi))';
else
    method.forward_model = 'V = T_{s,b}(x - dx_c - eta_s - A*sin(EO*theta + phi))';
end
method.objective_mode = S.objective_mode;
method.core_mask_mode = S.core_mask_mode;
method.query_guard_mode = S.query_guard_mode;
method.vp_bundle_mode = S.vp_bundle_mode;
method.vp_rank_score_mode = getfield_default_local( ...
    S, 'vp_rank_score_mode', 'linear_residual');
method.vp_fit_bundle_split = ['VP scans all EO on the same eligible main-pulse/query-safe bundle; ' ...
    'adaptive EO candidates enter full nonlinear waveform refinement'];
method.eo_low_gap_ratio_threshold = S.eo_low_gap_ratio_threshold;
if S.vp_include_sensor_eta
    method.vp_seed_model = 'q=-(V-T)/Tprime = dx_c + eta_s + a*sin(EO*theta) + b*cos(EO*theta)';
elseif ~isfield(S, 'active_eta_model') || strcmpi(S.active_eta_model, 'none')
    method.vp_seed_model = 'q=-(V-T)/Tprime = dx_c + a*sin(EO*theta) + b*cos(EO*theta)';
else
    method.vp_seed_model = 'q=-(V-T)/Tprime - eta_s = dx_c + a*sin(EO*theta) + b*cos(EO*theta)';
end
method.eo_selection = sprintf(['adaptive_full_fit: VP all EO, retain top-%d minimum, ' ...
    'always protect top-rank-%d, near-best gap <= %.3g, neighbors radius=%d, capped at %d full-fit EO'], ...
    S.vp_top_k_min, getfield_default_local(S, 'vp_top_rank_keep', S.vp_top_k_min), ...
    S.vp_gap_ratio_keep, S.neighbor_eo_radius, ...
    max(S.vp_top_k_max, getfield_default_local(S, 'vp_top_rank_keep', S.vp_top_k_max)));
method.eo_decision_score_mode = getfield_default_local( ...
    S, 'eo_decision_score_mode', 'weighted_voltage_rmse');
method.frequency_decision = ['single-window EO uses gradient-displacement VP seeds, ' ...
    'then selects EO by full-waveform voltage RMSE refinement; static inversion residual is not the formal EO decision score'];
method.frequency_output_policy = ['ResonanceSummary.reported_eo uses the joint EO summary first; ' ...
    'Trend mode is retained as a diagnostic fallback'];
method.parameter_constraints = sprintf(['A limit=%s, dx_c limit=%s, eta limit=%s, ' ...
    'eta regularization=%g, overshoot penalty=%g'], ...
    numeric_limit_text_local(S.amplitude_limit_mm), ...
    numeric_limit_text_local(S.dx_c_limit_mm), ...
    numeric_limit_text_local(S.sensor_eta_limit_mm), ...
    S.sensor_eta_reg_weight_v_per_mm, S.overshoot_penalty_weight);
method.template_domain_validity_protection = sprintf([ ...
    'invalid query penalty scale=%g; this protects template-domain coverage ' ...
    'and is not an EO prior or a physical parameter bound'], ...
    S.invalid_query_penalty_scale);
method.eta_model_policy = S.eta_model_policy;
method.eta_models = S.eta_models;
method.static_eta_source = getfield_default_local(S, 'static_eta_source', '');
method.sensor_eta_prior_mm = getfield_default_local(S, 'sensor_eta_prior_mm', []);
if isfield(S, 'active_eta_model') && strcmpi(S.active_eta_model, 'none')
    method.eta_gauge = sprintf('formal model is NoEta; eta_s is zero for all CH, reference CH%d retained for diagnostics', S.analysis_sensors(1));
else
    method.eta_gauge = sprintf('eta_s fixed to zero for reference CH%d; other CH use external static eta prior', S.analysis_sensors(1));
end
method.opr_events_per_revolution_policy = 'explicit dataset fact resolved from Template/Step02/OPRTable fallback';
method.opr_timing_reference_policy = 'Step04 Template and Step02 OPR reference time must use the same configured OPR reference';
method.angle_map_mode = S.angle_map_mode;
method.query_safe_policy = 'official route requires Step04 x_query_safe_domain; invalid query-safe domain is a template-contract error';
method.core_point_policy = 'finite & main_pulse & inside_domain & inside_query_guard';
method.frequency_policy = 'gradient-displacement VP all-EO scan, adaptive full nonlinear refinement';
method.parameter_bound_policy = 'formal A and dx_c bounds are enforced during objective evaluation and reporting; eta_s is not fitted in the main route';
method.phase_safe_policy = ['phase reference only selects phase-safe expanded points; ' ...
    'expanded pass reruns VP all-EO/adaptive full-waveform RMSE, so phase reference is not an EO prior'];
method.template_center_subtraction = true;
method.no_template_center_subtraction = false;
method.x_coordinate_policy = 'x_abs = OPRCenterStd mapped coordinate; x = x_abs - tpl.xc';
method.no_sensor_affine_voltage_projection = true;
method.no_within_pulse_uniform_downsampling = true;
method.core_mask_mode = S.core_mask_mode;
method.phase_safe_expansion = S.phase_safe_expansion;
method.two_stage_region_policy = ['full raw candidate waveform -> query-safe core region -> ' ...
    'optional phase-safe expanded region'];
end


function txt = numeric_limit_text_local(x)
if isfinite(x)
    txt = sprintf('%.6g', x);
else
    txt = 'Inf';
end
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


function validate_step05_template_contract_local(Template, cfg, S)
required_top = {'SensorBlade', 'Standard_Relative_Angles', ...
    'SchemaVersion', 'Final_Template_Point_Policy'};
for i = 1:numel(required_top)
    if ~isfield(Template, required_top{i})
        error('Step05:TemplateContract', ...
            'Step04 Template missing required field: %s.', required_top{i});
    end
end

if ~strcmp(string(Template.SchemaVersion), "LowSpeedOPRCenterStdTemplate/v1")
    error('Step05:TemplateContract', ...
        ['Step04 Template.SchemaVersion must be ' ...
         'LowSpeedOPRCenterStdTemplate/v1, but got %s.'], ...
        char(string(Template.SchemaVersion)));
end

if ~strcmp(string(Template.Final_Template_Point_Policy), ...
        "trusted_domain_recalibration")
    error('Step05:TemplateContract', ...
        ['Step04 Template.Final_Template_Point_Policy must be ' ...
         'trusted_domain_recalibration, but got %s.'], ...
        char(string(Template.Final_Template_Point_Policy)));
end

if ~isfield(Template, 'OPR_Events_Per_Revolution') && ...
        ~isfield(Template, 'opr_events_per_revolution')
    error('Step05:TemplateContract', ...
        'Step04 Template missing OPR events/revolution metadata.');
end

template_epr = resolve_step05_opr_events_per_revolution_local( ...
    cfg, Template, table(), struct());
if strcmp(string(cfg.dataset), "20251222") && template_epr ~= 6
    error('Step05:OPRSemantics', ...
        '20251222 should use 6 OPR events/rev, but Template has %g.', ...
        template_epr);
end
if ~isfield(Template, 'PointCloud_Settings') || ...
        ~isfield(Template.PointCloud_Settings, 'angle_map_mode') || ...
        ~strcmpi(Template.PointCloud_Settings.angle_map_mode, S.angle_map_mode)
    error('Step05:TemplateContract', ...
        'Step04 Template does not use the formal unequal-slot OPR angle mapping.');
end
if ~isfield(Template.PointCloud_Settings, 'slot_angle_deg') || ...
        max(abs(double(Template.PointCloud_Settings.slot_angle_deg(:).') - ...
        S.opr_slot_angle_deg)) > 1e-10
    error('Step05:TemplateContract', ...
        'Step04 Template and Step05 OPR6 slot-angle tables differ.');
end

angle_size = size(Template.Standard_Relative_Angles);
if numel(angle_size) < 2 || angle_size(2) < max(S.target_blades)
    error('Step05:TemplateContract', ...
        'Template.Standard_Relative_Angles does not cover requested sensors/blades.');
end

for sid = S.analysis_sensors
    for blade_id = S.target_blades
        get_standard_angle_local(Template, sid, blade_id);

        tpl = get_template_entry_local(Template, sid, blade_id);
        if isempty(tpl)
            error('Step05:TemplateContract', ...
                'Missing Step04 template for CH%d Blade%d.', sid, blade_id);
        end

        if ~isfield(tpl, 'x_domain') || numel(tpl.x_domain) ~= 2 || ...
                any(~isfinite(tpl.x_domain)) || tpl.x_domain(2) <= tpl.x_domain(1)
            error('Step05:TemplateContract', ...
                'Template CH%d Blade%d has invalid x_domain.', sid, blade_id);
        end

        xc = get_template_xcenter_local(tpl);
        if ~isfinite(xc)
            error('Step05:TemplateContract', ...
                'Template CH%d Blade%d has invalid center xc.', sid, blade_id);
        end

        if ~isfield(tpl, 'x_query_safe_domain') || ...
                numel(tpl.x_query_safe_domain) ~= 2 || ...
                any(~isfinite(tpl.x_query_safe_domain)) || ...
                tpl.x_query_safe_domain(2) <= tpl.x_query_safe_domain(1)
            error('Step05:TemplateContract', ...
                ['Official Step05 requires valid Step04 x_query_safe_domain ' ...
                 'for CH%d Blade%d when query_guard_mode=step04. ' ...
                 'Re-run unified Step04 before identification.'], ...
                sid, blade_id);
        end
    end
end
end


function epr = resolve_step05_opr_events_per_revolution_local( ...
    cfg, Template, OPRTable, step02_metadata)
epr = NaN;

if isstruct(Template)
    names = {'OPR_Events_Per_Revolution', 'opr_events_per_revolution', ...
        'OPR_Pulses_Per_Rev', 'opr_pulses_per_rev'};
    epr = first_scalar_field_value_local(Template, names);
end

if ~isfinite(epr) && isstruct(step02_metadata)
    names = {'OPR_Events_Per_Revolution', 'opr_events_per_revolution', ...
        'OPREventsPerRevolution', 'OPR_Pulses_Per_Rev', 'opr_pulses_per_rev'};
    epr = first_scalar_field_value_local(step02_metadata, names);
end

if ~isfinite(epr) && istable(OPRTable) && ...
        ismember('opr_events_per_revolution', OPRTable.Properties.VariableNames)
    v = OPRTable.opr_events_per_revolution;
    v = v(isfinite(v));
    if ~isempty(v)
        epr = mode(v);
    end
end

if ~isfinite(epr)
    switch string(cfg.dataset)
        case "20241106"
            epr = 1;
        case {"20250527", "20251222"}
            epr = 6;
        otherwise
            error('Cannot resolve OPR events per revolution for dataset %s.', ...
                cfg.dataset);
    end
end

if epr <= 0 || abs(epr - round(epr)) > 1e-9
    error('Invalid OPR events per revolution: %g.', epr);
end

epr = round(epr);
end


function v = first_scalar_field_value_local(s, names)
v = NaN;
for i = 1:numel(names)
    name = names{i};
    if isfield(s, name) && isnumeric(s.(name)) && ...
            isscalar(s.(name)) && isfinite(s.(name))
        v = s.(name);
        return;
    end
end
end

function opr_times = resolve_step05_opr_reference_times_local(OPRTable, cfg, Template)
if ismember('opr_reference_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_reference_time_s;
    source = "opr_reference_time_s";
    template_method = getfield_default_local(Template, 'OPR_Timing_Method', '');
    if ~isempty(template_method)
        method = lower(strtrim(string(template_method)));
        expects_rising = any(method == ["threshold_rising_edge", "rising_edge", "start_edge"]);
        if expects_rising && ismember('opr_center_time_s', OPRTable.Properties.VariableNames) && ...
                isequal(opr_times, OPRTable.opr_center_time_s)
            error('Step05:OPRReferenceMismatch', ...
                ['Template uses %s, but Step02 opr_reference_time_s equals center time. ' ...
                 'Re-run Step02 with cfg.opr_timing_reference=''rising_edge''.'], ...
                char(string(template_method)));
        end
    end
    return;
end

ref = lower(strtrim(string(getfield_default_local(cfg, ...
    'opr_timing_reference', 'rising_edge'))));
if any(ref == ["rising_edge", "start_edge", "start", "threshold_rising_edge"]) && ...
        ismember('opr_start_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_start_time_s;
    source = "opr_start_time_s";
elseif any(ref == ["center", "opr_center", "pulse_center", "multi_threshold_width_center"]) && ...
        ismember('opr_center_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_center_time_s;
    source = "opr_center_time_s";
elseif any(ref == ["falling_edge", "end_edge", "end"]) && ...
        ismember('opr_end_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_end_time_s;
    source = "opr_end_time_s";
elseif ismember('opr_start_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_start_time_s;
    source = "opr_start_time_s";
else
    opr_times = OPRTable.opr_center_time_s;
    source = "opr_center_time_s";
end

template_method = getfield_default_local(Template, 'OPR_Timing_Method', '');
if ~isempty(template_method)
    method = lower(strtrim(string(template_method)));
    expects_rising = any(method == ["threshold_rising_edge", "rising_edge", "start_edge"]);
    if expects_rising && source == "opr_center_time_s"
        error('Step05:OPRReferenceMismatch', ...
            ['Template uses %s, but Step02 OPRTable has no usable rising/reference time column. ' ...
             'Re-run Step02 with cfg.opr_timing_reference=''rising_edge''.'], ...
            char(string(template_method)));
    end
end
end


function [F_omega_deg_s, diagnostic] = build_dynamic_speed_interpolant_local( ...
    opr_times, opr_events_per_revolution)
opr_times = opr_times(:);
epr = opr_events_per_revolution;

rev_period = opr_times((epr + 1):end) - opr_times(1:(end - epr));

speed_time = 0.5 .* (opr_times((epr + 1):end) + opr_times(1:(end - epr)));
speed_deg_s = 360 ./ max(rev_period, eps);

valid = isfinite(speed_time) & isfinite(speed_deg_s) & speed_deg_s > 0;

if nnz(valid) < 2
    error('Not enough valid OPR periods for dynamic speed interpolation.');
end

F_omega_deg_s = griddedInterpolant( ...
    speed_time(valid), speed_deg_s(valid), ...
    'linear', 'nearest');

diagnostic = struct();
diagnostic.opr_events_per_revolution = epr;
diagnostic.speed_time_s = speed_time(valid);
diagnostic.speed_deg_s = speed_deg_s(valid);
diagnostic.rpm = speed_deg_s(valid) / 360 * 60;
diagnostic.rev_period_s = rev_period(valid);
diagnostic.median_rpm = median(diagnostic.rpm, 'omitnan');
end


function passes_by_sensor = build_analysis_passes_by_sensor_local(case_data, S, blade_id)
passes_by_sensor = struct();

if ~isfield(case_data, 'ObservationTable')
    error('20251222 Step05 case_data must contain ObservationTable.');
end

Obs = case_data.ObservationTable;

for sid = S.analysis_sensors
    mask = Obs.sensor_id == sid & ...
           Obs.blade_id == blade_id & ...
           Obs.is_valid & ...
           isfinite(Obs.arrival_time_s) & ...
           isfinite(Obs.pulse_start_time_s) & ...
           isfinite(Obs.pulse_end_time_s) & ...
           isfinite(Obs.source_file_id);

    rows = Obs(mask, :);

    [~, order] = sort(rows.arrival_time_s);
    rows = rows(order, :);

    % Use the same target-blade passes after the same analysis start time for
    % every sensor. Lap indexing must be reset after this filter.
    rows = rows(rows.arrival_time_s >= S.analysis_start_time_s, :);

    if isempty(rows)
        warning('No Blade%d passes for CH%d after %.6f s.', ...
            blade_id, sid, S.analysis_start_time_s);
        continue;
    end

    if isfinite(S.target_laps) && height(rows) > S.target_laps
        rows = rows(1:S.target_laps, :);
    end

    analysis_lap = (1:height(rows)).';

    T = table( ...
        analysis_lap, ...
        rows.pulse_start_time_s, ...
        rows.pulse_end_time_s, ...
        rows.arrival_time_s, ...
        rows.blade_id, ...
        rows.source_file_id, ...
        rows.rev_id, ...
        rows.pulse_index, ...
        'VariableNames', {'analysis_lap','start_s','end_s','arrival_s', ...
        'blade_id','source_file_id','rev_id','pulse_index'});

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
    Template, passes_by_sensor, cfg, S, blade_id, W, ...
    opr_events_per_revolution)

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
            p.start_s, p.end_s, p.arrival_s, p.source_file_id, cfg, S);

        if numel(t_seg) < 5
            continue;
        end

        % -------------------------------------------------------------
        % 2. OPRCenterStd mapping.
        %    Every sample in one pulse segment uses the OPR immediately
        %    before the blade passage time as angular reference.
        % -------------------------------------------------------------
        theta_std_deg = get_standard_angle_local(Template, sid, blade_id);

        [x_mm, theta_rad, prev_opr_idx, keep] = ...
            map_dynamic_times_to_oprcenterstd_local( ...
            t_seg, opr_times, F_omega_deg_s, ...
            theta_std_deg, cfg, p.arrival_s, opr_events_per_revolution);

        t_seg = t_seg(keep);
        v_seg = v_seg(keep);
        x_abs_mm = x_mm(keep);
        x_center_mm = get_template_xcenter_local(tpl);
        x_mm = x_abs_mm - x_center_mm;
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

        sensor_threshold = get_step05_sensor_threshold_local(cfg, sid, tpl, S);
        main_pulse_mask = build_main_pulse_mask_local( ...
            t_seg, v_seg, S, sensor_threshold);

        % -------------------------------------------------------------
        % 5. Query-safe shrink.
        %    Since final query coordinate is:
        %       x_query = x - dx_c - A*sin(...)
        %    the first/core pass must use a smaller region than x_domain.
        % -------------------------------------------------------------
        base_for_guard = finite_mask & inside_domain;
        if any(strcmpi(S.core_mask_mode, {'legacy_base', 'legacy_base_query_safe'}))
            base_for_guard = base_for_guard & dynamic_effective & main_pulse_mask;
        end

        [query_safe_domain, query_guard, query_guard_source] = ...
            resolve_query_safe_domain_local(x_mm, v_seg, tpl, S, base_for_guard);

        inside_guard = x_mm >= query_safe_domain(1) & ...
                       x_mm <= query_safe_domain(2);

        [fit_weight, weight_parts] = build_step05_fit_weight_local( ...
            t_seg, v_seg, x_mm, tpl, template_weight, query_guard, S);

        switch lower(strtrim(S.core_mask_mode))
            case 'legacy_base'
                % Historical soft-domain route: query guard contributes to
                % weights, but does not hard-remove valid pulse points.
                core_mask = finite_mask & ...
                            inside_domain & ...
                            dynamic_effective & ...
                            main_pulse_mask;

            case 'legacy_base_query_safe'
                % Historical base mask from the original direct-template
                % route, plus query-safe shrink.
                base_mask = finite_mask & ...
                            inside_domain & ...
                            dynamic_effective & ...
                            main_pulse_mask;

                core_mask = base_mask & inside_guard;

            case 'sg_threshold_mainpulse_trust'
                % SG Step5-compatible core contract:
                % threshold-isolated main pulse inside the Step04 template
                % trust domain. Query-safe shrink is a diagnostic here; it
                % must not change the EO objective relative to the SG route.
                core_mask = finite_mask & inside_domain & main_pulse_mask;

            otherwise
                % Main deterministic route for the paper:
                % complete per-pass main pulse inside the query-safe
                % shrunken trusted domain.
                core_mask = finite_mask & inside_domain & inside_guard & ...
                            main_pulse_mask;
        end

        n = numel(x_mm);

        bundle.x = [bundle.x; x_mm(:)]; %#ok<AGROW>
        bundle.x_abs = [bundle.x_abs; x_abs_mm(:)]; %#ok<AGROW>
        bundle.x_center = [bundle.x_center; repmat(x_center_mm, n, 1)]; %#ok<AGROW>
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
        bundle.edge_weight = [bundle.edge_weight; weight_parts.edge(:)]; %#ok<AGROW>
        bundle.domain_weight = [bundle.domain_weight; weight_parts.domain(:)]; %#ok<AGROW>
        bundle.query_weight = [bundle.query_weight; weight_parts.query(:)]; %#ok<AGROW>
        bundle.gradient_fit_weight = [bundle.gradient_fit_weight; weight_parts.gradient(:)]; %#ok<AGROW>
        bundle.fit_weight = [bundle.fit_weight; fit_weight(:)]; %#ok<AGROW>

        bundle.inside_domain_mask = [bundle.inside_domain_mask; inside_domain(:)]; %#ok<AGROW>
        bundle.inside_guard_mask = [bundle.inside_guard_mask; inside_guard(:)]; %#ok<AGROW>
        bundle.dynamic_effective_mask = [bundle.dynamic_effective_mask; dynamic_effective(:)]; %#ok<AGROW>
        bundle.main_pulse_mask = [bundle.main_pulse_mask; main_pulse_mask(:)]; %#ok<AGROW>
        bundle.finite_mask = [bundle.finite_mask; finite_mask(:)]; %#ok<AGROW>
        bundle.core_mask = [bundle.core_mask; core_mask(:)]; %#ok<AGROW>
        bundle.query_guard_mm = [bundle.query_guard_mm; repmat(query_guard, n, 1)]; %#ok<AGROW>
        bundle.query_guard_source = [bundle.query_guard_source; repmat(string(query_guard_source), n, 1)]; %#ok<AGROW>
        bundle.query_safe_left_mm = [bundle.query_safe_left_mm; repmat(query_safe_domain(1), n, 1)]; %#ok<AGROW>
        bundle.query_safe_right_mm = [bundle.query_safe_right_mm; repmat(query_safe_domain(2), n, 1)]; %#ok<AGROW>

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
bundle.x_abs = [];
bundle.x_center = [];
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
bundle.edge_weight = [];
bundle.domain_weight = [];
bundle.query_weight = [];
bundle.gradient_fit_weight = [];
bundle.fit_weight = [];

bundle.inside_domain_mask = [];
bundle.inside_guard_mask = [];
bundle.dynamic_effective_mask = [];
bundle.main_pulse_mask = [];
bundle.finite_mask = [];
bundle.core_mask = [];
bundle.query_guard_mm = [];
bundle.query_guard_source = strings(0, 1);
bundle.query_safe_left_mm = [];
bundle.query_safe_right_mm = [];

bundle.point_count = 0;
bundle.core_point_count = 0;
bundle.valid_pass_count = 0;
bundle.rot_freq_mean_hz = NaN;
bundle.rot_rpm_mean = NaN;
end


function [t_seg, v_seg] = read_dynamic_pulse_segment_local( ...
    case_dir, case_data, sid, start_s, end_s, arrival_s, source_file_id, cfg, S)

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
    case_dir, case_data, sid, t0, t1, source_file_id, cfg);
end


function [t, v] = read_raw_time_window_local(case_dir, case_data, sid, t0, t1, source_file_id, cfg)
t = [];
v = [];

if ~isfield(case_data, 'file_ids') || ~isfield(case_data, 'file_index_map')
    error('Step02 case_data must contain file_ids and file_index_map. Re-run improved Step02.');
end

candidate_file_indices = 1:numel(case_data.file_ids);
if nargin >= 6 && isfinite(source_file_id)
    idx0 = find(case_data.file_ids == source_file_id, 1, 'first');
    if ~isempty(idx0)
        candidate_file_indices = max(1, idx0 - 1):min(numel(case_data.file_ids), idx0 + 1);
    end
end

for iFile = candidate_file_indices
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
persistent raw_cache raw_cache_keys
if isempty(raw_cache)
    raw_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    raw_cache_keys = {};
end

cache_key = char(file);
if isKey(raw_cache, cache_key)
    raw = raw_cache(cache_key);
    return;
end

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

raw_cache(cache_key) = raw;
raw_cache_keys{end + 1} = cache_key; %#ok<AGROW>

max_cached_files = 12;
while numel(raw_cache_keys) > max_cached_files
    old_key = raw_cache_keys{1};
    raw_cache_keys(1) = [];
    if isKey(raw_cache, old_key)
        remove(raw_cache, old_key);
    end
end
end


function calibration = load_opr_slot_calibration_step05_local(cfg)
file = cfg.opr_slot_calibration_file;
if ~isfile(file)
    error('Step05:MissingOPRSlotCalibration', ...
        'Missing formal OPR6 slot calibration: %s', file);
end
S = load(file, 'SlotAngleCalibration');
if ~isfield(S, 'SlotAngleCalibration')
    error('Step05:InvalidOPRSlotCalibration', ...
        'Calibration file lacks SlotAngleCalibration: %s', file);
end
calibration = S.SlotAngleCalibration;
valid = isfield(calibration, 'SchemaVersion') && ...
    strcmpi(calibration.SchemaVersion, 'OPR6_UNEQUAL_SLOT_ANGLE_V1') && ...
    calibration.EventsPerRevolution == 6 && ...
    numel(calibration.IntervalAngleDeg) == 6 && ...
    all(isfinite(calibration.IntervalAngleDeg)) && ...
    abs(sum(calibration.IntervalAngleDeg) - 360) < 1e-9 && ...
    numel(calibration.AnchorAngleDeg) == 6 && ...
    numel(calibration.PulseWidthAngleDeg) == 6;
if ~valid
    error('Step05:InvalidOPRSlotCalibration', ...
        'Invalid formal OPR6 slot calibration contract: %s', file);
end
end


function [bestOffset, Audit] = align_dynamic_opr_slots_local( ...
    OPRTable, referencePulseWidthAngleDeg, epr)
required = {'event_in_rev','opr_width_s','rev_period_s'};
assert(all(ismember(required, OPRTable.Properties.VariableNames)), ...
    'OPRTable lacks fields required for physical slot alignment.');
eventInRev = double(OPRTable.event_in_rev(:));
widthAngleDeg = 360 .* double(OPRTable.opr_width_s(:)) ./ ...
    double(OPRTable.rev_period_s(:));
referencePulseWidthAngleDeg = double(referencePulseWidthAngleDeg(:).');
assert(numel(referencePulseWidthAngleDeg) == epr, ...
    'Reference OPR pulse-width fingerprint size mismatch.');

rows = repmat(struct('offset', NaN, 'score_deg', NaN, ...
    'score_margin_deg', NaN, 'observed_width_angle_deg', []), epr, 1);
for shift = 0:(epr - 1)
    physicalSlot = mod(eventInRev - 1 + shift, epr) + 1;
    observed = nan(1, epr);
    for slot = 1:epr
        m = physicalSlot == slot & isfinite(widthAngleDeg);
        observed(slot) = median(widthAngleDeg(m), 'omitnan');
    end
    rows(shift + 1).offset = shift;
    rows(shift + 1).score_deg = sqrt(mean( ...
        (observed - referencePulseWidthAngleDeg).^2, 'omitnan'));
    rows(shift + 1).observed_width_angle_deg = observed;
end
Audit = sortrows(struct2table(rows), 'score_deg', 'ascend');
Audit.score_margin_deg(:) = NaN;
if height(Audit) > 1
    Audit.score_margin_deg(1) = Audit.score_deg(2) - Audit.score_deg(1);
end
assert(isfinite(Audit.score_deg(1)) && Audit.score_deg(1) < 0.25 && ...
    Audit.score_margin_deg(1) > 0.05, ...
    'Unequal OPR physical-slot alignment is not sufficiently identifiable.');
bestOffset = Audit.offset(1);
if bestOffset > floor(epr / 2)
    bestOffset = bestOffset - epr;
end
Audit.offset(1) = bestOffset;
end


function [x_mm, theta_rad, prev_opr_idx, keep] = map_dynamic_times_to_oprcenterstd_local( ...
    t, opr_times, F_omega_deg_s, theta_std_deg, cfg, t_peak, ...
    opr_events_per_revolution)

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
idx_next = idx_ref + 1;
if idx_next > numel(opr_times)
    keep = false(size(t));
    return;
end
local_period_s = opr_times(idx_next) - t_ref;
if ~isfinite(local_period_s) || local_period_s <= 0
    keep = false(size(t));
    return;
end
event_slot = mod(idx_ref - 1 + cfg.opr_slot_index_offset, ...
    numel(cfg.opr_slot_angle_deg)) + 1;
interval_angle_deg = cfg.opr_slot_angle_deg(event_slot);
theta_actual_deg = interval_angle_deg .* (t - t_ref) ./ local_period_s;

theta_rel_deg = wrap_to_180_local(theta_actual_deg - theta_std_deg);

x_mm = theta_rel_deg .* (pi / 180) .* cfg.r_tip_mm;

theta_rad = map_time_to_rotor_phase_local( ...
    opr_times, t, opr_events_per_revolution, cfg);

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


function theta_rot = map_time_to_rotor_phase_local( ...
    opr_times, sample_times, opr_events_per_revolution, cfg)
sample_times = sample_times(:);
theta_rot = nan(size(sample_times));

epr = opr_events_per_revolution;

if numel(opr_times) <= epr
    return;
end

rev_anchor_times = opr_times(:);
event_index = (0:numel(rev_anchor_times)-1).';
event_slot = mod(event_index + cfg.opr_slot_index_offset, epr) + 1;
anchor_phase = deg2rad(cfg.opr_slot_anchor_deg(event_slot));
rev_phase = unwrap(anchor_phase(:));

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
interp_method = getfield_default_local(S, 'template_interp_method', 'pchip');

template_v = interp1( ...
    tpl.x_grid(:), ...
    tpl.v_grid(:), ...
    x, ...
    interp_method, ...
    NaN);

template_dv_dx = interp1( ...
    tpl.x_grid(:), ...
    tpl.dv_dx(:), ...
    x, ...
    interp_method, ...
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


function [w_total, parts] = build_step05_fit_weight_local(t, v, x, tpl, template_weight, query_guard_mm, S)
% Match the original full-wave direct-template weighting: dynamic edge
% emphasis times template reliability, with optional soft domain/query and
% gradient terms.  No point is selected by this function; it only scores
% points that already passed the Step05 mask.
t = t(:);
v = v(:);
x = x(:);
template_weight = template_weight(:);

if strcmpi(S.fit_weight_mode, 'sg_edge_only')
    w_edge = build_edge_weight_step05_local(t, v, S.weight_floor);
    w_total = max(S.weight_floor, w_edge);
    if max(w_total) > 0
        w_total = max(S.weight_floor, w_total ./ max(w_total));
    end
    parts = struct( ...
        'edge', w_edge, ...
        'domain', ones(size(x)), ...
        'query', ones(size(x)), ...
        'gradient', ones(size(x)));
    return;
end

if strcmpi(S.fit_weight_mode, 'sg_edge_template')
    w_edge = build_edge_weight_step05_local(t, v, S.weight_floor);
    w_template = template_weight;
    w_template(~isfinite(w_template)) = S.weight_floor;
    w_template = max(w_template, S.weight_floor);
    if max(w_template) > 0
        w_template = max(S.weight_floor, w_template ./ max(w_template));
    end

    w_total = max(S.weight_floor, w_edge .* w_template);
    if max(w_total) > 0
        w_total = max(S.weight_floor, w_total ./ max(w_total));
    end

    parts = struct( ...
        'edge', w_edge, ...
        'domain', ones(size(x)), ...
        'query', ones(size(x)), ...
        'gradient', ones(size(x)));
    return;
end

if ~strcmpi(S.fit_weight_mode, 'legacy_product')
    w_total = template_weight;
    w_total(~isfinite(w_total)) = S.weight_floor;
    w_total = max(w_total, S.weight_floor);
    parts = struct( ...
        'edge', ones(size(x)), ...
        'domain', ones(size(x)), ...
        'query', ones(size(x)), ...
        'gradient', ones(size(x)));
    return;
end

w_edge = build_edge_weight_step05_local(t, v, S.weight_floor);
w_template = template_weight;
w_template(~isfinite(w_template)) = S.weight_floor;
w_template = max(w_template, S.weight_floor);
w_domain = build_domain_soft_weight_step05_local( ...
    x, tpl.x_domain, S.domain_soft_margin_mm, S.weight_floor);
w_query = build_query_guard_soft_weight_step05_local( ...
    x, tpl.x_domain, query_guard_mm, S.domain_soft_margin_mm, S.weight_floor);
w_gradient = build_template_gradient_weight_step05_local( ...
    tpl, x, S.domain_selection_mode, S.weight_floor);

w_total = max(S.weight_floor, w_edge .* w_template .* w_domain .* w_query .* w_gradient);
if max(w_total) > 0
    w_total = max(S.weight_floor, w_total ./ max(w_total));
end

parts = struct( ...
    'edge', w_edge, ...
    'domain', w_domain, ...
    'query', w_query, ...
    'gradient', w_gradient);
end


function w_edge = build_edge_weight_step05_local(t, v, floor_w)
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


function w_domain = build_domain_soft_weight_step05_local(x, x_domain, margin_mm, floor_w)
if margin_mm <= 0
    w_domain = ones(size(x));
    return;
end

dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
ratio = min(max(dist_to_edge ./ margin_mm, 0), 1);
w_domain = floor_w + (1 - floor_w) .* ratio;
w_domain(~isfinite(w_domain)) = floor_w;
end


function w_query = build_query_guard_soft_weight_step05_local(x, x_domain, query_guard_mm, margin_mm, floor_w)
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


function w_gradient = build_template_gradient_weight_step05_local(tpl, x, domain_selection_mode, floor_w)
if ~strcmpi(domain_selection_mode, 'soft') || ...
        ~isfield(tpl, 'dv_dx') || isempty(tpl.dv_dx)
    w_gradient = ones(size(x));
    return;
end

g = abs(interp1(tpl.x_grid(:), tpl.dv_dx(:), x(:), 'linear', 0));
g_max = max(g, [], 'omitnan');
if ~isfinite(g_max) || g_max <= 0
    w_gradient = ones(size(x));
    return;
end
w_gradient = floor_w + (1 - floor_w) .* g ./ g_max;
w_gradient(~isfinite(w_gradient)) = floor_w;
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
% Match the original 20251222 direct-template route.  The effective dynamic
% waveform is the union of informative template-gradient, time-gradient, and
% peak-voltage regions, with a low-voltage guard and fallbacks.
v = v(:);
dVdt = dVdt(:);
template_dv_dx = template_dv_dx(:);

finite = isfinite(v) & isfinite(dVdt) & isfinite(template_dv_dx);
if ~any(finite)
    mask = false(size(v));
    return;
end

max_tpl_grad = max(abs(template_dv_dx(finite)), [], 'omitnan');

if ~isfinite(max_tpl_grad) || max_tpl_grad <= 0
    tpl_grad_mask = finite;
else
    tpl_grad_mask = abs(template_dv_dx) >= ...
        S.dynamic_template_gradient_min_ratio * max_tpl_grad;
end

max_time_grad = max(abs(dVdt(finite)), [], 'omitnan');

if ~isfinite(max_time_grad) || max_time_grad <= 0
    time_grad_mask = false(size(v));
else
    time_grad_mask = abs(dVdt) >= ...
        S.dynamic_time_gradient_min_ratio * max_time_grad;
end

v_finite = v(finite);
peak_q = min(max(S.dynamic_peak_quantile, 0), 100);
peak_level = prctile(v_finite, peak_q);

if isfield(tpl, 'threshold') && ~isempty(tpl.threshold) && isfinite(tpl.threshold)
    threshold = tpl.threshold;
elseif isfield(S, 'default_sensor_threshold') && ...
        ~isempty(S.default_sensor_threshold) && isfinite(S.default_sensor_threshold)
    threshold = S.default_sensor_threshold;
else
    threshold = -inf;
end

if ~isfinite(peak_level)
    peak_mask = false(size(v));
else
    peak_mask = v >= max(threshold, peak_level);
end

mask = finite & (tpl_grad_mask | time_grad_mask | peak_mask) & ...
       v >= 0.5 * threshold;

if nnz(mask) < 8 && isfinite(threshold)
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end


function mask = build_main_pulse_mask_local(t, v, S, threshold_override)
% Diagnostic mask only in the main query-safe route.
t = t(:); %#ok<NASGU>
v = v(:);

if nargin < 4 || isempty(threshold_override) || ~isfinite(threshold_override)
    threshold_override = NaN;
end

mask = true(size(v));

if ~strcmpi(S.pulse_mode, 'single') || numel(v) < 5
    return;
end

if isfield(S, 'main_pulse_mode') && strcmpi(S.main_pulse_mode, 'legacy_threshold')
    threshold = threshold_override;
    if ~isfinite(threshold)
        threshold = S.default_sensor_threshold;
    end
    if ~isfinite(threshold)
        threshold = 0.5;
    end
    mask = isolate_legacy_main_pulse_by_threshold_local(v, threshold);
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


function mask = isolate_legacy_main_pulse_by_threshold_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);

if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax - span):min(numel(v), imax + span)) = true;
    return;
end

jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];

if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
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
        peak_val = max(v(seg), [], 'omitnan');

        if peak_val > best_peak
            best_peak = peak_val;
            best_seg = iseg;
        end
    end

    mask(idx(starts(best_seg)):idx(ends(best_seg))) = true;
end
end


function threshold = get_step05_sensor_threshold_local(cfg, sid, tpl, S)
threshold = NaN;

if isfield(tpl, 'threshold') && ~isempty(tpl.threshold) && ...
        isfinite(tpl.threshold)
    threshold = tpl.threshold;
    return;
end

if isfield(cfg, 'sensor_thresholds')
    st = cfg.sensor_thresholds;
    if isa(st, 'containers.Map')
        key_candidates = {sid, double(sid), char(string(sid)), string(sid)};
        for i = 1:numel(key_candidates)
            try
                if isKey(st, key_candidates{i})
                    threshold = st(key_candidates{i});
                    break;
                end
            catch
            end
        end
    elseif isnumeric(st) && numel(st) >= sid && isfinite(st(sid))
        threshold = st(sid);
    end
end

if ~isfinite(threshold) && isfield(cfg, 'sensor_threshold_default') && ...
        isfinite(cfg.sensor_threshold_default)
    threshold = cfg.sensor_threshold_default;
end

if ~isfinite(threshold) && isfield(S, 'default_sensor_threshold') && ...
        isfinite(S.default_sensor_threshold)
    threshold = S.default_sensor_threshold;
end

if ~isfinite(threshold)
    threshold = 0.5;
end
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


function [query_safe_domain, query_guard, source] = resolve_query_safe_domain_local(x, v, tpl, S, base_mask)
% Resolve the core query-safe domain for Step05.
%
% The 20251222 foundation Step04 stores x_query_safe_domain after checking
% low-speed coverage, response, repeatability and continuity. The official
% Step05 uses that domain by default; adaptive/fixed guards remain available
% as explicit diagnostics and for parity checks with the 20250527 route.

mode = lower(strtrim(string(S.query_guard_mode)));

if any(mode == ["step04", "step04_override", "step04_query_safe", ...
        "template_override"]) && ...
        isfield(tpl, 'x_query_safe_domain') && ...
        numel(tpl.x_query_safe_domain) == 2 && ...
        all(isfinite(tpl.x_query_safe_domain)) && ...
        tpl.x_query_safe_domain(2) > tpl.x_query_safe_domain(1)

    query_safe_domain = tpl.x_query_safe_domain(:).';
    query_safe_domain(1) = max(query_safe_domain(1), tpl.x_domain(1));
    query_safe_domain(2) = min(query_safe_domain(2), tpl.x_domain(2));
    query_guard = min([query_safe_domain(1) - tpl.x_domain(1), ...
        tpl.x_domain(2) - query_safe_domain(2)]);
    query_guard = max(query_guard, 0);
    source = "step04_x_query_safe_domain";
    return;
end

query_guard = compute_query_guard_local(x, v, tpl, S, base_mask);
query_safe_domain = [tpl.x_domain(1) + query_guard, tpl.x_domain(2) - query_guard];
if query_safe_domain(2) <= query_safe_domain(1)
    query_safe_domain = [NaN NaN];
end
source = string(S.query_guard_mode);
end


function qg = compute_query_guard_local(x, v, tpl, S, base_mask)
% Resolve query guard for the first/core pass.
%
% The guard shrinks the template trusted domain because final fitting queries:
%   x_query = x - dx_c - A*sin(EO*theta + phi)
%
% Adaptive mode estimates apparent displacement by locally inverting the
% low-speed template voltage and taking a high quantile.

if nargin < 5 || isempty(base_mask)
    base_mask = isfinite(x) & isfinite(v);
end

if strcmpi(S.query_guard_mode, 'fixed')
    if isempty(S.query_guard_mm)
        qg = S.query_guard_min_mm;
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
% Local inverse of the non-parametric template.
% Match the verified 20251222 OPRCenterStd route: find voltage crossing
% points by linear interpolation and choose the branch closest to x_ref.
% A nearest-grid inverse is too coarse for this guard estimate and can
% inflate the apparent displacement on non-monotonic pulse shoulders.

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

    vv = v(i);
    diff_v = vg - vv;
    crossing_x = [];

    exact_idx = find(abs(diff_v) <= 1e-10);
    if ~isempty(exact_idx)
        crossing_x = xg(exact_idx);
    end

    for k = 1:numel(diff_v)-1
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k+1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k+1) > 0
            continue;
        end

        denom = vg(k+1) - vg(k);
        if abs(denom) < eps
            continue;
        end

        alpha = (vv - vg(k)) / denom;
        crossing_x(end+1, 1) = xg(k) + ...
            alpha * (xg(k+1) - xg(k)); %#ok<AGROW>
    end

    if isempty(crossing_x)
        [~, idx] = min(abs(diff_v));
        x_static(i) = xg(idx);
    else
        [~, idx] = min(abs(crossing_x - x_ref(i)));
        x_static(i) = crossing_x(idx);
    end
end

x_static = reshape(x_static, size(v));
end


function x_static = invert_bundle_template_static_x_local(bundle, S)
x_static = nan(size(bundle.x(:)));
if ~isfield(bundle, 'TemplateMini') || isempty(bundle.TemplateMini)
    return;
end
sensors = S.analysis_sensors(:).';
for sid = sensors
    idx_sensor = bundle.sensor_id(:) == sid;
    if ~any(idx_sensor)
        continue;
    end
    blades = unique(bundle.blade_id_vec(idx_sensor).');
    for bid = blades
        idx = idx_sensor & bundle.blade_id_vec(:) == bid;
        if ~any(idx)
            continue;
        end
        tpl = bundle_template_from_point_local(bundle, sid, bid);
        if isempty(tpl)
            continue;
        end
        mode = lower(strtrim(getfield_default_local( ...
            S, 'vp_template_inverse_mode', 'branch_monotone')));
        switch mode
            case {'branch_monotone', 'monotone_branch'}
                x_static(idx) = invert_template_branch_monotone_local( ...
                    tpl, bundle.v(idx), bundle.x(idx));
            case {'nearest_crossing', 'crossing'}
                x_static(idx) = invert_template_voltage_for_guard_local( ...
                    tpl, bundle.v(idx), bundle.x(idx));
            otherwise
                error('Unsupported S.vp_template_inverse_mode: %s', mode);
        end
    end
end
end


function x_static = invert_template_branch_monotone_local(tpl, v_obs, x_ref)
x_static = nan(size(v_obs(:)));
side = sign(x_ref(:));
side(side == 0) = 1;
left = side < 0;
right = side >= 0;
x_static(left) = inverse_template_one_side_monotone_local(tpl, v_obs(left), -1);
x_static(right) = inverse_template_one_side_monotone_local(tpl, v_obs(right), 1);
end


function xq = inverse_template_one_side_monotone_local(tpl, v_obs, side_sign)
xq = nan(size(v_obs(:)));
if isempty(v_obs)
    return;
end
xg = tpl.x_grid(:);
vg = tpl.v_grid(:);
valid = isfinite(xg) & isfinite(vg);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && ...
        all(isfinite(tpl.x_domain(1:2)))
    valid = valid & xg >= tpl.x_domain(1) & xg <= tpl.x_domain(2);
end
if side_sign < 0
    valid = valid & xg <= 0;
    [x_path, order] = sort(xg(valid), 'ascend');
else
    valid = valid & xg >= 0;
    [x_path, order] = sort(xg(valid), 'descend');
end
vg_side = vg(valid);
v_path = vg_side(order);
good = isfinite(x_path) & isfinite(v_path);
x_path = x_path(good);
v_path = v_path(good);
if numel(x_path) < 3
    return;
end
v_mono = cummax(v_path);
[v_unique, ia] = unique(v_mono, 'stable');
x_unique = x_path(ia);
if numel(v_unique) < 2 || v_unique(end) <= v_unique(1)
    [~, nearest_idx] = min(abs(v_path(:).' - v_obs(:)), [], 2);
    xq(:) = x_path(nearest_idx);
    return;
end
v_clip = min(max(v_obs(:), v_unique(1)), v_unique(end));
xq(:) = interp1(v_unique, x_unique, v_clip, 'linear', 'extrap');
end


function eta = get_fixed_joint_static_eta_local(S)
eta = zeros(1, numel(S.analysis_sensors));
if ~isfield(S, 'sensor_eta_prior_mm') || ...
        numel(S.sensor_eta_prior_mm) ~= numel(S.analysis_sensors)
    error('Fixed joint static eta prior is missing or has the wrong length.');
end
eta = double(S.sensor_eta_prior_mm(:).');
if any(~isfinite(eta))
    error('Fixed joint static eta prior contains non-finite values.');
end
eta(1) = 0;
end


function eta = estimate_static_sensor_eta_local(bundle, S)
eta = zeros(1, numel(S.analysis_sensors));

if isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta')
    eta = get_fixed_joint_static_eta_local(S);
    return;
end

if ~isfield(S, 'static_eta_enable') || ~S.static_eta_enable || isempty(bundle) || ...
        ~isstruct(bundle) || ~isfield(bundle, 'sensor_id')
    return;
end

for is = 1:numel(S.analysis_sensors)
    sid = S.analysis_sensors(is);
    mask = bundle.sensor_id(:) == sid & ...
           isfinite(bundle.x(:)) & ...
           isfinite(bundle.v(:)) & ...
           isfinite(bundle.template_v(:)) & ...
           isfinite(bundle.template_dv_dx(:)) & ...
           isfinite(bundle.fit_weight(:));

    if nnz(mask) < max(20, S.min_fit_points / max(numel(S.analysis_sensors), 1))
        continue;
    end

    dv = bundle.v(mask) - bundle.template_v(mask);
    grad = bundle.template_dv_dx(mask);
    strong = abs(grad) >= max(0.10 * max(abs(grad), [], 'omitnan'), eps);
    dx = -dv(strong) ./ grad(strong);
    w = bundle.fit_weight(mask);
    w = w(strong);
    good = isfinite(dx) & isfinite(w);

    if nnz(good) < 20
        continue;
    end

    eta(is) = weighted_median_local(dx(good), w(good));
end

eta = eta - eta(1);
eta(1) = 0;
eta(~isfinite(eta)) = 0;
end


function eta = get_static_sensor_eta_local(bundle, S)
eta = zeros(1, numel(S.analysis_sensors));
if isfield(S, 'active_eta_model') && strcmpi(S.active_eta_model, 'none')
    return;
end
if isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta')
    eta = get_fixed_joint_static_eta_local(S);
    return;
end
if isstruct(bundle) && isfield(bundle, 'static_sensor_eta') && ...
        numel(bundle.static_sensor_eta) == numel(S.analysis_sensors)
    eta = bundle.static_sensor_eta(:).';
end
eta(1) = 0;
eta(~isfinite(eta)) = 0;
end


function m = weighted_median_local(x, w)
x = x(:);
w = w(:);
good = isfinite(x) & isfinite(w) & w > 0;
x = x(good);
w = w(good);
if isempty(x)
    m = NaN;
    return;
end
[x, order] = sort(x);
w = w(order);
cw = cumsum(w);
target = 0.5 * cw(end);
idx = find(cw >= target, 1, 'first');
m = x(idx);
end


function eo_candidates = build_eo_candidates_local(rot_freq_mean_hz, freq_search_hz, eo_pad)
if ~isfinite(rot_freq_mean_hz) || rot_freq_mean_hz <= 0
    error('Invalid rot_freq_mean_hz for EO candidate generation.');
end

freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min = max(1, ceil(freq_lo / rot_freq_mean_hz));
eo_max = max(eo_min, floor(freq_hi / rot_freq_mean_hz));

if eo_max < eo_min
    eo_center = max(1, round(mean([freq_lo, freq_hi]) / rot_freq_mean_hz));
    eo_min = max(1, eo_center - max(0, eo_pad));
    eo_max = max(eo_min, eo_center + max(0, eo_pad));
end

eo_candidates = eo_min:eo_max;
freq_candidates = eo_candidates .* rot_freq_mean_hz;
inside_band = freq_candidates >= freq_lo & freq_candidates <= freq_hi;
eo_candidates = eo_candidates(inside_band);

if isempty(eo_candidates)
    error('No EO candidates remain inside [%.3f, %.3f] Hz at rot_freq %.6f Hz.', ...
        freq_lo, freq_hi, rot_freq_mean_hz);
end
end


function sub = select_bundle_points_local(bundle, mask, S)
mask = mask(:) & bundle.finite_mask(:);

sub = bundle;

fields = { ...
    'x', ...
    'x_abs', ...
    'x_center', ...
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
    'edge_weight', ...
    'domain_weight', ...
    'query_weight', ...
    'gradient_fit_weight', ...
    'fit_weight', ...
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

if ~isfield(sub, 'fit_weight') || isempty(sub.fit_weight)
    sub.fit_weight = sub.template_weight(:);
end
sub.fit_weight(~isfinite(sub.fit_weight)) = S.weight_floor;
sub.fit_weight = max(sub.fit_weight, S.weight_floor);
end


function [vp_bundle, info] = build_vp_screening_bundle_local( ...
    source_bundle, stage_mask, fallback_bundle, S, stage_label)
info = make_empty_vp_screening_info_local(stage_label);
mode = lower(strtrim(getfield_default_local(S, ...
    'vp_bundle_mode', 'legacy_screening_query_safe')));

if isempty(fallback_bundle) || ~isstruct(fallback_bundle) || ...
        ~isfield(fallback_bundle, 'x') || isempty(fallback_bundle.x)
    vp_bundle = fallback_bundle;
    info.fallback_reason = 'empty_fallback_bundle';
    return;
end

if isempty(source_bundle) || ~isstruct(source_bundle) || ...
        ~isfield(source_bundle, 'x') || isempty(source_bundle.x)
    vp_bundle = fallback_bundle;
    info.fallback_reason = 'empty_source_bundle';
    info.selected_point_count = get_bundle_point_count_local(vp_bundle);
    info.valid_pass_count = get_bundle_valid_pass_count_local(vp_bundle);
    return;
end

n = numel(source_bundle.x);

if isempty(stage_mask) || numel(stage_mask) ~= n
    vp_bundle = fallback_bundle;
    info.fallback_reason = 'invalid_stage_mask';
    info.selected_point_count = get_bundle_point_count_local(vp_bundle);
    info.valid_pass_count = get_bundle_valid_pass_count_local(vp_bundle);
    return;
end

mask = stage_mask(:);
info.stage_mask_point_count = nnz(mask);

switch mode
    case {'full_fit_bundle', 'fit_bundle'}
        info.policy = 'fit_bundle';
        vp_bundle = fallback_bundle;
        info.selected_point_count = get_bundle_point_count_local(vp_bundle);
        info.valid_pass_count = get_bundle_valid_pass_count_local(vp_bundle);
        return;

    case 'legacy_screening_query_safe'
        info.policy = 'legacy_screening_query_safe';
        if isfield(source_bundle, 'finite_mask') && ...
                numel(source_bundle.finite_mask) == n
            mask = mask & source_bundle.finite_mask(:);
        end
        if isfield(source_bundle, 'inside_domain_mask') && ...
                numel(source_bundle.inside_domain_mask) == n
            mask = mask & source_bundle.inside_domain_mask(:);
        else
            info.missing_inside_domain_mask = true;
        end
        if isfield(source_bundle, 'dynamic_effective_mask') && ...
                numel(source_bundle.dynamic_effective_mask) == n
            mask = mask & source_bundle.dynamic_effective_mask(:);
        else
            info.missing_dynamic_effective_mask = true;
        end
        if isfield(source_bundle, 'main_pulse_mask') && ...
                numel(source_bundle.main_pulse_mask) == n
            mask = mask & source_bundle.main_pulse_mask(:);
        else
            info.missing_main_pulse_mask = true;
        end

    otherwise
        error('Unsupported Step05 VP bundle mode: %s.', mode);
end

candidate = select_bundle_points_local(source_bundle, mask, S);
info.selected_point_count = candidate.point_count;
info.screened_point_count = candidate.point_count;
info.valid_pass_count = candidate.valid_pass_count;

if candidate.point_count >= S.min_fit_points && ...
        candidate.valid_pass_count >= S.min_valid_pass_count
    vp_bundle = candidate;
    info.fallback_reason = '';
else
    vp_bundle = fallback_bundle;
    info.fallback_reason = sprintf( ...
        'fallback_fit_bundle_after_%s_points_%d_passes_%d', ...
        info.policy, candidate.point_count, candidate.valid_pass_count);
    info.selected_point_count = get_bundle_point_count_local(vp_bundle);
    info.valid_pass_count = get_bundle_valid_pass_count_local(vp_bundle);
end
end


function info = make_empty_vp_screening_info_local(stage_label)
info = struct( ...
    'stage', string(stage_label), ...
    'policy', 'legacy_screening_query_safe', ...
    'fallback_reason', '', ...
    'stage_mask_point_count', NaN, ...
    'screened_point_count', NaN, ...
    'selected_point_count', NaN, ...
    'valid_pass_count', NaN, ...
    'missing_inside_domain_mask', false, ...
    'missing_dynamic_effective_mask', false, ...
    'missing_main_pulse_mask', false);
end


function seed_table = solve_vp_seed_eo_scan_local(bundle, eo_candidates, S)
% VP seed scan.
%
% Formal 20251222 route:
%   x_static = T^{-1}_{s,b}(V) on the observed pulse branch
%   y_obs    = x - x_static - eta_s
%   y_obs ~= dx_c + a*sin(EO*theta) + b*cos(EO*theta)
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
    'vp_linear_rmse', inf, ...
    'static_residual_rmse_mm', inf, ...
    'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, ...
    'point_count', bundle.point_count, ...
    'rank_score', inf);

seed_table = repmat(template, numel(eo_candidates), 1);
static_eta = get_static_sensor_eta_local(bundle, S);

x = bundle.x(:);
v = bundle.v(:);
theta = bundle.theta(:);
T0 = bundle.template_v(:);
Tp = bundle.template_dv_dx(:);
w = bundle.fit_weight(:);
sensor_id = bundle.sensor_id(:);

observation_mode = lower(strtrim(getfield_default_local( ...
    S, 'vp_seed_observation_mode', 'static_residual_template_inversion')));

switch observation_mode
    case {'static_residual_template_inversion', 'template_inversion', 'static_residual'}
        x_static = invert_bundle_template_static_x_local(bundle, S);
        eta_vec0 = sensor_eta_vector_local(sensor_id, sensors, static_eta);
        y_obs = x - x_static - eta_vec0(:);
        valid = isfinite(x) & isfinite(v) & isfinite(theta) & ...
                isfinite(w) & isfinite(y_obs) & isfinite(eta_vec0(:));
    case {'gradient_displacement', 'gradient_displacement_vp', ...
            'q_space', 'qspace', 'displacement_gradient'}
        eta_vec0 = sensor_eta_vector_local(sensor_id, sensors, static_eta);
        grad_abs = abs(Tp);
        grad_ref = prctile(grad_abs(isfinite(grad_abs)), ...
            getfield_default_local(S, 'vp_gradient_reference_quantile', 95));
        if ~isfinite(grad_ref) || grad_ref <= 0
            grad_ref = max(grad_abs(isfinite(grad_abs)), [], 'omitnan');
        end
        if ~isfinite(grad_ref) || grad_ref <= 0
            return;
        end
        grad_min = max(0, getfield_default_local(S, ...
            'vp_gradient_min_ratio', 0.10)) * grad_ref;
        valid = isfinite(v) & ...
                isfinite(theta) & ...
                isfinite(T0) & ...
                isfinite(Tp) & ...
                isfinite(w) & ...
                isfinite(eta_vec0(:)) & ...
                grad_abs >= max(grad_min, eps);
        y_obs = nan(size(v));
        y_obs(valid) = -(v(valid) - T0(valid)) ./ Tp(valid) - eta_vec0(valid);
    case {'linear_voltage_residual', 'voltage_gradient'}
        valid = isfinite(v) & ...
                isfinite(theta) & ...
                isfinite(T0) & ...
                isfinite(Tp) & ...
                isfinite(w) & ...
                abs(Tp) > eps;
        y_obs = v - T0;
    otherwise
        error('Unsupported S.vp_seed_observation_mode: %s', observation_mode);
end

if nnz(valid) < S.min_fit_points
    return;
end

x = x(valid); %#ok<NASGU>
v = v(valid);
theta = theta(valid);
T0 = T0(valid);
Tp = Tp(valid);
w = w(valid);
sensor_id = sensor_id(valid);
y_obs = y_obs(valid);

wq = w;
wq(~isfinite(wq)) = S.weight_floor;
wq = max(wq, S.weight_floor);
if any(strcmp(observation_mode, {'gradient_displacement', ...
        'gradient_displacement_vp', 'q_space', 'qspace', ...
        'displacement_gradient'}))
    grad_abs = abs(Tp);
    grad_ref = prctile(grad_abs(isfinite(grad_abs)), ...
        getfield_default_local(S, 'vp_gradient_reference_quantile', 95));
    if isfinite(grad_ref) && grad_ref > 0
        grad_power = getfield_default_local(S, 'vp_gradient_weight_power', 2);
        grad_weight = min(grad_abs ./ max(grad_ref, eps), 1) .^ ...
            max(0, grad_power);
        wq = max(S.weight_floor, wq .* grad_weight);
    end
end

for iEO = 1:numel(eo_candidates)
    EO = eo_candidates(iEO);

    s = sin(EO .* theta);
    c = cos(EO .* theta);

    % Columns:
    %   common dx_c
    %   sin coefficient a
    %   cos coefficient b
    switch observation_mode
        case {'static_residual_template_inversion', 'template_inversion', ...
                'static_residual', 'gradient_displacement', ...
                'gradient_displacement_vp', 'q_space', 'qspace', ...
                'displacement_gradient'}
            X = [ones(size(s)), s, c];
        otherwise
            X = [-Tp, -Tp .* s, -Tp .* c];
    end

    if isfield(S, 'vp_include_sensor_eta') && S.vp_include_sensor_eta
        for is = 2:n_sensors
            switch observation_mode
                case {'static_residual_template_inversion', 'template_inversion', ...
                        'static_residual', 'gradient_displacement', ...
                        'gradient_displacement_vp', 'q_space', 'qspace', ...
                        'displacement_gradient'}
                    X = [X, double(sensor_id == sensors(is))]; %#ok<AGROW>
                otherwise
                    X = [X, -Tp .* double(sensor_id == sensors(is))]; %#ok<AGROW>
            end
        end
    end

    good = all(isfinite(X), 2) & isfinite(y_obs) & isfinite(wq);

    if nnz(good) < max(5, size(X,2) + 1)
        continue;
    end

    Xg = X(good, :);
    yg = y_obs(good);
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
    linear_resid = yg - Xg * beta;
    linear_rmse = sqrt(sum(wq(good) .* (linear_resid .^ 2)) ./ ...
        max(sum(wq(good)), eps));

    eta = static_eta;
    if isfield(S, 'vp_include_sensor_eta') && S.vp_include_sensor_eta && ...
            n_sensors > 1 && numel(beta) >= 3 + n_sensors - 1
        eta(2:end) = static_eta(2:end) + beta(4:end).';
    end
    eta(1) = 0;

    A = hypot(a, b);
    phi = atan2(b, a);

    if ~isfinite(A)
        A = 0.05;
    end
    if ~isfinite(dx_c)
        dx_c = 0;
    end
    eta(~isfinite(eta)) = 0;
    if isfinite(S.amplitude_limit_mm)
        A = min(max(abs(A), 0), abs(S.amplitude_limit_mm));
    end
    dx_c = clamp_scalar_local(dx_c, S.dx_c_limit_mm);
    eta = clamp_vector_local(eta, S.sensor_eta_limit_mm);

    [wrmse, prmse] = evaluate_full_template_model_local( ...
        bundle, EO, A, phi, dx_c, eta, S);

    seed_table(iEO).EO = EO;
    seed_table(iEO).A = A;
    seed_table(iEO).phi = wrap_to_pi_local(phi);
    seed_table(iEO).dx_c = dx_c;
    seed_table(iEO).sensor_eta = eta;
    seed_table(iEO).vp_linear_rmse = linear_rmse;
    seed_table(iEO).static_residual_rmse_mm = linear_rmse;
    seed_table(iEO).weighted_voltage_rmse = wrmse;
    seed_table(iEO).plain_voltage_rmse = prmse;
    seed_table(iEO).point_count = bundle.point_count;
    switch lower(strtrim(getfield_default_local( ...
            S, 'vp_rank_score_mode', 'linear_residual')))
        case {'linear', 'linear_residual', 'vp_linear', ...
                'static_residual_mm', 'static_residual', 'template_inversion'}
            seed_table(iEO).rank_score = linear_rmse;
        case {'full_seed', 'full_model_seed', 'weighted_voltage_rmse'}
            seed_table(iEO).rank_score = wrmse;
        otherwise
            error('Unsupported S.vp_rank_score_mode: %s', ...
                getfield_default_local(S, 'vp_rank_score_mode', ''));
    end
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


function [eo_selected, info] = select_eo_candidates_adaptive_local(seed_table, S)
info = struct();
info.policy = 'adaptive_full_fit';
info.top_k_min = S.vp_top_k_min;
info.top_k_max = S.vp_top_k_max;
info.top_rank_keep = getfield_default_local(S, 'vp_top_rank_keep', S.vp_top_k_min);
info.gap_ratio_keep = S.vp_gap_ratio_keep;
info.include_neighbor_eo = S.include_neighbor_eo;
info.neighbor_eo_radius = S.neighbor_eo_radius;
info.SelectionTable = table();
info.selected_EO = [];

eo_selected = [];

if isempty(seed_table)
    return;
end

T = struct2table(seed_table);
if ~all(ismember({'EO', 'rank_score'}, T.Properties.VariableNames))
    return;
end

valid = isfinite(T.EO) & isfinite(T.rank_score);
T = T(valid, :);

if isempty(T)
    return;
end

T = sortrows(T, {'rank_score', 'EO'}, {'ascend', 'ascend'});

top_min = max(1, floor(S.vp_top_k_min));
top_max = max(top_min, floor(S.vp_top_k_max));
top_rank_keep = max(top_min, floor(getfield_default_local(S, ...
    'vp_top_rank_keep', top_min)));
top_max = max(top_max, top_rank_keep);
gap_thr = max(S.vp_gap_ratio_keep, 0);

best_score = T.rank_score(1);
gap_ratio = (T.rank_score - best_score) ./ max(abs(best_score), eps);
keep = false(height(T), 1);

keep(1:min(top_min, height(T))) = true;
keep(1:min(top_rank_keep, height(T))) = true;
keep = keep | gap_ratio <= gap_thr;

eo = T.EO(keep).';

if isfield(S, 'include_neighbor_eo') && S.include_neighbor_eo
    radius = max(0, floor(S.neighbor_eo_radius));
    eo_expand = [];
    for i = 1:numel(eo)
        eo_expand = [eo_expand, (eo(i) - radius):(eo(i) + radius)]; %#ok<AGROW>
    end
    eo = unique(eo_expand);
end

eo_all = T.EO(:).';
eo = eo(ismember(eo, eo_all));

[~, loc] = ismember(eo, T.EO);
loc = sort(loc(loc > 0));
loc = loc(1:min(numel(loc), top_max));

eo_selected = T.EO(loc).';

T.vp_rank = (1:height(T)).';
T.vp_gap_ratio = gap_ratio;
T.selected_for_full_fit = ismember(T.EO, eo_selected);

info.SelectionTable = T(:, {'EO', 'rank_score', ...
    'weighted_voltage_rmse', 'vp_rank', 'vp_gap_ratio', ...
    'selected_for_full_fit'});
info.selected_EO = eo_selected;
end


function starts = build_direct_fit_start_points_local(seed, bundle, S)
eta0 = seed.sensor_eta;
if isempty(eta0) || numel(eta0) ~= numel(S.analysis_sensors)
    eta0 = get_static_sensor_eta_local(bundle, S);
end
eta0(1) = 0;
eta0(~isfinite(eta0)) = 0;

A0 = abs(seed.A);
if ~isfinite(A0) || A0 <= 0
    A0 = 0.05;
end

phi0 = seed.phi;
if ~isfinite(phi0)
    phi0 = 0;
end
phi0 = wrap_to_pi_local(phi0);

dx0 = seed.dx_c;
if ~isfinite(dx0)
    dx0 = 0;
end

eta_static = get_static_sensor_eta_local(bundle, S);
eta_zero = eta_static;
eta_zero(1) = 0;

starts = [
    pack_params_local(A0, phi0, dx0, eta0, S)
    pack_params_local(A0, phi0, dx0, eta_zero, S)
    pack_params_local(A0, phi0 + pi, dx0, eta0, S)
    pack_params_local(A0, phi0, 0, eta0, S)
    pack_params_local(0.05, 0, 0, eta_zero, S)
    ];

if isfield(S, 'direct_multistart_enable') && ~S.direct_multistart_enable
    starts = starts(1, :);
end

% Remove exact duplicate rows while preserving order.
[~, ia] = unique(round(starts, 12), 'rows', 'stable');
starts = starts(sort(ia), :);

% Keep the current bundle argument in the signature intentionally: future
% data-driven starts can use bundle statistics without changing callers.
if isempty(bundle) %#ok<UNRCH>
end
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
    'decision_score', inf, ...
    'decision_score_mode', string(getfield_default_local( ...
        S, 'eo_decision_score_mode', 'weighted_voltage_rmse')), ...
    'static_residual_rmse_mm', inf, ...
    'seed_weighted_voltage_rmse', inf, ...
    'seed_plain_voltage_rmse', inf, ...
    'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, ...
    'point_count', bundle.point_count, ...
    'exitflag', NaN, ...
    'start_count', NaN, ...
    'best_start_index', NaN);

Candidate = repmat(candidate_template, numel(final_eos), 1);

for i = 1:numel(final_eos)
    EO = final_eos(i);

    seed = get_seed_for_eo_local(seed_table, EO, S);
    seed_static_residual_rmse_mm = getfield_default_local( ...
        seed, 'static_residual_rmse_mm', inf);
    seed_weighted_voltage_rmse = getfield_default_local( ...
        seed, 'weighted_voltage_rmse', inf);
    seed_plain_voltage_rmse = getfield_default_local( ...
        seed, 'plain_voltage_rmse', inf);

    obj = @(p) direct_template_objective_local(p, EO, bundle, S);
    starts = build_direct_fit_start_points_local(seed, bundle, S);
    start_scores = inf(size(starts, 1), 1);
    for istart = 1:size(starts, 1)
        start_scores(istart) = obj(starts(istart, :));
    end
    [~, start_order] = sort(start_scores, 'ascend', 'MissingPlacement', 'last');
    if isfield(S, 'direct_multistart_refine_top_n') && ...
            isfinite(S.direct_multistart_refine_top_n) && ...
            S.direct_multistart_refine_top_n > 0
        start_order = start_order(1:min(numel(start_order), S.direct_multistart_refine_top_n));
    end

    opts = optimset( ...
        'Display', 'off', ...
        'MaxIter', S.fminsearch_max_iter, ...
        'MaxFunEvals', S.fminsearch_max_fun, ...
        'TolX', 1e-7, ...
        'TolFun', 1e-9);

    p_best = starts(1, :);
    best_obj = inf;
    exitflag = -1;
    best_start_index = 1;

    for iorder = 1:numel(start_order)
        istart = start_order(iorder);
        p0 = starts(istart, :);
        try
            [p_try, obj_try, exitflag_try] = fminsearch(obj, p0, opts);
        catch
            p_try = p0;
            obj_try = obj(p0);
            exitflag_try = -1;
        end

        if isfinite(obj_try) && obj_try < best_obj
            p_best = p_try;
            best_obj = obj_try;
            exitflag = exitflag_try;
            best_start_index = istart;
        end
    end

    [A, phi, dx_c, eta] = unpack_params_for_bundle_local(p_best, S, bundle);

    [wrmse, prmse, ~, ~, objective_score] = evaluate_full_template_model_local( ...
        bundle, EO, A, phi, dx_c, eta, S);

    Candidate(i).EO = EO;
    Candidate(i).A = A;
    Candidate(i).phi = wrap_to_pi_local(phi);
    Candidate(i).dx_c = dx_c;
    Candidate(i).sensor_eta = eta;
    Candidate(i).objective_score = objective_score;
    Candidate(i).decision_score_mode = string(getfield_default_local( ...
        S, 'eo_decision_score_mode', 'weighted_voltage_rmse'));
    Candidate(i).static_residual_rmse_mm = seed_static_residual_rmse_mm;
    Candidate(i).seed_weighted_voltage_rmse = seed_weighted_voltage_rmse;
    Candidate(i).seed_plain_voltage_rmse = seed_plain_voltage_rmse;
    Candidate(i).weighted_voltage_rmse = wrmse;
    Candidate(i).plain_voltage_rmse = prmse;
    Candidate(i).point_count = bundle.point_count;
    Candidate(i).exitflag = exitflag;
    Candidate(i).start_count = size(starts, 1);
    Candidate(i).best_start_index = best_start_index;
    Candidate(i).decision_score = compute_candidate_decision_score_local( ...
        Candidate(i), S);
end

best_idx = select_best_candidate_index_local(Candidate, S);

best = Candidate(best_idx);

if ~isfinite(best.weighted_voltage_rmse)
    result = make_failed_result_local('all_candidates_failed');
    result.CandidateTable = make_sorted_candidate_table_local(Candidate);
    result.EOAmbiguity = summarize_eo_candidate_ambiguity_local(Candidate);
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
result.decision_score = best.decision_score;
result.decision_score_mode = best.decision_score_mode;
result.static_residual_rmse_mm = best.static_residual_rmse_mm;
result.seed_weighted_voltage_rmse = best.seed_weighted_voltage_rmse;
result.seed_plain_voltage_rmse = best.seed_plain_voltage_rmse;
result.weighted_voltage_rmse = wrmse;
result.plain_voltage_rmse = prmse;
result.point_count = bundle.point_count;
result.valid_segment_count = numel(unique(bundle.pass_id));
result.Coverage = coverage;
result.CandidateTable = make_sorted_candidate_table_local(Candidate);
result.EOAmbiguity = summarize_eo_candidate_ambiguity_local(Candidate);
result.V_pred_preview = Vpred(1:min(end, 200));

if S.reject_ambiguous_frequency && ...
        isfield(result.EOAmbiguity, 'rmse_gap_ratio') && ...
        isfinite(result.EOAmbiguity.rmse_gap_ratio) && ...
        result.EOAmbiguity.rmse_gap_ratio < max(S.eo_low_gap_ratio_threshold, S.unconstrained_unidentifiable_gap_ratio)
    result.status = 'unidentifiable_without_gauge';
    result.failure_reason = sprintf( ...
        ['EO candidate gap ratio %.4g is too small for an unconstrained direct-template ' ...
        'frequency decision; top candidates=%s'], ...
        result.EOAmbiguity.rmse_gap_ratio, result.EOAmbiguity.top_eo_list);
    result.EO_id = NaN;
    result.fn_id = NaN;
end
end


function T = make_sorted_candidate_table_local(Candidate)
if isempty(Candidate)
    T = table();
    return;
end
T = struct2table(Candidate);
decision_mode = "";
if any(strcmp(T.Properties.VariableNames, 'decision_score_mode'))
    mode_values = string(T.decision_score_mode);
    mode_values = mode_values(strlength(mode_values) > 0);
    if ~isempty(mode_values)
        decision_mode = lower(strtrim(mode_values(1)));
    end
end

if any(strcmp(T.Properties.VariableNames, 'decision_score')) && ...
        ~any(decision_mode == ["weighted_voltage_rmse", ...
        "voltage_rmse", "full_fit_rmse"])
    T = sortrows(T, {'decision_score', 'weighted_voltage_rmse', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
elseif any(strcmp(T.Properties.VariableNames, 'static_residual_rmse_mm'))
    T = sortrows(T, {'weighted_voltage_rmse', 'static_residual_rmse_mm', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
elseif any(strcmp(T.Properties.VariableNames, 'objective_score'))
    T = sortrows(T, {'objective_score', 'weighted_voltage_rmse', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
elseif any(strcmp(T.Properties.VariableNames, 'weighted_voltage_rmse'))
    T = sortrows(T, {'weighted_voltage_rmse', 'EO'}, {'ascend', 'ascend'});
end
end


function best_idx = select_best_candidate_index_local(Candidate, S)
scores = arrayfun(@(c) compute_candidate_decision_score_local(c, S), Candidate);
rmse = [Candidate.weighted_voltage_rmse];
eos = [Candidate.EO];
scores = scores(:).';
rmse = rmse(:).';
eos = eos(:).';

valid = isfinite(scores) & isfinite(rmse) & isfinite(eos);
if ~any(valid)
    scores = [Candidate.objective_score];
    valid = isfinite(scores) & isfinite(rmse) & isfinite(eos);
end

if ~any(valid)
    best_idx = 1;
    return;
end

idx = find(valid);
idx = idx(:);
score_col = scores(idx);
rmse_col = rmse(idx);
eo_col = eos(idx);
sort_matrix = [score_col(:), rmse_col(:), eo_col(:)];
[~, order] = sortrows(sort_matrix, [1 2 3]);
best_idx = idx(order(1));
end


function score = compute_candidate_decision_score_local(candidate, S)
mode = lower(strtrim(char(string(getfield_default_local( ...
    S, 'eo_decision_score_mode', 'weighted_voltage_rmse')))));

switch mode
    case {'static_residual_mm', 'static_residual', ...
            'template_inversion', 'vp_rank_score'}
        score = getfield_default_local(candidate, ...
            'static_residual_rmse_mm', inf);
    case {'weighted_voltage_rmse', 'voltage_rmse', 'full_fit_rmse'}
        score = getfield_default_local(candidate, ...
            'weighted_voltage_rmse', inf);
    case {'objective', 'objective_score'}
        score = getfield_default_local(candidate, ...
            'objective_score', inf);
    otherwise
        error('Unsupported S.eo_decision_score_mode: %s', mode);
end

if ~isfinite(score)
    score = inf;
end
end


function ambiguity = summarize_eo_candidate_ambiguity_local(Candidate)
ambiguity = struct( ...
    'candidate_count', numel(Candidate), ...
    'finite_candidate_count', 0, ...
    'best_eo', NaN, ...
    'best_rmse_v', NaN, ...
    'best_decision_score', NaN, ...
    'best_static_residual_rmse_mm', NaN, ...
    'second_best_eo', NaN, ...
    'second_best_rmse_v', NaN, ...
    'second_best_decision_score', NaN, ...
    'second_best_static_residual_rmse_mm', NaN, ...
    'rmse_gap_v', NaN, ...
    'rmse_gap_ratio', NaN, ...
    'decision_score_gap', NaN, ...
    'decision_score_gap_ratio', NaN, ...
    'top_eo_list', "", ...
    'top_rmse_list', "", ...
    'top_decision_score_list', "");

if isempty(Candidate)
    return;
end

eo = [Candidate.EO];
rmse = [Candidate.weighted_voltage_rmse];
eo = eo(:).';
rmse = rmse(:).';
if isfield(Candidate, 'decision_score')
    decision_score = [Candidate.decision_score];
else
    decision_score = rmse;
end
if isfield(Candidate, 'static_residual_rmse_mm')
    static_residual = [Candidate.static_residual_rmse_mm];
else
    static_residual = nan(size(rmse));
end
decision_score = decision_score(:).';
static_residual = static_residual(:).';

valid = isfinite(decision_score) & isfinite(rmse) & isfinite(eo);
ambiguity.finite_candidate_count = nnz(valid);

if ~any(valid)
    return;
end

eo_valid = eo(valid);
rmse_valid = rmse(valid);
decision_valid = decision_score(valid);
static_valid = static_residual(valid);
[~, order] = sortrows([decision_valid(:), rmse_valid(:), eo_valid(:)], [1 2 3]);
eo_sorted = eo_valid(order);
rmse_sorted = rmse_valid(order);
decision_sorted = decision_valid(order);
static_sorted = static_valid(order);

ambiguity.best_eo = eo_sorted(1);
ambiguity.best_rmse_v = rmse_sorted(1);
ambiguity.best_decision_score = decision_sorted(1);
ambiguity.best_static_residual_rmse_mm = static_sorted(1);

if numel(eo_sorted) >= 2
    ambiguity.second_best_eo = eo_sorted(2);
    ambiguity.second_best_rmse_v = rmse_sorted(2);
    ambiguity.second_best_decision_score = decision_sorted(2);
    ambiguity.second_best_static_residual_rmse_mm = static_sorted(2);
    ambiguity.rmse_gap_v = rmse_sorted(2) - rmse_sorted(1);
    ambiguity.rmse_gap_ratio = ambiguity.rmse_gap_v / max(rmse_sorted(1), eps);
    ambiguity.decision_score_gap = decision_sorted(2) - decision_sorted(1);
    ambiguity.decision_score_gap_ratio = ambiguity.decision_score_gap / ...
        max(abs(decision_sorted(1)), eps);
end

n_top = min(5, numel(eo_sorted));
ambiguity.top_eo_list = strjoin(compose('%d', eo_sorted(1:n_top)), '|');
ambiguity.top_rmse_list = strjoin(compose('%.6g', rmse_sorted(1:n_top)), '|');
ambiguity.top_decision_score_list = strjoin( ...
    compose('%.6g', decision_sorted(1:n_top)), '|');
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

seed.A = abs(seed.A);

if ~isfinite(seed.A)
    seed.A = 0.05;
end

seed.phi = wrap_to_pi_local(seed.phi);

if ~isfinite(seed.phi)
    seed.phi = 0;
end

if ~isfinite(seed.dx_c)
    seed.dx_c = 0;
end

if isempty(seed.sensor_eta) || numel(seed.sensor_eta) ~= numel(S.analysis_sensors)
    seed.sensor_eta = zeros(1, numel(S.analysis_sensors));
end

seed.sensor_eta(1) = 0;
seed.sensor_eta(~isfinite(seed.sensor_eta)) = 0;
end


function p = pack_params_local(A, phi, dx_c, eta, S)
eta = eta(:).';

if numel(eta) ~= numel(S.analysis_sensors)
    eta = zeros(1, numel(S.analysis_sensors));
end

eta(1) = 0;

% Parameter vector:
%   [A, phi, dx_c, eta_2, eta_3, ...]
if isfield(S, 'fit_sensor_eta') && S.fit_sensor_eta
    p = [A, phi, dx_c, eta(2:end)];
else
    p = [A, phi, dx_c];
end
p(~isfinite(p)) = 0;
end


function [A, phi, dx_c, eta] = unpack_params_local(p, S)
p = p(:).';

A = p(1);
phi = p(2);
dx_c = p(3);

eta = zeros(1, numel(S.analysis_sensors));

if isfield(S, 'fit_sensor_eta') && S.fit_sensor_eta && numel(p) > 3
    eta(2:end) = p(4:end);
end

% Raw unpacking only normalizes amplitude sign and phase.  The official
% bounded route applies physical clamps in unpack_params_for_bundle_local.
A = abs(A);
phi = wrap_to_pi_local(phi);

eta(1) = 0;
eta(~isfinite(eta)) = 0;
end


function [A, phi, dx_c, eta] = unpack_params_for_bundle_local(p, S, bundle)
[A, phi, dx_c, eta] = unpack_params_local(p, S);
if isfield(S, 'active_eta_model') && strcmpi(S.active_eta_model, 'none')
    eta = zeros(1, numel(S.analysis_sensors));
elseif ~isfield(S, 'fit_sensor_eta') || ~S.fit_sensor_eta
    eta = get_static_sensor_eta_local(bundle, S);
end

enforce_bounds = ...
    (isfield(S, 'enforce_parameter_bounds') && S.enforce_parameter_bounds) || ...
    (isfield(S, 'enforce_legacy_parameter_bounds') && ...
     S.enforce_legacy_parameter_bounds);

if enforce_bounds
    if isfinite(S.amplitude_limit_mm)
        A = min(abs(A), abs(S.amplitude_limit_mm));
    end
    dx_c = clamp_scalar_local(dx_c, S.dx_c_limit_mm);
    fixed_eta_model = isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta') && ...
        (~isfield(S, 'fit_sensor_eta') || ~S.fit_sensor_eta);
    if ~fixed_eta_model
        eta = clamp_vector_local(eta, ...
            resolve_sensor_eta_limit_step05_local(bundle, S));
    end
    if ~isempty(eta)
        eta(1) = 0;
    end
end
end


function value = direct_template_objective_local(p, EO, bundle, S)
[A, phi, dx_c, eta] = unpack_params_for_bundle_local(p, S, bundle);

[~, ~, ~, ~, objective_score] = evaluate_full_template_model_local( ...
    bundle, EO, A, phi, dx_c, eta, S);

if ~isfinite(objective_score)
    objective_score = 1e12;
end

value = objective_score;
end


function [wrmse, prmse, Vpred, coverage, objective_score] = evaluate_full_template_model_local( ...
    bundle, EO, A, phi, dx_c, eta, S)

sensors = S.analysis_sensors(:).';
extrapolation_mode = 'nan_outside';
if isfield(S, 'extrapolation_mode') && ~isempty(S.extrapolation_mode)
    extrapolation_mode = lower(strtrim(S.extrapolation_mode));
end
interp_method = getfield_default_local(S, 'template_interp_method', 'pchip');

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

        switch extrapolation_mode
            case 'clip_to_grid'
                xq_eval = min(max(xq(idx), min(tpl.x_grid)), max(tpl.x_grid));
            otherwise
                xq_eval = xq(idx);
        end

        Vpred(idx) = interp1( ...
            tpl.x_grid(:), ...
            tpl.v_grid(:), ...
            xq_eval, ...
            interp_method, ...
            NaN);
    end
end

valid = isfinite(Vpred) & ...
        isfinite(bundle.v(:)) & ...
        isfinite(bundle.fit_weight(:));

residual_sse = inf;
overshoot_penalty = inf;
eta_reg_penalty = inf;
invalid_penalty = inf;
objective_score = inf;

if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
else
    res = bundle.v(:) - Vpred(:);
    w = bundle.fit_weight(:);
    w(~isfinite(w)) = S.weight_floor;
    w = max(w, S.weight_floor);
    residual_sse = nansum(w(valid) .* res(valid).^2);
    invalid_fraction = 1 - nnz(valid) / max(bundle.point_count, 1);
    residual_scale = residual_sse / max(nnz(valid), 1);
    invalid_penalty = S.invalid_query_penalty_scale * ...
        bundle.point_count * invalid_fraction * max(residual_scale, eps);
    overshoot_penalty = S.overshoot_penalty_weight * ...
        nansum(w(valid) .* overshoot(valid).^2);
    eta_reg_penalty = bundle.point_count * ...
        (S.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta(:).^2)))^2;
    rmse_objective = residual_sse + overshoot_penalty + ...
        eta_reg_penalty + invalid_penalty;

    if isfield(S, 'objective_mode') && strcmpi(S.objective_mode, 'corr')
        corr_score = compute_weighted_sensor_correlation_score_local( ...
            bundle, Vpred, valid, w);
        objective_score = bundle.point_count * max(0, 1 - corr_score) + ...
            overshoot_penalty + eta_reg_penalty + invalid_penalty;
    else
        corr_score = NaN;
        objective_score = rmse_objective;
    end

    wrmse = sqrt(objective_score / max(bundle.point_count, 1));

    prmse = sqrt(nanmean(res(valid).^2));
end

coverage = struct();
coverage.clamp_fraction = mean(overshoot > 0, 'omitnan');
coverage.overshoot_rms_mm = sqrt(nanmean(overshoot(:).^2));
coverage.max_overshoot_mm = max(overshoot(:), [], 'omitnan');
coverage.residual_sse = residual_sse;
coverage.rmse_objective = rmse_objective;
coverage.correlation_score = corr_score;
coverage.overshoot_penalty = overshoot_penalty;
coverage.sensor_eta_reg_penalty = eta_reg_penalty;
coverage.invalid_penalty = invalid_penalty;
coverage.valid_query_fraction = nnz(valid) / max(bundle.point_count, 1);
coverage.extrapolation_mode = string(extrapolation_mode);
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
    wi = w(idx);
    good = isfinite(y) & isfinite(yp) & isfinite(wi) & wi > 0;
    y = y(good);
    yp = yp(good);
    wi = wi(good);
    if numel(y) < 10
        continue;
    end

    wi = wi ./ max(sum(wi), eps);
    my = sum(wi .* y);
    mp = sum(wi .* yp);
    yc = y - my;
    pc = yp - mp;
    denom = sqrt(sum(wi .* yc.^2) * sum(wi .* pc.^2));
    if denom <= eps
        continue;
    end
    scores(end + 1, 1) = sum(wi .* yc .* pc) / denom; %#ok<AGROW>
    weights(end + 1, 1) = numel(y); %#ok<AGROW>
end

if isempty(scores)
    score = -inf;
else
    weights = weights ./ sum(weights);
    score = sum(weights .* scores);
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


function eta_limit = resolve_sensor_eta_limit_step05_local(bundle, S)
eta_limit = abs(S.sensor_eta_limit_mm) * ones(1, numel(S.analysis_sensors));

if ~isfield(S, 'use_template_eta_limit') || ~S.use_template_eta_limit || ...
        ~isfield(bundle, 'TemplateMini') || isempty(bundle.TemplateMini)
    return;
end

for i = 1:min(numel(bundle.TemplateMini), numel(eta_limit))
    if isfield(bundle.TemplateMini(i), 'eta_limit_mm') && ...
            isfinite(bundle.TemplateMini(i).eta_limit_mm) && ...
            bundle.TemplateMini(i).eta_limit_mm > 0
        eta_limit(i) = min(eta_limit(i), abs(bundle.TemplateMini(i).eta_limit_mm));
    end
end
end


function [phase_ref, info] = resolve_phase_safe_reference_local( ...
    prev_result, core_result, seed_table, window_id, S)
mode = lower(strtrim(getfield_default_local(S, ...
    'phase_safe_reference_mode', 'prev_window')));

if strcmpi(mode, 'prev_window')
    [ok_prev, reason_prev] = is_prev_phase_reference_usable_local( ...
        prev_result, window_id, S);

    if ok_prev
        phase_ref = phase_reference_from_result_local(prev_result, S, ...
            'prev_window_final');
        info = struct('source', 'prev_window_final', ...
            'reason', 'prev_window_quality_pass', ...
            'window_id', prev_result.window_id);
        return;
    end

    fallback_mode = lower(strtrim(getfield_default_local(S, ...
        'phase_safe_fallback_mode', 'linear_vp')));
    [phase_ref, fallback_source] = phase_reference_from_fallback_local( ...
        fallback_mode, core_result, seed_table, S);
    info = struct('source', fallback_source, ...
        'reason', sprintf('%s_then_%s', reason_prev, fallback_mode), ...
        'window_id', window_id);
    return;
end

if strcmpi(mode, 'linear_vp')
    phase_ref = phase_reference_from_seed_local(seed_table, S);
    info = struct('source', 'linear_vp_seed', ...
        'reason', 'configured_linear_vp', ...
        'window_id', window_id);
else
    phase_ref = phase_reference_from_result_local(core_result, S, ...
        'current_core_final');
    info = struct('source', 'current_core_final', ...
        'reason', 'configured_core_current', ...
        'window_id', window_id);
end
end


function [phase_ref, source] = phase_reference_from_fallback_local( ...
    fallback_mode, core_result, seed_table, S)
switch lower(strtrim(fallback_mode))
    case 'linear_vp'
        phase_ref = phase_reference_from_seed_local(seed_table, S);
        source = 'linear_vp_seed';
    case 'final_core'
        phase_ref = phase_reference_from_result_local(core_result, S, ...
            'current_core_final');
        source = 'current_core_final';
    otherwise
        phase_ref = phase_reference_from_result_local(core_result, S, ...
            'current_core_final');
        source = 'current_core_final';
end
end


function phase_ref = phase_reference_from_seed_local(seed_table, S)
if isempty(seed_table)
    phase_ref = phase_reference_from_result_local(make_failed_result_local( ...
        'empty_seed_table'), S, 'linear_vp_seed_empty');
    return;
end

seed = seed_table(1);
phase_ref = struct();
phase_ref.window_id = NaN;
phase_ref.EO_id = seed.EO;
phase_ref.A_id = min(abs(seed.A), abs(S.amplitude_limit_mm));
phase_ref.phi_id_wrapped = wrap_to_pi_local(seed.phi);
phase_ref.dx_c_id = clamp_scalar_local(seed.dx_c, S.dx_c_limit_mm);
phase_ref.sensor_eta_id = zeros(1, numel(S.analysis_sensors));

if isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta') && ...
        isfield(S, 'sensor_eta_prior_mm') && ...
        numel(S.sensor_eta_prior_mm) == numel(S.analysis_sensors)
    phase_ref.sensor_eta_id = double(S.sensor_eta_prior_mm(:).');
elseif isfield(seed, 'sensor_eta') && ~isempty(seed.sensor_eta)
    n = min(numel(phase_ref.sensor_eta_id), numel(seed.sensor_eta));
    phase_ref.sensor_eta_id(1:n) = seed.sensor_eta(1:n);
    phase_ref.sensor_eta_id = clamp_vector_local( ...
        phase_ref.sensor_eta_id, S.sensor_eta_limit_mm);
end
phase_ref.sensor_eta_id(1) = 0;
phase_ref.weighted_voltage_rmse = seed.weighted_voltage_rmse;
phase_ref.reference_source = 'linear_vp_seed';
end


function phase_ref = phase_reference_from_result_local(result, S, source)
phase_ref = struct();
phase_ref.window_id = getfield_default_local(result, 'window_id', NaN);
phase_ref.EO_id = getfield_default_local(result, 'EO_id', NaN);
phase_ref.A_id = min(abs(getfield_default_local(result, 'A_id', NaN)), ...
    abs(S.amplitude_limit_mm));
phase_ref.phi_id_wrapped = wrap_to_pi_local( ...
    getfield_default_local(result, 'phi_id_wrapped', NaN));
phase_ref.dx_c_id = clamp_scalar_local( ...
    getfield_default_local(result, 'dx_c_id', NaN), S.dx_c_limit_mm);
phase_ref.sensor_eta_id = zeros(1, numel(S.analysis_sensors));

eta = getfield_default_local(result, 'sensor_eta_id', []);
if isfield(S, 'active_eta_model') && ...
        strcmpi(S.active_eta_model, 'fixed_joint_static_eta') && ...
        isfield(S, 'sensor_eta_prior_mm') && ...
        numel(S.sensor_eta_prior_mm) == numel(S.analysis_sensors)
    phase_ref.sensor_eta_id = double(S.sensor_eta_prior_mm(:).');
elseif ~isempty(eta)
    n = min(numel(phase_ref.sensor_eta_id), numel(eta));
    phase_ref.sensor_eta_id(1:n) = eta(1:n);
    phase_ref.sensor_eta_id = clamp_vector_local( ...
        phase_ref.sensor_eta_id, S.sensor_eta_limit_mm);
end
phase_ref.sensor_eta_id(1) = 0;
phase_ref.weighted_voltage_rmse = getfield_default_local( ...
    result, 'weighted_voltage_rmse', NaN);
phase_ref.reference_source = source;
end


function [ok, reason] = is_prev_phase_reference_usable_local( ...
    prev_result, window_id, S)
ok = false;

if isempty(prev_result) || ~isstruct(prev_result)
    reason = 'no_previous_window';
    return;
end

refresh_every = getfield_default_local(S, 'phase_safe_refresh_every', 0);
if refresh_every > 0 && mod(window_id - 1, refresh_every) == 0
    reason = 'scheduled_refresh';
    return;
end

required = {'window_id', 'EO_id', 'A_id', 'phi_id_wrapped', ...
    'dx_c_id', 'weighted_voltage_rmse', 'Coverage'};
for i = 1:numel(required)
    if ~isfield(prev_result, required{i})
        reason = ['previous_missing_', required{i}];
        return;
    end
end

vals = [prev_result.window_id, prev_result.EO_id, prev_result.A_id, ...
    prev_result.phi_id_wrapped, prev_result.dx_c_id, ...
    prev_result.weighted_voltage_rmse];
if any(~isfinite(vals))
    reason = 'previous_nonfinite_reference';
    return;
end

max_clamp = getfield_default_local(S, ...
    'phase_safe_prev_max_clamp_fraction', inf);
if isfield(prev_result.Coverage, 'clamp_fraction') && ...
        prev_result.Coverage.clamp_fraction > max_clamp
    reason = 'previous_clamped';
    return;
end

min_gap = getfield_default_local(S, ...
    'phase_safe_prev_min_final_gap_ratio', -inf);
gap_ratio = final_candidate_gap_ratio_step05_local(prev_result);
if isfinite(min_gap) && gap_ratio < min_gap
    reason = 'previous_final_gap_too_small';
    return;
end

ok = true;
reason = 'prev_window_quality_pass';
end


function gap_ratio = final_candidate_gap_ratio_step05_local(result)
gap_ratio = inf;
if ~isfield(result, 'CandidateTable') || isempty(result.CandidateTable) || ...
        height(result.CandidateTable) < 2 || ...
        ~ismember('weighted_voltage_rmse', result.CandidateTable.Properties.VariableNames)
    return;
end

rmse = result.CandidateTable.weighted_voltage_rmse;
rmse = sort(rmse(isfinite(rmse)), 'ascend');
if numel(rmse) < 2
    return;
end

gap_ratio = (rmse(2) - rmse(1)) / max(rmse(1), eps);
end


function y = clamp_scalar_local(x, limit_abs)
if ~isfinite(x)
    y = x;
    return;
end

if ~isfinite(limit_abs)
    y = x;
else
    limit_abs = abs(limit_abs);
    y = min(max(x, -limit_abs), limit_abs);
end
end


function y = clamp_vector_local(x, limit_abs)
y = x;
if isempty(limit_abs) || all(~isfinite(limit_abs))
    return;
end

limit_abs = abs(limit_abs(:).');
if numel(limit_abs) == 1
    limit_abs = repmat(limit_abs, size(y));
elseif numel(limit_abs) < numel(y)
    limit_abs(numel(limit_abs)+1:numel(y)) = limit_abs(end);
else
    limit_abs = limit_abs(1:numel(y));
end
limit_abs(~isfinite(limit_abs)) = inf;
y = min(max(y, -limit_abs), limit_abs);
end


function expanded_mask = build_phase_safe_expansion_mask_local(bundle, result, S)
% Second-stage deterministic phase-safe expansion.
%
% Use fitted vibration/offset parameters to compute:
%   x_query = x - dx_c - A*sin(EO*theta + phi)
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
if isfield(bundle, 'main_pulse_mask') && ...
        numel(bundle.main_pulse_mask) == numel(expanded_mask)
    expanded_mask = expanded_mask & bundle.main_pulse_mask(:);
end
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

    core_score = getfield_default_local(core_result, ...
        'decision_score', core_result.weighted_voltage_rmse);
    expanded_score = getfield_default_local(expanded_result, ...
        'decision_score', expanded_result.weighted_voltage_rmse);
    improve = core_score - expanded_score;

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
    'xc', NaN, ...
    'x_domain', [NaN NaN], ...
    'eta_limit_mm', NaN), ...
    numel(sensors), 1);

for i = 1:numel(sensors)
    tpl = get_template_entry_local(Template, sensors(i), blade_id);

    mini(i).sensor_id = sensors(i);
    mini(i).blade_id = blade_id;
    mini(i).x_grid = tpl.x_grid(:);
    mini(i).v_grid = tpl.v_grid(:);
    mini(i).dv_dx = tpl.dv_dx(:);
    mini(i).xc = get_template_xcenter_local(tpl);
    mini(i).x_domain = tpl.x_domain;
    if isfield(tpl, 'eta_limit_mm') && isfinite(tpl.eta_limit_mm)
        mini(i).eta_limit_mm = tpl.eta_limit_mm;
    end
end

bundle.TemplateMini = mini;
end


function tpl = get_template_entry_local(Template, sid, blade_id)
tpl = [];

if isfield(Template, 'SensorBlade')
    for i = 1:numel(Template.SensorBlade)
        if Template.SensorBlade(i).sensor_id == sid && ...
                Template.SensorBlade(i).blade_id == blade_id

            tpl = Template.SensorBlade(i);
            return;
        end
    end
end

% Compatibility path for the verified 20251222 legacy assembled template
% used only in override diagnostics. That template stores one target blade
% as Template.Sensor instead of the foundation SensorBlade array.
if isfield(Template, 'Sensor') && ...
        (~isfield(Template, 'TargetBlade') || Template.TargetBlade == blade_id)
    for i = 1:numel(Template.Sensor)
        if Template.Sensor(i).sensor_id == sid
            tpl = Template.Sensor(i);
            return;
        end
    end
end
end


function theta_std_deg = get_standard_angle_local(Template, sid, blade_id)
angles = Template.Standard_Relative_Angles;

if ndims(angles) ~= 2 || size(angles, 2) < blade_id
    error('Step05:TemplateContract', ...
        'Template.Standard_Relative_Angles does not cover Blade%d.', ...
        blade_id);
end

row = [];
if isfield(Template, 'Sensor_IDs') && ~isempty(Template.Sensor_IDs)
    sensor_ids = double(Template.Sensor_IDs(:));
    if size(angles, 1) >= max([sensor_ids; sid])
        row = sid;
    else
        row = find(sensor_ids == sid, 1, 'first');
    end
elseif size(angles, 1) >= sid
    row = sid;
end

if isempty(row) || row < 1 || row > size(angles, 1)
    error('Step05:TemplateContract', ...
        ['Template.Standard_Relative_Angles cannot map CH%d to a row; ' ...
         'provide direct sensor-indexed rows or Template.Sensor_IDs.'], sid);
end

theta_std_deg = angles(row, blade_id);
if ~isscalar(theta_std_deg) || ~isfinite(theta_std_deg)
    error('Step05:TemplateContract', ...
        'Template.Standard_Relative_Angles has invalid value for CH%d Blade%d.', ...
        sid, blade_id);
end
end


function xc = get_template_xcenter_local(tpl)
xc = 0;
if isstruct(tpl) && isfield(tpl, 'xc') && isfinite(tpl.xc)
    xc = tpl.xc;
elseif isstruct(tpl) && isfield(tpl, 'xc_mm') && isfinite(tpl.xc_mm)
    xc = tpl.xc_mm;
elseif isstruct(tpl) && isfield(tpl, 'x_center_mm') && isfinite(tpl.x_center_mm)
    xc = tpl.x_center_mm;
end
end


function row = make_empty_trend_row_local()
row = struct( ...
    'case_name', "", ...
    'blade_id', NaN, ...
    'sensor_tag', "", ...
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
    'decision_score', NaN, ...
    'static_residual_rmse_mm', NaN, ...
    'seed_weighted_voltage_rmse', NaN, ...
    'weighted_voltage_rmse', NaN, ...
    'plain_voltage_rmse', NaN, ...
    'second_best_EO_id', NaN, ...
    'second_best_rmse_v', NaN, ...
    'second_best_decision_score', NaN, ...
    'eo_decision_score_gap', NaN, ...
    'eo_decision_score_gap_ratio', NaN, ...
    'eo_rmse_gap_v', NaN, ...
    'eo_rmse_gap_ratio', NaN, ...
    'eo_candidate_count', NaN, ...
    'eo_finite_candidate_count', NaN, ...
    'eo_top_list', "", ...
    'core_EO_id', NaN, ...
    'expanded_EO_id', NaN, ...
    'phase_reference_EO_id', NaN, ...
    'final_EO_changed_after_expansion', NaN, ...
    'valid_segment_count', NaN, ...
    'valid_sensor_count', NaN, ...
    'valid_lap_count', NaN, ...
    'expected_valid_pass_count', NaN, ...
    'min_valid_pass_count', NaN, ...
    'min_valid_sensor_count', NaN, ...
    'min_valid_lap_count', NaN, ...
    'valid_pass_fraction', NaN, ...
    'point_count', NaN, ...
    'vp_point_count', NaN, ...
    'core_vp_point_count', NaN, ...
    'expanded_vp_point_count', NaN, ...
    'vp_screening_policy', "", ...
    'vp_screening_fallback', "", ...
    'candidate_point_count', NaN, ...
    'core_point_count', NaN, ...
    'expanded_point_count', NaN, ...
    'inside_domain_count', NaN, ...
    'inside_guard_count', NaN, ...
    'dynamic_effective_count', NaN, ...
    'main_pulse_count', NaN, ...
    'query_guard_source', "", ...
    'query_guard_mean_mm', NaN, ...
    'query_safe_left_min_mm', NaN, ...
    'query_safe_right_max_mm', NaN, ...
    'local_A_id', NaN, ...
    'local_EO_id', NaN, ...
    'local_fn_id', NaN, ...
    'local_phi_id_wrapped', NaN, ...
    'local_dx_c_id', NaN, ...
    'local_decision_score', NaN, ...
    'local_weighted_voltage_rmse', NaN, ...
    'local_plain_voltage_rmse', NaN, ...
    'local_status', "", ...
    'local_failure_reason', "", ...
    'joint_reported_EO_id', NaN, ...
    'reported_eo_source', "", ...
    'reported_refit_applied', NaN, ...
    'reported_refit_status', "", ...
    'reported_refit_failure_reason', "", ...
    'reported_vs_local_rmse_delta_v', NaN, ...
    'reported_vs_local_rmse_ratio', NaN, ...
    'fit_stage', "", ...
    'status', "", ...
    'failure_reason', "");
end


function Wres = make_empty_window_result_local()
Wres = struct( ...
    'window_id', NaN, ...
    'lap_range', [NaN NaN], ...
    'time_window_s', [NaN NaN], ...
    'Result', [], ...
    'LocalResult', [], ...
    'ReportedResult', [], ...
    'CoreResult', [], ...
    'ExpandedResult', [], ...
    'ExpansionInfo', [], ...
    'EOSelectionInfo', [], ...
    'CoreCoverage', [], ...
    'ExpandedCoverage', [], ...
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

    if isfield(bundle, 'Coverage') && isstruct(bundle.Coverage)
        trend_row = attach_coverage_to_trend_row_local( ...
            trend_row, bundle.Coverage);
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
result.decision_score = inf;
result.decision_score_mode = '';
result.static_residual_rmse_mm = inf;
result.seed_weighted_voltage_rmse = inf;
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
result.EOAmbiguity = struct();
result.BundleSummary = struct();
result.ReconstructionPreview = [];
end


function Summary = make_empty_joint_reported_refit_summary_local()
Summary = struct( ...
    'applied', false, ...
    'policy', '', ...
    'joint_eo', NaN, ...
    'joint_status', '', ...
    'window_count', 0, ...
    'refit_count', 0, ...
    'kept_local_count', 0, ...
    'failed_refit_count', 0, ...
    'skip_reason', '');
end


function [trend_rows, WindowResult, best, Summary] = apply_joint_reported_eo_refit_local( ...
    trend_rows, WindowResult, FinalBundleStore, FinalSeedTableStore, ...
    WindowPlan, JointEOSummary, best, S)

Summary = make_empty_joint_reported_refit_summary_local();
Summary.policy = getfield_default_local(S, ...
    'joint_reported_refit_policy', 'refit_each_window_to_joint_eo');
Summary.joint_eo = getfield_default_local(JointEOSummary, 'best_EO', NaN);
Summary.joint_status = char(string(getfield_default_local( ...
    JointEOSummary, 'status', 'empty')));

skip_status = string(getfield_default_local(S, ...
    'joint_reported_refit_skip_status', {}));
if isempty(skip_status)
    skip_status = strings(0, 1);
end

if ~isfinite(Summary.joint_eo)
    Summary.skip_reason = 'nonfinite_joint_eo';
    return;
end

if any(strcmpi(string(Summary.joint_status), skip_status))
    Summary.skip_reason = ['joint_status_', Summary.joint_status];
    return;
end

Summary.applied = true;
best = make_failed_result_local('no_joint_reported_window');

for iw = 1:numel(WindowResult)
    local_result = [];
    if isfield(WindowResult(iw), 'Result')
        local_result = WindowResult(iw).Result;
    end

    if isempty(local_result) || ~isstruct(local_result) || ...
            ~isfield(local_result, 'status') || ...
            ~strcmpi(local_result.status, 'ok')
        continue;
    end

    Summary.window_count = Summary.window_count + 1;

    bundle = [];
    seed_table = [];
    if iw <= numel(FinalBundleStore)
        bundle = FinalBundleStore{iw};
    end
    if iw <= numel(FinalSeedTableStore)
        seed_table = FinalSeedTableStore{iw};
    end

    refit_applied = false;
    reported = local_result;
    refit_status = 'kept_local_same_eo';
    refit_failure = '';

    if isempty(bundle) || ~isstruct(bundle) || isempty(seed_table)
        refit_status = 'failed';
        refit_failure = 'missing_final_bundle_or_seed_table';
        Summary.failed_refit_count = Summary.failed_refit_count + 1;
    elseif isfinite(local_result.EO_id) && local_result.EO_id == Summary.joint_eo
        Summary.kept_local_count = Summary.kept_local_count + 1;
    else
        refit_applied = true;
        reported_try = refine_direct_template_fit_local( ...
            bundle, seed_table, Summary.joint_eo, S);
        if isstruct(reported_try) && isfield(reported_try, 'status') && ...
                strcmpi(reported_try.status, 'ok')
            reported = reported_try;
            Summary.refit_count = Summary.refit_count + 1;
            refit_status = 'ok';
        else
            Summary.failed_refit_count = Summary.failed_refit_count + 1;
            refit_status = 'failed';
            refit_failure = getfield_default_local( ...
                reported_try, 'failure_reason', 'joint_refit_failed');
        end
    end

    if iw <= height(WindowPlan)
        W = WindowPlan(iw, :);
    else
        W = table();
    end

    reported = decorate_joint_reported_result_local( ...
        reported, local_result, bundle, seed_table, Summary.joint_eo, ...
        W, iw, S);

    if refit_applied && strcmpi(refit_status, 'failed')
        reported = local_result;
    end

    row = build_trend_row_local(iw, W, reported, bundle);
    row = attach_local_result_to_trend_row_local( ...
        row, local_result, reported, Summary.joint_eo, ...
        refit_applied, refit_status, refit_failure);
    trend_rows(iw) = row;

    WindowResult(iw).LocalResult = local_result;
    WindowResult(iw).ReportedResult = reported;
    WindowResult(iw).Result = reported;

    if strcmpi(reported.status, 'ok')
        score = getfield_default_local(reported, ...
            'decision_score', reported.weighted_voltage_rmse);
        best_score = getfield_default_local(best, ...
            'decision_score', best.weighted_voltage_rmse);
        if isfinite(score) && (~isfinite(best_score) || ...
                score < best_score || ...
                (abs(score - best_score) <= eps(max(abs(best_score), 1)) && ...
                reported.weighted_voltage_rmse < best.weighted_voltage_rmse))
            best = reported;
            best.window_id = iw;
            best.lap_range = [W.lap_start, W.lap_end];
            best.time_window_s = [W.time_start_s, W.time_end_s];
        end
    end
end

if isfield(best, 'status') && strcmpi(best.status, 'failed')
    best = make_failed_result_local('no_valid_joint_reported_window');
end
end


function reported = decorate_joint_reported_result_local( ...
    reported, local_result, bundle, seed_table, joint_eo, W, iw, S)

if isempty(reported) || ~isstruct(reported)
    reported = local_result;
end

copy_fields = {'fit_stage', 'CoreEO', 'ExpandedEO', 'PhaseReferenceEO', ...
    'FinalEOChangedAfterExpansion', 'CoreVPPointCount', 'CoreVPPolicy', ...
    'CoreVPFallback', 'ExpandedVPPointCount', 'ExpandedVPPolicy', ...
    'ExpandedVPFallback', 'VPPointCount', 'VPPolicy', 'VPFallback', ...
    'BundleSummary'};

for i = 1:numel(copy_fields)
    f = copy_fields{i};
    if isfield(local_result, f)
        reported.(f) = local_result.(f);
    end
end

reported.window_id = iw;
if istable(W) && ~isempty(W)
    reported.lap_range = [W.lap_start, W.lap_end];
    reported.time_window_s = [W.time_start_s, W.time_end_s];
end

reported.ReportedEOSelectionPolicy = 'joint_eo_summary_refit';
reported.ReportedEO_id = joint_eo;

if isfield(local_result, 'CandidateTable')
    reported.LocalCandidateTable = local_result.CandidateTable;
    reported.CandidateTable = local_result.CandidateTable;
end

if isfield(local_result, 'EOAmbiguity')
    reported.LocalEOAmbiguity = local_result.EOAmbiguity;
    reported.EOAmbiguity = local_result.EOAmbiguity;
end

if ~isfield(reported, 'VPSeedTable') || isempty(reported.VPSeedTable)
    if ~isempty(seed_table) && isstruct(seed_table)
        reported.VPSeedTable = struct2table(seed_table);
    elseif isfield(local_result, 'VPSeedTable')
        reported.VPSeedTable = local_result.VPSeedTable;
    end
end

if isstruct(bundle) && strcmpi(getfield_default_local(reported, 'status', ''), 'ok')
    reported.ReconstructionPreview = build_reconstruction_preview_local( ...
        bundle, reported, S.store_bundle_preview_points);
end
end


function row = attach_local_result_to_trend_row_local( ...
    row, local_result, reported, joint_eo, refit_applied, ...
    refit_status, refit_failure)

row.local_status = string(getfield_default_local(local_result, 'status', ''));
row.local_failure_reason = string(getfield_default_local( ...
    local_result, 'failure_reason', ''));
row.local_A_id = getfield_default_local(local_result, 'A_id', NaN);
row.local_EO_id = getfield_default_local(local_result, 'EO_id', NaN);
row.local_fn_id = getfield_default_local(local_result, 'fn_id', NaN);
row.local_phi_id_wrapped = getfield_default_local( ...
    local_result, 'phi_id_wrapped', NaN);
row.local_dx_c_id = getfield_default_local(local_result, 'dx_c_id', NaN);
row.local_decision_score = getfield_default_local( ...
    local_result, 'decision_score', NaN);
row.local_weighted_voltage_rmse = getfield_default_local( ...
    local_result, 'weighted_voltage_rmse', NaN);
row.local_plain_voltage_rmse = getfield_default_local( ...
    local_result, 'plain_voltage_rmse', NaN);

row.joint_reported_EO_id = joint_eo;
row.reported_eo_source = "joint_eo_refit";
row.reported_refit_applied = double(refit_applied);
row.reported_refit_status = string(refit_status);
row.reported_refit_failure_reason = string(refit_failure);

reported_rmse = getfield_default_local(reported, ...
    'weighted_voltage_rmse', NaN);
local_rmse = row.local_weighted_voltage_rmse;
row.reported_vs_local_rmse_delta_v = reported_rmse - local_rmse;
row.reported_vs_local_rmse_ratio = row.reported_vs_local_rmse_delta_v ./ ...
    max(abs(local_rmse), eps);
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

    if isfield(bundle, 'Coverage') && isstruct(bundle.Coverage)
        row = attach_coverage_to_trend_row_local(row, bundle.Coverage);
    end
end

if isempty(result) || ~isstruct(result)
    row.status = "failed";
    row.failure_reason = "empty_result";
    return;
end

row.status = string(result.status);
row.failure_reason = string(result.failure_reason);

if isfield(result, 'fit_stage')
    row.fit_stage = string(result.fit_stage);
end

if isfield(result, 'BundleSummary') && isstruct(result.BundleSummary)
    B = result.BundleSummary;
    row.candidate_point_count = getfield_default_local(B, 'candidate_point_count', NaN);
    row.inside_domain_count = getfield_default_local(B, 'inside_domain_count', NaN);
    row.inside_guard_count = getfield_default_local(B, 'inside_guard_count', NaN);
    row.dynamic_effective_count = getfield_default_local(B, 'dynamic_effective_count', NaN);
    row.main_pulse_count = getfield_default_local(B, 'main_pulse_count', NaN);
    row.query_guard_source = string(getfield_default_local(B, 'query_guard_source', ""));
    row.query_guard_mean_mm = getfield_default_local(B, 'query_guard_mean_mm', NaN);
    row.query_safe_left_min_mm = getfield_default_local(B, 'query_safe_left_min_mm', NaN);
    row.query_safe_right_max_mm = getfield_default_local(B, 'query_safe_right_max_mm', NaN);
end

if isfield(result, 'EOAmbiguity') && isstruct(result.EOAmbiguity)
    A = result.EOAmbiguity;
    row.second_best_EO_id = getfield_default_local(A, 'second_best_eo', NaN);
    row.second_best_rmse_v = getfield_default_local(A, 'second_best_rmse_v', NaN);
    row.second_best_decision_score = getfield_default_local(A, 'second_best_decision_score', NaN);
    row.eo_decision_score_gap = getfield_default_local(A, 'decision_score_gap', NaN);
    row.eo_decision_score_gap_ratio = getfield_default_local(A, 'decision_score_gap_ratio', NaN);
    row.eo_rmse_gap_v = getfield_default_local(A, 'rmse_gap_v', NaN);
    row.eo_rmse_gap_ratio = getfield_default_local(A, 'rmse_gap_ratio', NaN);
    row.eo_candidate_count = getfield_default_local(A, 'candidate_count', NaN);
    row.eo_finite_candidate_count = getfield_default_local(A, 'finite_candidate_count', NaN);
    row.eo_top_list = string(getfield_default_local(A, 'top_eo_list', ""));
end

row.core_EO_id = getfield_default_local(result, 'CoreEO', NaN);
row.expanded_EO_id = getfield_default_local(result, 'ExpandedEO', NaN);
row.phase_reference_EO_id = getfield_default_local(result, 'PhaseReferenceEO', NaN);
row.final_EO_changed_after_expansion = ...
    double(getfield_default_local(result, 'FinalEOChangedAfterExpansion', false));
row.vp_point_count = getfield_default_local(result, 'VPPointCount', NaN);
row.core_vp_point_count = getfield_default_local(result, 'CoreVPPointCount', NaN);
row.expanded_vp_point_count = getfield_default_local(result, 'ExpandedVPPointCount', NaN);
row.vp_screening_policy = string(getfield_default_local(result, 'VPPolicy', ""));
row.vp_screening_fallback = string(getfield_default_local(result, 'VPFallback', ""));

if strcmpi(result.status, 'ok')
    row.A_id = result.A_id;
    row.EO_id = result.EO_id;
    row.fn_id = result.fn_id;
    row.phi_id_wrapped = result.phi_id_wrapped;
    row.dx_c_id = result.dx_c_id;
    row.d0_id = result.d0_id;
    row.sensor_eta_max_abs_mm = result.sensor_eta_max_abs_mm;
    row.decision_score = getfield_default_local(result, 'decision_score', NaN);
    row.static_residual_rmse_mm = getfield_default_local(result, 'static_residual_rmse_mm', NaN);
    row.seed_weighted_voltage_rmse = getfield_default_local(result, 'seed_weighted_voltage_rmse', NaN);
    row.weighted_voltage_rmse = result.weighted_voltage_rmse;
    row.plain_voltage_rmse = result.plain_voltage_rmse;
    row.point_count = result.point_count;
    row.valid_segment_count = result.valid_segment_count;
    row.local_A_id = row.A_id;
    row.local_EO_id = row.EO_id;
    row.local_fn_id = row.fn_id;
    row.local_phi_id_wrapped = row.phi_id_wrapped;
    row.local_dx_c_id = row.dx_c_id;
    row.local_decision_score = row.decision_score;
    row.local_weighted_voltage_rmse = row.weighted_voltage_rmse;
    row.local_plain_voltage_rmse = row.plain_voltage_rmse;
    row.local_status = string(result.status);
    row.local_failure_reason = string(result.failure_reason);
    row.reported_eo_source = "realtime_window_local";
    row.reported_refit_applied = 0;

    if isfield(bundle, 'Coverage') && isstruct(bundle.Coverage)
        row = attach_coverage_to_trend_row_local(row, bundle.Coverage);
    end

    if isfield(result, 'core_point_count')
        row.core_point_count = result.core_point_count;
    end

    if isfield(result, 'expanded_point_count')
        row.expanded_point_count = result.expanded_point_count;
    end
end
end


function Coverage = build_bundle_coverage_local(bundle, S)
Coverage = struct();
Coverage.ok = false;
Coverage.failure_reason = 'empty_bundle';
Coverage.valid_pass_count = 0;
Coverage.valid_sensor_count = 0;
Coverage.valid_lap_count = 0;
Coverage.expected_valid_pass_count = S.window_laps * numel(S.analysis_sensors);
Coverage.min_valid_pass_count = getfield_default_local( ...
    S, 'min_valid_pass_count', 1);
Coverage.min_valid_sensor_count = getfield_default_local( ...
    S, 'min_valid_sensor_count', min(2, numel(S.analysis_sensors)));
Coverage.min_valid_lap_count = getfield_default_local( ...
    S, 'min_valid_lap_count', min(2, S.window_laps));
Coverage.valid_pass_fraction = 0;
Coverage.valid_pass_ids = [];
Coverage.valid_sensor_ids = [];
Coverage.valid_laps = [];

if isempty(bundle) || ~isstruct(bundle) || ~isfield(bundle, 'pass_id') || ...
        isempty(bundle.pass_id)
    return;
end

point_valid = true(size(bundle.pass_id(:)));
if isfield(bundle, 'finite_mask') && numel(bundle.finite_mask) == numel(point_valid)
    point_valid = point_valid & bundle.finite_mask(:);
end
if isfield(bundle, 'x') && numel(bundle.x) == numel(point_valid)
    point_valid = point_valid & isfinite(bundle.x(:));
end
if isfield(bundle, 'v') && numel(bundle.v) == numel(point_valid)
    point_valid = point_valid & isfinite(bundle.v(:));
end

Coverage.valid_pass_ids = unique(bundle.pass_id(point_valid & ...
    isfinite(bundle.pass_id(:))));
Coverage.valid_pass_count = numel(Coverage.valid_pass_ids);

if isfield(bundle, 'sensor_id') && numel(bundle.sensor_id) == numel(point_valid)
    Coverage.valid_sensor_ids = unique(bundle.sensor_id(point_valid & ...
        isfinite(bundle.sensor_id(:))));
    Coverage.valid_sensor_count = numel(Coverage.valid_sensor_ids);
end

if isfield(bundle, 'analysis_lap') && numel(bundle.analysis_lap) == numel(point_valid)
    Coverage.valid_laps = unique(bundle.analysis_lap(point_valid & ...
        isfinite(bundle.analysis_lap(:))));
    Coverage.valid_lap_count = numel(Coverage.valid_laps);
end

Coverage.valid_pass_fraction = Coverage.valid_pass_count / ...
    max(Coverage.expected_valid_pass_count, 1);

if Coverage.valid_pass_count < Coverage.min_valid_pass_count
    Coverage.failure_reason = 'too_few_valid_core_passes';
elseif Coverage.valid_sensor_count < Coverage.min_valid_sensor_count
    Coverage.failure_reason = 'too_few_valid_core_sensors';
elseif Coverage.valid_lap_count < Coverage.min_valid_lap_count
    Coverage.failure_reason = 'too_few_valid_core_laps';
else
    Coverage.ok = true;
    Coverage.failure_reason = '';
end
end


function row = attach_coverage_to_trend_row_local(row, Coverage)
if isempty(Coverage) || ~isstruct(Coverage)
    return;
end

row.valid_segment_count = getfield_default_local( ...
    Coverage, 'valid_pass_count', row.valid_segment_count);
row.valid_sensor_count = getfield_default_local( ...
    Coverage, 'valid_sensor_count', row.valid_sensor_count);
row.valid_lap_count = getfield_default_local( ...
    Coverage, 'valid_lap_count', row.valid_lap_count);
row.expected_valid_pass_count = getfield_default_local( ...
    Coverage, 'expected_valid_pass_count', row.expected_valid_pass_count);
row.min_valid_pass_count = getfield_default_local( ...
    Coverage, 'min_valid_pass_count', row.min_valid_pass_count);
row.min_valid_sensor_count = getfield_default_local( ...
    Coverage, 'min_valid_sensor_count', row.min_valid_sensor_count);
row.min_valid_lap_count = getfield_default_local( ...
    Coverage, 'min_valid_lap_count', row.min_valid_lap_count);
row.valid_pass_fraction = getfield_default_local( ...
    Coverage, 'valid_pass_fraction', row.valid_pass_fraction);
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

if isstruct(all_bundle) && isfield(all_bundle, 'query_guard_source') && ~isempty(all_bundle.query_guard_source)
    summary.query_guard_source = strjoin(unique(string(all_bundle.query_guard_source)).', '|');
else
    summary.query_guard_source = '';
end

if isstruct(all_bundle) && isfield(all_bundle, 'query_safe_left_mm') && ~isempty(all_bundle.query_safe_left_mm)
    summary.query_safe_left_min_mm = min(all_bundle.query_safe_left_mm, [], 'omitnan');
    summary.query_safe_right_max_mm = max(all_bundle.query_safe_right_mm, [], 'omitnan');
else
    summary.query_safe_left_min_mm = NaN;
    summary.query_safe_right_max_mm = NaN;
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
    elseif isstring(value)
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
S.overshoot_penalty_weight = 0;
S.sensor_eta_reg_weight_v_per_mm = 0;
S.invalid_query_penalty_scale = 5;
S.template_interp_method = 'pchip';
end


function Joint = build_joint_eo_summary_local(WindowResult, S)
if nargin < 2
    S = struct();
end

joint_gap_threshold = getfield_default_local(S, ...
    'joint_eo_gap_ratio_threshold', 0.03);
missing_window_penalty = getfield_default_local(S, ...
    'joint_missing_window_penalty', 0.50);
min_window_fraction = getfield_default_local(S, ...
    'joint_min_window_fraction', 0.50);

Joint = struct( ...
    'status', 'empty', ...
    'dominant_eo', NaN, ...
    'best_EO', NaN, ...
    'second_EO', NaN, ...
    'mean_objective', NaN, ...
    'second_best_eo', NaN, ...
    'second_best_mean_objective', NaN, ...
    'best_mean_decision_score', NaN, ...
    'second_mean_decision_score', NaN, ...
    'best_mean_weighted_rmse', NaN, ...
    'second_mean_weighted_rmse', NaN, ...
    'best_joint_score', NaN, ...
    'second_joint_score', NaN, ...
    'gap_ratio', NaN, ...
    'gap_ratio_threshold', joint_gap_threshold, ...
    'joint_missing_window_penalty', missing_window_penalty, ...
    'joint_min_window_fraction', min_window_fraction, ...
    'candidate_count', 0, ...
    'window_count', 0, ...
    'valid_window_count', 0, ...
    'candidate_table', table(), ...
    'EOTable', table());

if isempty(WindowResult)
    return;
end

eo_vec = [];
decision_vec = [];
weighted_rmse_vec = [];
objective_vec = [];
window_vec = [];

for iw = 1:numel(WindowResult)
    if ~isfield(WindowResult(iw), 'Result') || isempty(WindowResult(iw).Result) || ...
            ~isfield(WindowResult(iw).Result, 'CandidateTable')
        continue;
    end

    C = WindowResult(iw).Result.CandidateTable;
    if isempty(C) || ~istable(C) || ...
            ~ismember('EO', C.Properties.VariableNames)
        continue;
    end

    if ismember('weighted_voltage_rmse', C.Properties.VariableNames)
        score = C.weighted_voltage_rmse;
    elseif ismember('objective_score', C.Properties.VariableNames)
        score = C.objective_score;
    elseif ismember('decision_score', C.Properties.VariableNames)
        score = C.decision_score;
    elseif ismember('static_residual_rmse_mm', C.Properties.VariableNames)
        score = C.static_residual_rmse_mm;
    else
        continue;
    end

    if ismember('weighted_voltage_rmse', C.Properties.VariableNames)
        weighted_rmse = C.weighted_voltage_rmse;
    else
        weighted_rmse = score;
    end

    if ismember('objective_score', C.Properties.VariableNames)
        objective = C.objective_score;
    else
        objective = score;
    end

    ok = isfinite(C.EO) & isfinite(score) & isfinite(weighted_rmse);
    if ismember('status', C.Properties.VariableNames)
        ok = ok & strcmpi(string(C.status), "ok");
    end

    if ~any(ok)
        continue;
    end

    eo_vec = [eo_vec; C.EO(ok)]; %#ok<AGROW>
    decision_vec = [decision_vec; score(ok)]; %#ok<AGROW>
    weighted_rmse_vec = [weighted_rmse_vec; weighted_rmse(ok)]; %#ok<AGROW>
    objective_vec = [objective_vec; objective(ok)]; %#ok<AGROW>
    window_vec = [window_vec; repmat(iw, nnz(ok), 1)]; %#ok<AGROW>
end

if isempty(eo_vec)
    Joint.status = 'no_valid_candidate_table';
    return;
end

valid_window_count = numel(unique(window_vec));
eo_list = unique(eo_vec(:));

rows = repmat(struct( ...
    'EO', NaN, ...
    'window_count', 0, ...
    'mean_decision_score', NaN, ...
    'median_decision_score', NaN, ...
    'min_decision_score', NaN, ...
    'mean_weighted_rmse', NaN, ...
    'median_weighted_rmse', NaN, ...
    'min_weighted_rmse', NaN, ...
    'mean_objective', NaN, ...
    'median_objective', NaN, ...
    'std_objective', NaN, ...
    'window_fraction', NaN, ...
    'low_coverage_flag', false, ...
    'joint_score', NaN), numel(eo_list), 1);

for i = 1:numel(eo_list)
    eo = eo_list(i);
    m = eo_vec == eo;
    rows(i).EO = eo;
    rows(i).window_count = numel(unique(window_vec(m)));
    rows(i).mean_decision_score = mean(decision_vec(m), 'omitnan');
    rows(i).median_decision_score = median(decision_vec(m), 'omitnan');
    rows(i).min_decision_score = min(decision_vec(m), [], 'omitnan');
    rows(i).mean_weighted_rmse = mean(weighted_rmse_vec(m), 'omitnan');
    rows(i).median_weighted_rmse = median(weighted_rmse_vec(m), 'omitnan');
    rows(i).min_weighted_rmse = min(weighted_rmse_vec(m), [], 'omitnan');
    rows(i).mean_objective = mean(objective_vec(m), 'omitnan');
    rows(i).median_objective = median(objective_vec(m), 'omitnan');
    rows(i).std_objective = std(objective_vec(m), 0, 'omitnan');
    rows(i).window_fraction = rows(i).window_count / max(valid_window_count, 1);
    rows(i).low_coverage_flag = rows(i).window_fraction < min_window_fraction;
    rows(i).joint_score = rows(i).mean_decision_score .* ...
        (1 + missing_window_penalty .* (1 - rows(i).window_fraction));
end

T = struct2table(rows);
T = sortrows(T, { ...
    'low_coverage_flag', ...
    'joint_score', ...
    'median_decision_score', ...
    'EO'}, ...
    {'ascend', 'ascend', 'ascend', 'ascend'});

Joint.status = 'ok';
Joint.dominant_eo = T.EO(1);
Joint.best_EO = T.EO(1);
Joint.best_mean_decision_score = T.mean_decision_score(1);
Joint.best_mean_weighted_rmse = T.mean_weighted_rmse(1);
Joint.best_joint_score = T.joint_score(1);
Joint.mean_objective = T.mean_objective(1);
Joint.candidate_count = height(T);
Joint.window_count = valid_window_count;
Joint.valid_window_count = valid_window_count;
Joint.candidate_table = T;
Joint.EOTable = T;

if T.low_coverage_flag(1)
    Joint.status = 'low_window_coverage_ambiguous';
end

if height(T) >= 2
    Joint.second_EO = T.EO(2);
    Joint.second_best_eo = T.EO(2);
    Joint.second_mean_decision_score = T.mean_decision_score(2);
    Joint.second_mean_weighted_rmse = T.mean_weighted_rmse(2);
    Joint.second_best_mean_objective = T.mean_objective(2);
    Joint.second_joint_score = T.joint_score(2);
    Joint.gap_ratio = (T.joint_score(2) - T.joint_score(1)) ./ ...
        max(abs(T.joint_score(1)), eps);

    if Joint.gap_ratio < joint_gap_threshold && ...
            ~strcmpi(Joint.status, 'low_window_coverage_ambiguous')
        Joint.status = 'low_gap_ambiguous';
    end
end
end


function Summary = build_resonance_summary_local(Trend, JointEOSummary)
Summary = struct();
if nargin < 2 || isempty(JointEOSummary)
    JointEOSummary = struct();
end

if isempty(Trend) || height(Trend) == 0
    Summary.valid_window_count = 0;
    Summary.best_window_id = NaN;
    Summary.dominant_eo = NaN;
    Summary.trend_dominant_eo = NaN;
    Summary.joint_best_eo = NaN;
    Summary.joint_status = 'empty';
    Summary.joint_eo_status = "";
    Summary.joint_eo_gap_ratio = NaN;
    Summary.reported_eo = NaN;
    Summary.reported_freq_hz = NaN;
    Summary.reported_eo_source = 'none';
    Summary.mean_freq_hz = NaN;
    Summary.median_freq_hz = NaN;
    return;
end

valid = strcmpi(string(Trend.status), 'ok') & ...
        isfinite(Trend.weighted_voltage_rmse);

Summary.valid_window_count = nnz(valid);
Summary.joint_best_eo = getfield_default_local(JointEOSummary, 'best_EO', NaN);
Summary.joint_status = getfield_default_local(JointEOSummary, 'status', 'empty');
Summary.joint_eo_status = string(Summary.joint_status);
Summary.joint_eo_gap_ratio = getfield_default_local(JointEOSummary, 'gap_ratio', NaN);

if ~any(valid)
    Summary.best_window_id = NaN;
    Summary.dominant_eo = NaN;
    Summary.trend_dominant_eo = NaN;
    Summary.reported_eo = NaN;
    Summary.reported_freq_hz = NaN;
    Summary.reported_eo_source = 'no_valid_realtime_trend';
    Summary.mean_freq_hz = NaN;
    Summary.median_freq_hz = NaN;
    return;
end

valid_rows = find(valid);

[~, local_best] = min(Trend.weighted_voltage_rmse(valid));

Summary.best_window_id = Trend.window_id(valid_rows(local_best));

eos = Trend.EO_id(valid);

Summary.trend_dominant_eo = mode(eos);
Summary.dominant_eo = Summary.trend_dominant_eo;
Summary.mean_freq_hz = mean(Trend.fn_id(valid), 'omitnan');
Summary.median_freq_hz = median(Trend.fn_id(valid), 'omitnan');
Summary.reported_eo = Summary.trend_dominant_eo;
Summary.reported_eo_source = 'realtime_trend_mode';

Summary.reported_freq_hz = compute_reported_frequency_local( ...
    Trend, valid, Summary.reported_eo);
end


function f_hz = compute_reported_frequency_local(Trend, valid, reported_eo)
f_hz = NaN;

if isempty(Trend) || ~isfinite(reported_eo)
    return;
end

same_eo = valid & isfinite(Trend.EO_id) & Trend.EO_id == reported_eo;

if any(same_eo) && ismember('fn_id', Trend.Properties.VariableNames)
    f_hz = median(Trend.fn_id(same_eo), 'omitnan');
end

if isfinite(f_hz)
    return;
end

if ismember('rot_freq_mean_hz', Trend.Properties.VariableNames)
    rot_hz = median(Trend.rot_freq_mean_hz(valid), 'omitnan');
    if isfinite(rot_hz)
        f_hz = reported_eo * rot_hz;
    end
end
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
