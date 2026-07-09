clc; close all;

%STEP05_SINGLESYNC_DIRECTTEMPLATE_IDENTIFICATION_20241106
% VERSION: FOUNDATION_MAINPULSE_ADAPTIVE_FIXEDJOINTETA_20260707
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
%   5. The formal 20241106 foundation route uses a fixed residual static
%      sensor eta_s calibrated from dynamic NoEta residuals. NoEta/free-eta
%      and Step04-xc variants belong in comparison/diagnostic scripts.
%   6. Optional expanded bundle uses active-model phase-safe x_query:
%          x_query = x - dx_c - eta_s - A*sin(EO*theta + phi)
%      and only expands per-pass main pulse points whose x_query stays
%      inside the template domain.
%
% Forward model:
%   V_i ~= T_{sensor,blade}(x_i - dx_c - eta_s - A*sin(EO*theta_i + phi))
%
% Usage:
%   Run this script directly after Step02 and Step04.

cfg = BTTProjectConfig_20241106();

%% Step05 local settings
S = struct();
S.method_name = 'foundation_mainpulse_adaptive_fixed_joint_static_eta';
S.version = 'FOUNDATION_MAINPULSE_ADAPTIVE_FIXEDJOINTETA_20260707';
S.method_file_tag = 'FoundationMainPulseAdaptiveFixedJointEta';

S.show_plots = true;
S.save_figures = cfg.save_figures;
S.force_rebuild = false;
S.target_cases = cfg.dynamic_cases;
S.target_blades = cfg.step05_target_blades;
S.analysis_sensors = cfg.step05_analysis_sensors;
S.reference_sensor_id = S.analysis_sensors(1);

S.analysis_start_time_s = cfg.step05_analysis_start_time_s;
S.target_laps = cfg.step05_target_laps;
S.window_laps = cfg.step05_window_laps;
S.sliding_step_laps = cfg.step05_sliding_step_laps;
S.debug_max_windows = inf;

S.output_dir = cfg.step05_output_dir;
S.figure_dir = cfg.step05_figure_dir;

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

force_rebuild_env = strtrim(getenv('STEP05_FORCE_REBUILD'));
if ~isempty(force_rebuild_env)
    S.force_rebuild = any(strcmpi(force_rebuild_env, {'1', 'true', 'yes', 'on'}));
end

debug_max_windows_env = strtrim(getenv('STEP05_DEBUG_MAX_WINDOWS'));
if ~isempty(debug_max_windows_env)
    parsed_debug_windows = str2double(debug_max_windows_env);
    if isfinite(parsed_debug_windows) && parsed_debug_windows > 0
        S.debug_max_windows = parsed_debug_windows;
    end
end

target_blades_env = strtrim(getenv('STEP05_TARGET_BLADES'));
if ~isempty(target_blades_env)
    blade_tokens = regexp(target_blades_env, '[,;\s]+', 'split');
    parsed_blades = [];
    for iBladeToken = 1:numel(blade_tokens)
        if isempty(blade_tokens{iBladeToken})
            continue;
        end
        blade_value = str2double(blade_tokens{iBladeToken});
        if isfinite(blade_value)
            parsed_blades(end+1) = round(blade_value); %#ok<AGROW>
        end
    end
    if ~isempty(parsed_blades)
        S.target_blades = unique(parsed_blades, 'stable');
    end
end

%% Step05 identification settings
S.freq_search_hz = cfg.step05_freq_search_hz;
S.eo_pad = 2;

S.eo_candidate_policy = 'adaptive_full_fit';
S.vp_top_k_min = 3;
S.vp_top_k_max = 7;
S.vp_gap_ratio_keep = 0.05;
S.include_neighbor_eo = true;
S.neighbor_eo_radius = 1;
S.vp_include_sensor_eta = false;
S.vp_seed_observation_mode = 'gradient_displacement';
S.vp_gradient_min_ratio = 0.10;
S.vp_gradient_reference_quantile = 95;
S.vp_gradient_weight_power = 2;
S.fit_sensor_eta = false;
S.eta_model_policy = 'single';
S.active_eta_model = 'fixed_joint_static_eta';
S.eta_models = {S.active_eta_model};
S.eta_reference_sensor = S.analysis_sensors(1);
S.sensor_eta_mode = 'fixed_static';
S.static_eta_source = 'joint_dynamic_residual_eta_preview';
S.static_eta_file = '';
S.sensor_eta_prior_mm = zeros(1, numel(S.analysis_sensors));
S.joint_missing_window_penalty = 0.50;
S.joint_min_window_fraction = 0.50;
S.joint_eo_gap_ratio_threshold = 0.03;
S.eo_low_gap_ratio_threshold = 0.03;

S.amplitude_limit_mm = 0.50;
S.dx_c_limit_mm = 0.35;
S.sensor_eta_limit_mm = cfg.step05_sensor_eta_limit_mm;
S.sensor_eta_reg_weight_v_per_mm = cfg.step05_sensor_eta_reg_weight_v_per_mm;

S.dynamic_window_mode = get_cfg_text_local( ...
    cfg, 'step05_dynamic_window_mode', 'legacy_row_bounds');
S.pulse_pad_sec = get_cfg_scalar_local(cfg, 'step05_pulse_pad_sec', 0);
S.pulse_window_sec = 6e-4;
S.pulse_pad_fraction = 0.30;
S.pulse_min_pad_points = 50;

S.core_mask_mode = 'per_pass_main_pulse_query_safe';
S.domain_margin_mm = 0.02;
S.domain_selection_mode = 'hard';
S.domain_soft_margin_mm = 0.10;

S.query_guard_mode = 'step05_adaptive_from_x_domain';
S.query_guard_mm = [];
S.query_guard_min_mm = 0.12;
S.query_guard_max_mm = S.amplitude_limit_mm + S.dx_c_limit_mm + ...
    S.sensor_eta_limit_mm + 0.05;
S.query_guard_safety_mm = 0.05;
S.query_guard_quantile = 95;
S.sensor_query_guard_scale_table = [5 0.80; 7 0.75];

S.phase_safe_expansion = true;
S.phase_safe_margin_mm = get_cfg_scalar_local( ...
    cfg, 'step05_phase_safe_margin_mm', 0.02);
S.phase_safe_reference_mode = 'core_current';
S.choose_better_pass = true;
S.better_pass_rmse_tolerance_v = 0;

%% Step05 diagnostic, gate, and optimizer settings
S.pulse_mode = 'single';
S.default_sensor_threshold = 0.5;
S.dynamic_effective_mode = 'gradient';
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

S.fminsearch_max_iter = 800;
S.fminsearch_max_fun = 2000;
S.overshoot_penalty_weight = 100;
S.objective_mode = 'rmse';
S.method_preset = 'foundation_mainpulse_adaptive_fixed_joint_static_eta';
S.store_bundle_preview_points = 6000;
S.raw_file_cache_enable = true;
S.raw_file_cache_max_files = 12;

store_bundle_preview_points_env = strtrim(getenv('STEP05_STORE_BUNDLE_PREVIEW_POINTS'));
if ~isempty(store_bundle_preview_points_env)
    if any(strcmpi(store_bundle_preview_points_env, {'inf', 'infinite', 'all'}))
        S.store_bundle_preview_points = inf;
    else
        parsed_preview_points = str2double(store_bundle_preview_points_env);
        if isfinite(parsed_preview_points) && parsed_preview_points > 0
            S.store_bundle_preview_points = parsed_preview_points;
        else
            error('STEP05_STORE_BUNDLE_PREVIEW_POINTS must be positive, inf, infinite, or all.');
        end
    end
end

method_preset_env = strtrim(getenv('STEP05_METHOD_PRESET'));
if ~isempty(method_preset_env)
    error(['STEP05_METHOD_PRESET belongs to archived 20241106 Step05 routes. ' ...
        'The official foundation Step05 uses adaptive_full_fit with fixed joint static eta.']);
end

static_eta_file_env = strtrim(getenv('STEP05_STATIC_ETA_FILE'));
if ~isempty(static_eta_file_env)
    S.static_eta_file = static_eta_file_env;
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
fprintf('Eta policy: %s, active model=%s, reference CH%d.\n', ...
    S.eta_model_policy, S.active_eta_model, S.eta_reference_sensor);
fprintf('EO candidate policy: %s (VP all EO, adaptive full-fit set).\n', ...
    S.eo_candidate_policy);
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
        nan(numel(S.analysis_sensors), 1), ...
        strings(numel(S.analysis_sensors), 1), ...
        'VariableNames', {'sensor_id', 'eta_prior_mm', ...
        'eta_plan_iqr_mm', 'joint_static_eta_file'});
end

template_mode = get_template_source_mode_local(Template);
fprintf('EO selection: pure residual ranking; no preferred EO tie-break. Template source: %s\n', template_mode);

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


function [eta, info] = load_joint_static_eta_prior_local(cfg, S)
sensors = S.analysis_sensors(:).';
blade_id = S.target_blades(1);
sensor_tag = ['S', sprintf('%d', sensors)];

eta_file = '';
if isfield(S, 'static_eta_file') && ~isempty(S.static_eta_file)
    eta_file = char(S.static_eta_file);
end
if isempty(eta_file)
    eta_file = fullfile(cfg.output_root, 'step05_joint_static_eta_preview', ...
        sprintf('JointStaticEtaPreview_B%d_%s.mat', blade_id, sensor_tag));
end
if ~isfile(eta_file)
    error(['Formal Step05 now requires dynamic residual static eta. ' ...
        'Run Calibrate_Step05_JointStaticEta_FromPreview_20241106 first, ' ...
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
template_source = get_cfg_text_local(cfg, 'step05_template_source', 'step04');

if strcmpi(template_source, 'legacy_newflow')
    legacy_file = get_cfg_text_local( ...
        cfg, 'step05_legacy_template_library_file', '');

    if ~isempty(legacy_file) && isfile(legacy_file)
        Template = load_legacy_newflow_template_library_local(cfg, legacy_file);
        fprintf('Step05 template source: legacy NewFlow\n  %s\n', ...
            legacy_file);
        return;
    end

    warning('Legacy NewFlow template source requested, but file is missing: %s. Falling back to Step04 template.', ...
        legacy_file);
end

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
    error('No Step04 template file found under %s. Please run Step04_Build_LowSpeed_OPRCenterStd_Template_20241106.m first.', ...
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


function Template = load_legacy_newflow_template_library_local(cfg, legacy_file)
loaded = load(legacy_file, 'LowSpeedTemplateLibrary');

if ~isfield(loaded, 'LowSpeedTemplateLibrary')
    error('Legacy template file does not contain LowSpeedTemplateLibrary: %s', ...
        legacy_file);
end

L = loaded.LowSpeedTemplateLibrary;

if ~isfield(L, 'entry') || isempty(L.entry)
    error('Legacy LowSpeedTemplateLibrary has no entry records: %s', ...
        legacy_file);
end

SensorBlade = repmat(make_empty_legacy_template_entry_local(), 0, 1);

for ie = 1:numel(L.entry)
    blade_id = L.entry(ie).blade_id;
    bundle_file = char(L.entry(ie).template_file);

    if ~isfile(bundle_file)
        error('Legacy template bundle missing: %s', bundle_file);
    end

    loaded_bundle = load(bundle_file, 'TemplateBundle');

    if ~isfield(loaded_bundle, 'TemplateBundle') || ...
            ~isfield(loaded_bundle.TemplateBundle, 'SensorFiles')
        error('Invalid legacy template bundle: %s', bundle_file);
    end

    files = loaded_bundle.TemplateBundle.SensorFiles;

    for ir = 1:height(files)
        sensor_id = files.sensor_id(ir);
        template_file = char(files.template_file(ir));

        if ~isfile(template_file)
            error('Legacy sensor template missing: %s', template_file);
        end

        loaded_sensor = load(template_file);
        template_sensor = extract_legacy_template_sensor_local( ...
            loaded_sensor, sensor_id, template_file);

        entry = make_empty_legacy_template_entry_local();
        entry.sensor_id = sensor_id;
        entry.blade_id = blade_id;
        entry.x_grid = template_sensor.x_grid(:);
        entry.v_grid = template_sensor.v_grid(:);
        entry.v_grid_raw = template_sensor.v_grid(:);
        entry.v_grid_baseline_removed = template_sensor.v_grid(:);
        entry.dv_dx = template_sensor.dv_dx(:);
        entry.weight_grid = get_optional_vector_local( ...
            template_sensor, 'bin_weight', ones(size(entry.x_grid)));
        entry.count_grid = nan(size(entry.x_grid));
        entry.valid_grid_mask = isfinite(entry.x_grid) & ...
            isfinite(entry.v_grid) & isfinite(entry.dv_dx);
        entry.domain_effective_mask = entry.valid_grid_mask;
        entry.domain_mask = entry.x_grid >= template_sensor.x_domain(1) & ...
            entry.x_grid <= template_sensor.x_domain(2);
        entry.x_domain = template_sensor.x_domain(:).';
        entry.baseline = get_optional_scalar_local(template_sensor, 'baseline', NaN);
        entry.threshold = get_optional_scalar_local(template_sensor, 'threshold', 0.5);
        entry.xc = get_optional_scalar_local(template_sensor, 'xc', 0);
        entry.amplitude = max(entry.v_grid, [], 'omitnan') - ...
            min(entry.v_grid, [], 'omitnan');
        entry.pulse_count = get_optional_scalar_local(template_sensor, 'lap_count', NaN);
        entry.point_count = get_optional_scalar_local(template_sensor, 'point_count', NaN);
        entry.quality_status = 'legacy_newflow';
        entry.x_points_preview = [];
        entry.v_points_preview = [];
        SensorBlade(end+1, 1) = entry; %#ok<AGROW>
    end
end

Template = struct();
Template.Dataset = cfg.dataset;
Template.Case_Name = cfg.low_speed_case;
Template.CreatedBy = mfilename;
Template.CreatedOn = datestr(now, 31);
Template.LowSpeed_Template_Source_Mode = 'legacy_newflow_template_library';
Template.Legacy_Template_Library_File = legacy_file;
Template.Sensor_IDs = unique([SensorBlade.sensor_id], 'stable');
Template.Blades_Num = cfg.blades_num;
Template.R_Tip_mm = cfg.r_tip_mm;
Template.Standard_Relative_Angles = load_step05_standard_angles_local(cfg);
Template.SensorBlade = SensorBlade;
Template.Metadata = struct( ...
    'source', 'legacy_newflow', ...
    'legacy_template_library_file', legacy_file);
end


function mode = get_template_source_mode_local(Template)
mode = 'unknown';
if isfield(Template, 'LowSpeed_Template_Source_Mode') && ...
        ~isempty(Template.LowSpeed_Template_Source_Mode)
    mode = char(string(Template.LowSpeed_Template_Source_Mode));
elseif isfield(Template, 'Metadata') && isfield(Template.Metadata, 'source') && ...
        ~isempty(Template.Metadata.source)
    mode = char(string(Template.Metadata.source));
end
end


function entry = make_empty_legacy_template_entry_local()
entry = struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'x_grid', [], ...
    'v_grid', [], ...
    'v_grid_raw', [], ...
    'v_grid_baseline_removed', [], ...
    'dv_dx', [], ...
    'weight_grid', [], ...
    'count_grid', [], ...
    'valid_grid_mask', [], ...
    'domain_effective_mask', [], ...
    'domain_mask', [], ...
    'x_domain', [NaN NaN], ...
    'baseline', NaN, ...
    'threshold', NaN, ...
    'xc', 0, ...
    'amplitude', NaN, ...
    'pulse_count', NaN, ...
    'point_count', NaN, ...
    'quality_status', '', ...
    'x_points_preview', [], ...
    'v_points_preview', []);
end


function template_sensor = extract_legacy_template_sensor_local( ...
    loaded_sensor, sensor_id, template_file)
if isfield(loaded_sensor, 'Template')
    template = loaded_sensor.Template;
elseif isfield(loaded_sensor, 'TemplateSingle')
    template = loaded_sensor.TemplateSingle;
elseif isfield(loaded_sensor, 'sensor_template')
    template = loaded_sensor.sensor_template;
else
    error('Unsupported legacy sensor template file: %s', template_file);
end

if ~isfield(template, 'Sensor') || isempty(template.Sensor)
    error('Legacy template has no Sensor field: %s', template_file);
end

template_sensor = template.Sensor;

if numel(template_sensor) > 1
    template_sensor = template_sensor([template_sensor.sensor_id] == sensor_id);
end

if isempty(template_sensor)
    error('Legacy template file does not contain CH%d: %s', ...
        sensor_id, template_file);
end
end


function value = get_optional_vector_local(S, field, default_value)
if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
    value = S.(field)(:);
else
    value = default_value(:);
end
end


function value = get_optional_scalar_local(S, field, default_value)
if isstruct(S) && isfield(S, field) && ...
        isscalar(S.(field)) && isfinite(S.(field))
    value = S.(field);
else
    value = default_value;
end
end


function angles = load_step05_standard_angles_local(cfg)
angles = nan(max(cfg.sensor_ids), cfg.blades_num);
sensor_config_file = fullfile( ...
    cfg.step01_output_dir, 'Sensor_Config_20241106.mat');

if ~isfile(sensor_config_file)
    return;
end

loaded = load(sensor_config_file, 'Sensor_Config');

if ~isfield(loaded, 'Sensor_Config')
    return;
end

SC = loaded.Sensor_Config;

if isfield(SC, 'Standard_Relative_Angles_OPRCenter')
    A = SC.Standard_Relative_Angles_OPRCenter;
elseif isfield(SC, 'Standard_Relative_Angles')
    A = SC.Standard_Relative_Angles;
else
    return;
end

n1 = min(size(angles, 1), size(A, 1));
n2 = min(size(angles, 2), size(A, 2));
angles(1:n1, 1:n2) = A(1:n1, 1:n2);
end


function blade_results = process_step05_case_local(case_name, cfg, S, Template)
step02_case_dir = fullfile(cfg.step02_output_dir, case_name);

if ~isfolder(step02_case_dir) || ...
        ~isfile(fullfile(step02_case_dir, 'Step02_Dynamic_BTT_Extraction_20241106.mat'))

    error('Step02 products missing for %s. Please run Step02_Extract_Dynamic_BTT_20241106.m first.', case_name);
end

loaded = load(fullfile(step02_case_dir, ...
    'Step02_Dynamic_BTT_Extraction_20241106.mat'), ...
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

opr_times = load_dynamic_opr_times_local(step02_case_dir, cfg);
opr_events_per_revolution = resolve_step05_opr_events_per_revolution_local( ...
    cfg, Template, step02_metadata);
cfg.opr_pulses_per_rev = opr_events_per_revolution;

[F_omega_deg_s, speed_diag] = build_dynamic_speed_interpolant_local( ...
    opr_times, opr_events_per_revolution, cfg);
speed_diag.opr_events_per_revolution = opr_events_per_revolution;

blade_results = repmat(struct('blade_id', NaN, 'Result', []), ...
    numel(S.target_blades), 1);

for ib = 1:numel(S.target_blades)
    blade_id = S.target_blades(ib);

    fprintf('\n>>> [Step05][%s] Target blade %d\n', case_name, blade_id);

    out_dir = fullfile(S.output_dir, case_name);
    fig_case_dir = fullfile(S.figure_dir, case_name);

    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end

    if exist(fig_case_dir, 'dir') ~= 7
        mkdir(fig_case_dir);
    end

    sensor_tag = ['S', sprintf('%d', S.analysis_sensors)];
    method_file_tag = S.method_file_tag;

    result_file = fullfile(out_dir, sprintf( ...
        'Result_Step05_%s_B%d_%s_%s.mat', ...
        method_file_tag, blade_id, sensor_tag, cfg.dataset));

    trend_file = fullfile(out_dir, sprintf( ...
        'Trend_Step05_%s_B%d_%s_%s.csv', ...
        method_file_tag, blade_id, sensor_tag, cfg.dataset));

    if isfile(result_file) && isfile(trend_file) && ~S.force_rebuild
        D = load(result_file, 'Result');
        if isfield(D, 'Result') && can_reuse_step05_result_local(D.Result, S)
            Result = D.Result;
            fprintf('Reuse existing Step05 result:\n  %s\n', result_file);
            if S.show_plots || S.save_figures
                plot_step05_result_local(Result, S, fig_case_dir, S.show_plots);
            end
            blade_results(ib).blade_id = blade_id;
            blade_results(ib).Result = Result;
            continue;
        end
        fprintf('Existing Step05 result is stale for current eta prior or method; rebuilding:\n  %s\n', result_file);
    end

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

    seed_table_core = solve_vp_seed_eo_scan_local(core_bundle, eo_candidates, S);
    [final_eos_core, eo_selection_core] = ...
        select_eo_candidates_adaptive_local(seed_table_core, S);

    fit_core = refine_direct_template_fit_local( ...
        core_bundle, seed_table_core, final_eos_core, S);

    fit_core.fit_stage = 'core_query_safe';
    fit_core.core_point_count = core_bundle.point_count;
    fit_core.expanded_point_count = core_bundle.point_count;
    fit_core.EOSelectionInfo = eo_selection_core;

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
        expanded_coverage = build_bundle_coverage_local(expanded_bundle, S);
        expanded_bundle.Coverage = expanded_coverage;

        expand_info.applied = true;
        expand_info.expanded_point_count = expanded_bundle.point_count;
        expand_info.expanded_valid_pass_count = expanded_coverage.valid_pass_count;
        expand_info.expanded_valid_sensor_count = expanded_coverage.valid_sensor_count;
        expand_info.expanded_valid_lap_count = expanded_coverage.valid_lap_count;
        expand_info.reason = 'phase_safe_from_core_current';

        if expanded_bundle.point_count >= S.min_fit_points && ...
                expanded_coverage.ok

            seed_table_expanded = solve_vp_seed_eo_scan_local( ...
                expanded_bundle, eo_candidates, S);

            [final_eos_expanded, eo_selection_expanded] = ...
                select_eo_candidates_adaptive_local(seed_table_expanded, S);
            expand_info.ExpandedEOSelectionInfo = eo_selection_expanded;

            fit_expanded = refine_direct_template_fit_local( ...
                expanded_bundle, seed_table_expanded, final_eos_expanded, S);

            fit_expanded.fit_stage = 'phase_safe_expanded';
            fit_expanded.core_point_count = core_bundle.point_count;
            fit_expanded.expanded_point_count = expanded_bundle.point_count;
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
    if strcmpi(result.fit_stage, 'phase_safe_expanded') && ...
            isfield(expand_info, 'ExpandedEOSelectionInfo')
        result.EOSelectionInfo = expand_info.ExpandedEOSelectionInfo;
    else
        result.EOSelectionInfo = eo_selection_core;
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

Trend = struct2table(trend_rows, 'AsArray', true);

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
Result.JointEOSummary = build_joint_eo_summary_local(WindowResult, S);
Result.ResonanceSummary = build_resonance_summary_local( ...
    Trend, Result.JointEOSummary);
end


function tf = can_reuse_step05_result_local(Result, S)
tf = false;
if ~isstruct(Result) || ~isfield(Result, 'RunInfo') || ...
        ~isstruct(Result.RunInfo) || ~isfield(Result.RunInfo, 'step_settings')
    return;
end
oldS = Result.RunInfo.step_settings;
if ~isfield(oldS, 'method_name') || ~strcmpi(char(oldS.method_name), S.method_name)
    return;
end
if ~isfield(oldS, 'version') || ~strcmpi(char(oldS.version), S.version)
    return;
end
if ~isfield(oldS, 'sensor_eta_prior_mm') || ...
        numel(oldS.sensor_eta_prior_mm) ~= numel(S.sensor_eta_prior_mm)
    return;
end
delta = abs(double(oldS.sensor_eta_prior_mm(:).') - ...
    double(S.sensor_eta_prior_mm(:).'));
tf = all(delta < 1e-10);
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
method.forward_model = 'V = T_{s,b}(x - dx_c - eta_s - A*sin(EO*theta + phi))';
method.vp_seed_model = 'VP excludes eta_s; full fit uses fixed dynamic-residual eta_s prior';
method.eo_selection = sprintf(['adaptive_full_fit: VP all EO, retain top-%d minimum, ' ...
    'near-best gap <= %.3g, neighbors radius=%d, capped at %d full-fit EO'], ...
    S.vp_top_k_min, S.vp_gap_ratio_keep, S.neighbor_eo_radius, S.vp_top_k_max);
method.eta_model_policy = S.eta_model_policy;
method.eta_models = S.eta_models;
method.static_eta_source = getfield_default_local(S, 'static_eta_source', '');
method.sensor_eta_prior_mm = getfield_default_local(S, 'sensor_eta_prior_mm', []);
method.eta_gauge = sprintf('eta_s fixed to zero for reference CH%d; other CH use dynamic residual static eta prior', S.analysis_sensors(1));
method.no_template_center_subtraction = true;
method.no_sensor_affine_voltage_projection = true;
method.no_within_pulse_uniform_downsampling = true;
method.core_mask_mode = S.core_mask_mode;
method.phase_safe_expansion = S.phase_safe_expansion;
method.two_stage_region_policy = ['full raw candidate waveform -> query-safe core region -> ' ...
    'optional phase-safe expanded region'];
method.frequency_output_policy = ['ResonanceSummary.reported_eo uses the joint EO summary first; ' ...
    'Trend mode is retained as a diagnostic fallback'];
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


function epr = resolve_step05_opr_events_per_revolution_local(cfg, Template, step02_metadata)
values = [];

if isstruct(cfg)
    names = {'step05_opr_events_per_revolution', ...
        'step04_opr_events_per_revolution', 'opr_pulses_per_rev'};
    for i = 1:numel(names)
        if isfield(cfg, names{i}) && isnumeric(cfg.(names{i})) && ...
                isscalar(cfg.(names{i})) && isfinite(cfg.(names{i}))
            values(end+1) = cfg.(names{i}); %#ok<AGROW>
        end
    end
end

if isstruct(Template)
    names = {'OPR_Events_Per_Revolution', 'opr_events_per_revolution'};
    for i = 1:numel(names)
        if isfield(Template, names{i}) && isnumeric(Template.(names{i})) && ...
                isscalar(Template.(names{i})) && isfinite(Template.(names{i}))
            values(end+1) = Template.(names{i}); %#ok<AGROW>
        end
    end
end

if isstruct(step02_metadata)
    names = {'OPR_Events_Per_Revolution', 'opr_events_per_revolution', ...
        'opr_pulses_per_rev'};
    for i = 1:numel(names)
        if isfield(step02_metadata, names{i}) && ...
                isnumeric(step02_metadata.(names{i})) && ...
                isscalar(step02_metadata.(names{i})) && ...
                isfinite(step02_metadata.(names{i}))
            values(end+1) = step02_metadata.(names{i}); %#ok<AGROW>
        end
    end
end

values = values(isfinite(values) & values > 0);
if isempty(values)
    error('Cannot resolve Step05 OPR events/rev for 20241106.');
end

epr = round(values(1));

if strcmp(string(get_cfg_text_local(cfg, 'dataset', '')), "20241106") && epr ~= 1
    error('20241106 Step05 must use 1 OPR event/rev, got %g.', epr);
end
end


function [F_omega_deg_s, diagnostic] = build_dynamic_speed_interpolant_local( ...
    opr_times, opr_pulses_per_rev, cfg)
opr_times = opr_times(:);

omega_file = '';
if isstruct(cfg) && isfield(cfg, 'step05_legacy_pulse_dir') && ...
        ~isempty(cfg.step05_legacy_pulse_dir)
    omega_file = fullfile(cfg.step05_legacy_pulse_dir, 'omega.mat');
end

if isfile(omega_file)
    omega_vars = whos('-file', omega_file);

    if any(strcmp({omega_vars.name}, 'omega'))
        omega_loaded = load(omega_file, 'omega');
        omega = omega_loaded.omega;

        if ~isempty(omega) && size(omega, 2) >= 2
        valid_omega = isfinite(omega(:,1)) & isfinite(omega(:,2));

        if nnz(valid_omega) >= 2
            F_omega_deg_s = griddedInterpolant( ...
                omega(valid_omega, 1), omega(valid_omega, 2), ...
                'linear', 'nearest');

            diagnostic = struct();
            diagnostic.speed_time_s = omega(valid_omega, 1);
            diagnostic.speed_deg_s = omega(valid_omega, 2);
            diagnostic.rpm = omega(valid_omega, 2) / 360 * 60;
            diagnostic.source = 'legacy_omega_mat';
            diagnostic.omega_file = omega_file;
            return;
        end
        end
    end
end

idx0 = 1:(numel(opr_times) - opr_pulses_per_rev);
idx1 = idx0 + opr_pulses_per_rev;

rev_period = opr_times(idx1) - opr_times(idx0);
speed_time = 0.5 .* (opr_times(idx0) + opr_times(idx1));
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
diagnostic.source = 'opr_period';
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

        if isfield(tpl, 'xc') && isscalar(tpl.xc) && isfinite(tpl.xc)
            x_mm = x_mm - tpl.xc;
        end

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

        threshold = get_step05_sensor_threshold_local(tpl, S);

        dynamic_effective = build_dynamic_effective_mask_local( ...
            x_mm, v_seg, t_seg, template_dv_dx, threshold, S);

        main_pulse_mask = build_main_pulse_mask_local( ...
            v_seg, threshold, S);

        % -------------------------------------------------------------
        % 5. Query-safe shrink.
        %    Since final query coordinate is:
        %       x_query = x - dx_c - A*sin(...)
        %    the first/core pass must use a smaller region than x_domain.
        % -------------------------------------------------------------
        base_for_guard = finite_mask & inside_domain;

        query_guard = compute_query_guard_local( ...
            x_mm, v_seg, tpl, S, base_for_guard, sid);

        inside_guard = x_mm >= tpl.x_domain(1) + query_guard & ...
                       x_mm <= tpl.x_domain(2) - query_guard;

        legacy_fit_weight = build_legacy_fit_weight_local( ...
            t_seg, x_mm, v_seg, tpl, query_guard, S);

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
                % complete per-pass main pulse inside the query-safe
                % shrunken trusted domain.
                core_mask = finite_mask & inside_domain & inside_guard & ...
                            main_pulse_mask;
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
        bundle.legacy_fit_weight = [bundle.legacy_fit_weight; legacy_fit_weight(:)]; %#ok<AGROW>

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
bundle.legacy_fit_weight = [];

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
    % Match the 20241106 legacy NewFlow Step06 rule exactly:
    %   t_start = jilublade(row,1) - pulsePadSec
    %   t_end   = jilublade(row,2) + pulsePadSec
    % The default pulsePadSec is zero.
    t0 = start_s - S.pulse_pad_sec;
    t1 = end_s + S.pulse_pad_sec;
else
    % Original OPRCenterStd DynamicMap style:
    % fixed window centered at the target-blade passage time.
    t0 = arrival_s - S.pulse_window_sec;
    t1 = arrival_s + S.pulse_window_sec;
end

[t_seg, v_seg] = read_raw_time_window_local( ...
    case_dir, case_data, sid, t0, t1, cfg, S);
end


function [t, v] = read_raw_time_window_local(case_dir, case_data, sid, t0, t1, cfg, S)
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

    raw = load_raw_mat_file_local(file, S);

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


function raw = load_raw_mat_file_local(file, S)
persistent raw_cache raw_cache_keys
if nargin < 2
    S = struct();
end

cache_enable = getfield_default_local(S, 'raw_file_cache_enable', true);
cache_key = char(file);

if cache_enable
    if isempty(raw_cache)
        raw_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
        raw_cache_keys = {};
    end

    if isKey(raw_cache, cache_key)
        raw = raw_cache(cache_key);
        return;
    end
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

if cache_enable
    raw_cache(cache_key) = raw;
    raw_cache_keys{end + 1} = cache_key; %#ok<AGROW>

    max_cached_files = getfield_default_local(S, 'raw_file_cache_max_files', 12);
    while numel(raw_cache_keys) > max_cached_files
        old_key = raw_cache_keys{1};
        raw_cache_keys(1) = [];
        if isKey(raw_cache, old_key)
            remove(raw_cache, old_key);
        end
    end
end
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
    opr_times, t, cfg.opr_pulses_per_rev);

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


function theta_rot = map_time_to_rotor_phase_local(opr_times, sample_times, opr_pulses_per_rev)
sample_times = sample_times(:);
theta_rot = nan(size(sample_times));

if isempty(sample_times) || numel(opr_times) <= opr_pulses_per_rev
    return;
end

for i = 1:numel(sample_times)
    prev_idx = find(opr_times <= sample_times(i), 1, 'last');
    next_idx = prev_idx + opr_pulses_per_rev;

    if isempty(prev_idx) || next_idx > numel(opr_times)
        continue;
    end

    rev_dt = opr_times(next_idx) - opr_times(prev_idx);

    if rev_dt <= eps
        continue;
    end

    theta_rot(i) = 2 * pi * ...
        (sample_times(i) - opr_times(prev_idx)) / rev_dt;
end
end


function [template_v, template_dv_dx, template_weight, inside_domain] = ...
    evaluate_template_at_points_local(tpl, x, S)

x = x(:);

template_v = interp1( ...
    tpl.x_grid(:), ...
    tpl.v_grid(:), ...
    x, ...
    'pchip', ...
    NaN);

template_dv_dx = interp1( ...
    tpl.x_grid(:), ...
    tpl.dv_dx(:), ...
    x, ...
    'pchip', ...
    NaN);

if isfield(tpl, 'weight_grid') && ~isempty(tpl.weight_grid)
    template_weight = interp1( ...
        tpl.x_grid(:), ...
        tpl.weight_grid(:), ...
        x, ...
        'pchip', ...
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


function threshold = get_step05_sensor_threshold_local(tpl, S)
threshold = S.default_sensor_threshold;

if isstruct(tpl) && isfield(tpl, 'threshold') && ...
        isscalar(tpl.threshold) && isfinite(tpl.threshold)
    threshold = tpl.threshold;
end
end


function mask = build_dynamic_effective_mask_local( ...
    x, v, t, template_dv_dx, threshold, S)
% Match the 20241106 legacy Step06 effective-point selector.
x = x(:);
v = v(:);
t = t(:);
template_dv_dx = template_dv_dx(:);

finite = isfinite(x) & isfinite(v) & isfinite(t);

if ~strcmpi(S.dynamic_window_mode, 'legacy_row_bounds') || ...
        ~strcmpi(S.core_mask_mode, 'legacy_base_query_safe')
    finite = finite & isfinite(template_dv_dx);
end

if ~strcmpi(S.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end

g_tpl = zeros(size(x));
g_tpl(finite & isfinite(template_dv_dx)) = ...
    abs(template_dv_dx(finite & isfinite(template_dv_dx)));

g_tpl_max = max(g_tpl(finite), [], 'omitnan');

if ~isfinite(g_tpl_max) || g_tpl_max <= 0
    mask_tpl = finite;
else
    mask_tpl = g_tpl >= S.dynamic_template_gradient_min_ratio * g_tpl_max;
end

g_time = zeros(size(v));

if nnz(finite) >= 5 && range(t(finite)) > 0
    g_time(finite) = abs(gradient(v(finite), t(finite)));
end

g_time_max = max(g_time(finite), [], 'omitnan');

if ~isfinite(g_time_max) || g_time_max <= 0
    mask_time = false(size(v));
else
    mask_time = g_time >= S.dynamic_time_gradient_min_ratio * g_time_max;
end

if any(finite)
    peak_level = prctile(v(finite), ...
        min(max(S.dynamic_peak_quantile, 0), 100));
else
    peak_level = NaN;
end

if ~isfinite(peak_level)
    mask_peak = false(size(v));
else
    mask_peak = v >= max(threshold, peak_level);
end

mask = finite & (mask_tpl | mask_time | mask_peak) & ...
       v >= 0.5 * threshold;

if nnz(mask) < 8
    mask = finite & v >= threshold;
end

if nnz(mask) < 8
    mask = finite;
end
end


function mask = build_main_pulse_mask_local(v, threshold, S)
% Match the 20241106 legacy Step06 main-pulse selector.
v = v(:);
mask = true(size(v));

if ~strcmpi(S.pulse_mode, 'single') || numel(v) < 5
    return;
end

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
        peak_val = max(v(seg));

        if peak_val > best_peak
            best_peak = peak_val;
            best_seg = iseg;
        end
    end

    mask(idx(starts(best_seg)):idx(ends(best_seg))) = true;
end
end


function qg = compute_query_guard_local(x, v, tpl, S, base_mask, sensor_id)
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

if nargin < 6
    sensor_id = NaN;
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
qg = apply_sensor_query_guard_scale_local(qg, S, sensor_id);
qg = max(0, qg);
end


function qg = apply_sensor_query_guard_scale_local(qg, S, sensor_id)
if ~isfield(S, 'sensor_query_guard_scale_table') || ...
        isempty(S.sensor_query_guard_scale_table) || ~isfinite(sensor_id)
    return;
end

table_now = S.sensor_query_guard_scale_table;
idx = find(table_now(:,1) == sensor_id, 1);

if isempty(idx)
    return;
end

scale = table_now(idx, 2);

if isfinite(scale) && scale > 0
    qg = qg * scale;
end

if isfield(S, 'query_guard_min_mm') && isfinite(S.query_guard_min_mm)
    qg = max(qg, 0.5 * S.query_guard_min_mm);
end

if isfield(S, 'query_guard_max_mm') && isfinite(S.query_guard_max_mm)
    qg = min(qg, S.query_guard_max_mm);
end
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


function w_total = build_legacy_fit_weight_local(t, x, v, tpl, query_guard_mm, S)
% Match the legacy Step06 weighting route used for weighted RMSE ranking.
w_edge = build_edge_weight_local(t, v, S.weight_floor);
w_domain = build_domain_soft_weight_local( ...
    x, tpl.x_domain, S.domain_soft_margin_mm, S.weight_floor);
w_query = build_query_guard_soft_weight_local( ...
    x, [min(tpl.x_grid(:)), max(tpl.x_grid(:))], ...
    query_guard_mm, S.domain_soft_margin_mm, S.weight_floor);
w_gradient = build_template_gradient_weight_local( ...
    tpl, x, S.domain_selection_mode, S.weight_floor);
w_total = max(S.weight_floor, ...
    w_edge(:) .* w_domain(:) .* w_query(:) .* w_gradient(:));
if max(w_total) > 0
    w_total = max(S.weight_floor, w_total ./ max(w_total));
end
end


function w_edge = build_edge_weight_local(t, v, floor_w)
if numel(v) < 3 || range(t) <= 0
    w_edge = ones(size(v(:)));
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
    w_domain = ones(size(x(:)));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
ratio = min(max(dist_to_edge ./ margin_mm, 0), 1);
w_domain = floor_w + (1 - floor_w) .* ratio;
w_domain(~isfinite(w_domain)) = floor_w;
end


function w_query = build_query_guard_soft_weight_local(x, x_domain, query_guard_mm, margin_mm, floor_w)
if query_guard_mm <= 0 || margin_mm <= 0
    w_query = ones(size(x(:)));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
soft_start = max(query_guard_mm - margin_mm, 0);
ratio = min(max((dist_to_edge - soft_start) ./ max(margin_mm, eps), 0), 1);
w_query = floor_w + (1 - floor_w) .* ratio;
w_query(~isfinite(w_query)) = floor_w;
end


function w_gradient = build_template_gradient_weight_local(tpl, x, domain_selection_mode, floor_w)
if ~strcmpi(domain_selection_mode, 'soft') || ~isfield(tpl, 'dv_dx') || isempty(tpl.dv_dx)
    w_gradient = ones(size(x(:)));
    return;
end
g = abs(interp1(tpl.x_grid(:), tpl.dv_dx(:), x(:), 'linear', 0));
g_max = max(g, [], 'omitnan');
if ~isfinite(g_max) || g_max <= 0
    w_gradient = ones(size(x(:)));
    return;
end
w_gradient = floor_w + (1 - floor_w) .* g ./ g_max;
w_gradient(~isfinite(w_gradient)) = floor_w;
end


function eo_candidates = build_eo_candidates_local(rot_freq_mean_hz, freq_search_hz, eo_pad)
if ~isfinite(rot_freq_mean_hz) || rot_freq_mean_hz <= 0
    error('Invalid rot_freq_mean_hz for EO candidate generation.');
end

freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min = max(1, ceil(freq_lo / rot_freq_mean_hz) - max(0, floor(eo_pad)));
eo_max = max(eo_min, floor(freq_hi / rot_freq_mean_hz) + max(0, floor(eo_pad)));
eo_all = eo_min:eo_max;
freq_all = eo_all .* rot_freq_mean_hz;
eo_candidates = eo_all(freq_all >= freq_lo & freq_all <= freq_hi);
if isempty(eo_candidates)
    eo_center = max(1, round(mean([freq_lo, freq_hi]) / rot_freq_mean_hz));
    eo_candidates = max(1, eo_center - eo_pad):max(1, eo_center + eo_pad);
end
eo_candidates = unique(round(eo_candidates(:).'));
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
    'legacy_fit_weight', ...
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

if isfield(sub, 'legacy_fit_weight') && numel(sub.legacy_fit_weight) == sub.point_count
    sub.fit_weight = sub.legacy_fit_weight(:);
else
    sub.fit_weight = sub.template_weight(:);
end
sub.fit_weight(~isfinite(sub.fit_weight)) = S.weight_floor;
sub.fit_weight = max(sub.fit_weight, S.weight_floor);
if max(sub.fit_weight) > 0
    sub.fit_weight = max(S.weight_floor, sub.fit_weight ./ max(sub.fit_weight));
end
end
function seed_table = solve_vp_seed_eo_scan_local(bundle, eo_candidates, S)
% First-order VP seed scan using non-parametric template gradient.
%
% Approximation:
%   V - T(x) ~= -T'(x) * [dx_c + a*sin(EO*theta) + b*cos(EO*theta)]
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
    'linear_vp_rmse', inf, ...
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

grad_abs = abs(Tp);
grad_ref = prctile(grad_abs(isfinite(grad_abs)), S.vp_gradient_reference_quantile);
if ~isfinite(grad_ref) || grad_ref <= 0
    grad_ref = max(grad_abs(isfinite(grad_abs)), [], 'omitnan');
end
if ~isfinite(grad_ref) || grad_ref <= 0
    return;
end

static_eta = resolve_sensor_eta_prior_local(S);
eta_vec0 = zeros(size(sensor_id));
for is = 1:n_sensors
    eta_vec0(sensor_id == sensors(is)) = static_eta(is);
end

valid = isfinite(v) & ...
        isfinite(theta) & ...
        isfinite(T0) & ...
        isfinite(Tp) & ...
        isfinite(w) & ...
        isfinite(eta_vec0) & ...
        grad_abs >= max(S.vp_gradient_min_ratio * grad_ref, eps);

if nnz(valid) < S.min_fit_points
    return;
end

v = v(valid);
theta = theta(valid);
T0 = T0(valid);
Tp = Tp(valid);
w = w(valid);
sensor_id = sensor_id(valid);
eta_vec0 = eta_vec0(valid);

q_obs = -(v - T0) ./ Tp - eta_vec0;

wq = w .* min(abs(Tp) ./ max(grad_ref, eps), 1) .^ ...
    max(0, S.vp_gradient_weight_power);
wq(~isfinite(wq)) = S.weight_floor;
wq = max(wq, S.weight_floor);

for iEO = 1:numel(eo_candidates)
    EO = eo_candidates(iEO);

    s = sin(EO .* theta);
    c = cos(EO .* theta);

    basis = [ones(size(q_obs)), s, c];
    y = q_obs;
    good = all(isfinite(basis), 2) & isfinite(y) & isfinite(wq);

    if nnz(good) < max(5, size(basis, 2) + 1)
        continue;
    end

    Xg = basis(good, :);
    yg = y(good);
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
    A = hypot(a, b);
    phi = atan2(b, a);
    lin_res = y - basis * beta;

    A = min(max(abs(A), 0), S.amplitude_limit_mm);
    dx_c = max(min(dx_c, S.dx_c_limit_mm), -S.dx_c_limit_mm);

    [wrmse, prmse] = evaluate_full_template_model_local( ...
        bundle, EO, A, phi, dx_c, eta, S);

    seed_table(iEO).EO = EO;
    seed_table(iEO).A = A;
    seed_table(iEO).phi = wrap_to_pi_local(phi);
    seed_table(iEO).dx_c = dx_c;
    seed_table(iEO).sensor_eta = eta;
    seed_table(iEO).weighted_voltage_rmse = wrmse;
    seed_table(iEO).plain_voltage_rmse = prmse;
    seed_table(iEO).linear_vp_rmse = sqrt(nansum(wq(good) .* lin_res(good).^2) / ...
        max(nansum(wq(good)), eps));
    seed_table(iEO).point_count = bundle.point_count;
    seed_table(iEO).rank_score = wrmse;
end
end


function final_eos = select_top_eos_local(seed_table, top_k)
scores = [[seed_table.weighted_voltage_rmse].', ...
          [seed_table.linear_vp_rmse].', ...
          [seed_table.plain_voltage_rmse].', ...
          [seed_table.EO].'];
eos = [seed_table.EO];

valid = all(isfinite(scores), 2).' & isfinite(eos);

if ~any(valid)
    final_eos = [];
    return;
end

scores_valid = scores(valid, :);
eos_valid = eos(valid);

[~, order] = sortrows(scores_valid, [1, 2, 3, 4]);

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
gap_thr = max(S.vp_gap_ratio_keep, 0);

best_score = T.rank_score(1);
gap_ratio = (T.rank_score - best_score) ./ max(abs(best_score), eps);
keep = false(height(T), 1);

keep(1:min(top_min, height(T))) = true;
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
    seed.sensor_eta = zeros(1, numel(S.analysis_sensors));

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

if isfield(S, 'sensor_eta_mode') && strcmpi(S.sensor_eta_mode, 'fixed_static')
    seed.sensor_eta = resolve_sensor_eta_prior_local(S);
elseif ~isfield(S, 'fit_sensor_eta') || ~S.fit_sensor_eta
    seed.sensor_eta(:) = 0;
end

seed.sensor_eta = project_sensor_eta_local(seed.sensor_eta, S);
end


function p = pack_params_local(A, phi, dx_c, eta, S)
eta = eta(:).';

if numel(eta) ~= numel(S.analysis_sensors)
    eta = zeros(1, numel(S.analysis_sensors));
end

eta(1) = 0;

% Parameter vector:
%   [A, phi, dx_c] for the formal NoEta route.
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
elseif isfield(S, 'sensor_eta_mode') && strcmpi(S.sensor_eta_mode, 'fixed_static')
    eta = resolve_sensor_eta_prior_local(S);
end

% Soft projection for reported parameters.
% Objective also uses penalties on the unprojected values.
A = min(max(abs(A), 0), S.amplitude_limit_mm);
phi = wrap_to_pi_local(phi);
dx_c = max(min(dx_c, S.dx_c_limit_mm), -S.dx_c_limit_mm);

eta = project_sensor_eta_local(eta, S);
end


function eta = resolve_sensor_eta_prior_local(S)
eta = zeros(1, numel(S.analysis_sensors));
if isfield(S, 'sensor_eta_prior_mm') && ...
        isnumeric(S.sensor_eta_prior_mm) && ...
        numel(S.sensor_eta_prior_mm) == numel(S.analysis_sensors)
    eta = S.sensor_eta_prior_mm(:).';
end
eta(~isfinite(eta)) = 0;
eta(1) = 0;
end


function eta = project_sensor_eta_local(eta, S)
eta = eta(:).';
if numel(eta) ~= numel(S.analysis_sensors)
    eta = zeros(1, numel(S.analysis_sensors));
end
eta(~isfinite(eta)) = 0;
eta(1) = 0;

if isfield(S, 'sensor_eta_mode') && strcmpi(S.sensor_eta_mode, 'fixed_static')
    eta = resolve_sensor_eta_prior_local(S);
    return;
end

eta = max(min(eta, S.sensor_eta_limit_mm), ...
          -S.sensor_eta_limit_mm);
eta(1) = 0;
end


function value = direct_template_objective_local(p, EO, bundle, S)
[A, phi, dx_c, eta] = unpack_params_local(p, S);

[wrmse, ~, ~, coverage, objective_score] = evaluate_full_template_model_local( ...
    bundle, EO, A, phi, dx_c, eta, S);

if ~isfinite(objective_score)
    objective_score = 1e12;
end

value = objective_score;
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
            'pchip', ...
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

    residual_obj = nansum(w(valid) .* res(valid).^2);
    overshoot_obj = S.overshoot_penalty_weight * ...
        nansum(w(valid) .* overshoot(valid).^2);
    eta_reg_obj = bundle.point_count * ...
        (S.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta(:).^2))) ^ 2;
    missing_obj = 1e6 * nnz(~valid);

    rmse_objective = residual_obj + overshoot_obj + eta_reg_obj + missing_obj;
    if isfield(S, 'objective_mode') && strcmpi(S.objective_mode, 'corr')
        corr_score = compute_weighted_sensor_correlation_score_local( ...
            bundle, Vpred, valid, w);
        objective_score = max(0, 1 - corr_score) + ...
            (overshoot_obj + eta_reg_obj + missing_obj) / max(bundle.point_count, 1);
    else
        corr_score = NaN;
        objective_score = rmse_objective;
    end
    wrmse = sqrt(rmse_objective / max(bundle.point_count, 1));

    prmse = sqrt(nanmean(res(valid).^2));
end

coverage = struct();
coverage.clamp_fraction = mean(overshoot > 0, 'omitnan');
coverage.overshoot_rms_mm = sqrt(nanmean(overshoot(:).^2));
coverage.max_overshoot_mm = max(overshoot(:), [], 'omitnan');
coverage.correlation_score = corr_score;
if exist('objective_score', 'var')
    coverage.objective_score = objective_score;
    coverage.overshoot_penalty = overshoot_obj;
    coverage.sensor_eta_reg_penalty = eta_reg_obj;
else
    coverage.objective_score = inf;
coverage.overshoot_penalty = NaN;
    coverage.sensor_eta_reg_penalty = NaN;
end
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
    'valid_sensor_count', NaN, ...
    'valid_lap_count', NaN, ...
    'expected_valid_pass_count', NaN, ...
    'min_valid_pass_count', NaN, ...
    'min_valid_sensor_count', NaN, ...
    'min_valid_lap_count', NaN, ...
    'valid_pass_fraction', NaN, ...
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

    if isfield(bundle, 'Coverage') && isstruct(bundle.Coverage)
        row = attach_coverage_to_trend_row_local(row, bundle.Coverage);
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
S.sensor_eta_reg_weight_v_per_mm = 0.02;
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
rmse_vec = [];
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
    else
        continue;
    end

    if ismember('objective_score', C.Properties.VariableNames)
        objective = C.objective_score;
    else
        objective = score;
    end

    ok = isfinite(C.EO) & isfinite(score);
    if ~any(ok)
        continue;
    end

    eo_vec = [eo_vec; C.EO(ok)]; %#ok<AGROW>
    rmse_vec = [rmse_vec; score(ok)]; %#ok<AGROW>
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
    rows(i).mean_weighted_rmse = mean(rmse_vec(m), 'omitnan');
    rows(i).median_weighted_rmse = median(rmse_vec(m), 'omitnan');
    rows(i).min_weighted_rmse = min(rmse_vec(m), [], 'omitnan');
    rows(i).mean_objective = mean(objective_vec(m), 'omitnan');
    rows(i).median_objective = median(objective_vec(m), 'omitnan');
    rows(i).std_objective = std(objective_vec(m), 0, 'omitnan');
    rows(i).window_fraction = rows(i).window_count / max(valid_window_count, 1);
    rows(i).low_coverage_flag = rows(i).window_fraction < min_window_fraction;
    rows(i).joint_score = rows(i).mean_weighted_rmse .* ...
        (1 + missing_window_penalty .* (1 - rows(i).window_fraction));
end

T = struct2table(rows);
T = sortrows(T, { ...
    'low_coverage_flag', ...
    'joint_score', ...
    'median_weighted_rmse', ...
    'EO'}, ...
    {'ascend', 'ascend', 'ascend', 'ascend'});

Joint.status = 'ok';
Joint.dominant_eo = T.EO(1);
Joint.best_EO = T.EO(1);
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
    Summary.dominant_eo = Summary.joint_best_eo;
    Summary.trend_dominant_eo = NaN;
    Summary.reported_eo = Summary.joint_best_eo;
    Summary.reported_freq_hz = compute_reported_frequency_local( ...
        Trend, isfinite(Trend.rot_freq_mean_hz), Summary.reported_eo);
    Summary.reported_eo_source = 'joint_no_valid_trend';
    Summary.mean_freq_hz = Summary.reported_freq_hz;
    Summary.median_freq_hz = Summary.reported_freq_hz;
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

if isfinite(Summary.joint_best_eo)
    Summary.reported_eo = Summary.joint_best_eo;
    Summary.reported_eo_source = 'joint_eo_summary';
else
    Summary.reported_eo = Summary.trend_dominant_eo;
    Summary.reported_eo_source = 'trend_mode';
end

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

if any(ok)
    stem(ax, T.window_id(ok), ones(nnz(ok), 1), ...
        'filled', 'DisplayName', 'ok');
end

if any(~ok)
    stem(ax, T.window_id(~ok), zeros(nnz(~ok), 1), ...
        'filled', 'DisplayName', 'failed');
end

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


function value = get_cfg_text_local(cfg, field, default_value)
if isstruct(cfg) && isfield(cfg, field) && ~isempty(cfg.(field))
    value = char(cfg.(field));
else
    value = char(default_value);
end
end


function value = get_cfg_scalar_local(cfg, field, default_value)
if isstruct(cfg) && isfield(cfg, field) && ...
        isscalar(cfg.(field)) && isfinite(cfg.(field))
    value = cfg.(field);
else
    value = default_value;
end
end


function opr_times = load_dynamic_opr_times_local(step02_case_dir, cfg)
timing_source = get_cfg_text_local(cfg, 'step05_timing_source', 'step02');

if strcmpi(timing_source, 'legacy_newflow') && ...
        isfield(cfg, 'step05_legacy_pulse_dir') && ...
        ~isempty(cfg.step05_legacy_pulse_dir)
    legacy_file = fullfile(cfg.step05_legacy_pulse_dir, 'jiluOPR.mat');

    if isfile(legacy_file)
        jiluopr_file = legacy_file;
    else
        warning('Legacy timing source requested, but jiluOPR is missing: %s. Falling back to Step02.', ...
            legacy_file);
        jiluopr_file = fullfile(step02_case_dir, 'jiluOPR.mat');
    end
else
    jiluopr_file = fullfile(step02_case_dir, 'jiluOPR.mat');
end

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
