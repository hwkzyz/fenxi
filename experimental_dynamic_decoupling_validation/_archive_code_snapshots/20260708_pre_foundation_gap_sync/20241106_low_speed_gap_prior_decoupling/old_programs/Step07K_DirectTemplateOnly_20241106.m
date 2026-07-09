%% Step07K: 20241106 direct template identification
% This script now follows the same direct-template workflow used in
% 20250527 / 20251222 Step03:
%   1) load low-speed template + continuous high-speed DynamicMap windows
%   2) select effective waveform points inside each window
%   3) first-order template VP screens EO candidates
%   4) bounded full-waveform refinement solves [A, phi, dx_c, eta_s]
%   5) phase-safe expansion uses previous-window or current-window reference
%
% Only the data source is changed to the 20241106 Step06K DynamicMap.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
outDir = fullfile(routeDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07k_20241106_direct_only');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Parameters to tune
P = struct();
P.targetBlade = 1;
P.analysisSensors = [5 7];
P.selectedWindowId = 1;      % detail plot / exported summary window
P.debugMaxWindows = inf;     % set to 1..N to limit processing

P.directDynamicMapFile = '';
P.templateFile = '';

% Keep these aligned with 20250527 / 20251222 Step03 defaults.
P.freqSearchHz = [100 1000];
P.eoPad = 2;
P.topKEO = 3;
P.amplitudeLimitMm = 0.50;
P.dxLimitMm = 0.35;
P.sensorEtaLimitMm = 0;
P.sensorEtaRegWeightVPerMm = 0;

P.dynamicEffectiveMode = 'gradient';
P.dynamicTemplateGradientMinRatio = 0.08;
P.dynamicTimeGradientMinRatio = 0.15;
P.dynamicPeakQuantile = 85;
P.pulseMode = 'single';

P.domainSelectionMode = 'hard';
P.domainSoftMarginMm = 0;
P.queryGuardMM = [];
P.queryGuardMode = 'adaptive';
P.queryGuardQuantile = 95;
P.queryGuardSafetyMm = 0.05;
P.queryGuardMinMm = 0.12;
P.queryGuardMaxMm = [];
P.overshootPenaltyWeight = 100;

P.phaseSafeExpansion = true;
P.phaseSafeMarginMm = 0.03;
P.phaseSafeReferenceMode = 'prev_window';
P.phaseSafeFallbackMode = 'linear_vp';
P.phaseSafeRefreshEvery = 5;
P.phaseSafePrevMaxClampFraction = 1e-6;
P.phaseSafePrevMinFinalGapRatio = 0.005;

P.allEoWarmupWindows = 0;
P.vpGapRatioFallback = -inf;
P.vpLinearGapRatioFallback = -inf;
P.defaultSensorThresholdMv = 0.5;

P.targetBlade = round(parse_scalar_env_local('STEP07K_TARGET_BLADE', P.targetBlade));
P.analysisSensors = parse_int_env_local('STEP07K_ANALYSIS_SENSORS', P.analysisSensors);
P.selectedWindowId = round(parse_scalar_env_local('STEP07K_WINDOW_ID', P.selectedWindowId));
P.debugMaxWindows = parse_scalar_env_local('STEP07K_DEBUG_MAX_WINDOWS', P.debugMaxWindows);

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
if isempty(P.directDynamicMapFile)
    P.directDynamicMapFile = fullfile(outDir, sprintf( ...
        'Step06K_DirectDynamicMap_20241106_B%d_%s.mat', P.targetBlade, sensorTag));
end
if isempty(P.templateFile)
    P.templateFile = fullfile(outDir, sprintf( ...
        'Step06I_Generated_LowSpeed_Template_20241106_B%d_%s.mat', ...
        P.targetBlade, sensorTag));
end

if ~isfile(P.directDynamicMapFile)
    error('Direct DynamicMap file not found. Run Step06K first: %s', P.directDynamicMapFile);
end
if ~isfile(P.templateFile)
    error('Low-speed template file not found: %s', P.templateFile);
end

fprintf('\n=== Step07K 20241106 direct-template-only ===\n');
fprintf('Direct DynamicMap: %s\n', P.directDynamicMapFile);
fprintf('Template:          %s\n', P.templateFile);
fprintf('Blade B%d, sensors %s\n', P.targetBlade, mat2str(P.analysisSensors));
fprintf('This run copies the 20250527/20251222 Step03 direct logic.\n');
fprintf('Bounded parameters: A <= %.3f mm, |dx_c| <= %.3f mm, |eta_s| <= %.3f mm\n', ...
    P.amplitudeLimitMm, P.dxLimitMm, P.sensorEtaLimitMm);

S = load(P.templateFile, 'Template');
D = load(P.directDynamicMapFile, 'DynamicMap');

Template = clean_template_local(S.Template, P.analysisSensors);
DynamicMap = filter_dynamic_map_sensors_local(D.DynamicMap, P.analysisSensors, sensorTag);
coordinateCheck = check_template_dynamic_xcenter_local(Template, DynamicMap, P.analysisSensors, 1e-9);
source = make_source_info_local(DynamicMap, P);
cfg = build_cfg_local(DynamicMap, P);
method = build_method_local(P, cfg);

fprintf('Raw low-speed folder:  %s\n', source.lowDir);
fprintf('Raw high-speed folder: %s\n', source.highDir);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s)\n', ...
    coordinateCheck.max_abs_delta_mm, coordinateCheck.status);
fprintf('Selected waveform mode: %s / %s / phase-safe=%d\n', ...
    method.pulse_selection_mode, method.dynamic_effective_mode, method.phase_safe_expansion_enabled);

numWindows = numel(DynamicMap.Window);
if isfinite(cfg.debug_max_windows)
    numWindows = min(numWindows, max(1, floor(cfg.debug_max_windows)));
end
if P.selectedWindowId < 1 || P.selectedWindowId > numWindows
    error('STEP07K_WINDOW_ID=%d is outside processed windows 1-%d.', P.selectedWindowId, numWindows);
end

trendRows = repmat(struct( ...
    'window_id', NaN, ...
    'lap_start', NaN, ...
    'lap_end', NaN, ...
    'window_center_time', NaN, ...
    'rot_freq_mean_hz', NaN, ...
    'rot_rpm_mean', NaN, ...
    'A_id', NaN, ...
    'EO_id', NaN, ...
    'fn_id', NaN, ...
    'phi_id_wrapped', NaN, ...
    'dx_c_id', NaN, ...
    'sensor_eta_max_abs_mm', NaN, ...
    'weighted_voltage_rmse', NaN, ...
    'plain_voltage_rmse', NaN, ...
    'valid_segment_count', NaN, ...
    'point_count', NaN, ...
    'core_point_count', NaN, ...
    'added_point_count', NaN, ...
    'phase_safe_reference_source', '', ...
    'phase_safe_reference_reason', ''), numWindows, 1);

windowResults = repmat(struct( ...
    'window_id', NaN, ...
    'lap_range', [], ...
    'time_window', [], ...
    'bundle', [], ...
    'seed_table', struct([]), ...
    'CandidateTable', table(), ...
    'Result', []), numWindows, 1);

prevPhaseResult = [];
for iw = 1:numWindows
    [trendRows(iw), windowResults(iw)] = process_window_local( ...
        DynamicMap.Window(iw), Template, cfg, method, iw, prevPhaseResult);
    prevPhaseResult = windowResults(iw).Result;
    fprintf('Window %2d / %2d, laps %s, EO %2d, A %.5f mm, dx_c %.5f mm, WRMSE %.3f mV\n', ...
        iw, numWindows, mat2str(windowResults(iw).lap_range), ...
        windowResults(iw).Result.EO_id, windowResults(iw).Result.A_id, ...
        windowResults(iw).Result.dx_c_id, windowResults(iw).Result.weighted_voltage_rmse);
end

trendTable = struct2table(trendRows);
selected = windowResults(P.selectedWindowId);
selectedSummary = build_selected_summary_local(selected, source, P);
selectedEOScan = struct2table(selected.seed_table);
selectedEOScan = sortrows(selectedEOScan, {'weighted_voltage_rmse','linear_vp_rmse','EO'}, ...
    {'ascend','ascend','ascend'});

result = struct();
result.dataset = '20241106';
result.method = 'direct_template_step03_equivalent_20241106';
result.description = ['Copied direct-template workflow from 20250527/20251222 Step03; ' ...
    'only the data source is changed to 20241106 Step06K DynamicMap.'];
result.P = P;
result.cfg = cfg;
result.methodConfig = method;
result.source = source;
result.coordinateCheck = coordinateCheck;
result.Template = Template;
result.TrendTable = trendTable;
result.WindowResults = windowResults;
result.SelectedWindow = selected;
result.SelectedWindowSummary = selectedSummary;

base = sprintf('Step07K_DirectTemplateOnly_20241106_B%d_%s', P.targetBlade, sensorTag);
matFile = fullfile(outDir, [base, '.mat']);
trendCsv = fullfile(outDir, [base, '_WindowTrend.csv']);
seedCsv = fullfile(outDir, [base, '_EOScan.csv']);
summaryCsv = fullfile(outDir, [base, '_Summary.csv']);
figFile = fullfile(figDir, [base, sprintf('_W%d.png', P.selectedWindowId)]);

save(matFile, 'result', 'P', 'cfg', 'method', 'source', 'Template', ...
    'trendTable', 'windowResults', 'selectedSummary', '-v7.3');
writetable(trendTable, trendCsv);
writetable(selectedEOScan, seedCsv);
writetable(selectedSummary, summaryCsv);
plot_selected_window_local(selected, Template, figFile);

fprintf('\nSelected window %d summary:\n', P.selectedWindowId);
disp(selectedSummary);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    matFile, trendCsv, seedCsv, summaryCsv, figFile);

%% Local functions
function values = parse_int_env_local(name, defaultValues)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValues;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    values = defaultValues;
end
end

function value = parse_scalar_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
tmp = str2double(raw);
if isfinite(tmp)
    value = tmp;
else
    value = defaultValue;
end
end

function cfg = build_cfg_local(DynamicMap, P)
cfg = struct();
cfg.target_blades = P.targetBlade;
cfg.analysis_sensors = P.analysisSensors(:).';
cfg.analysis_start_time = get_source_scalar_local(DynamicMap.SourceSettings, 'analysisStartTimeSec', NaN);
cfg.target_laps = get_source_scalar_local(DynamicMap.SourceSettings, 'targetLaps', NaN);
cfg.analysis_win_size = get_source_scalar_local(DynamicMap.SourceSettings, 'windowLaps', NaN);
cfg.sliding_step = get_source_scalar_local(DynamicMap.SourceSettings, 'slidingStepLaps', NaN);
cfg.freq_search_hz = P.freqSearchHz;
cfg.eo_pad = P.eoPad;
cfg.amplitude_limit_mm = P.amplitudeLimitMm;
cfg.dx_c_limit_mm = P.dxLimitMm;
cfg.sensor_eta_limit_mm = P.sensorEtaLimitMm;
cfg.sensor_eta_reg_weight_v_per_mm = P.sensorEtaRegWeightVPerMm;
cfg.debug_max_windows = P.debugMaxWindows;
cfg.vp_top_k_eo = P.topKEO;
cfg.vp_all_eo_warmup_windows = P.allEoWarmupWindows;
cfg.vp_gap_ratio_fallback = P.vpGapRatioFallback;
cfg.vp_linear_gap_ratio_fallback = P.vpLinearGapRatioFallback;
end

function value = get_source_scalar_local(S, fieldName, fallback)
value = fallback;
if isstruct(S) && isfield(S, fieldName)
    tmp = S.(fieldName);
    if isnumeric(tmp) && isscalar(tmp) && isfinite(tmp)
        value = tmp;
    end
end
end

function method = build_method_local(P, cfg)
method = struct();
method.name = 'template_only_main_vp_top3_synchronous_waveform';
method.selection_rule = 'first_order_template_vp_top3_then_synchronous_waveform_rmse';
method.use_reference_eo_constraint = false;
method.use_vp = true;
method.vp_top_k_eo = cfg.vp_top_k_eo;
method.refine_params = {'A', 'phi', 'dx_c'};
method.template_forward = 'interp1_low_speed_template';
method.vp_seed_model = 'V_minus_Tx_equals_minus_Tprime_times_u';
method.final_waveform_objective = 'fixed_eo_direct_template_voltage_residual_without_sensor_affine_projection';
method.store_y_obs = false;
method.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
method.sensor_eta_reg_weight_v_per_mm = cfg.sensor_eta_reg_weight_v_per_mm;
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.query_guard_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm + method.coverage_safety_margin_mm;
if ~isempty(P.queryGuardMM)
    method.query_guard_mm = P.queryGuardMM;
end
method.query_guard_mode = lower(strtrim(P.queryGuardMode));
method.query_guard_quantile = P.queryGuardQuantile;
method.query_guard_safety_mm = P.queryGuardSafetyMm;
method.query_guard_min_mm = P.queryGuardMinMm;
method.query_guard_max_mm = method.query_guard_mm;
if ~isempty(P.queryGuardMaxMm)
    method.query_guard_max_mm = P.queryGuardMaxMm;
end
method.phase_safe_expansion_enabled = logical(P.phaseSafeExpansion);
method.phase_safe_margin_mm = max(0, P.phaseSafeMarginMm);
method.phase_safe_reference_mode = lower(strtrim(P.phaseSafeReferenceMode));
method.phase_safe_fallback_mode = lower(strtrim(P.phaseSafeFallbackMode));
method.phase_safe_refresh_every = max(0, floor(P.phaseSafeRefreshEvery));
method.phase_safe_prev_max_clamp_fraction = P.phaseSafePrevMaxClampFraction;
method.phase_safe_prev_min_final_gap_ratio = P.phaseSafePrevMinFinalGapRatio;
method.phase_safe_u_limit_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm;
method.domain_selection_mode = lower(strtrim(P.domainSelectionMode));
method.domain_soft_margin_mm = P.domainSoftMarginMm;
method.overshoot_penalty_weight = P.overshootPenaltyWeight;
method.coverage_mode = 'prefer_points_safe_for_all_bounded_u';
method.extrapolation_mode = 'clamp_to_template_edge';
method.default_sensor_threshold = P.defaultSensorThresholdMv;
method.pulse_selection_mode = lower(strtrim(P.pulseMode));
method.dynamic_effective_mode = lower(strtrim(P.dynamicEffectiveMode));
method.dynamic_template_gradient_min_ratio = P.dynamicTemplateGradientMinRatio;
method.dynamic_time_gradient_min_ratio = P.dynamicTimeGradientMinRatio;
method.dynamic_peak_quantile = P.dynamicPeakQuantile;
end

function source = make_source_info_local(DynamicMap, P)
source = struct();
source.lowDir = get_string_field_local(DynamicMap, 'lowDir', '<unknown low-speed raw folder>');
source.highDir = get_string_field_local(DynamicMap, 'highDir', '<unknown high-speed raw folder>');
source.templateFile = P.templateFile;
source.directDynamicMapFile = P.directDynamicMapFile;
source.lowDataRole = '20241106_2\\900 low-speed pulses used to build the template';
source.highDataRole = '20241106_2\\3000_3150 high-speed continuous pulses used to build DynamicMap windows';
end

function value = get_string_field_local(S, fieldName, fallback)
value = fallback;
if isstruct(S) && isfield(S, fieldName)
    tmp = S.(fieldName);
    if isstring(tmp) || ischar(tmp)
        value = char(tmp);
    end
end
end

function Template = clean_template_local(Template, analysisSensors)
available = [Template.Sensor.sensor_id];
idx = zeros(size(analysisSensors));
for is = 1:numel(analysisSensors)
    pos = find(available == analysisSensors(is), 1, 'first');
    if isempty(pos)
        error('Template does not contain sensor %d.', analysisSensors(is));
    end
    idx(is) = pos;
end
Template.Sensor = Template.Sensor(idx);
for is = 1:numel(Template.Sensor)
    Tsen = Template.Sensor(is);
    x = Tsen.x_grid(:);
    v = Tsen.v_grid(:);
    finite = isfinite(x) & isfinite(v);
    x = x(finite);
    v = v(finite);
    [x, order] = sort(x);
    v = v(order);
    [x, keep] = unique(x, 'stable');
    v = v(keep);
    if numel(x) < 20
        error('Template sensor %d has too few finite points.', Tsen.sensor_id);
    end
    Template.Sensor(is).x_grid = x;
    Template.Sensor(is).v_grid = v;
    Template.Sensor(is).dv_dx = gradient(v, x);
    Template.Sensor(is).x_domain = [min(x), max(x)];
    if ~isfield(Template.Sensor(is), 'bin_weight') || isempty(Template.Sensor(is).bin_weight)
        Template.Sensor(is).bin_weight = ones(size(x));
    else
        bw = Template.Sensor(is).bin_weight(:);
        bw = bw(finite);
        bw = bw(order);
        bw = bw(keep);
        if numel(bw) ~= numel(x)
            bw = ones(size(x));
        end
        Template.Sensor(is).bin_weight = bw;
    end
    if ~isfield(Template.Sensor(is), 'threshold') || isempty(Template.Sensor(is).threshold) || ...
            ~isfinite(Template.Sensor(is).threshold)
        Template.Sensor(is).threshold = 0.5;
    end
    if ~isfield(Template.Sensor(is), 'xc') || ~isfinite(Template.Sensor(is).xc)
        Template.Sensor(is).xc = 0;
    end
end
Template.SensorIDs = analysisSensors(:).';
Template.SensorTag = ['S', sprintf('%d', analysisSensors)];
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, analysisSensors, sensorTag)
for iw = 1:numel(DynamicMap.Window)
    available = [DynamicMap.Window(iw).Sensor.sensor_id];
    idx = zeros(size(analysisSensors));
    for is = 1:numel(analysisSensors)
        pos = find(available == analysisSensors(is), 1, 'first');
        if isempty(pos)
            error('DynamicMap window %d does not contain sensor %d.', iw, analysisSensors(is));
        end
        idx(is) = pos;
    end
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(idx);
end
DynamicMap.SensorIDs = analysisSensors(:).';
DynamicMap.SensorTag = sensorTag;
end

function check = check_template_dynamic_xcenter_local(Template, DynamicMap, analysisSensors, toleranceMm)
rows = repmat(struct('sensor_id', NaN, 'template_xc_mm', NaN, ...
    'dynamic_xc_mm', NaN, 'delta_mm', NaN), numel(analysisSensors), 1);
for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    rows(is).sensor_id = sid;
    rows(is).template_xc_mm = Tpl.xc;
    rows(is).dynamic_xc_mm = infer_dynamic_xcenter_local(DynamicMap, sid);
    rows(is).delta_mm = rows(is).template_xc_mm - rows(is).dynamic_xc_mm;
end
delta = [rows.delta_mm];
check = struct();
check.table = struct2table(rows);
check.tolerance_mm = toleranceMm;
check.max_abs_delta_mm = max(abs(delta), [], 'omitnan');
check.is_consistent = isfinite(check.max_abs_delta_mm) && check.max_abs_delta_mm <= toleranceMm;
if check.is_consistent
    check.status = 'consistent';
else
    check.status = 'mismatch';
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
    tmp = S(idx).x_abs(:) - S(idx).x_rel(:);
    xc = median(tmp(isfinite(tmp)), 'omitnan');
    return;
end
end

function [trendRow, windowResult] = process_window_local(Wmap, Template, cfg, method, wIdx, prevPhaseResult)
if nargin < 6
    prevPhaseResult = [];
end
bundle = build_template_observation_bundle_local(Wmap, Template, cfg.analysis_sensors, method);
eoCandidates = build_eo_candidates_local(bundle.rot_freq_mean_hz, cfg.freq_search_hz, cfg.eo_pad);
coreBundle = bundle;
phaseRefInfo = empty_phase_safe_reference_info_local();

if method.phase_safe_expansion_enabled
    if strcmpi(method.phase_safe_reference_mode, 'template_inverse')
        phaseRef = build_template_inverse_phase_safe_reference_local();
        phaseRefInfo = make_phase_safe_reference_info_local( ...
            phaseRef.reference_source, 'configured_template_inverse', NaN, false);
        expandedBundle = build_template_observation_bundle_local( ...
            Wmap, Template, cfg.analysis_sensors, method, phaseRef);
        if expandedBundle.point_count > coreBundle.point_count
            bundle = expandedBundle;
        end
        seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
        [finalEOCandidates, selectionInfo] = select_vp_adaptive_eo_local( ...
            seedTable, eoCandidates, method, cfg, wIdx);
        result = refine_template_waveform_fit_local(bundle, seedTable, finalEOCandidates, cfg);
        result.CorePass = summarize_phase_safe_reference_local(phaseRef, coreBundle, method);
        result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
            coreBundle, bundle, phaseRef, method);
    else
        if strcmpi(method.phase_safe_reference_mode, 'prev_window')
            [phaseRef, phaseRefInfo] = build_prev_window_phase_safe_reference_local( ...
                prevPhaseResult, wIdx, method, cfg);
            if isempty(phaseRef)
                [phaseRef, phaseRefInfo] = build_current_window_phase_safe_reference_local( ...
                    bundle, eoCandidates, cfg, method, wIdx, phaseRefInfo.reason);
            end
        else
            [phaseRef, phaseRefInfo] = build_current_window_phase_safe_reference_local( ...
                bundle, eoCandidates, cfg, method, wIdx, 'configured_current_window_reference');
        end
        expandedBundle = build_template_observation_bundle_local( ...
            Wmap, Template, cfg.analysis_sensors, method, phaseRef);
        if expandedBundle.point_count > coreBundle.point_count
            expandedSeedTable = solve_template_seed_eo_scan_local(expandedBundle, eoCandidates);
            [expandedEOCandidates, expandedSelectionInfo] = select_vp_adaptive_eo_local( ...
                expandedSeedTable, eoCandidates, method, cfg, wIdx);
            result = refine_template_waveform_fit_local( ...
                expandedBundle, expandedSeedTable, expandedEOCandidates, cfg);
            result.CorePass = summarize_phase_safe_reference_local(phaseRef, coreBundle, method);
            result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
                coreBundle, expandedBundle, phaseRef, method);
            seedTable = expandedSeedTable;
            finalEOCandidates = expandedEOCandidates;
            selectionInfo = expandedSelectionInfo;
            bundle = expandedBundle;
        else
            seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
            [finalEOCandidates, selectionInfo] = select_vp_adaptive_eo_local( ...
                seedTable, eoCandidates, method, cfg, wIdx);
            result = refine_template_waveform_fit_local(bundle, seedTable, finalEOCandidates, cfg);
            result.CorePass = summarize_phase_safe_reference_local(phaseRef, coreBundle, method);
            result.PhaseSafeExpansion = summarize_phase_safe_expansion_local( ...
                coreBundle, expandedBundle, phaseRef, method);
        end
    end
else
    seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
    [finalEOCandidates, selectionInfo] = select_vp_adaptive_eo_local( ...
        seedTable, eoCandidates, method, cfg, wIdx);
    result = refine_template_waveform_fit_local(bundle, seedTable, finalEOCandidates, cfg);
end

result.window_id = wIdx;
result.lap_range = Wmap.lap_range;
result.VPSeedTable = struct2table(seedTable);
result.VPSelectedEO = finalEOCandidates(:).';
result.VPSelectionInfo = selectionInfo;
result.PhaseSafeReferenceInfo = phaseRefInfo;

trendRow = struct( ...
    'window_id', wIdx, ...
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
    'sensor_eta_max_abs_mm', max(abs(result.sensor_eta_id), [], 'omitnan'), ...
    'weighted_voltage_rmse', result.weighted_voltage_rmse, ...
    'plain_voltage_rmse', result.plain_voltage_rmse, ...
    'valid_segment_count', result.valid_segment_count, ...
    'point_count', result.point_count, ...
    'core_point_count', coreBundle.point_count, ...
    'added_point_count', result.point_count - coreBundle.point_count, ...
    'phase_safe_reference_source', phaseRefInfo.source, ...
    'phase_safe_reference_reason', phaseRefInfo.reason);

windowResult = struct();
windowResult.window_id = wIdx;
windowResult.lap_range = Wmap.lap_range;
windowResult.time_window = Wmap.time_window;
windowResult.bundle = bundle;
windowResult.seed_table = seedTable;
windowResult.CandidateTable = result.CandidateTable;
windowResult.Result = result;
end

function phaseRef = build_template_inverse_phase_safe_reference_local()
phaseRef = struct();
phaseRef.EO_id = NaN;
phaseRef.A_id = NaN;
phaseRef.phi_id_wrapped = NaN;
phaseRef.dx_c_id = NaN;
phaseRef.sensor_eta_id = [];
phaseRef.weighted_voltage_rmse = NaN;
phaseRef.plain_voltage_rmse = NaN;
phaseRef.reference_source = 'template_inverse_apparent_query';
end

function phaseRef = build_phase_safe_reference_from_seed_local(seedTable, cfg)
seed = seedTable(1);
phaseRef = struct();
phaseRef.EO_id = seed.EO;
phaseRef.A_id = min(abs(seed.A), abs(cfg.amplitude_limit_mm));
phaseRef.phi_id_wrapped = wrap_to_pi_local(seed.phi);
phaseRef.dx_c_id = max(min(seed.dx_c, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phaseRef.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
phaseRef.weighted_voltage_rmse = seed.weighted_voltage_rmse;
phaseRef.plain_voltage_rmse = seed.plain_voltage_rmse;
phaseRef.reference_source = 'linear_vp_seed';
end

function info = empty_phase_safe_reference_info_local()
info = make_phase_safe_reference_info_local('', '', NaN, false);
end

function info = make_phase_safe_reference_info_local(source, reason, sourceWindowId, usedPreviousWindow)
info = struct();
info.source = source;
info.reason = reason;
info.source_window_id = sourceWindowId;
info.used_previous_window = usedPreviousWindow;
end

function [phaseRef, info] = build_current_window_phase_safe_reference_local( ...
    bundle, eoCandidates, cfg, method, wIdx, fallbackReason)
mode = method.phase_safe_reference_mode;
if strcmpi(mode, 'prev_window')
    mode = method.phase_safe_fallback_mode;
end

if strcmpi(mode, 'template_inverse')
    phaseRef = build_template_inverse_phase_safe_reference_local();
    reason = sprintf('%s_then_%s', fallbackReason, mode);
elseif strcmpi(mode, 'linear_vp')
    seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
    phaseRef = build_phase_safe_reference_from_seed_local(seedTable, cfg);
    reason = sprintf('%s_then_%s', fallbackReason, mode);
elseif strcmpi(mode, 'final_core')
    seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates);
    [coreEOCandidates, ~] = select_vp_adaptive_eo_local(seedTable, eoCandidates, method, cfg, wIdx);
    phaseRef = refine_template_waveform_fit_local(bundle, seedTable, coreEOCandidates, cfg);
    phaseRef.reference_source = 'current_core_final_waveform';
    reason = sprintf('%s_then_%s', fallbackReason, mode);
else
    error('Unsupported phase-safe current-window reference mode: %s.', mode);
end

info = make_phase_safe_reference_info_local(phaseRef.reference_source, reason, wIdx, false);
end

function [phaseRef, info] = build_prev_window_phase_safe_reference_local(prevResult, wIdx, method, cfg)
phaseRef = [];
[ok, reason] = is_prev_window_phase_reference_usable_local(prevResult, wIdx, method);
if ~ok
    info = make_phase_safe_reference_info_local('', reason, NaN, false);
    return;
end

phaseRef = struct();
phaseRef.EO_id = prevResult.EO_id;
phaseRef.A_id = min(abs(prevResult.A_id), abs(cfg.amplitude_limit_mm));
phaseRef.phi_id_wrapped = wrap_to_pi_local(prevResult.phi_id_wrapped);
phaseRef.dx_c_id = max(min(prevResult.dx_c_id, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phaseRef.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
if isfield(prevResult, 'sensor_eta_id') && ~isempty(prevResult.sensor_eta_id)
    n = min(numel(cfg.analysis_sensors), numel(prevResult.sensor_eta_id));
    phaseRef.sensor_eta_id(1:n) = prevResult.sensor_eta_id(1:n);
end
phaseRef.weighted_voltage_rmse = prevResult.weighted_voltage_rmse;
phaseRef.plain_voltage_rmse = prevResult.plain_voltage_rmse;
phaseRef.reference_source = 'prev_window_final';
phaseRef.source_window_id = prevResult.window_id;
info = make_phase_safe_reference_info_local(phaseRef.reference_source, ...
    'prev_window_quality_pass', prevResult.window_id, true);
end

function [ok, reason] = is_prev_window_phase_reference_usable_local(prevResult, wIdx, method)
ok = false;
if isempty(prevResult) || ~isstruct(prevResult)
    reason = 'no_previous_window';
    return;
end
if method.phase_safe_refresh_every > 0 && mod(wIdx - 1, method.phase_safe_refresh_every) == 0
    reason = 'scheduled_refresh';
    return;
end
required = {'window_id', 'EO_id', 'A_id', 'phi_id_wrapped', 'dx_c_id', 'weighted_voltage_rmse', 'Coverage'};
for i = 1:numel(required)
    if ~isfield(prevResult, required{i})
        reason = ['previous_missing_', required{i}];
        return;
    end
end
vals = [prevResult.window_id, prevResult.EO_id, prevResult.A_id, ...
    prevResult.phi_id_wrapped, prevResult.dx_c_id, prevResult.weighted_voltage_rmse];
if any(~isfinite(vals))
    reason = 'previous_nonfinite_reference';
    return;
end
if isfield(prevResult.Coverage, 'clamp_fraction') && ...
        prevResult.Coverage.clamp_fraction > method.phase_safe_prev_max_clamp_fraction
    reason = 'previous_clamped';
    return;
end
if isfield(prevResult, 'VPSelectionInfo') && isfield(prevResult.VPSelectionInfo, 'fallback_reason') && ...
        ~isempty(prevResult.VPSelectionInfo.fallback_reason)
    reason = 'previous_vp_fallback';
    return;
end
gapRatio = final_candidate_gap_ratio_local(prevResult);
if isfinite(method.phase_safe_prev_min_final_gap_ratio) && ...
        gapRatio < method.phase_safe_prev_min_final_gap_ratio
    reason = 'previous_final_gap_too_small';
    return;
end
ok = true;
reason = 'prev_window_quality_pass';
end

function gapRatio = final_candidate_gap_ratio_local(result)
gapRatio = inf;
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
gapRatio = (rmse(2) - rmse(1)) / max(rmse(1), eps);
end

function bundle = build_template_observation_bundle_local(Wmap, Template, analysisSensors, method, phaseRef)
if nargin < 5
    phaseRef = [];
end
usePhaseSafe = ~isempty(phaseRef) && method.phase_safe_expansion_enabled;

X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
Yobs = [];
F0all = [];
Fxall = [];
sensorIndexPt = [];
eventTime = [];
windowMeta = struct('sensor_id', {}, 'valid_points', {}, 'shape_rmse', {});

for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    vRaw = D.V(:);
    xRaw = D.x_rel(:);
    tRaw = D.t(:);
    thetaRaw = D.theta(:);
    if isfield(D, 'event_time')
        eventRaw = D.event_time(:);
    else
        eventRaw = NaN(size(tRaw));
    end

    f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw, 'pchip', NaN);
    fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw, 'pchip', NaN);

    if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
        threshold = Tpl.threshold;
    else
        threshold = method.default_sensor_threshold;
    end
    xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];

    if strcmpi(method.pulse_selection_mode, 'all')
        [maskPulse, pulseSegmentCount] = isolate_all_pulses_local(vRaw, threshold);
    else
        [maskPulse, pulseSegmentCount] = isolate_main_pulse_local(vRaw, threshold);
    end
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, method);
    maskDomain = xRaw >= Tpl.x_domain(1) + method.domain_margin_mm & ...
                 xRaw <= Tpl.x_domain(2) - method.domain_margin_mm;
    maskFinite = isfinite(f0) & isfinite(fx0) & isfinite(vRaw) & isfinite(thetaRaw);
    maskBase = maskPulse & maskEffective & maskDomain & maskFinite;

    queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, maskBase, method);
    maskQuerySafe = xRaw >= xDomain(1) + queryGuardMm & xRaw <= xDomain(2) - queryGuardMm;
    maskPhaseSafe = false(size(maskQuerySafe));
    if usePhaseSafe
        maskPhaseSafe = build_phase_safe_query_mask_local( ...
            Tpl, xRaw, vRaw, thetaRaw, xDomain, is, phaseRef, method);
        if nnz(maskBase & maskPhaseSafe) >= 8
            maskQuerySafe = maskPhaseSafe;
        end
    end

    mask = maskBase;
    if strcmpi(method.domain_selection_mode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
    if nnz(mask) < 8
        mask = maskEffective & maskDomain & maskFinite;
        if strcmpi(method.domain_selection_mode, 'hard')
            safeMask = mask & maskQuerySafe;
            if nnz(safeMask) >= 8
                mask = safeMask;
            end
        end
    end
    if nnz(mask) < 8
        continue;
    end

    tSel = tRaw(mask);
    xSel = xRaw(mask);
    vSel = vRaw(mask);
    thetaSel = thetaRaw(mask);
    eventSel = eventRaw(mask);
    f0Sel = f0(mask);
    fxSel = fx0(mask);
    if method.store_y_obs
        yObsSel = xSel - invert_template_voltage_local(Tpl, vSel, xSel);
    else
        yObsSel = [];
    end

    wEdge = build_edge_weight_local(tSel, vSel, method.weight_floor);
    if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
        wTemplate = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), xSel, 'linear', method.weight_floor);
    else
        wTemplate = ones(size(xSel));
    end
    wDomain = build_domain_soft_weight_local(xSel, xDomain, method.domain_soft_margin_mm, method.weight_floor);
    wQuery = build_query_guard_soft_weight_local( ...
        xSel, xDomain, queryGuardMm, method.domain_soft_margin_mm, method.weight_floor);
    wGradient = build_template_gradient_weight_local(Tpl, xSel, method.domain_selection_mode, method.weight_floor);
    wTotal = max(method.weight_floor, wEdge .* wTemplate .* wDomain .* wQuery .* wGradient);
    if max(wTotal) > 0
        wTotal = max(method.weight_floor, wTotal ./ max(wTotal));
    end

    vStatic = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xSel, 'pchip', NaN);
    shapeRmse = sqrt(mean((vSel - vStatic).^2, 'omitnan'));

    X = [X; xSel]; %#ok<AGROW>
    T = [T; tSel]; %#ok<AGROW>
    V = [V; vSel]; %#ok<AGROW>
    S = [S; repmat(sid, numel(tSel), 1)]; %#ok<AGROW>
    W = [W; wTotal]; %#ok<AGROW>
    Theta = [Theta; thetaSel]; %#ok<AGROW>
    eventTime = [eventTime; eventSel]; %#ok<AGROW>
    if method.store_y_obs
        Yobs = [Yobs; yObsSel]; %#ok<AGROW>
    end
    F0all = [F0all; f0Sel]; %#ok<AGROW>
    Fxall = [Fxall; fxSel]; %#ok<AGROW>
    sensorIndexPt = [sensorIndexPt; repmat(is, numel(tSel), 1)]; %#ok<AGROW>

    metaIdx = numel(windowMeta) + 1;
    windowMeta(metaIdx).sensor_id = sid;
    windowMeta(metaIdx).valid_points = numel(tSel);
    windowMeta(metaIdx).pulse_segment_count = pulseSegmentCount;
    windowMeta(metaIdx).dynamic_effective_points = nnz(maskEffective);
    windowMeta(metaIdx).query_safe_points = nnz(maskQuerySafe(mask));
    windowMeta(metaIdx).query_guard_mm = queryGuardMm;
    windowMeta(metaIdx).domain_selection_mode = method.domain_selection_mode;
    windowMeta(metaIdx).domain_soft_margin_mm = method.domain_soft_margin_mm;
    windowMeta(metaIdx).phase_safe_expansion = usePhaseSafe;
    windowMeta(metaIdx).phase_safe_points = nnz(maskBase & maskPhaseSafe);
    windowMeta(metaIdx).template_domain = xDomain;
    windowMeta(metaIdx).selected_x_range = [min(xSel), max(xSel)];
    windowMeta(metaIdx).shape_rmse = shapeRmse;
end

if isempty(T)
    error('No valid waveform points were constructed for window %d.', Wmap.window_id);
end

interpV = cell(numel(analysisSensors), 1);
xDomainBySensor = NaN(numel(analysisSensors), 2);
for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    interpV{is} = griddedInterpolant(Tpl.x_grid(:), Tpl.v_grid(:), 'pchip', 'none');
    xDomainBySensor(is, :) = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
end

bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.EventTime = eventTime;
bundle.Y_obs = Yobs;
bundle.F0 = F0all;
bundle.Fx = Fxall;
bundle.sensor_index = sensorIndexPt;
bundle.sensor_ids = analysisSensors(:).';
bundle.interp_v = interpV;
bundle.x_domain_by_sensor = xDomainBySensor;
bundle.query_guard_mm = max([windowMeta.query_guard_mm], [], 'omitnan');
bundle.overshoot_penalty_weight = method.overshoot_penalty_weight;
bundle.sensor_eta_limit_mm = method.sensor_eta_limit_mm;
bundle.sensor_eta_reg_weight_v_per_mm = method.sensor_eta_reg_weight_v_per_mm;
bundle.window_meta = windowMeta;
bundle.valid_segment_count = sum([windowMeta.pulse_segment_count]);
bundle.point_count = numel(T);
if usePhaseSafe
    bundle.selection_pass = 'phase_safe_expanded';
else
    bundle.selection_pass = 'core_hard_adaptive';
end
bundle.phase_safe_margin_mm = method.phase_safe_margin_mm;
bundle.time_window = [min(T), max(T)];
bundle.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
bundle.rot_rpm_mean = Wmap.rot_rpm_mean;
end

function maskSafe = build_phase_safe_query_mask_local(Tpl, x, v, theta, xDomain, sensorLocalIndex, phaseRef, method)
if isfield(phaseRef, 'reference_source') && strcmpi(phaseRef.reference_source, 'template_inverse_apparent_query')
    maskSafe = build_template_inverse_safe_mask_local(Tpl, x, v, xDomain, method);
    return;
end
A = phaseRef.A_id;
eo = phaseRef.EO_id;
phi = phaseRef.phi_id_wrapped;
dxc = phaseRef.dx_c_id;
eta = 0;
if isfield(phaseRef, 'sensor_eta_id') && numel(phaseRef.sensor_eta_id) >= sensorLocalIndex
    eta = phaseRef.sensor_eta_id(sensorLocalIndex);
end
uEst = A .* sin(eo .* theta(:) + phi);
xQuery = x(:) - dxc - eta - uEst;
margin = method.phase_safe_margin_mm;
maskSafe = xQuery >= xDomain(1) + margin & xQuery <= xDomain(2) - margin;
maskSafe = reshape(maskSafe, size(x));
end

function maskSafe = build_template_inverse_safe_mask_local(Tpl, x, v, xDomain, method)
xStat = invert_template_voltage_local(Tpl, v(:), x(:));
margin = method.phase_safe_margin_mm;
uApp = abs(x(:) - xStat(:));
maskSafe = isfinite(xStat) & ...
    xStat >= xDomain(1) + margin & xStat <= xDomain(2) - margin & ...
    uApp <= method.phase_safe_u_limit_mm + margin;
maskSafe = reshape(maskSafe, size(x));
end

function summary = summarize_phase_safe_reference_local(phaseRef, bundle, method)
summary = struct();
summary.EO_id = phaseRef.EO_id;
summary.A_id = phaseRef.A_id;
summary.dx_c_id = phaseRef.dx_c_id;
summary.phi_id_wrapped = phaseRef.phi_id_wrapped;
summary.weighted_voltage_rmse = phaseRef.weighted_voltage_rmse;
summary.plain_voltage_rmse = phaseRef.plain_voltage_rmse;
summary.point_count = bundle.point_count;
summary.selection_pass = bundle.selection_pass;
summary.reference_mode = method.phase_safe_reference_mode;
summary.reference_source = get_optional_field_local(phaseRef, 'reference_source', method.phase_safe_reference_mode);
end

function summary = summarize_phase_safe_expansion_local(coreBundle, expandedBundle, phaseRef, method)
summary = struct();
summary.enabled = method.phase_safe_expansion_enabled;
summary.reference_EO_id = phaseRef.EO_id;
summary.reference_A_id = phaseRef.A_id;
summary.reference_dx_c_id = phaseRef.dx_c_id;
summary.reference_phi_id_wrapped = phaseRef.phi_id_wrapped;
summary.reference_mode = method.phase_safe_reference_mode;
summary.margin_mm = method.phase_safe_margin_mm;
summary.core_point_count = coreBundle.point_count;
summary.expanded_point_count = expandedBundle.point_count;
summary.added_point_count = expandedBundle.point_count - coreBundle.point_count;
summary.added_fraction_of_core = summary.added_point_count / max(coreBundle.point_count, 1);
end

function value = get_optional_field_local(S, fieldName, defaultValue)
if isfield(S, fieldName)
    value = S.(fieldName);
else
    value = defaultValue;
end
end

function xStat = invert_template_voltage_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diffV = vGrid - vv;
    crossingX = [];
    exactIdx = find(abs(diffV) <= 1e-10);
    if ~isempty(exactIdx)
        crossingX = xGrid(exactIdx);
    end
    for k = 1:numel(diffV)-1
        if ~isfinite(diffV(k)) || ~isfinite(diffV(k+1))
            continue;
        end
        if diffV(k) == 0 || diffV(k) * diffV(k+1) > 0
            continue;
        end
        denom = vGrid(k+1) - vGrid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - vGrid(k)) / denom;
        crossingX(end+1, 1) = xGrid(k) + alpha * (xGrid(k+1) - xGrid(k)); %#ok<AGROW>
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function queryGuardMm = resolve_query_guard_mm_local(Tpl, x, v, baseMask, method)
queryGuardMm = method.query_guard_mm;
if ~strcmpi(method.query_guard_mode, 'adaptive') || nnz(baseMask) < 8
    return;
end
xSel = x(baseMask);
vSel = v(baseMask);
xStat = invert_template_voltage_local(Tpl, vSel, xSel);
uApp = abs(xSel(:) - xStat(:));
uApp = uApp(isfinite(uApp));
if isempty(uApp)
    return;
end
q = min(max(method.query_guard_quantile, 0), 100);
adaptiveGuard = prctile(uApp, q) + method.query_guard_safety_mm;
queryGuardMm = min(max(adaptiveGuard, method.query_guard_min_mm), method.query_guard_max_mm);
end

function eoCandidates = build_eo_candidates_local(rotFreqHz, freqSearchHz, eoPad)
if nargin < 3 || isempty(eoPad)
    eoPad = 0;
end
rotFreqHz = max(rotFreqHz, eps);
freqLo = min(freqSearchHz);
freqHi = max(freqSearchHz);
eoMinStrict = max(1, ceil(freqLo / rotFreqHz));
eoMaxStrict = max(eoMinStrict, floor(freqHi / rotFreqHz));
eoCandidates = eoMinStrict:eoMaxStrict;
if isempty(eoCandidates)
    eoCenter = max(1, round(mean([freqLo, freqHi]) / rotFreqHz));
    eoCandidates = max(1, eoCenter - eoPad):max(1, eoCenter + eoPad);
end
freqCandidates = eoCandidates .* rotFreqHz;
eoCandidates = eoCandidates(freqCandidates >= freqLo & freqCandidates <= freqHi);
if isempty(eoCandidates)
    error('No EO candidates remain inside [%.3f, %.3f] Hz at rot_freq %.6f Hz.', ...
        freqLo, freqHi, rotFreqHz);
end
end

function seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'weighted_voltage_rmse', inf, 'plain_voltage_rmse', inf, ...
    'linear_vp_rmse', inf), numel(eoCandidates), 1);
nSensor = numel(bundle.sensor_ids);

for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo * bundle.Theta);
    c1 = cos(eo * bundle.Theta);
    basis = [-bundle.Fx(:), -bundle.Fx(:) .* s1, -bundle.Fx(:) .* c1];
    wSqrt = sqrt(bundle.W(:));
    y = bundle.V(:) - bundle.F0(:);
    coeff = (basis .* wSqrt) \ (y .* wSqrt);
    dxc = coeff(1);
    eta = zeros(1, nSensor);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y - basis * coeff;
    [obj, plainRmse] = template_synchronous_objective_local([A, phi, dxc, eta], eo, bundle);

    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx_c = dxc;
    rows(i).sensor_eta = eta;
    rows(i).sensor_eta_max_abs = max(abs(eta), [], 'omitnan');
    rows(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    rows(i).plain_voltage_rmse = plainRmse;
    rows(i).linear_vp_rmse = sqrt(sum(bundle.W(:) .* linRes.^2) / max(sum(bundle.W(:)), eps));
end

score = [[rows.weighted_voltage_rmse].', [rows.linear_vp_rmse].', [rows.plain_voltage_rmse].', [rows.EO].'];
[~, order] = sortrows(score, [1, 2, 3]);
seedTable = rows(order);
end

function eoKeep = select_vp_topk_eo_local(seedTable, topK)
if nargin < 2 || isempty(topK) || ~isfinite(topK)
    eoKeep = unique(round([seedTable.EO]), 'stable');
    return;
end
topK = min(max(1, floor(topK)), numel(seedTable));
eoKeep = unique(round([seedTable(1:topK).EO]), 'stable');
end

function [eoKeep, info] = select_vp_adaptive_eo_local(seedTable, eoCandidates, method, cfg, windowId)
allEo = unique(round(eoCandidates(:).'));
topKEo = select_vp_topk_eo_local(seedTable, method.vp_top_k_eo);
eoKeep = intersect(allEo, topKEo, 'stable');

bestRmse = seedTable(1).weighted_voltage_rmse;
if numel(seedTable) >= 2
    secondRmse = seedTable(2).weighted_voltage_rmse;
else
    secondRmse = inf;
end
vpGapRatio = (secondRmse - bestRmse) / max(bestRmse, eps);
linearGapRatio = inf;
if numel(seedTable) >= 2
    linearGapRatio = (seedTable(2).linear_vp_rmse - seedTable(1).linear_vp_rmse) / ...
        max(seedTable(1).linear_vp_rmse, eps);
end

fallbackReason = '';
if windowId <= cfg.vp_all_eo_warmup_windows
    fallbackReason = 'warmup_all_eo';
elseif vpGapRatio < cfg.vp_gap_ratio_fallback
    fallbackReason = 'ambiguous_weighted_vp_gap';
elseif linearGapRatio < cfg.vp_linear_gap_ratio_fallback
    fallbackReason = 'ambiguous_linear_vp_gap';
end
if ~isempty(fallbackReason)
    eoKeep = allEo;
end

info = struct();
info.window_id = windowId;
info.all_eo = allEo;
info.top_k_eo = topKEo;
info.selected_eo = eoKeep;
info.best_seed_eo = seedTable(1).EO;
info.best_seed_weighted_rmse = bestRmse;
info.second_seed_weighted_rmse = secondRmse;
info.vp_gap_ratio = vpGapRatio;
info.linear_gap_ratio = linearGapRatio;
info.fallback_reason = fallbackReason;
end

function result = refine_template_waveform_fit_local(bundle, seedTable, eoCandidates, cfg)
ampLimitMm = abs(cfg.amplitude_limit_mm);
dxcLimitMm = abs(cfg.dx_c_limit_mm);
sensorEtaLimitMm = abs(cfg.sensor_eta_limit_mm);
eoCandidates = unique(round(eoCandidates(:).'));
nSensor = numel(bundle.sensor_ids);
candidate = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'fn_hz', NaN, 'objective', inf, 'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, 'V_pred', [], 'u_est', []), numel(eoCandidates), 1);
best = struct('objective', inf);

for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find([seedTable.EO] == eo, 1, 'first');
    if isempty(seedIdx)
        seed = seedTable(1);
    else
        seed = seedTable(seedIdx);
    end
    eta0 = zeros(1, nSensor);
    if isfield(seed, 'sensor_eta') && ~isempty(seed.sensor_eta)
        eta0(1:min(nSensor, numel(seed.sensor_eta))) = seed.sensor_eta(1:min(nSensor, numel(seed.sensor_eta)));
    end
    seedParams = [seed.A, seed.phi, seed.dx_c, eta0];
    fun = @(p) template_synchronous_objective_local(p, eo, bundle);
    opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
    pOpt = fminsearch(@(p) bounded_synchronous_objective_local( ...
        p, fun, ampLimitMm, dxcLimitMm, sensorEtaLimitMm, nSensor), seedParams, opts);
    pOpt(1) = min(abs(pOpt(1)), ampLimitMm);
    pOpt(2) = wrap_to_pi_local(pOpt(2));
    pOpt(3) = max(min(pOpt(3), dxcLimitMm), -dxcLimitMm);
    pOpt = bound_sensor_eta_params_local(pOpt, sensorEtaLimitMm, nSensor);
    [obj, plainRmse, vPred, uEst, coverage] = template_synchronous_objective_local(pOpt, eo, bundle);

    candidate(i).EO = eo;
    candidate(i).A = pOpt(1);
    candidate(i).phi = pOpt(2);
    candidate(i).dx_c = pOpt(3);
    candidate(i).sensor_eta = pOpt(4:3+nSensor);
    candidate(i).sensor_eta_max_abs = max(abs(candidate(i).sensor_eta), [], 'omitnan');
    candidate(i).fn_hz = eo * bundle.rot_freq_mean_hz;
    candidate(i).objective = obj;
    candidate(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    candidate(i).plain_voltage_rmse = plainRmse;
    candidate(i).clamp_fraction = coverage.clamp_fraction;
    candidate(i).max_query_overshoot_mm = coverage.max_query_overshoot_mm;
    candidate(i).overshoot_penalty = coverage.overshoot_penalty;
    candidate(i).query_guard_mm = bundle.query_guard_mm;
    candidate(i).V_pred = vPred;
    candidate(i).u_est = uEst;

    if obj < best.objective
        best = candidate(i);
    end
end

candidateTable = struct2table(rmfield(candidate, {'V_pred','u_est','sensor_eta'}));
candidateTable = sortrows(candidateTable, {'weighted_voltage_rmse','plain_voltage_rmse','EO'}, ...
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
result.CandidateTable = candidateTable;
result.Coverage = struct( ...
    'clamp_fraction', best.clamp_fraction, ...
    'max_query_overshoot_mm', best.max_query_overshoot_mm, ...
    'query_guard_mm', bundle.query_guard_mm, ...
    'overshoot_penalty', best.overshoot_penalty, ...
    'overshoot_penalty_weight', bundle.overshoot_penalty_weight);
result.metric_note = 'Copied Step03 direct-template bounded waveform fit with shared dx_c.';
end

function obj = bounded_synchronous_objective_local(p, fun, ampLimitMm, dxcLimitMm, sensorEtaLimitMm, nSensor)
pUse = p;
pUse(1) = min(abs(pUse(1)), ampLimitMm);
pUse(2) = wrap_to_pi_local(pUse(2));
pUse(3) = max(min(pUse(3), dxcLimitMm), -dxcLimitMm);
pUse = bound_sensor_eta_params_local(pUse, sensorEtaLimitMm, nSensor);
obj = fun(pUse);
end

function p = bound_sensor_eta_params_local(p, sensorEtaLimitMm, nSensor)
need = 3 + nSensor;
if numel(p) < need
    p(numel(p)+1:need) = 0;
end
p = p(1:need);
eta = p(4:end);
eta = max(min(eta, sensorEtaLimitMm), -sensorEtaLimitMm);
if ~isempty(eta)
    eta = eta - eta(1);
    eta = max(min(eta, sensorEtaLimitMm), -sensorEtaLimitMm);
end
p(4:end) = eta;
end

function [obj, plainRmse, vPred, uEst, coverage] = template_synchronous_objective_local(params, eo, bundle)
A = params(1);
phi = params(2);
dxc = params(3);
params = bound_sensor_eta_params_local(params, abs(bundle.sensor_eta_limit_mm), numel(bundle.sensor_ids));
eta = params(4:end).';
uEst = A .* sin(eo .* bundle.Theta + phi);
etaSample = eta(bundle.sensor_index(:));
xIn = bundle.X - dxc - etaSample - uEst;
vPred = NaN(size(bundle.V));
clamped = false(size(bundle.V));
overshoot = zeros(size(bundle.V));

for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    xLo = bundle.x_domain_by_sensor(is, 1);
    xHi = bundle.x_domain_by_sensor(is, 2);
    xEvalRaw = xIn(mask);
    xEval = min(max(xEvalRaw, xLo), xHi);
    clamped(mask) = xEvalRaw < xLo | xEvalRaw > xHi;
    overshoot(mask) = max(xLo - xEvalRaw, 0) + max(xEvalRaw - xHi, 0);
    vPred(mask) = bundle.interp_v{is}(xEval);
end

valid = isfinite(vPred);
res = bundle.V(valid) - vPred(valid);
overshootPenalty = bundle.overshoot_penalty_weight * sum(bundle.W(valid) .* (overshoot(valid).^2));
etaRegPenalty = bundle.point_count * (bundle.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta.^2))) ^ 2;
obj = sum(bundle.W(valid) .* (res.^2)) + overshootPenalty + etaRegPenalty + 1e6 * nnz(~valid);
plainRmse = sqrt(mean(res.^2));
coverage = struct();
coverage.clamp_fraction = nnz(clamped) / max(numel(clamped), 1);
coverage.max_query_overshoot_mm = max(overshoot, [], 'omitnan');
coverage.overshoot_penalty = overshootPenalty;
coverage.sensor_eta_reg_penalty = etaRegPenalty;
end

function [mask, segmentCount] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segmentCount = numel(starts);
end

function [mask, segmentCount] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segmentCount = 1;
    return;
end

segmentGap = idx(starts(2:end)) - idx(ends(1:end-1));
largeGapThreshold = max(10, round(0.02 * numel(v)));
pulseBreaks = find(segmentGap > largeGapThreshold);
groupStarts = [1; pulseBreaks + 1];
groupEnds = [pulseBreaks; numel(starts)];
for ig = 1:numel(groupStarts)
    segIds = groupStarts(ig):groupEnds(ig);
    bestSeg = segIds(1);
    bestPeak = -inf;
    for iseg = segIds
        seg = idx(starts(iseg):ends(iseg));
        peakVal = max(v(seg));
        if peakVal > bestPeak
            bestPeak = peakVal;
            bestSeg = iseg;
        end
    end
    mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end
segmentCount = numel(groupStarts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end

gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gTplMax = max(gTpl(finite), [], 'omitnan');
if ~isfinite(gTplMax) || gTplMax <= 0
    maskTpl = finite;
else
    maskTpl = gTpl >= method.dynamic_template_gradient_min_ratio * gTplMax;
end

gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if ~isfinite(gTimeMax) || gTimeMax <= 0
    maskTime = false(size(v));
else
    maskTime = gTime >= method.dynamic_time_gradient_min_ratio * gTimeMax;
end

peakQ = min(max(method.dynamic_peak_quantile, 0), 100);
peakLevel = prctile(v(finite), peakQ);
maskPeak = v >= max(threshold, peakLevel);

mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function wEdge = build_edge_weight_local(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_local(x, xDomain, marginMm, floorW)
if marginMm <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMm, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMm, marginMm, floorW)
if queryGuardMm <= 0 || marginMm <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMm - marginMm, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMm, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_local(Tpl, x, domainSelectionMode, floorW)
if ~strcmpi(domainSelectionMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function summary = build_selected_summary_local(selected, source, P)
bundle = selected.bundle;
result = selected.Result;
u = result.u_est(:);
valid = isfinite(u);
if any(valid)
    uMin = min(u(valid), [], 'omitnan');
    uMax = max(u(valid), [], 'omitnan');
    uStd = std(u(valid), 'omitnan');
else
    uMin = NaN;
    uMax = NaN;
    uStd = NaN;
end

summary = table( ...
    P.targetBlade, ...
    string(['S', sprintf('%d', bundle.sensor_ids)]), ...
    selected.window_id, ...
    string(mat2str(selected.lap_range)), ...
    selected.time_window(1), ...
    selected.time_window(2), ...
    string(source.lowDir), ...
    string(source.highDir), ...
    string(source.templateFile), ...
    string(source.directDynamicMapFile), ...
    string(selected.Result.PhaseSafeReferenceInfo.source), ...
    string(selected.Result.PhaseSafeReferenceInfo.reason), ...
    selected.Result.CorePass.point_count, ...
    selected.Result.point_count - selected.Result.CorePass.point_count, ...
    bundle.point_count, ...
    result.EO_id, ...
    result.fn_id, ...
    result.A_id, ...
    result.phi_id_wrapped, ...
    result.dx_c_id, ...
    max(abs(result.sensor_eta_id), [], 'omitnan'), ...
    result.weighted_voltage_rmse, ...
    result.plain_voltage_rmse, ...
    uMax - uMin, ...
    uStd, ...
    uMin, ...
    uMax, ...
    result.Coverage.clamp_fraction, ...
    result.Coverage.max_query_overshoot_mm, ...
    bundle.query_guard_mm, ...
    abs(abs(result.dx_c_id) - P.dxLimitMm) < 1e-4, ...
    abs(abs(result.A_id) - P.amplitudeLimitMm) < 1e-4, ...
    'VariableNames', {'blade','sensorTag','windowId','lapRange', ...
    'timeStartSec','timeEndSec','lowDir','highDir','templateFile','directDynamicMapFile', ...
    'phaseSafeReferenceSource','phaseSafeReferenceReason', ...
    'corePointCount','addedPointCount','pointCount', ...
    'EO','frequencyHz','amplitudeMm','phaseRad','dxMm','sensorEtaMaxAbsMm', ...
    'weightedRmseMv','plainRmseMv','uPeakToPeakMm','uStdMm','uMinMm','uMaxMm', ...
    'clampFraction','maxQueryOvershootMm','queryGuardMm','dxAtLimit','amplitudeAtLimit'});
end

function plot_selected_window_local(selected, Template, figFile)
bundle = selected.bundle;
result = selected.Result;
seedTable = struct2table(selected.seed_table);
seedTable = sortrows(seedTable, {'weighted_voltage_rmse','linear_vp_rmse','EO'}, ...
    {'ascend','ascend','ascend'});

colors = lines(numel(bundle.sensor_ids));
figure('Name', '20241106 direct template selected window', 'Color', 'w', ...
    'Position', [80, 60, 1500, 900], 'NumberTitle', 'off');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    plot_event_segments_local(bundle.T(mask), bundle.V(mask), bundle.EventTime(mask), ...
        colors(is, :), sprintf('P%d', bundle.sensor_ids(is)));
end
xlabel('Time (s)');
ylabel('Voltage (mV)');
title(sprintf('Selected waveform W%d, laps %s', selected.window_id, mat2str(selected.lap_range)));
legend('Location', 'best');

nexttile; hold on; grid on; box on;
for is = 1:numel(Template.Sensor)
    Tpl = Template.Sensor(is);
    plot(Tpl.x_grid, Tpl.v_grid, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('P%d', Tpl.sensor_id));
end
xlabel('x_{rel} (mm)');
ylabel('Template V (mV)');
title('Low-speed template');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
plot(seedTable.EO, seedTable.weighted_voltage_rmse, 'ko-', 'LineWidth', 1.1, ...
    'DisplayName', 'seed weighted RMSE');
xline(result.EO_id, 'r--', sprintf('best EO%d', result.EO_id), ...
    'LabelVerticalAlignment', 'bottom');
xlabel('EO');
ylabel('RMSE (mV)');
title('VP seed scan');

nexttile; hold on; grid on; box on;
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    plot_event_segments_local(bundle.T(mask), result.u_est(mask), bundle.EventTime(mask), ...
        colors(is, :), sprintf('P%d', bundle.sensor_ids(is)));
end
xlabel('Time (s)');
ylabel('u(t) (mm)');
title(sprintf('A %.4f mm, dx_c %.4f mm', result.A_id, result.dx_c_id));
legend('Location', 'best');

for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    nexttile; hold on; grid on; box on;
    mask = bundle.sensor_index == is;
    xRaw = bundle.X(mask);
    vObs = bundle.V(mask);
    xComp = xRaw - result.dx_c_id - result.u_est(mask);
    if ~isempty(result.sensor_eta_id)
        xComp = xComp - result.sensor_eta_id(is);
    end
    Tpl = Template.Sensor(is);
    scatter(xRaw, vObs, 6, [0.65 0.65 0.65], 'filled', ...
        'MarkerFaceAlpha', 0.18, 'DisplayName', 'raw x');
    plot(Tpl.x_grid, Tpl.v_grid, 'r-', 'LineWidth', 1.1, 'DisplayName', 'template');
    scatter(xComp, vObs, 6, colors(is, :), 'filled', ...
        'MarkerFaceAlpha', 0.18, 'DisplayName', 'compensated x');
    xlabel('x_{rel} (mm)');
    ylabel('Voltage (mV)');
    title(sprintf('P%d collapse, clamp %.4f', sid, result.Coverage.clamp_fraction));
    legend('Location', 'best');
end

sgtitle(sprintf('20241106 direct template, window %d, EO %d, WRMSE %.3f mV', ...
    selected.window_id, result.EO_id, result.weighted_voltage_rmse));
saveas(gcf, figFile);
end

function plot_event_segments_local(t, y, eventTime, color, labelText)
t = t(:);
y = y(:);
eventTime = eventTime(:);
if numel(eventTime) ~= numel(t) || all(~isfinite(eventTime))
    plot(t, y, '.', 'Color', color, 'MarkerSize', 4, 'DisplayName', labelText);
    return;
end
events = unique(eventTime(isfinite(eventTime)));
for i = 1:numel(events)
    mask = abs(eventTime - events(i)) < 1e-10;
    [tt, order] = sort(t(mask));
    yy = y(mask);
    yy = yy(order);
    if i == 1
        plot(tt, yy, '-', 'Color', color, 'LineWidth', 0.8, 'DisplayName', labelText);
    else
        plot(tt, yy, '-', 'Color', color, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
end
