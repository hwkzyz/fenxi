%% Step03B_Run_TemplateOnly_VP_Candidate_20250527
% Template-only identification with VP candidate retention and full
% waveform-template refinement.
%
% Same settings:
%   cfg.target_blades        = 1
%   cfg.analysis_sensors     = [1 3] by default; override with STEP03_ANALYSIS_SENSORS
%   cfg.analysis_start_time  = 1.5 s
%   cfg.target_laps          = 20
%   cfg.analysis_win_size    = 3 laps
%   cfg.sliding_step         = 1 lap
%
% Purpose:
%   This is not a best-window selection script. The steady-speed resonance
%   region is evaluated window by window. The frequency is identified by
%   the proposed method in each window, not fixed from the old SG result.
%
% Difference from Step03:
%   Variable projection only generates multiple frequency/initial-value
%   candidates. The final EO/frequency/amplitude/phase are selected after
%   full shifted-template waveform optimization.

clear; clc; close all;

%% Settings
route_dir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1,  6];
sensor_override = strtrim(getenv('STEP03_ANALYSIS_SENSORS'));
if ~isempty(sensor_override)
    parsed_sensors = sscanf(sensor_override, '%d').';
    if isempty(parsed_sensors)
        error('STEP03_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysis_sensors = parsed_sensors;
end
cfg.analysis_start_time = 1.5;
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];

method = struct();
method.ridge_factor = 1e-8;
method.sensor_eta_limit_mm = 0.03;
method.sweep_eta_limits_mm = [0, 0.01, 0.03, 0.06, 0.12];
method.eo_candidates = 3:24;
method.eo_rel_rmse_tolerance = 1.02;
method.eo_scan_offsets_hz = [-0.5, 0, 0.5];
method.local_scan_half_width_hz = [6, 0.8, 0.18];
method.local_scan_step_hz = [0.5, 0.1, 0.02];
method.full_refine_half_width_hz = 0.8;
method.eo_lock_window_count = 3;     % first few 3-rev windows only; no future-window lookahead
method.vp_top_k = 5;
method.lock_vp_top_k = 5;
method.lock_vp_rel_tolerance = 1.01;
method.lock_max_candidates = 3;
method.per_window_eo_candidates = 3;
method.local_vp_top_k = 1;
method.local_max_candidates = 1;
method.full_refine_top_m = 3;
method.lock_include_best_per_eo = true;
method.candidate_refine_levels = 3;
method.candidate_refine_min_amplitude_mm = 0.02;
method.candidate_phase_offsets_rad = 0;
method.candidate_amplitude_scales = 1.00;
method.candidate_d0_offsets_mm = 0;
method.validity_max_vp_rank = 5;
method.validity_max_full_rel_rmse = 1.05;
method.validity_min_valid_fraction = 0.90;
method.validity_max_sensor_rmse_ratio = 2.50;

sg_best_frequency_hz = 580.219489;
sg_best_amplitude_mm = 0.390398;
sg_best_phase_rad = -0.280964;
sg_best_d0_mm = 0.181145;
sg_best_eo = 14;

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
result_dir = fullfile(output_dir, 'identification');
figure_dir = fullfile(output_dir, 'figures');
if exist(result_dir, 'dir') ~= 7; mkdir(result_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_file = find_template_file_for_sensors_local(template_dir, cfg.target_blades, cfg.analysis_sensors);
dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, cfg.target_blades, cfg.analysis_sensors);
result_file = fullfile(result_dir, sprintf('Result_TemplateOnlyVPCandidateSliding_B%d_%s_20250527.mat', ...
    cfg.target_blades, sensor_tag));

loaded_template = load(template_file, 'Template');
loaded_dynamic = load(dynamic_file, 'DynamicMap');
Template = filter_template_sensors_local(loaded_template.Template, cfg.analysis_sensors, sensor_tag);
DynamicMap = filter_dynamic_map_sensors_local(loaded_dynamic.DynamicMap, cfg.analysis_sensors, sensor_tag);

fprintf('\n=== Step03B: template-only VP candidate + full waveform refinement ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Template source: %s\n', template_file);
fprintf('Dynamic map source: %s\n', dynamic_file);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('No gap variation is estimated in this route.\n');
fprintf('EO initialization: integer EO %d-%d with offsets %s Hz\n', ...
    min(method.eo_candidates), max(method.eo_candidates), mat2str(method.eo_scan_offsets_hz));
fprintf('Local frequency scan stages: half-width %s Hz, step %s Hz\n', ...
    mat2str(method.local_scan_half_width_hz), mat2str(method.local_scan_step_hz));
fprintf('Initial EO lock retains up to %d EO candidate(s); after EO lock one VP seed is refined locally.\n', ...
    method.lock_max_candidates);

if DynamicMap.TargetBlade ~= cfg.target_blades
    error('DynamicMap target blade (%d) does not match cfg.target_blades (%d).', ...
        DynamicMap.TargetBlade, cfg.target_blades);
end
if ~isequal(DynamicMap.SensorIDs(:).', cfg.analysis_sensors(:).')
    error('DynamicMap sensors %s do not match cfg.analysis_sensors %s.', ...
        mat2str(DynamicMap.SensorIDs), mat2str(cfg.analysis_sensors));
end
source_cfg = DynamicMap.SourceSettings;
if source_cfg.analysis_start_time ~= cfg.analysis_start_time || ...
        source_cfg.target_laps ~= cfg.target_laps || ...
        source_cfg.analysis_win_size ~= cfg.analysis_win_size || ...
        source_cfg.sliding_step ~= cfg.sliding_step
    error('DynamicMap SourceSettings do not match the Step03 cfg.');
end

num_windows = numel(DynamicMap.Window);
fprintf('DynamicMap contains %d sliding windows.\n', num_windows);
if isfield(DynamicMap, 'OldSGBestWindowID')
    fprintf('Old SG comparison window: %d, laps %s\n', ...
        DynamicMap.OldSGBestWindowID, mat2str(DynamicMap.OldSGBestLapRange));
end

%% Independent per-window EO candidate refinement
fprintf('\nRunning independent per-window EO candidate refinement...\n');
fprintf('Each window keeps %d EO candidate(s); final output is selected by full waveform objective.\n', ...
    method.per_window_eo_candidates);
global_eo = NaN;
global_eo_table = table();
num_eo_lock_windows = 0;
LockCandidateResult = [];

%% Run proposed identification for every sliding window
WindowResult = [];
Trend = table();
previous_solution = [];

for w_idx = 1:num_windows
    Wmap = DynamicMap.Window(w_idx);
    fprintf('\nWindow %02d/%02d, laps %s, time %.6f-%.6f s\n', ...
        w_idx, num_windows, mat2str(Wmap.lap_range), Wmap.time_window(1), Wmap.time_window(2));

    Rw = identify_one_window_top_eo_candidates_local(Wmap, Template, cfg.analysis_sensors, ...
        method, previous_solution);
    Rw.initial_eo = Rw.nominal_EO;
    previous_solution.theta = Rw.theta_best;
    previous_solution.t0 = Rw.t0;
    if w_idx == 1
        WindowResult = Rw;
    else
        WindowResult(w_idx) = Rw;
    end

    Trend = [Trend; table(w_idx, Wmap.lap_range(1), Wmap.lap_range(end), ...
        mean(Wmap.time_window), Wmap.rot_freq_mean_hz, Wmap.rot_rpm_mean, ...
        Rw.frequency_hz, Rw.EO, Rw.amplitude_mm, Rw.phase_rad, ...
        Rw.mean_x0_mm, Rw.weighted_voltage_rmse_V, Rw.point_count, ...
        Rw.Validity.is_valid, Rw.Validity.vp_rank, Rw.Validity.candidate_count, ...
        Rw.Validity.full_rel_rmse, Rw.Validity.valid_fraction, Rw.Validity.sensor_rmse_ratio, ...
        'VariableNames', {'window_id','lap_start','lap_end','window_center_time', ...
        'rot_freq_hz','rot_rpm','frequency_hz','EO','amplitude_mm','phase_rad', ...
        'mean_x0_mm','weighted_rmse_V','point_count','is_valid','vp_rank','candidate_count', ...
        'full_rel_rmse','valid_fraction','sensor_rmse_ratio'})]; %#ok<AGROW>

    fprintf('  Step03B: f %.6f Hz, EO %.4f, A %.6f mm, RMSE %.6f V, candidates %d, valid %d\n', ...
        Rw.frequency_hz, Rw.EO, Rw.amplitude_mm, Rw.weighted_voltage_rmse_V, ...
        height(Rw.CandidateTable), Rw.Validity.is_valid);
end

if isfield(DynamicMap, 'OldSGTrends')
    old_tr = DynamicMap.OldSGTrends;
    Trend.old_sg_frequency_hz = old_tr.Freq(:);
    Trend.old_sg_amplitude_mm = old_tr.Amp(:);
    Trend.old_sg_weighted_rmse_V = old_tr.WeightedRMSE(:);
end

Result = struct();
Result.Method = 'template_only_vp_candidate_full_waveform_refinement_sliding_windows';
Result.AnalysisSettings = cfg;
Result.MethodSettings = method;
Result.DynamicMapFile = dynamic_file;
Result.TargetBlade = cfg.target_blades;
Result.SensorIDs = cfg.analysis_sensors;
Result.SensorTag = sensor_tag;
Result.GlobalEO = global_eo;
Result.InitialEO = global_eo;
Result.GlobalEOSelectionTable = global_eo_table;
Result.EOLock.Mode = 'independent_per_window_eo_candidates_no_global_lock';
Result.EOLock.WindowCount = num_eo_lock_windows;
Result.EOLock.WindowIDs = [];
Result.EOLock.LapRanges = [];
Result.EOLock.LockCandidateResult = LockCandidateResult;
Result.Trend = Trend;
Result.WindowResult = WindowResult;
if isfield(DynamicMap, 'OldSGBestWindowID')
    Result.SG.frequency_hz = sg_best_frequency_hz;
    Result.SG.amplitude_mm = sg_best_amplitude_mm;
    Result.SG.phase_rad = sg_best_phase_rad;
    Result.SG.d0_mm = sg_best_d0_mm;
    Result.SG.EO = sg_best_eo;
    Result.SG.best_window_id = DynamicMap.OldSGBestWindowID;
    Result.SG.best_lap_range = DynamicMap.OldSGBestLapRange;
end

save(result_file, 'Result', '-v7.3');

fprintf('\n=== Proposed sliding-window resonance identification ===\n');
fprintf('Each 3-rev window was identified independently from retained EO candidates.\n');
disp(Trend(:, {'window_id','lap_start','lap_end','frequency_hz','EO','amplitude_mm','weighted_rmse_V'}));

if isfield(DynamicMap, 'OldSGBestWindowID')
    fprintf('\nOld SG comparison:\n');
    fprintf('  comparison window = %d, laps %s\n', DynamicMap.OldSGBestWindowID, mat2str(DynamicMap.OldSGBestLapRange));
    fprintf('  f  = %.6f Hz, A = %.6f mm, EO = %d\n', ...
        sg_best_frequency_hz, sg_best_amplitude_mm, sg_best_eo);
else
    fprintf('\nOld SG comparison skipped for sensor tag %s.\n', sensor_tag);
end

fprintf('Saved result: %s\n', result_file);

%% Visualization: sliding-window trends
fig1 = figure('Name', 'Step03 proposed sliding-window trends', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.08 0.12 0.78 0.68]);
tiledlayout(fig1, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(Trend.window_id, Trend.frequency_hz, 'o-r', 'LineWidth', 1.3); hold on;
if isfield(DynamicMap, 'OldSGBestWindowID')
    yline(sg_best_frequency_hz, 'b:', 'SG best f', 'LineWidth', 1.1);
end
xlabel('Window id'); ylabel('Frequency (Hz)');
title('Frequency trend');
box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');

nexttile;
plot(Trend.window_id, Trend.amplitude_mm, 'o-r', 'LineWidth', 1.3); hold on;
if isfield(DynamicMap, 'OldSGBestWindowID')
    yline(sg_best_amplitude_mm, 'b:', 'SG best A', 'LineWidth', 1.1);
end
xlabel('Window id'); ylabel('Amplitude (mm)');
title('Amplitude trend');
box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');

nexttile;
plot(Trend.window_id, Trend.weighted_rmse_V, 'o-k', 'LineWidth', 1.3); hold on;
if ismember('old_sg_weighted_rmse_V', Trend.Properties.VariableNames)
    plot(Trend.window_id, Trend.old_sg_weighted_rmse_V, 's-b', 'LineWidth', 1.1);
    legend('proposed', 'old SG', 'Location', 'best');
else
    legend('proposed', 'Location', 'best');
end
xlabel('Window id'); ylabel('Weighted RMSE (V)');
title('Window-by-window residual');
box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');

exportgraphics(fig1, fullfile(figure_dir, sprintf('Step03_ProposedSlidingWindowTrends_B%d_%s.png', ...
    cfg.target_blades, sensor_tag)), 'Resolution', 300);

%% Local calculation blocks
function template_file = find_template_file_for_sensors_local(template_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', ...
    target_blade, sensor_tag));
if isfile(template_file)
    return;
end

files = dir(fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_S*_20250527.mat', target_blade)));
for k = 1:numel(files)
    candidate_file = fullfile(files(k).folder, files(k).name);
    loaded = load(candidate_file, 'Template');
    if isfield(loaded.Template, 'SensorIDs') && all(ismember(analysis_sensors, loaded.Template.SensorIDs))
        template_file = candidate_file;
        return;
    end
end
error('Template file for sensors %s was not found. Run Step01 with these sensors or a superset.', mat2str(analysis_sensors));
end

function dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows_20250527.mat', ...
    target_blade, sensor_tag));
if isfile(dynamic_file)
    return;
end

files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_B%d_S*_SlidingWindows_20250527.mat', target_blade)));
for k = 1:numel(files)
    candidate_file = fullfile(files(k).folder, files(k).name);
    loaded = load(candidate_file, 'DynamicMap');
    if isfield(loaded.DynamicMap, 'SensorIDs') && all(ismember(analysis_sensors, loaded.DynamicMap.SensorIDs))
        dynamic_file = candidate_file;
        return;
    end
end
error('DynamicMap file for sensors %s was not found. Run Step02 with these sensors or a superset.', mat2str(analysis_sensors));
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

function R = identify_one_window_candidates_local(Wmap, Template, analysis_sensors, method, global_eo, previous_solution, include_best_per_eo)
[t_all, x_all, V_all, W_all, S_all, F0_all, Fx_all] = pack_window_data_local(Wmap, Template, analysis_sensors);
t0 = min(t_all);
t_rel = t_all - t0;
W_all = W_all ./ max(W_all);
sw_all = sqrt(W_all);

V_norm = NaN(size(V_all));
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    idx = S_all == sid;
    H_aff = [ones(nnz(idx), 1), F0_all(idx)];
    beta_aff = (H_aff .* sw_all(idx)) \ (V_all(idx) .* sw_all(idx));
    gain0 = beta_aff(2);
    if abs(gain0) < 1e-6
        gain0 = 1;
    end
    V_norm(idx) = (V_all(idx) - beta_aff(1)) ./ gain0;
end
DV = V_norm - F0_all;

if isempty(global_eo)
    eo_list = method.eo_candidates(:).';
else
    eo_list = global_eo(:).';
end

vp_rows = table();
for eo = eo_list
    f_center = eo * Wmap.rot_freq_mean_hz;
    [freq_trace, rmse_trace, beta_trace] = adaptive_first_order_scan_local( ...
        f_center, t_rel, DV, Fx_all, x_all, V_all, W_all, S_all, Template, analysis_sensors, method);
    n = numel(freq_trace);
    vp_rows = [vp_rows; table(repmat(Wmap.window_id, n, 1), repmat(eo, n, 1), ...
        freq_trace(:), rmse_trace(:), beta_trace, ...
        'VariableNames', {'window_id','nominal_EO','vp_frequency_hz','vp_rmse','vp_beta'})]; %#ok<AGROW>
end

CandidateVP = retain_vp_candidates_local(vp_rows, method, include_best_per_eo || numel(eo_list) > 1, Wmap, global_eo, previous_solution);
num_candidates = height(CandidateVP);
CandidateDetail = repmat(empty_candidate_detail_local(), num_candidates, 1);
refined_frequency_hz = NaN(num_candidates, 1);
refined_EO = NaN(num_candidates, 1);
amplitude_mm = NaN(num_candidates, 1);
phase_rad = NaN(num_candidates, 1);
mean_x0_mm = NaN(num_candidates, 1);
weighted_rmse_V = NaN(num_candidates, 1);

for ic = 1:num_candidates
    beta = CandidateVP.vp_beta(ic, :).';
    x0_init = beta(1:numel(analysis_sensors));
    d0_init = mean(x0_init, 'omitnan');
    A_init = max(method.candidate_refine_min_amplitude_mm, hypot(beta(end-1), beta(end)));
    phi_init = atan2(beta(end), beta(end-1));

    theta_init = [CandidateVP.vp_frequency_hz(ic), A_init, phi_init, d0_init, ...
        zeros(1, numel(analysis_sensors))];
    if method.full_refine_top_m > 0 && ic <= method.full_refine_top_m
        [theta_best, best_rmse, best_pred, best_sensor_fit] = refine_full_template_candidate_local( ...
            CandidateVP.vp_frequency_hz(ic), A_init, phi_init, d0_init, t0, t_rel, ...
            x_all, V_all, W_all, S_all, Template, analysis_sensors, method, previous_solution);
    elseif method.full_refine_top_m < 0
        theta_best = theta_init;
        best_rmse = CandidateVP.vp_rmse(ic);
        best_pred = NaN(size(V_all));
        best_sensor_fit = table();
    else
        theta_best = theta_init;
        [best_rmse, best_pred, best_sensor_fit] = full_template_objective_local( ...
            theta_best, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors);
    end

    detail = empty_candidate_detail_local();
    detail.index = ic;
    detail.nominal_EO = CandidateVP.nominal_EO(ic);
    detail.vp_frequency_hz = CandidateVP.vp_frequency_hz(ic);
    detail.vp_rmse = CandidateVP.vp_rmse(ic);
    detail.vp_beta = beta(:).';
    detail.theta_best = theta_best;
    detail.weighted_voltage_rmse_V = best_rmse;
    detail.V_pred = best_pred;
    detail.sensor_fit = best_sensor_fit;
    detail.u_vibration = theta_best(2) * sin(2*pi*theta_best(1)*t_rel + theta_best(3));
    detail.u_total = compute_u_total_local(theta_best, t_rel, S_all, analysis_sensors);
    CandidateDetail(ic) = detail;

    refined_frequency_hz(ic) = theta_best(1);
    refined_EO(ic) = theta_best(1) / Wmap.rot_freq_mean_hz;
    amplitude_mm(ic) = theta_best(2);
    phase_rad(ic) = theta_best(3);
    mean_x0_mm(ic) = theta_best(4);
    weighted_rmse_V(ic) = best_rmse;
end

CandidateTable = table((1:num_candidates).', CandidateVP.nominal_EO, CandidateVP.vp_frequency_hz, ...
    CandidateVP.vp_rmse, CandidateVP.protected_candidate, refined_frequency_hz, refined_EO, amplitude_mm, phase_rad, ...
    mean_x0_mm, weighted_rmse_V, ...
    'VariableNames', {'candidate_id','nominal_EO','vp_frequency_hz','vp_rmse','protected_candidate', ...
    'frequency_hz','EO','amplitude_mm','phase_rad','mean_x0_mm','weighted_rmse_V'});
CandidateTable = sortrows(CandidateTable, 'weighted_rmse_V', 'ascend');
CandidateDetail = CandidateDetail(CandidateTable.candidate_id);
CandidateTable.candidate_id = (1:height(CandidateTable)).';
for ic = 1:numel(CandidateDetail)
    CandidateDetail(ic).index = ic;
end

R = base_window_result_local(Wmap, t_all, t0, t_rel, x_all, V_all, W_all, S_all);
R.CandidateTable = CandidateTable;
R.CandidateDetail = CandidateDetail;
R.VPTrace = vp_rows;
if method.full_refine_top_m < 0
    R.selected_candidate_id = NaN;
    R.nominal_EO = NaN;
    R.Validity = struct();
else
    R = apply_candidate_to_result_local(R, 1, Wmap.rot_freq_mean_hz, analysis_sensors);
    R.nominal_EO = R.CandidateTable.nominal_EO(R.selected_candidate_id);
    R.Validity = evaluate_vibration_validity_local(R, method);
end
end

function CandidateVP = retain_vp_candidates_local(vp_rows, method, include_best_per_eo, ~, ~, ~)
vp_rows = sortrows(vp_rows, 'vp_rmse', 'ascend');
vp_rows.protected_candidate = false(height(vp_rows), 1);
if include_best_per_eo && method.lock_include_best_per_eo
    top_k = method.lock_vp_top_k;
    eo_best = table();
    eo_values = unique(vp_rows.nominal_EO(:).');
    for eo = eo_values
        idx = find(vp_rows.nominal_EO == eo);
        [~, best_local] = min(vp_rows.vp_rmse(idx));
        eo_best = [eo_best; vp_rows(idx(best_local), :)]; %#ok<AGROW>
    end
    eo_best = sortrows(eo_best, 'vp_rmse', 'ascend');
    min_rmse = min(eo_best.vp_rmse);
    protect = eo_best.vp_rmse <= method.lock_vp_rel_tolerance * max(min_rmse, eps);
    keep = false(height(eo_best), 1);
    keep(1:min(top_k, height(eo_best))) = true;
    keep = keep | protect;
    eo_best.protected_candidate = protect;
    CandidateVP = cap_candidate_table_local(eo_best, keep, method.lock_max_candidates);
else
    top_k = method.local_vp_top_k;
    keep = false(height(vp_rows), 1);
    keep(1:min(top_k, height(vp_rows))) = true;
    CandidateVP = cap_candidate_table_local(vp_rows, keep, method.local_max_candidates);
end

freq_key = round(CandidateVP.vp_frequency_hz * 100) / 100;
[~, unique_idx] = unique([CandidateVP.nominal_EO, freq_key], 'rows', 'stable');
CandidateVP = CandidateVP(unique_idx, :);
end

function CandidateVP = cap_candidate_table_local(rows, keep, max_candidates)
protected_idx = find(keep & rows.protected_candidate);
regular_idx = find(keep & ~rows.protected_candidate);
selected_idx = [protected_idx; regular_idx];
if numel(selected_idx) > max_candidates
    selected_idx = selected_idx(1:max_candidates);
end
CandidateVP = rows(selected_idx, :);
CandidateVP = sortrows(CandidateVP, 'vp_rmse', 'ascend');
end

function R = identify_one_window_top_eo_candidates_local(Wmap, Template, analysis_sensors, method, previous_solution)
C = integer_eo_shifted_template_scan_one_window_local(Wmap, Template, analysis_sensors, method);
eo_values = unique(C.eo_candidates(:).');
eo_rmse = nan(numel(eo_values), 1);
for ie = 1:numel(eo_values)
    idx = C.eo_candidates == eo_values(ie);
    eo_rmse(ie) = min(C.rmse(idx));
end
[~, order] = sort(eo_rmse, 'ascend');
top_eo = eo_values(order(1:min(method.per_window_eo_candidates, numel(order))));
R = identify_one_window_candidates_local(Wmap, Template, analysis_sensors, ...
    method, top_eo, previous_solution, false);
R.CoarseEOScan = C;
end

function R = select_candidate_by_locked_eo_local(R, locked_eo, rot_freq_mean_hz, method)
idx = find(R.CandidateTable.nominal_EO == locked_eo, 1);
if isempty(idx)
    [~, idx] = min(abs(R.CandidateTable.EO - locked_eo));
    R.locked_eo_missing_from_candidates = true;
else
    R.locked_eo_missing_from_candidates = false;
end
R = apply_candidate_to_result_local(R, idx, rot_freq_mean_hz, []);
R.selected_locked_EO = locked_eo;
R.nominal_EO = R.CandidateTable.nominal_EO(idx);
R.Validity = evaluate_vibration_validity_local(R, method);
end

function C = integer_eo_shifted_template_scan_one_window_local(Wmap, Template, analysis_sensors, method)
[t_all, x_all, V_all, W_all, S_all, F0_all, Fx_all] = pack_window_data_local(Wmap, Template, analysis_sensors);
t_rel = t_all - min(t_all);
W_all = W_all ./ max(W_all);
sw_all = sqrt(W_all);

V_norm = NaN(size(V_all));
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    idx = S_all == sid;
    H_aff = [ones(nnz(idx), 1), F0_all(idx)];
    beta_aff = (H_aff .* sw_all(idx)) \ (V_all(idx) .* sw_all(idx));
    gain0 = beta_aff(2);
    if abs(gain0) < 1e-6
        gain0 = 1;
    end
    V_norm(idx) = (V_all(idx) - beta_aff(1)) ./ gain0;
end
DV = V_norm - F0_all;

eo_candidates = method.eo_candidates(:);
offsets = method.eo_scan_offsets_hz(:).';
freq_grid = [];
eo_for_freq = [];
for ie = 1:numel(eo_candidates)
    freq_grid = [freq_grid; eo_candidates(ie) * Wmap.rot_freq_mean_hz + offsets(:)]; %#ok<AGROW>
    eo_for_freq = [eo_for_freq; repmat(eo_candidates(ie), numel(offsets), 1)]; %#ok<AGROW>
end

rmse = NaN(numel(freq_grid), 1);
for jf = 1:numel(freq_grid)
    f = freq_grid(jf);
    beta = solve_first_order_beta_local(f, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors, method.ridge_factor);
    rmse(jf) = evaluate_shifted_template_rmse_local( ...
        beta, f, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors);
end

C = struct();
C.window_id = Wmap.window_id;
C.eo_candidates = eo_for_freq;
C.freq_grid = freq_grid;
C.rmse = rmse;
[C.best_rmse, best_idx] = min(rmse);
C.best_frequency_hz = freq_grid(best_idx);
end

function [global_eo, eo_table] = select_global_eo_from_coarse_scans_local(CoarseScan, Window, method)
eo_candidates = method.eo_candidates(:);
support_count = zeros(numel(eo_candidates), 1);
median_rel_rmse = nan(numel(eo_candidates), 1);
mean_rel_rmse = nan(numel(eo_candidates), 1);
score = nan(numel(eo_candidates), 1);
for ie = 1:numel(eo_candidates)
    eo = eo_candidates(ie);
    rel = nan(numel(CoarseScan), 1);
    for iw = 1:numel(CoarseScan)
        fg = CoarseScan(iw).freq_grid(:);
        rm = CoarseScan(iw).rmse(:);
        eo_id = CoarseScan(iw).eo_candidates(:);
        f_target = eo * Window(iw).rot_freq_mean_hz;
        idx = find(eo_id == eo);
        [best_eo_rmse, best_local_idx] = min(rm(idx));
        f_err = abs(fg(idx(best_local_idx)) - f_target);
        rel(iw) = best_eo_rmse / max(min(rm), eps) + 1e-3 * f_err;
    end
    support_count(ie) = sum(rel <= method.eo_rel_rmse_tolerance);
    median_rel_rmse(ie) = median(rel, 'omitnan');
    mean_rel_rmse(ie) = mean(rel, 'omitnan');
    score(ie) = median_rel_rmse(ie) + 0.02 * (numel(CoarseScan) - support_count(ie));
end
eo_table = table(eo_candidates, support_count, median_rel_rmse, mean_rel_rmse, score, ...
    'VariableNames', {'EO','support_count','median_rel_rmse','mean_rel_rmse','score'});
eo_table = sortrows(eo_table, {'support_count','median_rel_rmse'}, {'descend','ascend'});
global_eo = eo_table.EO(1);
end

function [sequence_eo, eo_table] = select_sequence_eo_from_candidate_scans_local(CandidateScanResult, method)
eo_candidates = method.eo_candidates(:);
support_count = zeros(numel(eo_candidates), 1);
median_vp_rank = nan(numel(eo_candidates), 1);
mean_vp_rank = nan(numel(eo_candidates), 1);
mean_vp_rel_rmse = nan(numel(eo_candidates), 1);
score = nan(numel(eo_candidates), 1);
for ie = 1:numel(eo_candidates)
    eo = eo_candidates(ie);
    ranks = nan(numel(CandidateScanResult), 1);
    rel_rmse = nan(numel(CandidateScanResult), 1);
    for iw = 1:numel(CandidateScanResult)
        T = CandidateScanResult(iw).CandidateTable;
        idx = find(T.nominal_EO == eo, 1);
        if isempty(idx)
            continue;
        end
        [~, order] = sort(T.vp_rmse, 'ascend');
        ranks(iw) = find(order == idx, 1);
        rel_rmse(iw) = T.vp_rmse(idx) / max(min(T.vp_rmse), eps);
    end
    support_count(ie) = sum(isfinite(ranks));
    median_vp_rank(ie) = median(ranks, 'omitnan');
    mean_vp_rank(ie) = mean(ranks, 'omitnan');
    mean_vp_rel_rmse(ie) = mean(rel_rmse, 'omitnan');
    score(ie) = mean_vp_rel_rmse(ie) + 0.08 * max(median_vp_rank(ie) - 1, 0) + ...
        0.03 * (numel(CandidateScanResult) - support_count(ie));
end
eo_table = table(eo_candidates, support_count, median_vp_rank, mean_vp_rank, mean_vp_rel_rmse, score, ...
    'VariableNames', {'EO','support_count','median_vp_rank','mean_vp_rank','mean_vp_rel_rmse','score'});
eo_table = sortrows(eo_table, {'support_count','score'}, {'descend','ascend'});
sequence_eo = eo_table.EO(1);
end

function R = refine_one_window_for_sequence_eo_local(Rscan, Wmap, Template, analysis_sensors, method, sequence_eo, previous_solution)
T = Rscan.CandidateTable;
idx = find(T.nominal_EO == sequence_eo, 1);
if isempty(idx)
    idx = 1;
end
[t_all, x_all, V_all, W_all, S_all] = deal(Rscan.t, Rscan.x_rel, Rscan.V_meas, Rscan.W, Rscan.sensor_id);
t0 = Rscan.t0;
t_rel = Rscan.t_rel;
beta = Rscan.CandidateDetail(idx).vp_beta(:);
d0_init = mean(beta(1:numel(analysis_sensors)), 'omitnan');
A_init = max(method.candidate_refine_min_amplitude_mm, hypot(beta(end-1), beta(end)));
phi_init = atan2(beta(end), beta(end-1));
[theta_best, best_rmse, best_pred, best_sensor_fit] = refine_full_template_candidate_local( ...
    T.vp_frequency_hz(idx), A_init, phi_init, d0_init, t0, t_rel, ...
    x_all, V_all, W_all, S_all, Template, analysis_sensors, method, previous_solution);

R = Rscan;
detail = R.CandidateDetail(idx);
detail.theta_best = theta_best;
detail.weighted_voltage_rmse_V = best_rmse;
detail.V_pred = best_pred;
detail.sensor_fit = best_sensor_fit;
detail.u_vibration = theta_best(2) * sin(2*pi*theta_best(1)*t_rel + theta_best(3));
detail.u_total = compute_u_total_local(theta_best, t_rel, S_all, analysis_sensors);
R.CandidateDetail(idx) = detail;
R.CandidateTable.frequency_hz(idx) = theta_best(1);
R.CandidateTable.EO(idx) = theta_best(1) / Wmap.rot_freq_mean_hz;
R.CandidateTable.amplitude_mm(idx) = theta_best(2);
R.CandidateTable.phase_rad(idx) = theta_best(3);
R.CandidateTable.mean_x0_mm(idx) = theta_best(4);
R.CandidateTable.weighted_rmse_V(idx) = best_rmse;
R = apply_candidate_to_result_local(R, idx, Wmap.rot_freq_mean_hz, analysis_sensors);
R.selected_locked_EO = sequence_eo;
R.nominal_EO = sequence_eo;
R.Validity = evaluate_vibration_validity_local(R, method);
end

function [global_eo, eo_table] = select_global_eo_from_refined_candidates_local(LockCandidateResult, Window, method)
eo_candidates = method.eo_candidates(:);
support_count = zeros(numel(eo_candidates), 1);
median_rel_rmse = nan(numel(eo_candidates), 1);
mean_rel_rmse = nan(numel(eo_candidates), 1);
score = nan(numel(eo_candidates), 1);
for ie = 1:numel(eo_candidates)
    eo = eo_candidates(ie);
    rel = nan(numel(LockCandidateResult), 1);
    for iw = 1:numel(LockCandidateResult)
        T = LockCandidateResult(iw).CandidateTable;
        idx = find(T.nominal_EO == eo);
        if isempty(idx)
            continue;
        end
        [best_eo_rmse, best_local_idx] = min(T.vp_rmse(idx));
        f_target = eo * Window(iw).rot_freq_mean_hz;
        f_err = abs(T.vp_frequency_hz(idx(best_local_idx)) - f_target);
        rel(iw) = best_eo_rmse / max(min(T.vp_rmse), eps) + 1e-3 * f_err;
    end
    support_count(ie) = sum(rel <= method.eo_rel_rmse_tolerance);
    median_rel_rmse(ie) = median(rel, 'omitnan');
    mean_rel_rmse(ie) = mean(rel, 'omitnan');
    score(ie) = median_rel_rmse(ie) + 0.02 * (numel(LockCandidateResult) - support_count(ie));
end
eo_table = table(eo_candidates, support_count, median_rel_rmse, mean_rel_rmse, score, ...
    'VariableNames', {'EO','support_count','median_rel_rmse','mean_rel_rmse','score'});
eo_table = sortrows(eo_table, {'support_count','median_rel_rmse'}, {'descend','ascend'});
global_eo = eo_table.EO(1);
end

function R = select_best_refined_candidate_for_eo_local(Rall, global_eo, rot_freq_mean_hz, method)
T = Rall.CandidateTable;
idx = find(T.nominal_EO == global_eo);
if isempty(idx)
    [~, idx] = min(abs(T.EO - global_eo));
end
[~, best_local] = min(T.weighted_rmse_V(idx));
R = apply_candidate_to_result_local(Rall, idx(best_local), rot_freq_mean_hz, []);
R.selected_locked_EO = global_eo;
R.Validity = evaluate_vibration_validity_local(R, method);
end

function R = base_window_result_local(Wmap, t_all, t0, t_rel, x_all, V_all, W_all, S_all)
R = struct();
R.window_id = Wmap.window_id;
R.lap_range = Wmap.lap_range;
R.time_window = Wmap.time_window;
R.point_count = numel(t_all);
R.t = t_all;
R.t0 = t0;
R.t_rel = t_rel;
R.x_rel = x_all;
R.V_meas = V_all;
R.W = W_all;
R.sensor_id = S_all;
R.frequency_hz = NaN;
R.EO = NaN;
R.amplitude_mm = NaN;
R.phase_rad = NaN;
R.mean_x0_mm = NaN;
R.sensor_x0_mm = [];
R.weighted_voltage_rmse_V = NaN;
R.sensor_fit = table();
R.V_pred = NaN(size(V_all));
R.u_total = NaN(size(t_all));
R.u_vibration = NaN(size(t_all));
R.theta_best = [];
R.selected_candidate_id = NaN;
R.selected_locked_EO = NaN;
R.CandidateTable = table();
R.CandidateDetail = empty_candidate_detail_local();
R.VPTrace = table();
R.Validity = struct();
end

function R = apply_candidate_to_result_local(R, candidate_idx, rot_freq_mean_hz, analysis_sensors)
detail = R.CandidateDetail(candidate_idx);
theta_best = detail.theta_best;
R.frequency_hz = theta_best(1);
R.EO = theta_best(1) / rot_freq_mean_hz;
R.amplitude_mm = theta_best(2);
R.phase_rad = theta_best(3);
R.mean_x0_mm = theta_best(4);
R.theta_best = theta_best;
R.weighted_voltage_rmse_V = detail.weighted_voltage_rmse_V;
R.sensor_fit = detail.sensor_fit;
R.V_pred = detail.V_pred;
R.u_total = detail.u_total;
R.u_vibration = detail.u_vibration;
R.selected_candidate_id = candidate_idx;
eta = theta_best(5:end);
if isempty(analysis_sensors)
    analysis_sensors = unique(R.sensor_id(:).', 'stable');
end
R.sensor_x0_mm = theta_best(4) + eta(:).';
for is = 1:numel(analysis_sensors)
    R.sensor_fit.x0_mm(R.sensor_fit.sensor_id == analysis_sensors(is)) = R.sensor_x0_mm(is);
end
end

function Validity = evaluate_vibration_validity_local(R, method)
T = R.CandidateTable;
selected = R.selected_candidate_id;
best_vp = min(T.vp_rmse);
best_full = min(T.weighted_rmse_V);
Validity = struct();
Validity.candidate_count = height(T);
Validity.vp_rank = find(sort(T.vp_rmse, 'ascend') == T.vp_rmse(selected), 1, 'first');
Validity.full_rank = find(sort(T.weighted_rmse_V, 'ascend') == T.weighted_rmse_V(selected), 1, 'first');
Validity.vp_rel_rmse = T.vp_rmse(selected) / max(best_vp, eps);
Validity.full_rel_rmse = T.weighted_rmse_V(selected) / max(best_full, eps);
Validity.valid_fraction = mean(isfinite(R.V_pred));
if isempty(R.sensor_fit) || ~ismember('rmse', R.sensor_fit.Properties.VariableNames)
    Validity.sensor_rmse_ratio = NaN;
else
    sensor_rmse = R.sensor_fit.rmse;
    Validity.sensor_rmse_ratio = max(sensor_rmse) / max(min(sensor_rmse), eps);
end
Validity.is_valid = Validity.vp_rank <= method.validity_max_vp_rank && ...
    Validity.valid_fraction >= method.validity_min_valid_fraction && ...
    (isnan(Validity.sensor_rmse_ratio) || Validity.sensor_rmse_ratio <= method.validity_max_sensor_rmse_ratio);
if isfield(R, 'locked_eo_missing_from_candidates') && R.locked_eo_missing_from_candidates
    Validity.is_valid = false;
end
end

function detail = empty_candidate_detail_local()
detail = struct('index', NaN, 'nominal_EO', NaN, 'vp_frequency_hz', NaN, ...
    'vp_rmse', NaN, 'vp_beta', [], 'theta_best', [], ...
    'weighted_voltage_rmse_V', NaN, 'V_pred', [], 'sensor_fit', table(), ...
    'u_total', [], 'u_vibration', []);
end

function u_total = compute_u_total_local(theta, t_rel, S_all, analysis_sensors)
f = theta(1);
A = theta(2);
phi = theta(3);
d0 = theta(4);
eta = theta(5:end);
u_total = d0 + A * sin(2*pi*f*t_rel + phi);
for is = 1:numel(analysis_sensors)
    u_total(S_all == analysis_sensors(is)) = u_total(S_all == analysis_sensors(is)) + eta(is);
end
end

function [theta_best, best_rmse, best_pred, best_sensor_fit] = refine_full_template_candidate_local( ...
    f0, A0, phi0, d00, t0, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors, method, previous_solution) %#ok<INUSD>
eta_limit = method.sensor_eta_limit_mm;
f_half_width = method.full_refine_half_width_hz;
theta_lb = [f0 - f_half_width, 0.00, -pi, -0.60, -eta_limit * ones(1, numel(analysis_sensors))];
theta_ub = [f0 + f_half_width, 0.80,  pi,  0.60,  eta_limit * ones(1, numel(analysis_sensors))];

phi_starts = atan2(sin(phi0 + method.candidate_phase_offsets_rad), ...
    cos(phi0 + method.candidate_phase_offsets_rad));
A_starts = max(method.candidate_refine_min_amplitude_mm, A0 * method.candidate_amplitude_scales);
d0_starts = d00 + method.candidate_d0_offsets_mm;
if ~isempty(previous_solution)
    prev_theta = previous_solution.theta;
    prev_phi = atan2(sin(prev_theta(3) + 2*pi*prev_theta(1)*(t0 - previous_solution.t0)), ...
        cos(prev_theta(3) + 2*pi*prev_theta(1)*(t0 - previous_solution.t0)));
    phi_starts = unique(round([phi_starts, prev_phi], 8), 'stable');
    A_starts = unique(round([A_starts, prev_theta(2)], 8), 'stable');
    d0_starts = unique(round([d0_starts, prev_theta(4)], 8), 'stable');
end

theta_best = [f0, max(A0, method.candidate_refine_min_amplitude_mm), phi0, d00, zeros(1, numel(analysis_sensors))];
theta_best = min(max(theta_best, theta_lb), theta_ub);
best_rmse = inf;
best_pred = NaN(size(V_all));
best_sensor_fit = table();
for ia = 1:numel(A_starts)
    for ip = 1:numel(phi_starts)
        for id0 = 1:numel(d0_starts)
            theta_try = [f0, A_starts(ia), phi_starts(ip), d0_starts(id0), zeros(1, numel(analysis_sensors))];
            theta_try = min(max(theta_try, theta_lb), theta_ub);
            [rmse_try, pred_try, sensor_fit_try] = full_template_objective_local( ...
                theta_try, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors);
            if rmse_try < best_rmse
                best_rmse = rmse_try;
                theta_best = theta_try;
                best_pred = pred_try;
                best_sensor_fit = sensor_fit_try;
            end
        end
    end
end

step = [0.35, 0.035, 0.20, 0.035, 0.02 * ones(1, numel(analysis_sensors))];
for level = 1:method.candidate_refine_levels
    improved = true;
    while improved
        improved = false;
        candidate_set = theta_best;
        for j = 1:numel(theta_best)
            cand_p = theta_best; cand_p(j) = cand_p(j) + step(j);
            cand_m = theta_best; cand_m(j) = cand_m(j) - step(j);
            candidate_set = [candidate_set; cand_p; cand_m]; %#ok<AGROW>
        end
        for ic = 1:size(candidate_set, 1)
            theta_try = min(max(candidate_set(ic, :), theta_lb), theta_ub);
            theta_try(3) = atan2(sin(theta_try(3)), cos(theta_try(3)));
            [rmse_try, pred_try, sensor_fit_try] = full_template_objective_local( ...
                theta_try, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors);
            if rmse_try + 1e-10 < best_rmse
                best_rmse = rmse_try;
                theta_best = theta_try;
                best_pred = pred_try;
                best_sensor_fit = sensor_fit_try;
                improved = true;
            end
        end
    end
    step = step * 0.5;
end
end

function [freq_trace, rmse_trace, beta_trace, best_f, best_beta] = adaptive_first_order_scan_local( ...
    f_center, t_rel, DV, Fx_all, x_all, V_all, W_all, S_all, Template, analysis_sensors, method)
freq_trace = [];
rmse_trace = [];
beta_trace = [];
center = f_center;
for stage = 1:numel(method.local_scan_half_width_hz)
    half_width = method.local_scan_half_width_hz(stage);
    step_hz = method.local_scan_step_hz(stage);
    freq_grid = (center - half_width):step_hz:(center + half_width);
    stage_rmse = NaN(numel(freq_grid), 1);
    stage_beta = NaN(numel(freq_grid), numel(analysis_sensors) + 2);
    for jf = 1:numel(freq_grid)
        f = freq_grid(jf);
        beta = solve_first_order_beta_local(f, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors, method.ridge_factor);
        stage_beta(jf, :) = beta(:).';
        stage_rmse(jf) = evaluate_first_order_projection_rmse_local( ...
            beta, f, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors);
    end
    freq_trace = [freq_trace; freq_grid(:)]; %#ok<AGROW>
    rmse_trace = [rmse_trace; stage_rmse(:)]; %#ok<AGROW>
    beta_trace = [beta_trace; stage_beta]; %#ok<AGROW>
    [~, best_stage_idx] = min(stage_rmse);
    center = freq_grid(best_stage_idx);
end

[~, best_idx] = min(rmse_trace);
best_f = freq_trace(best_idx);
best_beta = beta_trace(best_idx, :).';
end

function [t_all, x_all, V_all, W_all, S_all, F0_all, Fx_all] = pack_window_data_local(Wmap, Template, analysis_sensors)
t_all = [];
x_all = [];
V_all = [];
W_all = [];
S_all = [];
F0_all = [];
Fx_all = [];
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
    F0 = interp1(Tpl.x_grid, Tpl.v_grid, D.x_rel, 'pchip', NaN);
    Fx = interp1(Tpl.x_grid, Tpl.dv_dx, D.x_rel, 'pchip', NaN);
    keep = isfinite(F0) & isfinite(Fx) & isfinite(D.V) & isfinite(D.t);
    t_all = [t_all; D.t(keep)]; %#ok<AGROW>
    x_all = [x_all; D.x_rel(keep)]; %#ok<AGROW>
    V_all = [V_all; D.V(keep)]; %#ok<AGROW>
    if isempty(D.W)
        W_all = [W_all; ones(nnz(keep), 1)]; %#ok<AGROW>
    else
        W_all = [W_all; max(D.W(keep), 0.05)]; %#ok<AGROW>
    end
    S_all = [S_all; repmat(sid, nnz(keep), 1)]; %#ok<AGROW>
    F0_all = [F0_all; F0(keep)]; %#ok<AGROW>
    Fx_all = [Fx_all; Fx(keep)]; %#ok<AGROW>
end
end

function beta = solve_first_order_beta_local(f, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors, ridge_factor)
s1 = sin(2*pi*f*t_rel);
c1 = cos(2*pi*f*t_rel);
H_x0 = zeros(numel(t_rel), numel(analysis_sensors));
for is = 1:numel(analysis_sensors)
    H_x0(:, is) = -Fx_all .* (S_all == analysis_sensors(is));
end
H = [H_x0, -Fx_all .* s1, -Fx_all .* c1];
sw = sqrt(W_all);
Hw = H .* sw;
yw = DV .* sw;
G = Hw.' * Hw;
rhs = Hw.' * yw;
beta = (G + ridge_factor * max(trace(G), eps) * eye(size(G))) \ rhs;
end

function rmse = evaluate_first_order_projection_rmse_local(beta, f, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors)
s1 = sin(2*pi*f*t_rel);
c1 = cos(2*pi*f*t_rel);
H_x0 = zeros(numel(t_rel), numel(analysis_sensors));
for is = 1:numel(analysis_sensors)
    H_x0(:, is) = -Fx_all .* (S_all == analysis_sensors(is));
end
H = [H_x0, -Fx_all .* s1, -Fx_all .* c1];
residual = DV - H * beta(:);
rmse = sqrt(sum(W_all .* residual.^2) / sum(W_all));
end

function rmse = evaluate_shifted_template_rmse_local(beta, f, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors)
s1 = sin(2*pi*f*t_rel);
c1 = cos(2*pi*f*t_rel);
x0_sample = zeros(size(t_rel));
for is = 1:numel(analysis_sensors)
    x0_sample(S_all == analysis_sensors(is)) = beta(is);
end
u = x0_sample + beta(end-1) * s1 + beta(end) * c1;
r_all = [];
w_all = [];
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    idx = S_all == sid;
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    F_shift = interp1(Tpl.x_grid, Tpl.v_grid, x_all(idx) - u(idx), 'pchip', NaN);
    keep = isfinite(F_shift);
    if nnz(keep) < 20
        continue;
    end
    H_aff = [ones(nnz(keep), 1), F_shift(keep)];
    ww = W_all(idx);
    VV = V_all(idx);
    beta_aff = (H_aff .* sqrt(ww(keep))) \ (VV(keep) .* sqrt(ww(keep)));
    rr = VV(keep) - H_aff * beta_aff;
    r_all = [r_all; rr]; %#ok<AGROW>
    w_all = [w_all; ww(keep)]; %#ok<AGROW>
end
rmse = sqrt(sum(w_all .* r_all.^2) / sum(w_all));
end

function [rmse, V_pred, sensor_rows] = full_template_objective_local(theta, t_rel, x_all, V_all, W_all, S_all, Template, analysis_sensors)
f = theta(1);
A = theta(2);
phi = theta(3);
d0 = theta(4);
eta = theta(5:end);
u_vib = A * sin(2*pi*f*t_rel + phi);
u_total = d0 * ones(size(t_rel)) + u_vib;
for is = 1:numel(analysis_sensors)
    u_total(S_all == analysis_sensors(is)) = u_total(S_all == analysis_sensors(is)) + eta(is);
end

V_pred = NaN(size(V_all));
r_all = [];
w_all = [];
sensor_rows = table();
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    idx = S_all == sid;
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    F_shift = interp1(Tpl.x_grid, Tpl.v_grid, x_all(idx) - u_total(idx), 'pchip', NaN);
    keep = isfinite(F_shift);
    if nnz(keep) < 20
        continue;
    end
    H_aff = [ones(nnz(keep), 1), F_shift(keep)];
    ww = W_all(idx);
    VV = V_all(idx);
    beta_aff = (H_aff .* sqrt(ww(keep))) \ (VV(keep) .* sqrt(ww(keep)));
    pred_local = H_aff * beta_aff;
    rr = VV(keep) - pred_local;
    local_indices = find(idx);
    V_pred(local_indices(keep)) = pred_local;
    r_all = [r_all; rr]; %#ok<AGROW>
    w_all = [w_all; ww(keep)]; %#ok<AGROW>
    sensor_rows = [sensor_rows; table(sid, beta_aff(2), beta_aff(1), ...
        sqrt(sum(ww(keep).*rr.^2)/sum(ww(keep))), nnz(keep), d0 + eta(is), ...
        'VariableNames', {'sensor_id','gain','bias','rmse','point_count','x0_mm'})]; %#ok<AGROW>
end
rmse = sqrt(sum(w_all .* r_all.^2) / sum(w_all));
end
