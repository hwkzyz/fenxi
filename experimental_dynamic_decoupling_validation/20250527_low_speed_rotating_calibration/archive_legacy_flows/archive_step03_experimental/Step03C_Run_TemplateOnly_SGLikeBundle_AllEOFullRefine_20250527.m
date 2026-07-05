%% Step03B_Run_TemplateOnly_AllEOFullRefine_20250527
% Template-only identification with full waveform refinement for every EO.
%
% Same settings:
%   cfg.target_blades        = 1
%   cfg.analysis_sensors     = [1 6] by default; override with STEP03_ANALYSIS_SENSORS
%   cfg.analysis_start_time  = read from DynamicMap.SourceSettings
%   cfg.target_laps          = 20
%   cfg.analysis_win_size    = 3 laps
%   cfg.sliding_step         = 1 lap
%
% Purpose:
%   Each sliding window evaluates all integer EO candidates with the same
%   full low-speed-template waveform objective. The selected EO/frequency is
%   decided only by the post-refinement weighted RMSE in that window.
%
% Difference from Step03:
%   No reference EO override and no first-three-window EO lock are used.
%   Variable projection is used only to generate an initial frequency,
%   amplitude, phase and offset for each EO candidate.

clear; clc; close all;

%% Settings
route_dir = fileparts(mfilename('fullpath'));
validation_root = fileparts(route_dir);
workspace_root = fileparts(validation_root);

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1, 6];
sensor_override = strtrim(getenv('STEP03_ANALYSIS_SENSORS'));
if ~isempty(sensor_override)
    parsed_sensors = sscanf(sensor_override, '%d').';
    if isempty(parsed_sensors)
        error('STEP03_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysis_sensors = parsed_sensors;
end
cfg.analysis_start_time = NaN;       % overwritten from DynamicMap.SourceSettings
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
cfg.debug_max_windows = parse_positive_integer_env_local('STEP03B_MAX_WINDOWS', inf);
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];

method = struct();
method.ridge_factor = 1e-8;
method.sensor_eta_limit_mm = 0.03;
method.sweep_eta_limits_mm = [0, 0.01, 0.03, 0.06, 0.12];
method.eo_candidates = 3:24;
method.local_scan_half_width_hz = [6, 0.8, 0.18];
method.local_scan_step_hz = [0.5, 0.1, 0.02];
method.full_refine_half_width_hz = 4.0;
method.full_refine_frequency_bound = 'integer_eo_center';
method.selection_rule = 'all_integer_eo_full_waveform_rmse';
method.use_sg_like_dynamic_bundle = true;
method.dynamic_weight_floor = 0.05;
method.template_affine_mode = 'refit_each_theta';

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
result_file = fullfile(result_dir, sprintf('Result_TemplateOnlySGLikeBundleAllEOFullRefine_B%d_%s_20250527.mat', ...
    cfg.target_blades, sensor_tag));

loaded_template = load(template_file, 'Template');
loaded_dynamic = load(dynamic_file, 'DynamicMap');
Template = filter_template_sensors_local(loaded_template.Template, cfg.analysis_sensors, sensor_tag);
DynamicMap = filter_dynamic_map_sensors_local(loaded_dynamic.DynamicMap, cfg.analysis_sensors, sensor_tag);
source_cfg = DynamicMap.SourceSettings;
cfg.analysis_start_time = source_cfg.analysis_start_time;
cfg.target_laps = source_cfg.target_laps;
cfg.analysis_win_size = source_cfg.analysis_win_size;
cfg.sliding_step = source_cfg.sliding_step;

fprintf('\n=== Step03C: template-only SG-like-bundle all-EO full-refine identification ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Template source: %s\n', template_file);
fprintf('Dynamic map source: %s\n', dynamic_file);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('No gap variation is estimated in this route.\n');
fprintf('EO candidates: integer EO %d-%d; every EO enters full waveform refinement.\n', ...
    min(method.eo_candidates), max(method.eo_candidates));
fprintf('Dynamic samples: main-pulse/trust-domain/template-slope weighting enabled.\n');
fprintf('Template affine gain/bias mode: %s\n', method.template_affine_mode);
fprintf('Local frequency scan stages: half-width %s Hz, step %s Hz\n', ...
    mat2str(method.local_scan_half_width_hz), mat2str(method.local_scan_step_hz));
fprintf('Full waveform refinement frequency bound: each EO nominal frequency +/- %.2f Hz\n', ...
    method.full_refine_half_width_hz);

if DynamicMap.TargetBlade ~= cfg.target_blades
    error('DynamicMap target blade (%d) does not match cfg.target_blades (%d).', ...
        DynamicMap.TargetBlade, cfg.target_blades);
end
if ~isequal(DynamicMap.SensorIDs(:).', cfg.analysis_sensors(:).')
    error('DynamicMap sensors %s do not match cfg.analysis_sensors %s.', ...
        mat2str(DynamicMap.SensorIDs), mat2str(cfg.analysis_sensors));
end
if source_cfg.analysis_start_time ~= cfg.analysis_start_time || ...
        source_cfg.target_laps ~= cfg.target_laps || ...
        source_cfg.analysis_win_size ~= cfg.analysis_win_size || ...
        source_cfg.sliding_step ~= cfg.sliding_step
    error('DynamicMap SourceSettings do not match the Step03 cfg.');
end

num_windows = numel(DynamicMap.Window);
if isfinite(cfg.debug_max_windows)
    num_windows = min(num_windows, cfg.debug_max_windows);
end
fprintf('DynamicMap contains %d sliding windows.\n', num_windows);
if isfield(DynamicMap, 'OldSGBestWindowID')
    fprintf('Old SG comparison window: %d, laps %s\n', ...
        DynamicMap.OldSGBestWindowID, mat2str(DynamicMap.OldSGBestLapRange));
end

%% Run proposed identification for every sliding window
WindowResultCell = cell(num_windows, 1);
Trend = table(nan(num_windows,1), nan(num_windows,1), nan(num_windows,1), ...
    nan(num_windows,1), nan(num_windows,1), nan(num_windows,1), ...
    nan(num_windows,1), nan(num_windows,1), nan(num_windows,1), ...
    nan(num_windows,1), nan(num_windows,1), nan(num_windows,1), ...
    nan(num_windows,1), ...
    'VariableNames', {'window_id','lap_start','lap_end','window_center_time', ...
    'rot_freq_hz','rot_rpm','frequency_hz','EO','amplitude_mm','phase_rad', ...
    'mean_x0_mm','weighted_rmse_V','point_count'});

for w_idx = 1:num_windows
    Wmap = DynamicMap.Window(w_idx);
    fprintf('\nWindow %02d/%02d, laps %s, time %.6f-%.6f s\n', ...
        w_idx, num_windows, mat2str(Wmap.lap_range), Wmap.time_window(1), Wmap.time_window(2));

    Rw = identify_one_window_all_eo_local(Wmap, Template, cfg.analysis_sensors, method);
    WindowResultCell{w_idx} = Rw;

    Trend{w_idx, :} = [w_idx, Wmap.lap_range(1), Wmap.lap_range(end), ...
        mean(Wmap.time_window), Wmap.rot_freq_mean_hz, Wmap.rot_rpm_mean, ...
        Rw.frequency_hz, Rw.EO, Rw.amplitude_mm, Rw.phase_rad, ...
        Rw.mean_x0_mm, Rw.weighted_voltage_rmse_V, Rw.point_count];

    fprintf('  all-EO full refine: f %.6f Hz, EO %.4f, A %.6f mm, RMSE %.6f V\n', ...
        Rw.frequency_hz, Rw.EO, Rw.amplitude_mm, Rw.weighted_voltage_rmse_V);
end
WindowResult = vertcat(WindowResultCell{:});

if isfield(DynamicMap, 'OldSGTrends')
    old_tr = DynamicMap.OldSGTrends;
    n_old = min(num_windows, numel(old_tr.Freq));
    Trend.old_sg_frequency_hz = nan(height(Trend), 1);
    Trend.old_sg_amplitude_mm = nan(height(Trend), 1);
    Trend.old_sg_weighted_rmse_V = nan(height(Trend), 1);
    Trend.old_sg_frequency_hz(1:n_old) = old_tr.Freq(1:n_old);
    Trend.old_sg_amplitude_mm(1:n_old) = old_tr.Amp(1:n_old);
    Trend.old_sg_weighted_rmse_V(1:n_old) = old_tr.WeightedRMSE(1:n_old);
end

Result = struct();
Result.Method = 'template_only_sg_like_bundle_all_integer_eo_full_waveform_refine';
Result.AnalysisSettings = cfg;
Result.MethodSettings = method;
Result.DynamicMapFile = dynamic_file;
Result.TargetBlade = cfg.target_blades;
Result.SensorIDs = cfg.analysis_sensors;
Result.SensorTag = sensor_tag;
Result.GlobalEO = mode(round(Trend.EO));
Result.InitialEO = NaN;
Result.GlobalEOSelectionTable = build_eo_summary_from_window_candidates_local(WindowResult, method);
Result.EOLock.Mode = 'none';
Result.EOLock.Source = 'all EO candidates refined in each window; selection by final waveform RMSE';
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
fprintf('No EO lock or reference EO was used; each 3-rev window selected by final waveform RMSE.\n');
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
fig1 = figure('Name', 'Step03C SG-like bundle all-EO full-refine trends', 'Color', 'w', ...
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

exportgraphics(fig1, fullfile(figure_dir, sprintf('Step03C_SGLikeBundleAllEOFullRefineTrends_B%d_%s.png', ...
    cfg.target_blades, sensor_tag)), 'Resolution', 300);

%% Local calculation blocks
function value = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    value = default_value;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 1 || fix(value) ~= value
    error('%s must be a positive integer.', name);
end
end

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

function R = identify_one_window_all_eo_local(Wmap, Template, analysis_sensors, method)
[t_all, x_all, V_all, W_all, S_all, F0_all, Fx_all] = pack_window_data_local(Wmap, Template, analysis_sensors, method);
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
SensorData = build_sensor_data_local(x_all, V_all, W_all, S_all, Template, analysis_sensors, method);

eo_candidates = method.eo_candidates(:);
candidate = repmat(struct( ...
    'EO_integer', NaN, 'frequency_hz', NaN, 'EO_refined', NaN, ...
    'amplitude_mm', NaN, 'phase_rad', NaN, 'mean_x0_mm', NaN, ...
    'weighted_rmse_V', inf, 'first_order_frequency_hz', NaN, ...
    'first_order_rmse_V', inf, 'theta_best', [], 'V_pred', [], ...
    'sensor_fit', table(), 'freq_grid_fine', [], 'rmse_fine', []), numel(eo_candidates), 1);

for ie = 1:numel(eo_candidates)
    eo = eo_candidates(ie);
    f_center = eo * Wmap.rot_freq_mean_hz;
    [freq_grid_fine, rmse_fine, ~, best_f, best_beta] = adaptive_first_order_scan_local( ...
        f_center, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors, method, SensorData);

    best_x0 = best_beta(1:numel(analysis_sensors));
    best_d0 = mean(best_x0, 'omitnan');
    best_A = hypot(best_beta(end-1), best_beta(end));
    best_phi = atan2(best_beta(end), best_beta(end-1));

    [theta_best, best_rmse, best_pred, best_sensor_fit] = refine_full_template_local( ...
        best_f, best_A, best_phi, best_d0, t0, t_rel, SensorData, analysis_sensors, method, [], f_center);

    candidate(ie).EO_integer = eo;
    candidate(ie).frequency_hz = theta_best(1);
    candidate(ie).EO_refined = theta_best(1) / Wmap.rot_freq_mean_hz;
    candidate(ie).amplitude_mm = theta_best(2);
    candidate(ie).phase_rad = theta_best(3);
    candidate(ie).mean_x0_mm = theta_best(4);
    candidate(ie).weighted_rmse_V = best_rmse;
    candidate(ie).first_order_frequency_hz = best_f;
    candidate(ie).first_order_rmse_V = min(rmse_fine);
    candidate(ie).theta_best = theta_best;
    candidate(ie).V_pred = best_pred;
    candidate(ie).sensor_fit = best_sensor_fit;
    candidate(ie).freq_grid_fine = freq_grid_fine;
    candidate(ie).rmse_fine = rmse_fine;
end

candidate_table = struct2table(rmfield(candidate, {'theta_best','V_pred','sensor_fit','freq_grid_fine','rmse_fine'}));
candidate_table = sortrows(candidate_table, {'weighted_rmse_V','first_order_rmse_V','EO_integer'}, {'ascend','ascend','ascend'});
[~, best_original_idx] = min([candidate.weighted_rmse_V]);
best = candidate(best_original_idx);
theta_best = best.theta_best;
best_f = theta_best(1);
best_A = theta_best(2);
best_phi = theta_best(3);
best_d0 = theta_best(4);
best_eta = theta_best(5:end);
best_x0 = best_d0 + best_eta;
best_u_vib = best_A * sin(2*pi*best_f*t_rel + best_phi);
best_u = best_d0 * ones(size(t_all)) + best_u_vib;
for js = 1:numel(analysis_sensors)
    best_u(S_all == analysis_sensors(js)) = best_u(S_all == analysis_sensors(js)) + best_eta(js);
end

best_sensor_fit = best.sensor_fit;
for is = 1:numel(analysis_sensors)
    best_sensor_fit.x0_mm(best_sensor_fit.sensor_id == analysis_sensors(is)) = best_x0(is);
end

R = struct();
R.window_id = Wmap.window_id;
R.lap_range = Wmap.lap_range;
R.time_window = Wmap.time_window;
R.frequency_hz = best_f;
R.EO = best_f / Wmap.rot_freq_mean_hz;
R.EO_integer = best.EO_integer;
R.amplitude_mm = best_A;
R.phase_rad = best_phi;
R.mean_x0_mm = best_d0;
R.sensor_x0_mm = best_x0(:).';
R.weighted_voltage_rmse_V = best.weighted_rmse_V;
R.sensor_fit = best_sensor_fit;
R.point_count = numel(t_all);
R.t = t_all;
R.t0 = t0;
R.t_rel = t_rel;
R.x_rel = x_all;
R.V_meas = V_all;
R.V_pred = best.V_pred;
R.u_total = best_u;
R.u_vibration = best_u_vib;
R.sensor_id = S_all;
R.global_eo = best.EO_integer;
R.theta_best = theta_best;
R.CandidateTable = candidate_table;
R.Candidates = candidate;
R.freq_grid_coarse = [];
R.rmse_coarse = [];
R.freq_grid_fine = best.freq_grid_fine;
R.rmse_fine = best.rmse_fine;
end

function eo_table = build_eo_summary_from_window_candidates_local(WindowResult, method)
eo_candidates = method.eo_candidates(:);
selected_eo = round([WindowResult.EO_integer]).';
win_count = numel(WindowResult);
selected_count = zeros(numel(eo_candidates), 1);
median_rank = nan(numel(eo_candidates), 1);
median_rmse = nan(numel(eo_candidates), 1);
mean_rmse = nan(numel(eo_candidates), 1);
for ie = 1:numel(eo_candidates)
    eo = eo_candidates(ie);
    selected_count(ie) = sum(selected_eo == eo);
    ranks = nan(win_count, 1);
    rmse = nan(win_count, 1);
    for iw = 1:win_count
        T = WindowResult(iw).CandidateTable;
        idx = find(T.EO_integer == eo, 1, 'first');
        if ~isempty(idx)
            ranks(iw) = idx;
            rmse(iw) = T.weighted_rmse_V(idx);
        end
    end
    median_rank(ie) = median(ranks, 'omitnan');
    median_rmse(ie) = median(rmse, 'omitnan');
    mean_rmse(ie) = mean(rmse, 'omitnan');
end
eo_table = table(eo_candidates, selected_count, median_rank, median_rmse, mean_rmse, ...
    'VariableNames', {'EO','selected_count','median_final_rank','median_weighted_rmse_V','mean_weighted_rmse_V'});
eo_table = sortrows(eo_table, {'selected_count','median_final_rank','median_weighted_rmse_V'}, ...
    {'descend','ascend','ascend'});
end

function [freq_trace, rmse_trace, beta_trace, best_f, best_beta] = adaptive_first_order_scan_local( ...
    f_center, t_rel, DV, Fx_all, W_all, S_all, analysis_sensors, method, SensorData)
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
        stage_rmse(jf) = evaluate_shifted_template_rmse_fast_local(beta, f, t_rel, SensorData, method);
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

function [t_all, x_all, V_all, W_all, S_all, F0_all, Fx_all] = pack_window_data_local(Wmap, Template, analysis_sensors, method)
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
    if isfield(method, 'use_sg_like_dynamic_bundle') && method.use_sg_like_dynamic_bundle
        keep = keep & build_template_dynamic_keep_mask_local(D, Tpl);
    end
    t_all = [t_all; D.t(keep)]; %#ok<AGROW>
    x_all = [x_all; D.x_rel(keep)]; %#ok<AGROW>
    V_all = [V_all; D.V(keep)]; %#ok<AGROW>
    W_dynamic = build_template_dynamic_weight_local(D, Tpl, method);
    W_all = [W_all; W_dynamic(keep)]; %#ok<AGROW>
    S_all = [S_all; repmat(sid, nnz(keep), 1)]; %#ok<AGROW>
    F0_all = [F0_all; F0(keep)]; %#ok<AGROW>
    Fx_all = [Fx_all; Fx(keep)]; %#ok<AGROW>
end
end

function keep = build_template_dynamic_keep_mask_local(D, Tpl)
v = D.V(:);
x = D.x_rel(:);
keep_pulse = isolate_main_pulse_local(v, Tpl.threshold);
domain_margin = 0.02;
keep_domain = x >= Tpl.x_domain(1) + domain_margin & x <= Tpl.x_domain(2) - domain_margin;
keep = keep_pulse & keep_domain;
if nnz(keep) < 20
    keep = keep_domain;
end
end

function W = build_template_dynamic_weight_local(D, Tpl, method)
floor_w = method.dynamic_weight_floor;
v = D.V(:);
x = D.x_rel(:);
w_edge = build_edge_weight_local(D.t(:), v, floor_w);
slope = abs(interp1(Tpl.x_grid, Tpl.dv_dx, x, 'linear', 0));
if max(slope) > 0
    w_slope = floor_w + (1 - floor_w) * slope ./ max(slope);
else
    w_slope = ones(size(slope));
end
if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
    w_template = interp1(Tpl.x_grid, Tpl.bin_weight, x, 'linear', floor_w);
else
    w_template = ones(size(x));
end
W = max(floor_w, w_edge .* w_slope .* w_template);
if max(W) > 0
    W = max(floor_w, W ./ max(W));
end
end

function mask = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
[~, idx_peak] = max(v);
best_seg = 1;
for i = 1:numel(starts)
    seg = idx(starts(i):ends(i));
    if any(seg == idx_peak)
        best_seg = i;
        break;
    end
end
seg = idx(starts(best_seg):ends(best_seg));
mask(seg) = true;
end

function w_edge = build_edge_weight_local(t, v, floor_w)
if numel(v) < 3 || range(t) <= 0
    w_edge = ones(size(v));
    return;
end
span = min(11, 2 * floor((numel(v) - 1) / 2) + 1);
if span >= 3
    v_smooth = movmedian(v(:), span);
else
    v_smooth = v(:);
end
dv = abs(gradient(v_smooth, t(:)));
if max(dv) > 0
    w_edge = floor_w + (1 - floor_w) * dv ./ max(dv);
else
    w_edge = ones(size(dv));
end
w_edge = max(floor_w, w_edge);
end

function SensorData = build_sensor_data_local(x_all, V_all, W_all, S_all, Template, analysis_sensors, method)
SensorData = repmat(struct( ...
    'sensor_id', NaN, 'idx', [], 'x', [], 'V', [], 'W', [], 'sqrtW', [], ...
    'interp_v', [], 'local_indices', [], 'gain_ref', NaN, 'bias_ref', NaN), numel(analysis_sensors), 1);
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    idx = S_all == sid;
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    SensorData(is).sensor_id = sid;
    SensorData(is).idx = idx;
    SensorData(is).x = x_all(idx);
    SensorData(is).V = V_all(idx);
    SensorData(is).W = W_all(idx);
    SensorData(is).sqrtW = sqrt(W_all(idx));
    SensorData(is).interp_v = griddedInterpolant(Tpl.x_grid(:), Tpl.v_grid(:), 'pchip', 'none');
    SensorData(is).local_indices = find(idx);
    F0 = SensorData(is).interp_v(SensorData(is).x);
    keep = isfinite(F0);
    if nnz(keep) >= 20
        H_aff = [ones(nnz(keep), 1), F0(keep)];
        sw = SensorData(is).sqrtW(keep);
        beta_aff = (H_aff .* sw) \ (SensorData(is).V(keep) .* sw);
        SensorData(is).bias_ref = beta_aff(1);
        SensorData(is).gain_ref = beta_aff(2);
    else
        SensorData(is).bias_ref = 0;
        SensorData(is).gain_ref = 1;
    end
    if ~isfinite(SensorData(is).gain_ref) || abs(SensorData(is).gain_ref) < 1e-6
        SensorData(is).gain_ref = 1;
    end
    if ~isfinite(SensorData(is).bias_ref)
        SensorData(is).bias_ref = 0;
    end
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

function rmse = evaluate_shifted_template_rmse_fast_local(beta, f, t_rel, SensorData, method)
s1 = sin(2*pi*f*t_rel);
c1 = cos(2*pi*f*t_rel);
x0_sample = zeros(size(t_rel));
for is = 1:numel(SensorData)
    x0_sample(SensorData(is).idx) = beta(is);
end
u = x0_sample + beta(end-1) * s1 + beta(end) * c1;
weighted_sse = 0;
weight_sum = 0;
for is = 1:numel(SensorData)
    SD = SensorData(is);
    F_shift = SD.interp_v(SD.x - u(SD.idx));
    keep = isfinite(F_shift);
    if nnz(keep) < 20
        continue;
    end
    if strcmpi(method.template_affine_mode, 'fixed_unshifted_per_window')
        pred_local = SD.bias_ref + SD.gain_ref * F_shift(keep);
    else
        H_aff = [ones(nnz(keep), 1), F_shift(keep)];
        sw = SD.sqrtW(keep);
        beta_aff = (H_aff .* sw) \ (SD.V(keep) .* sw);
        pred_local = H_aff * beta_aff;
    end
    rr = SD.V(keep) - pred_local;
    weighted_sse = weighted_sse + sum(SD.W(keep) .* rr.^2);
    weight_sum = weight_sum + sum(SD.W(keep));
end
rmse = sqrt(weighted_sse / weight_sum);
end

function [theta_best, best_rmse, best_pred, best_sensor_fit] = refine_full_template_local( ...
    f0, A0, phi0, d00, t0, t_rel, SensorData, analysis_sensors, method, previous_solution, f_bound_center)
eta_limit = method.sensor_eta_limit_mm;
f_half_width = method.full_refine_half_width_hz;
if nargin < 11 || isempty(f_bound_center) || ~isfinite(f_bound_center)
    f_bound_center = f0;
end
f0 = min(max(f0, f_bound_center - f_half_width), f_bound_center + f_half_width);
theta_best = [f0, max(A0, 0.05), phi0, d00, zeros(1, numel(analysis_sensors))];
theta_lb = [f_bound_center - f_half_width, 0.00, -pi, -0.60, -eta_limit * ones(1, numel(analysis_sensors))];
theta_ub = [f_bound_center + f_half_width, 0.80,  pi,  0.60,  eta_limit * ones(1, numel(analysis_sensors))];

eta_starts = zeros(1, numel(analysis_sensors));
if isempty(previous_solution)
    phi_starts = linspace(-pi, pi, 17);
    A_starts = [0.10, 0.18, 0.25, 0.32, 0.40, 0.50];
    d0_starts = [d00 - 0.10, d00, d00 + 0.10];
else
    prev_theta = previous_solution.theta;
    prev_phi = atan2(sin(prev_theta(3) + 2*pi*prev_theta(1)*(t0 - previous_solution.t0)), ...
        cos(prev_theta(3) + 2*pi*prev_theta(1)*(t0 - previous_solution.t0)));
    phi_starts = atan2(sin([phi0, phi0 - pi/4, phi0 + pi/4, ...
        prev_phi, prev_phi - pi/6, prev_phi + pi/6]), ...
        cos([phi0, phi0 - pi/4, phi0 + pi/4, prev_phi, prev_phi - pi/6, prev_phi + pi/6]));
    A_starts = [max(0.05, A0), max(0.05, 0.90*prev_theta(2)), ...
        max(0.05, prev_theta(2)), min(0.80, 1.10*prev_theta(2)), 0.38];
    A_starts = unique(round(A_starts, 6), 'stable');
    d0_starts = [d00, prev_theta(4), 0.5*(d00 + prev_theta(4))];
    if numel(prev_theta) >= 4 + numel(analysis_sensors)
        eta_starts = [eta_starts; prev_theta(5:end)];
    end
end

best_rmse = inf;
for ia = 1:numel(A_starts)
    for ip = 1:numel(phi_starts)
        for id0 = 1:numel(d0_starts)
            for ie = 1:size(eta_starts, 1)
                theta_try = [f0, A_starts(ia), phi_starts(ip), d0_starts(id0), eta_starts(ie, :)];
                theta_try = min(max(theta_try, theta_lb), theta_ub);
                rmse_try = full_template_rmse_fast_local(theta_try, t_rel, SensorData, analysis_sensors, method);
                if rmse_try < best_rmse
                    best_rmse = rmse_try;
                    theta_best = theta_try;
                end
            end
        end
    end
end

step = [0.25, 0.04, 0.25, 0.04, 0.03 * ones(1, numel(analysis_sensors))];
for level = 1:7
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
            rmse_try = full_template_rmse_fast_local(theta_try, t_rel, SensorData, analysis_sensors, method);
            if rmse_try + 1e-10 < best_rmse
                best_rmse = rmse_try;
                theta_best = theta_try;
                improved = true;
            end
        end
    end
    step = step * 0.5;
end
[~, best_pred, best_sensor_fit] = full_template_objective_local(theta_best, t_rel, SensorData, analysis_sensors, method);
end

function rmse = full_template_rmse_fast_local(theta, t_rel, SensorData, analysis_sensors, method)
f = theta(1);
A = theta(2);
phi = theta(3);
d0 = theta(4);
eta = theta(5:end);
u_vib = A * sin(2*pi*f*t_rel + phi);
u_total = d0 * ones(size(t_rel)) + u_vib;
for is = 1:numel(analysis_sensors)
    u_total(SensorData(is).idx) = u_total(SensorData(is).idx) + eta(is);
end

weighted_sse = 0;
weight_sum = 0;
for is = 1:numel(analysis_sensors)
    SD = SensorData(is);
    F_shift = SD.interp_v(SD.x - u_total(SD.idx));
    keep = isfinite(F_shift);
    if nnz(keep) < 20
        continue;
    end
    if strcmpi(method.template_affine_mode, 'fixed_unshifted_per_window')
        pred_local = SD.bias_ref + SD.gain_ref * F_shift(keep);
    else
        H_aff = [ones(nnz(keep), 1), F_shift(keep)];
        sw = SD.sqrtW(keep);
        beta_aff = (H_aff .* sw) \ (SD.V(keep) .* sw);
        pred_local = H_aff * beta_aff;
    end
    rr = SD.V(keep) - pred_local;
    weighted_sse = weighted_sse + sum(SD.W(keep) .* rr.^2);
    weight_sum = weight_sum + sum(SD.W(keep));
end
rmse = sqrt(weighted_sse / weight_sum);
end

function [rmse, V_pred, sensor_rows] = full_template_objective_local(theta, t_rel, SensorData, analysis_sensors, method)
f = theta(1);
A = theta(2);
phi = theta(3);
d0 = theta(4);
eta = theta(5:end);
u_vib = A * sin(2*pi*f*t_rel + phi);
u_total = d0 * ones(size(t_rel)) + u_vib;
for is = 1:numel(analysis_sensors)
    u_total(SensorData(is).idx) = u_total(SensorData(is).idx) + eta(is);
end

total_points = sum(arrayfun(@(s) numel(s.V), SensorData));
V_pred = NaN(total_points, 1);
weighted_sse = 0;
weight_sum = 0;
rows = repmat(struct('sensor_id', NaN, 'gain', NaN, 'bias', NaN, ...
    'rmse', NaN, 'point_count', NaN, 'x0_mm', NaN), numel(analysis_sensors), 1);
for is = 1:numel(analysis_sensors)
    SD = SensorData(is);
    F_shift = SD.interp_v(SD.x - u_total(SD.idx));
    keep = isfinite(F_shift);
    if nnz(keep) < 20
        continue;
    end
    if strcmpi(method.template_affine_mode, 'fixed_unshifted_per_window')
        pred_local = SD.bias_ref + SD.gain_ref * F_shift(keep);
        gain_local = SD.gain_ref;
        bias_local = SD.bias_ref;
    else
        H_aff = [ones(nnz(keep), 1), F_shift(keep)];
        sw = SD.sqrtW(keep);
        beta_aff = (H_aff .* sw) \ (SD.V(keep) .* sw);
        pred_local = H_aff * beta_aff;
        gain_local = beta_aff(2);
        bias_local = beta_aff(1);
    end
    rr = SD.V(keep) - pred_local;
    V_pred(SD.local_indices(keep)) = pred_local;
    local_sse = sum(SD.W(keep) .* rr.^2);
    local_weight = sum(SD.W(keep));
    weighted_sse = weighted_sse + local_sse;
    weight_sum = weight_sum + local_weight;
    rows(is).sensor_id = SD.sensor_id;
    rows(is).gain = gain_local;
    rows(is).bias = bias_local;
    rows(is).rmse = sqrt(local_sse / local_weight);
    rows(is).point_count = nnz(keep);
    rows(is).x0_mm = d0 + eta(is);
end
rmse = sqrt(weighted_sse / weight_sum);
sensor_rows = struct2table(rows);
sensor_rows = sensor_rows(isfinite(sensor_rows.sensor_id), :);
end

