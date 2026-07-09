function Summary = Diagnose_Step05_SGWindow_CurrentTemplateInversion_20251222()
%DIAGNOSE_STEP05_SGWINDOW_CURRENTTEMPLATEINVERSION_20251222
% Differential diagnostic for the 20251222 Step04/Step05 failure.
%
% This script deliberately uses the old SG-success Step5 dynamic observation
% bundle as the common high-speed data source, then replaces only the static
% calibration model:
%   1) old SG calibration, as a positive control;
%   2) current Step04 non-parametric template, branch-monotone inverse;
%   3) current Step04 non-parametric template, nearest-crossing inverse.
%
% If the old SG control ranks EO14 but the current Step04 inverse does not,
% the remaining error is in the static template / inversion contract rather
% than in high-speed OPR row selection.

cfg = BTTProjectConfig_20251222();

sg_root = resolve_sg_success_root_local(fileparts(mfilename('fullpath')));
if isempty(sg_root) || ~isfolder(sg_root)
    error(['Cannot locate the SG-success route. Set SG_SUCCESS_ROOT to the ' ...
        'folder containing SGCalib_20251222_B1_S1.mat.']);
end
sg_path_cleanup = add_sg_paths_for_whole_diagnostic_local(sg_root); %#ok<NASGU>

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_sgwindow_current_template_inversion_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

blade_id = read_env_number_local('STEP05_DIAG_TARGET_BLADE', 1);
sensors = read_env_vector_local('STEP05_DIAG_SENSORS', cfg.sensor_ids);
start_time = read_env_number_local('STEP05_DIAG_START_TIME', 50.0);
target_laps = read_env_number_local('STEP05_DIAG_TARGET_LAPS', 20);
window_laps = read_env_number_local('STEP05_DIAG_WINDOW_LAPS', 3);
max_windows = read_env_number_local('STEP05_DIAG_MAX_WINDOWS', inf);

Template = load_current_step04_template_local(cfg);
[experiment_case, SGCalib, old_cfg] = build_sg_success_case_local( ...
    sg_root, blade_id, sensors, start_time, target_laps, window_laps);

num_windows = floor((old_cfg.target_laps - old_cfg.analysis_win_size) / ...
    old_cfg.sliding_step) + 1;
num_windows = min(num_windows, max_windows);
if ~isfinite(num_windows)
    num_windows = floor((old_cfg.target_laps - old_cfg.analysis_win_size) / ...
        old_cfg.sliding_step) + 1;
end

all_candidates = table();
summary_rows = repmat(make_empty_summary_row_local(), num_windows, 1);

for iw = 1:num_windows
    lap_start = 1 + (iw - 1) * old_cfg.sliding_step;
    lap_end = lap_start + old_cfg.analysis_win_size - 1;
    lap_range = lap_start:lap_end;

    bundle = build_experiment_observation_bundle_single_sync_20251222( ...
        experiment_case, SGCalib, lap_range, sensors);
    eo_candidates = build_eo_candidates_local( ...
        bundle.rot_freq_mean_hz, old_cfg.freq_search_hz, old_cfg.eo_pad);

    Tsg = scan_sg_calibration_eo_local(bundle, SGCalib, eo_candidates, old_cfg);
    Tmono = scan_current_template_eo_local( ...
        bundle, Template, eo_candidates, sensors, blade_id, ...
        old_cfg.weight_floor, 'branch_monotone_inverse');
    Tcross = scan_current_template_eo_local( ...
        bundle, Template, eo_candidates, sensors, blade_id, ...
        old_cfg.weight_floor, 'nearest_crossing_inverse');

    Tsg.model = repmat("sg_success_control", height(Tsg), 1);
    Tmono.model = repmat("current_step04_branch_monotone_inverse", height(Tmono), 1);
    Tcross.model = repmat("current_step04_nearest_crossing_inverse", height(Tcross), 1);

    T = [Tsg; Tmono; Tcross];
    T.window_id = repmat(iw, height(T), 1);
    T.lap_start = repmat(lap_start, height(T), 1);
    T.lap_end = repmat(lap_end, height(T), 1);
    T = movevars(T, {'window_id', 'lap_start', 'lap_end', 'model'}, 'Before', 1);
    all_candidates = [all_candidates; T]; %#ok<AGROW>

    summary_rows(iw) = summarize_window_local(iw, lap_start, lap_end, bundle, T);
end

CandidateTable = all_candidates;
WindowSummary = struct2table(summary_rows);
TemplateDiagnostics = diagnose_current_template_branches_local(Template, sensors, blade_id);

candidate_csv = fullfile(out_dir, ...
    sprintf('SGWindow_CurrentTemplate_Candidates_B%d_S%s_20251222.csv', ...
    blade_id, join_sensor_tag_local(sensors)));
summary_csv = fullfile(out_dir, ...
    sprintf('SGWindow_CurrentTemplate_WindowSummary_B%d_S%s_20251222.csv', ...
    blade_id, join_sensor_tag_local(sensors)));
template_csv = fullfile(out_dir, ...
    sprintf('SGWindow_CurrentTemplate_BranchDiagnostics_B%d_S%s_20251222.csv', ...
    blade_id, join_sensor_tag_local(sensors)));
mat_file = fullfile(out_dir, ...
    sprintf('SGWindow_CurrentTemplate_Diagnostic_B%d_S%s_20251222.mat', ...
    blade_id, join_sensor_tag_local(sensors)));

writetable(CandidateTable, candidate_csv);
writetable(WindowSummary, summary_csv);
writetable(TemplateDiagnostics, template_csv);

Summary = struct();
Summary.SGRoot = sg_root;
Summary.OutputDir = out_dir;
Summary.WindowSummary = WindowSummary;
Summary.CandidateTable = CandidateTable;
Summary.TemplateDiagnostics = TemplateDiagnostics;
save(mat_file, 'Summary', '-v7.3');

fprintf('\n=== Step05 SG-window/current-template inversion diagnostic ===\n');
fprintf('SG route:\n  %s\n', sg_root);
fprintf('Output:\n  %s\n', out_dir);
fprintf('Windows: %d, B%d, S%s, start %.3f s\n', ...
    height(WindowSummary), blade_id, join_sensor_tag_local(sensors), start_time);
disp(WindowSummary);
print_sequence_local(WindowSummary);
end


function [experiment_case, Calib, cfg] = build_sg_success_case_local( ...
    sg_root, blade_id, sensors, start_time, target_laps, window_laps)
cfg = build_single_sync_experiment_config_20251222();
cfg.target_blade = blade_id;
cfg.target_blades = blade_id;
cfg.analysis_sensors = sensors(:).';
cfg.analysis_start_time = start_time;
cfg.target_laps = target_laps;
cfg.analysis_win_size = window_laps;
cfg.sliding_step = 1;
cfg.show_figures = false;
cfg.step5_plot_case = false;
cfg.step5_plot_speed = false;
cfg.step5_plot_window_trends = false;
cfg.step5_plot_best_window = false;
cfg.step5_plot_seed_scan = false;
cfg.step5_save_figures = false;
cfg.save_experiment_case = false;
cfg.step5_case_mode = 'memory';

experiment_case = prepare_identification_case_single_sync_20251222( ...
    cfg, false, false);

sensor_structs = cell(numel(sensors), 1);
for i = 1:numel(sensors)
    sid = sensors(i);
    f = fullfile(sg_root, sprintf('SGCalib_20251222_B%d_S%d.mat', blade_id, sid));
    if ~isfile(f)
        error('Missing SG calibration file: %s', f);
    end
    loaded = load(f, 'Calib');
    sensor_structs{i} = loaded.Calib.Sensor;
end


function cleanup_obj = add_sg_paths_local(sg_root)
addpath(sg_root, '-begin');
old_scripts = fullfile(sg_root, '旧文件夹', 'old_scripts');
if isfolder(old_scripts)
    addpath(old_scripts, '-begin');
    cleanup_obj = onCleanup(@() cleanup_sg_paths_local(sg_root, old_scripts));
else
    cleanup_obj = onCleanup(@() rmpath(sg_root));
end
end


function cleanup_sg_paths_local(sg_root, old_scripts)
if isfolder(old_scripts)
    rmpath(old_scripts);
end
rmpath(sg_root);
end

Calib = struct();
Calib.Sensor = vertcat(sensor_structs{:});
Calib.SensorIDs = sensors(:).';
if isfield(loaded.Calib, 'ReferenceFrame')
    Calib.ReferenceFrame = loaded.Calib.ReferenceFrame;
else
    Calib.ReferenceFrame = struct();
end
Calib.ReferenceFrame.sensor_ids = sensors(:).';
Calib.ReferenceFrame.standard_relative_angles = [Calib.Sensor.alpha_ref];
end


function T = scan_sg_calibration_eo_local(bundle, Calib, eo_candidates, cfg)
sensors = [Calib.Sensor.sensor_id];
x_stat = nan(size(bundle.X(:)));
W = bundle.W(:);
Vpred_fun = @(yfit) evaluate_sg_voltage_rmse_local(bundle, Calib, yfit, W);

for i = 1:numel(sensors)
    sid = sensors(i);
    C = Calib.Sensor(i);
    idx = bundle.S(:) == sid;
    if ~any(idx)
        continue;
    end
    x_centered = bundle.X(idx) - C.xc;
    v_above = max(bundle.V(idx) - C.baseline, 1e-6);
    v_clip = min(max(v_above ./ max(C.B, 1e-6), 1e-6), 0.999999);
    sgn = sign(x_centered);
    sgn(sgn == 0) = 1;
    x_stat(idx) = C.xc + sgn .* C.w .* (-log(v_clip)).^(1 ./ C.n);
end

y_obs = bundle.X(:) - x_stat(:);
T = scan_y_observation_eo_local(bundle.Theta(:), y_obs, W, eo_candidates, Vpred_fun);
end


function T = scan_current_template_eo_local( ...
    bundle, Template, eo_candidates, sensors, blade_id, weight_floor, inverse_mode)

Xabs = bundle.X(:);
V = bundle.V(:);
sensor_id = bundle.S(:);
theta = bundle.Theta(:);
W = bundle.W(:);
W(~isfinite(W)) = weight_floor;
W = max(W, weight_floor);

x_rel = nan(size(Xabs));
x_static_rel = nan(size(Xabs));

for sid = sensors(:).'
    idx = sensor_id == sid;
    if ~any(idx)
        continue;
    end
    tpl = get_current_template_entry_local(Template, sid, blade_id);
    if isempty(tpl)
        continue;
    end
    xc = get_template_xc_local(tpl);
    xr = Xabs(idx) - xc;
    x_rel(idx) = xr;
    switch lower(strtrim(inverse_mode))
        case 'branch_monotone_inverse'
            x_static_rel(idx) = invert_template_branch_monotone_local(tpl, V(idx), xr);
        case 'nearest_crossing_inverse'
            x_static_rel(idx) = invert_template_nearest_crossing_local(tpl, V(idx), xr);
        otherwise
            error('Unsupported inverse mode: %s', inverse_mode);
    end
end

y_obs = x_rel - x_static_rel;
Vpred_fun = @(yfit) evaluate_current_template_voltage_rmse_local( ...
    Template, sensors, blade_id, Xabs, sensor_id, V, yfit, W);
T = scan_y_observation_eo_local(theta, y_obs, W, eo_candidates, Vpred_fun);
end


function T = scan_y_observation_eo_local(theta, y_obs, W, eo_candidates, Vpred_fun)
rows = repmat(struct( ...
    'EO', NaN, 'rank', NaN, 'A_mm', NaN, 'phi_rad', NaN, 'd0_mm', NaN, ...
    'weighted_y_rmse_mm', inf, 'weighted_voltage_rmse_v', inf, ...
    'plain_voltage_rmse_v', inf, 'point_count', NaN, ...
    'valid_query_fraction', NaN), numel(eo_candidates), 1);

base_valid = isfinite(theta) & isfinite(y_obs) & isfinite(W) & W > 0;
for i = 1:numel(eo_candidates)
    EO = eo_candidates(i);
    X = [sin(EO .* theta), cos(EO .* theta), ones(size(theta))];
    valid = base_valid & all(isfinite(X), 2);
    if nnz(valid) < max(20, size(X, 2) + 2)
        continue;
    end
    wg = sqrt(W(valid));
    beta = (X(valid, :) .* wg) \ (y_obs(valid) .* wg);
    if numel(beta) < 3 || any(~isfinite(beta))
        continue;
    end
    y_fit = X * beta;
    y_res = y_obs(valid) - X(valid, :) * beta;
    y_rmse = sqrt(sum(W(valid) .* y_res.^2) ./ max(sum(W(valid)), eps));
    [wrmse, prmse, valid_fraction] = Vpred_fun(y_fit);

    rows(i).EO = EO;
    rows(i).A_mm = hypot(beta(1), beta(2));
    rows(i).phi_rad = wrap_to_pi_local(atan2(beta(2), beta(1)));
    rows(i).d0_mm = beta(3);
    rows(i).weighted_y_rmse_mm = y_rmse;
    rows(i).weighted_voltage_rmse_v = wrmse;
    rows(i).plain_voltage_rmse_v = prmse;
    rows(i).point_count = nnz(valid);
    rows(i).valid_query_fraction = valid_fraction;
end

T = struct2table(rows);
T = T(isfinite(T.EO) & isfinite(T.weighted_voltage_rmse_v), :);
if isempty(T)
    return;
end
T = sortrows(T, {'weighted_voltage_rmse_v', 'weighted_y_rmse_mm', 'EO'}, ...
    {'ascend', 'ascend', 'ascend'});
T.rank = (1:height(T)).';
end


function [wrmse, prmse, valid_fraction] = evaluate_sg_voltage_rmse_local(bundle, Calib, yfit, W)
Vpred = nan(size(bundle.V(:)));
for i = 1:numel(Calib.Sensor)
    C = Calib.Sensor(i);
    idx = bundle.S(:) == C.sensor_id;
    if ~any(idx)
        continue;
    end
    xq = bundle.X(idx) - yfit(idx);
    Vpred(idx) = C.B .* exp(-abs((xq - C.xc) ./ C.w).^C.n) + C.baseline;
end
[wrmse, prmse, valid_fraction] = weighted_voltage_rmse_local(bundle.V(:), Vpred, W);
end


function [wrmse, prmse, valid_fraction] = evaluate_current_template_voltage_rmse_local( ...
    Template, sensors, blade_id, Xabs, sensor_id, V, yfit, W)
Vpred = nan(size(V(:)));
for sid = sensors(:).'
    idx = sensor_id == sid;
    if ~any(idx)
        continue;
    end
    tpl = get_current_template_entry_local(Template, sid, blade_id);
    if isempty(tpl)
        continue;
    end
    xq = Xabs(idx) - get_template_xc_local(tpl) - yfit(idx);
    Vpred(idx) = interp1(tpl.x_grid(:), tpl.v_grid(:), xq, 'pchip', NaN);
end
[wrmse, prmse, valid_fraction] = weighted_voltage_rmse_local(V, Vpred, W);
end


function [wrmse, prmse, valid_fraction] = weighted_voltage_rmse_local(V, Vpred, W)
valid = isfinite(V) & isfinite(Vpred) & isfinite(W) & W > 0;
valid_fraction = nnz(valid) / max(numel(V), 1);
if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
    return;
end
res = V(:) - Vpred(:);
wrmse = sqrt(sum(W(valid) .* res(valid).^2) ./ max(sum(W(valid)), eps));
prmse = sqrt(mean(res(valid).^2));
end


function x_static = invert_template_branch_monotone_local(tpl, v_obs, x_ref)
x_static = nan(size(v_obs(:)));
side = sign(x_ref(:));
side(side == 0) = 1;
left = side < 0;
right = side >= 0;
x_static(left) = inverse_one_side_monotone_local(tpl, v_obs(left), -1);
x_static(right) = inverse_one_side_monotone_local(tpl, v_obs(right), 1);
end


function xq = inverse_one_side_monotone_local(tpl, v_obs, side_sign)
xq = nan(size(v_obs(:)));
if isempty(v_obs)
    return;
end
xg = tpl.x_grid(:);
vg = tpl.v_grid(:);
valid = isfinite(xg) & isfinite(vg);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && all(isfinite(tpl.x_domain(1:2)))
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


function x_static = invert_template_nearest_crossing_local(tpl, v_obs, x_ref)
xg = tpl.x_grid(:);
vg = tpl.v_grid(:);
valid_grid = isfinite(xg) & isfinite(vg);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && all(isfinite(tpl.x_domain(1:2)))
    valid_grid = valid_grid & xg >= tpl.x_domain(1) & xg <= tpl.x_domain(2);
end
xg = xg(valid_grid);
vg = vg(valid_grid);
x_static = nan(size(v_obs(:)));
for i = 1:numel(v_obs)
    if ~isfinite(v_obs(i)) || ~isfinite(x_ref(i)) || isempty(xg)
        continue;
    end
    diff_v = vg - v_obs(i);
    crossing_x = [];
    exact_idx = find(abs(diff_v) <= 1e-10);
    if ~isempty(exact_idx)
        crossing_x = xg(exact_idx);
    end
    for k = 1:(numel(diff_v) - 1)
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k + 1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k + 1) > 0
            continue;
        end
        denom = vg(k + 1) - vg(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (v_obs(i) - vg(k)) / denom;
        crossing_x(end + 1, 1) = xg(k) + alpha * (xg(k + 1) - xg(k)); %#ok<AGROW>
    end
    if isempty(crossing_x)
        [~, idx] = min(abs(diff_v));
        x_static(i) = xg(idx);
    else
        [~, idx] = min(abs(crossing_x - x_ref(i)));
        x_static(i) = crossing_x(idx);
    end
end
end


function T = diagnose_current_template_branches_local(Template, sensors, blade_id)
rows = repmat(struct( ...
    'sensor_id', NaN, 'blade_id', blade_id, 'xc_mm', NaN, ...
    'x_domain_left_mm', NaN, 'x_domain_right_mm', NaN, ...
    'left_point_count', NaN, 'right_point_count', NaN, ...
    'left_violation_fraction', NaN, 'right_violation_fraction', NaN, ...
    'left_voltage_span_v', NaN, 'right_voltage_span_v', NaN), numel(sensors), 1);
for i = 1:numel(sensors)
    sid = sensors(i);
    tpl = get_current_template_entry_local(Template, sid, blade_id);
    rows(i).sensor_id = sid;
    if isempty(tpl)
        continue;
    end
    rows(i).xc_mm = get_template_xc_local(tpl);
    rows(i).x_domain_left_mm = tpl.x_domain(1);
    rows(i).x_domain_right_mm = tpl.x_domain(2);
    [rows(i).left_point_count, rows(i).left_violation_fraction, rows(i).left_voltage_span_v] = ...
        branch_monotonicity_local(tpl, -1);
    [rows(i).right_point_count, rows(i).right_violation_fraction, rows(i).right_voltage_span_v] = ...
        branch_monotonicity_local(tpl, 1);
end
T = struct2table(rows);
end


function [n, frac, span_v] = branch_monotonicity_local(tpl, side_sign)
x = tpl.x_grid(:);
v = tpl.v_grid(:);
valid = isfinite(x) & isfinite(v);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && all(isfinite(tpl.x_domain(1:2)))
    valid = valid & x >= tpl.x_domain(1) & x <= tpl.x_domain(2);
end
if side_sign < 0
    valid = valid & x <= 0;
    [~, order] = sort(x(valid), 'ascend');
else
    valid = valid & x >= 0;
    [~, order] = sort(x(valid), 'descend');
end
vv = v(valid);
vv = vv(order);
n = numel(vv);
if n < 3
    frac = NaN;
    span_v = NaN;
    return;
end
dv = diff(vv);
tol = max(1e-4, 0.002 * range(vv));
frac = nnz(dv < -tol) / numel(dv);
span_v = max(vv, [], 'omitnan') - min(vv, [], 'omitnan');
end


function row = summarize_window_local(iw, lap_start, lap_end, bundle, T)
row = make_empty_summary_row_local();
row.window_id = iw;
row.lap_start = lap_start;
row.lap_end = lap_end;
row.point_count = bundle.point_count;
row.rot_freq_mean_hz = bundle.rot_freq_mean_hz;
models = ["sg_success_control", ...
    "current_step04_branch_monotone_inverse", ...
    "current_step04_nearest_crossing_inverse"];
for model = models
    M = T(strcmpi(T.model, model), :);
    if isempty(M)
        continue;
    end
    top_list = strjoin(compose('%d', M.EO(1:min(5, height(M))).'), "|");
    idx14 = find(M.EO == 14, 1, 'first');
    eo14_rank = NaN;
    eo14_rmse = NaN;
    if ~isempty(idx14)
        eo14_rank = M.rank(idx14);
        eo14_rmse = M.weighted_voltage_rmse_v(idx14);
    end
    switch model
        case "sg_success_control"
            row.sg_eo = M.EO(1);
            row.sg_rmse_v = M.weighted_voltage_rmse_v(1);
            row.sg_eo14_rank = eo14_rank;
            row.sg_eo14_rmse_v = eo14_rmse;
            row.sg_top_list = string(top_list);
        case "current_step04_branch_monotone_inverse"
            row.current_mono_eo = M.EO(1);
            row.current_mono_rmse_v = M.weighted_voltage_rmse_v(1);
            row.current_mono_y_rmse_mm = M.weighted_y_rmse_mm(1);
            row.current_mono_eo14_rank = eo14_rank;
            row.current_mono_eo14_rmse_v = eo14_rmse;
            row.current_mono_top_list = string(top_list);
        case "current_step04_nearest_crossing_inverse"
            row.current_cross_eo = M.EO(1);
            row.current_cross_rmse_v = M.weighted_voltage_rmse_v(1);
            row.current_cross_y_rmse_mm = M.weighted_y_rmse_mm(1);
            row.current_cross_eo14_rank = eo14_rank;
            row.current_cross_eo14_rmse_v = eo14_rmse;
            row.current_cross_top_list = string(top_list);
    end
end
end


function row = make_empty_summary_row_local()
row = struct( ...
    'window_id', NaN, 'lap_start', NaN, 'lap_end', NaN, ...
    'point_count', NaN, 'rot_freq_mean_hz', NaN, ...
    'sg_eo', NaN, 'sg_rmse_v', NaN, 'sg_eo14_rank', NaN, ...
    'sg_eo14_rmse_v', NaN, 'sg_top_list', "", ...
    'current_mono_eo', NaN, 'current_mono_rmse_v', NaN, ...
    'current_mono_y_rmse_mm', NaN, 'current_mono_eo14_rank', NaN, ...
    'current_mono_eo14_rmse_v', NaN, 'current_mono_top_list', "", ...
    'current_cross_eo', NaN, 'current_cross_rmse_v', NaN, ...
    'current_cross_y_rmse_mm', NaN, 'current_cross_eo14_rank', NaN, ...
    'current_cross_eo14_rmse_v', NaN, 'current_cross_top_list', "");
end


function eo_candidates = build_eo_candidates_local(rot_freq_hz, freq_search_hz, eo_pad)
rot_freq_hz = max(rot_freq_hz, eps);
eo_min = max(1, ceil(min(freq_search_hz) / rot_freq_hz) - max(0, eo_pad));
eo_max = max(eo_min, floor(max(freq_search_hz) / rot_freq_hz) + max(0, eo_pad));
eo_candidates = eo_min:eo_max;
eo_candidates = eo_candidates(eo_candidates > 0);
end


function Template = load_current_step04_template_local(cfg)
override_file = strtrim(getenv('STEP05_DIAG_TEMPLATE_FILE'));
if ~isempty(override_file)
    f = override_file;
else
    f = fullfile(cfg.step04_output_dir, ...
        'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');
end
if ~isfile(f)
    error('Missing current Step04 template: %s', f);
end
loaded = load(f, 'Template');
Template = loaded.Template;
end


function tpl = get_current_template_entry_local(Template, sid, blade_id)
tpl = [];
if ~isfield(Template, 'SensorBlade')
    return;
end
for i = 1:numel(Template.SensorBlade)
    entry = Template.SensorBlade(i);
    if double(entry.sensor_id) == double(sid) && double(entry.blade_id) == double(blade_id)
        tpl = entry;
        return;
    end
end
end


function xc = get_template_xc_local(tpl)
xc = 0;
if isfield(tpl, 'xc') && isfinite(tpl.xc)
    xc = tpl.xc;
elseif isfield(tpl, 'xc_mm') && isfinite(tpl.xc_mm)
    xc = tpl.xc_mm;
end
end


function sg_root = resolve_sg_success_root_local(route_dir)
sg_root = strtrim(getenv('SG_SUCCESS_ROOT'));
if ~isempty(sg_root) && isfolder(sg_root)
    return;
end
validation_root = fileparts(route_dir);
repo_root = fileparts(validation_root);
common_root = fileparts(repo_root);
candidate = fullfile(common_root, ...
    '7超高斯模型-权重-瞬态', '程序', '直叶片验证', '实验验证', ...
    '超高斯', '20251222适配-3号传感间隙变化-改');
if isfolder(candidate)
    sg_root = candidate;
else
    sg_root = '';
end
end


function cleanup_obj = add_sg_paths_for_whole_diagnostic_local(sg_root)
path_text = genpath(sg_root);
addpath(path_text, '-begin');
cleanup_obj = onCleanup(@() rmpath(path_text));
end


function value = read_env_number_local(name, default_value)
txt = strtrim(getenv(name));
if isempty(txt)
    value = default_value;
    return;
end
if any(strcmpi(txt, {'inf', 'all'}))
    value = inf;
    return;
end
value = str2double(txt);
if ~isfinite(value) && ~isinf(value)
    value = default_value;
end
end


function values = read_env_vector_local(name, default_values)
txt = strtrim(getenv(name));
if isempty(txt)
    values = default_values(:).';
    return;
end
parts = regexp(txt, '[,;\s]+', 'split');
values = cellfun(@str2double, parts);
values = values(isfinite(values));
if isempty(values)
    values = default_values(:).';
else
    values = unique(values(:).', 'stable');
end
end


function tag = join_sensor_tag_local(sensors)
tag = sprintf('%d', sensors(:).');
end


function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end


function print_sequence_local(W)
if isempty(W)
    return;
end
fprintf('\nEO sequences:\n');
fprintf('  SG control:          %s\n', strjoin(compose('%d', W.sg_eo.'), "|"));
fprintf('  Current branch inv:  %s\n', strjoin(compose('%d', W.current_mono_eo.'), "|"));
fprintf('  Current crossing inv:%s\n', strjoin(compose('%d', W.current_cross_eo.'), "|"));
fprintf('\nEO14 median ranks: SG %.2f, branch %.2f, crossing %.2f\n', ...
    median(W.sg_eo14_rank, 'omitnan'), ...
    median(W.current_mono_eo14_rank, 'omitnan'), ...
    median(W.current_cross_eo14_rank, 'omitnan'));
end
