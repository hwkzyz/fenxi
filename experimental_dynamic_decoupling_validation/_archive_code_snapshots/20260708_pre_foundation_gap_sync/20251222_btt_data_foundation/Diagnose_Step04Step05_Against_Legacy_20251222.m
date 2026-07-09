%% Diagnose_Step04Step05_Against_Legacy_20251222
% Compare the new foundation low-speed template and Step05 candidates
% against the established 20251222 legacy main-route result. This script is
% diagnostic only; it does not feed legacy data into the foundation pipeline.

clear; clc;

this_file = mfilename('fullpath');
foundation_dir = fileparts(this_file);
repo_root = fileparts(foundation_dir);
addpath(foundation_dir);

cfg = BTTDataConfig_20251222();

diag_dir = fullfile(cfg.output_root, 'diagnostics');
if exist(diag_dir, 'dir') ~= 7
    mkdir(diag_dir);
end

blade_id = 1;
sensor_ids = [1 2 3];
sensor_tag = sprintf('S%s', strjoin(string(sensor_ids), ''));

legacy_result_file = fullfile(repo_root, '20251222_low_speed_rotating_calibration', ...
    'output', 'identification', ...
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_20251222.mat', ...
    blade_id, sensor_tag));

foundation_result_file = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_corr_joint_unconstrained_w3', ...
    cfg.dynamic_cases{1}, ...
    sprintf('Result_Step05_SingleSyncDirectTemplate_B%d_%s_%s.mat', ...
    blade_id, sensor_tag, cfg.dataset));

foundation_template_file = fullfile(cfg.output_root, 'step04_low_speed_template', ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', sensor_tag, cfg.dataset));

if ~isfile(legacy_result_file)
    error('Legacy result not found: %s', legacy_result_file);
end
if ~isfile(foundation_result_file)
    error('Foundation result not found: %s', foundation_result_file);
end
if ~isfile(foundation_template_file)
    error('Foundation template not found: %s', foundation_template_file);
end

legacy_loaded = load(legacy_result_file, 'Result');
legacy = legacy_loaded.Result;

foundation_loaded = load(foundation_result_file, 'Result');
foundation = foundation_loaded.Result;

template_loaded = load(foundation_template_file, 'Template');
Template = template_loaded.Template;

legacy_template_file = legacy.TemplateFile;
if ~isfile(legacy_template_file)
    error('Legacy template referenced by legacy result is not readable: %s', legacy_template_file);
end
legacy_tpl_loaded = load(legacy_template_file);
legacy_tpl = pick_legacy_template_struct_local(legacy_tpl_loaded);

fprintf('\n=== Diagnose Step04/Step05 against legacy route ===\n');
fprintf('Legacy result:     %s\n', legacy_result_file);
fprintf('Foundation result: %s\n', foundation_result_file);
fprintf('Legacy template:   %s\n', legacy_template_file);
fprintf('Foundation tpl:    %s\n', foundation_template_file);

window_compare = build_window_comparison_local(legacy, foundation);
candidate_compare = build_candidate_comparison_local(legacy, foundation);
template_compare = build_template_comparison_local(legacy_tpl, Template, sensor_ids, blade_id);

window_file = fullfile(diag_dir, sprintf('Step05_LegacyVsFoundation_WindowCompare_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset));
candidate_file = fullfile(diag_dir, sprintf('Step05_LegacyVsFoundation_CandidateCompare_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset));
template_file = fullfile(diag_dir, sprintf('Step04_LegacyVsFoundation_TemplateCompare_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset));

writetable(window_compare, window_file);
writetable(candidate_compare, candidate_file);
writetable(template_compare, template_file);

fprintf('\nSaved diagnostics:\n  %s\n  %s\n  %s\n', window_file, candidate_file, template_file);

disp(' ');
disp('Window comparison preview:');
disp(window_compare);
disp(' ');
disp('Candidate comparison preview:');
disp(candidate_compare(1:min(18, height(candidate_compare)), :));
disp(' ');
disp('Template comparison:');
disp(template_compare);

function T = build_window_comparison_local(legacy, foundation)
legacy_trend = legacy.Trend;
foundation_trend = foundation.Trend;
n = min(3, min(height(legacy_trend), height(foundation_trend)));
rows = repmat(struct(), n, 1);
for i = 1:n
    rows(i).window_id = i;
    rows(i).legacy_lap_start = legacy_trend.lap_start(i);
    rows(i).legacy_lap_end = legacy_trend.lap_end(i);
    rows(i).foundation_lap_start = foundation_trend.lap_start(i);
    rows(i).foundation_lap_end = foundation_trend.lap_end(i);
    rows(i).legacy_center_time_s = get_table_value_local(legacy_trend, i, ...
        {'window_center_time', 'window_center_time_s'});
    rows(i).foundation_center_time_s = get_table_value_local(foundation_trend, i, ...
        {'window_center_time_s', 'window_center_time'});
    rows(i).center_time_delta_s = rows(i).foundation_center_time_s - rows(i).legacy_center_time_s;
    rows(i).legacy_rot_freq_hz = legacy_trend.rot_freq_mean_hz(i);
    rows(i).foundation_rot_freq_hz = foundation_trend.rot_freq_mean_hz(i);
    rows(i).legacy_EO = legacy_trend.EO_id(i);
    rows(i).foundation_EO = foundation_trend.EO_id(i);
    rows(i).legacy_fn_hz = legacy_trend.fn_id(i);
    rows(i).foundation_fn_hz = foundation_trend.fn_id(i);
    rows(i).legacy_points = legacy_trend.point_count(i);
    rows(i).foundation_points = foundation_trend.point_count(i);
    rows(i).foundation_candidate_points = foundation_trend.candidate_point_count(i);
    rows(i).legacy_weighted_rmse = legacy_trend.weighted_voltage_rmse(i);
    rows(i).foundation_weighted_metric = foundation_trend.weighted_voltage_rmse(i);
    rows(i).foundation_status = string(foundation_trend.status(i));
    rows(i).foundation_top_list = string(foundation_trend.eo_top_list(i));
end
T = struct2table(rows);
end

function T = build_candidate_comparison_local(legacy, foundation)
rows = repmat(struct('window_id', NaN, 'source', "", 'rank', NaN, 'EO', NaN, ...
    'fn_hz', NaN, 'A_mm', NaN, 'phi_rad', NaN, 'dx_c_mm', NaN, ...
    'eta_max_abs_mm', NaN, 'objective', NaN, 'weighted_metric', NaN, ...
    'plain_rmse_v', NaN, 'point_count', NaN), 0, 1);

n = min(3, min(numel(legacy.WindowResult), numel(foundation.WindowResult)));
for iw = 1:n
    old_ct = legacy.WindowResult(iw).CandidateTable;
    new_ct = foundation.WindowResult(iw).Result.CandidateTable;
    rows = append_candidate_rows_local(rows, iw, "legacy", old_ct, 8);
    rows = append_candidate_rows_local(rows, iw, "foundation_corr_unconstrained", new_ct, 8);
end
T = struct2table(rows);
end

function rows = append_candidate_rows_local(rows, iw, source, CT, top_n)
if isempty(CT) || height(CT) == 0
    return;
end
CT = sort_candidate_table_for_display_local(CT);
n = min(top_n, height(CT));
for ir = 1:n
    row = struct();
    row.window_id = iw;
    row.source = source;
    row.rank = ir;
    row.EO = get_table_value_local(CT, ir, {'EO', 'EO_id'});
    row.fn_hz = get_table_value_local(CT, ir, {'fn_hz', 'fn_id'});
    row.A_mm = get_table_value_local(CT, ir, {'A', 'A_id'});
    row.phi_rad = get_table_value_local(CT, ir, {'phi', 'phi_id_wrapped'});
    row.dx_c_mm = get_table_value_local(CT, ir, {'dx_c', 'dx_c_id'});
    row.eta_max_abs_mm = get_table_value_local(CT, ir, {'sensor_eta_max_abs', 'sensor_eta_max_abs_mm'});
    row.objective = get_table_value_local(CT, ir, {'objective', 'objective_score'});
    row.weighted_metric = get_table_value_local(CT, ir, {'weighted_voltage_rmse'});
    row.plain_rmse_v = get_table_value_local(CT, ir, {'plain_voltage_rmse'});
    row.point_count = get_table_value_local(CT, ir, {'point_count'});
    rows(end + 1, 1) = row; %#ok<AGROW>
end
end

function CT = sort_candidate_table_for_display_local(CT)
if ~istable(CT) || height(CT) == 0
    return;
end
vars = CT.Properties.VariableNames;
if any(strcmp(vars, 'objective_score'))
    CT = sortrows(CT, {'objective_score', 'weighted_voltage_rmse', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
elseif any(strcmp(vars, 'objective'))
    CT = sortrows(CT, {'objective', 'weighted_voltage_rmse', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
elseif any(strcmp(vars, 'weighted_voltage_rmse'))
    CT = sortrows(CT, {'weighted_voltage_rmse', 'EO'}, {'ascend', 'ascend'});
end
end

function T = build_template_comparison_local(legacy_tpl, Template, sensor_ids, blade_id)
rows = repmat(struct(), numel(sensor_ids), 1);
for is = 1:numel(sensor_ids)
    sid = sensor_ids(is);
    old_entry = get_legacy_template_entry_local(legacy_tpl, sid);
    new_entry = get_foundation_template_entry_local(Template, sid, blade_id);
    [curve_rmse, curve_corr, x_left, x_right, old_n, new_n] = compare_template_curves_local(old_entry, new_entry);
    rows(is).sensor_id = sid;
    rows(is).blade_id = blade_id;
    rows(is).legacy_xc_mm = get_struct_value_local(old_entry, {'xc_mm', 'x_center_mm', 'xc'}, NaN);
    rows(is).foundation_xc_mm = get_struct_value_local(new_entry, {'xc_mm', 'x_center_mm', 'xc'}, NaN);
    rows(is).xc_delta_mm = rows(is).foundation_xc_mm - rows(is).legacy_xc_mm;
    rows(is).legacy_domain_left_mm = get_domain_value_local(old_entry, 1);
    rows(is).legacy_domain_right_mm = get_domain_value_local(old_entry, 2);
    rows(is).foundation_domain_left_mm = get_domain_value_local(new_entry, 1);
    rows(is).foundation_domain_right_mm = get_domain_value_local(new_entry, 2);
    rows(is).overlap_left_mm = x_left;
    rows(is).overlap_right_mm = x_right;
    rows(is).legacy_grid_points = old_n;
    rows(is).foundation_grid_points = new_n;
    rows(is).overlap_curve_rmse_v = curve_rmse;
    rows(is).overlap_curve_corr = curve_corr;
end
T = struct2table(rows);
end

function tpl = pick_legacy_template_struct_local(S)
names = fieldnames(S);
preferred = {'Template', 'Tpl', 'TemplateLibrary', 'TemplateData'};
for i = 1:numel(preferred)
    if isfield(S, preferred{i})
        tpl = S.(preferred{i});
        return;
    end
end
tpl = S.(names{1});
end

function entry = get_foundation_template_entry_local(Template, sid, blade_id)
entry = struct();
if isfield(Template, 'SensorBlade') && ~isempty(Template.SensorBlade)
    SB = Template.SensorBlade;
    if isstruct(SB)
        for i = 1:numel(SB)
            if isfield(SB(i), 'sensor_id') && isfield(SB(i), 'blade_id') && ...
                    SB(i).sensor_id == sid && SB(i).blade_id == blade_id
                entry = SB(i);
                return;
            end
        end
        if ~isvector(SB) && size(SB, 1) >= sid && size(SB, 2) >= blade_id
            entry = SB(sid, blade_id);
            return;
        end
    end
end
end

function entry = get_legacy_template_entry_local(legacy_tpl, sid)
entry = struct();
if isstruct(legacy_tpl)
    if isfield(legacy_tpl, 'SensorTemplate')
        ST = legacy_tpl.SensorTemplate;
        if numel(ST) >= sid
            entry = ST(sid);
            return;
        end
    end
    if isfield(legacy_tpl, 'Sensor')
        ST = legacy_tpl.Sensor;
        for i = 1:numel(ST)
            if isfield(ST(i), 'sensor_id') && ST(i).sensor_id == sid
                entry = ST(i);
                return;
            end
        end
        if numel(ST) >= sid
            entry = ST(sid);
            return;
        end
    end
    if isfield(legacy_tpl, 'Template')
        T = legacy_tpl.Template;
        if isstruct(T) && numel(T) >= sid
            entry = T(sid);
            return;
        end
    end
    if isfield(legacy_tpl, 'Entries')
        E = legacy_tpl.Entries;
        if numel(E) >= sid
            entry = E(sid);
            return;
        end
    end
    if isfield(legacy_tpl, 'Templates')
        T = legacy_tpl.Templates;
        if numel(T) >= sid
            entry = T(sid);
            return;
        end
    end
    if isfield(legacy_tpl, 'x_grid') && isfield(legacy_tpl, 'v_grid')
        entry = legacy_tpl;
        return;
    end
end
end

function [curve_rmse, curve_corr, x_left, x_right, old_n, new_n] = compare_template_curves_local(old_entry, new_entry)
[x_old, v_old] = get_curve_local(old_entry);
[x_new, v_new] = get_curve_local(new_entry);
old_n = numel(x_old);
new_n = numel(x_new);
curve_rmse = NaN;
curve_corr = NaN;
x_left = NaN;
x_right = NaN;
if old_n < 5 || new_n < 5
    return;
end
x_left = max(min(x_old), min(x_new));
x_right = min(max(x_old), max(x_new));
if x_right <= x_left
    return;
end
xq = linspace(x_left, x_right, 600).';
vo = interp1(x_old, v_old, xq, 'linear', NaN);
vn = interp1(x_new, v_new, xq, 'linear', NaN);
valid = isfinite(vo) & isfinite(vn);
if nnz(valid) < 20
    return;
end
dv = vo(valid) - vn(valid);
curve_rmse = sqrt(mean(dv .^ 2, 'omitnan'));
vo = vo(valid) - mean(vo(valid), 'omitnan');
vn = vn(valid) - mean(vn(valid), 'omitnan');
den = sqrt(sum(vo .^ 2, 'omitnan') * sum(vn .^ 2, 'omitnan'));
if den > 0
    curve_corr = sum(vo .* vn, 'omitnan') / den;
end
end

function [x, v] = get_curve_local(entry)
x = [];
v = [];
if isempty(entry) || ~isstruct(entry)
    return;
end
x_names = {'x_grid', 'X_grid', 'x', 'X'};
v_names = {'v_grid', 'V_grid', 'v_template', 'V_template', 'v'};
for i = 1:numel(x_names)
    if isfield(entry, x_names{i})
        x = entry.(x_names{i});
        break;
    end
end
for i = 1:numel(v_names)
    if isfield(entry, v_names{i})
        v = entry.(v_names{i});
        break;
    end
end
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
[x, ia] = unique(x, 'stable');
v = v(ia);
end

function val = get_table_value_local(T, row_idx, names)
val = NaN;
for i = 1:numel(names)
    name = names{i};
    if istable(T) && any(strcmp(T.Properties.VariableNames, name))
        data = T.(name);
        if row_idx <= numel(data)
            val = data(row_idx);
            if iscell(val)
                val = val{1};
            end
            return;
        end
    end
end
end

function val = get_struct_value_local(S, names, default_val)
val = default_val;
if isempty(S) || ~isstruct(S)
    return;
end
for i = 1:numel(names)
    if isfield(S, names{i}) && ~isempty(S.(names{i}))
        tmp = S.(names{i});
        if isnumeric(tmp)
            val = tmp(1);
        end
        return;
    end
end
end

function val = get_domain_value_local(S, idx)
val = NaN;
if isempty(S) || ~isstruct(S)
    return;
end
names = {'x_domain', 'x_query_safe_domain', 'domain', 'x_range'};
for i = 1:numel(names)
    if isfield(S, names{i})
        d = S.(names{i});
        if isnumeric(d) && numel(d) >= idx
            val = d(idx);
            return;
        end
    end
end
end
