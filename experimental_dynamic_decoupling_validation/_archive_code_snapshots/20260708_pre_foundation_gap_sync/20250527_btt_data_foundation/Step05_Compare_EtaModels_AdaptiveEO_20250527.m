clc; close all;

%STEP05_COMPARE_ETAMODELS_ADAPTIVEEO_20250527
% Diagnostic comparison wrapper for the 20250527 foundation Step05.
%
% The formal Step05 entry point runs one active eta model. This script runs
% the same main program for NoEta and bounded window-eta, then writes compact
% comparison CSV files. Keep this diagnostic route separate from the main
% identification entry point.

cfg = BTTProjectConfig_20250527();

case_name = '20250526_2500-3500_t400';
blade_id = 1;
sensor_tag = ['S', sprintf('%d', cfg.sensor_ids)];

models = {'none', 'window_bounded'};
tags = {'eta_none', 'eta_window_bounded'};

env_names = { ...
    'STEP05_ETA_MODEL', ...
    'STEP05_OUTPUT_TAG', ...
    'STEP05_SHOW_PLOTS', ...
    'STEP05_SAVE_FIGURES', ...
    'STEP05_FULL_EO_SCAN_DIAGNOSTIC_ENABLE'};
old_env = capture_env_local(env_names);
cleanup_obj = onCleanup(@() restore_env_local(env_names, old_env)); %#ok<NASGU>

if isempty(strtrim(getenv('STEP05_SHOW_PLOTS')))
    setenv('STEP05_SHOW_PLOTS', '0');
end

if isempty(strtrim(getenv('STEP05_SAVE_FIGURES')))
    setenv('STEP05_SAVE_FIGURES', '0');
end

if isempty(strtrim(getenv('STEP05_FULL_EO_SCAN_DIAGNOSTIC_ENABLE')))
    setenv('STEP05_FULL_EO_SCAN_DIAGNOSTIC_ENABLE', '1');
end

fprintf('\n=== Step05 eta-model comparison: 20250527 ===\n');
fprintf('Case: %s | Blade %d | Sensors %s\n', case_name, blade_id, sensor_tag);

result_files = strings(numel(models), 1);
trend_files = strings(numel(models), 1);

for im = 1:numel(models)
    fprintf('\n--- Running model: %s ---\n', models{im});

    setenv('STEP05_ETA_MODEL', models{im});
    setenv('STEP05_OUTPUT_TAG', tags{im});

    Step05_SingleSync_DirectTemplate_Identification_20250527;

    out_dir = fullfile(cfg.output_root, ...
        ['step05_single_sync_direct_template_', tags{im}], case_name);

    result_files(im) = fullfile(out_dir, sprintf( ...
        'Result_Step05_SingleSyncDirectTemplate_B%d_%s_%s.mat', ...
        blade_id, sensor_tag, cfg.dataset));

    trend_files(im) = fullfile(out_dir, sprintf( ...
        'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
        blade_id, sensor_tag, cfg.dataset));
end

compare_dir = fullfile(cfg.output_root, 'step05_eta_model_compare', case_name);
if ~exist(compare_dir, 'dir')
    mkdir(compare_dir);
end

[WindowCompare, Summary, PairedCompare] = build_eta_model_compare_tables_local( ...
    models, result_files, trend_files);

window_csv = fullfile(compare_dir, sprintf( ...
    'Step05_EtaModelCompare_Window_%s_B%d_%s.csv', ...
    cfg.dataset, blade_id, sensor_tag));

summary_csv = fullfile(compare_dir, sprintf( ...
    'Step05_EtaModelCompare_Summary_%s_B%d_%s.csv', ...
    cfg.dataset, blade_id, sensor_tag));

paired_csv = fullfile(compare_dir, sprintf( ...
    'Step05_EtaModelCompare_Paired_%s_B%d_%s.csv', ...
    cfg.dataset, blade_id, sensor_tag));

writetable(WindowCompare, window_csv);
writetable(Summary, summary_csv);
writetable(PairedCompare, paired_csv);

fprintf('\nSaved eta-model window comparison:\n  %s\n', window_csv);
fprintf('Saved eta-model summary comparison:\n  %s\n', summary_csv);
fprintf('Saved eta-model paired comparison:\n  %s\n', paired_csv);


function old_env = capture_env_local(names)
old_env = cell(size(names));
for i = 1:numel(names)
    old_env{i} = getenv(names{i});
end
end


function restore_env_local(names, old_env)
for i = 1:numel(names)
    setenv(names{i}, old_env{i});
end
end


function [WindowCompare, Summary, PairedCompare] = build_eta_model_compare_tables_local( ...
    models, result_files, trend_files)

WindowCompare = table();
Summary = table();

for im = 1:numel(models)
    model_name = string(models{im});

    if ~isfile(trend_files(im))
        warning('Missing trend file for model %s: %s', ...
            model_name, trend_files(im));
        continue;
    end

    T = readtable(trend_files(im), 'TextType', 'string');
    T.run_eta_model = repmat(model_name, height(T), 1);
    T.result_file = repmat(result_files(im), height(T), 1);
    T.trend_file = repmat(trend_files(im), height(T), 1);

    WindowCompare = [WindowCompare; T]; %#ok<AGROW>

    ok = strcmpi(T.status, "ok") & isfinite(T.EO_id);

    if any(ok)
        unique_eo = join(string(unique(T.EO_id(ok)).'), ",");
        median_eo = median(T.EO_id(ok), 'omitnan');
        median_fn_hz = median(T.fn_id(ok), 'omitnan');
        median_A_mm = median(T.A_id(ok), 'omitnan');
        median_dx_c_mm = median(T.dx_c_id(ok), 'omitnan');
        median_eta_max_abs_mm = median(T.sensor_eta_max_abs_mm(ok), 'omitnan');
        median_weighted_rmse_v = median(T.weighted_voltage_rmse(ok), 'omitnan');
    else
        unique_eo = "";
        median_eo = NaN;
        median_fn_hz = NaN;
        median_A_mm = NaN;
        median_dx_c_mm = NaN;
        median_eta_max_abs_mm = NaN;
        median_weighted_rmse_v = NaN;
    end

    Srow = table();
    Srow.run_eta_model = model_name;
    Srow.window_count = height(T);
    Srow.ok_window_count = nnz(ok);
    Srow.unique_eo = unique_eo;
    Srow.median_eo = median_eo;
    Srow.median_fn_hz = median_fn_hz;
    Srow.median_A_mm = median_A_mm;
    Srow.median_dx_c_mm = median_dx_c_mm;
    Srow.median_eta_max_abs_mm = median_eta_max_abs_mm;
    Srow.median_weighted_rmse_v = median_weighted_rmse_v;
    Srow.result_file = result_files(im);

    Summary = [Summary; Srow]; %#ok<AGROW>
end

PairedCompare = build_paired_eta_compare_local(models, result_files, trend_files);
end


function PairedCompare = build_paired_eta_compare_local(models, result_files, trend_files)
PairedCompare = table();

model_names = string(models);
idx_none = find(strcmpi(model_names, "none"), 1, 'first');
idx_eta = find(strcmpi(model_names, "window_bounded"), 1, 'first');

if isempty(idx_none) || isempty(idx_eta) || ...
        ~isfile(trend_files(idx_none)) || ~isfile(trend_files(idx_eta))
    return;
end

Tn = readtable(trend_files(idx_none), 'TextType', 'string');
Tw = readtable(trend_files(idx_eta), 'TextType', 'string');

if ~ismember('window_id', Tn.Properties.VariableNames) || ...
        ~ismember('window_id', Tw.Properties.VariableNames)
    return;
end

common_windows = intersect(Tn.window_id, Tw.window_id);
if isempty(common_windows)
    return;
end

settings = read_compare_settings_local(result_files(idx_eta), result_files(idx_none));
diag_none = collect_full_scan_summary_local(result_files(idx_none), "none");
diag_eta = collect_full_scan_summary_local(result_files(idx_eta), "window_bounded");

rows = repmat(make_empty_paired_row_local(), numel(common_windows), 1);

for i = 1:numel(common_windows)
    wid = common_windows(i);
    rn = Tn(Tn.window_id == wid, :);
    rw = Tw(Tw.window_id == wid, :);

    if height(rn) > 1
        rn = rn(1, :);
    end
    if height(rw) > 1
        rw = rw(1, :);
    end

    rows(i).window_id = wid;
    rows(i).status_none = get_table_string_local(rn, 'status', "");
    rows(i).status_window_bounded = get_table_string_local(rw, 'status', "");
    rows(i).EO_none = get_table_numeric_local(rn, 'EO_id', NaN);
    rows(i).EO_window_bounded = get_table_numeric_local(rw, 'EO_id', NaN);
    rows(i).same_EO = isequaln(rows(i).EO_none, rows(i).EO_window_bounded) && ...
        isfinite(rows(i).EO_none);
    rows(i).RMSE_none_v = get_table_numeric_local(rn, 'weighted_voltage_rmse', NaN);
    rows(i).RMSE_window_bounded_v = get_table_numeric_local(rw, 'weighted_voltage_rmse', NaN);
    rows(i).RMSE_improvement_v = rows(i).RMSE_none_v - rows(i).RMSE_window_bounded_v;
    rows(i).RMSE_improvement_ratio = rows(i).RMSE_improvement_v / ...
        max(abs(rows(i).RMSE_none_v), eps);
    rows(i).eta_max_abs_mm = get_table_numeric_local(rw, 'sensor_eta_max_abs_mm', NaN);
    rows(i).eta_at_bound = isfinite(rows(i).eta_max_abs_mm) && ...
        rows(i).eta_max_abs_mm >= settings.eta_boundary_fraction * ...
        max(abs(settings.sensor_eta_limit_mm), eps);

    dn = lookup_diag_row_local(diag_none, wid);
    de = lookup_diag_row_local(diag_eta, wid);

    rows(i).full_scan_applied_none = get_table_logical_local(dn, 'full_scan_applied', false);
    rows(i).full_scan_best_EO_none = get_table_numeric_local(dn, 'full_scan_best_EO', NaN);
    rows(i).full_scan_gap_ratio_none = get_table_numeric_local(dn, 'full_scan_gap_ratio', NaN);
    rows(i).adaptive_missed_full_scan_best_none = ...
        isfinite(rows(i).full_scan_best_EO_none) && ...
        isfinite(rows(i).EO_none) && ...
        rows(i).full_scan_best_EO_none ~= rows(i).EO_none;

    rows(i).full_scan_applied_window_bounded = ...
        get_table_logical_local(de, 'full_scan_applied', false);
    rows(i).full_scan_best_EO_window_bounded = ...
        get_table_numeric_local(de, 'full_scan_best_EO', NaN);
    rows(i).full_scan_gap_ratio_window_bounded = ...
        get_table_numeric_local(de, 'full_scan_gap_ratio', NaN);
    rows(i).adaptive_missed_full_scan_best_window_bounded = ...
        isfinite(rows(i).full_scan_best_EO_window_bounded) && ...
        isfinite(rows(i).EO_window_bounded) && ...
        rows(i).full_scan_best_EO_window_bounded ~= rows(i).EO_window_bounded;

    rows(i).preferred_model_by_window = decide_paired_preference_local( ...
        rows(i), settings);
end

PairedCompare = struct2table(rows);
end


function row = make_empty_paired_row_local()
row = struct( ...
    'window_id', NaN, ...
    'status_none', "", ...
    'status_window_bounded', "", ...
    'EO_none', NaN, ...
    'EO_window_bounded', NaN, ...
    'same_EO', false, ...
    'RMSE_none_v', NaN, ...
    'RMSE_window_bounded_v', NaN, ...
    'RMSE_improvement_v', NaN, ...
    'RMSE_improvement_ratio', NaN, ...
    'eta_max_abs_mm', NaN, ...
    'eta_at_bound', false, ...
    'full_scan_applied_none', false, ...
    'full_scan_best_EO_none', NaN, ...
    'full_scan_gap_ratio_none', NaN, ...
    'adaptive_missed_full_scan_best_none', false, ...
    'full_scan_applied_window_bounded', false, ...
    'full_scan_best_EO_window_bounded', NaN, ...
    'full_scan_gap_ratio_window_bounded', NaN, ...
    'adaptive_missed_full_scan_best_window_bounded', false, ...
    'preferred_model_by_window', "");
end


function settings = read_compare_settings_local(primary_result_file, fallback_result_file)
settings = struct();
settings.sensor_eta_limit_mm = 0.20;
settings.eta_boundary_fraction = 0.90;
settings.eta_accept_rmse_improvement_ratio = 0.03;
settings.eta_accept_rmse_improvement_v = 0.003;

files = [primary_result_file, fallback_result_file];

for i = 1:numel(files)
    if ~isfile(files(i))
        continue;
    end

    loaded = load(files(i), 'Result');
    if ~isfield(loaded, 'Result') || ~isfield(loaded.Result, 'AnalysisSettings')
        continue;
    end

    S = loaded.Result.AnalysisSettings;
    settings.sensor_eta_limit_mm = get_struct_numeric_local( ...
        S, 'sensor_eta_limit_mm', settings.sensor_eta_limit_mm);
    settings.eta_boundary_fraction = get_struct_numeric_local( ...
        S, 'eta_boundary_fraction', settings.eta_boundary_fraction);
    settings.eta_accept_rmse_improvement_ratio = get_struct_numeric_local( ...
        S, 'eta_accept_rmse_improvement_ratio', ...
        settings.eta_accept_rmse_improvement_ratio);
    settings.eta_accept_rmse_improvement_v = get_struct_numeric_local( ...
        S, 'eta_accept_rmse_improvement_v', settings.eta_accept_rmse_improvement_v);
    return;
end
end


function diag_table = collect_full_scan_summary_local(result_file, eta_model)
diag_table = table();

if ~isfile(result_file)
    return;
end

loaded = load(result_file, 'Result');

if ~isfield(loaded, 'Result') || ~isfield(loaded.Result, 'WindowResult')
    return;
end

WR = loaded.Result.WindowResult;
rows = repmat(make_empty_full_scan_row_local(), numel(WR), 1);

for i = 1:numel(WR)
    rows(i).window_id = get_struct_numeric_local(WR(i), 'window_id', i);

    D = [];
    if isfield(WR(i), 'FullScanDiagnosticFinal')
        D = WR(i).FullScanDiagnosticFinal;
    elseif isfield(WR(i), 'Result') && isfield(WR(i).Result, 'FullScanDiagnosticFinal')
        D = WR(i).Result.FullScanDiagnosticFinal;
    end

    if isempty(D) || ~isstruct(D)
        continue;
    end

    rows(i).full_scan_applied = get_struct_logical_local(D, 'applied', false);
    rows(i).trigger_reason = string(get_struct_string_local(D, 'trigger_reason', ""));

    if ~rows(i).full_scan_applied
        continue;
    end

    switch lower(char(eta_model))
        case {'none', 'noeta', 'eta_none'}
            model_result = get_struct_field_local(D, 'NoEta', []);
        case {'window_bounded', 'bounded', 'witheta'}
            model_result = get_struct_field_local(D, 'WindowEta', []);
        otherwise
            model_result = [];
    end

    [best_eo, best_rmse, gap_ratio] = best_eo_from_full_scan_result_local(model_result);
    rows(i).full_scan_best_EO = best_eo;
    rows(i).full_scan_best_rmse_v = best_rmse;
    rows(i).full_scan_gap_ratio = gap_ratio;
end

diag_table = struct2table(rows);
end


function row = make_empty_full_scan_row_local()
row = struct( ...
    'window_id', NaN, ...
    'full_scan_applied', false, ...
    'trigger_reason', "", ...
    'full_scan_best_EO', NaN, ...
    'full_scan_best_rmse_v', NaN, ...
    'full_scan_gap_ratio', NaN);
end


function [best_eo, best_rmse, gap_ratio] = best_eo_from_full_scan_result_local(model_result)
best_eo = NaN;
best_rmse = NaN;
gap_ratio = NaN;

if isempty(model_result) || ~isstruct(model_result) || ...
        ~isfield(model_result, 'CandidateTable') || ...
        isempty(model_result.CandidateTable) || ...
        ~istable(model_result.CandidateTable)
    return;
end

T = model_result.CandidateTable;

if ~all(ismember({'EO', 'weighted_voltage_rmse'}, T.Properties.VariableNames))
    return;
end

ok = isfinite(T.EO) & isfinite(T.weighted_voltage_rmse);

if ismember('status', T.Properties.VariableNames)
    ok = ok & strcmpi(string(T.status), "ok");
end

T = T(ok, :);

if isempty(T)
    return;
end

T = sortrows(T, {'weighted_voltage_rmse', 'EO'}, {'ascend', 'ascend'});
best_eo = T.EO(1);
best_rmse = T.weighted_voltage_rmse(1);

if height(T) >= 2
    gap_ratio = (T.weighted_voltage_rmse(2) - T.weighted_voltage_rmse(1)) / ...
        max(abs(T.weighted_voltage_rmse(1)), eps);
end
end


function decision = decide_paired_preference_local(row, settings)
none_ok = strcmpi(row.status_none, "ok") && isfinite(row.RMSE_none_v);
eta_ok = strcmpi(row.status_window_bounded, "ok") && isfinite(row.RMSE_window_bounded_v);

if ~none_ok && eta_ok
    decision = "window_bounded_only_noeta_failed";
elseif none_ok && ~eta_ok
    decision = "none_only_eta_failed";
elseif ~none_ok && ~eta_ok
    decision = "neither_ok";
elseif row.adaptive_missed_full_scan_best_none || ...
        row.adaptive_missed_full_scan_best_window_bounded
    decision = "diagnostic_full_scan_disagrees";
elseif ~row.same_EO
    decision = "none_primary_eta_changes_eo";
elseif row.eta_at_bound
    decision = "none_primary_eta_at_bound";
elseif row.RMSE_improvement_v >= settings.eta_accept_rmse_improvement_v && ...
        row.RMSE_improvement_ratio >= settings.eta_accept_rmse_improvement_ratio
    decision = "window_bounded_supported_same_eo";
else
    decision = "none_primary_eta_improvement_small";
end
end


function sub = lookup_diag_row_local(T, window_id)
sub = table();
if isempty(T) || ~istable(T) || ~ismember('window_id', T.Properties.VariableNames)
    return;
end

sub = T(T.window_id == window_id, :);
if height(sub) > 1
    sub = sub(1, :);
end
end


function value = get_table_numeric_local(T, field_name, default_value)
value = default_value;
if isempty(T) || ~istable(T) || ~ismember(field_name, T.Properties.VariableNames)
    return;
end

candidate = T.(field_name);
if ~isempty(candidate) && isnumeric(candidate) && isfinite(candidate(1))
    value = candidate(1);
end
end


function value = get_table_string_local(T, field_name, default_value)
value = string(default_value);
if isempty(T) || ~istable(T) || ~ismember(field_name, T.Properties.VariableNames)
    return;
end

candidate = T.(field_name);
if ~isempty(candidate)
    value = string(candidate(1));
end
end


function value = get_table_logical_local(T, field_name, default_value)
value = default_value;
if isempty(T) || ~istable(T) || ~ismember(field_name, T.Properties.VariableNames)
    return;
end

candidate = T.(field_name);
if ~isempty(candidate)
    value = logical(candidate(1));
end
end


function value = get_struct_numeric_local(S, field_name, default_value)
value = default_value;
if isstruct(S) && isfield(S, field_name) && ...
        isnumeric(S.(field_name)) && isscalar(S.(field_name)) && ...
        isfinite(S.(field_name))
    value = S.(field_name);
end
end


function value = get_struct_logical_local(S, field_name, default_value)
value = default_value;
if isstruct(S) && isfield(S, field_name) && isscalar(S.(field_name))
    value = logical(S.(field_name));
end
end


function value = get_struct_string_local(S, field_name, default_value)
value = string(default_value);
if isstruct(S) && isfield(S, field_name)
    value = string(S.(field_name));
end
end


function value = get_struct_field_local(S, field_name, default_value)
value = default_value;
if isstruct(S) && isfield(S, field_name)
    value = S.(field_name);
end
end
