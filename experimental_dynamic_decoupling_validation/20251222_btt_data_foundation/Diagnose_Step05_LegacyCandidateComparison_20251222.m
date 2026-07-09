%% Diagnose_Step05_LegacyCandidateComparison_20251222
% Compare legacy 20251222 B1/S123 direct-template identification against
% foundation Step05 variants at the candidate and waveform-bundle levels.

clear; clc;

this_file = mfilename('fullpath');
foundation_dir = fileparts(this_file);
repo_root = fileparts(foundation_dir);
cfg = BTTProjectConfig_20251222();

diag_dir = fullfile(cfg.output_root, 'diagnostics');
fig_dir = fullfile(cfg.figure_root, 'diagnostics');
if exist(diag_dir, 'dir') ~= 7
    mkdir(diag_dir);
end
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

legacy_dir = fullfile(repo_root, '20251222_low_speed_rotating_calibration');
legacy_result_file = fullfile(legacy_dir, 'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_20251222.mat');
legacy_template_file = fullfile(legacy_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_OPRCenterStd_20251222.mat');
legacy_dynamic_map_file = fullfile(legacy_dir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat');

foundation_template_file = fullfile(cfg.output_root, 'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');

variant = [
    struct('tag', "legacy_noeta", 'kind', "legacy", 'file', legacy_result_file)
    struct('tag', "foundation_raw_current", 'kind', "foundation", 'file', fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_c22xcsmk', cfg.dynamic_cases{1}, ...
        'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat'))
    struct('tag', "foundation_dmap_current", 'kind', "foundation", 'file', fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_c22dmapst4smk', cfg.dynamic_cases{1}, ...
        'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat'))
    struct('tag', "foundation_dmap_legacytpl", 'kind', "foundation", 'file', fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_oldtpl_xcheck', cfg.dynamic_cases{1}, ...
        'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat'))
    struct('tag', "foundation_dmap_legacytpl_noeta", 'kind', "foundation", 'file', fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_oldtpl_noeta', cfg.dynamic_cases{1}, ...
        'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat'))
    struct('tag', "foundation_dmap_legacytpl_lmask", 'kind', "foundation", 'file', fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_oldtpl_lmask', cfg.dynamic_cases{1}, ...
        'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat'))
    ];

fprintf('\n=== Step05 legacy candidate comparison ===\n');
fprintf('Legacy result: %s\n', legacy_result_file);
fprintf('Legacy template: %s\n', legacy_template_file);
fprintf('Legacy DynamicMap: %s\n', legacy_dynamic_map_file);

TrendSummary = build_trend_summary_local(variant);
SeedRanks = build_seed_rank_table_local(variant, 1:3);
BundleSummary = build_bundle_summary_local(variant, 1:3);
TemplateCompare = build_template_compare_local( ...
    legacy_template_file, foundation_template_file, [1 2 3], 1);

trend_file = fullfile(diag_dir, ...
    'Step05_LegacyCandidateCompare_TrendSummary_B1_S123_20251222.csv');
seed_file = fullfile(diag_dir, ...
    'Step05_LegacyCandidateCompare_SeedRanks_B1_S123_20251222.csv');
bundle_file = fullfile(diag_dir, ...
    'Step05_LegacyCandidateCompare_BundleSummary_B1_S123_20251222.csv');
template_file = fullfile(diag_dir, ...
    'Step05_LegacyCandidateCompare_TemplateCompare_B1_S123_20251222.csv');

writetable(TrendSummary, trend_file);
writetable(SeedRanks, seed_file);
writetable(BundleSummary, bundle_file);
writetable(TemplateCompare, template_file);

plot_comparison_local(TrendSummary, SeedRanks, fig_dir);

fprintf('\nSaved diagnostics:\n  %s\n  %s\n  %s\n  %s\n', ...
    trend_file, seed_file, bundle_file, template_file);
disp(TrendSummary);
disp(TemplateCompare);


function T = build_trend_summary_local(variant)
rows = repmat(struct( ...
    'tag', "", 'kind', "", 'file_exists', false, 'window_count', 0, ...
    'ok_count', 0, 'eo14_count', 0, 'eo_sequence', "", ...
    'eo_mode', NaN, 'mean_freq_hz', NaN, 'mean_rmse_v', NaN, ...
    'mean_point_count', NaN, 'mean_dx_mm', NaN, 'mean_amp_mm', NaN, ...
    'status_sequence', ""), numel(variant), 1);

for i = 1:numel(variant)
    rows(i).tag = variant(i).tag;
    rows(i).kind = variant(i).kind;
    rows(i).file_exists = isfile(variant(i).file);
    if ~rows(i).file_exists
        continue;
    end
    R = load_result_local(variant(i).file);
    Tr = R.Trend;
    n = height(Tr);
    rows(i).window_count = n;
    status = get_table_string_local(Tr, 'status', repmat("ok", n, 1));
    eo = get_table_numeric_local(Tr, 'EO_id');
    ok = isfinite(eo) & strcmpi(status, "ok");
    rows(i).ok_count = nnz(ok);
    rows(i).eo14_count = nnz(ok & eo == 14);
    rows(i).eo_sequence = strjoin(string(eo(:).'), "|");
    rows(i).status_sequence = strjoin(status(:).', "|");
    if any(ok)
        rows(i).eo_mode = mode(eo(ok));
        rows(i).mean_freq_hz = mean(get_table_numeric_local(Tr, 'fn_id'), 'omitnan');
        rows(i).mean_rmse_v = mean(get_table_numeric_local(Tr, 'weighted_voltage_rmse'), 'omitnan');
        rows(i).mean_point_count = mean(get_table_numeric_local(Tr, 'point_count'), 'omitnan');
        rows(i).mean_dx_mm = mean(get_table_numeric_local(Tr, 'dx_c_id'), 'omitnan');
        rows(i).mean_amp_mm = mean(get_table_numeric_local(Tr, 'A_id'), 'omitnan');
    end
end
T = struct2table(rows);
end


function T = build_seed_rank_table_local(variant, windows)
rows = repmat(struct('tag', "", 'kind', "", 'window_id', NaN, ...
    'rank', NaN, 'EO', NaN, 'weighted_voltage_rmse', NaN, ...
    'plain_voltage_rmse', NaN, 'A_mm', NaN, 'dx_c_mm', NaN, ...
    'point_count', NaN), 0, 1);

for i = 1:numel(variant)
    if ~isfile(variant(i).file)
        continue;
    end
    R = load_result_local(variant(i).file);
    for iw = windows
        if iw > numel(R.WindowResult)
            continue;
        end
        S = get_seed_table_local(R.WindowResult(iw), variant(i).kind);
        if isempty(S)
            continue;
        end
        S = normalize_seed_table_local(S);
        n = min(12, height(S));
        for ir = 1:n
            row = struct();
            row.tag = variant(i).tag;
            row.kind = variant(i).kind;
            row.window_id = iw;
            row.rank = ir;
            row.EO = table_value_local(S, ir, {'EO', 'EO_id'});
            row.weighted_voltage_rmse = table_value_local(S, ir, {'weighted_voltage_rmse'});
            row.plain_voltage_rmse = table_value_local(S, ir, {'plain_voltage_rmse'});
            row.A_mm = table_value_local(S, ir, {'A', 'A_id'});
            row.dx_c_mm = table_value_local(S, ir, {'dx_c', 'dx_c_id'});
            row.point_count = table_value_local(S, ir, {'point_count'});
            rows(end + 1, 1) = row; %#ok<AGROW>
        end
    end
end
T = struct2table(rows);
end


function T = build_bundle_summary_local(variant, windows)
rows = repmat(struct('tag', "", 'kind', "", 'window_id', NaN, ...
    'point_count', NaN, 'valid_segment_count', NaN, ...
    'query_guard_mm', NaN, 'weight_mean', NaN, ...
    'ch1_points', NaN, 'ch2_points', NaN, 'ch3_points', NaN), 0, 1);
for i = 1:numel(variant)
    if ~isfile(variant(i).file)
        continue;
    end
    R = load_result_local(variant(i).file);
    for iw = windows
        if iw > numel(R.WindowResult)
            continue;
        end
        B = get_bundle_like_local(R.WindowResult(iw), variant(i).kind);
        if isempty(B)
            continue;
        end
        row = struct();
        row.tag = variant(i).tag;
        row.kind = variant(i).kind;
        row.window_id = iw;
        [sid, w, qg, pc, vseg] = unpack_bundle_fields_local(B, variant(i).kind);
        row.point_count = pc;
        row.valid_segment_count = vseg;
        row.query_guard_mm = qg;
        row.weight_mean = mean(w, 'omitnan');
        row.ch1_points = nnz(sid == 1);
        row.ch2_points = nnz(sid == 2);
        row.ch3_points = nnz(sid == 3);
        rows(end + 1, 1) = row; %#ok<AGROW>
    end
end
T = struct2table(rows);
end


function T = build_template_compare_local(legacy_template_file, foundation_template_file, sensors, blade_id)
L = load(legacy_template_file, 'Template');
F = load(foundation_template_file, 'Template');
rows = repmat(struct('sensor_id', NaN, 'blade_id', blade_id, ...
    'legacy_xc_mm', NaN, 'foundation_xc_mm', NaN, 'xc_delta_mm', NaN, ...
    'overlap_left_mm', NaN, 'overlap_right_mm', NaN, ...
    'curve_rmse_v', NaN, 'curve_corr', NaN, ...
    'legacy_domain_width_mm', NaN, 'foundation_domain_width_mm', NaN), numel(sensors), 1);
for i = 1:numel(sensors)
    sid = sensors(i);
    old = L.Template.Sensor([L.Template.Sensor.sensor_id] == sid);
    idx = find(arrayfun(@(e) e.sensor_id == sid && e.blade_id == blade_id, ...
        F.Template.SensorBlade), 1, 'first');
    new = F.Template.SensorBlade(idx);
    left = max(min(old.x_grid), min(new.x_grid));
    right = min(max(old.x_grid), max(new.x_grid));
    xq = linspace(left, right, 1000);
    vo = interp1(old.x_grid(:), old.v_grid(:), xq, 'pchip', NaN);
    vn = interp1(new.x_grid(:), new.v_grid(:), xq, 'pchip', NaN);
    rows(i).sensor_id = sid;
    rows(i).legacy_xc_mm = old.xc;
    rows(i).foundation_xc_mm = getfield_default_local(new, 'xc', NaN);
    rows(i).xc_delta_mm = rows(i).foundation_xc_mm - rows(i).legacy_xc_mm;
    rows(i).overlap_left_mm = left;
    rows(i).overlap_right_mm = right;
    rows(i).curve_rmse_v = sqrt(mean((vo - vn).^2, 'omitnan'));
    rows(i).curve_corr = corr(vo(:), vn(:), 'Rows', 'complete');
    rows(i).legacy_domain_width_mm = diff(old.x_domain);
    rows(i).foundation_domain_width_mm = diff(new.x_domain);
end
T = struct2table(rows);
end


function plot_comparison_local(TrendSummary, SeedRanks, fig_dir)
fig = figure('Color', 'w', 'Position', [80, 80, 1180, 720], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(categorical(TrendSummary.tag), TrendSummary.eo14_count ./ max(TrendSummary.window_count, 1));
ylabel('EO14 fraction');
title('Expected EO coverage');
ylim([0 1]);
grid on;
xtickangle(30);

nexttile;
bar(categorical(TrendSummary.tag), TrendSummary.mean_rmse_v);
ylabel('mean weighted RMSE (V)');
title('Residual');
grid on;
xtickangle(30);

nexttile;
bar(categorical(TrendSummary.tag), TrendSummary.mean_point_count);
ylabel('mean points');
title('Fitting-point count');
grid on;
xtickangle(30);

nexttile;
mask = SeedRanks.window_id == 1 & ismember(SeedRanks.EO, [8 13 14 20 21]);
S = SeedRanks(mask, :);
if ~isempty(S)
    scatter(categorical(S.tag), S.weighted_voltage_rmse, 34, S.EO, 'filled');
    cb = colorbar;
    cb.Label.String = 'EO';
end
ylabel('W1 seed weighted RMSE (V)');
title('Window 1 seed landscape');
grid on;
xtickangle(30);

out_png = fullfile(fig_dir, ...
    'Step05_LegacyCandidateCompare_B1_S123_20251222.png');
out_pdf = fullfile(fig_dir, ...
    'Step05_LegacyCandidateCompare_B1_S123_20251222.pdf');
exportgraphics(fig, out_png, 'Resolution', 220);
exportgraphics(fig, out_pdf, 'ContentType', 'vector');
close(fig);
end


function R = load_result_local(file)
S = load(file, 'Result');
if isfield(S, 'Result')
    R = S.Result;
else
    names = fieldnames(S);
    R = S.(names{1});
end
end


function T = get_seed_table_local(W, kind)
T = table();
if strcmpi(kind, "legacy")
    if isfield(W, 'seed_table')
        T = W.seed_table;
    end
else
    if isfield(W, 'SeedTable')
        T = W.SeedTable;
    elseif isfield(W, 'Result') && isfield(W.Result, 'VPSeedTable')
        T = W.Result.VPSeedTable;
    end
end
end


function T = normalize_seed_table_local(T)
if isstruct(T)
    T = struct2table(T);
end
if isempty(T)
    return;
end
if any(strcmp(T.Properties.VariableNames, 'weighted_voltage_rmse'))
    T = sortrows(T, {'weighted_voltage_rmse', 'plain_voltage_rmse', 'EO'}, ...
        {'ascend', 'ascend', 'ascend'});
end
end


function B = get_bundle_like_local(W, kind)
B = [];
if strcmpi(kind, "legacy")
    if isfield(W, 'bundle')
        B = W.bundle;
    end
else
    if isfield(W, 'CoreBundlePreview')
        B = W.CoreBundlePreview;
    elseif isfield(W, 'BundlePreview')
        B = W.BundlePreview;
    end
end
end


function [sid, w, qg, pc, vseg] = unpack_bundle_fields_local(B, kind)
if strcmpi(kind, "legacy")
    sid = B.S(:);
    w = B.W(:);
    qg = getfield_default_local(B, 'query_guard_mm', NaN);
    pc = getfield_default_local(B, 'point_count', numel(sid));
    vseg = getfield_default_local(B, 'valid_segment_count', NaN);
else
    sid = B.sensor_id(:);
    w = B.fit_weight(:);
    qg = mean(getfield_default_local(B, 'query_guard_mm', NaN), 'omitnan');
    pc = getfield_default_local(B, 'point_count', numel(sid));
    vseg = getfield_default_local(B, 'valid_pass_count', NaN);
end
end


function v = get_table_numeric_local(T, name)
if any(strcmp(T.Properties.VariableNames, name))
    v = T.(name);
else
    v = NaN(height(T), 1);
end
end


function s = get_table_string_local(T, name, default_value)
if any(strcmp(T.Properties.VariableNames, name))
    s = string(T.(name));
else
    s = default_value;
end
end


function v = table_value_local(T, row, names)
v = NaN;
for i = 1:numel(names)
    if any(strcmp(T.Properties.VariableNames, names{i}))
        raw = T.(names{i});
        if isnumeric(raw) || islogical(raw)
            v = raw(row);
        end
        return;
    end
end
end


function v = getfield_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
    v = S.(field_name);
else
    v = default_value;
end
end
