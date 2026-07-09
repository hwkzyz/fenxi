clc; clear; close all;

%COMPARE_STEP05_METHODRESULTS_ACROSSFOUNDATIONS
% Summarize existing Step05 identification outputs across the three
% foundation folders. This script does not rerun identification. It keeps
% formal long-window outputs separate from smoke/diagnostic outputs.

root_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(root_dir, 'output', 'cross_foundation_step05_compare');
fig_dir = fullfile(root_dir, 'output', 'figures', 'cross_foundation_step05_compare');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

specs = build_result_specs_local(root_dir);
rows = table();
missing = table();

for i = 1:numel(specs)
    spec = specs(i);
    file = fullfile(root_dir, spec.relative_file);
    if exist(file, 'file') ~= 2
        missing = [missing; struct2table(struct( ...
            'dataset', string(spec.dataset), ...
            'method_tag', string(spec.method_tag), ...
            'result_class', string(spec.result_class), ...
            'relative_file', string(spec.relative_file), ...
            'reason', "missing_file"))]; %#ok<AGROW>
        continue;
    end

    T = readtable(file, 'TextType', 'string');
    if spec.kind == "trend"
        R = summarize_trend_local(T, spec, file);
    else
        R = summarize_summary_local(T, spec, file);
    end
    rows = [rows; R]; %#ok<AGROW>
end

rows = sortrows(rows, {'dataset', 'result_class', 'method_tag', 'blade_id'});

summary_file = fullfile(out_dir, 'Step05_MethodComparison_CurrentOutputs.csv');
writetable(rows, summary_file);

missing_file = fullfile(out_dir, 'Step05_MethodComparison_MissingOutputs.csv');
writetable(missing, missing_file);

fprintf('\n=== Step05 method comparison from existing outputs ===\n');
disp(rows(:, {'dataset','method_tag','result_class','blade_id','sensor_tag', ...
    'window_count','eo_report','freq_report_hz','median_A_mm','median_rmse_v', ...
    'joint_status','interpretation'}));
fprintf('Saved: %s\n', summary_file);
if height(missing) > 0
    fprintf('Missing expected outputs: %s\n', missing_file);
end

plot_frequency_compare_local(rows, fig_dir);

function specs = build_result_specs_local(root_dir) %#ok<INUSD>
specs = repmat(struct( ...
    'dataset', "", ...
    'method_tag', "", ...
    'result_class', "", ...
    'kind', "", ...
    'blade_id', NaN, ...
    'sensor_tag', "", ...
    'relative_file', ""), 0, 1);

specs(end+1) = make_spec_local( ...
    "20241106", "vp_top3_bounded_rmse", "formal_default", "trend", 4, "S257", ...
    fullfile('20241106_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template', '3000_3150', ...
    'Trend_Step05_SingleSyncDirectTemplate_B4_S257_20241106.csv'));

specs(end+1) = make_spec_local( ...
    "20250527", "vp_top3_bounded_rmse", "formal_default", "trend", 1, "S136", ...
    fullfile('20250527_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template', '20250526_2500-3500_t400', ...
    'Trend_Step05_SingleSyncDirectTemplate_B1_S136_20250527.csv'));

% 20251222 formal long-window joint runs. B1/B2/B3 are the blades currently
% audited against the old 20251222 route.
specs(end+1) = make_spec_local( ...
    "20251222", "all_eo_unconstrained_corr_joint", "formal_long_window", ...
    "summary", 1, "S123", ...
    fullfile('20251222_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template_official_default_unconstrained_w18', ...
    'Step05_AllBlade_Summary_20251222.csv'));

specs(end+1) = make_spec_local( ...
    "20251222", "all_eo_unconstrained_corr_joint", "formal_long_window", ...
    "summary", 2, "S123", ...
    fullfile('20251222_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template_b2w18', ...
    'Step05_AllBlade_Summary_20251222.csv'));

specs(end+1) = make_spec_local( ...
    "20251222", "all_eo_unconstrained_corr_joint", "formal_long_window", ...
    "summary", 3, "S123", ...
    fullfile('20251222_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template_b3w18', ...
    'Step05_AllBlade_Summary_20251222.csv'));

% Existing 20251222 smoke scan. This is useful as a contrast, but it is not
% a formal frequency decision because it only uses three windows and lacks
% the long-window joint ambiguity gate.
specs(end+1) = make_spec_local( ...
    "20251222", "quick_smoke_existing", "diagnostic_smoke", ...
    "summary", NaN, "S123", ...
    fullfile('20251222_btt_data_foundation', 'output', ...
    'step05_single_sync_direct_template_smoke_allblades_w3_v2', ...
    'Step05_AllBlade_Summary_20251222.csv'));
end

function spec = make_spec_local(dataset, method_tag, result_class, kind, ...
    blade_id, sensor_tag, relative_file)
spec = struct();
spec.dataset = string(dataset);
spec.method_tag = string(method_tag);
spec.result_class = string(result_class);
spec.kind = string(kind);
spec.blade_id = blade_id;
spec.sensor_tag = string(sensor_tag);
spec.relative_file = char(relative_file);
end

function R = summarize_trend_local(T, spec, file)
ok = true(height(T), 1);
if any(strcmp(T.Properties.VariableNames, 'status'))
    ok = strcmpi(string(T.status), "ok");
end
if any(strcmp(T.Properties.VariableNames, 'EO_id'))
    ok = ok & isfinite(T.EO_id);
end

Tok = T(ok, :);
eo = NaN;
freq = NaN;
median_A = NaN;
median_rmse = NaN;
eo_sequence = "";
if height(Tok) > 0
    if any(strcmp(Tok.Properties.VariableNames, 'EO_id'))
        eo_values = Tok.EO_id(isfinite(Tok.EO_id));
        eo = mode(eo_values);
        eo_sequence = strjoin(string(unique(eo_values(:).', 'stable')), "|");
    end
    if any(strcmp(Tok.Properties.VariableNames, 'fn_id'))
        freq = mean(Tok.fn_id, 'omitnan');
    end
    if any(strcmp(Tok.Properties.VariableNames, 'A_id'))
        median_A = median(Tok.A_id, 'omitnan');
    end
    if any(strcmp(Tok.Properties.VariableNames, 'weighted_voltage_rmse'))
        median_rmse = median(Tok.weighted_voltage_rmse, 'omitnan');
    end
end

R = result_row_local(spec.dataset, spec.method_tag, spec.result_class, ...
    spec.blade_id, spec.sensor_tag, height(Tok), height(T), eo, eo_sequence, ...
    freq, median_A, median_rmse, "", "", file);
end

function R = summarize_summary_local(T, spec, file)
if isfinite(spec.blade_id) && any(strcmp(T.Properties.VariableNames, 'blade_id'))
    T = T(T.blade_id == spec.blade_id, :);
end

R = table();
for i = 1:height(T)
    row = T(i, :);
    blade_id = get_table_scalar_local(row, 'blade_id', spec.blade_id);
    sensor_tag = get_table_string_local(row, 'sensor_tag', spec.sensor_tag);
    windows = get_table_scalar_local(row, 'valid_window_count', NaN);
    total_windows = get_table_scalar_local(row, 'total_window_count', NaN);
    eo = get_table_scalar_local(row, 'dominant_eo', NaN);
    joint_eo = get_table_scalar_local(row, 'joint_eo', NaN);
    if ~isfinite(eo) && isfinite(joint_eo)
        eo = joint_eo;
    end
    if isfinite(joint_eo)
        eo_sequence = "joint:" + string(joint_eo);
    else
        eo_sequence = "";
    end

    freq = get_table_scalar_local(row, 'mean_freq_hz', NaN);
    joint_freq = get_table_scalar_local(row, 'joint_freq_hz', NaN);
    if ~isfinite(freq) && isfinite(joint_freq)
        freq = joint_freq;
    end
    median_A = get_table_scalar_local(row, 'median_A_mm', NaN);
    median_rmse = get_table_scalar_local(row, 'median_rmse_v', NaN);
    joint_status = get_table_string_local(row, 'joint_eo_status', "");

    Ri = result_row_local(spec.dataset, spec.method_tag, spec.result_class, ...
        blade_id, sensor_tag, windows, total_windows, eo, eo_sequence, ...
        freq, median_A, median_rmse, joint_status, "", file);
    R = [R; Ri]; %#ok<AGROW>
end
end

function R = result_row_local(dataset, method_tag, result_class, blade_id, ...
    sensor_tag, window_count, total_window_count, eo_report, eo_sequence, ...
    freq_report_hz, median_A_mm, median_rmse_v, joint_status, interpretation, file)
if strlength(interpretation) == 0
    interpretation = interpret_result_local(result_class, window_count, ...
        total_window_count, eo_report, freq_report_hz, joint_status);
end

R = table( ...
    string(dataset), string(method_tag), string(result_class), ...
    double(blade_id), string(sensor_tag), double(window_count), ...
    double(total_window_count), double(eo_report), string(eo_sequence), ...
    double(freq_report_hz), double(median_A_mm), double(median_rmse_v), ...
    string(joint_status), string(interpretation), string(file), ...
    'VariableNames', {'dataset','method_tag','result_class','blade_id', ...
    'sensor_tag','window_count','total_window_count','eo_report', ...
    'eo_sequence','freq_report_hz','median_A_mm','median_rmse_v', ...
    'joint_status','interpretation','source_file'});
end

function txt = interpret_result_local(result_class, window_count, ...
    total_window_count, eo_report, freq_report_hz, joint_status)
if result_class == "diagnostic_smoke"
    txt = "diagnostic only: short-window smoke result";
elseif strlength(joint_status) > 0 && strcmpi(joint_status, "ambiguous")
    txt = "ambiguous: candidate EO reported for diagnostics, final frequency should be NaN";
elseif strlength(joint_status) > 0 && strcmpi(joint_status, "ok")
    txt = "accepted by long-window joint EO gate";
elseif window_count > 0 && isfinite(eo_report) && isfinite(freq_report_hz)
    txt = "accepted by per-window VP top-k direct-template route";
elseif total_window_count > 0
    txt = "no accepted per-window fit";
else
    txt = "no usable result";
end
end

function x = get_table_scalar_local(T, name, default_value)
if any(strcmp(T.Properties.VariableNames, name))
    x = T.(name)(1);
else
    x = default_value;
end
if iscell(x)
    x = x{1};
end
if isstring(x) || ischar(x)
    x = str2double(x);
end
if isempty(x)
    x = default_value;
end
end

function s = get_table_string_local(T, name, default_value)
if any(strcmp(T.Properties.VariableNames, name))
    s = string(T.(name)(1));
else
    s = string(default_value);
end
if ismissing(s)
    s = string(default_value);
end
end

function plot_frequency_compare_local(rows, fig_dir)
if height(rows) == 0
    return;
end
plot_rows = rows(isfinite(rows.freq_report_hz), :);
if height(plot_rows) == 0
    return;
end

labels = plot_rows.dataset + " B" + string(plot_rows.blade_id) + newline + ...
    plot_rows.method_tag;
fig = figure('Name', 'Step05 cross-foundation frequency comparison', ...
    'Color', 'w', 'Position', [100 100 1300 520]);
bar(plot_rows.freq_report_hz, 'FaceColor', [0.20 0.45 0.75]);
grid on;
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
xtickangle(30);
ylabel('Reported frequency (Hz)');
title('Step05 identification results from existing outputs');
for i = 1:height(plot_rows)
    text(i, plot_rows.freq_report_hz(i), sprintf(' EO %.0f', plot_rows.eo_report(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8);
end
exportgraphics(fig, fullfile(fig_dir, 'Step05_MethodComparison_Frequency.png'), ...
    'Resolution', 200);
savefig(fig, fullfile(fig_dir, 'Step05_MethodComparison_Frequency.fig'));
end
