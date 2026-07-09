clc; clear; close all;

%RUN_STEP05_METHODMATRIX_SELECTEDBLADES
% Run three Step05 identification presets on the three selected
% blade/sensor objects:
%   20241106 B4/S257
%   20250527 B1/S136
%   20251222 B1/S123

root_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(root_dir, 'output', 'step05_method_matrix_selected');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

targets = [
    make_target_local("20241106", "20241106_btt_data_foundation", ...
        "Step05_SingleSync_DirectTemplate_Identification_20241106", ...
        "3000_3150", 4, "S257")
    make_target_local("20250527", "20250527_btt_data_foundation", ...
        "Step05_SingleSync_DirectTemplate_Identification_20250527", ...
        "20250526_2500-3500_t400", 1, "S136")
    make_target_local("20251222", "20251222_btt_data_foundation", ...
        "Step05_SingleSync_DirectTemplate_Identification_20251222", ...
        "1000_2500_3500", 1, "S123")
    ];

methods = [
    make_method_local("vp_top3_bounded_rmse", "VP top-3 + bounded RMSE")
    make_method_local("all_eo_unconstrained_rmse", "all EO + unconstrained RMSE")
    make_method_local("all_eo_unconstrained_corr", "all EO + unconstrained corr")
    ];

run_rows = table();
for it = 1:numel(targets)
    target = targets(it);
    folder = fullfile(root_dir, char(target.folder));

    for im = 1:numel(methods)
        method = methods(im);
        tag = sprintf('matrix_%s_%s', target.dataset, method.name);
        tag = regexprep(tag, '[^\w\-]', '_');
        if strlength(tag) > 24
            tag = extractBefore(tag, 25);
        end

        fprintf('\n=== Running %s %s with %s ===\n', ...
            target.dataset, target.object_label, method.name);

        old_dir = pwd;
        cleanup = onCleanup(@() restore_env_and_dir_local(old_dir));
        cd(folder);
        setenv('STEP05_METHOD_PRESET', char(method.name));
        setenv('STEP05_OUTPUT_TAG', char(tag));
        setenv('STEP05_SHOW_PLOTS', '0');
        setenv('STEP05_SAVE_FIGURES', '0');
        setenv('STEP05_FORCE_REBUILD', '1');

        if target.dataset == "20251222"
            setenv('STEP05_TARGET_BLADES', char(string(target.blade_id)));
            setenv('STEP05_TARGET_LAPS', '20');
            setenv('STEP05_WINDOW_LAPS', '3');
            setenv('STEP05_SLIDING_STEP_LAPS', '1');
        end

        status = "ok";
        message = "";
        try
            run(char(target.script));
        catch ME
            status = "failed";
            message = string(ME.message);
            fprintf(2, 'FAILED: %s\n', ME.message);
        end

        clear cleanup;

        run_rows = [run_rows; table( ...
            target.dataset, target.object_label, method.name, method.description, ...
            string(tag), status, message, ...
            'VariableNames', {'dataset','object_label','method_tag', ...
            'method_description','output_tag','run_status','message'})]; %#ok<AGROW>
    end
end

run_log_file = fullfile(out_dir, 'Step05_MethodMatrix_RunLog.csv');
writetable(run_rows, run_log_file);
fprintf('\nRun log saved: %s\n', run_log_file);

summary = summarize_matrix_outputs_local(root_dir, targets, methods);
summary_file = fullfile(out_dir, 'Step05_MethodMatrix_SelectedBlades_Summary.csv');
writetable(summary, summary_file);
fprintf('Summary saved: %s\n', summary_file);
disp(summary);

function target = make_target_local(dataset, folder, script, case_name, blade_id, sensor_tag)
target = struct();
target.dataset = string(dataset);
target.folder = string(folder);
target.script = string(script);
target.case_name = string(case_name);
target.blade_id = blade_id;
target.sensor_tag = string(sensor_tag);
target.object_label = sprintf('%s B%d/%s', dataset, blade_id, sensor_tag);
end

function method = make_method_local(name, description)
method = struct();
method.name = string(name);
method.description = string(description);
end

function restore_env_and_dir_local(old_dir)
cd(old_dir);
setenv('STEP05_METHOD_PRESET', '');
setenv('STEP05_OUTPUT_TAG', '');
setenv('STEP05_SHOW_PLOTS', '');
setenv('STEP05_SAVE_FIGURES', '');
setenv('STEP05_FORCE_REBUILD', '');
setenv('STEP05_TARGET_BLADES', '');
setenv('STEP05_TARGET_LAPS', '');
setenv('STEP05_WINDOW_LAPS', '');
setenv('STEP05_SLIDING_STEP_LAPS', '');
end

function summary = summarize_matrix_outputs_local(root_dir, targets, methods)
summary = table();
for it = 1:numel(targets)
    target = targets(it);
    for im = 1:numel(methods)
        method = methods(im);
        tag = sprintf('matrix_%s_%s', target.dataset, method.name);
        tag = regexprep(tag, '[^\w\-]', '_');
        if strlength(tag) > 24
            tag = extractBefore(tag, 25);
        end
        step05_dir = fullfile(root_dir, char(target.folder), 'output', ...
            ['step05_single_sync_direct_template_', char(tag)], ...
            char(target.case_name));
        trend_file = fullfile(step05_dir, sprintf( ...
            'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
            target.blade_id, target.sensor_tag, target.dataset));

        if exist(trend_file, 'file') ~= 2
            row = make_summary_row_local(target, method, tag, "missing", ...
                NaN, NaN, NaN, NaN, NaN, NaN, "", trend_file);
        else
            T = readtable(trend_file, 'TextType', 'string');
            row = summarize_trend_file_local(T, target, method, tag, trend_file);
        end
        summary = [summary; row]; %#ok<AGROW>
    end
end
end

function row = summarize_trend_file_local(T, target, method, tag, trend_file)
status = "ok";
if any(strcmp(T.Properties.VariableNames, 'status'))
    ok = strcmpi(string(T.status), "ok");
else
    ok = true(height(T), 1);
end
if any(strcmp(T.Properties.VariableNames, 'EO_id'))
    ok = ok & isfinite(T.EO_id);
end

Tok = T(ok, :);
if height(Tok) == 0
    reason = "";
    if any(strcmp(T.Properties.VariableNames, 'failure_reason')) && height(T) > 0
        reason = strjoin(unique(string(T.failure_reason)), " | ");
    end
    row = make_summary_row_local(target, method, tag, "no_ok_window", ...
        height(Tok), height(T), NaN, NaN, NaN, NaN, reason, trend_file);
    return;
end

eo = mode(Tok.EO_id);
freq = mean(Tok.fn_id, 'omitnan');
median_freq = median(Tok.fn_id, 'omitnan');
median_A = median(Tok.A_id, 'omitnan');
median_rmse = median(Tok.weighted_voltage_rmse, 'omitnan');
eo_unique = strjoin(string(unique(Tok.EO_id(:).', 'stable')), "|");
status = classify_result_local(Tok, eo_unique);

row = make_summary_row_local(target, method, tag, status, ...
    height(Tok), height(T), eo, freq, median_freq, median_A, ...
    median_rmse, eo_unique, trend_file);
end

function status = classify_result_local(Tok, eo_unique)
unique_count = numel(strsplit(char(eo_unique), '|'));
if unique_count <= 1
    status = "stable_single_eo";
elseif unique_count <= 3
    status = "mixed_eo";
else
    status = "unstable_many_eo";
end
if any(strcmp(Tok.Properties.VariableNames, 'eo_rmse_gap_ratio'))
    gap = median(Tok.eo_rmse_gap_ratio, 'omitnan');
    if isfinite(gap) && gap < 0.03
        status = status + "_low_gap";
    end
end
end

function row = make_summary_row_local(target, method, tag, result_status, ...
    ok_windows, total_windows, eo, mean_freq, median_freq, median_A, ...
    median_rmse, note, source_file)
row = table( ...
    target.dataset, target.object_label, target.blade_id, target.sensor_tag, ...
    method.name, method.description, string(tag), string(result_status), ...
    double(ok_windows), double(total_windows), double(eo), double(mean_freq), ...
    double(median_freq), double(median_A), double(median_rmse), ...
    string(note), string(source_file), ...
    'VariableNames', {'dataset','object_label','blade_id','sensor_tag', ...
    'method_tag','method_description','output_tag','result_status', ...
    'ok_window_count','total_window_count','eo_mode','mean_freq_hz', ...
    'median_freq_hz','median_A_mm','median_rmse_v','note','source_file'});
end
