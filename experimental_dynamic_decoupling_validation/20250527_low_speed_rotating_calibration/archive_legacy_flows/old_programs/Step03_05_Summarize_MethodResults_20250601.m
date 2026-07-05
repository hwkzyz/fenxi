% Summarize current Step03 low-speed-template-only method results.
%
% This script does not rerun identification and does not change any Step03
% algorithm. It only reads the current result MAT files plus the measured
% runtime CSV, then exports compact method summaries and complete per-window
% vibration parameters.

clear; clc;

root_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(root_dir, 'output', 'identification');
runtime_csv = fullfile(result_dir, 'Step03_Method_RunSummary_20250601.csv');

methods = {
    'Step03_01', 'DirectLowSpeedWaveform',      '无加速；全 EO 完整波形优化';
    'Step03_02', 'TemplateInverseSeedTopK',     '模板反查种子；top-K 后完整波形优化';
    'Step03_03', 'FirstOrderVPTopK',            '一阶模板 VP；top-K 后完整波形优化';
    'Step03_04', 'FirstOrderVPAdaptive',        '一阶模板 VP；严格 top-3；同步频率；优化 dx_c'
};
sensors = {'S13', 'S136'};

if exist(runtime_csv, 'file') ~= 2
    error('Runtime summary CSV not found: %s', runtime_csv);
end
runtime_tbl = readtable(runtime_csv, 'TextType', 'string');

summary_rows = table();
detail_rows = table();

for im = 1:size(methods, 1)
    method_id = string(methods{im, 1});
    method_label = string(methods{im, 2});
    method_note = string(methods{im, 3});

    for is = 1:numel(sensors)
        sensor_tag = string(sensors{is});
        mat_file = fullfile(result_dir, sprintf('Result_%s_%s_B1_%s_20250527.mat', ...
            method_id, method_label, sensor_tag));
        if exist(mat_file, 'file') ~= 2
            error('Result MAT file not found: %s', mat_file);
        end

        S = load(mat_file, 'Result');
        R = S.Result;
        T = R.Trend;
        if ismember('dx_c_id', T.Properties.VariableNames)
            dx_c_vals = T.dx_c_id;
        else
            dx_c_vals = T.d0_id;
        end
        nwin = height(T);

        candidate_counts = nan(nwin, 1);
        selected_counts = nan(nwin, 1);
        fallback_flags = false(nwin, 1);
        for iw = 1:nwin
            W = R.WindowResult(iw);
            if isfield(W, 'Result') && isfield(W.Result, 'CandidateTable')
                candidate_counts(iw) = height(W.Result.CandidateTable);
            end
            if isfield(W, 'Result') && isfield(W.Result, 'VPSelectedEO')
                selected_counts(iw) = numel(W.Result.VPSelectedEO);
            else
                selected_counts(iw) = candidate_counts(iw);
            end
            if isfield(W, 'Result') && isfield(W.Result, 'VPSelectionInfo')
                info = W.Result.VPSelectionInfo;
                if isfield(info, 'mode')
                    fallback_flags(iw) = contains(string(info.mode), "fallback") || contains(string(info.mode), "warmup");
                elseif isfield(info, 'fallback_reason')
                    fallback_flags(iw) = strlength(string(info.fallback_reason)) > 0;
                end
            end
        end

        rt_mask = runtime_tbl.method == method_id & runtime_tbl.sensor == sensor_tag;
        if nnz(rt_mask) ~= 1
            error('Runtime row not found or duplicated for %s %s.', method_id, sensor_tag);
        end
        rt = runtime_tbl(rt_mask, :);

        summary_rows = [summary_rows; table( ...
            method_id, sensor_tag, method_note, nwin, rt.elapsed_s, rt.sec_per_window, ...
            mode(T.EO_id), all(T.EO_id == 14), mean(T.fn_id), median(T.fn_id), std(T.fn_id), ...
            mean(T.A_id), median(T.A_id), std(T.A_id), ...
            mean(T.phi_id_wrapped), median(T.phi_id_wrapped), std(T.phi_id_wrapped), ...
            mean(dx_c_vals), median(dx_c_vals), std(dx_c_vals), ...
            mean(T.weighted_voltage_rmse), median(T.weighted_voltage_rmse), std(T.weighted_voltage_rmse), ...
            mean(candidate_counts, 'omitnan'), median(candidate_counts, 'omitnan'), ...
            mean(selected_counts, 'omitnan'), median(selected_counts, 'omitnan'), nnz(fallback_flags), ...
            'VariableNames', {'method','sensor','strategy','windows','elapsed_s','sec_per_window', ...
            'dominant_eo','all_windows_eo14','mean_f_hz','median_f_hz','std_f_hz', ...
            'mean_A_mm','median_A_mm','std_A_mm', ...
            'mean_phi_rad','median_phi_rad','std_phi_rad', ...
            'mean_dx_c_mm','median_dx_c_mm','std_dx_c_mm', ...
            'mean_rmse_V','median_rmse_V','std_rmse_V', ...
            'mean_refined_candidates','median_refined_candidates', ...
            'mean_selected_candidates','median_selected_candidates','fallback_or_warmup_windows'})];

        detail_rows = [detail_rows; table( ...
            repmat(method_id, nwin, 1), repmat(sensor_tag, nwin, 1), T.window_id, ...
            T.lap_start, T.lap_end, T.rot_freq_mean_hz, T.EO_id, T.fn_id, ...
            T.A_id, T.phi_id_wrapped, dx_c_vals, T.weighted_voltage_rmse, ...
            T.plain_voltage_rmse, T.valid_segment_count, T.point_count, ...
            candidate_counts, selected_counts, fallback_flags, ...
            'VariableNames', {'method','sensor','window_id','lap_start','lap_end', ...
            'rot_freq_mean_hz','EO','f_hz','A_mm','phi_rad','dx_c_mm', ...
            'weighted_rmse_V','plain_rmse_V','valid_segment_count','point_count', ...
            'refined_candidate_count','selected_candidate_count','fallback_or_warmup'})];
    end
end

summary_csv = fullfile(result_dir, 'Step03_Method_ResultSummary_20250601.csv');
detail_csv = fullfile(result_dir, 'Step03_Method_WindowParameters_20250601.csv');
writetable(summary_rows, summary_csv);
writetable(detail_rows, detail_csv);

fprintf('Wrote method summary: %s\n', summary_csv);
fprintf('Wrote per-window parameters: %s\n', detail_csv);
disp(summary_rows(:, {'method','sensor','windows','elapsed_s','sec_per_window', ...
    'dominant_eo','all_windows_eo14','mean_f_hz','mean_A_mm','mean_phi_rad', ...
    'mean_dx_c_mm','mean_rmse_V','mean_refined_candidates','fallback_or_warmup_windows'}));
