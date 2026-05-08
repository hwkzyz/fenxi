function aggregateTable = aggregate_multi_case_table(detailTable)
%aggregate_multi_case_table  Aggregate multi-case identification errors.

methodList = unique(string(detailTable.method), 'stable');
rows = struct('method', {}, 'mean_gap_error_mm', {}, 'mean_f1_error_Hz', {}, ...
    'mean_f2_error_Hz', {}, 'mean_A1_error_mm', {}, 'mean_A2_error_mm', {}, ...
    'mean_freq_error_Hz', {}, 'mean_amp_error_mm', {}, 'mean_rmse', {}, ...
    'mean_eta_g', {}, 'mean_elapsed_s', {});

for im = 1:numel(methodList)
    m = methodList(im);
    idx = string(detailTable.method) == m;
    rows(end+1).method = m; %#ok<AGROW>
    rows(end).mean_gap_error_mm = mean(detailTable.gap_error_mm(idx));
    rows(end).mean_f1_error_Hz = mean(detailTable.f1_error_Hz(idx));
    rows(end).mean_f2_error_Hz = mean(detailTable.f2_error_Hz(idx));
    rows(end).mean_A1_error_mm = mean(detailTable.A1_error_mm(idx));
    rows(end).mean_A2_error_mm = mean(detailTable.A2_error_mm(idx));
    rows(end).mean_freq_error_Hz = mean(detailTable.mean_freq_error_Hz(idx));
    rows(end).mean_amp_error_mm = mean(detailTable.mean_amp_error_mm(idx));
    rows(end).mean_rmse = mean(detailTable.rmse(idx));
    rows(end).mean_eta_g = mean(detailTable.eta_g(idx));
    rows(end).mean_elapsed_s = mean(detailTable.elapsed_s(idx));
end

aggregateTable = struct2table(rows);
end
