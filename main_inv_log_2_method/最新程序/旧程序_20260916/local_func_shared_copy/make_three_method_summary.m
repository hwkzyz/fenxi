function T = make_three_method_summary(R, truth)
%make_three_method_summary  Build one-row summary for three-method comparison.

[fTrue, order] = sort(truth.f(:).');
ATrue = truth.A(order);
deltaJ = get_field_or_nan(R, 'deltaJ');
rhoJ = get_field_or_nan(R, 'rhoJ');
T = table(R.method, R.g_used, abs(R.g_used - truth.g_high), ...
    R.f_id(1), R.f_id(2), R.A_id(1), R.A_id(2), ...
    mean(abs(R.f_id - fTrue)), mean(abs(R.A_id - ATrue)), ...
    R.rmse, R.eta_g, deltaJ, rhoJ, R.elapsed_s, ...
    'VariableNames', {'method','gap_used_mm','gap_error_mm', ...
    'f1_Hz','f2_Hz','A1_mm','A2_mm','mean_freq_error_Hz', ...
    'mean_amp_error_mm','rmse','eta_g','deltaJ','rhoJ','elapsed_s'});
end

function val = get_field_or_nan(S, name)
if isfield(S, name)
    val = S.(name);
else
    val = NaN;
end
end
