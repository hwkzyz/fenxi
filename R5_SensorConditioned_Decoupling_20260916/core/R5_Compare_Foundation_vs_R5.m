function comparison = R5_Compare_Foundation_vs_R5(r5ResultFile, outputBase)
%R5_COMPARE_FOUNDATION_VS_R5 Compare fixed-gap baseline with R5 per window.
% The comparison is descriptive. It never changes R5 acceptance.
if nargin < 1 || isempty(r5ResultFile)
    error('R5:MissingResult','An R5 result MAT file is required.');
end
if nargin < 2 || isempty(outputBase)
    outputBase = fullfile(fileparts(r5ResultFile), 'foundation_vs_r5_comparison');
end
S = load(r5ResultFile,'out');
assert(isfield(S,'out') && isfield(S.out,'rows'), ...
    'R5:ResultSchema','The result must contain out.rows.');
rows = S.out.rows;
n = numel(rows);
window_id = nan(n,1); EO_foundation = nan(n,1); f_foundation_hz = nan(n,1);
A_foundation_mm = nan(n,1); dx_foundation_mm = nan(n,1);
rmse_foundation_mv = nan(n,1); EO_R5 = nan(n,1); f_R5_hz = nan(n,1);
A_R5_mm = nan(n,1); dx_R5_mm = nan(n,1); rmse_R5_mv = nan(n,1);
rmse_reduction_mv = nan(n,1); rmse_reduction_fraction = nan(n,1);
delta_gap_R5_mm = cell(n,1); support_fraction_R5 = nan(n,1);
quality_class = strings(n,1); replay_status = strings(n,1);
for k = 1:n
    r = rows(k);
    window_id(k) = field_num(r,'window_id');
    EO_foundation(k) = field_num(r,'foundation_baseline_EO');
    f_foundation_hz(k) = field_num(r,'foundation_baseline_frequency_hz');
    A_foundation_mm(k) = field_num(r,'foundation_baseline_amplitude_mm');
    dx_foundation_mm(k) = field_num(r,'foundation_baseline_dx_mm');
    rmse_foundation_mv(k) = field_num(r,'baseline_rmse_mv');
    EO_R5(k) = field_num(r,'EO'); f_R5_hz(k) = field_num(r,'frequency_hz');
    A_R5_mm(k) = field_num(r,'amplitude_mm'); dx_R5_mm(k) = field_num(r,'dx_mm');
    rmse_R5_mv(k) = field_num(r,'rmse_mv');
    delta_gap_R5_mm{k} = field_value(r,'delta_gap_mm',[]);
    support_fraction_R5(k) = field_num(r,'support_fraction');
    quality_class(k) = string(field_value(r,'quality_class',''));
    replay_status(k) = string(field_value(r,'foundation_replay_status',''));
    if isfinite(rmse_foundation_mv(k)) && isfinite(rmse_R5_mv(k))
        rmse_reduction_mv(k) = rmse_foundation_mv(k) - rmse_R5_mv(k);
        rmse_reduction_fraction(k) = rmse_reduction_mv(k) / max(rmse_foundation_mv(k),eps);
    end
end
comparison = struct();
comparison.schema = 'R5_FOUNDATION_VS_R5_COMPARISON_V1';
comparison.sourceResult = r5ResultFile;
comparison.baselineDefinition = 'frozen_static_state_with_delta_g_zero';
comparison.proposedDefinition = 'frozen_static_state_with_sensorwise_delta_g';
comparison.replayIsDiagnosticOnly = true;
comparison.table = table(window_id,EO_foundation,f_foundation_hz,A_foundation_mm, ...
    dx_foundation_mm,rmse_foundation_mv,EO_R5,f_R5_hz,A_R5_mm,dx_R5_mm, ...
    rmse_R5_mv,rmse_reduction_mv,rmse_reduction_fraction,delta_gap_R5_mm, ...
    support_fraction_R5,quality_class,replay_status);
comparison.summary = struct();
comparison.summary.windowCount = n;
comparison.summary.meanFoundationRmseMv = mean(rmse_foundation_mv,'omitnan');
comparison.summary.meanR5RmseMv = mean(rmse_R5_mv,'omitnan');
comparison.summary.meanRmseReductionMv = mean(rmse_reduction_mv,'omitnan');
comparison.summary.meanRmseReductionFraction = mean(rmse_reduction_fraction,'omitnan');
comparison.summary.meanFoundationFrequencyHz = mean(f_foundation_hz,'omitnan');
comparison.summary.meanR5FrequencyHz = mean(f_R5_hz,'omitnan');
comparison.summary.meanFoundationAmplitudeMm = mean(A_foundation_mm,'omitnan');
comparison.summary.meanR5AmplitudeMm = mean(A_R5_mm,'omitnan');
comparison.summary.r5ImprovesWindowCount = nnz(rmse_reduction_mv > 0);
comparison.summary.r5WorsensWindowCount = nnz(rmse_reduction_mv < 0);
comparison.summary.r5EqualWindowCount = nnz(rmse_reduction_mv == 0);
if ~exist(fileparts(outputBase),'dir'), mkdir(fileparts(outputBase)); end
save([outputBase '.mat'],'comparison','-v7.3');
writetable(comparison.table,[outputBase '.csv']);
fprintf('Foundation vs R5 comparison written: %s.csv\n',outputBase);
end

function v = field_num(s,name)
v = NaN;
if isfield(s,name) && isnumeric(s.(name)) && isscalar(s.(name)), v = double(s.(name)); end
end
function v = field_value(s,name,d)
v = d;
if isfield(s,name), v = s.(name); end
end
