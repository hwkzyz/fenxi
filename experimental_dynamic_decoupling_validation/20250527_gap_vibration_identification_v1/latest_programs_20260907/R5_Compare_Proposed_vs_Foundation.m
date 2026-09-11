function Comparison = R5_Compare_Proposed_vs_Foundation(resultFile)
%R5_COMPARE_PROPOSED_VS_FOUNDATION Formal aligned method comparison.
% Foundation is the previous direct low-speed-template method without
% dynamic clearance variation. Proposed R5 is the frozen PC1-anchor,
% sensor-wise dynamic-clearance, full-wave joint identification method.

if nargin < 1 || isempty(resultFile) || exist(resultFile,'file') ~= 2
    error('R5:ComparisonInput','Pass the exact formal Main06 result file.');
end
S = load(resultFile,'Result');
assert(isfield(S,'Result'),'R5:ComparisonInput','MAT file lacks Result.');
R = S.Result;
assert(isfield(R,'method') && ...
    strcmp(string(R.method),'R5_PC1_AnchorSurface_SensorWiseDg_FullWave'), ...
    'R5:ComparisonMethod','Input is not the frozen proposed R5 method.');
assert(isfield(R,'foundationStep05File') && exist(R.foundationStep05File,'file')==2, ...
    'R5:ComparisonFoundation','Bound Foundation result is missing.');
F = load(R.foundationStep05File,'Result');
assert(isfield(F,'Result') && isfield(F.Result,'Trend'), ...
    'R5:ComparisonFoundation','Foundation MAT lacks Result.Trend.');

P = R.Trend;
B = F.Result.Trend;
requiredP = {'window_id','lap_start','lap_end','gap_EO','gap_frequency_hz', ...
    'gap_amplitude_mm','gap_dx_mm','gap_rmse_mV'};
requiredB = {'window_id','lap_start','lap_end','EO_id','fn_id','A_id', ...
    'dx_c_id','plain_voltage_rmse'};
assert(all(ismember(requiredP,P.Properties.VariableNames)), ...
    'R5:ComparisonSchema','Proposed result lacks comparison fields.');
assert(all(ismember(requiredB,B.Properties.VariableNames)), ...
    'R5:ComparisonSchema','Foundation result lacks comparison fields.');
assert(height(P)==height(B) && isequal(double(P.window_id),double(B.window_id)) && ...
    isequal(double(P.lap_start),double(B.lap_start)) && ...
    isequal(double(P.lap_end),double(B.lap_end)), ...
    'R5:ComparisonAlignment','Foundation and proposed windows are not exactly aligned.');

window_id = double(P.window_id);
lap_start = double(P.lap_start);
lap_end = double(P.lap_end);
foundation_EO = double(B.EO_id);
proposed_EO = double(P.gap_EO);
foundation_frequency_hz = double(B.fn_id);
proposed_frequency_hz = double(P.gap_frequency_hz);
foundation_amplitude_mm = double(B.A_id);
proposed_amplitude_mm = double(P.gap_amplitude_mm);
foundation_dx_mm = double(B.dx_c_id);
proposed_dx_mm = double(P.gap_dx_mm);
foundation_rmse_mV = 1000*double(B.plain_voltage_rmse);
proposed_rmse_mV = double(P.gap_rmse_mV);
rmse_reduction_mV = foundation_rmse_mV-proposed_rmse_mV;
rmse_reduction_percent = 100*rmse_reduction_mV./foundation_rmse_mV;
eo_agreement = foundation_EO==proposed_EO;

WindowComparison = table(window_id,lap_start,lap_end,foundation_EO,proposed_EO, ...
    foundation_frequency_hz,proposed_frequency_hz,foundation_amplitude_mm, ...
    proposed_amplitude_mm,foundation_dx_mm,proposed_dx_mm,foundation_rmse_mV, ...
    proposed_rmse_mV,rmse_reduction_mV,rmse_reduction_percent,eo_agreement);

method = ["Foundation_direct_low_speed_template_no_dynamic_clearance"; ...
    "R5_PC1_AnchorSurface_SensorWiseDg_FullWave"];
dominant_EO = [mode_finite_local(foundation_EO);mode_finite_local(proposed_EO)];
EO_dominant_fraction = [mean(foundation_EO==dominant_EO(1),'omitnan'); ...
    mean(proposed_EO==dominant_EO(2),'omitnan')];
mean_frequency_hz = [mean(foundation_frequency_hz,'omitnan');mean(proposed_frequency_hz,'omitnan')];
std_frequency_hz = [std(foundation_frequency_hz,'omitnan');std(proposed_frequency_hz,'omitnan')];
mean_amplitude_mm = [mean(foundation_amplitude_mm,'omitnan');mean(proposed_amplitude_mm,'omitnan')];
mean_rmse_mV = [mean(foundation_rmse_mV,'omitnan');mean(proposed_rmse_mV,'omitnan')];
median_rmse_mV = [median(foundation_rmse_mV,'omitnan');median(proposed_rmse_mV,'omitnan')];
Summary = table(method,dominant_EO,EO_dominant_fraction,mean_frequency_hz, ...
    std_frequency_hz,mean_amplitude_mm,mean_rmse_mV,median_rmse_mV);

outDir = fullfile(fileparts(resultFile),'comparison');
if exist(outDir,'dir')~=7, mkdir(outDir); end
blade = double(R.cfg.targetBlade);
sensorTag = ['S',sprintf('%d',double(R.cfg.analysisSensors))];
regionTag = string(R.flowConfig.identification.resonanceRegionShortTag);
stem = sprintf('R5_FormalMethodComparison_B%d_%s_%s',blade,sensorTag,regionTag);
windowCsv = fullfile(outDir,[stem '_windows.csv']);
summaryCsv = fullfile(outDir,[stem '_summary.csv']);
matFile = fullfile(outDir,[stem '.mat']);
writetable(WindowComparison,windowCsv);
writetable(Summary,summaryCsv);
Comparison = struct('resultFile',resultFile,'foundationFile',R.foundationStep05File, ...
    'windowComparison',WindowComparison,'summary',Summary, ...
    'windowCsv',windowCsv,'summaryCsv',summaryCsv,'matFile',matFile, ...
    'eoAgreementFraction',mean(eo_agreement,'omitnan'), ...
    'meanRmseReductionMv',mean(rmse_reduction_mV,'omitnan'));
save(matFile,'Comparison','-v7.3');
disp(Summary);
fprintf('Formal aligned comparison saved:\n  %s\n  %s\n',windowCsv,summaryCsv);
end

function value = mode_finite_local(x)
x = x(isfinite(x));
if isempty(x), value=NaN; else, value=mode(x); end
end
