function Result = Run_13_AdaptiveNullDistributionAudit()
%RUN_13_ADAPTIVENULLDISTRIBUTIONAUDIT Summarize formal adaptive null fits.

root=fileparts(mfilename('fullpath'));P=R3_Protocol();
S=load(fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat'),'Result');
g0=S.Result.C0.pathCal.g0;
T=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
T=T(T.method=="adaptive"&T.A_true_mm==0,:);
assert(height(T)==numel(P.gMainMm)*75,'Unexpected adaptive-null row count.');
T.expected_delta_coordinate_mm=g0+(T.g_truth_mm-P.gReferenceMm);
T.coordinate_error_mm=T.g_model_coordinate_mm-T.expected_delta_coordinate_mm;
gaps=unique(T.g_truth_mm);Gap=table();EO=table();
for i=1:numel(gaps)
    Q=T(T.g_truth_mm==gaps(i),:);assert(height(Q)==75);
    Gap=[Gap;table(gaps(i),height(Q),sum(Q.false_positive),mean(Q.false_positive),...
        median(Q.A_est_mm),max(Q.A_est_mm),median(Q.coordinate_error_mm),...
        min(Q.coordinate_error_mm),max(Q.coordinate_error_mm),...
        'VariableNames',{'g_truth_mm','count','false_positive_count','P_FP','median_A_est_mm',...
        'max_A_est_mm','median_coordinate_error_mm','min_coordinate_error_mm','max_coordinate_error_mm'})]; %#ok<AGROW>
    eos=unique(Q.eo_est);
    for j=1:numel(eos)
        E=Q(Q.eo_est==eos(j),:);
        EO=[EO;table(gaps(i),eos(j),height(E),mean(E.false_positive),...
            'VariableNames',{'g_truth_mm','eo_est','count','P_FP'})]; %#ok<AGROW>
    end
end
outDir=fullfile(root,'output','13_adaptive_null_distribution_audit');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Gap,fullfile(outDir,'adaptive_null_by_gap.csv'));
writetable(EO,fullfile(outDir,'adaptive_null_eo_distribution.csv'));
writetable(T,fullfile(outDir,'adaptive_null_rows_annotated.csv'));
Result=struct('Gap',Gap,'EO',EO,'Detail',T,'outputDir',outDir);
save(fullfile(outDir,'adaptive_null_distribution.mat'),'Result');
fprintf('ADAPTIVE_NULL_DISTRIBUTION_OK rows=%d\n',height(T));
end
