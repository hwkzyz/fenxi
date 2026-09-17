function Result = Run_10b_FixedFailureAttribution()
%RUN_10B_FIXEDFAILUREATTRIBUTION Attribute fixed-method errors in formal rows.

root=fileparts(mfilename('fullpath'));
T=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
T=T(T.method=="fixed" & T.A_true_mm>0,:);
T.eo_fail=T.eo_est~=T.eo_true;
T.amplitude_fail=abs(T.amplitude_error_mm)>T.amplitude_tolerance_mm;
T.phase_fail=T.phase_error_rad>deg2rad(10);
T.failure_pattern=strings(height(T),1);
for i=1:height(T)
    labels=["EO" "amplitude" "phase"];
    labels=labels([T.eo_fail(i) T.amplitude_fail(i) T.phase_fail(i)]);
    T.failure_pattern(i)=strjoin(labels," + ");
end
Summary=groupsummary(T,{'g_truth_mm','failure_pattern'},'numel','eo_est');
Summary.Properties.VariableNames{'GroupCount'}='count';
Summary.Properties.VariableNames{'numel_eo_est'}=[];
GapTotals=groupsummary(T,'g_truth_mm','numel','eo_est');
Summary.fraction=zeros(height(Summary),1);
for i=1:height(Summary)
    Summary.fraction(i)=Summary.count(i)/GapTotals.numel_eo_est(GapTotals.g_truth_mm==Summary.g_truth_mm(i));
end
Overall=groupsummary(T,'failure_pattern','numel','eo_est');
Overall.Properties.VariableNames{'GroupCount'}='count';
Overall.Properties.VariableNames{'numel_eo_est'}=[];
Overall.fraction=Overall.count/height(T);
outDir=fullfile(root,'output','10_fixed_failure_attribution');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Summary,fullfile(outDir,'fixed_failure_by_gap.csv'));
writetable(Overall,fullfile(outDir,'fixed_failure_overall.csv'));
writetable(T,fullfile(outDir,'fixed_nonzero_rows_annotated.csv'));
Result=struct('Summary',Summary,'Overall',Overall,'rowCount',height(T),'outputDir',outDir);
save(fullfile(outDir,'fixed_failure_attribution.mat'),'Result');
fprintf('FIXED_FAILURE_ATTRIBUTION_OK rows=%d\n',height(T));
end
