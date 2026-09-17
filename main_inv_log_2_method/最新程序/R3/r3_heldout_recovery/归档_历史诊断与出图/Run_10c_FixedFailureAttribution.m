function Result = Run_10c_FixedFailureAttribution()
%RUN_10C_FIXEDFAILUREATTRIBUTION Attribute fixed-method errors in formal rows.

root=fileparts(mfilename('fullpath'));
T=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
T=T(T.method=="fixed" & T.A_true_mm>0,:);
T.eo_fail=T.eo_est~=T.eo_true;
T.amplitude_fail=abs(T.amplitude_error_mm)>T.amplitude_tolerance_mm;
T.phase_fail=T.phase_error_rad>deg2rad(10);
T.failure_pattern=repmat("none",height(T),1);
for i=1:height(T)
    labels=["EO" "amplitude" "phase"];
    labels=labels([T.eo_fail(i) T.amplitude_fail(i) T.phase_fail(i)]);
    if ~isempty(labels),T.failure_pattern(i)=strjoin(labels," + ");end
end
patterns=unique(T.failure_pattern);gaps=unique(T.g_truth_mm);
byGap=table();
for i=1:numel(gaps)
    idxGap=T.g_truth_mm==gaps(i);
    for j=1:numel(patterns)
        count=sum(idxGap&T.failure_pattern==patterns(j));
        if count>0
            byGap=[byGap;table(gaps(i),patterns(j),count,count/sum(idxGap),...
                'VariableNames',{'g_truth_mm','failure_pattern','count','fraction'})]; %#ok<AGROW>
        end
    end
end
overall=table();
for j=1:numel(patterns)
    count=sum(T.failure_pattern==patterns(j));
    overall=[overall;table(patterns(j),count,count/height(T),...
        'VariableNames',{'failure_pattern','count','fraction'})]; %#ok<AGROW>
end
outDir=fullfile(root,'output','10_fixed_failure_attribution');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(byGap,fullfile(outDir,'fixed_failure_by_gap.csv'));
writetable(overall,fullfile(outDir,'fixed_failure_overall.csv'));
writetable(T,fullfile(outDir,'fixed_nonzero_rows_annotated.csv'));
Result=struct('byGap',byGap,'overall',overall,'rowCount',height(T),'outputDir',outDir);
save(fullfile(outDir,'fixed_failure_attribution.mat'),'Result');
fprintf('FIXED_FAILURE_ATTRIBUTION_OK rows=%d\n',height(T));
end
