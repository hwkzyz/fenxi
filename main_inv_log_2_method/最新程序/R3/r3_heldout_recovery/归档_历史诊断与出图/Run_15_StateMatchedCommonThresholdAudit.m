function Result = Run_15_StateMatchedCommonThresholdAudit()
%RUN_15_STATEMATCHEDCOMMONTHRESHOLDAUDIT Apply common thresholds to SM-delta.

root=fileparts(mfilename('fullpath'));gaps=[0.5 0.7 0.8 0.9 1.0 1.2];All=table();
Formal=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
for g=gaps
    if g==0.8
        D=readtable(fullfile(root,'output','08_state_matched_coordinate_diagnostic',...
            'reference_coordinate_detail.csv'),'TextType','string');
    else
        tag=sprintf('g%04d',round(1000*g));
        D=readtable(fullfile(root,'output','09_state_matched_offreference_diagnostic',...
            [tag '_detail.csv']),'TextType','string');
    end
    D=D(D.coordinate_mode=="SM_delta"&D.A_true_mm==0,:);assert(height(D)==75);
    for i=1:height(D)
        F=Formal(Formal.g_truth_mm==g&Formal.A_true_mm==0&Formal.high_seed==D.high_seed(i),:);
        assert(height(F)==3);
        D.threshold_fixed_mm(i)=F.detection_threshold_mm(F.method=="fixed");
        D.threshold_adaptive_mm(i)=F.detection_threshold_mm(F.method=="adaptive");
        D.threshold_state_matched_mm(i)=F.detection_threshold_mm(F.method=="state_matched");
    end
    D.fp_fixed_threshold=D.A_est_mm>D.threshold_fixed_mm;
    D.fp_adaptive_threshold=D.A_est_mm>D.threshold_adaptive_mm;
    D.fp_state_matched_threshold=D.A_est_mm>D.threshold_state_matched_mm;
    All=[All;D]; %#ok<AGROW>
end
Summary=table();
for g=gaps
    Q=All(All.g_truth_mm==g,:);
    Summary=[Summary;table(g,height(Q),mean(Q.fp_fixed_threshold),mean(Q.fp_adaptive_threshold),...
        mean(Q.fp_state_matched_threshold),median(Q.A_est_mm),max(Q.A_est_mm),...
        'VariableNames',{'g_truth_mm','count','PFP_fixed_threshold','PFP_adaptive_threshold',...
        'PFP_state_matched_threshold','median_A_est_mm','max_A_est_mm'})]; %#ok<AGROW>
end
outDir=fullfile(root,'output','15_state_matched_common_threshold_audit');if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(All,fullfile(outDir,'state_matched_delta_null_common_threshold_detail.csv'));
writetable(Summary,fullfile(outDir,'state_matched_delta_null_common_threshold_summary.csv'));
Result=struct('Detail',All,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,'state_matched_common_threshold_audit.mat'),'Result');
fprintf('STATE_MATCHED_COMMON_THRESHOLD_AUDIT_OK rows=%d\n',height(All));
end
