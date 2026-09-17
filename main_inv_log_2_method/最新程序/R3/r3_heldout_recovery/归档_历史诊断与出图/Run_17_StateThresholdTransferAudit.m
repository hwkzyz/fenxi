function Result = Run_17_StateThresholdTransferAudit()
%RUN_17_STATETHRESHOLDTRANSFERAUDIT Quantify gap-dependent null thresholds.

root=fileparts(mfilename('fullpath'));P=R3_Protocol();
Formal=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
T=readtable(fullfile(root,'output','16_corrected_state_matched_benchmark',...
    'corrected_state_matched_detail.csv'),'TextType','string');
T=T(T.A_true_mm==0,:);assert(height(T)==450);
gaps=unique(T.g_truth_mm);rows=repmat(empty_row(),numel(gaps),1);
for i=1:numel(gaps)
    g=gaps(i);Q=T(T.g_truth_mm==g,:);F=Formal(Formal.method=="fixed"&Formal.g_truth_mm==g&Formal.A_true_mm==0,:);
    A=Formal(Formal.method=="adaptive"&Formal.g_truth_mm==g&Formal.A_true_mm==0,:);
    SM=Formal(Formal.method=="state_matched"&Formal.g_truth_mm==g&Formal.A_true_mm==0,:);
    assert(height(Q)==75&&height(F)==1&&height(A)==1&&height(SM)==1);
    q95=quantile(Q.A_est_mm,0.95);q99=quantile(Q.A_est_mm,0.99);
    r=empty_row();r.g_truth_mm=g;r.count=height(Q);r.q95_sm_delta_mm=q95;r.q99_sm_delta_mm=q99;
    r.formal_fixed_threshold_mm=F.detection_threshold_mm;r.formal_adaptive_threshold_mm=A.detection_threshold_mm;
    r.formal_state_matched_threshold_mm=SM.detection_threshold_mm;
    r.PFP_at_formal_fixed=mean(Q.A_est_mm>r.formal_fixed_threshold_mm);
    r.PFP_at_formal_adaptive=mean(Q.A_est_mm>r.formal_adaptive_threshold_mm);
    r.PFP_at_formal_state_matched=mean(Q.A_est_mm>r.formal_state_matched_threshold_mm);
    r.required_threshold_ratio_to_fixed=q95/r.formal_fixed_threshold_mm;
    r.required_threshold_ratio_to_adaptive=q95/r.formal_adaptive_threshold_mm;
    rows(i)=r;
end
Summary=struct2table(rows);
Overall=table(mean(Summary.PFP_at_formal_fixed),mean(Summary.PFP_at_formal_adaptive),...
    mean(Summary.PFP_at_formal_state_matched),quantile(T.A_est_mm,0.95),...
    quantile(T.A_est_mm,0.99),...
    'VariableNames',{'mean_gap_PFP_fixed','mean_gap_PFP_adaptive','mean_gap_PFP_state_matched',...
    'pooled_q95_sm_delta_mm','pooled_q99_sm_delta_mm'});
outDir=fullfile(root,'output','17_state_threshold_transfer_audit');if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Summary,fullfile(outDir,'state_threshold_transfer_by_gap.csv'));
writetable(Overall,fullfile(outDir,'state_threshold_transfer_overall.csv'));
Result=struct('Summary',Summary,'Overall',Overall,'outputDir',outDir);
save(fullfile(outDir,'state_threshold_transfer_audit.mat'),'Result');
fprintf('STATE_THRESHOLD_TRANSFER_AUDIT_OK gaps=%d rows=%d\n',height(Summary),height(T));
end

function r=empty_row()
r=struct('g_truth_mm',NaN,'count',NaN,'q95_sm_delta_mm',NaN,'q99_sm_delta_mm',NaN,...
    'formal_fixed_threshold_mm',NaN,'formal_adaptive_threshold_mm',NaN,...
    'formal_state_matched_threshold_mm',NaN,'PFP_at_formal_fixed',NaN,...
    'PFP_at_formal_adaptive',NaN,'PFP_at_formal_state_matched',NaN,...
    'required_threshold_ratio_to_fixed',NaN,'required_threshold_ratio_to_adaptive',NaN);
end
