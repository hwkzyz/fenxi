function Result = Run_11_StateMatchedCoordinateAudit()
%RUN_11_STATEMATCHEDCOORDINATEAUDIT Merge and audit all coordinate replays.

root=fileparts(mfilename('fullpath'));
gaps=[0.5 0.7 0.8 0.9 1.0 1.2];
All=table();detailRows=0;
for i=1:numel(gaps)
    g=gaps(i);
    if g==0.8
        outDir=fullfile(root,'output','08_state_matched_coordinate_diagnostic');
        S=readtable(fullfile(outDir,'reference_coordinate_summary.csv'),'TextType','string');
        D=readtable(fullfile(outDir,'reference_coordinate_detail.csv'),'TextType','string');
    else
        outDir=fullfile(root,'output','09_state_matched_offreference_diagnostic');
        tag=sprintf('g%04d',round(1000*g));
        S=readtable(fullfile(outDir,[tag '_summary.csv']),'TextType','string');
        D=readtable(fullfile(outDir,[tag '_detail.csv']),'TextType','string');
    end
    assert(height(S)==14 && all(S.GroupCount==75),'Unexpected summary cardinality at g=%.3f.',g);
    assert(height(D)==1050,'Unexpected detail cardinality at g=%.3f.',g);
    S=S(:,{'coordinate_mode','A_true_mm','GroupCount','P_vib','P_FP','P_reproduces_formal_sm'});
    S.g_truth_mm=repmat(g,height(S),1);
    S=movevars(S,'g_truth_mm','Before','coordinate_mode');
    All=[All;S]; %#ok<AGROW>
    detailRows=detailRows+height(D);
end
assert(height(All)==84 && detailRows==6300,'Coordinate audit is incomplete.');
assert(all(All.P_reproduces_formal_sm(All.coordinate_mode=="SM_abs")==1),...
    'At least one old-coordinate replay does not reproduce the formal state-matched row.');

modes=["SM_abs" "SM_delta"];
GapSummary=table();Overall=table();
for i=1:numel(gaps)
    for j=1:numel(modes)
        Q=All(All.g_truth_mm==gaps(i)&All.coordinate_mode==modes(j),:);
        nonzero=Q.A_true_mm>0;
        GapSummary=[GapSummary;table(gaps(i),modes(j),sum(Q.GroupCount(nonzero)),...
            sum(Q.P_vib(nonzero).*Q.GroupCount(nonzero))/sum(Q.GroupCount(nonzero)),...
            Q.P_FP(Q.A_true_mm==0),...
            'VariableNames',{'g_truth_mm','coordinate_mode','nonzero_count','P_vib_pooled','P_FP'})]; %#ok<AGROW>
    end
end
for j=1:numel(modes)
    Q=GapSummary(GapSummary.coordinate_mode==modes(j),:);
    Overall=[Overall;table(modes(j),sum(Q.nonzero_count),...
        sum(Q.P_vib_pooled.*Q.nonzero_count)/sum(Q.nonzero_count),...
        mean(Q.P_FP),'VariableNames',{'coordinate_mode','nonzero_count','P_vib_pooled','P_FP_mean_across_gaps'})]; %#ok<AGROW>
end
outDir=fullfile(root,'output','11_state_matched_coordinate_audit');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(All,fullfile(outDir,'all_cell_summary.csv'));
writetable(GapSummary,fullfile(outDir,'gap_summary.csv'));
writetable(Overall,fullfile(outDir,'overall_summary.csv'));
Result=struct('All',All,'GapSummary',GapSummary,'Overall',Overall,'detailRows',detailRows,'outputDir',outDir);
save(fullfile(outDir,'state_matched_coordinate_audit.mat'),'Result');
fprintf('STATE_MATCHED_COORDINATE_AUDIT_OK detail_rows=%d summary_rows=%d\n',detailRows,height(All));
end
