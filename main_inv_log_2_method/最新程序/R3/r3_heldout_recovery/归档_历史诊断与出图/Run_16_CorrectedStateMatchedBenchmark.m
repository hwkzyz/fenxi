function Result = Run_16_CorrectedStateMatchedBenchmark()
%RUN_16_CORRECTEDSTATEMATCHEDBENCHMARK Assemble corrected oracle benchmark.

root=fileparts(mfilename('fullpath'));gaps=[0.5 0.7 0.8 0.9 1.0 1.2];Detail=table();
Formal=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
fixedThreshold=unique(Formal.detection_threshold_mm(Formal.method=="fixed"));
assert(numel(fixedThreshold)==1&&isfinite(fixedThreshold));
for g=gaps
    if g==0.8
        D=readtable(fullfile(root,'output','08_state_matched_coordinate_diagnostic',...
            'reference_coordinate_detail.csv'),'TextType','string');
    else
        tag=sprintf('g%04d',round(1000*g));
        D=readtable(fullfile(root,'output','09_state_matched_offreference_diagnostic',...
            [tag '_detail.csv']),'TextType','string');
    end
    D=D(D.coordinate_mode=="SM_delta",:);assert(height(D)==525);
    D.corrected_detection_threshold_mm=repmat(fixedThreshold,height(D),1);
    D.corrected_false_positive=D.A_true_mm==0&D.A_est_mm>fixedThreshold;
    D=D(:,{'coordinate_mode','g_truth_mm','A_true_mm','replicate_index','high_seed','g_query_mm',...
        'eo_est','A_est_mm','phase_error_rad','rmse_V','vib_success',...
        'corrected_detection_threshold_mm','corrected_false_positive'});
    Detail=[Detail;D]; %#ok<AGROW>
end
assert(height(Detail)==3150);
Cell=table();Gap=table();
for g=gaps
    G=Detail(Detail.g_truth_mm==g,:);assert(height(G)==525);
    amps=unique(G.A_true_mm);
    for A=amps.'
        Q=G(G.A_true_mm==A,:);assert(height(Q)==75);
        Cell=[Cell;table(g,A,height(Q),mean(Q.vib_success),mean(Q.corrected_false_positive),...
            'VariableNames',{'g_truth_mm','A_true_mm','count','P_vib','P_FP'})]; %#ok<AGROW>
    end
    N=G(G.A_true_mm>0,:);Z=G(G.A_true_mm==0,:);
    Gap=[Gap;table(g,height(N),mean(N.vib_success),mean(Z.corrected_false_positive),...
        'VariableNames',{'g_truth_mm','nonzero_count','P_vib_pooled','P_FP'})]; %#ok<AGROW>
end
assert(height(Cell)==42&&all(Cell.count==75));
Overall=table(height(Detail(Detail.A_true_mm>0,:)),mean(Detail.vib_success(Detail.A_true_mm>0)),...
    mean(Detail.corrected_false_positive(Detail.A_true_mm==0)),fixedThreshold,...
    'VariableNames',{'nonzero_count','P_vib_pooled','P_FP_pooled','corrected_threshold_mm'});
outDir=fullfile(root,'output','16_corrected_state_matched_benchmark');if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'corrected_state_matched_detail.csv'));
writetable(Cell,fullfile(outDir,'corrected_state_matched_cell_summary.csv'));
writetable(Gap,fullfile(outDir,'corrected_state_matched_gap_summary.csv'));
writetable(Overall,fullfile(outDir,'corrected_state_matched_overall_summary.csv'));
Result=struct('Detail',Detail,'Cell',Cell,'Gap',Gap,'Overall',Overall,'outputDir',outDir,...
    'threshold_provenance',"Corrected state match equals fixed at reference g=0.8, so it shares the fixed null threshold.");
save(fullfile(outDir,'corrected_state_matched_benchmark.mat'),'Result','-v7.3');
fprintf('CORRECTED_STATE_MATCHED_BENCHMARK_OK rows=%d cells=%d\n',height(Detail),height(Cell));
end
