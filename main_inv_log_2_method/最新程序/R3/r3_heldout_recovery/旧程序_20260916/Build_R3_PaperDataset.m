function Result = Build_R3_PaperDataset()
%BUILD_R3_PAPERDATASET Freeze the paper-facing R3 data convention.
root=fileparts(mfilename('fullpath')); addpath(fullfile(root,'local_func')); P=R3_Protocol();
formalDir=fullfile(root,'output','05_main_15db'); correctedDir=fullfile(root,'output','16_corrected_state_matched_benchmark');
Formal=readtable(fullfile(formalDir,'main_detail_merged.csv'),'TextType','string');
Corrected=readtable(fullfile(correctedDir,'corrected_state_matched_detail.csv'),'TextType','string');
assert(height(Formal)==9450&&height(Corrected)==3150&&all(Corrected.coordinate_mode=="SM_delta"));
[Summary,~]=Summarize_R3(Formal,P); Summary=Summary(Summary.method=="fixed"|Summary.method=="adaptive",:);
Summary=[Summary;corrected_summary(Corrected,Formal,P)]; Summary=sortrows(Summary,{'method','g_truth_mm','A_true_mm'}); assert(height(Summary)==126);
Paired=readtable(fullfile(formalDir,'main_paired_summary.csv'),'TextType','string'); assert(height(Paired)==36);
[FailureOverall,FailureByGap,FailureDetail]=fixed_failure_20deg(Formal,P);
outDir=fullfile(root,'output','18_paper_data_sync'); if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Summary,fullfile(outDir,'r3_paper_summary.csv')); writetable(Paired,fullfile(outDir,'r3_paper_paired_summary.csv'));
writetable(FailureOverall,fullfile(outDir,'fixed_failure_composition_20deg_overall.csv'));
writetable(FailureByGap,fullfile(outDir,'fixed_failure_composition_20deg_by_gap.csv'));
writetable(FailureDetail,fullfile(outDir,'fixed_failure_rows_20deg.csv'));
Result=struct('P',P,'Summary',Summary,'Paired',Paired,'CorrectedStateMatchedDetail',Corrected,...
 'FailureOverall',FailureOverall,'FailureByGap',FailureByGap,'outputDir',outDir,...
 'state_matched_status',"corrected reference-anchored coordinate only",...
 'replicate_design',"75 uniformly spaced phase values, each paired one-to-one with a deterministic high-speed noise seed");
save(fullfile(outDir,'r3_paper_result.mat'),'Result','-v7.3');
fprintf('R3_PAPER_DATA_SYNC_OK summary=%d paired=%d corrected_sm=%d\n',height(Summary),height(Paired),height(Corrected));
end

function S=corrected_summary(C,Formal,P)
gaps=sort(unique(C.g_truth_mm)).'; amps=sort(unique(C.A_true_mm)).'; rows=repmat(empty_summary_row(),numel(gaps)*numel(amps),1); k=0;
for g=gaps
 for A=amps
  Q=C(C.g_truth_mm==g&C.A_true_mm==A,:); assert(height(Q)==75);
  truth=Formal(Formal.method=="fixed"&Formal.g_truth_mm==g&Formal.A_true_mm==A,{'high_seed','A_truth_mapped_mm'});
  M=innerjoin(Q,truth,'Keys','high_seed'); assert(height(M)==75); k=k+1; r=empty_summary_row();
  r.method="state_matched_corrected"; r.g_truth_mm=g; r.delta_gap_mm=g-P.gReferenceMm; r.A_true_mm=A; r.n=height(M);
  r.median_amplitude_error_mm=median(M.A_est_mm-M.A_truth_mapped_mm,'omitnan'); r.median_gap_error_mm=0;
  if A>0
   n=height(M); nv=nnz(M.vib_success); r.P_vib=nv/n; [r.P_vib_wilson_low,r.P_vib_wilson_high]=r3_wilson_interval(nv,n,P.wilsonAlpha);
   r.P_joint=r.P_vib; r.P_joint_wilson_low=r.P_vib_wilson_low; r.P_joint_wilson_high=r.P_vib_wilson_high;
  else
   n=height(M); nf=nnz(M.corrected_false_positive); r.P_FP=nf/n; [r.P_FP_wilson_low,r.P_FP_wilson_high]=r3_wilson_interval(nf,n,P.wilsonAlpha);
  end
  rows(k)=r;
 end
end
S=struct2table(rows);
end

function [Overall,ByGap,D]=fixed_failure_20deg(Formal,P)
D=Formal(Formal.method=="fixed"&Formal.A_true_mm>0,:); D.eo_fail=D.eo_est~=D.eo_true;
D.amplitude_fail=abs(D.amplitude_error_mm)>D.amplitude_tolerance_mm; D.phase_fail=D.phase_error_rad>P.phaseToleranceRad; D.failure_pattern=strings(height(D),1);
for i=1:height(D)
 if ~D.eo_fail(i)&&~D.amplitude_fail(i)&&~D.phase_fail(i),D.failure_pattern(i)="success";continue;end
 labels=strings(0,1); if D.eo_fail(i),labels(end+1)="EO";end; if D.amplitude_fail(i),labels(end+1)="amplitude";end; if D.phase_fail(i),labels(end+1)="phase";end
 D.failure_pattern(i)=strjoin(labels," + ");
end
Overall=composition_table(D,NaN); gaps=sort(unique(D.g_truth_mm)); ByGap=table();
for i=1:numel(gaps),ByGap=[ByGap;composition_table(D(D.g_truth_mm==gaps(i),:),gaps(i))];end %#ok<AGROW>
assert(sum(Overall.count)==height(D));
end

function O=composition_table(D,g)
patterns=sort(unique(D.failure_pattern)); O=table(); failures=nnz(D.failure_pattern~="success");
for p=patterns.'
 count=nnz(D.failure_pattern==p); fracFail=NaN; if p~="success",fracFail=count/failures;end
 O=[O;table(g,p,count,count/height(D),fracFail,'VariableNames',{'g_truth_mm','failure_pattern','count','fraction_of_all','fraction_of_failures'})]; %#ok<AGROW>
end
end

function r=empty_summary_row()
r=struct('method',"",'g_truth_mm',NaN,'delta_gap_mm',NaN,'A_true_mm',NaN,'n',NaN,'P_vib',NaN,'P_vib_wilson_low',NaN,'P_vib_wilson_high',NaN,'P_FP',NaN,'P_FP_wilson_low',NaN,'P_FP_wilson_high',NaN,'P_joint',NaN,'P_joint_wilson_low',NaN,'P_joint_wilson_high',NaN,'median_amplitude_error_mm',NaN,'median_gap_error_mm',NaN);
end
