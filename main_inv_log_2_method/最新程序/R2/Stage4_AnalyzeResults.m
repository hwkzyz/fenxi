function T = Stage4_AnalyzeResults(caseIDs)
%STAGE4_ANALYZERESULTS 汇总已完成的 R2 辨识结果，不重新生成波形。
if nargin<1, caseIDs=1:6; end
here=fileparts(mfilename('fullpath')); outDir=fullfile(here,'r2_profiled_recovery','output','stages'); rows=repmat(struct('case_id',NaN,'eo_true',NaN,'eo_hat',NaN,'gap_true_mm',NaN,'gap_hat_mm',NaN,'A_true_mm',NaN,'A_hat_mm',NaN,'rmse_V',NaN,'status',""),0,1);
for k=caseIDs
 f=fullfile(outDir,sprintf('R2_identification_case_%02d.mat',k)); if ~exist(f,'file'), warning('缺少结果，跳过 case %d',k); continue; end
 S=load(f,'estimate'); b=S.estimate.backend.best; t=S.estimate.truth;
 rows(end+1)=struct('case_id',k,'eo_true',t.eo,'eo_hat',b.eo,'gap_true_mm',t.gHigh,'gap_hat_mm',b.gap,'A_true_mm',t.A,'A_hat_mm',b.A,'rmse_V',b.fullWaveRmseV,'status',S.estimate.backend.status); %#ok<AGROW>
end
if isempty(rows), T=table(); return; end
T=struct2table(rows); writetable(T,fullfile(outDir,'R2_stage_summary.csv')); disp(T);
end
