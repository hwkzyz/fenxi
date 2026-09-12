function out=R5_Run_V1Increment_Ablation_20260912(oldFile,sidecar,template,foundation,outputFile)
here=fileparts(mfilename('fullpath')); addpath(here); addpath(fullfile(here,'common'));
q=load(oldFile,'Result'); w=q.Result.WindowResult; cfg=struct(); cfg.dataset='ablation'; cfg.frequency=struct('candidateTopK',3,'maxOffsetHz',2); cfg.r5=struct('useV1DynamicIncrement',true,'auditAllEO',true,'useNestedFoundationAnchor',false);
ad=fullfile(tempdir,'r5_v1inc_adapter.mat'); R5_Adapt_OldGapOnlyWindows_20260911(oldFile,ad);
% The adapter is the foundation input because it carries CoreBundlePreview.
out=R5_Run_SensorConditionedWindowIdentification(cfg,@R5_Build_FromCorrectedLibrary,sidecar,template,ad,outputFile);
end
