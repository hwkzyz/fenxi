function R = Run_LowIncrementConsistencyAudit()
%RUN_LOWINCREMENTCONSISTENCYAUDIT Minimal audit of the frozen path contract.
root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
g0=median(lib.gapTrain); x=lib.xGrid(:);
path=make_low_increment_path_model(lib,g0);
delta=.10;
Tlow=path.templateLow(:);
y0=eval_path_increment_template(path,0,x);
y=eval_path_increment_template(path,delta,x);
% Independent reference evaluation of the documented equation.
Rbase=eval_gap_template(lib,g0,x);
Rquery=eval_gap_template(lib,g0+delta,x);
yRef=Tlow+(Rquery-Rbase);
R=struct('g0_mm',g0,'delta_g_mm',delta,...
    'zero_increment_max_error_V',max(abs(y0-Tlow)),...
    'formula_max_error_V',max(abs(y-yRef)),...
    'public_absolute_query_max_error_V',max(abs(eval_gap_template(path,g0+delta,x)-y)),...
    'pass',max(abs(y0-Tlow))<1e-10 && max(abs(y-yRef))<1e-10 && ...
    max(abs(eval_gap_template(path,g0+delta,x)-y))<1e-4);
out=fullfile(root,'output','low_increment_consistency_audit');
if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(R),fullfile(out,'low_increment_consistency_audit.csv'));
save(fullfile(out,'low_increment_consistency_audit.mat'),'R');
disp(struct2table(R));
end
