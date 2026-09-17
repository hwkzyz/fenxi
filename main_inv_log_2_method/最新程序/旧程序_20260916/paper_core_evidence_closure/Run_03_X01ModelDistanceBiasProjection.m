function Result=Run_03_X01ModelDistanceBiasProjection()
% Connect constrained model distance, directed discrepancy and full margin.
root=fileparts(mfilename('fullpath'));programDir=fileparts(root);
addpath(programDir,'-begin');addpath(fullfile(programDir,'local_func'),'-begin');
addpath(fullfile(programDir,'identifiability_mechanism_analysis','local_func'),'-begin');
src=fullfile(programDir,'low_high_gap_fast_bridge','output','x01_deterministic_bias_audit','audit.mat');
if ~isfile(src),error('route32:MissingX01','Run X01 deterministic bias audit first.');end
S=load(src,'map','cal','c','cfg');map=S.map;c=S.c;cfg=S.cfg;formal=S.cal.templateModel;
ctx=load_inv_log_2_project_context(programDir);tr=load_fixed_trust_domain(programDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
reference=identity_model(lib,c.gLow,cfg);trueEO=c.eoTrue;wrongEO=[10 18];
refT=fixed_fit(map,reference,cfg,trueEO,c.gHigh);refW=fixed_fit(map,reference,cfg,wrongEO,c.gHigh);
forT=fixed_fit(map,formal,cfg,trueEO,c.gHigh);forW=fixed_fit(map,formal,cfg,wrongEO,c.gHigh);
valid=all(isfinite([refT.VFit refW.VFit forT.VFit forW.VFit map.V_a(:)]),2);
y=map.V_a(valid);d=refW.VFit(valid)-refT.VFit(valid);Dnorm=norm(d);h=d/max(Dnorm,eps);
bTrue=forT.VFit(valid)-refT.VFit(valid);bWrong=forW.VFit(valid)-refW.VFit(valid);
signedProj=dot(bTrue,h)/sqrt(nnz(valid));beta=abs(signedProj);D=Dnorm/sqrt(nnz(valid));
Mref=refW.rmse-refT.rmse;Mformal=forW.rmse-forT.rmse;
% Fixed-reference first-order diagnostic: ||d+b||^2-||b||^2.
localShiftRatio=2*dot(bTrue,d)/max(dot(d,d),eps);
Summary=table(D,beta,signedProj,beta/max(D,eps),localShiftRatio,...
    norm(bTrue)/sqrt(nnz(valid)),norm(bWrong)/sqrt(nnz(valid)),...
    refT.rmse,refW.rmse,Mref,forT.rmse,forW.rmse,Mformal,Mref>0&&Mformal<0,...
    'VariableNames',{'reference_model_distance_rms_V','bias_projection_abs_rms_V',...
    'bias_projection_signed_rms_V','bias_to_distance_ratio','local_squared_margin_shift_ratio',...
    'true_manifold_discrepancy_rms_V','wrong_manifold_discrepancy_rms_V',...
    'reference_true_rmse_V','reference_wrong_rmse_V','reference_margin_V',...
    'formal_true_rmse_V','formal_wrong_rmse_V','formal_margin_V','full_margin_flip'});
Wave=table(y,refT.VFit(valid),refW.VFit(valid),forT.VFit(valid),forW.VFit(valid),d,bTrue,bWrong,...
    'VariableNames',{'observation_V','reference_true_fit_V','reference_wrong_fit_V',...
    'formal_true_fit_V','formal_wrong_fit_V','reference_discriminant_V',...
    'true_model_discrepancy_V','wrong_model_discrepancy_V'});
out=fullfile(root,'output','03_x01_model_distance_bias_projection');if ~exist(out,'dir'),mkdir(out);end
writetable(Summary,fullfile(out,'summary.csv'));writetable(Wave,fullfile(out,'waveform_components.csv'));
save(fullfile(out,'x01_distance_bias.mat'),'Summary','Wave','refT','refW','forT','forW','-v7.3');
write_report(out,Summary);disp(Summary);Result=struct('Summary',Summary,'outputDir',out);
end
function model=identity_model(lib,g,cfg)
T=repmat(eval_gap_template(lib,g,lib.xGrid),1,numel(cfg.alpha_k));
pc=struct('g0',g,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));
model=build_path_template_model(lib,T,pc);model=rmfield(model,'pathCache');
end
function fit=fixed_fit(map,model,cfg,pair,g)
q=cfg;q.dualSyncCandidateEOPairs=pair;q.funnelFullRefineCount=1;q.dualSyncMaxIter=600;
state=struct('gHat',model.pathCal.g0,'g_low_hat',model.pathCal.g0,'dx0',0);
fit=run_dual_sync_voltage_vp_funnel(map,model,q,state).fit;
end
function write_report(out,S)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# X01 model-distance and directed-bias closure\n\n');
fprintf(fid,'- Reference constrained model distance: %.9g V RMS\n',S.reference_model_distance_rms_V);
fprintf(fid,'- Absolute true-model discrepancy projection: %.9g V RMS\n',S.bias_projection_abs_rms_V);
fprintf(fid,'- Bias/distance diagnostic ratio: %.9g\n',S.bias_to_distance_ratio);
fprintf(fid,'- Reference wrong-minus-true margin: %.9g V\n',S.reference_margin_V);
fprintf(fid,'- Formal wrong-minus-true margin: %.9g V\n',S.formal_margin_V);
fprintf(fid,'- Full optimized margin flip: %d\n\n',S.full_margin_flip);
fprintf(fid,'The projection ratio is a local diagnostic. The optimized residual-margin sign is the authoritative model-selection evidence.\n');fclose(fid);
end
