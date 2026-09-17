function Result=Run_DualSyncStructuredRegression(snrDb,numRevs)
if nargin<1,snrDb=5;end
if nargin<2,numRevs=8;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=numRevs;cfg.snrDb=snrDb;
eoCases={[10 26],[14 24],[18 28]};gapPairs=[.2 .8;.8 .5;.5 .2];A=[.25 .15];phi=[pi/4 -pi/3];
rows=repmat(row0(),numel(eoCases),1);fRot=cfg.RPM_high/60;
for i=1:numel(eoCases)
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gapPairs(i,1),z),cfg,lib.domain,snrDb,206000+i);
    lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
    fTrue=eoCases{i}*fRot;
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gapPairs(i,2),z),...
        cfg.RPM_high,numRevs,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,A,fTrue,phi,'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(206100+i,'twister');d.V_cap=d.V_clean+sig/10^(snrDb/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticCase=tic;fit=run_dual_sync_voltage_vp(map,lib,cfg,state);elapsed=toc(ticCase);
    err=max(abs(fit.f_id-fTrue));rows(i)=struct('eo1_true',eoCases{i}(1),'eo2_true',eoCases{i}(2),...
        'eo1_est',fit.eo_id(1),'eo2_est',fit.eo_id(2),'max_frequency_error',err,...
        'g_true',gapPairs(i,2),'g_est',fit.g_used,'rmse',fit.rmse,'elapsed_s',elapsed,...
        'used_samples',fit.used_sample_count,'success',isequal(fit.eo_id,eoCases{i})&&abs(fit.g_used-gapPairs(i,2))<=.05);
end
Audit=struct2table(rows);out=fullfile(root,'output',sprintf('dual_sync_structured_%gdb_%grev',snrDb,numRevs));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'dual_sync_structured.csv'));save(fullfile(out,'dual_sync_structured.mat'),'Audit');disp(Audit);
if ~all(Audit.success),error('Structured dual-synchronous regression failed.');end
Result=struct('Audit',Audit,'outputDir',out);
end
function r=row0()
r=struct('eo1_true',NaN,'eo2_true',NaN,'eo1_est',NaN,'eo2_est',NaN,...
    'max_frequency_error',NaN,'g_true',NaN,'g_est',NaN,'rmse',NaN,...
    'elapsed_s',NaN,'used_samples',NaN,'success',false);
end
