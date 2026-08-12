function Result=Run_SyncAsyncStructuredRegression(snrDb,numRevs)
if nargin<1,snrDb=5;end
if nargin<2,numRevs=8;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=numRevs;cfg.snrDb=snrDb;
eo=[10 14 18];fa=[733 1217 1383];gapPairs=[.2 .8;.8 .5;.5 .2];A=[.25 .15];phi=[pi/4 -pi/3];
rows=repmat(row0(),numel(eo),1);fRot=cfg.RPM_high/60;
for i=1:numel(eo)
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gapPairs(i,1),z),cfg,lib.domain,snrDb,207000+i);
    lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
    fTrue=[eo(i)*fRot,fa(i)];
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gapPairs(i,2),z),...
        cfg.RPM_high,numRevs,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,A,fTrue,phi,'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(207100+i,'twister');d.V_cap=d.V_clean+sig/10^(snrDb/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticCase=tic;fit=run_sync_async_voltage_vp(map,lib,cfg,state);elapsed=toc(ticCase);
    err=max(abs(fit.f_id-fTrue));rows(i)=struct('eo_true',eo(i),'eo_est',fit.sync_eo,...
        'async_true',fa(i),'async_est',fit.async_frequency_hz,'max_frequency_error',err,...
        'g_true',gapPairs(i,2),'g_est',fit.g_used,'rmse',fit.rmse,'elapsed_s',elapsed,...
        'used_samples',fit.used_sample_count,'success',fit.sync_eo==eo(i)&&...
        abs(fit.async_frequency_hz-fa(i))<=1&&abs(fit.g_used-gapPairs(i,2))<=.05);
end
Audit=struct2table(rows);out=fullfile(root,'output',sprintf('sync_async_structured_%gdb_%grev',snrDb,numRevs));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'sync_async_structured.csv'));save(fullfile(out,'sync_async_structured.mat'),'Audit');disp(Audit);
if ~all(Audit.success),error('Structured synchronous-asynchronous regression failed.');end
Result=struct('Audit',Audit,'outputDir',out);
end
function r=row0()
r=struct('eo_true',NaN,'eo_est',NaN,'async_true',NaN,'async_est',NaN,...
    'max_frequency_error',NaN,'g_true',NaN,'g_est',NaN,'rmse',NaN,...
    'elapsed_s',NaN,'used_samples',NaN,'success',false);
end
