function Result=Run_SingleSyncStructuredRegression(snrDb,numRevs)
%RUN_SINGLESYNCSTRUCTUREDREGRESSION Integer-EO voltage-domain regression.
if nargin<1,snrDb=5;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.snrDb=snrDb;
if nargin>=2&&~isempty(numRevs),cfg.NumRevs_high=numRevs;end
cfg.route30SingleFrequencyRangeHz=[300 1500];cfg.route30UseAllTrustedSamples=true;
cfg.syncSingleKeepPerGap=25;cfg.syncSingleRefineCount=12;
cfg.syncSingleIteratedReplayCount=2000;
eoList=[10 14 18];AList=[.10 .15 .25];gapPairs=[.2 .8;.8 .5;.5 .2];
rows=repmat(row0(),numel(eoList),1);fRot=cfg.RPM_high/60;
for i=1:numel(eoList)
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gapPairs(i,1),z),...
        cfg,lib.domain,cfg.snrDb,204000+i);
    lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);
    state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
    fTrue=eoList(i)*fRot;
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gapPairs(i,2),z),...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
        AList(i),fTrue,pi/4,'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(204100+i,'twister');
    d.V_cap=d.V_clean+sig/10^(cfg.snrDb/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    uTime=AList(i)*sin(2*pi*fTrue*map.t_v(:)+pi/4);
    B=[sin(eoList(i)*map.theta_v(:)),cos(eoList(i)*map.theta_v(:))];bc=B\uTime;
    uTheta=B*bc;VTheta=eval_gap_template(lib,gapPairs(i,2),map.x_v(:)-uTheta);
    thetaRmse=sqrt(mean((map.V_a(:)-VTheta).^2));thetaUErr=max(abs(uTime-uTheta));
    ticCase=tic;fit=run_single_sync_voltage_vp(map,lib,cfg,state);elapsed=toc(ticCase);
    rows(i)=struct('eo_true',eoList(i),'eo_est',fit.eo_id,'f_true',fTrue,...
        'f_est',fit.f_id,'A_true',AList(i),'A_est',fit.A_id,...
        'g_true',gapPairs(i,2),'g_est',fit.g_used,'rmse',fit.rmse,...
        'elapsed_s',elapsed,'used_samples',fit.used_sample_count,...
        'theta_model_rmse',thetaRmse,'theta_u_error',thetaUErr,...
        'success',fit.eo_id==eoList(i)&&abs(fit.g_used-gapPairs(i,2))<=.05);
end
Audit=struct2table(rows);out=fullfile(root,'output',sprintf('single_sync_structured_%gdb_%grev',snrDb,cfg.NumRevs_high));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'single_sync_structured.csv'));save(fullfile(out,'single_sync_structured.mat'),'Audit');disp(Audit);
if ~all(Audit.success),error('Structured synchronous single regression failed.');end
Result=struct('Audit',Audit,'outputDir',out);
end

function r=row0()
r=struct('eo_true',NaN,'eo_est',NaN,'f_true',NaN,'f_est',NaN,...
    'A_true',NaN,'A_est',NaN,'g_true',NaN,'g_est',NaN,'rmse',NaN,...
    'elapsed_s',NaN,'used_samples',NaN,'theta_model_rmse',NaN,...
    'theta_u_error',NaN,'success',false);
end
