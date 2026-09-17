function Result=Run_SingleAsyncStructuredRegression(snrDb,numRevs)
%RUN_SINGLEASYNCSTRUCTUREDREGRESSION Known asynchronous single-frequency audit.
if nargin<1,snrDb=5;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.snrDb=snrDb;
if nargin>=2&&~isempty(numRevs),cfg.NumRevs_high=numRevs;end
cfg.route30SingleFrequencyRangeHz=[300 1500];cfg.route30FrequencyStepHz=5;
fList=[531 733 917];AList=[.10 .15 .25];gapPairs=[.2 .8;.8 .5;.5 .2];
rows=repmat(row0(),numel(fList),1);
for i=1:numel(fList)
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gapPairs(i,1),z),...
        cfg,lib.domain,snrDb,205000+i);
    lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);
    state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gapPairs(i,2),z),...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
        AList(i),fList(i),pi/4,'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(205100+i,'twister');
    d.V_cap=d.V_clean+sig/10^(snrDb/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticCase=tic;fit=run_single_async_voltage_vp(map,lib,cfg,state);elapsed=toc(ticCase);
    rows(i)=struct('f_true',fList(i),'f_est',fit.f_id,'frequency_error',abs(fit.f_id-fList(i)),...
        'A_true',AList(i),'A_est',fit.A_id,'g_true',gapPairs(i,2),'g_est',fit.g_used,...
        'rmse',fit.rmse,'elapsed_s',elapsed,'used_samples',fit.used_sample_count,...
        'success',abs(fit.f_id-fList(i))<=1&&abs(fit.g_used-gapPairs(i,2))<=.05);
end
Audit=struct2table(rows);out=fullfile(root,'output',sprintf('single_async_structured_%gdb_%grev',snrDb,cfg.NumRevs_high));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'single_async_structured.csv'));save(fullfile(out,'single_async_structured.mat'),'Audit');disp(Audit);
if ~all(Audit.success),error('Structured asynchronous single regression failed.');end
Result=struct('Audit',Audit,'outputDir',out);
end

function r=row0()
r=struct('f_true',NaN,'f_est',NaN,'frequency_error',NaN,'A_true',NaN,'A_est',NaN,...
    'g_true',NaN,'g_est',NaN,'rmse',NaN,'elapsed_s',NaN,'used_samples',NaN,'success',false);
end
