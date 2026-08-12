function Result=Run_SynchronousOrderFastBenchmark()
%RUN_SYNCHRONOUSORDERFASTBENCHMARK EO-range search versus general dual route.
% The structured route searches every integer EO in the configured band.
% Truth is used only after fitting to calculate benchmark errors.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=8;
cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="absolute";cfg0.route30UseAllTrustedSamples=true;
fRot=cfg0.RPM_high/60;fTrue=[700 1200];eoTrue=fTrue/fRot;
Atrue=[.25 .15];phiTrue=[pi/4 -pi/3];gLow=.8;gHigh=.5;
snrList=[25 15 5];seed=20260809;methods=["general_fast","integer_eo_search"];
rows=repmat(row0(),0,1);
for isnr=1:numel(snrList)
    cfg=cfg0;cfg.snrDb=snrList(isnr);
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
        cfg,lib.domain,cfg.snrDb,seed+1000*isnr);
    lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);
    state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
    high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
        Atrue,fTrue,phiTrue,'noise_ratio');
    sig=rms(high.V_clean-min(high.V_clean));rng(seed+100000+1000*isnr,'twister');
    high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
    map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    for im=1:numel(methods)
        q=cfg;
        if methods(im)=="integer_eo_search"
            q.route30FrequencyStructureMode="dual_sync_sync";
            q.route30StructuredFastBudgets=true;
            q.dualSyncCandidateEOPairs=[]; % enumerate EO; never inject truth
        else
            q.route30FrequencyStructureMode="general";
            q.route30UseFastDualCandidates=true;
        end
        tc=tic;fit=run_inv_log_2_main_method(map,lib,q,state);elapsed=toc(tc);
        [fEst,ord]=sort(fit.f_id(:).');aEst=fit.A_id(ord);
        r=row0();r.method=methods(im);r.snr_db=cfg.snrDb;r.seed=seed;
        r.eo_grid_count=floor(1500/fRot)-ceil(300/fRot)+1;
        r.unordered_eo_pair_count=nchoosek(r.eo_grid_count,2);
        r.f1_est=fEst(1);r.f2_est=fEst(2);
        r.frequency_error_hz=max(abs(fEst-fTrue));
        r.amplitude_error_mm=max(abs(aEst-Atrue));r.gap_error_mm=fit.g_used-gHigh;
        r.rmse_V=fit.rmse;r.elapsed_s=elapsed;r.success=fit.model_order==2&&...
            r.frequency_error_hz<=2&&r.amplitude_error_mm<=.06&&abs(r.gap_error_mm)<=.05;
        r.candidate_count=getv(fit,'candidate_count',getv(fit,'retained_candidate_count',NaN));
        r.coarse_time_s=getv(fit,'coarse_time_s',NaN);
        r.replay_time_s=getv(fit,'replay_time_s',NaN);
        r.refine_time_s=getv(fit,'refine_time_s',NaN);
        rows(end+1,1)=r; %#ok<AGROW>
        fprintf('%s %g dB: f=%.3f,%.3f err=%.3g t=%.2f success=%d\n',...
            methods(im),cfg.snrDb,fEst(1),fEst(2),r.frequency_error_hz,elapsed,r.success);
    end
end
Audit=struct2table(rows);Summary=groupsummary(Audit,'method',{'mean','sum'},...
    {'success','elapsed_s','frequency_error_hz','amplitude_error_mm','gap_error_mm'});
out=fullfile(root,'output','synchronous_order_fast_benchmark');if~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'case_metrics.csv'));writetable(Summary,fullfile(out,'summary.csv'));
save(fullfile(out,'benchmark.mat'),'Audit','Summary','fTrue','eoTrue','-v7.3');
Result=struct('Audit',Audit,'Summary',Summary,'outputDir',out);disp(Summary);
end
function v=getv(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
function r=row0()
r=struct('method',"",'snr_db',NaN,'seed',NaN,'eo_grid_count',NaN,...
    'unordered_eo_pair_count',NaN,'f1_est',NaN,'f2_est',NaN,...
    'frequency_error_hz',NaN,'amplitude_error_mm',NaN,'gap_error_mm',NaN,...
    'rmse_V',NaN,'elapsed_s',NaN,'candidate_count',NaN,'coarse_time_s',NaN,...
    'replay_time_s',NaN,'refine_time_s',NaN,'success',false);
end
