function Detail=Run_DualSyncFunnelV1Check(snrDb,seeds,shortPairCount,runReference)
%RUN_DUALSYNCFUNNELV1CHECK Same-data reference versus Funnel V1 check.
if nargin<1||isempty(snrDb),snrDb=25;end
if nargin<2||isempty(seeds),seeds=1;end
if nargin<3||isempty(shortPairCount),shortPairCount=12;end
if nargin<4||isempty(runReference),runReference=true;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','dual_sync_funnel_v1');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.route30LowTemplateMethod="adaptive_sg";
cfg0.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',...
    [.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg0.route30FrequencyStructureMode="dual_sync_sync";cfg0.route30StructuredFastBudgets=false;
cfg0.route30DualFrequencyRangeHz=[500 1300];cfg0.route30GapHalfWidthMm=.20;
cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;
cfg0.dualSyncIteratedReplayCount=2100;cfg0.dualSyncRefineCount=50;
cfg0.dualSyncTopPairCount=40;cfg0.dualSyncStartsPerTopPair=2;
cfg0.funnelTopResidualPairCount=max(1,shortPairCount-2);cfg0.funnelShortGnPairCount=shortPairCount;
cfg0.funnelJointGnSteps=2;cfg0.funnelFullRefineCount=3;
cfg0.structuredComponentAmplitudeFloorMm=.05;cfg0.returnDualSyncDiagnostics=true;
fRot=cfg0.RPM_high/60;eoTrue=[10 26];ATrue=[.25 .15];
cfg0.f_true=eoTrue*fRot;cfg0.A_true=ATrue;cfg0.phi_true=[pi/4 -pi/3];
gLow=.5;deltaTrue=.1;truthT=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg0.alpha_k));
truthCal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',1:numel(cfg0.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);
sigma=reference_noise_std(truthModel,deltaTrue,cfg0,lib.domain,snrDb);
if runReference,strategies=["coarse_top_m_multistart","funnel_v1"];else,strategies="funnel_v1";end
rows=repmat(empty_row(),numel(seeds)*numel(strategies),1);ir=0;
checkpoint=fullfile(outDir,sprintf('checkpoint_mappednoise_v2_%gdB_short%d_ref%d.csv',...
    snrDb,shortPairCount,runReference));
if exist(checkpoint,'file')
    old=readtable(checkpoint,'TextType','string');
    keep=ismember(old.seed,seeds)&ismember(old.strategy,strategies);
    oldRows=table2struct(old(keep,:));ir=min(numel(oldRows),numel(rows));
    rows(1:ir)=oldRows(1:ir);fprintf('Resuming from %d completed fits.\n',ir);
end
allDone=true;
for seed=seeds(:).'
    for is=1:numel(strategies)
        allDone=allDone&&is_done(rows,ir,seed,strategies(is));
    end
end
if allDone
    Detail=struct2table(rows(1:ir));disp(Detail);return;
end
% One low-speed calibration per SNR operating condition.  It is reused by
% every high-speed noise seed so the Monte Carlo loop isolates high-speed
% identification variability and does not repeat the offline calibration.
lowSeed=991001;
low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg0,lib.domain,...
    sigma,lowSeed,'fixed_std');tc=tic;
cal=calibrate_inv_log_2_low_speed(low,lib,cfg0);calTime=toc(tc);
% The deterministic high-speed truth and spatial geometry are also common
% to all Monte Carlo seeds.  Generate/map them once; each seed changes only
% the additive fixed-standard-deviation voltage noise.
highClean=simulate_highspeed_from_low_increment(truthModel,deltaTrue,cfg0,Inf,...
    'fixed_std',1);
cleanMap=map_highspeed_to_space(highClean,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
for seed=seeds(:).'
    done=true;for is=1:numel(strategies),done=done&&is_done(rows,ir,seed,strategies(is));end
    if done,continue;end
    highSeed=992000+seed;
    rng(highSeed,'twister');highMap=cleanMap;
    highMap.V_a=cleanMap.V_a+sigma*randn(size(cleanMap.V_a));
    highMap.snr_db_equiv=snrDb;
    for is=1:numel(strategies)
        if is_done(rows,ir,seed,strategies(is)),continue;end
        cfg=cfg0;cfg.dualSyncSearchStrategy=strategies(is);tc=tic;
        fit=run_inv_log_2_high_with_calibration(cal,highMap,cfg);elapsed=toc(tc);
        ir=ir+1;rows(ir)=pack_row(seed,snrDb,sigma,strategies(is),fit,cal,calTime,...
            elapsed,eoTrue,ATrue,deltaTrue,lowSeed,highSeed);
        writetable(struct2table(rows(1:ir)),checkpoint);
        fprintf('%s seed=%d: EO=(%d,%d), dA=%.4g mm, ddg=%.4g mm, %.2f s\n',...
            strategies(is),seed,fit.eo_id(1),fit.eo_id(2),rows(ir).max_amplitude_error_mm,...
            rows(ir).delta_gap_error_mm,elapsed);
    end
end
Detail=struct2table(rows(1:ir));writetable(Detail,fullfile(outDir,'detail.csv'));
save(fullfile(outDir,'dual_sync_funnel_v1.mat'),'Detail','cfg0','eoTrue','ATrue',...
    'deltaTrue','snrDb','seeds','sigma','shortPairCount','runReference','-v7.3');disp(Detail);
end

function tf=is_done(rows,n,seed,strategy)
tf=false;
for i=1:n
    if rows(i).seed==seed&&string(rows(i).strategy)==strategy,tf=true;return;end
end
end

function r=empty_row()
r=struct('seed',NaN,'snr_db',NaN,'noise_std_V',NaN,'low_seed',NaN,'high_seed',NaN,...
    'strategy',"",'eo1_est',NaN,'eo2_est',NaN,'eo_correct',false,...
    'A1_est_mm',NaN,'A2_est_mm',NaN,'max_amplitude_error_mm',NaN,...
    'delta_gap_est_mm',NaN,'delta_gap_error_mm',NaN,'rmse_V',NaN,...
    'calibration_time_s',NaN,'fit_time_s',NaN,'coarse_time_s',NaN,...
    'level1_time_s',NaN,'level2_time_s',NaN,'full_refine_time_s',NaN,...
    'level1_candidate_count',NaN,'short_gn_pair_count',NaN,'full_refine_count',NaN,...
    'true_eo_rank_gn1',NaN,'true_eo_in_short_set',false,'true_eo_in_top3',false,...
    'full_iterations_total',NaN,'full_function_evaluations_total',NaN,'success',false);
end

function r=pack_row(seed,snr,sigma,strategy,fit,cal,calTime,elapsed,eoTrue,ATrue,dgTrue,lowSeed,highSeed)
r=empty_row();r.seed=seed;r.snr_db=snr;r.noise_std_V=sigma;r.low_seed=lowSeed;r.high_seed=highSeed;
r.strategy=strategy;r.eo1_est=fit.eo_id(1);r.eo2_est=fit.eo_id(2);r.eo_correct=all(fit.eo_id==eoTrue);
[~,ord]=sort(fit.eo_id);A=fit.A_id(ord);r.A1_est_mm=A(1);r.A2_est_mm=A(2);
r.max_amplitude_error_mm=max(abs(A-ATrue));r.delta_gap_est_mm=fit.g_used-cal.pathCal.g0;
r.delta_gap_error_mm=abs(r.delta_gap_est_mm-dgTrue);r.rmse_V=fit.rmse;
r.calibration_time_s=calTime;r.fit_time_s=elapsed;r.coarse_time_s=get_field(fit,'coarse_time_s',NaN);
r.level1_time_s=get_field(fit,'exact_gn1_time_s',get_field(fit,'replay_time_s',NaN));
r.level2_time_s=get_field(fit,'short_gn_time_s',NaN);r.full_refine_time_s=get_field(fit,'refine_time_s',NaN);
r.level1_candidate_count=get_field(fit,'level1_candidate_count',get_field(fit,'iterated_replay_count',NaN));
r.short_gn_pair_count=get_field(fit,'short_gn_pair_count',NaN);
r.full_refine_count=get_field(fit,'nonlinear_refine_count',NaN);
if strategy=="funnel_v1"
    eo=fit.diagnostic_pair_eo;score=fit.diagnostic_pair_gn1_score;[~,q]=sort(score);
    loc=find(all(eo(q,:)==eoTrue,2),1);if ~isempty(loc),r.true_eo_rank_gn1=loc;end
    r.true_eo_in_short_set=any(all(fit.diagnostic_short_eo==eoTrue,2));
    r.true_eo_in_top3=any(all(fit.diagnostic_full_start_eo==eoTrue,2));
    r.full_iterations_total=sum(fit.full_refine_iterations,'omitnan');
    r.full_function_evaluations_total=sum(fit.full_refine_function_evaluations,'omitnan');
end
r.success=r.eo_correct&&r.max_amplitude_error_mm<=.01&&r.delta_gap_error_mm<=.01;
end

function sigma=reference_noise_std(model,delta,cfg,domain,snr)
clean=simulate_highspeed_from_low_increment(model,delta,cfg,Inf,'snr_db',1);
q=cfg;q.A_true=[];q.f_true=[];q.phi_true=[];
stat=simulate_highspeed_from_low_increment(model,delta,q,Inf,'snr_db',1);
a=map_highspeed_to_space(clean,cfg.alpha_k,cfg.R_tip,domain,.02);
b=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,domain,.02);
sigma=rms(a.V_a-b.V_a)/10^(snr/20);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
