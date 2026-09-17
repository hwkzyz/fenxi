function Probe=Run_TrueEOMultistartProbe(snrDb,smoothSpan,seed)
%RUN_TRUEEOMULTISTARTPROBE Resolve which gap starts reach the true-EO basin.
if nargin<1,snrDb=5;end;if nargin<2,smoothSpan=9;end;if nargin<3,seed=1;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment";cfg.A_true=[.25 .15];fRot=cfg.RPM_high/60;eoTrue=[10 26];cfg.f_true=eoTrue*fRot;cfg.phi_true=[pi/4 -pi/3];
cfg.route30DualFrequencyRangeHz=eoTrue*fRot;cfg.dualSyncCandidateEOPairs=eoTrue;cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;
cfg.dualSyncSearchStrategy="balanced_per_pair";cfg.dualSyncBalancedStartsPerPair=11;cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
g0=.5;delta=.1;truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg.alpha_k));truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,snrDb,941000+seed);
high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,snrDb,'snr_db',942000+seed);map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'mean',struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',smoothSpan,'minCoverage',.90));
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,~]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);
state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;fit=run_dual_sync_voltage_vp(map,model,cfg,state);
Probe=table(fit.diagnostic_refine_start_g,fit.diagnostic_refine_start_dx,fit.diagnostic_refined_rmse,...
    fit.diagnostic_refined_eligible,'VariableNames',{'start_g_mm','start_dx_mm','refined_rmse_V','eligible'});
[~,ord]=sort(Probe.start_g_mm);Probe=Probe(ord,:);out=fullfile(root,'output','eo_search_strategy_study');if ~exist(out,'dir'),mkdir(out);end
tag=sprintf('%gdb_seed%d_SG%d',snrDb,seed,smoothSpan);writetable(Probe,fullfile(out,['true_eo_multistart_probe_',tag,'.csv']));
save(fullfile(out,['true_eo_multistart_probe_',tag,'.mat']),'Probe','fit','cfg','-v7.3');disp(Probe);
end
