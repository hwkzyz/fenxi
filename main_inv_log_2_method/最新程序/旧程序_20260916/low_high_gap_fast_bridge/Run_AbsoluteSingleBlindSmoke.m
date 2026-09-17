function R = Run_AbsoluteSingleBlindSmoke()
%RUN_ABSOLUTESINGLEBLINDSMOKE Reproduce the original absolute-template route.
% This is an A/B control for Run_LowIncrementSingleBlindSmoke.  The high-speed
% waveform and inverse model both use the absolute response surface R(g,x).

root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2); cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=20; cfg.NumRevs_high=8;
cfg.fs=base.cfgCase.fs; cfg.snrDb=Inf;
eoTrue=11; gLow=.50; gHigh=.70; ATrue=.20; phiTrue=.4;
fTrue=eoTrue*(cfg.RPM_high/60);

% Low-speed data is used only to provide the static reference gap.
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
    cfg,lib.domain,Inf,951001,'fixed_std');
% Original absolute forward model; no path template or low_increment calls.
high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
    cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
    ATrue,fTrue,phiTrue,'noise_ratio');
high.V_cap=high.V_clean;
map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);

% Route-30 single scan: absolute voltage VP candidates followed by exact
% bounded [g,dx,a,b,f] full-waveform refinement.
cfg.route30ForwardModel="absolute";
cfg.route30HierarchicalOrderSelection=true;
cfg.route30SingleGrid=5:5:1500;
cfg.route30SingleNuisanceHz=1500;
cfg.route30SingleCandidateCount=40;
cfg.route30SingleGapHalfWidthMm=.12;
cfg.route30DxHalfWidthMm=.50;
cfg.route30AmplitudeUpperMm=.80;
cfg.route30SingleRefineHalfWidthHz=20;
cfg.route30SingleMaxIter=180;
cfg.route30UseFastDualCandidates=false;
cfg.route30EnableRhoFallback=false;
cfg.route30FrequencyStructureMode="general";
ticCase=tic;
fit=run_inv_log_2_low_high_main(lowData,map,lib,cfg);
elapsed=toc(ticCase);

uTrue=ATrue*sin(2*pi*fTrue*map.t_v+phiTrue);
truthPred=eval_gap_template(lib,gHigh,map.x_v-uTrue);
ok=isfinite(truthPred)&isfinite(map.V_a);
truthRmse=sqrt(mean((map.V_a(ok)-truthPred(ok)).^2));
R=struct('g_low_true_mm',gLow,'g_high_true_mm',gHigh,...
    'g_hat_mm',fit.g_used,'eo_true',eoTrue,'eo_hat',fit.f_id/(cfg.RPM_high/60),...
    'frequency_true_hz',fTrue,'frequency_hat_hz',fit.f_id,...
    'amplitude_true_mm',ATrue,'amplitude_hat_mm',fit.A_id,...
    'phase_true_rad',phiTrue,'phase_hat_rad',fit.phi_id,...
    'rmse_V',fit.rmse,'truth_rmse_V',truthRmse,'elapsed_s',elapsed,...
    'method',"absolute_route30_single");
out=fullfile(root,'output','absolute_single_blind_smoke');
if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(R),fullfile(out,'absolute_single_blind_smoke.csv'));
save(fullfile(out,'absolute_single_blind_smoke.mat'),'R','fit');
disp(struct2table(R));
end
