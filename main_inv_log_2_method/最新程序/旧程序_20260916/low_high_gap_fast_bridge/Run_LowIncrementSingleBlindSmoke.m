function R = Run_LowIncrementSingleBlindSmoke()
%RUN_LOWINCREMENTSINGLEBLINDSMOKE Single-frequency blind backend smoke test.
root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2); cfg=base.cfgCase;
% Keep the smoke test representative of the formal simulation contract:
% enough low-speed revolutions to build the frozen template and eight
% high-speed revolutions for waveform identification.
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=20; cfg.NumRevs_high=8;
cfg.fs=base.cfgCase.fs; eoTrue=11; cfg.A_true=.20; cfg.f_true=eoTrue*(cfg.RPM_high/60); cfg.phi_true=.4;
gLow=.50; gHigh=.70;
% Generate the low-speed reference and freeze the path.
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg,lib.domain,Inf,951001,'fixed_std');
low=build_low_speed_templates_adaptive(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'cv_sg',struct());
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false);
[pathCal,~]=calibrate_low_speed_path(lib,low,opts);
path=build_path_template_model(lib,low.templateLow,pathCal);
% Diagnostic switch: the exact pointwise evaluator is retained for a
% controlled comparison against the cached production evaluator.
useExactEvaluator = false;
if useExactEvaluator, path.pathCache = []; end
high=simulate_highspeed_from_low_increment(path,gHigh-pathCal.g0,cfg,Inf,'fixed_std',951002);
map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
% diagnostic: retain all sensors in the formal smoke test
bundle=struct('t',map.t_v,'x',map.x_v,'V',map.V_a,'theta',map.theta_v,...
    'sensorId',map.S_v,'turnId',map.rev_v,'Wvp',ones(size(map.V_a)),...
    'Wfull',ones(size(map.V_a)));
response=struct('evalF',@(dg,xq,sid)eval_path_increment_template(path,dg,xq,sid),...
    'evalFx',@(dg,xq,sid)eval_gap_derivative(path,dg,xq,sid));
uTruth=cfg.A_true*sin(eoTrue*map.theta_v+cfg.phi_true);
truthPred=eval_path_increment_template(path,gHigh-pathCal.g0,map.x_v-uTruth,map.S_v);
truthValid=isfinite(truthPred)&isfinite(map.V_a);
truthRmse=sqrt(mean((map.V_a(truthValid)-truthPred(truthValid)).^2));
% Explicit support audit: samples outside the displaced spatial support are
% baseline samples in the simulator and must not be interpreted as response
% surface observations by the inverse model.
supportTruth=(map.x_v-uTruth)>=lib.domain(1)&(map.x_v-uTruth)<=lib.domain(2);
truthSupportFraction=mean(supportTruth);
rotHz=cfg.RPM_high/60; dgGrid=(min(lib.gapTrain)-pathCal.g0):.02:(max(lib.gapTrain)-pathCal.g0);
bc=struct('gapBounds',[min(dgGrid),max(dgGrid)],'gapGrid',dgGrid,...
    'eoGrid',1:15,'topK',3,'rotHz',rotHz,'minGradient',1e-8,'minSamples',12,...
    'dxBound',.20,'ampCoeffBound',.8,'deltaFBound',2,'maxIter',60,'ambiguityMargin',.001,...
    'refineAmplitudeGrid',0:.025:.30);
fit=run_unified_blind_dynamic_backend(bundle,response,bc);
% Diagnostic only: evaluate the nearest truth-independent refinement-grid
% point for the true simulated case. This value is not used for selection.
R=struct('g_low_true_mm',gLow,'g_low_cal_mm',pathCal.g0,...
    'g_high_true_mm',gHigh,'delta_g_true_mm',gHigh-gLow,...
    'delta_g_hat_mm',fit.best.gap,'g_high_hat_mm',pathCal.g0+fit.best.gap,...
    'eo_true',eoTrue,'eo_hat',fit.best.eo,'frequency_true_hz',cfg.f_true,...
    'frequency_hat_hz',(fit.best.eo+fit.best.deltaF/rotHz)*rotHz,...
    'amplitude_true_mm',cfg.A_true,'amplitude_hat_mm',fit.best.A,...
    'rmse_V',fit.best.fullWaveRmseV,'truth_rmse_V',truthRmse,...
    'truth_support_fraction',truthSupportFraction,'status',fit.status,...
    'pass',isfinite(fit.best.fullWaveRmseV)&&abs(fit.best.gap-(gHigh-gLow))<.03&&fit.best.eo==11);
out=fullfile(root,'output','low_increment_single_blind_smoke');
if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(R),fullfile(out,'low_increment_single_blind_smoke.csv'));
save(fullfile(out,'low_increment_single_blind_smoke.mat'),'R','fit'); disp(struct2table(R));
fig=figure('Visible','off','Color','w','Position',[100 100 1000 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; scatter(map.x_v, uTruth, 8, double(supportTruth), 'filled'); box on; xlabel('x (mm)'); ylabel('u truth (mm)'); title('Displaced support mask'); colorbar;
nexttile; plot(map.x_v(truthValid),map.V_a(truthValid),'k.','MarkerSize',4); hold on; plot(map.x_v(truthValid),truthPred(truthValid),'r.','MarkerSize',4); box on; xlabel('x (mm)'); ylabel('V'); legend('observation','truth replay'); title(sprintf('Truth replay RMSE %.3g V',truthRmse));
nexttile; histogram(map.x_v-uTruth,40); hold on; xline(lib.domain(1),'r--'); xline(lib.domain(2),'r--'); box on; xlabel('x-u (mm)'); title('Displaced coordinate support');
nexttile; plot(map.t_v,map.V_a-truthPred,'b.','MarkerSize',4); yline(0,'k-'); box on; xlabel('t (s)'); ylabel('V residual'); title(sprintf('valid %.1f%%',100*truthSupportFraction));
exportgraphics(fig,fullfile(out,'low_increment_truth_support_diagnostic.png'),'Resolution',180); close(fig);
% Objective profile at the true EO/gap (diagnostic only; truth is not passed
% to the estimator).  This distinguishes a genuine amplitude-identifiability
% problem from a local optimizer problem.
Agrid=linspace(0,0.30,121); sseA=nan(size(Agrid));
for ia=1:numel(Agrid)
    uu=Agrid(ia)*sin(eoTrue*map.theta_v+cfg.phi_true);
    pp=eval_path_increment_template(path,gHigh-pathCal.g0,map.x_v-uu,map.S_v);
    ok=isfinite(pp)&isfinite(map.V_a); sseA(ia)=mean((map.V_a(ok)-pp(ok)).^2);
end
fig=figure('Visible','off','Color','w','Position',[100 100 650 450]);
plot(Agrid,sqrt(sseA),'LineWidth',1.4); hold on; xline(cfg.A_true,'k--','true A'); xline(fit.best.A,'r--','estimated A');
box on; xlabel('Amplitude (mm)'); ylabel('RMSE (V)'); title('Fixed EO/gap amplitude objective profile');
exportgraphics(fig,fullfile(out,'low_increment_amplitude_profile.png'),'Resolution',180); close(fig);
end
