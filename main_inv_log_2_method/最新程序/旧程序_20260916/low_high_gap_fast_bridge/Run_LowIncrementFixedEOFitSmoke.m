function R=Run_LowIncrementFixedEOFitSmoke()
% Targeted diagnostic for the fixed-EO amplitude re-fit.
root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2); cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=20; cfg.NumRevs_high=8; cfg.fs=base.cfgCase.fs;
eoTrue=11; gLow=.50; gHigh=.70; cfg.A_true=.20; cfg.f_true=eoTrue*(cfg.RPM_high/60); cfg.phi_true=.4;
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg,lib.domain,Inf,951001,'fixed_std');
low=build_low_speed_templates_adaptive(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'cv_sg',struct());
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false);
[pathCal,~]=calibrate_low_speed_path(lib,low,opts); path=build_path_template_model(lib,low.templateLow,pathCal);
high=simulate_highspeed_from_low_increment(path,gHigh-pathCal.g0,cfg,Inf,'fixed_std',951002);
map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
bundle=struct('t',map.t_v,'x',map.x_v,'V',map.V_a,'theta',map.theta_v,'sensorId',map.S_v,'turnId',map.rev_v,'Wvp',ones(size(map.V_a)),'Wfull',ones(size(map.V_a)));
response=struct('evalF',@(dg,xq,sid)eval_path_increment_template(path,dg,xq,sid),'evalFx',@(dg,xq,sid)eval_gap_derivative(path,dg,xq,sid));
rotHz=cfg.RPM_high/60; bc=struct('gapBounds',[.10,.30],'gapGrid',.10:.02:.30,'eoGrid',eoTrue,'topK',1,'rotHz',rotHz,'minGradient',1e-8,'minSamples',12,'dxBound',.20,'ampCoeffBound',.8,'deltaFBound',2,'maxIter',40,'ambiguityMargin',.001,'refineAmplitudeGrid',0:.025:.30);
fit=run_unified_blind_dynamic_backend(bundle,response,bc);
best=fit.best; R=struct('g_low_true_mm',gLow,'g_low_cal_mm',pathCal.g0,'g_high_true_mm',gHigh,'delta_g_true_mm',gHigh-gLow,'delta_g_hat_mm',best.gap,'eo_true',eoTrue,'eo_hat',best.eo,'amplitude_true_mm',cfg.A_true,'amplitude_hat_mm',best.A,'phase_true_rad',cfg.phi_true,'phase_hat_rad',best.phi,'frequency_true_hz',cfg.f_true,'frequency_hat_hz',(best.eo+best.deltaF/rotHz)*rotHz,'rmse_V',best.fullWaveRmseV,'truth_rmse_V',sqrt(mean((map.V_a-eval_path_increment_template(path,gHigh-pathCal.g0,map.x_v-cfg.A_true*sin(eoTrue*map.theta_v+cfg.phi_true),map.S_v)).^2)),'method',"fixed_eo_refit_targeted");
out=fullfile(root,'output','low_increment_fixed_eo_fit_smoke'); if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(R),fullfile(out,'low_increment_fixed_eo_fit_smoke.csv')); save(fullfile(out,'low_increment_fixed_eo_fit_smoke.mat'),'R','fit'); disp(struct2table(R));
end
