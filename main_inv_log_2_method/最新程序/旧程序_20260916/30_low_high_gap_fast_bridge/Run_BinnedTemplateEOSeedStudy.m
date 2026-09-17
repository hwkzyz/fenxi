function [Detail,Summary]=Run_BinnedTemplateEOSeedStudy(numSeeds,lowSnrDb,highSnrDb,numLowRevs,numHighRevs,methodList,smoothSpan)
%RUN_BINNEDTEMPLATEEOSEEDSTUDY Full EO search for the three binned methods.
% Noise-only repeatability study with identical low/high realizations within
% each seed. EO remains unknown to the solver.
if nargin<1||isempty(numSeeds),numSeeds=5;end
if nargin<2||isempty(lowSnrDb),lowSnrDb=15;end
if nargin<3||isempty(highSnrDb),highSnrDb=15;end
if nargin<4||isempty(numLowRevs),numLowRevs=20;end
if nargin<5||isempty(numHighRevs),numHighRevs=8;end
if nargin<6||isempty(methodList),methodList={'mean','median','huber'};end
if nargin<7||isempty(smoothSpan),smoothSpan=9;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=numLowRevs;cfg0.NumRevs_high=numHighRevs;cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];cfg0.route30DualFrequencyRangeHz=[min(cfg0.f_true) max(cfg0.f_true)];cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;cfg0.dualSyncIteratedReplayCount=200;cfg0.dualSyncRefineCount=20;cfg0.structuredComponentAmplitudeFloorMm=.05;
cfg0.returnDualSyncDiagnostics=true;
g0=.50;delta=.10;truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg0.alpha_k));truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg0.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);
methods=cellstr(methodList);empty=empty_row();Detail=repmat(empty,numSeeds*numel(methods),1);ir=0;
for iseed=1:numSeeds
    cfg=cfg0;lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,lowSnrDb,941000+iseed);high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,highSnrDb,'snr_db',942000+iseed);highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    for im=1:numel(methods)
        low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,methods{im},struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',smoothSpan,'minCoverage',.90));opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);ir=ir+1;Detail(ir)=pack_row(iseed,methods{im},low,truthT,lib.xGrid,pc,fit,g0,delta,cfg.A_true,eoTrue,cal);
    end
end
Tdetail=struct2table(Detail);u=unique(string(Tdetail.method),'stable');Summary=table(u,'VariableNames',{'method'});vars={'template_rmse_V','weighted_derivative_rmse_V_per_mm','g0_error_mm','delta_gap_error_mm','A1_error_mm','A2_error_mm','eo_error','success','engineering_success','identification_confident'};for iv=1:numel(vars),v=vars{iv};med=zeros(numel(u),1);iq=zeros(numel(u),1);for iu=1:numel(u),q=Tdetail.(v)(strcmp(string(Tdetail.method),u(iu)));med(iu)=median(q);iq(iu)=iqr(q);end;Summary.([v,'_median'])=med;Summary.([v,'_iqr'])=iq;end
outDir=fullfile(root,'output','binned_template_eo_seed_study');if ~exist(outDir,'dir'),mkdir(outDir);end;tag=sprintf('low_%gdb_high_%gdb_%dLrevs_%dHrevs_%dseeds_%dmethods_SG%d',lowSnrDb,highSnrDb,numLowRevs,numHighRevs,numSeeds,numel(methods),smoothSpan);writetable(struct2table(Detail),fullfile(outDir,['detail_',tag,'.csv']));writetable(Summary,fullfile(outDir,['summary_',tag,'.csv']));save(fullfile(outDir,['study_',tag,'.mat']),'Detail','Summary','truthT','methods','smoothSpan');disp(Summary);
end

function r=empty_row()
r=struct('seed',NaN,'method','','template_rmse_V',NaN,'derivative_rmse_V_per_mm',NaN,'weighted_derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,'g0_error_mm',NaN,'delta_gap_error_mm',NaN,'eo1_est',NaN,'eo2_est',NaN,'eo_error',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,'success',false,'engineering_success',false,'identification_confident',false,'noise_normalized_margin',NaN,'competitor_eo1',NaN,'competitor_eo2',NaN,'true_pair_bank_rank',NaN,'true_pair_refined',false,'true_pair_best_refined_rmse_V',NaN,'rmse_V',NaN,'calibration_rmse_V',NaN);
end

function r=pack_row(seed,name,low,Ttrue,x,pc,fit,g0,delta,Atrue,eotrue,cal)
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);de=D-Dt;w=Dt.^2;eo=fit.eo_id(:);ae=fit.A_id(:);r=empty_row();r.seed=seed;r.method=name;r.template_rmse_V=sqrt(mean((low.templateBySensor-Ttrue).^2,'all'));r.derivative_rmse_V_per_mm=sqrt(mean(de.^2,'all'));r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*de.^2,'all')/sum(w,'all'));r.derivative_correlation=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));r.g0_error_mm=abs(pc.g0-g0);r.delta_gap_error_mm=abs((fit.g_used-pc.g0)-delta);r.eo1_est=eo(1);r.eo2_est=eo(2);r.eo_error=max(abs(eo-eotrue(:)));r.A1_error_mm=abs(ae(1)-Atrue(1));r.A2_error_mm=abs(ae(2)-Atrue(2));r.A1_est_mm=ae(1);r.A2_est_mm=ae(2);r.success=r.eo_error==0&&r.delta_gap_error_mm<=.01&&max([r.A1_error_mm,r.A2_error_mm])<=.01;r.engineering_success=r.eo_error==0&&r.delta_gap_error_mm<=.01&&all([r.A1_error_mm/Atrue(1),r.A2_error_mm/Atrue(2)]<=.10);if isfield(fit,'identification_confident'),r.identification_confident=logical(fit.identification_confident);end;if isfield(fit,'noise_normalized_margin'),r.noise_normalized_margin=fit.noise_normalized_margin;end;if isfield(fit,'competing_frequency_or_order')&&numel(fit.competing_frequency_or_order)==2,fRot=fit.f_id(1)/fit.eo_id(1);ceo=sort(fit.competing_frequency_or_order(:)/fRot);r.competitor_eo1=ceo(1);r.competitor_eo2=ceo(2);end;if isfield(fit,'diagnostic_bank_eo'),idx=find(all(fit.diagnostic_bank_eo==eotrue,2),1);if ~isempty(idx),r.true_pair_bank_rank=idx;end;end;if isfield(fit,'diagnostic_refined_eo'),idx=all(fit.diagnostic_refined_eo==eotrue,2);r.true_pair_refined=any(idx);if any(idx),r.true_pair_best_refined_rmse_V=min(fit.diagnostic_refined_rmse(idx));end;end;r.rmse_V=fit.rmse;r.calibration_rmse_V=cal.rmse_V;
end

function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
