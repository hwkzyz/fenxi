function Result=Run_LowIncrementDualSyncOrderCheck(highSnrDb,fixedKnownOrder,lowSnrDb,useRobustTemplate,numLowRevs)
%RUN_LOWINCREMENTDUALSYNCORDERCHECK EO-constrained dual-frequency closure.
if nargin<1,highSnrDb=Inf;end
if nargin<2,fixedKnownOrder=false;end
if nargin<3,lowSnrDb=highSnrDb;end
if nargin<4,useRobustTemplate=false;end
if nargin<5,numLowRevs=8;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=numLowRevs;cfg.NumRevs_high=8;
cfg.snrDb=highSnrDb;cfg.route30ForwardModel="low_increment";
fRot=cfg.RPM_high/60;eoTrue=[10 26];cfg.A_true=[.25 .15];
cfg.f_true=eoTrue*fRot;cfg.phi_true=[pi/4 -pi/3];
cfg.route30DualFrequencyRangeHz=[min(cfg.f_true) max(cfg.f_true)];
if fixedKnownOrder,cfg.dualSyncCandidateEOPairs=eoTrue;end
cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;
cfg.dualSyncKeepPerGap=80;cfg.dualSyncIteratedReplayCount=300;
cfg.dualSyncRefineCount=30;cfg.structuredComponentAmplitudeFloorMm=.05;
gLow=.50;delta=.10;gHigh=gLow+delta;
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
    cfg,lib.domain,lowSnrDb,941001);
if useRobustTemplate
    low=aggregate_low_speed_template_robust(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,...
        struct('minCoverage',.90,'shiftGrid',0,'smoothSpan',1));
    low.templateBySensor=repmat(low.templateLow,1,numel(cfg.alpha_k));
    low.sensorIds=1:numel(cfg.alpha_k);
else
    % Mirror the actual Step04 route: each sensor owns an independently
    % constructed, coverage-gated and smoothed low-speed template.
    low=build_low_speed_templates_step04(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,...
        struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));
end
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],...
    'fitGainOffset',false,'g0Init',gLow);
[pathCal,calFit]=calibrate_low_speed_path(lib,low,opts);
pathCal.sensorIds=low.sensorIds;
pathModel=build_path_template_model(lib,low.templateBySensor,pathCal);
% The low-speed acquisition and its reconstruction are not part of the
% high-speed physical truth.  Generate high-speed data from the noise-free
% rotor reference, then invert it with the estimated Step04 templates.
trueCal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',low.sensorIds);
trueTemplate=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(low.sensorIds));
truthModel=build_path_template_model(lib,trueTemplate,trueCal);
high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,highSnrDb,'snr_db',942001);
map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);
state.gHat=pathCal.g0;state.g_low_hat=pathCal.g0;state.dx0=0;
ticCase=tic;fit=run_dual_sync_voltage_vp(map,pathModel,cfg,state);elapsed=toc(ticCase);
eoErr=max(abs(fit.eo_id-eoTrue));
Result=struct('high_snr_db',highSnrDb,'low_snr_db',lowSnrDb,...
    'fixed_known_order',fixedKnownOrder,'robust_template',useRobustTemplate,...
    'template_method',ternary(useRobustTemplate,"legacy_robust_diagnostic","step04_sensorwise_binned_smoothed"),...
    'low_template_revolutions',numLowRevs,...
    'eo_true',string(mat2str(eoTrue)),...
    'eo_est',string(mat2str(fit.eo_id)),'max_eo_error',eoErr,...
    'g_low_true',gLow,'g_low_hat',pathCal.g0,'g_high_true',gHigh,...
    'g_high_hat',fit.g_used,'gap_error_mm',abs(fit.g_used-gHigh),...
    'A_true',string(mat2str(cfg.A_true)),'A_est',string(mat2str(fit.A_id)),...
    'rmse',fit.rmse,'candidate_count',fit.candidate_count,...
    'component_floor_mm',fit.component_amplitude_floor_mm,...
    'calibration_rmse_V',calFit.rmse_V,...
    'template_rmse_V',sqrt(mean((low.templateBySensor-repmat(trueTemplate(:,1),1,size(low.templateBySensor,2))).^2,'all')),...
    'template_derivative_rmse_V_per_mm',template_derivative_rmse(low.templateBySensor,trueTemplate,lib.xGrid),...
    'template_smooth_span',get_low_field(low,'smooth_span',NaN),...
    'elapsed_s',elapsed,...
    'success',eoErr==0&&abs(fit.g_used-gHigh)<=.05&&min(fit.A_id)>=.05);
modeTag=ternary(fixedKnownOrder,'fixed','enumerated');
out=fullfile(root,'output',sprintf('low_increment_dual_sync_order_high_%s_low_%s_%s',...
    snr_tag(highSnrDb),snr_tag(lowSnrDb),modeTag));
if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(Result),fullfile(out,'result.csv'));save(fullfile(out,'result.mat'),'Result','fit','pathCal','calFit');disp(struct2table(Result));
end

function value=get_low_field(low,name,defaultValue)
if isfield(low,name)&&~isempty(low.(name)),value=low.(name);else,value=defaultValue;end
end

function e=template_derivative_rmse(T,Ttrue,x)
d=median(diff(x(:)));D=gradient(T,d);Dt=gradient(Ttrue,d);
e=sqrt(mean((D-Dt).^2,'all'));
end

function s=snr_tag(x)
if isinf(x),s='noiseless';else,s=strrep(sprintf('%gdb',x),'-','m');end
end

function v=ternary(c,a,b)
if c,v=a;else,v=b;end
end
