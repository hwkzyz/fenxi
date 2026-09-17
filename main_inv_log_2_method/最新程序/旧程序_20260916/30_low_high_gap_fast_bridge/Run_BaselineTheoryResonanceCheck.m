function Result=Run_BaselineTheoryResonanceCheck(snrDb,A,delta)
%RUN_BASELINETHEORYRESONANCECHECK Single-frequency low-increment closure test.
if nargin<1,snrDb=Inf;end
if nargin<2,A=.20;end
if nargin<3,delta=.10;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=8;
cfg.snrDb=snrDb;cfg.route30ForwardModel="low_increment";
cfg.route30SingleFrequencyRangeHz=[500 900];cfg.route30UseAllTrustedSamples=true;
cfg.A_true=A;cfg.f_true=700;cfg.phi_true=pi/4;
gLow=.50;gHigh=gLow+delta;
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
    cfg,lib.domain,snrDb,931001);
low=aggregate_low_speed_template(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],...
    'fitGainOffset',false,'g0Init',gLow);
[pathCal,calFit]=calibrate_low_speed_path(lib,low,opts);
pathModel=build_path_template_model(lib,low.templateLow,pathCal);
high=simulate_highspeed_from_low_increment(pathModel,delta,cfg,snrDb,'snr_db',932001);
map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);
state.gHat=pathCal.g0;state.g_low_hat=pathCal.g0;state.dx0=0;
ticCase=tic;fit=run_single_sync_voltage_vp(map,pathModel,cfg,state);elapsed=toc(ticCase);
Result=struct('snr_db',snrDb,'f_true',cfg.f_true,'f_est',fit.f_id,...
    'frequency_error_hz',abs(fit.f_id-cfg.f_true),'g_low_true',gLow,...
    'g_low_hat',pathCal.g0,'g_high_true',gHigh,'g_high_hat',fit.g_used,...
    'gap_error_mm',abs(fit.g_used-gHigh),'A_true',cfg.A_true,...
    'A_est',fit.A_id,'phase_true',cfg.phi_true,'phase_est',fit.phi_id,...
    'rmse',fit.rmse,'identification_status',get_status(fit),...
    'calibration_rmse_V',calFit.rmse_V,'elapsed_s',elapsed);
out=fullfile(root,'output',sprintf('baseline_theory_resonance_%s',snr_tag(snrDb)));
if ~exist(out,'dir'),mkdir(out);end
writetable(struct2table(Result),fullfile(out,'result.csv'));save(fullfile(out,'result.mat'),'Result');disp(struct2table(Result));
end

function s=get_status(f)
if isfield(f,'identification_status'),s=f.identification_status;else,s="single_candidate";end
end

function s=snr_tag(x)
if isinf(x),s='noiseless';else,s=strrep(sprintf('%gdb',x),'-','m');end
end
