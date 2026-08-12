function [Detail,Summary,CalibrationAudit]=Run_AdaptiveSGMainRegression(snrList,seeds)
%RUN_ADAPTIVESGMAINREGRESSION Main-entry regression with reused calibration.
if nargin<1||isempty(snrList),snrList=[25 15 5];end
if nargin<2||isempty(seeds),seeds=1:3;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.route30LowTemplateMethod="adaptive_sg";
cfg0.route30BaseSamples=900;cfg0.route30EnableRhoFallback=true;
gLow=.5;cases=struct('name',{'single_531_g08','dual_700_1200_g02'},...
    'gHigh',{.8,.2},'f',{531,[700 1200]},'A',{.25,[.25 .15]},...
    'phi',{pi/4,[pi/4 -pi/3]});
nCal=numel(snrList)*numel(seeds);n=nCal*numel(cases);
Detail=repmat(row0(),n,1);CalibrationAudit=repmat(cal0(),nCal,1);ir=0;ic=0;
outDir=fullfile(root,'output','adaptive_sg_main_regression');if ~exist(outDir,'dir'),mkdir(outDir);end
for is=1:numel(snrList)
    snrDb=snrList(is);
    for js=1:numel(seeds)
        seed=seeds(js);cfg=cfg0;cfg.snrDb=snrDb;
        lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
            cfg,lib.domain,snrDb,941000+seed);
        tc=tic;C=calibrate_inv_log_2_low_speed(lowData,lib,cfg);calTime=toc(tc);
        ic=ic+1;CalibrationAudit(ic)=pack_cal(snrDb,seed,C,calTime,gLow,lib);
        for k=1:numel(cases)
            cfgCase=cfg;cfgCase.f_true=cases(k).f;cfgCase.A_true=cases(k).A;
            cfgCase.phi_true=cases(k).phi;
            high=simulate_highspeed_from_low_increment(C.templateModel,...
                cases(k).gHigh-C.pathCal.g0,cfgCase,snrDb,'snr_db',...
                942000+100*k+seed);
            hm=map_highspeed_to_space(high,cfgCase.alpha_k,cfgCase.R_tip,lib.domain,.02);
            th=tic;fit=run_inv_log_2_high_with_calibration(C,hm,cfgCase);highTime=toc(th);
            ir=ir+1;Detail(ir)=pack_row(snrDb,seed,cases(k),fit,C,calTime,highTime);
            writetable(struct2table(Detail(1:ir)),fullfile(outDir,'checkpoint_detail.csv'));
            save(fullfile(outDir,'checkpoint.mat'),'Detail','CalibrationAudit','ir','ic');
        end
    end
end
T=struct2table(Detail);Ctab=struct2table(CalibrationAudit);Summary=make_summary(T);
writetable(T,fullfile(outDir,'detail.csv'));writetable(Summary,fullfile(outDir,'summary.csv'));
writetable(Ctab,fullfile(outDir,'calibration_audit.csv'));
save(fullfile(outDir,'regression.mat'),'Detail','Summary','CalibrationAudit','snrList','seeds');
disp(Summary);
end
function r=row0()
r=struct('snr_db',NaN,'seed',NaN,'case_name','', 'calibration_key','',...
    'selected_window_mm',NaN,'selected_span_bins',NaN,'calibration_time_s',NaN,...
    'high_fit_time_s',NaN,'true_order',NaN,'estimated_order',NaN,...
    'max_frequency_error_hz',NaN,'g_true_mm',NaN,'g_est_mm',NaN,...
    'g_error_mm',NaN,'max_amplitude_error_mm',NaN,'rmse_V',NaN,...
    'model_order_correct',false,'frequency_correct',false,'gap_correct',false,...
    'success',false);
end
function r=pack_row(snr,seed,c,fit,C,tc,th)
f0=sort(c.f(:).');fh=sort(fit.f_id(:).');r=row0();r.snr_db=snr;r.seed=seed;
r.case_name=c.name;r.calibration_key=sprintf('%gdb_seed%d',snr,seed);
r.selected_window_mm=C.low_speed_selected_window_mm;
r.selected_span_bins=C.low_speed_selected_span_bins;r.calibration_time_s=tc;
r.high_fit_time_s=th;r.true_order=numel(f0);r.estimated_order=fit.model_order;
if numel(fh)==numel(f0),r.max_frequency_error_hz=max(abs(fh-f0));end
r.g_true_mm=c.gHigh;r.g_est_mm=fit.g_used;r.g_error_mm=abs(fit.g_used-c.gHigh);
if numel(fit.A_id)==numel(c.A),r.max_amplitude_error_mm=max(abs(sort(fit.A_id)-sort(c.A)));end
r.rmse_V=fit.rmse;r.model_order_correct=fit.model_order==numel(f0);
r.frequency_correct=r.max_frequency_error_hz<=2;r.gap_correct=r.g_error_mm<=.05;
r.success=r.model_order_correct&&r.frequency_correct&&r.gap_correct;
end
function r=cal0()
r=struct('snr_db',NaN,'seed',NaN,'selected_window_mm',NaN,...
    'selected_span_bins',NaN,'selection_at_boundary',false,'g0_error_mm',NaN,...
    'template_rmse_V',NaN,'calibration_time_s',NaN,'reused_high_record_count',2);
end
function r=pack_cal(snr,seed,C,t,g0,lib)
r=cal0();r.snr_db=snr;r.seed=seed;r.selected_window_mm=C.low_speed_selected_window_mm;
r.selected_span_bins=C.low_speed_selected_span_bins;cv=C.low_speed_template_cv;
r.selection_at_boundary=r.selected_span_bins==cv.span_bins(1)||r.selected_span_bins==cv.span_bins(end);
r.g0_error_mm=abs(C.pathCal.g0-g0);truth=eval_gap_template(lib,g0,lib.xGrid);
r.template_rmse_V=sqrt(mean((C.lowTemplate.templateLow-truth).^2));r.calibration_time_s=t;
end
function S=make_summary(T)
u=unique(T(:,{'snr_db','case_name'}),'rows','stable');S=u;
vars={'max_frequency_error_hz','g_error_mm','max_amplitude_error_mm','rmse_V','high_fit_time_s'};
S.n=zeros(height(S),1);S.success_rate=zeros(height(S),1);
for v=vars,S.([v{1},'_median'])=zeros(height(S),1);S.([v{1},'_iqr'])=zeros(height(S),1);end
for i=1:height(S)
    idx=T.snr_db==S.snr_db(i)&strcmp(T.case_name,S.case_name{i});S.n(i)=nnz(idx);
    S.success_rate(i)=mean(T.success(idx));for v=vars,q=T.(v{1})(idx);
        S.([v{1},'_median'])(i)=median(q);S.([v{1},'_iqr'])(i)=iqr(q);end
end
end
