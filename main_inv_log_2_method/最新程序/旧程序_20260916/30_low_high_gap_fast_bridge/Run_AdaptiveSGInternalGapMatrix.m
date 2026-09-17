function [Detail,Summary,CalibrationAudit]=Run_AdaptiveSGInternalGapMatrix(snrList,seeds)
%RUN_ADAPTIVESGINTERNALGAPMATRIX Dual-frequency internal-gap regression.
if nargin<1||isempty(snrList),snrList=[15 5];end
if nargin<2||isempty(seeds),seeds=1:3;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.route30LowTemplateMethod="adaptive_sg";
% Freeze the validated full unordered-pair search. A separate fast-candidate
% branch is under active development and is not part of this gap study.
cfg0.route30UseFastDualCandidates=false;
cfg0.f_true=[700 1200];cfg0.A_true=[.25 .15];cfg0.phi_true=[pi/4 -pi/3];
gapPairs=[.5 .4;.5 .6;.5 .8;.8 .5];uLow=unique(gapPairs(:,1),'stable');
n=numel(snrList)*numel(seeds)*size(gapPairs,1);nCal=numel(snrList)*numel(seeds)*numel(uLow);
Detail=repmat(row0(),n,1);CalibrationAudit=repmat(cal0(),nCal,1);ir=0;ic=0;
out=fullfile(root,'output','adaptive_sg_internal_gap_matrix');if ~exist(out,'dir'),mkdir(out);end
for snr=snrList
    for seed=seeds
        cfg=cfg0;cfg.snrDb=snr;
        for il=1:numel(uLow)
            gLow=uLow(il);lowSeed=941000+1000*(il-1)+seed;
            low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
                cfg,lib.domain,snr,lowSeed);
            tc=tic;C=calibrate_inv_log_2_low_speed(low,lib,cfg);calTime=toc(tc);
            ic=ic+1;CalibrationAudit(ic)=pack_cal(snr,seed,gLow,C,calTime);
            pairIdx=find(gapPairs(:,1)==gLow).';
            for ip=pairIdx
                gHigh=gapPairs(ip,2);high=simulate_highspeed_from_low_increment(...
                    C.templateModel,gHigh-C.pathCal.g0,cfg,snr,'snr_db',945000+100*ip+seed);
                hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
                th=tic;
                try
                    fit=run_inv_log_2_high_with_calibration(C,hm,cfg);fitTime=toc(th);
                    row=pack_row(snr,seed,gLow,gHigh,C,fit,calTime,fitTime);
                catch ME
                    fitTime=toc(th);row=pack_error(snr,seed,gLow,gHigh,C,calTime,fitTime,ME);
                end
                ir=ir+1;Detail(ir)=row;
                writetable(struct2table(Detail(1:ir)),fullfile(out,'checkpoint_detail.csv'));
                save(fullfile(out,'checkpoint.mat'),'Detail','CalibrationAudit','ir','ic');
            end
        end
    end
end
T=struct2table(Detail);Ctab=struct2table(CalibrationAudit);Summary=make_summary(T);
writetable(T,fullfile(out,'detail.csv'));writetable(Summary,fullfile(out,'summary.csv'));
writetable(Ctab,fullfile(out,'calibration_audit.csv'));
save(fullfile(out,'internal_gap_matrix.mat'),'Detail','Summary','CalibrationAudit','gapPairs','snrList','seeds');
disp(Summary);
end
function r=row0()
r=struct('snr_db',NaN,'seed',NaN,'g_low_mm',NaN,'g_high_mm',NaN,...
    'delta_gap_mm',NaN,'calibration_key','', 'window_mm',NaN,'span_bins',NaN,...
    'calibration_time_s',NaN,'fit_time_s',NaN,'f1_est_hz',NaN,'f2_est_hz',NaN,...
    'max_frequency_error_hz',NaN,'g_est_mm',NaN,'g_error_mm',NaN,...
    'max_amplitude_error_mm',NaN,'rmse_V',NaN,'error_identifier',"",...
    'error_message',"", 'success',false);
end
function r=pack_row(snr,seed,gLow,gHigh,C,fit,tc,tf)
r=row0();f=sort(fit.f_id(:));r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;
r.g_high_mm=gHigh;r.delta_gap_mm=gHigh-gLow;
r.calibration_key=sprintf('%gdb_seed%d_g%g',snr,seed,gLow);
r.window_mm=C.low_speed_selected_window_mm;r.span_bins=C.low_speed_selected_span_bins;
r.calibration_time_s=tc;r.fit_time_s=tf;if numel(f)==2,r.f1_est_hz=f(1);r.f2_est_hz=f(2);r.max_frequency_error_hz=max(abs(f-[700;1200]));end
r.g_est_mm=fit.g_used;r.g_error_mm=abs(fit.g_used-gHigh);
if numel(fit.A_id)==2,r.max_amplitude_error_mm=max(abs(sort(fit.A_id)-[.15 .25]));end
r.rmse_V=fit.rmse;r.success=fit.model_order==2&&r.max_frequency_error_hz<=2&&r.g_error_mm<=.05;
end
function r=pack_error(snr,seed,gLow,gHigh,C,tc,tf,ME)
r=row0();r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;r.g_high_mm=gHigh;
r.delta_gap_mm=gHigh-gLow;r.calibration_key=sprintf('%gdb_seed%d_g%g',snr,seed,gLow);
r.window_mm=C.low_speed_selected_window_mm;r.span_bins=C.low_speed_selected_span_bins;
r.calibration_time_s=tc;r.fit_time_s=tf;r.error_identifier=string(ME.identifier);
r.error_message=string(ME.message);
end
function r=cal0()
r=struct('snr_db',NaN,'seed',NaN,'g_low_mm',NaN,'window_mm',NaN,...
    'span_bins',NaN,'g0_error_mm',NaN,'selection_at_boundary',false,...
    'calibration_time_s',NaN,'reused_record_count',NaN);
end
function r=pack_cal(snr,seed,gLow,C,t)
r=cal0();r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;
r.window_mm=C.low_speed_selected_window_mm;r.span_bins=C.low_speed_selected_span_bins;
r.g0_error_mm=abs(C.pathCal.g0-gLow);cv=C.low_speed_template_cv;
r.selection_at_boundary=r.span_bins==cv.span_bins(1)||r.span_bins==cv.span_bins(end);
r.calibration_time_s=t;r.reused_record_count=1+2*(gLow==.5);
end
function S=make_summary(T)
S=unique(T(:,{'snr_db','g_low_mm','g_high_mm','delta_gap_mm'}),'rows','stable');
vars={'max_frequency_error_hz','g_error_mm','max_amplitude_error_mm','rmse_V','fit_time_s'};
S.n=zeros(height(S),1);S.success_rate=zeros(height(S),1);
for v=vars,S.([v{1},'_median'])=zeros(height(S),1);S.([v{1},'_max'])=zeros(height(S),1);end
for i=1:height(S)
    q=T.snr_db==S.snr_db(i)&T.g_low_mm==S.g_low_mm(i)&T.g_high_mm==S.g_high_mm(i);
    S.n(i)=nnz(q);S.success_rate(i)=mean(T.success(q));for v=vars,x=T.(v{1})(q);
        S.([v{1},'_median'])(i)=median(x);S.([v{1},'_max'])(i)=max(x);end
end
end
