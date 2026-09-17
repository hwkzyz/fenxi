function [Detail,Summary,CalibrationAudit]=Run_NoiseProtocolGapMatrix(snrList,seeds)
%RUN_NOISEPROTOCOLGAPMATRIX Compare absolute-noise and vibration-SNR protocols.
% The SG parameter is still selected from low-speed records only.  For the
% fixed-absolute protocol, one voltage-noise standard deviation is calibrated
% at g_ref=0.5 mm and then held fixed for every low/high clearance pair.
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
cfg0.route30UseFastDualCandidates=false;
cfg0.f_true=[700 1200];cfg0.A_true=[.25 .15];cfg0.phi_true=[pi/4 -pi/3];
gapPairs=[.5 .4;.5 .6;.5 .8;.8 .5];uLow=unique(gapPairs(:,1),'stable');
protocols=["fixed_absolute","fixed_vibration_snr"];
sigmaRef=reference_vibration_noise(lib,cfg0,snrList,.5);
nCal=numel(snrList)*numel(seeds)*numel(uLow);
n=numel(protocols)*numel(snrList)*numel(seeds)*size(gapPairs,1);
Detail=repmat(row0(),n,1);CalibrationAudit=repmat(cal0(),nCal,1);ir=0;ic=0;
out=fullfile(root,'output','noise_protocol_gap_matrix');if ~exist(out,'dir'),mkdir(out);end
for isnr=1:numel(snrList)
    snr=snrList(isnr);sigmaAbs=sigmaRef(isnr);
    for seed=seeds
        cfg=cfg0;cfg.snrDb=snr;
        for il=1:numel(uLow)
            gLow=uLow(il);lowSeed=961000+1000*(il-1)+seed;
            low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
                cfg,lib.domain,Inf,lowSeed);
            low=inject_noise(low,sigmaAbs,lowSeed,"fixed_absolute",snr);
            tc=tic;C=calibrate_inv_log_2_low_speed(low,lib,cfg);calTime=toc(tc);
            ic=ic+1;CalibrationAudit(ic)=pack_cal(snr,seed,gLow,C,calTime,sigmaAbs);
            pairIdx=find(gapPairs(:,1)==gLow).';
            for ip=pairIdx
                gHigh=gapPairs(ip,2);delta=gHigh-C.pathCal.g0;
                clean=simulate_highspeed_from_low_increment(C.templateModel,delta,cfg,Inf,'snr_db',1);
                sigmaVib=vibration_noise_std(clean,C.templateModel,delta,cfg,lib.domain,snr);
                for im=1:numel(protocols)
                    protocol=protocols(im);if protocol=="fixed_absolute",sigma=sigmaAbs;else,sigma=sigmaVib;end
                    highSeed=965000+1000*im+100*ip+seed;
                    high=inject_noise(clean,sigma,highSeed,protocol,snr);
                    hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
                    th=tic;
                    try
                        fit=run_inv_log_2_high_with_calibration(C,hm,cfg);fitTime=toc(th);
                        row=pack_row(protocol,snr,seed,gLow,gHigh,C,fit,calTime,fitTime,sigma,sigmaAbs,sigmaVib);
                    catch ME
                        fitTime=toc(th);row=pack_error(protocol,snr,seed,gLow,gHigh,C,calTime,fitTime,sigma,sigmaAbs,sigmaVib,ME);
                    end
                    ir=ir+1;Detail(ir)=row;
                    writetable(struct2table(Detail(1:ir)),fullfile(out,'checkpoint_detail.csv'));
                    save(fullfile(out,'checkpoint.mat'),'Detail','CalibrationAudit','ir','ic','sigmaRef');
                end
            end
        end
    end
end
T=struct2table(Detail);Ctab=struct2table(CalibrationAudit);Summary=make_summary(T);
writetable(T,fullfile(out,'detail.csv'));writetable(Summary,fullfile(out,'summary.csv'));
writetable(Ctab,fullfile(out,'calibration_audit.csv'));
ReferenceNoise=table(snrList(:),sigmaRef(:),'VariableNames',{'reference_vibration_snr_db','fixed_noise_std_V'});
writetable(ReferenceNoise,fullfile(out,'reference_noise.csv'));
save(fullfile(out,'noise_protocol_gap_matrix.mat'),'Detail','Summary','CalibrationAudit',...
    'ReferenceNoise','gapPairs','snrList','seeds');disp(Summary);
end

function sigma=reference_vibration_noise(lib,cfg,snrList,gRef)
vib=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gRef,z),...
    cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
    cfg.A_true,cfg.f_true,cfg.phi_true,'noise_ratio');
stat=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gRef,z),...
    cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
    [],[],[],'noise_ratio');
mv=map_highspeed_to_space(vib,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
ms=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
vibRms=rms(mv.V_a-ms.V_a);sigma=vibRms./10.^(snrList/20);
end

function sigma=vibration_noise_std(clean,pathModel,delta,cfg,domain,snr)
cfgStatic=cfg;cfgStatic.A_true=[];cfgStatic.f_true=[];cfgStatic.phi_true=[];
stat=simulate_highspeed_from_low_increment(pathModel,delta,cfgStatic,Inf,'snr_db',1);
mv=map_highspeed_to_space(clean,cfg.alpha_k,cfg.R_tip,domain,.02);
ms=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,domain,.02);
sigma=rms(mv.V_a-ms.V_a)/10^(snr/20);
end

function data=inject_noise(data,sigma,seed,mode,snr)
rng(seed,'twister');data.V_cap=data.V_clean+sigma*randn(size(data.V_clean));
data.noise_std=sigma;data.noise_mode=char(mode);data.noise_level=sigma;
data.snr_db_equiv=snr;data.random_seed=seed;
end

function r=row0()
r=struct('protocol',"",'snr_db',NaN,'seed',NaN,'g_low_mm',NaN,'g_high_mm',NaN,...
    'delta_gap_mm',NaN,'window_mm',NaN,'span_bins',NaN,'noise_std_V',NaN,...
    'reference_noise_std_V',NaN,'case_vibration_noise_std_V',NaN,...
    'calibration_time_s',NaN,'fit_time_s',NaN,'f1_est_hz',NaN,'f2_est_hz',NaN,...
    'max_frequency_error_hz',NaN,'g_est_mm',NaN,'g_error_mm',NaN,...
    'max_amplitude_error_mm',NaN,'rmse_V',NaN,'error_identifier',"",...
    'error_message',"",'success',false);
end
function r=pack_row(protocol,snr,seed,gLow,gHigh,C,fit,tc,tf,sigma,sigmaAbs,sigmaVib)
r=row0();f=sort(fit.f_id(:));r.protocol=protocol;r.snr_db=snr;r.seed=seed;
r.g_low_mm=gLow;r.g_high_mm=gHigh;r.delta_gap_mm=gHigh-gLow;
r.window_mm=C.low_speed_selected_window_mm;r.span_bins=C.low_speed_selected_span_bins;
r.noise_std_V=sigma;r.reference_noise_std_V=sigmaAbs;r.case_vibration_noise_std_V=sigmaVib;
r.calibration_time_s=tc;r.fit_time_s=tf;if numel(f)==2,r.f1_est_hz=f(1);r.f2_est_hz=f(2);r.max_frequency_error_hz=max(abs(f-[700;1200]));end
r.g_est_mm=fit.g_used;r.g_error_mm=abs(fit.g_used-gHigh);
if numel(fit.A_id)==2,r.max_amplitude_error_mm=max(abs(sort(fit.A_id)-[.15 .25]));end
r.rmse_V=fit.rmse;r.success=fit.model_order==2&&r.max_frequency_error_hz<=2&&r.g_error_mm<=.05;
end
function r=pack_error(protocol,snr,seed,gLow,gHigh,C,tc,tf,sigma,sigmaAbs,sigmaVib,ME)
r=row0();r.protocol=protocol;r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;r.g_high_mm=gHigh;
r.delta_gap_mm=gHigh-gLow;r.window_mm=C.low_speed_selected_window_mm;r.span_bins=C.low_speed_selected_span_bins;
r.noise_std_V=sigma;r.reference_noise_std_V=sigmaAbs;r.case_vibration_noise_std_V=sigmaVib;
r.calibration_time_s=tc;r.fit_time_s=tf;r.error_identifier=string(ME.identifier);r.error_message=string(ME.message);
end
function r=cal0()
r=struct('snr_db',NaN,'seed',NaN,'g_low_mm',NaN,'window_mm',NaN,'span_bins',NaN,...
    'low_noise_std_V',NaN,'g0_error_mm',NaN,'selection_at_boundary',false,...
    'calibration_time_s',NaN,'reused_record_count',NaN);
end
function r=pack_cal(snr,seed,gLow,C,t,sigma)
r=cal0();r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;r.window_mm=C.low_speed_selected_window_mm;
r.span_bins=C.low_speed_selected_span_bins;r.low_noise_std_V=sigma;r.g0_error_mm=abs(C.pathCal.g0-gLow);
cv=C.low_speed_template_cv;r.selection_at_boundary=r.span_bins==cv.span_bins(1)||r.span_bins==cv.span_bins(end);
r.calibration_time_s=t;r.reused_record_count=1+2*(gLow==.5);
end
function S=make_summary(T)
S=unique(T(:,{'protocol','snr_db','g_low_mm','g_high_mm','delta_gap_mm'}),'rows','stable');
vars={'max_frequency_error_hz','g_error_mm','max_amplitude_error_mm','rmse_V','fit_time_s','noise_std_V'};
S.n=zeros(height(S),1);S.success_rate=zeros(height(S),1);
for v=vars,S.([v{1},'_median'])=zeros(height(S),1);S.([v{1},'_max'])=zeros(height(S),1);end
for i=1:height(S)
    q=T.protocol==S.protocol(i)&T.snr_db==S.snr_db(i)&T.g_low_mm==S.g_low_mm(i)&T.g_high_mm==S.g_high_mm(i);
    S.n(i)=nnz(q);S.success_rate(i)=mean(T.success(q));for v=vars,x=T.(v{1})(q);
        S.([v{1},'_median'])(i)=median(x);S.([v{1},'_max'])(i)=max(x);end
end
end
