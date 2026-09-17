function [Detail,Summary]=Run_NoiseSourceAttribution(snrList,seeds)
%RUN_NOISESOURCEATTRIBUTION Separate low-template and high-record noise effects.
% Runs only the two mixed cells omitted by the legacy/controlled matrices:
% legacy low + fixed high, and fixed low + legacy high.
if nargin<1||isempty(snrList),snrList=[15 5];end
if nargin<2||isempty(seeds),seeds=1:3;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
b=make_inv_log_2_demo_case(ctx,.2);cfg0=b.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);
cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;cfg0.route30ForwardModel="low_increment";
cfg0.route30LowTemplateMethod="adaptive_sg";cfg0.route30UseFastDualCandidates=false;
cfg0.f_true=[700 1200];cfg0.A_true=[.25 .15];cfg0.phi_true=[pi/4 -pi/3];
pairs=[.5 .4;.8 .5];modes=["legacy_low_fixed_high","fixed_low_legacy_high"];
ref=readtable(fullfile(root,'output','noise_protocol_gap_matrix','reference_noise.csv'));
n=numel(snrList)*numel(seeds)*size(pairs,1)*numel(modes);Detail=repmat(row0(),n,1);ir=0;
out=fullfile(root,'output','noise_source_attribution');if ~exist(out,'dir'),mkdir(out);end
for snr=snrList
    hit=find(abs(ref.reference_vibration_snr_db-snr)<eps,1);if isempty(hit),error('No reference sigma for %g dB.',snr);end
    sigma=ref.fixed_noise_std_V(hit);
    for seed=seeds
        cfg=cfg0;cfg.snrDb=snr;
        for ip=1:size(pairs,1)
            gLow=pairs(ip,1);gHigh=pairs(ip,2);
            for im=1:numel(modes)
                if modes(im)=="legacy_low_fixed_high"
                    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg,lib.domain,snr,971000+100*ip+seed);
                else
                    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg,lib.domain,sigma,972000+100*ip+seed,'fixed_std');
                end
                C=calibrate_inv_log_2_low_speed(low,lib,cfg);delta=gHigh-C.pathCal.g0;
                if modes(im)=="legacy_low_fixed_high"
                    high=simulate_highspeed_from_low_increment(C.templateModel,delta,cfg,sigma,'fixed_std',973000+1000*im+100*ip+seed);
                else
                    high=simulate_highspeed_from_low_increment(C.templateModel,delta,cfg,snr,'snr_db',973000+1000*im+100*ip+seed);
                end
                hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);tt=tic;
                fit=run_inv_log_2_high_with_calibration(C,hm,cfg);elapsed=toc(tt);ir=ir+1;
                Detail(ir)=pack(modes(im),snr,seed,gLow,gHigh,C,low.noise_std,high,fit,elapsed);
                writetable(struct2table(Detail(1:ir)),fullfile(out,'checkpoint_detail.csv'));
            end
        end
    end
end
T=struct2table(Detail);Summary=groupsummary(T,{'mode','snr_db','g_low_mm','g_high_mm'},...
    {'mean','median','max'},{'success','max_frequency_error_hz','g_error_mm','noise_std_V'});
writetable(T,fullfile(out,'detail.csv'));writetable(Summary,fullfile(out,'summary.csv'));
save(fullfile(out,'noise_source_attribution.mat'),'Detail','Summary','snrList','seeds');disp(Summary);
end
function r=row0()
r=struct('mode',"",'snr_db',NaN,'seed',NaN,'g_low_mm',NaN,'g_high_mm',NaN,...
    'low_window_mm',NaN,'low_noise_std_V',NaN,'noise_std_V',NaN,'fit_time_s',NaN,...
    'max_frequency_error_hz',NaN,'g_error_mm',NaN,'max_amplitude_error_mm',NaN,'success',false);
end
function r=pack(mode,snr,seed,gLow,gHigh,C,lowNoise,high,fit,elapsed)
r=row0();f=sort(fit.f_id(:));r.mode=mode;r.snr_db=snr;r.seed=seed;r.g_low_mm=gLow;r.g_high_mm=gHigh;
r.low_window_mm=C.low_speed_selected_window_mm;r.low_noise_std_V=lowNoise;r.noise_std_V=high.noise_std;
r.fit_time_s=elapsed;if numel(f)==2,r.max_frequency_error_hz=max(abs(f-[700;1200]));end
r.g_error_mm=abs(fit.g_used-gHigh);if numel(fit.A_id)==2,r.max_amplitude_error_mm=max(abs(sort(fit.A_id)-[.15 .25]));end
r.success=fit.model_order==2&&r.max_frequency_error_hz<=2&&r.g_error_mm<=.05;
end
