function Result=Run_ModelOrderAmplitudeRegression()
%RUN_MODELORDERAMPLITUDEREGRESSION Amplitude-aware single/dual order audit.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;
cfg.route30ForwardModel="absolute";cfg.route30UseAllTrustedSamples=true;
cfg.route30HierarchicalOrderSelection=true;cfg.route30EnableRhoFallback=false;
gLow=.8;gHigh=.5;fTrue=[700 1200];phi=[pi/4 -pi/3];
snrList=[20 5];weakAmplitude=[.15 .08 .04];
rows=repmat(row0(),numel(snrList)*numel(weakAmplitude),1);at=0;
for snrDb=snrList
    cfg.snrDb=snrDb;
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
        cfg,lib.domain,snrDb,20260820+snrDb);
    for A2=weakAmplitude
        at=at+1;A=[.25 A2];
        d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
            cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
            A,fTrue,phi,'noise_ratio');
        sig=rms(d.V_clean-min(d.V_clean));rng(20260900+at,'twister');
        d.V_cap=d.V_clean+sig/10^(snrDb/20)*randn(size(d.V_clean));
        map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        fit=run_inv_log_2_low_high_main(low,map,lib,cfg);test=fit.second_frequency_test;
        rows(at).snr_db=snrDb;rows(at).A1_mm=A(1);rows(at).A2_mm=A(2);
        rows(at).amplitude_ratio=A2/A(1);rows(at).estimated_order=fit.model_order;
        rows(at).detected_frequency_hz=test.second_frequency_hz;
        rows(at).second_delta_bic=test.delta_bic;rows(at).second_amplitude_z=test.amplitude_z;
        rows(at).strong_second_frequency=test.strong_second_frequency;
        rows(at).dual_search_skipped=fit.dual_search_skipped;
        rows(at).used_samples=fit.used_sample_count;rows(at).elapsed_s=fit.elapsed_s;
        if fit.model_order==2
            rows(at).max_frequency_error_hz=max(abs(sort(fit.f_id)-fTrue));
        end
    end
end
Audit=struct2table(rows);out=fullfile(root,'output','model_order_amplitude');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'model_order_amplitude.csv'));
save(fullfile(out,'model_order_amplitude.mat'),'Audit');disp(Audit);
Result=struct('Audit',Audit,'outputDir',out);
end

function r=row0()
r=struct('snr_db',NaN,'A1_mm',NaN,'A2_mm',NaN,'amplitude_ratio',NaN,...
    'estimated_order',NaN,'detected_frequency_hz',NaN,'second_delta_bic',NaN,...
    'second_amplitude_z',NaN,'strong_second_frequency',false,...
    'dual_search_skipped',false,'max_frequency_error_hz',NaN,...
    'used_samples',NaN,'elapsed_s',NaN);
end
