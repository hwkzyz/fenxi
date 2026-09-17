function Result = Run_FormalNoiseOffGridRegression(snrList,seedList,useInformationSelection)
%RUN_FORMALNOISEOFFGRIDREGRESSION Off-grid and SNR regression for main route.
if nargin<1||isempty(snrList),snrList=[20 10 5];end
if nargin<2||isempty(seedList),seedList=1;end
if nargin<3,useInformationSelection=false;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=8;
cfg0.route30BaseSamples=900;cfg0.route30ForwardModel="absolute";
cfg0.route30UseInformationSelection=logical(useInformationSelection);
freqCases={[531 1300],[733 1217],[917 1383]};
gapPairs=[.2 .8;.8 .2;.5 .8];A=[.25 .15];phi=[pi/4 -pi/3];
rows=repmat(row0(),0,1);
for isnr=1:numel(snrList)
    for is=1:numel(seedList)
        for ic=1:numel(freqCases)
            cfg=cfg0;cfg.snrDb=snrList(isnr);seed=seedList(is);fTrue=freqCases{ic};
            gLow=gapPairs(ic,1);gHigh=gapPairs(ic,2);
            lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
                cfg,lib.domain,cfg.snrDb,710000+10000*seed+100*isnr+ic);
            lowAgg=aggregate_low_speed_template(lowData,cfg.alpha_k,cfg.R_tip,...
                lib.domain,lib.xGrid);
            lowState=estimate_highspeed_static_gap_raw(lowAgg.mapped,lib,cfg);
            high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
                cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
                A,fTrue,phi,'noise_ratio');
            sig=rms(high.V_clean-min(high.V_clean));
            rng(720000+10000*seed+100*isnr+ic,'twister');
            high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
            hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
            ticCase=tic;fit=run_low_high_gap_fast_vp_method(hm,lib,cfg,lowState);elapsed=toc(ticCase);
            fEst=sort(fit.f_id(:).');r=row0();r.snr_db=cfg.snrDb;r.seed=seed;
            r.case_name=sprintf('%d+%d',fTrue);r.f1_true=fTrue(1);r.f2_true=fTrue(2);
            r.f1_est=fEst(1);r.f2_est=fEst(2);r.max_frequency_error=max(abs(fEst-fTrue));
            r.g_true=gHigh;r.g_est=fit.g_used;r.gap_error=fit.g_used-gHigh;
            r.rmse=fit.rmse;r.elapsed_s=elapsed;r.fallback_used=fit.fallback_used;
            r.used_samples=fit.used_sample_count;r.success=fit.model_order==2&&...
                r.max_frequency_error<=2&&abs(r.gap_error)<=.05;
            rows(end+1,1)=r; %#ok<AGROW>
            fprintf('%g dB seed %d %s: f=%g,%g err=%g t=%g fallback=%d success=%d\n',...
                r.snr_db,seed,r.case_name,r.f1_est,r.f2_est,r.max_frequency_error,...
                elapsed,r.fallback_used,r.success);
        end
    end
end
Audit=struct2table(rows);out=fullfile(root,'output','formal_noise_offgrid');
if~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'formal_noise_offgrid.csv'));
save(fullfile(out,'formal_noise_offgrid.mat'),'Audit');
Result=struct('Audit',Audit,'outputDir',out,'pass',all(Audit.success));disp(Audit);
end

function r=row0()
r=struct('snr_db',NaN,'seed',NaN,'case_name',"",'f1_true',NaN,'f2_true',NaN,...
    'f1_est',NaN,'f2_est',NaN,'max_frequency_error',NaN,'g_true',NaN,...
    'g_est',NaN,'gap_error',NaN,'rmse',NaN,'elapsed_s',NaN,...
    'fallback_used',false,'used_samples',NaN,'success',false);
end
