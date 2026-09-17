function Result=Run_ModelOrderBoundaryAudit()
%RUN_MODELORDERBOUNDARYAUDIT Probe order-selection limits without deleting data.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.route30ForwardModel="absolute";
cfg.route30UseAllTrustedSamples=true;cfg.route30HierarchicalOrderSelection=true;
cfg.route30EnableRhoFallback=false;gLow=.8;gHigh=.5;phi=[pi/4 -pi/3];
cases={};
% Frequency separation boundary at a visible second amplitude.
for snr=[10 5]
    for sep=[20 50 100]
        cases{end+1}=struct('name',sprintf('dual_sep_%g_snr_%g',sep,snr),...
            'snr',snr,'f',[700 700+sep],'A',[.25 .15]); %#ok<AGROW>
    end
end
% Weak-component boundary at a well-separated frequency.
for snr=[5 0]
    cases{end+1}=struct('name',sprintf('dual_weak_snr_%g',snr),...
        'snr',snr,'f',[700 1200],'A',[.25 .02]); %#ok<AGROW>
end
% False-dual alarm check for a single component.
for snr=[5 0]
    cases{end+1}=struct('name',sprintf('single_snr_%g',snr),...
        'snr',snr,'f',531,'A',.25); %#ok<AGROW>
end
rows=repmat(row0(),numel(cases),1);
for k=1:numel(cases)
    c=cases{k};cfg.snrDb=c.snr;
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
        cfg,lib.domain,c.snr,203000+k);
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
        c.A,c.f,phi(1:numel(c.f)),'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(203100+k,'twister');
    d.V_cap=d.V_clean+sig/10^(c.snr/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticCase=tic;fit=run_inv_log_2_low_high_main(low,map,lib,cfg);elapsed=toc(ticCase);
    test=fit.second_frequency_test;rows(k).name=c.name;rows(k).snr_db=c.snr;
    rows(k).true_order=numel(c.f);rows(k).estimated_order=fit.model_order;
    rows(k).true_f1=c.f(1);
    if numel(c.f)>1,rows(k).true_f2=c.f(2);end
    rows(k).estimated_f1=fit.f_id(1);
    if numel(fit.f_id)>1,rows(k).estimated_f2=fit.f_id(2);end
    rows(k).frequency_error=frequency_error(fit.f_id,c.f);
    rows(k).detected_f2=test.second_frequency_hz;rows(k).delta_bic=test.delta_bic;
    rows(k).amplitude_z=test.amplitude_z;rows(k).strong_second=test.strong_second_frequency;
    rows(k).dual_skipped=fit.dual_search_skipped;rows(k).elapsed_s=elapsed;
    rows(k).used_samples=fit.used_sample_count;
    fprintf('%s: order %g->%g, err %.4g Hz, dBIC %.3g, z %.3g, %.2f s\n',...
        c.name,numel(c.f),fit.model_order,rows(k).frequency_error,...
        test.delta_bic,test.amplitude_z,elapsed);
end
Audit=struct2table(rows);out=fullfile(root,'output','model_order_boundary');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'model_order_boundary.csv'));save(fullfile(out,'model_order_boundary.mat'),'Audit');disp(Audit);
Result=struct('Audit',Audit,'outputDir',out);
end

function r=row0()
r=struct('name',"",'snr_db',NaN,'true_order',NaN,'estimated_order',NaN,...
 'true_f1',NaN,'true_f2',NaN,'estimated_f1',NaN,'estimated_f2',NaN,...
 'frequency_error',NaN,'detected_f2',NaN,'delta_bic',NaN,'amplitude_z',NaN,...
 'strong_second',false,'dual_skipped',false,'elapsed_s',NaN,'used_samples',NaN);
end

function e=frequency_error(fEst,fTrue)
if numel(fEst)~=numel(fTrue),e=NaN;return;end
e=max(abs(sort(fEst(:).')-sort(fTrue(:).')));
end
