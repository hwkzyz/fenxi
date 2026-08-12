function Result=Run_StructuredRandomizedRegression(casesPerMode,seedBase)
%RUN_STRUCTUREDRANDOMIZEDREGRESSION Randomized four-mode acceptance for robust5.
if nargin<1,casesPerMode=2;end
if nargin<2,seedBase=230000;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);
cfg=configure_structured_identification(base.cfgCase,"experiment_current");
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=8;cfg.snrDb=5;
fRot=cfg.RPM_high/60;rng(seedBase,'twister');cases=make_cases(casesPerMode,fRot);
rows=repmat(row0(),numel(cases),1);
for i=1:numel(cases)
    c=cases(i);
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),...
        cfg,lib.domain,5,seedBase+10*i);
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,c.gHigh,z),...
        cfg.RPM_high,8,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,c.A,c.f,c.phi,'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(seedBase+10*i+1,'twister');
    d.V_cap=d.V_clean+sig/10^(5/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticFit=tic;fit=run_inv_log_2_structured_main(low,map,lib,cfg,c.mode);elapsed=toc(ticFit);
    fTrue=sort(c.f(:).');fEst=sort(fit.f_id(:).');ferr=max(abs(fTrue-fEst));
    parameterSuccess=ferr<=1&&abs(fit.g_used-c.gHigh)<=.05;
    falseOfficial=fit.identification_acceptable&&~parameterSuccess;
    rows(i)=struct('case_id',c.id,'mode',c.mode,'frequency_true',string(mat2str(fTrue,7)),...
        'frequency_est',string(mat2str(fEst,7)),'frequency_error',ferr,...
        'amplitude_true',string(mat2str(c.A,7)),'amplitude_est',string(mat2str(fit.A_id,7)),...
        'g_low',c.gLow,'g_true',c.gHigh,'g_est',fit.g_used,'rmse',fit.rmse,...
        'identification_status',fit.identification_status,...
        'identification_confident',fit.identification_confident,...
        'observation_action',fit.observation_action,'parameter_success',parameterSuccess,...
        'false_official_identification',falseOfficial,'used_samples',fit.used_sample_count,...
        'elapsed_s',elapsed);
end
Audit=struct2table(rows);out=fullfile(root,'output',sprintf('structured_randomized_n%d_seed%d',casesPerMode,seedBase));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'structured_randomized.csv'));
save(fullfile(out,'structured_randomized.mat'),'Audit','cases','seedBase','casesPerMode');disp(Audit);
if any(Audit.false_official_identification),error('A wrong randomized result passed the official reporting gate.');end
Result=struct('Audit',Audit,'successCount',nnz(Audit.parameter_success),...
    'officialCount',nnz(Audit.identification_confident),'outputDir',out);
end

function C=make_cases(n,fRot)
modes=["single_sync","single_async","dual_sync_sync","dual_sync_async"];
gList=[.2 .5 .8 1.1];C=repmat(case0(),numel(modes)*n,1);at=0;
for im=1:numel(modes)
    for j=1:n
        at=at+1;mode=modes(im);gLow=gList(randi(numel(gList)));gHigh=gList(randi(numel(gList)));
        switch mode
            case "single_sync"
                eo=randi([6 30]);f=eo*fRot;A=.1+.2*rand;phi=2*pi*rand;
            case "single_async"
                f=off_order_frequency(fRot,[]);A=.1+.2*rand;phi=2*pi*rand;
            case "dual_sync_sync"
                eo=sort(randperm(25,2)+5);f=eo*fRot;A=.1+.2*rand(1,2);phi=2*pi*rand(1,2);
            otherwise
                eo=randi([6 30]);fs=eo*fRot;fa=off_order_frequency(fRot,fs);
                f=[fs fa];A=.1+.2*rand(1,2);phi=2*pi*rand(1,2);
        end
        C(at)=struct('id',sprintf('%s_%02d',mode,j),'mode',mode,'gLow',gLow,...
            'gHigh',gHigh,'A',A,'f',f,'phi',phi);
    end
end
end

function f=off_order_frequency(fRot,avoid)
while true
    f=300+1200*rand;
    if abs(f/fRot-round(f/fRot))*fRot>=5&&...
            (isempty(avoid)||min(abs(f-avoid))>=20),return;end
end
end

function c=case0()
c=struct('id',"",'mode',"",'gLow',NaN,'gHigh',NaN,'A',NaN,'f',NaN,'phi',NaN);
end

function r=row0()
r=struct('case_id',"",'mode',"",'frequency_true',"",'frequency_est',"",...
    'frequency_error',NaN,'amplitude_true',"",'amplitude_est',"",...
    'g_low',NaN,'g_true',NaN,'g_est',NaN,'rmse',NaN,'identification_status',"",...
    'identification_confident',false,'observation_action',"",...
    'parameter_success',false,'false_official_identification',false,...
    'used_samples',NaN,'elapsed_s',NaN);
end
