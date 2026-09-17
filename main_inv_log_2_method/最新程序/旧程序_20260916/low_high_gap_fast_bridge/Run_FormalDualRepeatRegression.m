function Result = Run_FormalDualRepeatRegression(seedList)
%RUN_FORMALDUALREPEATREGRESSION Repeat shared-gap dual identification.
% The purpose is to test basin recovery across independent high-speed noise
% realizations without changing the production entry or its search budget.
if nargin < 1 || isempty(seedList), seedList = 1:3; end
root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2); cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=8; cfg.snrDb=20;
cfg.route30BaseSamples=900; cfg.route30ForwardModel="absolute";
cfg.route30EnableRhoFallback=false;
freqCases={[500 1300],[700 1200],[900 1400]};
ampCases=repmat([.25 .15],numel(freqCases),1);
phiCases=repmat([pi/4 -pi/3],numel(freqCases),1);
gapPairs=[.2 .8;.8 .2;.5 .8];
rows=repmat(row0(),0,1);
for is=1:numel(seedList)
    seed=seedList(is);
    for ic=1:numel(freqCases)
        fTrue=freqCases{ic}; gLow=gapPairs(ic,1); gHigh=gapPairs(ic,2);
        lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
            cfg,lib.domain,cfg.snrDb,880000+10000*seed+ic);
        lowAgg=aggregate_low_speed_template(lowData,cfg.alpha_k,cfg.R_tip,...
            lib.domain,lib.xGrid);
        lowState=estimate_highspeed_static_gap_raw(lowAgg.mapped,lib,cfg);
        high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
            cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
            ampCases(ic,:),fTrue,phiCases(ic,:),'noise_ratio');
        sig=rms(high.V_clean-min(high.V_clean));
        rng(990000+10000*seed+ic,'twister');
        high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
        hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        ticCase=tic; fit=run_low_high_gap_fast_vp_method(hm,lib,cfg,lowState); elapsed=toc(ticCase);
        fEst=sort(fit.f_id(:).');
        r=row0(); r.seed=seed; r.case_name=sprintf('%d+%d',fTrue);
        r.f1_true=fTrue(1); r.f2_true=fTrue(2); r.f1_est=fEst(1); r.f2_est=fEst(2);
        r.g_true=gHigh; r.g_est=fit.g_used; r.rmse=fit.rmse; r.elapsed_s=elapsed;
        r.max_frequency_error=max(abs(fEst-fTrue));
        r.success=fit.model_order==2 && r.max_frequency_error<=2 && abs(r.g_est-gHigh)<=.05;
        rows(end+1,1)=r; %#ok<AGROW>
        singleAnchor=NaN;
        if isfield(fit,'singleFit'), singleAnchor=fit.singleFit.f; end
        bankNearest=NaN;
        replayRank=NaN;
        if isfield(fit,'candidate_bank') && ~isempty(fit.candidate_bank)
            dd=arrayfun(@(q)max(abs(sort(q.f(:).')-fTrue)),fit.candidate_bank);
            bankNearest=min(dd);
        end
        if isfield(fit,'full_replay_frequencies')
            dd=max(abs(fit.full_replay_frequencies-fTrue),[],2);
            [~,truthIdx]=min(dd);
            replayRank=find(fit.full_replay_order==truthIdx,1,'first');
        end
        fprintf('seed %d case %s: single=%g f_est=%g,%g rmse=%g success=%d\n',...
            seed,r.case_name,singleAnchor,r.f1_est,r.f2_est,r.rmse,r.success);
        fprintf('  candidate nearest frequency error=%g Hz, replay rank=%g\n',...
            bankNearest,replayRank);
    end
end
Audit=struct2table(rows); out=fullfile(root,'output','formal_dual_repeat');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'formal_dual_repeat.csv')); save(fullfile(out,'formal_dual_repeat.mat'),'Audit');
Result=struct('Audit',Audit,'outputDir',out,'pass',all(Audit.success)); disp(Audit);
end

function r=row0()
r=struct('seed',NaN,'case_name',"",'f1_true',NaN,'f2_true',NaN,...
    'f1_est',NaN,'f2_est',NaN,'max_frequency_error',NaN,...
    'g_true',NaN,'g_est',NaN,'rmse',NaN,'elapsed_s',NaN,'success',false);
end
