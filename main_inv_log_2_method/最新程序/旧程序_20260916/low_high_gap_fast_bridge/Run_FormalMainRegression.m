function Result = Run_FormalMainRegression(useInformationSelection)
%RUN_FORMALMAINREGRESSION Unified-entry regression for the shared-gap main.
if nargin<1,useInformationSelection=false;end

root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir); tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN, ...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2); cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=8; cfg.snrDb=20;
cfg.route30BaseSamples=900; cfg.route30ForwardModel="absolute";
cfg.route30UseInformationSelection=logical(useInformationSelection);
cfg.route30EnableRhoFallback=false;
freqCases={[531],[733],[917],[500 1300],[700 1200],[900 1400]};
ampCases={[.25],[.25],[.25],[.25 .15],[.25 .15],[.25 .15]};
phiCases={[pi/4],[pi/4],[pi/4],[pi/4 -pi/3],[pi/4 -pi/3],[pi/4 -pi/3]};
gapPairs=[.2 .2; .8 .5; .5 .8; .8 .5; .2 .8; .8 .2];
rows=repmat(row0(),numel(freqCases),1);
for i=1:numel(freqCases)
    gLow=gapPairs(i,1); gHigh=gapPairs(i,2);
    lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z), ...
        cfg,lib.domain,cfg.snrDb,880000+i);
    high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z), ...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain, ...
        ampCases{i},freqCases{i},phiCases{i},'noise_ratio');
    sig=rms(high.V_clean-min(high.V_clean)); rng(881000+i,'twister');
    high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
    hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticCase=tic; fit=run_inv_log_2_low_high_main(lowData,hm,lib,cfg); elapsed=toc(ticCase);
    fEst=sort(fit.f_id(:).'); fTrue=sort(freqCases{i});
    fprintf('case %s: f_est=%g,%g rmse=%g\n',mat2str(fTrue),fEst(1),fEst(end),fit.rmse);
    rows(i).case_name=sprintf('%g',fTrue(1));
    if numel(fTrue)>1, rows(i).case_name=sprintf('%g+%g',fTrue); end
    rows(i).true_order=numel(fTrue); rows(i).estimated_order=fit.model_order;
    rows(i).f1_true=fTrue(1); rows(i).f1_est=fEst(1);
    if numel(fTrue)>1
        rows(i).f2_true=fTrue(2); rows(i).f2_est=fEst(2);
        [rows(i).candidate_nearest_error, rows(i).candidate_replay_rank] = ...
            audit_candidate_pool(fit, hm, lib, fTrue);
    end
    rows(i).g_true=gHigh; rows(i).g_est=fit.g_used;
    uTrue=zeros(size(hm.t_v));
    for j=1:numel(fTrue)
        uTrue=uTrue+ampCases{i}(j)*sin(2*pi*fTrue(j)*hm.t_v(:)+phiCases{i}(j));
    end
    VTrue=eval_gap_template(lib,gHigh,hm.x_v(:)-uTrue);
    rows(i).true_parameter_rmse=sqrt(mean((hm.V_a(:)-VTrue).^2));
    rows(i).max_frequency_error=max(abs(fEst-fTrue)); rows(i).rmse=fit.rmse;
    rows(i).elapsed_s=elapsed; rows(i).success=fit.model_order==numel(fTrue) && ...
        abs(fit.g_used-gHigh)<=.05 && rows(i).max_frequency_error<=2;
    rows(i).used_samples=fit.used_sample_count;
    rows(i).second_frequency_hz=fit.second_frequency_test.second_frequency_hz;
    rows(i).second_delta_bic=fit.second_frequency_test.delta_bic;
    rows(i).second_amplitude_z=fit.second_frequency_test.amplitude_z;
    rows(i).strong_second_frequency=fit.second_frequency_test.strong_second_frequency;
    if isfield(fit.information_selection,'information_min_eigenvalue')
        rows(i).information_min_eigenvalue=fit.information_selection.information_min_eigenvalue;
    end
end
Audit=struct2table(rows); out=fullfile(root,'output','formal_main_regression');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'formal_main_regression.csv')); save(fullfile(out,'formal_main_regression.mat'),'Audit');
disp(Audit);
Result=struct('Audit',Audit,'outputDir',out,'pass',all(Audit.success));
if ~Result.pass
    warning('invlog2:RegressionDiagnosticFailures', ...
        'Some cases remain local-frequency-competition diagnostics; see the CSV.');
end
end

function r=row0()
r=struct('case_name',"",'true_order',NaN,'estimated_order',NaN,...
    'f1_true',NaN,'f2_true',NaN,'f1_est',NaN,'f2_est',NaN,...
    'candidate_nearest_error',NaN,'candidate_replay_rank',NaN,...
    'true_parameter_rmse',NaN,'g_true',NaN,'g_est',NaN,...
    'max_frequency_error',NaN,'rmse',NaN,'elapsed_s',NaN,'used_samples',NaN,...
    'second_frequency_hz',NaN,'second_delta_bic',NaN,'second_amplitude_z',NaN,...
    'strong_second_frequency',false,...
    'information_min_eigenvalue',NaN,'success',false);
end

function [nearestError, replayRank] = audit_candidate_pool(fit, map, lib, fTrue)
nearestError=NaN; replayRank=NaN;
if ~isfield(fit,'candidate_bank') || isempty(fit.candidate_bank), return; end
bank=fit.candidate_bank; n=numel(bank); distance=inf(n,1); score=inf(n,1);
t=map.t_v(:); x=map.x_v(:); V=map.V_a(:);
for k=1:n
    distance(k)=max(abs(sort(bank(k).f(:).')-fTrue));
    V0=eval_gap_template(lib,bank(k).g,x-bank(k).dx-fit_u(bank(k).p,t));
    r=V-V0; good=isfinite(r);
    if any(good), score(k)=mean(r(good).^2); end
end
[nearestError,nearestIdx]=min(distance);
[~,ord]=sort(score,'ascend'); replayRank=find(ord==nearestIdx,1,'first');
end
