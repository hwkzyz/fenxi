function Result=Run_FastDualCandidateBenchmark(methods,snrList,caseIndices)
%RUN_FASTDUALCANDIDATEBENCHMARK Fixed-seed fast/full dual-frequency audit.
% Truth enters only below in simulation and post-fit metrics.  Candidate
% generation receives map, library, configuration, and calibrated gap only.
if nargin<1||isempty(methods),methods=["fast","full"];end
methods=string(methods);
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=8;
cfg0.route30ForwardModel="absolute";cfg0.route30UseAllTrustedSamples=true;
cfg0.route30EnableRhoFallback=true;cfg0.route30FrequencyStepHz=5;
if nargin<2||isempty(snrList),snrList=[25 15 5];end
seed=20260809;
cases=[struct('name',"grid_700_1200_boundary",'f',[700 1200],'gLow',.8,'gHigh',.2),...
    struct('name',"offgrid_733_1217",'f',[733 1217],'gLow',.2,'gHigh',.8),...
    struct('name',"offgrid_917p5_1382p5",'f',[917.5 1382.5],'gLow',.5,'gHigh',.8)];
if nargin>=3&&~isempty(caseIndices),cases=cases(caseIndices);end
Atrue=[.25 .15];phiTrue=[pi/4 -pi/3];
rows=repmat(row0(),0,1);
for isnr=1:numel(snrList)
    for ic=1:numel(cases)
        c=cases(ic);cfg=cfg0;cfg.snrDb=snrList(isnr);
        low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),...
            cfg,lib.domain,cfg.snrDb,seed+1000*isnr+10*ic);
        agg=aggregate_low_speed_template(low,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);
        lowState=estimate_highspeed_static_gap_raw(agg.mapped,lib,cfg);
        high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,c.gHigh,z),...
            cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
            Atrue,c.f,phiTrue,'noise_ratio');
        sig=rms(high.V_clean-min(high.V_clean));
        rng(seed+100000+1000*isnr+10*ic,'twister');
        high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
        hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        for im=1:numel(methods)
            q=cfg;q.route30UseFastDualCandidates=methods(im)=="fast";
            tc=tic;fit=run_low_high_gap_fast_vp_method(hm,lib,q,lowState);wall=toc(tc);
            r=row0();r.case_name=c.name;r.method=methods(im);r.snr_db=cfg.snrDb;
            r.seed=seed;r.f1_true=c.f(1);r.f2_true=c.f(2);r.g_low_true=c.gLow;
            r.g_high_true=c.gHigh;r.g_low_hat=lowState.gHat;r.g_est=fit.g_used;
            [fEst,ord]=sort(fit.f_id(:).');r.f1_est=fEst(1);r.f2_est=fEst(2);
            aEst=fit.A_id(ord);pEst=fit.phi_id(ord);
            r.a1_est=aEst(1);r.a2_est=aEst(2);
            r.frequency_error_hz=max(abs(fEst-c.f));
            r.amplitude_error_mm=max(abs(aEst-Atrue));
            r.phase_error_rad=max(abs(angle(exp(1i*(pEst-phiTrue)))));
            r.gap_error_mm=fit.g_used-c.gHigh;r.rmse_V=fit.rmse;r.elapsed_s=wall;
            r.model_order=fit.model_order;r.success=fit.model_order==2&&...
                r.frequency_error_hz<=2&&r.amplitude_error_mm<=.06&&abs(r.gap_error_mm)<=.05;
            r.fallback_used=getv(fit,'fast_candidate_fallback',false);
            fc=getv(fit,'fast_candidate',struct());
            r.generator_confident=getv(fc,'confident',false);
            r.generator_pair_count=getv(fc,'candidate_pair_count',0);
            r.full_pair_count=getv(fc,'full_unordered_pair_count',28203);
            cp=getv(fc,'candidate_pairs_hz',zeros(0,2));
            gp=getv(fc,'generated_candidate_pairs_hz',cp);
            r.generated_candidate_pair_count=size(gp,1);
            [r.generated_nearest_error_hz,r.generated_nearest_rank]=pair_error_rank(gp,c.f);
            r.nearest_candidate_error_hz=pair_error(cp,c.f);
            r.candidate_recalled=r.nearest_candidate_error_hz<=20;
            r.final_bank_recalled=bank_recalled(getv(fit,'candidate_bank',[]),c.f,5);
            r.profile_pair_evaluations=getv(fit,'frequency_pair_evaluations',NaN);
            pt=getv(fit,'pipeline_timing',struct());
            r.single_order_s=getv(pt,'single_order_s',NaN);
            r.second_test_s=getv(pt,'second_frequency_test_s',NaN);
            r.candidate_generation_s=getv(pt,'candidate_generation_s',0);
            r.dual_gap_profile_s=getv(pt,'dual_gap_profile_s',NaN);
            r.candidate_augmentation_s=getv(pt,'candidate_augmentation_s',NaN);
            r.full_replay_refine_s=getv(pt,'full_replay_refine_s',NaN);
            r.vp_precompute_s=getv(fit,'time_precomp_s',NaN);
            r.vp_screen_s=getv(fit,'time_screen_s',NaN);
            r.wave_refine_s=getv(fit,'time_refine_wave_s',NaN);
            r.joint_refine_s=getv(fit,'time_joint_s',NaN);
            r.boundary_case=abs(c.gHigh-min(lib.gapTrain))<1e-12;
            rows(end+1,1)=r; %#ok<AGROW>
            fprintf('%s %g dB %s: err=%.3g Hz gap=%.3g mm amp=%.3g mm t=%.2f s recall=%d fallback=%d\n',...
                methods(im),cfg.snrDb,c.name,r.frequency_error_hz,r.gap_error_mm,...
                r.amplitude_error_mm,wall,r.candidate_recalled,r.fallback_used);
        end
    end
end
Audit=struct2table(rows);Summary=groupsummary(Audit,'method',...
    {'mean','sum'},{'success','candidate_recalled','elapsed_s','fallback_used'});
tag=strjoin(methods,'_')+sprintf('_%dcases_%dsnrs',numel(cases),numel(snrList));
out=fullfile(root,'output',char("fast_dual_candidate_benchmark_"+tag));
if~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'case_metrics.csv'));
writetable(Summary,fullfile(out,'method_summary.csv'));
save(fullfile(out,'benchmark.mat'),'Audit','Summary','cases','snrList','seed','-v7.3');
Result=struct('Audit',Audit,'Summary',Summary,'outputDir',out);disp(Summary);
end

function tf=pair_recalled(pairs,truth,tol)
if isempty(pairs),tf=false;return;end
tf=any(max(abs(sort(pairs,2)-sort(truth)),[],2)<=tol);
end
function e=pair_error(pairs,truth)
if isempty(pairs),e=Inf;return;end
e=min(max(abs(sort(pairs,2)-sort(truth)),[],2));
end
function [e,k]=pair_error_rank(pairs,truth)
if isempty(pairs),e=Inf;k=NaN;return;end
[e,k]=min(max(abs(sort(pairs,2)-sort(truth)),[],2));
end
function tf=bank_recalled(bank,truth,tol)
if isempty(bank),tf=false;return;end
p=vertcat(bank.f);tf=pair_recalled(p,truth,tol);
end
function v=getv(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
function r=row0()
r=struct('case_name',"",'method',"",'snr_db',NaN,'seed',NaN,...
    'f1_true',NaN,'f2_true',NaN,'f1_est',NaN,'f2_est',NaN,...
    'g_low_true',NaN,'g_high_true',NaN,'g_low_hat',NaN,'g_est',NaN,...
    'a1_est',NaN,'a2_est',NaN,'frequency_error_hz',NaN,...
    'amplitude_error_mm',NaN,'phase_error_rad',NaN,'gap_error_mm',NaN,...
    'rmse_V',NaN,'elapsed_s',NaN,'model_order',NaN,'success',false,...
    'fallback_used',false,'generator_confident',false,'generator_pair_count',NaN,...
    'full_pair_count',NaN,'nearest_candidate_error_hz',NaN,...
    'generated_candidate_pair_count',NaN,'generated_nearest_error_hz',NaN,...
    'generated_nearest_rank',NaN,...
    'candidate_recalled',false,'final_bank_recalled',false,...
    'profile_pair_evaluations',NaN,'single_order_s',NaN,'second_test_s',NaN,...
    'candidate_generation_s',NaN,'dual_gap_profile_s',NaN,...
    'candidate_augmentation_s',NaN,'full_replay_refine_s',NaN,'boundary_case',false);
r.vp_precompute_s=NaN;r.vp_screen_s=NaN;r.wave_refine_s=NaN;r.joint_refine_s=NaN;
end
