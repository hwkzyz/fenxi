function Result=Run_FastDualGeneratorAudit()
%RUN_FASTDUALGENERATORAUDIT Candidate-only audit; truth is evaluation-only.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=8;
cfg0.f1Grid=300:5:1500;cfg0.f2Grid=cfg0.f1Grid;
cfg0.route30SingleGrid=5:5:1500;
cfg0.mainVpGapHalfWidth=.70;cfg0.route30FastCandidateGapAnchors=5;
cfg0.route30FastCandidatePairsPerGap=16;cfg0.route30FastCandidateMaxPairs=640;
cfg0.route30FastCandidateSingleAnchorCount=8;
snrList=[25 15 5];seed=20260809;
cases=[struct('name',"grid_700_1200_boundary",'f',[700 1200],'gLow',.8,'gHigh',.2),...
    struct('name',"offgrid_733_1217",'f',[733 1217],'gLow',.2,'gHigh',.8),...
    struct('name',"offgrid_917p5_1382p5",'f',[917.5 1382.5],'gLow',.5,'gHigh',.8)];
A=[.25 .15];phi=[pi/4 -pi/3];rows=repmat(row0(),0,1);
for isnr=1:numel(snrList)
    for ic=1:numel(cases)
        c=cases(ic);cfg=cfg0;cfg.snrDb=snrList(isnr);
        low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),...
            cfg,lib.domain,cfg.snrDb,seed+1000*isnr+10*ic);
        agg=aggregate_low_speed_template(low,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);
        st=estimate_highspeed_static_gap_raw(agg.mapped,lib,cfg);
        high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,c.gHigh,z),...
            cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
            A,c.f,phi,'noise_ratio');sig=rms(high.V_clean-min(high.V_clean));
        rng(seed+100000+1000*isnr+10*ic,'twister');
        high.V_cap=high.V_clean+sig/10^(cfg.snrDb/20)*randn(size(high.V_clean));
        hm=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        single=fit_single_vp_main(hm,lib,cfg,st);
        anchors=[single.f_id,[single.candidates.f_id]];
        [pairs,info]=generate_fast_dual_frequency_pairs(hm,lib,cfg,st,anchors);
        e=max(abs(sort(pairs,2)-sort(c.f)),[],2);r=row0();r.case_name=c.name;
        r.snr_db=cfg.snrDb;r.seed=seed;r.g_low_hat=st.gHat;r.g_high_true=c.gHigh;
        r.candidate_count=size(pairs,1);r.full_pair_count=info.full_unordered_pair_count;
        r.nearest_pair_error_hz=min(e);r.grid_recalled=r.nearest_pair_error_hz<=2.5;
        r.refinement_basin_recalled=r.nearest_pair_error_hz<=20;
        r.confident=info.confident;r.fallback=info.fallback_to_full;r.elapsed_s=info.elapsed_s;
        r.boundary_case=abs(c.gHigh-min(lib.gapTrain))<1e-12;rows(end+1,1)=r; %#ok<AGROW>
    end
end
Audit=struct2table(rows);out=fullfile(root,'output','fast_dual_generator_audit');
if~exist(out,'dir'),mkdir(out);end;writetable(Audit,fullfile(out,'generator_metrics.csv'));
save(fullfile(out,'generator_metrics.mat'),'Audit','-v7.3');Result=struct('Audit',Audit,'outputDir',out);disp(Audit);
end
function r=row0()
r=struct('case_name',"",'snr_db',NaN,'seed',NaN,'g_low_hat',NaN,...
    'g_high_true',NaN,'candidate_count',NaN,'full_pair_count',NaN,...
    'nearest_pair_error_hz',NaN,'grid_recalled',false,...
    'refinement_basin_recalled',false,'confident',false,'fallback',false,...
    'elapsed_s',NaN,'boundary_case',false);
end
