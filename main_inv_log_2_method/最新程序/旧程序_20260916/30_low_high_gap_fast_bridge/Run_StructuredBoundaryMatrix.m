function Result=Run_StructuredBoundaryMatrix(seedOffsets,stressMode,caseIds,numRevs,snrDb,layoutMode)
%RUN_STRUCTUREDBOUNDARYMATRIX Stress the four solvers at range and separation boundaries.
if nargin<1,seedOffsets=0;end
if nargin<2,stressMode=false;end
if nargin<3,caseIds=[];end
if nargin<4||isempty(numRevs),numRevs=8;end
if nargin<5||isempty(snrDb),snrDb=5;end
if nargin<6||isempty(layoutMode),layoutMode="current";end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=configure_structured_identification(base.cfgCase,"experiment_current");
switch string(layoutMode)
    case "uniform3",cfg.alpha_k=(0:2)*2*pi/3;
    case "uniform4",cfg.alpha_k=(0:3)*2*pi/4;
    case "uniform5",cfg.alpha_k=(0:4)*2*pi/5;
    case "golden3"
        golden=(1+sqrt(5))/2;cfg.alpha_k=sort(mod((0:2)*2*pi/golden,2*pi));
    case "golden5"
        golden=(1+sqrt(5))/2;cfg.alpha_k=sort(mod((0:4)*2*pi/golden,2*pi));
    case "optimized5"
        cfg.alpha_k=[0 .973466 2.42106 2.90560 4.06032];
end
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=numRevs;cfg.snrDb=snrDb;
cfg.structuredMinRevolutions=8;cfg.route30SingleFrequencyRangeHz=[300 1500];
cfg.route30DualFrequencyRangeHz=[300 1500];fRot=cfg.RPM_high/60;
if stressMode==1
    cfg.dualSyncKeepPerGap=100;cfg.dualSyncIteratedReplayCount=2100;cfg.dualSyncRefineCount=50;
    cfg.syncAsyncKeepPerGap=100;cfg.syncAsyncIteratedReplayCount=2100;cfg.syncAsyncRefineCount=50;
elseif stressMode==2
    cfg.dualSyncKeepPerGap=300;cfg.dualSyncIteratedReplayCount=6000;cfg.dualSyncRefineCount=100;
    cfg.syncAsyncKeepPerGap=300;cfg.syncAsyncIteratedReplayCount=6000;cfg.syncAsyncRefineCount=100;
end
allCases=build_cases(fRot);for j=1:numel(allCases),allCases(j).seedIndex=j;end
cases=allCases;
if ~isempty(caseIds),cases=cases(ismember([cases.id],string(caseIds)));end
rows=repmat(row0(),numel(cases)*numel(seedOffsets),1);at=0;
for iseed=1:numel(seedOffsets)
    seedOffset=seedOffsets(iseed);
    for i=1:numel(cases)
        at=at+1;c=cases(i);nComp=numel(c.f);phi=c.phi(1:nComp);
        low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),...
            cfg,lib.domain,snrDb,220000+100*c.seedIndex+seedOffset);
        d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,c.gHigh,z),...
            cfg.RPM_high,numRevs,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,c.A,c.f,phi,'noise_ratio');
        sig=rms(d.V_clean-min(d.V_clean));rng(220050+100*c.seedIndex+seedOffset,'twister');
        d.V_cap=d.V_clean+sig/10^(snrDb/20)*randn(size(d.V_clean));
        map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        ticCase=tic;fit=run_inv_log_2_structured_main(low,map,lib,cfg,c.mode);elapsed=toc(ticCase);
        syncOnlyEo=NaN;syncAverageEo=NaN;
        if c.mode=="dual_sync_async"
            lowMap=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);
            state=estimate_highspeed_static_gap_raw(lowMap,lib,cfg);
            qSync=run_single_sync_voltage_vp(map,lib,cfg,state);syncOnlyEo=qSync.eo_id;
            avgMap=average_map_by_revolution(map);
            qAvg=run_single_sync_voltage_vp(avgMap,lib,cfg,state);syncAverageEo=qAvg.eo_id;
        end
        fEst=sort(fit.f_id(:).');fTrue=sort(c.f(:).');ferr=max(abs(fEst-fTrue));
        modeOk=mode_success(c,fit,fRot);gapOk=abs(fit.g_used-c.gHigh)<=.05;
        rows(at)=struct('case_id',c.id,'mode',c.mode,'seed_offset',seedOffset,...
            'frequency_true',string(mat2str(fTrue,7)),'frequency_est',string(mat2str(fEst,7)),...
            'max_frequency_error',ferr,'amplitude_min',min(c.A),'g_true',c.gHigh,...
            'amplitude_est',string(mat2str(fit.A_id,7)),'amplitude_est_min',min(fit.A_id),...
            'sync_only_eo',syncOnlyEo,...
            'sync_average_eo',syncAverageEo,...
            'g_est',fit.g_used,'rmse',fit.rmse,'frequency_margin',fit.frequency_margin,...
            'noise_normalized_margin',fit.noise_normalized_margin,...
            'identification_confident',fit.identification_confident,...
            'identification_status',fit.identification_status,...
            'requires_more_observations',fit.requires_more_observations,...
            'observation_action',fit.observation_action,...
            'competitor',string(mat2str(fit.competing_frequency_or_order,7)),...
            'used_samples',fit.used_sample_count,'elapsed_s',elapsed,...
            'success',modeOk&&gapOk&&ferr<=1);
    end
end
Audit=struct2table(rows);tag="structured_boundary_matrix";
if numel(cases)~=numel(allCases)||stressMode~=0||numRevs~=8||snrDb~=5||string(layoutMode)~="current"
    tag=tag+sprintf('_%dcases_%drev_%gdb_stress%d_%s',numel(cases),numRevs,snrDb,stressMode,layoutMode);
end
out=fullfile(root,'output',char(tag));
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'structured_boundary_matrix.csv'));
save(fullfile(out,'structured_boundary_matrix.mat'),'Audit','cases','seedOffsets');disp(Audit);
if any(~Audit.success&Audit.identification_confident)
    error('A wrong boundary solution passed the confidence gate.');
end
Result=struct('Audit',Audit,'outputDir',out,'successCount',nnz(Audit.success));
end

function cases=build_cases(fRot)
phi=[.13 2.21];
cases=[...
 c("ss_low", "single_sync",.2,.8,.10,6*fRot,6,phi),...
 c("ss_mid", "single_sync",.8,.5,.25,17*fRot,17,phi),...
 c("ss_high","single_sync",.5,.2,.10,30*fRot,30,phi),...
 c("sa_low", "single_async",.2,.8,.10,303,NaN,phi),...
 c("sa_mid", "single_async",.8,.5,.25,947,NaN,phi),...
 c("sa_high","single_async",.5,.2,.10,1497,NaN,phi),...
 c("dss_adj_low", "dual_sync_sync",.2,.8,[.10 .10],[6 7]*fRot,[6 7],phi),...
 c("dss_wide",    "dual_sync_sync",.8,.5,[.25 .10],[10 30]*fRot,[10 30],phi),...
 c("dss_adj_high","dual_sync_sync",.5,.2,[.10 .25],[29 30]*fRot,[29 30],phi),...
 c("dsa_close",  "dual_sync_async",.2,.8,[.25 .15],[10*fRot 560],10,phi),...
 c("dsa_wide",   "dual_sync_async",.8,.5,[.25 .10],[6*fRot 1493],6,phi),...
 c("dsa_reverse","dual_sync_async",.5,.2,[.10 .25],[29*fRot 317],29,phi)];
end

function q=c(id,mode,gLow,gHigh,A,f,eo,phi)
q=struct('id',string(id),'mode',string(mode),'gLow',gLow,'gHigh',gHigh,...
    'A',A,'f',f,'eo',eo,'phi',phi);
end

function ok=mode_success(c,fit,fRot)
switch c.mode
    case "single_sync",ok=fit.eo_id==c.eo;
    case "single_async",ok=abs(fit.f_id-c.f)<=1;
    case "dual_sync_sync",ok=isequal(fit.eo_id,c.eo);
    otherwise
        ok=fit.sync_eo==c.eo&&abs(fit.async_frequency_hz-c.f(2))<=1&&...
            abs(fit.f_id(1)-c.eo*fRot)<=eps(c.eo*fRot);
end
end

function r=row0()
r=struct('case_id',"",'mode',"",'seed_offset',NaN,'frequency_true',"",...
    'frequency_est',"",'max_frequency_error',NaN,'amplitude_min',NaN,...
    'amplitude_est',"",'amplitude_est_min',NaN,...
    'sync_only_eo',NaN,...
    'sync_average_eo',NaN,...
    'g_true',NaN,'g_est',NaN,'rmse',NaN,'frequency_margin',NaN,...
    'noise_normalized_margin',NaN,'identification_confident',false,...
    'identification_status',"",'requires_more_observations',false,...
    'observation_action',"",...
    'competitor',"",'used_samples',NaN,...
    'elapsed_s',NaN,'success',false);
end

function out=average_map_by_revolution(in)
% Coherent average by sensor and within-passage index; asynchronous content cancels.
out=in;sensors=unique(in.S_v(:)).';revs=unique(in.rev_v(:)).';
groups=cell(numel(sensors),numel(revs));nMin=inf;
for is=1:numel(sensors)
    for ir=1:numel(revs)
        q=find(in.S_v==sensors(is)&in.rev_v==revs(ir));groups{is,ir}=q;nMin=min(nMin,numel(q));
    end
end
if ~isfinite(nMin)||nMin<8,out=in;return;end
t=zeros(nMin*numel(sensors),1);x=t;V=t;theta=t;S=t;
for is=1:numel(sensors)
    Q=zeros(nMin,numel(revs));T=Q;X=Q;Th=Q;
    for ir=1:numel(revs)
        q=groups{is,ir}(1:nMin);Q(:,ir)=in.V_a(q);T(:,ir)=in.t_v(q);
        X(:,ir)=in.x_v(q);Th(:,ir)=in.theta_v(q)-2*pi*(revs(ir)-1);
    end
    at=((is-1)*nMin+(1:nMin)).';V(at)=mean(Q,2);t(at)=mean(T,2);x(at)=mean(X,2);
    theta(at)=mean(Th,2);S(at)=sensors(is);
end
out.t_v=t;out.V_a=V;out.x_v=x;out.theta_v=theta;out.S_v=S;
out.rev_v=ones(size(t));out.selectionMode="coherent_revolution_average";
end
