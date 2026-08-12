function R=run_inv_log_2_unknown_structure(highMap,templateLib,cfg,staticState)
%RUN_INV_LOG_2_UNKNOWN_STRUCTURE Hierarchical unknown-structure router.
% The router first compares the two single-component models.  Dual models
% are evaluated only when the best single model cannot explain the record at
% the supplied voltage-noise level.  Expensive dual budgets are used only
% for a low-quality or close fast comparison.
if nargin<4,error('highMap, templateLib, cfg, and staticState are required.');end
noiseStd=get_field(cfg,'unknownStructureNoiseStdV',NaN);
noiseFactor=get_field(cfg,'unknownStructureNoiseRmseFactor',1.5);
bicTie=get_field(cfg,'unknownStructureBicTieThreshold',10);

singleModes=["single_sync","single_async"];
single=cell(1,2);singleBic=nan(1,2);singleTime=nan(1,2);
for i=1:2
    q=cfg;q.route30StructuredFastBudgets=true;
    tc=tic;single{i}=run_inv_log_2_structured_main(staticState,highMap,templateLib,q,singleModes(i));
    singleTime(i)=toc(tc);singleBic(i)=model_bic(single{i},numel(highMap.V_a));
end
[~,iSingle]=min(singleBic);bestSingle=single{iSingle};
singleAdequate=residual_adequate(bestSingle.rmse,noiseStd,noiseFactor)&&...
    basic_quality(bestSingle,staticState,cfg);

if singleAdequate
    R=bestSingle;
    R.structure_selected=singleModes(iSingle);
    R.structure_stage="single";R.structure_budget="fast";
    R.structure_candidates=pack_candidates(singleModes,single,singleBic,singleTime);
    R.structure_bic_margin=other_margin(singleBic,iSingle);
    R.structure_router_reason="single_model_noise_compatible";
    R.structure_router_elapsed_s=sum(singleTime);
    return
end

dualModes=["dual_sync_sync","dual_sync_async"];
dual=cell(1,2);dualBic=nan(1,2);dualTime=nan(1,2);
for i=1:2
    q=fast_dual_cfg(cfg,dualModes(i));
    tc=tic;dual{i}=run_inv_log_2_structured_main(staticState,highMap,templateLib,q,dualModes(i));
    dualTime(i)=toc(tc);dualBic(i)=model_bic(dual{i},numel(highMap.V_a));
end
% The coarse global EO scan is an inexpensive screening stage, but a low
% residual alone cannot reveal a wrong combination-order solution.  Always
% certify the synchronous two-component candidate with the Top-M x
% multi-start search before comparing structures.
q=standard_dual_cfg(cfg,"dual_sync_sync");
tc=tic;candidate=run_inv_log_2_structured_main(staticState,highMap,templateLib,q,"dual_sync_sync");
dualTime(1)=dualTime(1)+toc(tc);candidateBic=model_bic(candidate,numel(highMap.V_a));
if candidateBic<=dualBic(1),dual{1}=candidate;dualBic(1)=candidateBic;end
[~,iDual]=min(dualBic);dualMargin=other_margin(dualBic,iDual);
fastQuality=residual_adequate(dual{iDual}.rmse,noiseStd,noiseFactor)&&...
    basic_quality(dual{iDual},staticState,cfg);

rerun=find(dualBic-min(dualBic)<=bicTie);
if ~fastQuality||dualMargin<=bicTie
    if isempty(rerun),rerun=iDual;end
    for j=rerun
        q=standard_dual_cfg(cfg,dualModes(j));
        tc=tic;candidate=run_inv_log_2_structured_main(staticState,highMap,templateLib,q,dualModes(j));
        dualTime(j)=dualTime(j)+toc(tc);candidateBic=model_bic(candidate,numel(highMap.V_a));
        if candidateBic<=dualBic(j),dual{j}=candidate;dualBic(j)=candidateBic;end
    end
    [~,iDual]=min(dualBic);dualMargin=other_margin(dualBic,iDual);
    budget="standard";
else
    budget="certified_sync";
end

allModes=[singleModes dualModes];allFits=[single dual];allBic=[singleBic dualBic];allTime=[singleTime dualTime];
[~,iBest]=min(allBic);R=allFits{iBest};
R.structure_selected=allModes(iBest);
if R.model_order==1
    R.structure_stage="single";R.structure_budget="fast";
else
    R.structure_stage="dual";R.structure_budget=budget;
end
R.structure_candidates=pack_candidates(allModes,allFits,allBic,allTime);
R.structure_bic_margin=other_margin(allBic,iBest);
R.structure_router_reason="single_inadequate_compare_dual";
R.structure_router_elapsed_s=sum(allTime);
R.structure_noise_std_V=noiseStd;
R.structure_noise_rmse_factor=noiseFactor;
R.structure_dual_fast_bic_margin=dualMargin;
end

function q=fast_dual_cfg(cfg,mode)
q=cfg;q.route30StructuredFastBudgets=false;
if mode=="dual_sync_sync"
    q.dualSyncSearchStrategy="global_budget";
    q.dualSyncGapCount=get_field(cfg,'unknownFastDualSyncGapCount',15);
    q.dualSyncKeepPerGap=get_field(cfg,'unknownFastDualSyncKeepPerGap',40);
    q.dualSyncIteratedReplayCount=get_field(cfg,'unknownFastDualSyncReplayCount',400);
    q.dualSyncRefineCount=get_field(cfg,'unknownFastDualSyncRefineCount',12);
else
    q.syncAsyncGapCount=get_field(cfg,'unknownFastSyncAsyncGapCount',11);
    q.syncAsyncFrequencyStepHz=get_field(cfg,'unknownFastSyncAsyncStepHz',10);
    q.syncAsyncKeepPerGap=get_field(cfg,'unknownFastSyncAsyncKeepPerGap',40);
    q.syncAsyncIteratedReplayCount=get_field(cfg,'unknownFastSyncAsyncReplayCount',400);
    q.syncAsyncRefineCount=get_field(cfg,'unknownFastSyncAsyncRefineCount',12);
end
end

function q=standard_dual_cfg(cfg,mode)
q=cfg;q.route30StructuredFastBudgets=false;
if mode=="dual_sync_sync"
    q.dualSyncSearchStrategy="coarse_top_m_multistart";
    q.dualSyncGapCount=get_field(cfg,'unknownStandardDualSyncGapCount',11);
    q.dualSyncKeepPerGap=get_field(cfg,'unknownStandardDualSyncKeepPerGap',80);
    q.dualSyncIteratedReplayCount=get_field(cfg,'unknownStandardDualSyncReplayCount',2100);
    q.dualSyncRefineCount=get_field(cfg,'unknownStandardDualSyncRefineCount',50);
    q.dualSyncTopPairCount=get_field(cfg,'dualSyncTopPairCount',40);
    q.dualSyncStartsPerTopPair=get_field(cfg,'dualSyncStartsPerTopPair',2);
else
    q.syncAsyncGapCount=get_field(cfg,'unknownStandardSyncAsyncGapCount',15);
    q.syncAsyncFrequencyStepHz=get_field(cfg,'unknownStandardSyncAsyncStepHz',10);
    q.syncAsyncKeepPerGap=get_field(cfg,'unknownStandardSyncAsyncKeepPerGap',100);
    q.syncAsyncIteratedReplayCount=get_field(cfg,'unknownStandardSyncAsyncReplayCount',2100);
    q.syncAsyncRefineCount=get_field(cfg,'unknownStandardSyncAsyncRefineCount',50);
end
end

function b=model_bic(fit,n)
k=parameter_count(string(fit.mode));rss=max(n*fit.rmse^2,realmin);
b=n*log(rss/n)+k*log(n);
end

function k=parameter_count(mode)
switch mode
    case "single_sync",k=4;
    case "single_async",k=5;
    case "dual_sync_sync",k=6;
    case "dual_sync_async",k=7;
    otherwise,k=NaN;
end
end

function ok=residual_adequate(rmse,sigma,factor)
ok=isfinite(sigma)&&sigma>0&&isfinite(rmse)&&rmse<=factor*sigma;
end

function ok=basic_quality(fit,state,cfg)
half=get_field(cfg,'route30GapHalfWidthMm',.70);
gapBoundary=abs(fit.g_used-state.gHat)>=.97*half;
upper=get_field(cfg,'route30AmplitudeUpperMm',1.5);
amplitudeBoundary=any(~isfinite(fit.A_id))||any(fit.A_id>=.98*upper);
ok=isfinite(fit.rmse)&&~gapBoundary&&~amplitudeBoundary;
end

function C=pack_candidates(modes,fits,bic,time)
C=repmat(struct('mode',"",'model_order',NaN,'rmse_V',NaN,'bic',NaN,...
    'elapsed_s',NaN,'g_used_mm',NaN,'frequencies_hz',[],'amplitudes_mm',[]),numel(modes),1);
for i=1:numel(modes)
    C(i)=struct('mode',modes(i),'model_order',fits{i}.model_order,'rmse_V',fits{i}.rmse,...
        'bic',bic(i),'elapsed_s',time(i),'g_used_mm',fits{i}.g_used,...
        'frequencies_hz',fits{i}.f_id,'amplitudes_mm',fits{i}.A_id);
end
end

function m=other_margin(v,iBest)
q=v;q(iBest)=Inf;m=min(q)-v(iBest);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
