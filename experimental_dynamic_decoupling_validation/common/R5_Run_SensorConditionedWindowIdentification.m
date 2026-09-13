function out = R5_Run_SensorConditionedWindowIdentification(cfg, modelBuilder, sidecarFile, templateFile, foundationFile, outputFile)
% Run the formal sensor-conditioned window fit. No delta-mu/delta-tau state exists.
if nargin<6 || isempty(outputFile), outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat'); end
sidecarSourceFile=sidecarFile; foundationSourceFile=foundationFile;
sidecarFile = R5_StageMatForMatlab(sidecarFile,'r5_sidecar');
foundationFile = R5_StageMatForMatlab(foundationFile,'r5_foundation');
S=load(sidecarFile,'SensorConditionedLibrary'); L=S.SensorConditionedLibrary;
assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
Q=load(foundationFile,'Result'); R=Q.Result; W=R.WindowResult;
try
    models=modelBuilder(sidecarFile,templateFile,cfg.case.targetBlade,W(1).bundle.CorrectedGapLibrary);
catch ME
    if ~strcmp(ME.identifier,'MATLAB:TooManyInputs'), rethrow(ME); end
    models=modelBuilder(sidecarFile,templateFile,cfg.case.targetBlade);
end
ids=[models.sensorId]; assert(isequal(sort(ids),sort(cfg.case.analysisSensors)),'R5:SensorRoleMismatch');
staticAudit=[];
if isfield(cfg,'r5') && isfield(cfg.r5,'requireStaticAudit') && cfg.r5.requireStaticAudit
    staticAudit=R5_Audit_StaticObservationContract(models, ...
        R5_MakeStaticAuditQueries(models),1e-8);
    assert(staticAudit.pass,'R5:StaticContractAuditFailed', ...
        'Static observation contract audit failed before dynamic fitting.');
end
% Do not silently turn a gap-aware run into a no-gap fit when a builder
% forgets to declare the sensor role.  This was the cause of a misleading
% 20250527 result with all delta_gap values fixed at zero.
if isfield(cfg.case,'gapSensors')
    declaredGap=arrayfun(@(m) isfield(m,'isGapSensor') && m.isGapSensor, models);
    declaredGap=declaredGap(:).';
    expectedGap=ismember(ids,double(cfg.case.gapSensors));
    expectedGap=expectedGap(:).';
    assert(isequal(declaredGap,expectedGap), ...
        'R5:SensorRoleDeclaration', ...
        'Builder gap-role flags do not match cfg.case.gapSensors.');
end
fit_bounds(cfg);
rows=repmat(empty_row(),numel(W),1);
for iw=1:numel(W)
    B=W(iw).CoreBundlePreview; bundle=make_bundle(B,ids,cfg,W(iw));
    replay=foundation_replay_audit(bundle,models,cfg);
    rows(iw).window_id=W(iw).window_id; rows(iw).lap_range=W(iw).lap_range;
    seedTable=W(iw).SeedTable;
    % Prefer the legacy VP table when present.  AllEOSeedTableCore belongs
    % to the Foundation preparation stage and can have a different point
    % set/weighting than the frozen gap_only VP scan.
    useAllCore=isempty(seedTable);
    if useAllCore && isfield(W(iw),'AllEOSeedTableCore') && ~isempty(W(iw).AllEOSeedTableCore)
        seedTable=W(iw).AllEOSeedTableCore;
    end
    if isempty(seedTable)
        rows(iw).status='no_seed'; continue;
    end
expectedEO=NaN;
    if isfield(R,'Trend') && istable(R.Trend) && ismember('EO_id',R.Trend.Properties.VariableNames)
        kEO=find(double(R.Trend.window_id)==double(W(iw).window_id),1);
        if ~isempty(kEO), expectedEO=double(R.Trend.EO_id(kEO)); end
    end
fit=fit_window(bundle,models,seedTable,expectedEO,cfg,rows(iw).window_id);
    fit.foundation_replay_available=replay.available;
    fit.foundation_replay_rmse_mv=replay.rmse_mv;
    fit.foundation_replay_max_abs_mv=replay.max_abs_mv;
    fit.foundation_replay_pass=replay.pass;
    if replay.available && ~replay.pass && strcmp(fit.status,'pass')
        fit.status='foundation_replay_fail';
    end
    fit.quality_class=classify_fit(fit);
    fit.window_id=rows(iw).window_id; fit.lap_range=rows(iw).lap_range; rows(iw)=fit;
end

function q=classify_fit(r)
% Post-fit label used for reporting; it never changes the optimizer winner.
if strcmp(r.status,'foundation_replay_fail') || ~r.support_pass
    q='rejected';
elseif r.hit_boundary || r.a_bound_hit || r.dx_bound_hit || r.dg_bound_hit
    q='diagnostic_only';
elseif isfinite(r.foundation_replay_available) && r.foundation_replay_available && ...
        ~r.foundation_replay_pass
    q='rejected';
else
    q='usable';
end
end

out=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_RESULT_V1','dataset',cfg.dataset,...
    'targetBlade',cfg.case.targetBlade,'sensorIds',ids,'rows',rows,...
    'gapSensorIds',double(cfg.case.gapSensors(:).'),...
    'gapBaselineModel',getfield_or_text(cfg,'r5','gapBaselineModel','sensor_conditioned_scalar'),...
    'referenceEO',double(getfield_or_nan(cfg,'audit','referenceEO')),...
    'referenceEOTolerance',NaN,...
    'lockReferenceEO',logical(isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO),...
    'fitBounds',fit_bounds(cfg),...
    'forbiddenDynamicParameters',{{'deltaMu','deltaTau'}},...
    'useNestedFoundationAnchor',logical(isfield(cfg,'r5') && isfield(cfg.r5,'useNestedFoundationAnchor') && cfg.r5.useNestedFoundationAnchor),...
    'useV1DynamicIncrement',logical(isfield(cfg,'r5') && isfield(cfg.r5,'useV1DynamicIncrement') && cfg.r5.useV1DynamicIncrement),...
    'foundationFile',foundationSourceFile,'stagedFoundationFile',foundationFile, ...
    'sidecarFile',sidecarSourceFile,'stagedSidecarFile',sidecarFile, ...
    'staticAudit',staticAudit);
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'out','-v7.3');
end

function v=getfield_or_nan(s,group,name)
v=NaN;
if isfield(s,group) && isstruct(s.(group)) && isfield(s.(group),name), v=s.(group).(name); end
end

function v=getfield_or_text(s,group,name,d)
v=d;
if isfield(s,group) && isstruct(s.(group)) && isfield(s.(group),name)
    v=s.(group).(name);
end
end

function B=make_bundle(P,ids,cfg,W)
mask=logical(P.core_mask(:)) & isfinite(P.x(:)) & isfinite(P.v(:)) & isfinite(P.theta(:)) & isfinite(P.sensor_id(:));
[ok,si]=ismember(double(P.sensor_id(:)),double(ids)); mask=mask&ok;
B.X=double(P.x(mask)); B.T=double(P.t(mask)); B.V=1000*double(P.v(mask)); B.Theta=double(P.theta(mask)); B.sensorIndex=si(mask);
B.W=max(double(P.fit_weight(mask)),1e-6); B.sensorIds=double(ids(:).'); B.rotFreqMeanHz=double(P.rot_freq_mean_hz);
B.useV1LowTemplate=isfield(P,'v1_template') && ~isempty(P.v1_template);
B.useV1DynamicIncrement=isfield(cfg,'r5') && isfield(cfg.r5,'useV1DynamicIncrement') && logical(cfg.r5.useV1DynamicIncrement);
if B.useV1LowTemplate
    B.v1Template=P.v1_template; B.v1SensorInfo=[];
    if isfield(P,'v1_sensor_info'), B.v1SensorInfo=P.v1_sensor_info; end
    B.v1ResponseSurface=[]; B.v1SensorCalibration=[]; B.v1G0BySensor=[];
    if isfield(P,'v1_response_surface'), B.v1ResponseSurface=P.v1_response_surface; end
    if isfield(P,'v1_sensor_calibration'), B.v1SensorCalibration=P.v1_sensor_calibration; end
    if isfield(P,'v1_g0_by_sensor'), B.v1G0BySensor=double(P.v1_g0_by_sensor(:)); end
end
% Optional legacy nested anchor.  Foundation files from older routes may
% carry the already fitted high-speed prediction in the preview.  Use it
% only when it is point-aligned with the current preview; otherwise the
% formal low-speed-template contract remains active.
B.useNestedFoundationAnchor=false; B.anchorV=nan(size(B.V));
B.anchorEO=NaN; B.anchorF=NaN; B.anchorA=NaN; B.anchorPhi=NaN; B.anchorDx=NaN;
useNested=isfield(cfg,'r5') && isfield(cfg.r5,'useNestedFoundationAnchor') && cfg.r5.useNestedFoundationAnchor;
needReplay=isfield(cfg,'r5') && isfield(cfg.r5,'auditFoundationReplay') && cfg.r5.auditFoundationReplay;
if useNested || needReplay
    names={'v_pred','V_pred','VPred','vPred'}; av=[]; AP=[];
    if isfield(W,'Result') && isfield(W.Result,'ReconstructionPreview'), AP=W.Result.ReconstructionPreview; end
    for k=1:numel(names)
        if isfield(P,names{k}), av=double(P.(names{k})); break; end
    end
    if ~isempty(AP) && isfield(AP,'v_pred') && numel(AP.v_pred)>0
        av=map_anchor_preview(P,AP,ids);
    end
    if ~isempty(av) && numel(av)==numel(P.v)
        av=1000*av(:); av=av(mask);
        if all(isfinite(av)), B.anchorV=av; B.useNestedFoundationAnchor=true; end
    end
    if isfield(W,'CoreResult')
        C=W.CoreResult;
        if isfield(C,'EO_id'), B.anchorEO=double(C.EO_id); end
        if isfield(C,'fn_id'), B.anchorF=double(C.fn_id); end
        if isfield(C,'A_id'), B.anchorA=double(C.A_id); end
        if isfield(C,'phi_id_wrapped'), B.anchorPhi=double(C.phi_id_wrapped); end
        if isfield(C,'dx_c_id'), B.anchorDx=double(C.dx_c_id); end
    end
end
assert(numel(B.X)>20 && isfinite(B.rotFreqMeanHz),'R5:InvalidBundle');
end

function av=map_anchor_preview(P,AP,ids)
% Match the frozen Foundation prediction to the current core points.  The
% previews are sorted independently, so positional assignment is invalid.
n=numel(P.v); av=nan(n,1); sid=double(P.sensor_id(:)); pass=double(P.pass_id(:));
th=double(P.theta(:)); tt=double(P.t(:)); as=double(AP.sensor_id(:)); ap=double(AP.pass_id(:));
ath=double(AP.theta(:)); at=double(AP.t(:)); val=double(AP.v_pred(:));
for k=1:n
    q=(as==sid(k))&(ap==pass(k));
    if isfield(P,'blade_id_vec') && isfield(AP,'blade_id_vec')
        q=q&(double(AP.blade_id_vec(:))==double(P.blade_id_vec(k)));
    end
    j=find(q);
    if isempty(j), continue; end
    % Time and phase identify a pulse; the small x-coordinate differences
    % between routes are deliberately not used for the match.
    d=(at(j)-tt(k)).^2 + (ath(j)-th(k)).^2;
    [~,ii]=min(d); jj=j(ii);
    if isfinite(val(jj)), av(k)=val(jj); end
end
end

function r=fit_window(B,M,T,expectedEO,cfg,windowId)
if nargin<6, windowId=NaN; end
assert(~isempty(T),'R5:NoSeed');
if istable(T), T=table2struct(T); end
[B.FixedMask,supportInfo]=common_support_mask(B,M,cfg);
minFrac=0.98;
if isfield(cfg,'r5') && isfield(cfg.r5,'minSupportFraction')
    minFrac=double(cfg.r5.minSupportFraction);
end
supportPass=nnz(B.FixedMask)>=max(20,ceil(minFrac*numel(B.FixedMask)));
if ~supportPass
    r=empty_row(); r.status='support_fail'; r.support_pass=false;
    r.support_valid_count=nnz(B.FixedMask); r.support_total_count=numel(B.FixedMask);
    r.support_fraction=nnz(B.FixedMask)/max(numel(B.FixedMask),1);
    r.support_info=supportInfo; return;
end
% The legacy route retains several EO/phase candidates before the final
% continuous fit.  A single seed can trap lsqnonlin in a poor local basin,
% especially when gap and vibration terms partially compensate each other.
screen=screen_all_eo(B,M,T,cfg);
screenRmse=double([screen.r5_screen_rmse]); screenRmse(~isfinite(screenRmse))=Inf;
% The frozen V1 contract selects the EO candidates from its VP waveform
% scan, before the gap-aware full-wave refinement.  R5's linearized dg
% score is retained as a diagnostic only; using it to replace the V1
% candidate ranking can create EO switches that are absent in gap_only.
legacyRmse=NaN(size(screenRmse));
for ik=1:numel(T)
    if isfield(T(ik),'plain_voltage_rmse')
        legacyRmse(ik)=double(T(ik).plain_voltage_rmse);
    elseif isfield(T(ik),'fallbackRmseMv')
        legacyRmse(ik)=double(T(ik).fallbackRmseMv);
    elseif isfield(T(ik),'weighted_voltage_rmse')
        legacyRmse(ik)=double(T(ik).weighted_voltage_rmse);
    end
end
if nnz(isfinite(legacyRmse))>=max(1,ceil(0.8*numel(legacyRmse)))
    screenRmse=legacyRmse;
    for ik=1:numel(screen), screen(ik).candidate_screen_mode='legacy_vp'; end
else
    for ik=1:numel(screen), screen(ik).candidate_screen_mode='r5_linearized_dg'; end
end
% Keep the complete EO candidate set in play.  The Foundation trend is a
% useful diagnostic/initialization hint, but it is not an observation-level
% constraint: locking the dynamic fit to expectedEO can retain a wrong EO
% whenever the transferred gap increment changes the voltage ranking.  This
% was a major difference from the legacy full-wave route, which compared
% several EO candidates after the gap-aware refinement.
ordAll=1:numel(T);
[~,ord0]=sort(screenRmse(ordAll),'ascend'); ordAll=ordAll(ord0);
nKeep=min(candidate_keep_count(cfg),numel(ordAll));
ord=ordAll(1:nKeep);
for jj=1:numel(ord), screen(ord(jj)).selected_for_full_fit=true; end
if isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO
    eoAll=double([T.EO]); referenceEO=getfield_or_nan(cfg,'audit','referenceEO');
    assert(isfinite(referenceEO),'R5:ReferenceEOUnset', ...
        'A locked EO requires cfg.audit.referenceEO.');
    match=find(eoAll==referenceEO);
    assert(~isempty(match),'R5:ReferenceEOMissing','Locked reference EO is absent from the Foundation seeds.');
    ord=match(:).';
end
eoLocked=isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO;
% expectedEO is retained only as a provenance field by the caller.  It must
% not be injected into the candidate set: final EO selection is data-driven.
best=[]; bestRmse=Inf; candidates=repmat(empty_row(),0,1);
for kk=1:numel(ord)
    cand=fit_window_from_seed(B,M,T(ord(kk)),cfg);
    cand.support_pass=supportPass; cand.support_valid_count=nnz(B.FixedMask);
    cand.support_total_count=numel(B.FixedMask); cand.support_fraction=nnz(B.FixedMask)/numel(B.FixedMask);
    cand.support_info=supportInfo;
    if isfinite(cand.rmse_mv), candidates(end+1)=cand; end %#ok<AGROW>
    if isfinite(cand.rmse_mv) && cand.rmse_mv<bestRmse
        best=cand; bestRmse=cand.rmse_mv;
    end
end
if isempty(best), best=fit_window_from_seed(B,M,T(ord(1)),cfg); end
auditAll=logical(isfield(cfg,'r5') && isfield(cfg.r5,'auditAllEO') && cfg.r5.auditAllEO);
if auditAll && isfield(cfg.r5,'auditAllEOWindowIds') && ~isempty(cfg.r5.auditAllEOWindowIds)
    auditAll=ismember(double(windowId),double(cfg.r5.auditAllEOWindowIds));
end
allCandidates=repmat(empty_row(),0,1); allBest=Inf; allBestEO=NaN;
if auditAll
    % Diagnostic only: refit every EO with the same fixed support, bounds,
    % frequency refinement and multistart settings.  The production winner
    % above remains Top-K; this table measures whether Top-K omitted the
    % global full-wave winner.
    allCandidates=candidates;
    already=double([T(ord).EO]);
    for kk=1:numel(ordAll)
        if ismember(double(T(ordAll(kk)).EO),already), continue; end
        cand=fit_window_from_seed(B,M,T(ordAll(kk)),cfg);
        cand.support_pass=supportPass; cand.support_valid_count=nnz(B.FixedMask);
        cand.support_total_count=numel(B.FixedMask); cand.support_fraction=nnz(B.FixedMask)/numel(B.FixedMask);
        cand.support_info=supportInfo;
        if isfinite(cand.rmse_mv), allCandidates(end+1)=cand; end %#ok<AGROW>
    end
    if ~isempty(allCandidates)
        [allBest,ib]=min([allCandidates.rmse_mv]); allBestEO=allCandidates(ib).EO;
    end
end
best.eo_screen_table=screen;
best.candidate_table=candidates;
% Report near-equivalent EO solutions instead of silently treating a tiny
% SSE difference as a physically decisive order choice.
if ~isempty(candidates) && isfinite(bestRmse)
    near=double([candidates.rmse_mv]) <= bestRmse*(1+0.03);
    best.eo_near_tie=double([candidates(near).EO]);
    best.eo_near_tie_rmse_mv=double([candidates(near).rmse_mv]);
    best.eo_ambiguity=nnz(near)>1;
else
    best.eo_near_tie=[]; best.eo_near_tie_rmse_mv=[]; best.eo_ambiguity=false;
end
best.audit_all_eo_table=allCandidates;
best.audit_all_eo_best_rmse_mv=allBest;
best.audit_all_eo_best_EO=allBestEO;
if auditAll && isfinite(allBest)
    best.topk_recall_global=ismember(allBestEO,double([T(ord).EO]));
else
    best.topk_recall_global=NaN;
end
r=best;
end

function n=candidate_keep_count(cfg)
n=3;
if isfield(cfg,'frequency')
    if isfield(cfg.frequency,'candidateTopK') && isfinite(cfg.frequency.candidateTopK)
        n=max(1,round(double(cfg.frequency.candidateTopK)));
    elseif isfield(cfg.frequency,'vpTopK') && isfinite(cfg.frequency.vpTopK)
        n=max(1,round(double(cfg.frequency.vpTopK)));
    end
end

end

function screen=screen_all_eo(B,M,T,cfg)
% Fast R5-owned EO screening. Foundation supplies only A/phi/dx seeds;
% the score below already includes a linearized sensor-wise dg correction.
n=numel(T); screen=repmat(struct('EO',NaN,'f_seed_hz',NaN,'A_seed',NaN, ...
    'phi_seed',NaN,'dx_seed',NaN,'foundation_seed_rmse',NaN, ...
    'dg_projected',[],'r5_screen_rmse',Inf,'screen_rank',NaN, ...
    'selected_for_full_fit',false,'candidate_screen_mode',''),n,1);
[~,~,limDg]=fit_bounds_from_cfg(); dgHard=hard_gap_limits(M,limDg);
for k=1:n
    s=T(k); screen(k).EO=double(s.EO); screen(k).f_seed_hz=double(s.EO)*B.rotFreqMeanHz;
    screen(k).A_seed=double(s.A); screen(k).phi_seed=double(s.phi); screen(k).dx_seed=double(s.dx_c);
    if isfield(s,'plain_voltage_rmse'), screen(k).foundation_seed_rmse=double(s.plain_voltage_rmse); end
    z=[screen(k).A_seed screen(k).phi_seed screen(k).dx_seed zeros(1,numel(M)) 0];
    B.EO=screen(k).EO; p0=prediction(z,screen(k).f_seed_hz,B,M);
    fitMask=isfield(B,'FixedMask') & B.FixedMask & isfinite(B.V) & isfinite(p0);
    dg=zeros(1,numel(M));
    for ii=1:numel(M)
        if ~M(ii).isGapSensor || dgHard(ii)<=0, continue; end
        h=min(1e-3,0.25*dgHard(ii)); ze=z; ze(3+ii)=h;
        pe=prediction(ze,screen(k).f_seed_hz,B,M); q=fitMask & B.sensorIndex==ii & isfinite(pe);
        if nnz(q)>=4
            d=(pe(q)-p0(q))/h; rr=B.V(q)-p0(q); den=sum(d.^2);
            if isfinite(den) && den>eps, dg(ii)=min(max(sum(d.*rr)/den,-dgHard(ii)),dgHard(ii)); end
        end
    end
    z(4:3+numel(M))=dg; pp=prediction(z,screen(k).f_seed_hz,B,M);
    ok=fitMask & isfinite(pp);
    if nnz(ok)>=20, screen(k).r5_screen_rmse=sqrt(mean((pp(ok)-B.V(ok)).^2)); end
    screen(k).dg_projected=dg;
end
[~,ord]=sort([screen.r5_screen_rmse],'ascend');
for k=1:numel(ord), screen(ord(k)).screen_rank=k; end
end

function r=fit_window_from_seed(B,M,s,cfg)
n=numel(M); B.EO=double(s.EO); f0=double(s.EO)*B.rotFreqMeanHz; [limA,limDx,limDg]=fit_bounds_from_cfg();
% Keep a parameter slot for every analysed channel so the common vibration
% state is constrained by direct/template channels as well.  Direct-only
% channels have no clearance degree of freedom; their dg slot is fixed at 0.
isGap=arrayfun(@(m) isfield(m,'isGapSensor') && m.isGapSensor, M);
dgHard=hard_gap_limits(M,limDg);
dgTrust=local_gap_projection_limits(B,M,double(s.A),double(s.phi),double(s.dx_c),f0,dgHard,cfg);
dgTrust=dgTrust(:).'; isGap=logical(isGap(:).');
dgLo=-dgHard.*double(isGap); dgHi=dgHard.*double(isGap);
dfHalf=2.0;
if isfield(cfg,'frequency') && isfield(cfg.frequency,'refineEachCandidate') && cfg.frequency.refineEachCandidate
    if isfield(cfg.frequency,'refineHalfWidthHz') && isfinite(cfg.frequency.refineHalfWidthHz)
        dfHalf=double(cfg.frequency.refineHalfWidthHz);
    end
end
z0=[double(s.A),double(s.phi),double(s.dx_c),zeros(1,n),0];
lb=[0,-pi,-limDx,dgLo,-dfHalf]; ub=[limA,pi,limDx,dgHi,dfHalf]; z0=min(max(z0,lb),ub);
% Keep the Foundation vibration seed when requested.  The R5 zero-dg
% model is intentionally a low-speed template, so re-fitting vibration
% against it before the gap increment is introduced can absorb a static
% high-speed offset into amplitude/dx.  The old V1 route used the frozen
% Foundation solution for exactly this reason.
useSeedOnly=isfield(cfg,'r5') && isfield(cfg.r5,'useFoundationVibrationSeedOnly') && cfg.r5.useFoundationVibrationSeedOnly;
if ~useSeedOnly
    vibFun=@(q)resid([q(:).',zeros(1,n)],f0,B,M);
    vibOpts=optimoptions('lsqnonlin','Display','off','MaxIterations',80,'MaxFunctionEvaluations',1200);
    [qV,~,~,~]=safe_lsqnonlin(vibFun,z0(1:3),lb(1:3),ub(1:3),vibOpts);
    z0(1:3)=qV(:).';
end
% Obtain a local first-order dg estimate at the legacy vibration seed.
% Starting all dg values at zero can leave the joint optimizer in the
% Foundation basin when clearance and vibration terms partially compensate.
p0=prediction(z0,f0,B,M); epsDg=min(1e-3,0.02*max(dgHard));
for ii=1:n
    if ~isGap(ii), continue; end
    h=min(epsDg,0.25*dgHard(ii));
    if h<=0, continue; end
    ze=z0; ze(3+ii)=h; pe=prediction(ze,f0,B,M);
    q=(B.sensorIndex==ii) & B.FixedMask & isfinite(p0) & isfinite(pe) & isfinite(B.V);
    if nnz(q)>=4
        d=(pe(q)-p0(q))/h; rr=B.V(q)-p0(q); den=sum(d.^2);
        if isfinite(den) && den>eps
            z0(3+ii)=sum(d.*rr)/den;
        end
    end
end
if isfield(cfg,'r5') && isfield(cfg.r5,'useLocalGapProjection') && cfg.r5.useLocalGapProjection
    z0(4:3+n)=min(max(z0(4:3+n),-dgTrust),dgTrust);
end
z0=min(max(z0,lb),ub);
fun=@(z)resid(z,f0+z(end),B,M); opts=optimoptions('lsqnonlin','Display','off','MaxIterations',120,'MaxFunctionEvaluations',3000);
starts=z0;
if isfield(cfg,'r5') && isfield(cfg.r5,'multiStartCount')
    ns=max(1,round(double(cfg.r5.multiStartCount)));
    starts=repmat(z0,ns,1);
    if ns>=2, starts(2,1)=0.8*z0(1); end
    if ns>=3, starts(3,1)=min(1.2*z0(1),limA); end
    if ns>=4, starts(4,2)=wrap_phase(z0(2)+pi/4); end
    if ns>=5, starts(5,2)=wrap_phase(z0(2)-pi/4); end
    starts=min(max(starts,lb),ub);
end
bestZ=z0; bestRes=[]; bestCost=Inf; bestExit=NaN; startRows=repmat(struct('start_id',NaN,'A',NaN,'phi',NaN,'dx',NaN,'rmse_mv',NaN,'exitflag',NaN),size(starts,1),1);
for istart=1:size(starts,1)
    [zz,~,rr,ef]=safe_lsqnonlin(fun,starts(istart,:),lb,ub,opts);
    cc=norm(rr);
    startRows(istart)=struct('start_id',istart,'A',zz(1),'phi',zz(2),'dx',zz(3),'rmse_mv',sqrt(mean(rr.^2)),'exitflag',ef);
    if isfinite(cc) && cc<bestCost, bestZ=zz; bestRes=rr; bestCost=cc; bestExit=ef; end
end
z=bestZ; res=bestRes; exitflag=bestExit; fFit=f0+z(end); pred=prediction(z,fFit,B,M);
fitMask=B.FixedMask & isfinite(B.V); invalid=~isfinite(pred(fitMask));
r=empty_row(); r.status='pass'; r.EO=double(s.EO); r.frequency_hz=fFit; r.amplitude_mm=z(1); r.phase_rad=z(2); r.dx_mm=z(3); r.delta_gap_mm=z(4:3+n).';
if any(invalid), r.status='invalid_prediction'; r.rmse_mv=Inf; r.exitflag=exitflag; r.start_table=startRows; return; end
r.support_pass=true;
r.anchor_rmse_mv=NaN; r.anchor_max_abs_mv=NaN;
if B.useNestedFoundationAnchor && all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx]))
    za=[B.anchorA B.anchorPhi B.anchorDx zeros(1,n)]; pa=prediction(za,B.anchorEO*B.rotFreqMeanHz,B,M); qa=isfinite(pa)&isfinite(B.anchorV)&B.FixedMask;
    if any(qa), r.anchor_rmse_mv=sqrt(mean((pa(qa)-B.anchorV(qa)).^2)); r.anchor_max_abs_mv=max(abs(pa(qa)-B.anchorV(qa))); end
end
r.rmse_mv=sqrt(mean((pred(fitMask)-B.V(fitMask)).^2)); r.sse_mv2=sum((pred(fitMask)-B.V(fitMask)).^2); sr=nan(1,numel(M)); for ii=1:numel(M), q=B.sensorIndex==ii & fitMask; sr(ii)=sqrt(mean((pred(q)-B.V(q)).^2)); end
zBase=z; zBase(4:3+n)=0; predBase=prediction(zBase,fFit,B,M); baseOk=isfinite(predBase)&isfinite(B.V);
baseOk=baseOk & B.FixedMask;
r.baseline_rmse_mv=sqrt(mean((predBase(baseOk)-B.V(baseOk)).^2));
r.increment_rms_mv=sqrt(mean((pred(fitMask)-predBase(fitMask)).^2));
r.sensor_baseline_rmse_mv=nan(1,numel(M)); r.sensor_increment_rms_mv=nan(1,numel(M));
for ii=1:numel(M)
    q=B.sensorIndex==ii & baseOk; if any(q), r.sensor_baseline_rmse_mv(ii)=sqrt(mean((predBase(q)-B.V(q)).^2)); end
    q=B.sensorIndex==ii & fitMask; if any(q), r.sensor_increment_rms_mv(ii)=sqrt(mean((pred(q)-predBase(q)).^2)); end
end
r.sensor_rmse_mv=sr; r.exitflag=exitflag; r.residual_norm=norm(res);
% Treat numerical solutions within the optimizer's practical resolution of
% a bound as boundary hits; exact equality alone misses saturated fits.
tolBound=1e-5;
r.delta_gap_limit_mm=dgHard(:);
r.hit_boundary=any(isGap & abs(abs(z(4:3+n))-dgHard)<tolBound)| ...
    abs(abs(z(3))-limDx)<tolBound|abs(z(1)-limA)<tolBound;
r.a_bound_hit=abs(z(1)-limA)<tolBound; r.dx_bound_hit=abs(abs(z(3))-limDx)<tolBound;
r.dg_bound_hit=any(isGap & abs(abs(z(4:3+n))-dgHard)<tolBound);
r.multi_start_count=size(starts,1);
r.start_table=startRows;
end

function [z,cost,res,exitflag]=safe_lsqnonlin(fun,z0,lb,ub,opts)
try
    [z,cost,res,exitflag]=lsqnonlin(fun,z0,lb,ub,opts);
    if isempty(res) || ~isvector(res), error('R5:ResidualShape','Invalid residual shape'); end
catch
    obj=@(q)sum(fun(min(max(q,lb),ub)).^2,'omitnan');
    fo=optimset('Display','off','MaxIter',2000,'MaxFunEvals',5000);
    z=fminsearch(obj,z0,fo); z=min(max(z,lb),ub); res=fun(z); cost=norm(res); exitflag=-2;
end
end

function [limA,limDx,limDg]=fit_bounds_from_cfg()
% Kept as a local override point for the generic worker; the public cfg is
% copied into a persistent run-level variable before window processing.
global R5_SC_FIT_BOUNDS
if isempty(R5_SC_FIT_BOUNDS), R5_SC_FIT_BOUNDS=[0.50 0.35 0.25]; end
limA=R5_SC_FIT_BOUNDS(1); limDx=R5_SC_FIT_BOUNDS(2); limDg=R5_SC_FIT_BOUNDS(3);
end

function b=fit_bounds(cfg)
global R5_SC_FIT_BOUNDS
if isfield(cfg,'r5') && isfield(cfg.r5,'amplitudeLimitMm'), a=cfg.r5.amplitudeLimitMm; else, a=.5; end
if isfield(cfg,'r5') && isfield(cfg.r5,'dxLimitMm'), d=cfg.r5.dxLimitMm; else, d=.35; end
if isfield(cfg,'r5') && isfield(cfg.r5,'deltaGapLimitMm'), g=cfg.r5.deltaGapLimitMm; else, g=.25; end
assert(all(isfinite([a d g])) && all([a d g]>0),'R5:InvalidFitBounds');
R5_SC_FIT_BOUNDS=[a d g]; b=struct('amplitudeLimitMm',a,'dxLimitMm',d,'deltaGapLimitMm',g);
end

function e=resid(z,f,B,M)
p=prediction(z,f,B,M); good=isfinite(p)&isfinite(B.V);
if isfield(B,'FixedMask'), good=good & B.FixedMask; end
% lsqnonlin requires an invariant residual vector.  Keep all observation
% arrays explicitly column-shaped; otherwise a row/column mismatch in a
% template interpolant can make the residual length change between calls.
p=p(:); vv=B.V(:); good=isfinite(p)&isfinite(vv);
if isfield(B,'FixedMask'), good=good & logical(B.FixedMask(:)); end
% A malformed interpolant must be treated as an invalid candidate, not as a
% variable-length residual.  Preserve the observation length expected by
% lsqnonlin and apply the existing large invalid-point penalty.
N=numel(vv);
if numel(p)~=N
    pp=nan(N,1); gg=false(N,1); k=min(numel(p),N);
    pp(1:k)=p(1:k); gg(1:k)=good(1:k); p=pp; good=gg;
end
e=p-vv;
% Keep a fixed residual length. Invalid support is a candidate failure, not
% an opportunity to silently remove difficult observations from the fit.
if isfield(B,'FixedMask')
    fm=logical(B.FixedMask(:));
    e(~good & fm)=1e6;
    e(~fm)=0;
else
    e(~good)=1e6;
end
if isempty(e), e=1e6; end
e=e(:);
end

function p=wrap_phase(p)
p=mod(p+pi,2*pi)-pi;
end

function p=prediction(z,f,B,M)
[p,~]=prediction_components(z,f,B,M);
p=p(:);
end

function [p,c]=prediction_components(z,f,B,M)
[p,c]=R5_Evaluate_CanonicalForward(z,f,B,M);
end

function a=foundation_replay_audit(B,M,cfg)
a=struct('available',false,'rmse_mv',NaN,'max_abs_mv',NaN,'pass',false);
need=isfield(cfg,'r5') && isfield(cfg.r5,'auditFoundationReplay') && cfg.r5.auditFoundationReplay;
if ~need || ~all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx])) || ...
        ~any(isfinite(B.anchorV)), return; end
B0=B; B0.EO=B.anchorEO; B0.useNestedFoundationAnchor=false;
z0=[B.anchorA B.anchorPhi B.anchorDx zeros(1,numel(M))];
p0=R5_Evaluate_CanonicalForward(z0,B.anchorF,B0,M);
[fixedMask,~]=common_support_mask(B,M,cfg);
q=fixedMask & isfinite(B.anchorV) & isfinite(p0);
if nnz(q)<20, return; end
d=p0(q)-B.anchorV(q); a.available=true;
a.rmse_mv=sqrt(mean(d.^2)); a.max_abs_mv=max(abs(d));
tol=1e-6;
if isfield(cfg.r5,'foundationReplayToleranceMv'), tol=double(cfg.r5.foundationReplayToleranceMv); end
a.pass=isfinite(a.rmse_mv) && a.rmse_mv<=tol;
end

function r=empty_row()
r=struct('window_id',NaN,'lap_range',[NaN NaN],'status','unprocessed','quality_class','unprocessed','EO',NaN,'frequency_hz',NaN,'amplitude_mm',NaN,'phase_rad',NaN,'dx_mm',NaN,'delta_gap_mm',NaN,'delta_gap_limit_mm',NaN,'rmse_mv',NaN,'sse_mv2',NaN,'baseline_rmse_mv',NaN,'increment_rms_mv',NaN,'anchor_rmse_mv',NaN,'anchor_max_abs_mv',NaN,'foundation_replay_available',false,'foundation_replay_rmse_mv',NaN,'foundation_replay_max_abs_mv',NaN,'foundation_replay_pass',false,'sensor_rmse_mv',NaN,'sensor_baseline_rmse_mv',NaN,'sensor_increment_rms_mv',NaN,'exitflag',NaN,'residual_norm',NaN,'hit_boundary',false,'support_pass',false,'support_valid_count',NaN,'support_total_count',NaN,'support_fraction',NaN,'support_info',struct(),'a_bound_hit',false,'dx_bound_hit',false,'dg_bound_hit',false,'multi_start_count',NaN,'start_table',struct([]),'eo_screen_table',struct([]),'candidate_table',struct([]),'eo_near_tie',[],'eo_near_tie_rmse_mv',[],'eo_ambiguity',false,'audit_all_eo_table',struct([]),'audit_all_eo_best_rmse_mv',NaN,'audit_all_eo_best_EO',NaN,'topk_recall_global',NaN,'candidate_screen_mode','');
end

function limits=local_gap_projection_limits(B,M,A,phi,dx,f0,globalLimit,cfg)
limits=double(globalLimit(:)).';
if numel(limits)~=numel(M), error('R5:GapLimitShape','Gap limits must have one value per model.'); end
if ~(isfield(cfg,'r5') && isfield(cfg.r5,'useLocalGapProjection') && cfg.r5.useLocalGapProjection), return; end
alpha=cfg_num(cfg.r5,'gapProjectionAlpha',2.5); minLim=cfg_num(cfg.r5,'gapProjectionMinLimitMm',0.02);
hG=cfg_num(cfg.r5,'gapProjectionStepMm',1e-3); floorS=cfg_num(cfg.r5,'gapProjectionSensitivityFloorMvPerMm',0.05);
n=numel(M); z=[A phi dx zeros(1,n)]; p0=prediction(z,f0,B,M); valid=isfinite(p0)&isfinite(B.V);
if isfield(B,'FixedMask'), valid=valid & B.FixedMask; end
if nnz(valid)<max(8,n+3), return; end
Jnu=zeros(nnz(valid),3); Jg=zeros(nnz(valid),n); rv=B.V(valid)-p0(valid);
hs=[1e-4 1e-4 1e-4];
for k=1:3
    zp=z; zm=z; zp(k)=zp(k)+hs(k); zm(k)=zm(k)-hs(k);
    pp=prediction(zp,f0,B,M); pm=prediction(zm,f0,B,M); Jnu(:,k)=(pp(valid)-pm(valid))/(2*hs(k));
end
for is=1:n
    if ~M(is).isGapSensor, continue; end
    zp=z; zm=z; zp(3+is)=hG; zm(3+is)=-hG;
    pp=prediction(zp,f0,B,M); pm=prediction(zm,f0,B,M); Jg(:,is)=(pp(valid)-pm(valid))/(2*hG);
end
if any(~isfinite(Jnu(:))) || any(~isfinite(Jg(:))), return; end
Q=orth(Jnu); rPerp=rv-Q*(Q'*rv); G=Jg-Q*(Q'*Jg);
for is=1:n
    if ~M(is).isGapSensor, limits(is)=0; continue; end
    other=setdiff(1:n,is); Qo=orth(G(:,other)); q=G(:,is)-Qo*(Qo'*G(:,is)); qn=norm(q);
    if ~(isfinite(qn)&&qn>floorS), continue; end
    epsPerp=norm(rPerp-Qo*(Qo'*rPerp)); raw=alpha*epsPerp/qn;
    if isfinite(raw), limits(is)=min(limits(is),max(minLim,raw)); end
end
end

function v=cfg_num(s,name,d)
v=d; if isfield(s,name) && isfinite(s.(name)), v=double(s.(name)); end
end

function [mask,info]=common_support_mask(B,M,cfg)
[limA,limDx,limDg]=fit_bounds_from_cfg(); dgHard=hard_gap_limits(M,limDg); n=numel(B.X); mask=true(n,1);
for i=1:numel(M)
    q=B.sensorIndex==i; if ~any(q), continue; end
    xlo=B.X(q)-limA-limDx; xhi=B.X(q)+limA+limDx;
    % In the V1→R5 bridge the low-speed observation curve is the V1
    % direct/template curve.  Do not gate it with the R5 template domain:
    % the two files can have slightly different end points even though the
    % observed V1 points are valid.  For a direct channel there is no R5
    % response-surface query at all, so a response-domain gate would also
    % reject otherwise valid V1 points for no physical reason.
    if ~(isfield(B,'useV1LowTemplate') && B.useV1LowTemplate) && ...
            isfield(M(i),'templateDomainMm') && all(isfinite(M(i).templateDomainMm))
        d=M(i).templateDomainMm; mask(q)=mask(q)&xlo>=d(1)&xhi<=d(2);
    end
    if M(i).isGapSensor && isfield(M(i),'responseDomainMm') && all(isfinite(M(i).responseDomainMm)) && isfield(M(i),'coord')
        d=M(i).responseDomainMm; c=M(i).coord; xx1=c.responseScale*(xlo-c.responseTauMm-c.responseOffsetMm); xx2=c.responseScale*(xhi-c.responseTauMm-c.responseOffsetMm);
        mask(q)=mask(q)&min(xx1,xx2)>=d(1)&max(xx1,xx2)<=d(2);
    end
    if M(i).isGapSensor && isfield(M(i),'gapDomainMm') && all(isfinite(M(i).gapDomainMm))
        g=M(i).gapReferenceMm; d=M(i).gapDomainMm;
        gapOK=g-dgHard(i)>=d(1) && g+dgHard(i)<=d(2); mask(q)=mask(q)&gapOK;
    end
end
info=struct('limA_mm',limA,'limDx_mm',limDx,'limDg_mm',limDg, ...
    'dgHard_mm',dgHard,'validCount',nnz(mask),'totalCount',n,'fraction',nnz(mask)/max(n,1));
end

function lim=hard_gap_limits(M,globalLimit)
lim=zeros(1,numel(M));
for i=1:numel(M)
    if ~M(i).isGapSensor, continue; end
    lim(i)=globalLimit;
    if isfield(M(i),'gapDomainMm') && all(isfinite(M(i).gapDomainMm))
        lim(i)=min(lim(i),M(i).gapReferenceMm-M(i).gapDomainMm(1));
        lim(i)=min(lim(i),M(i).gapDomainMm(2)-M(i).gapReferenceMm);
    end
    lim(i)=max(0,lim(i));
end
end
