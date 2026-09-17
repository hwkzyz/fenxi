function out = R5_Run_SensorConditionedWindowIdentification(cfg, modelBuilder, sidecarFile, templateFile, foundationFile, outputFile)
% Run the formal sensor-conditioned window fit. No delta-mu/delta-tau state exists.
if nargin<6 || isempty(outputFile), outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat'); end
sidecarSourceFile=sidecarFile; foundationSourceFile=foundationFile;
sidecarFile = R5_StageMatForMatlab(sidecarFile,'r5_sidecar');
foundationFile = R5_StageMatForMatlab(foundationFile,'r5_foundation');
S=load(sidecarFile,'SensorConditionedLibrary'); L=S.SensorConditionedLibrary;
assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
models=modelBuilder(sidecarFile,templateFile,cfg.case.targetBlade);
ids=[models.sensorId]; assert(isequal(sort(ids),sort(cfg.case.analysisSensors)),'R5:SensorRoleMismatch');
staticAudit=[];
if isfield(cfg,'r5') && isfield(cfg.r5,'requireStaticAudit') && cfg.r5.requireStaticAudit
    minFullSupport=0.80;
    if isfield(cfg.r5,'minimumFullTemplateSupportFraction'), minFullSupport=cfg.r5.minimumFullTemplateSupportFraction; end
    allowGapExtrapolation=false;
    if isfield(cfg.r5,'allowPositiveGapExtrapolation')
        allowGapExtrapolation=logical(cfg.r5.allowPositiveGapExtrapolation);
    end
    staticAudit=R5_Audit_StaticObservationContract(models, ...
        R5_MakeStaticAuditQueries(models),1e-8,minFullSupport,allowGapExtrapolation);
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
Q=load(foundationFile,'Result'); R=Q.Result; W=R.WindowResult; rows=repmat(empty_row(),numel(W),1);
for iw=1:numel(W)
    B=W(iw).CoreBundlePreview; bundle=make_bundle(B,ids,cfg,W(iw));
    replay=foundation_replay_audit(bundle,models,cfg);
    rows(iw).window_id=W(iw).window_id; rows(iw).lap_range=W(iw).lap_range;
    seedTable=W(iw).SeedTable;
    if isfield(W(iw),'AllEOSeedTableCore') && ~isempty(W(iw).AllEOSeedTableCore)
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
    % A result outside the declared physical search band is not a valid
    % dynamic identification.  Keep it in the audit output, but reject it
    % instead of allowing a wrong 1x/subharmonic branch to look usable when
    % the foundation bundle and frozen sidecar are inconsistent.
    if isfield(cfg,'frequency') && isfield(cfg.frequency,'searchHz') && ...
            numel(cfg.frequency.searchHz)==2 && isfinite(fit.frequency_hz)
        band=double(cfg.frequency.searchHz(:).');
        if fit.frequency_hz < min(band) || fit.frequency_hz > max(band)
            fit.status='frequency_out_of_search_band';
            fit.quality_class='rejected';
            fit.support_pass=false;
        end
    end
    fit.foundation_replay_available=replay.available;
    fit.foundation_replay_status=replay.status;
    fit.foundation_replay_rmse_mv=replay.rmse_mv;
    fit.foundation_replay_max_abs_mv=replay.max_abs_mv;
    fit.foundation_replay_pass=replay.pass;
    fit.foundation_replay_point_count=replay.point_count;
    fit.foundation_replay_sensor_count=replay.sensor_count;
    fit.foundation_replay_coordinate_pass=replay.coordinate_pass;
    fit.foundation_replay_time_phase_pass=replay.time_phase_pass;
    % Foundation replay is a baseline diagnostic only.  It must never
    % change the R5 fit status or quality classification.
    fit.foundation_replay_used_for_acceptance=false;
    fit.foundation_replay_role='diagnostic_only';
    fit.foundation_baseline_EO=bundle.anchorEO;
    fit.foundation_baseline_frequency_hz=bundle.anchorF;
    fit.foundation_baseline_amplitude_mm=bundle.anchorA;
    fit.foundation_baseline_dx_mm=bundle.anchorDx;
    fit.foundation_baseline_phase_rad=bundle.anchorPhi;
    fit.quality_class=classify_fit(fit);
    fit.window_id=rows(iw).window_id; fit.lap_range=rows(iw).lap_range;
    % Keep the output contract stable when a fit exits through an early
    % diagnostic branch or gains a newly added metric.
    if ~isfield(rows,'residual_reduction_ratio')
        [rows.residual_reduction_ratio] = deal(NaN);
        [rows.fit_improves_baseline] = deal(false);
    end
    fit.residual_reduction_ratio = 1 - fit.rmse_mv / max(fit.baseline_rmse_mv,eps);
    fit.fit_improves_baseline = isfinite(fit.residual_reduction_ratio) && fit.residual_reduction_ratio > 0;
    rows(iw)=fit;
end

function q=classify_fit(r)
% Classify the R5 fit from its own contract; Foundation replay is diagnostic.
if ~r.support_pass || ~strcmp(r.status,'pass')
    q='rejected';
elseif isfield(r,'eo_consistency_pass') && ~r.eo_consistency_pass
    q='diagnostic_only';
elseif r.hit_boundary || r.a_bound_hit || r.dx_bound_hit || r.dg_bound_hit
    q='diagnostic_only';
elseif isfield(r,'amplitude_ambiguous') && r.amplitude_ambiguous
    q='diagnostic_only';
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
    'foundationFile',foundationSourceFile,'stagedFoundationFile',foundationFile, ...
    'sidecarFile',sidecarSourceFile,'stagedSidecarFile',sidecarFile, ...
    'methodContract',contract_snapshot(cfg), ...
    'staticAudit',staticAudit);
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'out','-v7.3');
end

function c=contract_snapshot(cfg)
% Store the decision-relevant blind contract with every result.  This makes
% historical comparisons auditable when two files point to the same inputs.
c=struct();
c.schema='R5_BLIND_CONTRACT_SNAPSHOT_V1';
c.targetBlade=double(cfg.case.targetBlade);
c.analysisSensors=double(cfg.case.analysisSensors(:).');
c.gapSensors=double(cfg.case.gapSensors(:).');
c.amplitudeLimitMm=getr5num(cfg,'amplitudeLimitMm',NaN);
c.dxLimitMm=getr5num(cfg,'dxLimitMm',NaN);
c.deltaGapLimitMm=getr5num(cfg,'deltaGapLimitMm',NaN);
c.useStaticGapSlopeInDynamic=getr5logical(cfg,'useStaticGapSlopeInDynamic',false);
c.useLocalGapProjection=getr5logical(cfg,'useLocalGapProjection',false);
c.useFoundationVibrationSeedOnly=getr5logical(cfg,'useFoundationVibrationSeedOnly',false);
c.useNestedFoundationAnchor=getr5logical(cfg,'useNestedFoundationAnchor',false);
c.useTargetXOffsetInDynamic=getr5logical(cfg,'useTargetXOffsetInDynamic',false);
c.lockReferenceEO=getr5logical(cfg,'lockReferenceEO',false);
c.referenceEO=NaN;
if isfield(cfg.r5,'referenceEO'), c.referenceEO=double(cfg.r5.referenceEO); end
c.candidateTopK=NaN;
if isfield(cfg,'frequency') && isfield(cfg.frequency,'candidateTopK'), c.candidateTopK=double(cfg.frequency.candidateTopK); end
c.minSupportFraction=getr5num(cfg,'minSupportFraction',NaN);
c.requireStaticAudit=getr5logical(cfg,'requireStaticAudit',false);
c.auditFoundationReplay=getr5logical(cfg,'auditFoundationReplay',false);
c.dynamicParameters={'EO','frequency','amplitude','phase','dx','delta_gap'};
c.forbiddenDynamicParameters={'deltaMu','deltaTau','dynamic_tilt','registration_drift'};
end

function v=getr5num(cfg,name,d)
v=d; if isfield(cfg,'r5') && isfield(cfg.r5,name) && isnumeric(cfg.r5.(name)) && isscalar(cfg.r5.(name)), v=double(cfg.r5.(name)); end
end
function v=getr5logical(cfg,name,d)
v=d; if isfield(cfg,'r5') && isfield(cfg.r5,name) && isscalar(cfg.r5.(name)), v=logical(cfg.r5.(name)); end
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
        if isfield(C,'V_pred_preview') && numel(C.V_pred_preview)==numel(P.v)
            av=double(C.V_pred_preview(:));
        end
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
% When the seed/VP stage is the declared EO localization stage, retain its
% best EO for the final waveform refinement.  Jointly re-selecting a lower
% EO after allowing dg_s and A to move can trade frequency for clearance and
% produce a smaller RMSE with the wrong physical order.  Keep Top-K fitting
% available as an audit mode, but make the selection policy explicit.
if isfield(cfg,'r5') && isfield(cfg.r5,'finalEOSelection') && ...
        strcmpi(string(cfg.r5.finalEOSelection),'screen_top1')
    ord=ordAll(1);
end
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
best.audit_all_eo_table=allCandidates;
best.audit_all_eo_best_rmse_mv=allBest;
best.audit_all_eo_best_EO=allBestEO;
if auditAll && isfinite(allBest)
    best.topk_recall_global=ismember(allBestEO,double([T(ord).EO]));
else
    best.topk_recall_global=NaN;
end
% A full-wave refit is allowed to improve the seed, but a switch away from
% the best VP/seed EO is evidence of EO--vibration compensation rather than a
% clean identification.  Make that discrepancy explicit and gate the status
% when the case requests an EO consistency audit.
best.eo_screen_top1=NaN;
best.eo_screen_top1_rmse_mv=NaN;
best.eo_final_vs_screen_top1=false;
best.eo_consistency_pass=true;
if ~isempty(screen)
    [~,isTop]=min([screen.r5_screen_rmse]);
    best.eo_screen_top1=double(screen(isTop).EO);
    best.eo_screen_top1_rmse_mv=double(screen(isTop).r5_screen_rmse);
    best.eo_final_vs_screen_top1=double(best.EO)==best.eo_screen_top1;
    gate=isfield(cfg,'r5') && isfield(cfg.r5,'requireEOConsistency') && cfg.r5.requireEOConsistency;
    if gate
        best.eo_consistency_pass=best.eo_final_vs_screen_top1;
        if ~best.eo_consistency_pass && strcmpi(best.status,'pass')
            best.status='eo_ambiguous';
            best.quality_class='ambiguous';
        end
    end
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
    'selected_for_full_fit',false),n,1);
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
    [qV,~,~,~]=lsqnonlin(vibFun,z0(1:3),lb(1:3),ub(1:3),vibOpts);
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
    [zz,~,rr,ef]=lsqnonlin(fun,starts(istart,:),lb,ub,opts);
    cc=norm(rr);
    startRows(istart)=struct('start_id',istart,'A',zz(1),'phi',zz(2),'dx',zz(3),'rmse_mv',sqrt(mean(rr.^2)),'exitflag',ef);
    if isfinite(cc) && cc<bestCost, bestZ=zz; bestRes=rr; bestCost=cc; bestExit=ef; end
end
z=bestZ; res=bestRes; exitflag=bestExit; fFit=f0+z(end); pred=prediction(z,fFit,B,M);
fitMask=B.FixedMask & isfinite(B.V); invalid=~isfinite(pred(fitMask));
r=empty_row(); r.status='pass'; r.EO=double(s.EO); r.frequency_hz=fFit; r.amplitude_mm=z(1); r.phase_rad=z(2); r.dx_mm=z(3); r.delta_gap_mm=z(4:3+n).';
r.VPred=pred;
if any(invalid), r.status='invalid_prediction'; r.rmse_mv=Inf; r.exitflag=exitflag; r.start_table=startRows; return; end
r.support_pass=true;
r.anchor_rmse_mv=NaN; r.anchor_max_abs_mv=NaN;
if B.useNestedFoundationAnchor && all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx]))
    za=[B.anchorA B.anchorPhi B.anchorDx zeros(1,n)]; pa=prediction(za,B.anchorEO*B.rotFreqMeanHz,B,M); qa=isfinite(pa)&isfinite(B.anchorV)&B.FixedMask;
    if any(qa), r.anchor_rmse_mv=sqrt(mean((pa(qa)-B.anchorV(qa)).^2)); r.anchor_max_abs_mv=max(abs(pa(qa)-B.anchorV(qa))); end
end
r.rmse_mv=sqrt(mean((pred(fitMask)-B.V(fitMask)).^2)); r.sse_mv2=sum((pred(fitMask)-B.V(fitMask)).^2); sr=nan(1,numel(M)); for ii=1:numel(M), q=B.sensorIndex==ii & fitMask; sr(ii)=sqrt(mean((pred(q)-B.V(q)).^2)); end
zBase=z; zBase(4:3+n)=0; predBase=prediction(zBase,fFit,B,M); baseOk=isfinite(predBase)&isfinite(B.V);
r.VPredNoGap=predBase;
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
% Multistart agreement is a post-fit identifiability diagnostic only.  It
% never changes the winner or enters the voltage objective.
aa=double([startRows.A]); cc=double([startRows.rmse_mv]);
okA=isfinite(aa)&isfinite(cc);
if any(okA)
    r.multistart_amplitude_sd_mm=std(aa(okA));
    r.multistart_amplitude_range_mm=max(aa(okA))-min(aa(okA));
    bestC=min(cc(okA));
    near=okA & cc <= bestC + max(0.10,0.01*max(bestC,eps));
    r.amplitude_ambiguous=nnz(near)>=2 && (max(aa(near))-min(aa(near)))>0.01;
else
    r.multistart_amplitude_sd_mm=NaN; r.multistart_amplitude_range_mm=NaN; r.amplitude_ambiguous=false;
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
e=p-B.V;
% Keep a fixed residual length. Invalid support is a candidate failure, not
% an opportunity to silently remove difficult observations from the fit.
if isfield(B,'FixedMask')
    e(~good & B.FixedMask)=1e6;
    e(~B.FixedMask)=0;
else
    e(~good)=1e6;
end
if isempty(e), e=1e6; end
end

function p=wrap_phase(p)
p=mod(p+pi,2*pi)-pi;
end

function p=prediction(z,f,B,M)
[p,~]=prediction_components(z,f,B,M);
end

function [p,c]=prediction_components(z,f,B,M)
[p,c]=R5_Evaluate_CanonicalForward(z,f,B,M);
end

function a=foundation_replay_audit(B,M,cfg)
a=struct('available',false,'status','replay_unavailable','rmse_mv',NaN, ...
    'max_abs_mv',NaN,'pass',false,'point_count',0,'sensor_count',0, ...
    'coordinate_pass',false,'time_phase_pass',false);
need=isfield(cfg,'r5') && isfield(cfg.r5,'auditFoundationReplay') && cfg.r5.auditFoundationReplay;
if ~need, return; end
if ~all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx])) || ...
        ~any(isfinite(B.anchorV)), return; end
B0=B; B0.EO=B.anchorEO; B0.useNestedFoundationAnchor=false;
z0=[B.anchorA B.anchorPhi B.anchorDx zeros(1,numel(M))];
p0=R5_Evaluate_CanonicalForward(z0,B.anchorF,B0,M);
[fixedMask,~]=common_support_mask(B,M,cfg);
q=fixedMask & isfinite(B.anchorV) & isfinite(p0);
a.point_count=nnz(q); a.sensor_count=nnz(unique(B.sensorIndex(q)));
a.coordinate_pass=isfield(B,'X') && isfield(B,'sensorIndex') && ...
    numel(B.X)==numel(B.anchorV) && numel(B.sensorIndex)==numel(B.X);
a.time_phase_pass=isfinite(B.anchorF) && B.anchorF>0 && ...
    isfinite(B.anchorEO) && B.anchorEO>0;
if nnz(q)<20 || ~a.coordinate_pass || ~a.time_phase_pass, return; end
d=p0(q)-B.anchorV(q); a.available=true; a.status='replay_pass';
a.rmse_mv=sqrt(mean(d.^2)); a.max_abs_mv=max(abs(d));
tol=1e-6;
if isfield(cfg.r5,'foundationReplayToleranceMv'), tol=double(cfg.r5.foundationReplayToleranceMv); end
a.pass=isfinite(a.rmse_mv) && a.rmse_mv<=tol;
if ~a.pass, a.status='replay_fail'; end
end

function r=empty_row()
r=struct('window_id',NaN,'lap_range',[NaN NaN],'status','unprocessed','quality_class','unprocessed','EO',NaN,'frequency_hz',NaN,'amplitude_mm',NaN,'phase_rad',NaN,'dx_mm',NaN,'delta_gap_mm',NaN,'delta_gap_limit_mm',NaN,'rmse_mv',NaN,'sse_mv2',NaN,'baseline_rmse_mv',NaN,'increment_rms_mv',NaN,'anchor_rmse_mv',NaN,'anchor_max_abs_mv',NaN,'foundation_replay_available',false,'foundation_replay_status','replay_unavailable','foundation_replay_rmse_mv',NaN,'foundation_replay_max_abs_mv',NaN,'foundation_replay_pass',false,'foundation_replay_point_count',0,'foundation_replay_sensor_count',0,'foundation_replay_coordinate_pass',false,'foundation_replay_time_phase_pass',false,'foundation_replay_used_for_acceptance',false,'foundation_replay_role','diagnostic_only','foundation_baseline_EO',NaN,'foundation_baseline_frequency_hz',NaN,'foundation_baseline_amplitude_mm',NaN,'foundation_baseline_dx_mm',NaN,'foundation_baseline_phase_rad',NaN,'sensor_rmse_mv',NaN,'sensor_baseline_rmse_mv',NaN,'sensor_increment_rms_mv',NaN,'exitflag',NaN,'residual_norm',NaN,'hit_boundary',false,'support_pass',false,'support_valid_count',NaN,'support_total_count',NaN,'support_fraction',NaN,'support_info',struct(),'a_bound_hit',false,'dx_bound_hit',false,'dg_bound_hit',false,'multi_start_count',NaN,'multistart_amplitude_sd_mm',NaN,'multistart_amplitude_range_mm',NaN,'amplitude_ambiguous',false,'start_table',struct([]),'eo_screen_table',struct([]),'candidate_table',struct([]),'audit_all_eo_table',struct([]),'audit_all_eo_best_rmse_mv',NaN,'audit_all_eo_best_EO',NaN,'topk_recall_global',NaN,'eo_screen_top1',NaN,'eo_screen_top1_rmse_mv',NaN,'eo_final_vs_screen_top1',false,'eo_consistency_pass',true,'VPred',[],'VPredNoGap',[]);
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
    if isfield(M(i),'templateDomainMm') && all(isfinite(M(i).templateDomainMm))
        d=M(i).templateDomainMm; mask(q)=mask(q)&xlo>=d(1)&xhi<=d(2);
    end
    if isfield(M(i),'responseDomainMm') && all(isfinite(M(i).responseDomainMm)) && isfield(M(i),'coord')
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
