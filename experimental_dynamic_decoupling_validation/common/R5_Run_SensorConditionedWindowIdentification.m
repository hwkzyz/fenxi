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
% Do not silently turn a gap-aware run into a no-gap fit when a builder
% forgets to declare the sensor role.  This was the cause of a misleading
% 20250527 result with all delta_gap values fixed at zero.
if isfield(cfg.case,'gapSensors')
    declaredGap=arrayfun(@(m) isfield(m,'isGapSensor') && m.isGapSensor, models);
    expectedGap=ismember(ids,double(cfg.case.gapSensors));
    assert(isequal(declaredGap,expectedGap), ...
        'R5:SensorRoleDeclaration', ...
        'Builder gap-role flags do not match cfg.case.gapSensors.');
end
fit_bounds(cfg);
Q=load(foundationFile,'Result'); R=Q.Result; W=R.WindowResult; rows=repmat(empty_row(),numel(W),1);
for iw=1:numel(W)
    B=W(iw).CoreBundlePreview; bundle=make_bundle(B,ids,cfg,W(iw));
    rows(iw).window_id=W(iw).window_id; rows(iw).lap_range=W(iw).lap_range;
    if isempty(W(iw).SeedTable)
        rows(iw).status='no_seed'; continue;
    end
expectedEO=NaN;
    if isfield(R,'Trend') && istable(R.Trend) && ismember('EO_id',R.Trend.Properties.VariableNames)
        kEO=find(double(R.Trend.window_id)==double(W(iw).window_id),1);
        if ~isempty(kEO), expectedEO=double(R.Trend.EO_id(kEO)); end
    end
fit=fit_window(bundle,models,W(iw).SeedTable,expectedEO,cfg);
    fit.window_id=rows(iw).window_id; fit.lap_range=rows(iw).lap_range; rows(iw)=fit;
end

out=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_RESULT_V1','dataset',cfg.dataset,...
    'targetBlade',cfg.case.targetBlade,'sensorIds',ids,'rows',rows,...
    'gapSensorIds',double(cfg.case.gapSensors(:).'),...
    'staticGapSlopeInDynamic',logical(isfield(cfg,'r5') && isfield(cfg.r5,'useStaticGapSlopeInDynamic') && cfg.r5.useStaticGapSlopeInDynamic),...
    'referenceEO',double(getfield_or_nan(cfg,'r5','referenceEO')),...
    'referenceEOTolerance',double(getfield_or_nan(cfg,'r5','referenceEOTolerance')),...
    'lockReferenceEO',logical(isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO),...
    'fitBounds',fit_bounds(cfg),...
    'forbiddenDynamicParameters',{{'deltaMu','deltaTau'}},...
    'useNestedFoundationAnchor',logical(isfield(cfg,'r5') && isfield(cfg.r5,'useNestedFoundationAnchor') && cfg.r5.useNestedFoundationAnchor),...
    'foundationFile',foundationSourceFile,'stagedFoundationFile',foundationFile, ...
    'sidecarFile',sidecarSourceFile,'stagedSidecarFile',sidecarFile);
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'out','-v7.3');
end

function v=getfield_or_nan(s,group,name)
v=NaN;
if isfield(s,group) && isstruct(s.(group)) && isfield(s.(group),name), v=s.(group).(name); end
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
B.anchorEO=NaN; B.anchorA=NaN; B.anchorPhi=NaN; B.anchorDx=NaN;
useNested=isfield(cfg,'r5') && isfield(cfg.r5,'useNestedFoundationAnchor') && cfg.r5.useNestedFoundationAnchor;
if useNested
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

function r=fit_window(B,M,T,expectedEO,cfg)
assert(~isempty(T),'R5:NoSeed');
% The legacy route retains several EO/phase candidates before the final
% continuous fit.  A single seed can trap lsqnonlin in a poor local basin,
% especially when gap and vibration terms partially compensate each other.
rmseSeed=double([T.plain_voltage_rmse]); rmseSeed(~isfinite(rmseSeed))=Inf;
% Keep the complete EO candidate set in play.  The Foundation trend is a
% useful diagnostic/initialization hint, but it is not an observation-level
% constraint: locking the dynamic fit to expectedEO can retain a wrong EO
% whenever the transferred gap increment changes the voltage ranking.  This
% was a major difference from the legacy full-wave route, which compared
% several EO candidates after the gap-aware refinement.
ordAll=1:numel(T);
[~,ord0]=sort(rmseSeed(ordAll),'ascend'); ordAll=ordAll(ord0);
nKeep=min(candidate_keep_count(cfg),numel(ordAll));
ord=ordAll(1:nKeep);
if isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO
    eoAll=double([T.EO]); match=find(eoAll==double(cfg.r5.referenceEO));
    assert(~isempty(match),'R5:ReferenceEOMissing','Locked reference EO is absent from the Foundation seeds.');
    ord=match(:).';
end
eoLocked=isfield(cfg,'r5') && isfield(cfg.r5,'lockReferenceEO') && cfg.r5.lockReferenceEO;
% expectedEO is retained only as a provenance field by the caller.  It must
% not be injected into the candidate set: final EO selection is data-driven.
best=[]; bestRmse=Inf; candidates=repmat(empty_row(),0,1);
for kk=1:numel(ord)
    cand=fit_window_from_seed(B,M,T(ord(kk)),cfg);
    if isfinite(cand.rmse_mv), candidates(end+1)=cand; end %#ok<AGROW>
    if isfinite(cand.rmse_mv) && cand.rmse_mv<bestRmse
        best=cand; bestRmse=cand.rmse_mv;
    end
end
if isempty(best), best=fit_window_from_seed(B,M,T(ord(1)),cfg); end
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

function r=fit_window_from_seed(B,M,s,cfg)
n=numel(M); B.EO=double(s.EO); f0=double(s.EO)*B.rotFreqMeanHz; [limA,limDx,limDg]=fit_bounds_from_cfg();
% Keep a parameter slot for every analysed channel so the common vibration
% state is constrained by direct/template channels as well.  Direct-only
% channels have no clearance degree of freedom; their dg slot is fixed at 0.
isGap=arrayfun(@(m) isfield(m,'isGapSensor') && m.isGapSensor, M);
dgLimit=local_gap_projection_limits(B,M,double(s.A),double(s.phi),double(s.dx_c),f0,limDg,cfg);
dgLimit=dgLimit(:).'; isGap=logical(isGap(:).');
dgLo=-dgLimit.*double(isGap); dgHi=dgLimit.*double(isGap);
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
p0=prediction(z0,f0,B,M); epsDg=min(1e-3,0.02*limDg);
for ii=1:n
    if ~isGap(ii), continue; end
    ze=z0; ze(3+ii)=epsDg; pe=prediction(ze,f0,B,M);
    q=(B.sensorIndex==ii) & isfinite(p0) & isfinite(pe) & isfinite(B.V);
    if nnz(q)>=4
        d=(pe(q)-p0(q))/epsDg; rr=B.V(q)-p0(q); den=sum(d.^2);
        if isfinite(den) && den>eps
            z0(3+ii)=sum(d.*rr)/den;
        end
    end
end
z0=min(max(z0,lb),ub);
fun=@(z)resid(z,f0+z(end),B,M); opts=optimoptions('lsqnonlin','Display','off','MaxIterations',120,'MaxFunctionEvaluations',3000);
[z,~,res,exitflag]=lsqnonlin(fun,z0,lb,ub,opts); pred=prediction(z,f0,B,M); ok=isfinite(pred)&isfinite(B.V);
fFit=f0+z(end); pred=prediction(z,fFit,B,M); ok=isfinite(pred)&isfinite(B.V);
r=empty_row(); r.status='pass'; r.EO=double(s.EO); r.frequency_hz=fFit; r.amplitude_mm=z(1); r.phase_rad=z(2); r.dx_mm=z(3); r.delta_gap_mm=z(4:3+n).';
r.anchor_rmse_mv=NaN; r.anchor_max_abs_mv=NaN;
if B.useNestedFoundationAnchor && all(isfinite([B.anchorEO B.anchorA B.anchorPhi B.anchorDx]))
    za=[B.anchorA B.anchorPhi B.anchorDx zeros(1,n)]; pa=prediction(za,B.anchorEO*B.rotFreqMeanHz,B,M); qa=isfinite(pa)&isfinite(B.anchorV);
    if any(qa), r.anchor_rmse_mv=sqrt(mean((pa(qa)-B.anchorV(qa)).^2)); r.anchor_max_abs_mv=max(abs(pa(qa)-B.anchorV(qa))); end
end
r.rmse_mv=sqrt(mean((pred(ok)-B.V(ok)).^2)); sr=nan(1,numel(M)); for ii=1:numel(M), q=B.sensorIndex==ii & ok; sr(ii)=sqrt(mean((pred(q)-B.V(q)).^2)); end
zBase=z; zBase(4:3+n)=0; predBase=prediction(zBase,fFit,B,M); baseOk=isfinite(predBase)&isfinite(B.V);
r.baseline_rmse_mv=sqrt(mean((predBase(baseOk)-B.V(baseOk)).^2));
r.increment_rms_mv=sqrt(mean((pred(ok)-predBase(ok)).^2));
r.sensor_baseline_rmse_mv=nan(1,numel(M)); r.sensor_increment_rms_mv=nan(1,numel(M));
for ii=1:numel(M)
    q=B.sensorIndex==ii & baseOk; if any(q), r.sensor_baseline_rmse_mv(ii)=sqrt(mean((predBase(q)-B.V(q)).^2)); end
    q=B.sensorIndex==ii & ok; if any(q), r.sensor_increment_rms_mv(ii)=sqrt(mean((pred(q)-predBase(q)).^2)); end
end
r.sensor_rmse_mv=sr; r.exitflag=exitflag; r.residual_norm=norm(res);
% Treat numerical solutions within the optimizer's practical resolution of
% a bound as boundary hits; exact equality alone misses saturated fits.
tolBound=1e-5;
r.delta_gap_limit_mm=dgLimit(:);
r.hit_boundary=any(isGap & abs(abs(z(4:3+n))-dgLimit)<tolBound)| ...
    abs(abs(z(3))-limDx)<tolBound|abs(z(1)-limA)<tolBound;
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
p=prediction(z,f,B,M); ok=isfinite(p)&isfinite(B.V); e=p(ok)-B.V(ok); if isempty(e), e=1e9; end
end

function p=prediction(z,f,B,M)
[p,~]=prediction_components(z,f,B,M);
end

function [p,c]=prediction_components(z,f,B,M)
% Preserve the old gap_only phase law: integer EO multiplies the measured
% accumulated rotor angle, while the continuous frequency refinement is a
% time-domain detuning term.  Scaling theta by f/rotFreqMeanHz is not
% equivalent when the speed varies within a window and can bias A/dx.
eo=double(B.EO); if ~isfinite(eo), eo=f/B.rotFreqMeanHz; end
phase=eo*B.Theta+z(2)+2*pi*(f-eo*B.rotFreqMeanHz)*(B.T-mean(B.T));
u=z(1)*sin(phase); p=nan(size(B.V));
c=struct('foundation',nan(size(B.V)),'noGapCurrent',nan(size(B.V)),...
    'noGapBase',nan(size(B.V)),'gapIncrement',nan(size(B.V)),...
    'xCurrent',B.X-z(3)-u,'xBase',nan(size(B.V)),'uCurrent',u);
% The nested anchor is a frozen Foundation waveform.  The R5 model adds
% only the change caused by the current vibration state and dg_s, matching
% the old V1 contract while keeping dynamic tilt/registration frozen.
for i=1:numel(M)
    q=B.sensorIndex==i; [p(q),d]=M(i).evaluate(z(3+i),c.xCurrent(q));
    c.noGapCurrent(q)=d.noGapMv; c.gapIncrement(q)=d.gapIncrementMv;
end
if B.useNestedFoundationAnchor
    % Anchor correction is applied outside the sensor loop so direct and
    % gap channels share the same frozen baseline.  The zero-dg operator
    % supplies the R5 observation at the current state; the anchor itself
    % is retained as the absolute baseline.
    q=isfinite(B.anchorV) & isfinite(p);
    p(q)=B.anchorV(q) + p(q) - local_zero_gap_prediction(z,f,B,M,q);
end
end

function p0=local_zero_gap_prediction(z,f,B,M,qKeep)
z0=z; z0(4:end)=0; z0(1:3)=z(1:3);
if B.useNestedFoundationAnchor && all(isfinite([B.anchorA B.anchorPhi B.anchorDx]))
    z0(1:3)=[B.anchorA B.anchorPhi B.anchorDx];
end
eo=double(B.EO); if ~isfinite(eo), eo=f/B.rotFreqMeanHz; end
phase0=eo*B.Theta+z0(2)+2*pi*(f-eo*B.rotFreqMeanHz)*(B.T-mean(B.T));
u0=z0(1)*sin(phase0); p0=nan(size(B.V));
for i=1:numel(M)
    q=(B.sensorIndex==i)&qKeep; [p0(q),~]=M(i).evaluate(0,B.X(q)-z0(3)-u0(q));
end
end

function r=empty_row()
r=struct('window_id',NaN,'lap_range',[NaN NaN],'status','unprocessed','EO',NaN,'frequency_hz',NaN,'amplitude_mm',NaN,'phase_rad',NaN,'dx_mm',NaN,'delta_gap_mm',NaN,'delta_gap_limit_mm',NaN,'rmse_mv',NaN,'baseline_rmse_mv',NaN,'increment_rms_mv',NaN,'anchor_rmse_mv',NaN,'anchor_max_abs_mv',NaN,'sensor_rmse_mv',NaN,'sensor_baseline_rmse_mv',NaN,'sensor_increment_rms_mv',NaN,'exitflag',NaN,'residual_norm',NaN,'hit_boundary',false);
end

function limits=local_gap_projection_limits(B,M,A,phi,dx,f0,globalLimit,cfg)
limits=repmat(globalLimit,numel(M),1);
if ~(isfield(cfg,'r5') && isfield(cfg.r5,'useLocalGapProjection') && cfg.r5.useLocalGapProjection), return; end
alpha=cfg_num(cfg.r5,'gapProjectionAlpha',2.5); minLim=cfg_num(cfg.r5,'gapProjectionMinLimitMm',0.02);
hG=cfg_num(cfg.r5,'gapProjectionStepMm',1e-3); floorS=cfg_num(cfg.r5,'gapProjectionSensitivityFloorMvPerMm',0.05);
n=numel(M); z=[A phi dx zeros(1,n)]; p0=prediction(z,f0,B,M); valid=isfinite(p0)&isfinite(B.V);
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
    if isfinite(raw), limits(is)=min(globalLimit,max(minLim,raw)); end
end
end

function v=cfg_num(s,name,d)
v=d; if isfield(s,name) && isfinite(s.(name)), v=double(s.(name)); end
end
