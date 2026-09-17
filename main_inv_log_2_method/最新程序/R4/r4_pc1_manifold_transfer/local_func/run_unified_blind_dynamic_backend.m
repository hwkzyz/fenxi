function R = run_unified_blind_dynamic_backend(bundle, response, cfg)
%RUN_UNIFIED_BLIND_DYNAMIC_BACKEND Unified blind gap/vibration backend.
% The response adapter is the only R2/R3/R4-specific object.  The estimator
% never receives truth parameters.  VP is used for candidate generation;
% final ranking is data-only full-wave voltage SSE.

if ~isfield(cfg,'gapGrid'), cfg.gapGrid=cfg.gapBounds(1):.05:cfg.gapBounds(2); end
if ~isfield(cfg,'eoGrid'), cfg.eoGrid=1:30; end
if ~isfield(cfg,'topK'), cfg.topK=3; end
if ~isfield(cfg,'minGradient'), cfg.minGradient=0; end
if ~isfield(cfg,'minSamples'), cfg.minSamples=20; end
if ~isfield(cfg,'dxBound'), cfg.dxBound=.2; end
if ~isfield(cfg,'ampCoeffBound'), cfg.ampCoeffBound=1.5; end
if ~isfield(cfg,'deltaFBound'), cfg.deltaFBound=2; end
if ~isfield(cfg,'maxIter'), cfg.maxIter=300; end
if ~isfield(cfg,'refineAmplitudeScales'), cfg.refineAmplitudeScales=[1 .5 2]; end
if ~isfield(cfg,'refinePhaseOffsets'), cfg.refinePhaseOffsets=[0 pi/2 pi -pi/2]; end
if ~isfield(cfg,'refineAmplitudeGrid'), cfg.refineAmplitudeGrid=[]; end
% VP can substantially underestimate A when the gap increment absorbs part
% of the waveform.  Always include a modest, truth-independent amplitude
% grid unless the caller supplied one explicitly; this mirrors the
% experiment's vibration-only prefit followed by full-wave refinement.
if isempty(cfg.refineAmplitudeGrid), cfg.refineAmplitudeGrid=0:.05:min(.50,cfg.ampCoeffBound); end
if ~isfield(cfg,'ambiguityMargin'), cfg.ambiguityMargin=.01; end
if ~isfield(cfg,'profileTiming'), cfg.profileTiming=false; end
tAll = tic;
required = {'t','x','V','theta','sensorId','turnId'};
for k=1:numel(required), assert(isfield(bundle,required{k}), 'Missing bundle.%s.',required{k}); end
assert(isfield(response,'evalF') && isfield(response,'evalFx'), 'Response must provide evalF/evalFx.');
t=bundle.t(:); x=bundle.x(:); V=bundle.V(:); th=bundle.theta(:);
sid=bundle.sensorId(:); turn=bundle.turnId(:); %#ok<NASGU>
if isfield(bundle,'Wvp') && ~isempty(bundle.Wvp), Wvp=bundle.Wvp(:); else, Wvp=ones(size(V)); end
if isfield(bundle,'Wfull') && ~isempty(bundle.Wfull), Wfull=bundle.Wfull(:); else, Wfull=ones(size(V)); end
assert(numel(t)==numel(V) && numel(x)==numel(V) && numel(th)==numel(V));
assert(all(isfinite(Wvp)) && all(Wvp>=0) && any(Wvp>0));
assert(all(isfinite(Wfull)) && all(Wfull>=0) && any(Wfull>0));
% Enforce the fixed-path query support before VP or nonlinear fitting.  The
% simulator leaves samples at baseline when the displaced coordinate exits
% the trusted domain; treating those as voltage residuals biases amplitude.

gapGrid=cfg.gapGrid(:); eoGrid=cfg.eoGrid(:); topK=min(cfg.topK,numel(gapGrid)*numel(eoGrid));
rotHz=cfg.rotHz; tCenter=mean(t);
rows=repmat(candidate_row(),0,1);
tScan=tic;
for ig=1:numel(gapGrid)
    g=gapGrid(ig); F0=call_response(response,'evalF',g,x,sid); Fx=call_response(response,'evalFx',g,x,sid);
    valid=isfinite(F0)&isfinite(Fx)&isfinite(V)&(abs(Fx)>=cfg.minGradient);
    w=Wvp; w(~valid)=0; q=-(V-F0); q(valid)=q(valid)./Fx(valid); q(~valid)=0;
    for ie=1:numel(eoGrid)
        eo=eoGrid(ie); s=sin(eo*th); c=cos(eo*th); X=[ones(size(th)),s,c];
        ww=sqrt(w); A=X.*ww; b=q.*ww; keep=w>0;
        if nnz(keep)<cfg.minSamples || rank(A(keep,:))<3, continue; end
        beta=A(keep,:)\b(keep); qhat=X*beta; vpRmse=sqrt(sum(w.*(q-qhat).^2)/max(sum(w),eps));
        a=beta(2); bb=beta(3); seed=[g,beta(1),a,bb,0];
        % Keep the scan genuinely variable-projection: the linear
        % displacement RMSE ranks the full grid.  Expensive voltage
        % evaluation is deferred to Top-K only (as in the experimental VP
        % backend), instead of evaluating every gap/EO pair.
        rows(end+1)=pack_candidate(g,eo,beta,vpRmse,NaN); %#ok<AGROW>
    end
end
if isempty(rows), R=struct('status',"failed",'candidates',rows); return; end
tScanElapsed=toc(tScan);
% The gradient model is a seed generator, not the final physical model.
% The VP coefficients are only a frequency/EO localization device.  The
% experimental procedure does not rank candidates by the VP amplitude
% replay, because that replay can inherit a biased low amplitude.  Retain
% candidates by the VP displacement score, then re-estimate gap, dx and the
% sine/cosine coefficients from the complete voltage waveform below.
[~,ordSeed]=sort([rows.vpDispRmseMm],'ascend');
rows=rows(ordSeed(1:min(topK,numel(ordSeed))));
for k=1:numel(rows)
    q=rows(k); seed=[q.gap,q.d,q.a,q.b,0];
    rSeed=full_residual_theta(seed,q.eo,t,th,x,V,sid,response,Wfull,tCenter);
    rows(k).seedVoltageRmseV=sqrt(mean(rSeed.^2));
    rows(k).seedVoltageRmseMv=1e3*rows(k).seedVoltageRmseV;
end
refined=repmat(refined_row(),0,1);
tRefine=tic;
for k=1:numel(rows)
    q=rows(k); eo=q.eo; p0=[q.gap,q.d,q.a,q.b,0];
    lb=[cfg.gapBounds(1),-cfg.dxBound,-cfg.ampCoeffBound,-cfg.ampCoeffBound,-cfg.deltaFBound];
    ub=[cfg.gapBounds(2), cfg.dxBound, cfg.ampCoeffBound, cfg.ampCoeffBound, cfg.deltaFBound];
    % The VP coefficients are a local linear seed.  A single nonlinear
    % start can remain at a low-amplitude basin when the response surface is
    % strongly nonlinear.  Use deterministic, truth-independent amplitude
    % rescalings around the VP seed and retain the best voltage fit.
    scales=cfg.refineAmplitudeScales(:).';
    amp0=hypot(p0(3),p0(4));
    if isempty(cfg.refineAmplitudeGrid)
        ampStarts=amp0*scales;
    else
        ampStarts=unique([amp0*scales, cfg.refineAmplitudeGrid(:).']);
    end
    % Stage 1: keep the VP frequency candidate fixed and re-estimate the
    % physical clearance, offset, and sine/cosine coefficients together.
    % This removes the VP displacement amplitude from the nonlinear start;
    % the coefficients are re-fit against the actual voltage waveform.
    lb4=lb(1:4); ub4=ub(1:4); bestFixedSse=inf; pFixed=[];
    % Refit the physical amplitude from several deterministic magnitude and
    % phase starts.  This is the important distinction from using the VP
    % amplitude as an initial-value constraint: each start is free to move
    % to the voltage-optimal A/phi basin.
    amp0=hypot(p0(3),p0(4));
    ampStarts=amp0*scales;
    if ~isempty(cfg.refineAmplitudeGrid), ampStarts=[ampStarts,cfg.refineAmplitudeGrid(:).']; end
    ampStarts=unique(max(0,min(cfg.ampCoeffBound,ampStarts)));
    fixedStarts=zeros(numel(ampStarts)*max(1,numel(cfg.refinePhaseOffsets)),4);
    is=0;
    for ia=1:numel(ampStarts)
        for ip=1:max(1,numel(cfg.refinePhaseOffsets))
            is=is+1;
            ph=atan2(p0(4),p0(3))+cfg.refinePhaseOffsets(min(ip,numel(cfg.refinePhaseOffsets)));
            fixedStarts(is,:)=[p0(1),p0(2),ampStarts(ia)*cos(ph),ampStarts(ia)*sin(ph)];
        end
    end
    % Experimental R5 performs a vibration-only nonlinear prefit before
    % introducing dynamic gap increments.  Reproduce that staging here:
    % first refine [dx,a,b] at the VP gap, then use the result as an
    % additional full-wave start.  This prevents dg from explaining the
    % vibration amplitude before the vibration state has been established.
    q0=p0(2:4);
    vibFun=@(q) fixed_eo_residual([p0(1),q(:).'],eo,t,th,x,V,sid,response,Wfull,tCenter);
    [qV,~]=solve_lsq_bounded(vibFun,q0,lb4(2:4),ub4(2:4),cfg.maxIter,1e-9,1e-9);
    fixedStarts(end+1,:)=[p0(1),qV(:).'];
    bestFixedSse=inf; pFixed=[];
    for ist=1:size(fixedStarts,1)
        pStart=min(max(fixedStarts(ist,:),lb4),ub4);
        [p4Try,~]=solve_lsq_bounded(@(z)fixed_eo_residual(z,eo,t,th,x,V,sid,response,Wfull,tCenter),...
            pStart,lb4,ub4,cfg.maxIter,1e-9,1e-9);
        r4=fixed_eo_residual(p4Try,eo,t,th,x,V,sid,response,Wfull,tCenter);
        sse4=sum(r4.^2);
        if isfinite(sse4) && sse4<bestFixedSse, bestFixedSse=sse4; pFixed=p4Try; end
    end
    if isempty(pFixed), pFixed=p0(1:4); end
    pFixed=min(max(pFixed,lb4),ub4);
    % Stage 2: retain the VP EO but allow a small frequency correction.
    bestSse=inf; p=[]; info=[];
    pFixed5=[pFixed,0];
    releaseStarts=[pFixed5; p0];
    for ir=1:size(releaseStarts,1)
        pStart=releaseStarts(ir,:);
        [pTry,infoTry]=solve_lsq_bounded(@(z)full_residual_theta(z,eo,t,th,x,V,sid,response,Wfull,tCenter),pStart,lb,ub,cfg.maxIter,1e-9,1e-9);
        rTry=full_residual_theta(pTry,eo,t,th,x,V,sid,response,Wfull,tCenter);
        sseTry=sum(rTry.^2);
        if isfinite(sseTry) && sseTry<bestSse
            bestSse=sseTry; p=pTry; info=infoTry;
        end
    end
    if isempty(p), p=p0; info=struct('status',"no_finite_refinement"); end
    rr=full_residual_theta(p,eo,t,th,x,V,sid,response,Wfull,tCenter); sse=sum(rr.^2); A=hypot(p(3),p(4)); phi=atan2(p(4),p(3));
    refined(end+1)=struct('gap',p(1),'dx',p(2),'A',A,'phi',phi,'deltaF',p(5),'eo',eo,'fullWaveSseV2',sse,'fullWaveRmseV',sqrt(mean(rr.^2)), ...
        'fullWaveSseMv2',1e6*sse,'fullWaveRmseMv',1e3*sqrt(mean(rr.^2)), ...
        'vpDispRmseMm',q.vpDispRmseMm,'seedVoltageRmseV',q.seedVoltageRmseV, ...
        'seedVoltageRmseMv',1e3*q.seedVoltageRmseV,'solveInfo',info); %#ok<AGROW>
end
[~,ord]=sort([refined.fullWaveSseMv2]); refined=refined(ord);
tRefineElapsed=toc(tRefine);
margin=NaN; if numel(refined)>1, margin=(refined(2).fullWaveSseMv2-refined(1).fullWaveSseMv2)/max(refined(1).fullWaveSseMv2,eps); end
best=refined(1); status="identified"; if ~isfinite(margin)||margin<cfg.ambiguityMargin, status="ambiguous_joint"; end
R=struct('status',status,'best',best,'secondBestMargin',margin,'candidates',rows,'refined',refined,'supportCount',nnz(Wfull>0),'method',"gradient_VP_seed_then_full_wave");
R.timing=struct('scan_s',tScanElapsed,'refine_s',tRefineElapsed,'total_s',toc(tAll), ...
    'nGap',numel(gapGrid),'nEO',numel(eoGrid),'nScanCandidates',numel(gapGrid)*numel(eoGrid), ...
    'nTopK',numel(rows),'nSeedReplay',numel(ordSeed));
end

function r=full_residual_theta(p,eo,t,theta,x,V,sid,response,W,tCenter)
phase=eo*theta+2*pi*p(5)*(t-tCenter); u=p(3)*sin(phase)+p(4)*cos(phase);
pred=call_response(response,'evalF',p(1),x-p(2)-u,sid); r=sqrt(W).*(V-pred);
r(~isfinite(r))=10*max(std(V),1e-3);
end

function r=fixed_eo_residual(p,eo,t,theta,x,V,sid,response,W,tCenter)
% Residual with the VP-selected EO held fixed; p=[gap,dx,a,b].
phase=eo*theta; u=p(3)*sin(phase)+p(4)*cos(phase);
pred=call_response(response,'evalF',p(1),x-p(2)-u,sid); r=sqrt(W).*(V-pred);
r(~isfinite(r))=10*max(std(V),1e-3);
end

function y=call_response(response,name,g,x,sid)
% Sensor-conditioned response is optional for legacy single-sensor adapters.
f=response.(name);
try
    y=f(g,x,sid);
catch
    y=f(g,x);
end
end

function r=pack_candidate(g,eo,beta,vp,seed)
r=struct('gap',g,'eo',eo,'d',beta(1),'a',beta(2),'b',beta(3),'vpDispRmseMm',vp,'seedVoltageRmseV',seed,'seedVoltageRmseMv',1e3*seed);
end
function r=candidate_row(), r=struct('gap',NaN,'eo',NaN,'d',NaN,'a',NaN,'b',NaN,'vpDispRmseMm',NaN,'seedVoltageRmseV',NaN,'seedVoltageRmseMv',NaN); end
function r=refined_row(), r=struct('gap',NaN,'dx',NaN,'A',NaN,'phi',NaN,'deltaF',NaN,'eo',NaN,'fullWaveSseV2',NaN,'fullWaveRmseV',NaN,'fullWaveSseMv2',NaN,'fullWaveRmseMv',NaN,'vpDispRmseMm',NaN,'seedVoltageRmseV',NaN,'seedVoltageRmseMv',NaN,'solveInfo',[]); end
