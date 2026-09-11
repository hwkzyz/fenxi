function fit = optimize_experimental_full_waveform(bundle, candidate, model, cfg)
%OPTIMIZE_EXPERIMENTAL_FULL_WAVEFORM Diagnostic full F_s optimization.
% Candidate supplies integer EO, A, phi and dx. The optimized variables are
% A, phi, dx and one dg per sensor; optional frequency refinement is bounded
% to candidate frequency +/- cfg.frequencyRefineHalfWidthHz.
V = bundle.V(:); X = bundle.X(:); Theta = bundle.Theta(:);
sensorIndex = bundle.sensorIndex(:); nSensor = numel(bundle.sensorIds);
W = bundle.W(:); if isempty(W), W = ones(size(V)); end
rotHz = bundle.rotFreqMeanHz;
eo = round(candidate.EO);
freq0 = eo * rotHz;
ampLim = local_cfg(cfg,'amplitudeLimitMm',0.50);
dxLim = local_cfg(cfg,'dxLimitMm',0.35);
dgLim = local_cfg(cfg,'deltaGapLimitMm',0.25);
regW = local_cfg(cfg,'staticRegWeightMv',0);
regScale = max(local_cfg(cfg,'deltaGapPriorScaleMm',0.10),eps);
lb = [0,-pi,-dxLim,-dgLim*ones(1,nSensor)];
ub = [ampLim,pi,dxLim,dgLim*ones(1,nSensor)];
x0 = [candidate.A,candidate.phi,candidate.dx,zeros(1,nSensor)];
x0 = min(max(x0,lb),ub);
obj = @(z, fHz) local_residual(z,fHz,bundle,model,cfg,regW,regScale);
opts = optimoptions('lsqnonlin','Display','off','MaxIterations',120, ...
    'MaxFunctionEvaluations',4000,'FunctionTolerance',1e-10,'StepTolerance',1e-10);
[zInt,~,resInt,exitInt] = lsqnonlin(@(z)obj(z,freq0),x0,lb,ub,opts);
fit = local_pack(zInt,eo,freq0,resInt,exitInt,bundle,model,cfg);
halfW = local_cfg(cfg,'frequencyRefineHalfWidthHz',2.0);
if local_cfg(cfg,'continuousFrequencyRefine',true) && isfinite(freq0)
    lbF = max(300,freq0-halfW); ubF = min(1000,freq0+halfW);
    if ubF > lbF
        z0 = [zInt,freq0]; lb2=[lb,lbF]; ub2=[ub,ubF];
        [z2,~,res2,exit2] = lsqnonlin(@(q)obj(q(1:end-1),q(end)),z0,lb2,ub2,opts);
        fitContinuous = local_pack(z2(1:end-1),eo,z2(end),res2,exit2,bundle,model,cfg);
        fitContinuous.integerPlainRmseMv = fit.plainRmseMv;
        fit = fitContinuous;
    end
end
fit.integerEO = eo; fit.integerFrequencyHz = freq0;
fit.frequencyDeltaFromEOHz = fit.frequencyHz - freq0;
end

function r = local_residual(z,fHz,bundle,model,cfg,regW,regScale)
A=z(1); phi=z(2); dx=z(3); dg=z(4:end);
u=A*sin((fHz/bundle.rotFreqMeanHz)*bundle.Theta(:)+phi);
vPred=nan(size(bundle.V));
for is=1:numel(bundle.sensorIds)
    mask=bundle.sensorIndex(:)==is;
    [v,~]=model(is).evaluate(dg(is),bundle.X(mask)-dx-u(mask));
    vPred(mask)=v(:);
end
valid=isfinite(vPred)&isfinite(bundle.V(:));
r=(vPred(valid)-bundle.V(valid));
if regW>0
    r=[r; sqrt(numel(bundle.V))*regW*(dg(:)/regScale)/sqrt(max(numel(dg),1))]; %#ok<AGROW>
end
if isempty(r), r=1e9; end
end

function fit=local_pack(z,eo,fHz,residual,exitflag,bundle,model,cfg)
A=z(1); phi=z(2); dx=z(3); dg=z(4:end);
u=A*sin((fHz/bundle.rotFreqMeanHz)*bundle.Theta(:)+phi);
vPred=nan(size(bundle.V));
for is=1:numel(bundle.sensorIds)
    mask=bundle.sensorIndex(:)==is;
    [v,~]=model(is).evaluate(dg(is),bundle.X(mask)-dx-u(mask)); vPred(mask)=v(:);
end
valid=isfinite(vPred)&isfinite(bundle.V(:)); plain=sqrt(mean((vPred(valid)-bundle.V(valid)).^2));
fit=struct('EO',eo,'freqHz',fHz,'frequencyHz',fHz,'amplitudeMm',A,'phaseRad',phi, ...
    'dxMm',dx,'deltaGapMm',dg(:),'plainRmseMv',plain,'weightedRmseMv',plain, ...
    'residualNorm',norm(residual), 'exitflag',exitflag,'VPred',vPred,'uMm',u, ...
    'hitAmplitudeBoundary',abs(A-local_cfg(cfg,'amplitudeLimitMm',0.50))<1e-8, ...
    'hitDxBoundary',abs(abs(dx)-local_cfg(cfg,'dxLimitMm',0.35))<1e-8, ...
    'hitGapBoundary',any(abs(abs(dg)-local_cfg(cfg,'deltaGapLimitMm',0.25))<1e-8));
end

function v=local_cfg(cfg,name,default)
if isfield(cfg,name)&&~isempty(cfg.(name))&&isfinite(cfg.(name)),v=cfg.(name);else,v=default;end
end
