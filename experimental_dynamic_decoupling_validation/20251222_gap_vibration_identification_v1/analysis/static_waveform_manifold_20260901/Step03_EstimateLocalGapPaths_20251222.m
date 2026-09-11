%% 03: Estimate each blade's local gap path from the B2 reference surface
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir,'results');
F = load(fullfile(outDir,'StaticWaveformFamily_20251222.mat'),'StaticWaveformFamily'); F=F.StaticWaveformFamily;
U = load(fullfile(outDir,'BladeWaveformSurface_20251222.mat'),'BladeWaveformSurface'); R=U.BladeWaveformSurface;
x=F.xiAxisMm(:); gaps=F.nominalGapMm(:).'; blades=F.bladeIds(:).'; nX=numel(x); nG=numel(gaps); nB=numel(blades); g0=R.g0Mm;
Y=R.observedWaveformsReferenceMv; C=R.coeffReference; b2=find(blades==2,1); if isempty(b2), error('B2 is required.'); end
cfg.deltaBoundsMm=[-0.30 0.30]; cfg.qBoundsMmPerMm=[-0.08 0.08]; cfg.lambdaQ=1e3;
out=nan(nB,8); pathX=nan(nX,nB); recon=nan(size(Y));
for ib=1:nB
    if ib==b2, out(ib,:)=[2 0 0 0 Inf 0 1 0]; pathX(:,ib)=0; recon(:,ib,:)=Y(:,ib,:); continue; end
    target=squeeze(Y(:,ib,:)); starts=[0 0; -.05 0; .05 0; 0 -.02; 0 .02]; best=[Inf 0 0];
    for is=1:size(starts,1)
        z=fminsearch(@(p)pathObjective(p,target,C(:,:,b2),x,gaps,g0,cfg),starts(is,:),optimset('Display','off','MaxIter',500,'MaxFunEvals',1000));
        z(1)=min(max(z(1),cfg.deltaBoundsMm(1)),cfg.deltaBoundsMm(2)); z(2)=min(max(z(2),cfg.qBoundsMmPerMm(1)),cfg.qBoundsMmPerMm(2));
        v=pathObjective(z,target,C(:,:,b2),x,gaps,g0,cfg); if v<best(1), best=[v z]; end
    end
    delta=best(2); q=best(3); path=delta+q*(x-median(x)); pathX(:,ib)=path;
    pred=evalSurface(C(:,:,b2),x,gaps+path,g0); recon(:,ib,:)=pred;
    valid=isfinite(target)&isfinite(pred); rmse=sqrt(mean((target(valid)-pred(valid)).^2));
    support=min(min(gaps)+min(path), max(gaps)-max(path)); atBound=(abs(delta)>=max(abs(cfg.deltaBoundsMm))-1e-9)||(abs(q)>=max(abs(cfg.qBoundsMmPerMm))-1e-9);
    conf=double(~atBound && support>=0); out(ib,:)=[blades(ib),delta,q,rmse,support,atBound,conf,best(1)];
end
BladeGapPathModel=struct('method','B2_surface_constrained_local_path','referenceBlade',2,'bladeIds',blades,'nominalGapMm',gaps,'xAxisMm',x,'deltaBoundsMm',cfg.deltaBoundsMm,'qBoundsMmPerMm',cfg.qBoundsMmPerMm,'deltaGMm',out(:,2),'qGMmPerMm',out(:,3),'gapPathOffsetMm',pathX,'reconstructionRMSEMv',out(:,4),'supportMarginMm',out(:,5),'atParameterBoundary',logical(out(:,6)),'confidenceGate',logical(out(:,7)),'reconstructedWaveformsMv',recon);
save(fullfile(outDir,'BladeGapPathModel_20251222.mat'),'BladeGapPathModel','-v7.3');
writematrix(out,fullfile(outDir,'03_BladeGapPathSummary_20251222.csv'));
fprintf('Path estimation saved. Accepted non-B2 blades: %d/%d.\n',nnz(out(:,7))-1,nB-1);

function val=pathObjective(p,target,C,x,gaps,g0,cfg)
path=p(1)+p(2)*(x-median(x)); G=gaps+path; if min(G(:))<=0, val=1e12; return; end
pred=evalSurface(C,x,G,g0); ok=isfinite(target)&isfinite(pred); if nnz(ok)<20, val=1e12; return; end
r=target(ok)-pred(ok); val=mean(r.^2)+cfg.lambdaQ*p(2)^2;
end
function Y=evalSurface(C,x,G,g0)
nX=numel(x); nG=size(G,2); Y=nan(nX,nG);
for ix=1:nX
    c=C(ix,:); g=G(ix,:); Y(ix,:)=c(1)+c(2)./g+c(3).*log(g./g0);
end
end
