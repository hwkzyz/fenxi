%% Main06A: estimate blade-relative static state on Blade-2 direct support
% Does not modify or rebuild the physical response library.  The reference
% library is Blade 2 at the 16 nominal setpoints; other blades only provide
% relative offset/gradient diagnostics on a fixed direct-overlap subset.
clear; clc;
thisDir=fileparts(mfilename('fullpath'));
step05=fullfile(thisDir,'outputs','Step05_Response_Surface_20250527.mat');
if ~isfile(step05), error('Run revised Main05 first: %s',step05); end
S=load(step05,'responseSurface'); R=S.responseSurface;
outDir=fullfile(thisDir,'outputs');
bladeIds=R.bladeIds(:).'; refId=2; refIdx=find(bladeIds==refId,1);
x=R.xGrid(:); gaps=R.trueGapMm(:); Y=R.waveforms;
gLo=min(gaps); gHi=max(gaps); x0=0;
% Fixed direct-overlap anchor set: selected before fitting, independent of fit.
anchorPool=gaps<=1.70 & gaps>=0.90;
anchorMask=false(size(gaps)); anchorIdx=find(anchorPool); anchorMask(anchorIdx(1:2:end))=true;
anchorGaps=gaps(anchorMask);
fitMask=R.effectiveWindow(:) & isfinite(x);
fprintf('Direct reference: Blade %d, gaps %.3f..%.3f mm; anchor states=%d\n',refId,gLo,gHi,nnz(anchorMask));
% Fitted response coefficients are those learned from the fixed Blade-2 library.
coef=R.coeff; gBasis0=R.g0Mm;
states=repmat(struct('bladeId',NaN,'deltaGapMm',NaN,'muGapPerXMm',NaN,'fitRmseMv',NaN,'holdoutRmseMv',NaN,'valid',false,'anchorGapsMm',[]),1,numel(bladeIds));
for ib=1:numel(bladeIds)
    b=bladeIds(ib); states(ib).bladeId=b; states(ib).anchorGapsMm=anchorGaps(:).';
    if b==refId, states(ib).deltaGapMm=0; states(ib).muGapPerXMm=0; states(ib).fitRmseMv=0; states(ib).holdoutRmseMv=0; states(ib).valid=true; continue; end
    p0=[0,0]; opts=optimset('Display','off','MaxIter',2000,'MaxFunEvals',5000,'TolX',1e-7,'TolFun',1e-6);
    obj=@(p) objective_local(p,Y(:,ib,:),x,gaps,coef,gBasis0,anchorMask,fitMask);
    p=fminsearch(obj,p0,opts); states(ib).deltaGapMm=p(1); states(ib).muGapPerXMm=p(2);
    states(ib).fitRmseMv=sqrt(obj(p));
    holdMask=anchorPool & ~anchorMask;
    states(ib).holdoutRmseMv=sqrt(objective_local(p,Y(:,ib,:),x,gaps,coef,gBasis0,holdMask,fitMask));
    states(ib).valid=isfinite(states(ib).fitRmseMv) && states(ib).fitRmseMv<200 && abs(p(1))<0.8 && abs(p(2))<0.1;
    fprintf('  Blade %d: delta %.5f mm, mu %.6f, fit %.2f mV, holdout %.2f mV, valid=%d\n',b,p(1),p(2),states(ib).fitRmseMv,states(ib).holdoutRmseMv,states(ib).valid);
end
BladeCondition=states; metadata=struct('referenceBladeId',refId,'anchorGapSetMm',anchorGaps(:).','directGapDomainMm',[gLo gHi],'note','Static overlap-only relative blade-state estimates; no vibration and no response-surface rewrite.');
save(fullfile(outDir,'BladeRelativeState_20250527_Blade2Anchor_20260831.mat'),'BladeCondition','metadata','-v7');
T=struct2table(states); writetable(T,fullfile(outDir,'BladeRelativeState_20250527_Blade2Anchor_20260831.csv'));
fprintf('Saved revised blade-state outputs in %s\n',outDir);

function J=objective_local(p,Yb,x,gaps,coef,g0,maskG,maskX)
delta=p(1); mu=p(2); gMat=gaps(maskG(:)).'+delta+mu*(x-x(ceil(numel(x)/2)));
if any(gMat(:)<min(gaps)-1e-8 | gMat(:)>max(gaps)+1e-8), J=1e12; return; end
ix=find(maskX); ig=find(maskG); err=[];
Yall=squeeze(Yb);
for jj=1:numel(ig)
    gq=gMat(ix,jj); F=evalF_local(gq,x(ix),x,coef,g0); y=Yall(ix,ig(jj));
    if any(~isfinite(y))||any(~isfinite(F)), J=1e12; return; end
    err=[err; y-F]; %#ok<AGROW>
end
J=mean(err.^2);
end

function F=evalF_local(g,x,xFull,coef,g0)
% coefficient rows are on the full response x-grid
xx=xFull(:);
c=zeros(numel(x),3); for k=1:3, ok=isfinite(coef(:,k)); c(:,k)=interp1(xx(ok),coef(ok,k),x,'linear','extrap'); end
F=c(:,1)+c(:,2)./g+c(:,3).*log(g./g0);
end
