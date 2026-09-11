%% Strict V3 comparison: held blade contributes only its low-state waveform
% The target waveform and every other gap state of the held blade are hidden.
% This separates an oracle base-gap test from deployable localization+transfer.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
U=load(fullfile(outDir,'StaticWaveformFamily_20251222.mat')); F=U.StaticWaveformFamily;
Y=F.waveformsB2FixedPeakAlignedMv; x=F.xB2OprAxisMm(:); gaps=F.nominalGapMm(:); blades=F.bladeIds(:).';
[nX,nB,nG]=size(Y); gGrid=linspace(min(gaps),max(gaps),281); rows={};

methodNames={ ...
    'R1_pooled_increment', ...
    'R2_oracle_basegap_conditioned', ...
    'R3raw_plus_R2', ...
    'R3raw_plus_R4', ...
    'R3gain_plus_R4', ...
    'R3nearest_surface', ...
    'R3_top2_surface_ensemble', ...
    'R3_top3_surface_ensemble', ...
    'R3_convex_surface_ensemble', ...
    'R3_shifted_nearest_surface', ...
    'R3_peak_aligned_nearest_surface', ...
    'R3_affine_x_nearest_surface', ...
    'R3_weighted_shift_nearest_g1', ...
    'R3_weighted_shift_nearest_g2', ...
    'R3_affine_gain_offset_shifted_surface', ...
    'R3_centered_shift_nearest_surface', ...
    'R3_shape_normalized_shift_nearest_surface', ...
    'R3_regularized_affine_shift_nearest_surface', ...
    'R3_peak_feature_shift_nearest_surface', ...
    'R3_continuous_affine_shifted_surface'};
predCase=nan(nX,nB*(nG-1),numel(methodNames)); lowCase=nan(nX,nB*(nG-1)); targetCase=nan(nX,nB*(nG-1));
caseMeta=nan(nB*(nG-1),3); caseIdx=0;

for ib=1:nB
    trainB=setdiff(1:nB,ib);
    % Training-blade surfaces used only for localization / nearest surface.
    Ctrain=cell(numel(trainB),1);
    for kb=1:numel(trainB), Ctrain{kb}=fitSurface(squeeze(Y(:,trainB(kb),:)),gaps); end
    % Shared basis increment coefficients, learned without the held blade.
    [C1,C2]=fitSharedIncrementBasis(Y(:,trainB,:),gaps);
    % Nonparametric base-gap-conditioned derivative and unconditioned pool.
    [dPooled,dByBase,dBase]=incrementLibraries(Y(:,trainB,:),gaps);

    for il=1:nG-1
        it=il+1; y0=squeeze(Y(:,ib,il)); yt=squeeze(Y(:,ib,it)); dg=gaps(it)-gaps(il);
        caseIdx=caseIdx+1; lowCase(:,caseIdx)=y0; targetCase(:,caseIdx)=yt; caseMeta(caseIdx,:)=[blades(ib),gaps(il),gaps(it)];
        [gRaw,bRaw,~,rawLocRmse]=locateOnTrainingSurfaces(y0,Ctrain,trainB,gGrid,false);
        [gGain,bGain,aGain,gainLocRmse]=locateOnTrainingSurfaces(y0,Ctrain,trainB,gGrid,true);

        d1=dPooled.*dg;
        % dByBase is a secant derivative stored at the transition midpoint.
        d2=interpColumns(dByBase,dBase,gaps(il)+0.5*dg).*dg;
        d3=interpColumns(dByBase,dBase,gRaw+0.5*dg).*dg;
        d4=basisIncrement(C1,C2,gRaw,gRaw+dg);
        d5=basisIncrement(C1,C2,gGain,gGain+dg);
        kRaw=find(trainB==bRaw,1); d6=surface(Ctrain{kRaw},gRaw+dg)-surface(Ctrain{kRaw},gRaw);
        [d7,g7,e7]=manifoldEnsembleIncrement(y0,Ctrain,gGrid,dg,2);
        [d8,g8,e8]=manifoldEnsembleIncrement(y0,Ctrain,gGrid,dg,3);
        [d9,g9,e9]=manifoldConvexIncrement(y0,Ctrain,gGrid,dg);
        [d10,g10,b10,e10]=shiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d11,g11,b11,e11]=peakAlignedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d12,g12,b12,e12]=affineXNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d13,g13,b13,e13]=weightedShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg,1);
        [d14,g14,b14,e14]=weightedShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg,2);
        [d15,g15,b15,a15,e15]=affineShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d16,g16,b16,e16]=centeredShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d17,g17,b17,e17]=shapeNormalizedShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d18,g18,b18,a18,e18]=regularizedAffineShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d19,g19,b19,e19]=peakFeatureShiftedNearestIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        [d20,g20,b20,a20,e20]=continuousAffineShiftedIncrement(y0,x,Ctrain,trainB,gGrid,dg);
        D={d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20};
        locGap=[NaN,gaps(il),gRaw,gRaw,gGain,gRaw,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20]; locBlade=[NaN,NaN,bRaw,bRaw,bGain,bRaw,NaN,NaN,NaN,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20];
        locGain=[NaN,NaN,1,1,aGain,1,1,1,1,1,1,1,1,1,a15,1,1,a18,1,a20]; locRmse=[NaN,NaN,rawLocRmse,rawLocRmse,gainLocRmse,rawLocRmse,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20];
        for im=1:numel(D)
            d=D{im}; ok=isfinite(y0)&isfinite(yt)&isfinite(d); if nnz(ok)<20, continue, end
            predCase(:,caseIdx,im)=y0+d;
            e=(y0(ok)+d(ok))-yt(ok); e0=y0(ok)-yt(ok);
            rmse=sqrt(mean(e.^2)); base=sqrt(mean(e0.^2)); rel=100*rmse/max(std(yt(ok)),eps);
            rows(end+1,:)={blades(ib),gaps(il),gaps(it),methodNames{im},rmse,base,rel,rmse/max(base,eps), ... %#ok<AGROW>
                locGap(im),locGap(im)-gaps(il),locBlade(im),locGain(im),locRmse(im)};
        end
    end
end

T=cell2table(rows,'VariableNames',{'heldBlade','lowGapMm','targetGapMm','method','transitionRMSEmV', ...
    'zeroIncrementRMSEmV','relativeRMSEpct','transferToBaselineRatio','estimatedBaseGapMm', ...
    'baseGapErrorMm','selectedTrainingBlade','localizationGain','localizationRMSEmV'});
writetable(T,fullfile(outDir,'StrictV3_MethodComparison_20251222.csv'));

summaryRows=cell(0,9);
for im=1:numel(methodNames)
    z=T(strcmp(T.method,methodNames{im}),:);
    summaryRows(end+1,:)={methodNames{im},height(z),mean(z.transitionRMSEmV,'omitnan'), ... %#ok<AGROW>
        percentile(z.relativeRMSEpct,90),percentile(z.transferToBaselineRatio,90), ...
        median(abs(z.baseGapErrorMm),'omitnan'),percentile(abs(z.baseGapErrorMm),90), ...
        mean(z.transferToBaselineRatio<1,'omitnan'),mean(z.transferToBaselineRatio<0.5,'omitnan')};
end
S=cell2table(summaryRows,'VariableNames',{'method','N','meanRMSEmV','P90RelativeRMSEpct', ...
    'P90TransferToBaselineRatio','medianAbsBaseGapErrorMm','P90AbsBaseGapErrorMm', ...
    'fractionBetterThanZeroIncrement','fractionBelowHalfBaseline'});
S.P50TransferRMSEmV=nan(height(S),1); S.P90TransferRMSEmV=nan(height(S),1);
for im=1:height(S), z=T(strcmp(T.method,S.method{im}),:); S.P50TransferRMSEmV(im)=percentile(z.transitionRMSEmV,50); S.P90TransferRMSEmV(im)=percentile(z.transitionRMSEmV,90); end
writetable(S,fullfile(outDir,'StrictV3_MethodComparison_Summary_20251222.csv')); disp(S);
save(fullfile(outDir,'StrictV3_WaveformCases_20251222.mat'),'predCase','lowCase','targetCase','caseMeta','methodNames','x','-v7.3');

function [dPool,dByBase,dBase]=incrementLibraries(Y,gaps)
    [nX,nB,nG]=size(Y); deriv=[]; dByBase=nan(nX,nG-1); dBase=nan(nG-1,1);
    for i=1:nG-1
        dg=gaps(i+1)-gaps(i); q=nan(nX,nB);
        for b=1:nB, q(:,b)=(squeeze(Y(:,b,i+1))-squeeze(Y(:,b,i)))./dg; end
        dByBase(:,i)=median(q,2,'omitnan'); dBase(i)=0.5*(gaps(i)+gaps(i+1)); deriv=[deriv,q]; %#ok<AGROW>
    end
    dPool=median(deriv,2,'omitnan');
end

function [C1,C2]=fitSharedIncrementBasis(Y,gaps)
    [nX,nB,nG]=size(Y); C1=nan(nX,1); C2=nan(nX,1); A=[]; D=nan(nB*(nG-1),nX); r=0;
    for b=1:nB
        for i=1:nG-1
            r=r+1; A(r,:)=[1/gaps(i+1)-1/gaps(i),log(gaps(i+1)/gaps(i))];
            D(r,:)=squeeze(Y(:,b,i+1)-Y(:,b,i)).';
        end
    end
    for ix=1:nX, ok=isfinite(D(:,ix)); if nnz(ok)>=4, c=A(ok,:)\D(ok,ix); C1(ix)=c(1); C2(ix)=c(2); end, end
end

function [gBest,bBest,aBest,eBest]=locateOnTrainingSurfaces(y,Ctrain,trainB,gGrid,fitGain)
    gBest=NaN; bBest=NaN; aBest=NaN; eBest=Inf;
    for kb=1:numel(Ctrain)
        for ig=1:numel(gGrid)
            p=surface(Ctrain{kb},gGrid(ig)); ok=isfinite(y)&isfinite(p); if nnz(ok)<20, continue, end
            a=1; if fitGain, a=max(0,(p(ok)'*y(ok))/(p(ok)'*p(ok))); end
            e=sqrt(mean((y(ok)-a*p(ok)).^2));
            if e<eBest, eBest=e; gBest=gGrid(ig); bBest=trainB(kb); aBest=a; end
        end
    end
end

function [d,gWeighted,eWeighted]=manifoldEnsembleIncrement(y,Ctrain,gGrid,dg,kKeep)
    n=numel(Ctrain); bestE=inf(n,1); bestG=nan(n,1); bestD=cell(n,1);
    for kb=1:n
        for ig=1:numel(gGrid)
            p=surface(Ctrain{kb},gGrid(ig)); ok=isfinite(y)&isfinite(p); if nnz(ok)<20, continue, end
            e=sqrt(mean((y(ok)-p(ok)).^2));
            if e<bestE(kb), bestE(kb)=e; bestG(kb)=gGrid(ig); bestD{kb}=surface(Ctrain{kb},gGrid(ig)+dg)-p; end
        end
    end
    [~,ord]=sort(bestE); ord=ord(1:min(kKeep,n)); w=1./max(bestE(ord).^2,eps); w=w/sum(w);
    d=zeros(size(y)); for j=1:numel(ord), d=d+w(j).*bestD{ord(j)}; end
    gWeighted=sum(w.*bestG(ord)); eWeighted=sum(w.*bestE(ord));
end

function [d,gWeighted,eFit]=manifoldConvexIncrement(y,Ctrain,gGrid,dg)
    n=numel(Ctrain); bestE=inf(n,1); bestG=nan(n,1); P=nan(numel(y),n); D=nan(numel(y),n);
    for kb=1:n
        for ig=1:numel(gGrid)
            p=surface(Ctrain{kb},gGrid(ig)); ok=isfinite(y)&isfinite(p); if nnz(ok)<20, continue, end
            e=sqrt(mean((y(ok)-p(ok)).^2));
            if e<bestE(kb), bestE(kb)=e; bestG(kb)=gGrid(ig); P(:,kb)=p; D(:,kb)=surface(Ctrain{kb},gGrid(ig)+dg)-p; end
        end
    end
    ok=isfinite(y)&all(isfinite(P),2)&all(isfinite(D),2);
    if nnz(ok)<20, d=nan(size(y)); gWeighted=NaN; eFit=NaN; return, end
    scale=max(norm(P(ok,:),'fro'),1); lambda=1e3*scale/max(sqrt(nnz(ok)),1);
    w=lsqnonneg([P(ok,:);lambda*ones(1,n)],[y(ok);lambda]);
    if sum(w)<=eps, w=1./max(bestE.^2,eps); end; w=w/sum(w);
    d=D*w; gWeighted=sum(w.*bestG); eFit=sqrt(mean((y(ok)-P(ok,:)*w).^2));
end

function [dBest,gBest,bBest,eBest]=shiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
    tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf;
    for kb=1:numel(Ctrain)
        for ig=1:numel(gGrid)
            p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
            for jt=1:numel(tauGrid)
                p0s=interp1(x,p0,x-tauGrid(jt),'pchip',NaN); ok=isfinite(y)&isfinite(p0s); if nnz(ok)<20, continue, end
                e=sqrt(mean((y(ok)-p0s(ok)).^2));
                if e<eBest
                    p1s=interp1(x,p1,x-tauGrid(jt),'pchip',NaN); dBest=p1s-p0s;
                    gBest=gGrid(ig); bBest=trainB(kb); eBest=e;
                end
            end
        end
    end
end

function [dBest,gBest,bBest,eBest]=peakAlignedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
    dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf; ys=smoothdata(y,'sgolay',31); [~,iy]=max(ys); xpy=x(iy);
    for kb=1:numel(Ctrain)
        for ig=1:numel(gGrid)
            p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
            ps=smoothdata(p0,'sgolay',31); [~,ip]=max(ps); tau=xpy-x(ip);
            p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s); if nnz(ok)<20, continue, end
            e=sqrt(mean((y(ok)-p0s(ok)).^2));
            if e<eBest, p1s=interp1(x,p1,x-tau,'pchip',NaN); dBest=p1s-p0s; gBest=gGrid(ig); bBest=trainB(kb); eBest=e; end
        end
    end
end

function [dBest,gBest,bBest,eBest]=affineXNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
    tauGrid=-0.50:0.10:0.50; kGrid=0.90:0.05:1.10;
    dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf;
    for kb=1:numel(Ctrain)
        for ig=1:numel(gGrid)
            p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
            for ik=1:numel(kGrid)
                for jt=1:numel(tauGrid)
                    xq=kGrid(ik).*(x-tauGrid(jt)); p0s=interp1(x,p0,xq,'pchip',NaN);
                    ok=isfinite(y)&isfinite(p0s); if nnz(ok)<20, continue, end
                    e=sqrt(mean((y(ok)-p0s(ok)).^2));
                    if e<eBest, p1s=interp1(x,p1,xq,'pchip',NaN); dBest=p1s-p0s; gBest=gGrid(ig); bBest=trainB(kb); eBest=e; end
                end
            end
        end
    end
end

function [dBest,gBest,bBest,eBest]=weightedShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg,gamma)
    tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf;
    yy=max(y,0); w=(yy./max(max(yy),eps)).^gamma; w=w./max(sum(w),eps);
    for kb=1:numel(Ctrain)
        for ig=1:numel(gGrid)
            p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
            for jt=1:numel(tauGrid)
                p0s=interp1(x,p0,x-tauGrid(jt),'pchip',NaN); ok=isfinite(y)&isfinite(p0s)&isfinite(w); if nnz(ok)<20, continue, end
                ww=w(ok); ww=ww./max(sum(ww),eps); e=sqrt(sum(ww.*(y(ok)-p0s(ok)).^2));
                if e<eBest, p1s=interp1(x,p1,x-tauGrid(jt),'pchip',NaN); dBest=p1s-p0s; gBest=gGrid(ig); bBest=trainB(kb); eBest=e; end
            end
        end
    end
end

function [dBest,gBest,bBest,aBest,eBest]=affineShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Static localization may use an affine nuisance map for blade-to-blade
% absolute voltage differences. The fitted gain scales the increment; the
% fitted offset is not carried into the predicted transition.
tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; aBest=NaN; eBest=Inf;
for kb=1:numel(Ctrain)
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
        for jt=1:numel(tauGrid)
            tau=tauGrid(jt); p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s);
            if nnz(ok)<20, continue, end
            H=[p0s(ok),ones(nnz(ok),1)]; ab=H\y(ok); a=max(ab(1),0); c=ab(2); %#ok<NASGU>
            e=sqrt(mean((y(ok)-a*p0s(ok)-c).^2));
            if e<eBest
                p1s=interp1(x,p1,x-tau,'pchip',NaN); dBest=a*(p1s-p0s);
                gBest=gGrid(ig); bBest=trainB(kb); aBest=a; eBest=e;
            end
        end
    end
end
end

function [dBest,gBest,bBest,eBest]=centeredShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Remove only the constant level used for static localization. The increment
% itself remains on the measured voltage scale and is not re-centered.
tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf; tauBest=0;
for kb=1:numel(Ctrain)
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
        for jt=1:numel(tauGrid)
            tau=tauGrid(jt); p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s);
            if nnz(ok)<20, continue, end
            yc=y(ok)-median(y(ok)); pc=p0s(ok)-median(p0s(ok)); e=sqrt(mean((yc-pc).^2));
            if e<eBest
                p1s=interp1(x,p1,x-tau,'pchip',NaN); dBest=p1s-p0s;
                gBest=gGrid(ig); bBest=trainB(kb); eBest=e;
            end
        end
    end
end
end

function [dBest,gBest,bBest,eBest]=shapeNormalizedShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Use shape only for blind localization; preserve the donor voltage scale in
% the predicted increment. This separates blade amplitude from gap state.
tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf;
yc=y-median(y(isfinite(y))); ys=max(prctile(y(isfinite(y)),95)-prctile(y(isfinite(y)),5),eps); yn=yc/ys;
for kb=1:numel(Ctrain)
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
        pc=p0-median(p0(isfinite(p0))); ps=max(prctile(p0(isfinite(p0)),95)-prctile(p0(isfinite(p0)),5),eps); pn=pc/ps;
        for jt=1:numel(tauGrid)
            tau=tauGrid(jt); pns=interp1(x,pn,x-tau,'pchip',NaN); ok=isfinite(yn)&isfinite(pns);
            if nnz(ok)<20, continue, end
            e=sqrt(mean((yn(ok)-pns(ok)).^2));
            if e<eBest
                p0s=interp1(x,p0,x-tau,'pchip',NaN); p1s=interp1(x,p1,x-tau,'pchip',NaN);
                dBest=p1s-p0s; gBest=gGrid(ig); bBest=trainB(kb); eBest=e;
            end
        end
    end
end
end

function [dBest,gBest,bBest,aBest,eBest]=regularizedAffineShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Fit a weakly regularized gain/offset for localization so amplitude mismatch
% is allowed without letting it replace the gap coordinate.
tauGrid=-0.50:0.05:0.50; lambda=0.10;
dBest=nan(size(y)); gBest=NaN; bBest=NaN; aBest=NaN; eBest=Inf;
ys=y(isfinite(y)); scale=max(prctile(ys,95)-prctile(ys,5),eps);
for kb=1:numel(Ctrain)
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
        for jt=1:numel(tauGrid)
            tau=tauGrid(jt); p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s);
            if nnz(ok)<20, continue, end
            H=[p0s(ok),ones(nnz(ok),1)]; ab=H\y(ok); a=max(ab(1),0); c=ab(2);
            r=y(ok)-a*p0s(ok)-c; score=mean(r.^2)+lambda*scale^2*((a-1)^2+(c/scale)^2);
            if score<eBest
                p1s=interp1(x,p1,x-tau,'pchip',NaN); dBest=a*(p1s-p0s);
                gBest=gGrid(ig); bBest=trainB(kb); aBest=a; eBest=score;
            end
        end
    end
end
eBest=sqrt(max(eBest,0));
end

function [dBest,gBest,bBest,eBest]=peakFeatureShiftedNearestIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Add low-state peak and tail amplitude as weak localization features.
tauGrid=-0.50:0.05:0.50; dBest=nan(size(y)); gBest=NaN; bBest=NaN; eBest=Inf;
ys=smoothdata(y,'sgolay',31); [~,iy]=max(ys); xpkY=x(iy);
tailY=abs(x-xpkY)>2.5; ampY=max(ys)-median(ys(tailY));
for kb=1:numel(Ctrain)
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig)); p1=surface(Ctrain{kb},gGrid(ig)+dg);
        ps=smoothdata(p0,'sgolay',31); [~,ip]=max(ps); xpkP=x(ip);
        tailP=abs(x-xpkP)>2.5; ampP=max(ps)-median(ps(tailP));
        for jt=1:numel(tauGrid)
            tau=tauGrid(jt); p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s);
            if nnz(ok)<20, continue, end
            r=sqrt(mean((y(ok)-p0s(ok)).^2));
            score=r+25*abs((xpkY-xpkP)-tau)+0.05*abs(ampY-ampP);
            if score<eBest
                p1s=interp1(x,p1,x-tau,'pchip',NaN); dBest=p1s-p0s;
                gBest=gGrid(ig); bBest=trainB(kb); eBest=score;
            end
        end
    end
end
end

function [dBest,gBest,bBest,aBest,eBest]=continuousAffineShiftedIncrement(y,x,Ctrain,trainB,gGrid,dg)
% Refine blind localization continuously in (base gap, spatial shift).
% Gain and offset are nuisance terms for localization; only gain-scaled
% waveform increment is transferred to the held blade.
tauBounds=[-0.50,0.50]; gBounds=[min(gGrid),max(gGrid)];
dBest=nan(size(y)); gBest=NaN; bBest=NaN; aBest=NaN; eBest=Inf;
for kb=1:numel(Ctrain)
    % Coarse seeds keep the local refinement deterministic and bounded.
    bestSeed=[NaN,NaN]; bestSeedE=Inf;
    for ig=1:numel(gGrid)
        p0=surface(Ctrain{kb},gGrid(ig));
        for tau=-0.50:0.10:0.50
            p0s=interp1(x,p0,x-tau,'pchip',NaN); ok=isfinite(y)&isfinite(p0s);
            if nnz(ok)<20, continue, end
            ab=[p0s(ok),ones(nnz(ok),1)]\y(ok); a=max(ab(1),0); c=ab(2);
            ee=sqrt(mean((y(ok)-a*p0s(ok)-c).^2));
            if ee<bestSeedE, bestSeedE=ee; bestSeed=[gGrid(ig),tau]; end
        end
    end
    if ~isfinite(bestSeedE), continue, end
    obj=@(z)localAffineObjective(z,y,x,Ctrain{kb},gBounds,tauBounds);
    opts=optimset('Display','off','MaxIter',250,'TolX',1e-5,'TolFun',1e-4);
    z=fminsearch(obj,bestSeed,opts); z(1)=min(max(z(1),gBounds(1)),gBounds(2)); z(2)=min(max(z(2),tauBounds(1)),tauBounds(2));
    [ee,a,p0s,p1s]=localAffineFit(z,y,x,Ctrain{kb},dg);
    if ee<eBest
        dBest=a*(p1s-p0s); gBest=z(1); bBest=trainB(kb); aBest=a; eBest=ee;
    end
end
end

function v=localAffineObjective(z,y,x,C,gBounds,tauBounds)
if z(1)<gBounds(1)||z(1)>gBounds(2)||z(2)<tauBounds(1)||z(2)>tauBounds(2), v=1e9+1e6*sum(max([gBounds(1)-z(1),z(1)-gBounds(2),tauBounds(1)-z(2),z(2)-tauBounds(2)],0).^2); return, end
[v,~,~,~]=localAffineFit(z,y,x,C,0);
end

function [ee,a,p0s,p1s]=localAffineFit(z,y,x,C,dg)
p0=surface(C,z(1)); p1=surface(C,z(1)+dg);
p0s=interp1(x,p0,x-z(2),'pchip',NaN); p1s=interp1(x,p1,x-z(2),'pchip',NaN);
ok=isfinite(y)&isfinite(p0s); if nnz(ok)<20, ee=1e9; a=0; return, end
ab=[p0s(ok),ones(nnz(ok),1)]\y(ok); a=max(ab(1),0); c=ab(2); ee=sqrt(mean((y(ok)-a*p0s(ok)-c).^2));
end

function C=fitSurface(Y,gaps)
    C=nan(size(Y,1),3); A=[ones(numel(gaps),1),1./gaps,log(gaps)];
    for ix=1:size(Y,1), y=Y(ix,:).'; ok=isfinite(y); if nnz(ok)>=4, C(ix,:)=A(ok,:)\y(ok); end, end
end
function y=surface(C,g), y=C(:,1)+C(:,2)./g+C(:,3).*log(g); end
function d=basisIncrement(C1,C2,g1,g2), d=C1.*(1/g2-1/g1)+C2.*log(g2/g1); end
function y=interpColumns(V,g,x), y=nan(size(V,1),1); for i=1:size(V,1), ok=isfinite(V(i,:)); if nnz(ok)>=2, y(i)=interp1(g(ok),V(i,ok),x,'linear','extrap'); end, end, end
function p=percentile(x,q), x=sort(x(isfinite(x))); if isempty(x),p=NaN;return,end; k=1+(numel(x)-1)*q/100; lo=floor(k); hi=ceil(k); if lo==hi,p=x(lo);else,p=x(lo)+(k-lo)*(x(hi)-x(lo));end,end
