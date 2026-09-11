%% Compare proposed increment routes on identical held-out transitions
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
U=load(fullfile(outDir,'StaticWaveformFamily_20251222.mat')); F=U.StaticWaveformFamily;
Y=F.waveformsB2FixedPeakAlignedMv; x=F.xB2OprAxisMm(:); gaps=F.nominalGapMm(:); blades=F.bladeIds(:).';
[nX,nB,nG]=size(Y); rows=[];
for ib=1:nB
    trainB=setdiff(1:nB,ib);
    for il=1:nG-1
        it=il+1; y0=squeeze(Y(:,ib,il)); yt=squeeze(Y(:,ib,it));
        dAll=[]; dSame=[];
        for jb=trainB
            yb0=squeeze(Y(:,jb,il)); ybt=squeeze(Y(:,jb,it));
            d= ybt-yb0; dSame=[dSame,d]; %#ok<AGROW>
            for lk=1:nG-1
                dAll=[dAll,squeeze(Y(:,jb,lk+1))-squeeze(Y(:,jb,lk))]; %#ok<AGROW>
            end
        end
        % Route 1: pooled increment, ignoring base gap.
        d1=nanmedian(dAll,2);
        % Route 2: base-gap-conditioned increment. With equal nominal gaps,
        % this is the cross-blade median for the same transition.
        d2=nanmedian(dSame,2);
        % Route 4: target-blade surface difference, with target gap held out.
        C=fitLeaveGap(squeeze(Y(:,ib,:)),gaps,it,g0Value());
        d4=surface(C,gaps(it),g0Value())-surface(C,gaps(il),g0Value());
        for im=1:3
            if im==1, d=d1; elseif im==2, d=d2; else, d=d4; end
            ok=isfinite(y0)&isfinite(yt)&isfinite(d); if nnz(ok)<20, continue, end
            yh=y0(ok)+d(ok); e=yh-yt(ok); e0=y0(ok)-yt(ok);
            rows(end+1,:)=[blades(ib),gaps(il),gaps(it),im,sqrt(mean(e.^2)),sqrt(mean(e0.^2)), ...
                100*sqrt(mean(e.^2))/max(std(yt(ok)),eps),sqrt(mean(e.^2))/max(sqrt(mean(e0.^2)),eps)]; %#ok<AGROW>
        end
    end
end
T=array2table(rows,'VariableNames',{'heldBlade','lowGapMm','targetGapMm','route','transitionRMSEmV','zeroIncrementRMSEmV','relativeRMSEpct','transferToBaselineRatio'});
writetable(T,fullfile(outDir,'IncrementRouteComparison_V1V3_20251222.csv'));
S=groupsummary(T,'route',{'mean','median','std'}, {'transitionRMSEmV','relativeRMSEpct','transferToBaselineRatio'});
writetable(S,fullfile(outDir,'IncrementRouteComparison_V1V3_Summary_20251222.csv')); disp(S);
for r=1:3
    z=T(T.route==r,:); fprintf('Route %d: N=%d, P90 rel RMSE=%.3f%%, P90 transfer/baseline=%.3f, mean RMSE=%.3f mV\n',r,height(z),percentile(z.relativeRMSEpct,90),percentile(z.transferToBaselineRatio,90),mean(z.transitionRMSEmV,'omitnan'));
end
function C=fitLeaveGap(Y,gaps,held,g0)
    C=nan(size(Y,1),3); tr=setdiff(1:numel(gaps),held);
    for ix=1:size(Y,1)
        y=Y(ix,tr); ok=isfinite(y); if nnz(ok)>=4
            A=[ones(nnz(ok),1),1./gaps(tr(ok)),log(gaps(tr(ok))/g0)]; C(ix,:)=A\y(ok).';
        end
    end
end
function y=surface(C,g,g0), y=C(:,1)+C(:,2)./g+C(:,3).*log(g/g0); end
function g0=g0Value(), g0=1; end
function p=percentile(x,q), x=sort(x(isfinite(x))); if isempty(x),p=NaN;return,end; k=1+(numel(x)-1)*q/100; lo=floor(k); hi=ceil(k); if lo==hi,p=x(lo);else,p=x(lo)+(k-lo)*(x(hi)-x(lo));end,end
