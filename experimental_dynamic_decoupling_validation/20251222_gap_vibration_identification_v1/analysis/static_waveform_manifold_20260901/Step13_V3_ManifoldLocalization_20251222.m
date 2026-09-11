%% Route 3: hidden static-state manifold localization
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
U=load(fullfile(outDir,'StaticWaveformFamily_20251222.mat')); F=U.StaticWaveformFamily;
Y=F.waveformsB2FixedPeakAlignedMv; x=F.xB2OprAxisMm(:); gaps=F.nominalGapMm(:); blades=F.bladeIds(:).';
[~,nB,nG]=size(Y); rows=[]; gGrid=linspace(min(gaps),max(gaps),281);
for ib=1:nB
    for ih=1:nG
        C=fitLeaveGap(squeeze(Y(:,ib,:)),gaps,ih,1);
        yt=squeeze(Y(:,ib,ih)); best=[Inf NaN NaN];
        for kg=1:numel(gGrid)
            yp=surface(C,gGrid(kg),1); ok=isfinite(yt)&isfinite(yp);
            if nnz(ok)<20, continue, end
            e=sqrt(mean((yp(ok)-yt(ok)).^2)); if e<best(1), best=[e,gGrid(kg),100*e/max(std(yt(ok)),eps)]; end
        end
        rows(end+1,:)=[blades(ib),gaps(ih),best(2),best(2)-gaps(ih),best(1),best(3)]; %#ok<AGROW>
    end
end
T=array2table(rows,'VariableNames',{'blade','trueGapMm','estimatedGapMm','gapErrorMm','localizationRMSEmV','relativeRMSEpct'});
writetable(T,fullfile(outDir,'V3_ManifoldLocalization_20251222.csv'));
S=groupsummary(T,'blade',{'mean','median','std'},{'gapErrorMm','localizationRMSEmV','relativeRMSEpct'});
writetable(S,fullfile(outDir,'V3_ManifoldLocalization_Summary_20251222.csv')); disp(S);
fprintf('Route 3 localization: N=%d, median |gap error|=%.5f mm, P90 |gap error|=%.5f mm, P90 rel RMSE=%.3f%%\n',height(T),median(abs(T.gapErrorMm)),percentile(abs(T.gapErrorMm),90),percentile(T.relativeRMSEpct,90));
function C=fitLeaveGap(Y,gaps,held,g0), C=nan(size(Y,1),3); tr=setdiff(1:numel(gaps),held); for ix=1:size(Y,1), y=Y(ix,tr); ok=isfinite(y); if nnz(ok)>=4, A=[ones(nnz(ok),1),1./gaps(tr(ok)),log(gaps(tr(ok))/g0)]; C(ix,:)=A\y(ok).'; end,end,end
function y=surface(C,g,g0), y=C(:,1)+C(:,2)./g+C(:,3).*log(g/g0); end
function p=percentile(x,q), x=sort(x(isfinite(x))); if isempty(x),p=NaN;return,end; k=1+(numel(x)-1)*q/100; lo=floor(k); hi=ceil(k); if lo==hi,p=x(lo);else,p=x(lo)+(k-lo)*(x(hi)-x(lo));end,end
