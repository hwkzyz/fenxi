function T = V3_LeaveOneBladeTransfer_Framework_20251222(trainFcn,predictFcn,methodId)
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results'); U=load(fullfile(outDir,'StaticWaveformFamily_20251222.mat')); F=U.StaticWaveformFamily;
Y=F.waveformsB2FixedPeakAlignedMv; gaps=F.nominalGapMm(:); blades=F.bladeIds(:).'; [~,nB,nG]=size(Y); rows=[];
for ih=1:nB, tr=setdiff(1:nB,ih); model=trainFcn(Y(:,tr,:),gaps); for il=1:nG-1, it=il+1; y0=squeeze(Y(:,ih,il)); yt=squeeze(Y(:,ih,it)); d=predictFcn(model,gaps(il),gaps(it),y0); yh=y0+d; ok=isfinite(yh)&isfinite(yt)&isfinite(y0); if nnz(ok)<20,continue,end; e=yh(ok)-yt(ok); e0=y0(ok)-yt(ok); rows(end+1,:)=[blades(ih),gaps(il),gaps(it),sqrt(mean(e.^2)),sqrt(mean(e0.^2)),sqrt(mean(e.^2))/max(sqrt(mean(e0.^2)),eps),mean(d(ok))]; end, end
T=array2table(rows,'VariableNames',{'heldOutBlade','lowGapMm','targetGapMm','incrementTransferRMSEmV','zeroIncrementRMSEmV','relativeToBaseline','meanPredictedIncrementmV'}); T.methodId=repmat(string(methodId),height(T),1);
end
