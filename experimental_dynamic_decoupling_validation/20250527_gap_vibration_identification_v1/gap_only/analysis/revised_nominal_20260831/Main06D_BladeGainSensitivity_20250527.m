%% Main06D: sensitivity of blade-relative state to a free blade gain
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'outputs');
A=load(fullfile(outDir,'Step05_Response_Surface_20250527.mat'),'responseSurface'); R=A.responseSurface;
x=R.xGrid(:); gaps=R.trueGapMm(:); coef=R.coeff; g0=R.g0Mm; W=R.effectiveWindow(:); idx=find(W); refId=2;
pool=gaps<=1.70; ii=find(pool); fitMask=false(size(gaps)); fitMask(ii(1:2:end))=true; holdMask=pool&~fitMask;
rows={};
for ib=1:numel(R.bladeIds)
 b=R.bladeIds(ib); Yall=squeeze(R.waveforms(:,ib,:));
 if b==refId, continue; end
 obj=@(p) obj3(p,Yall,x,idx,gaps,coef,g0,fitMask);
 p=fminsearch(obj,[0 0 1],optimset('Display','off','MaxIter',3000,'MaxFunEvals',8000));
 rmse=sqrt(obj(p)); holdrmse=sqrt(obj3(p,Yall,x,idx,gaps,coef,g0,holdMask));
 rows{end+1,1}=table(b,p(1),p(2),p(3),rmse,holdrmse,'VariableNames',{'bladeId','deltaGainFreeMm','muGainFreeMmPerX','relativeBladeGain','fitRmseMv','holdoutRmseMv'}); %#ok<AGROW>
 fprintf('Blade %d: delta %.5f, mu %.6f, c %.4f, fit %.2f, hold %.2f mV\n',b,p(1),p(2),p(3),rmse,holdrmse);
end
T=vertcat(rows{:}); writetable(T,fullfile(outDir,'BladeGainSensitivity_20250527_Blade2Anchor_20260831.csv'));

function J=obj3(p,Yall,x,idx,gaps,coef,g0,maskG)
 d=p(1); m=p(2); cGain=p(3); if cGain<=0, J=1e12; return; end
 gMat=gaps(maskG).'+d+m*(x-x(ceil(numel(x)/2))); if any(gMat(:)<min(gaps)-1e-8|gMat(:)>max(gaps)+1e-8),J=1e12;return;end
 ii=find(maskG); e=[];
 for j=1:numel(ii)
  gq=gMat(idx,j); F=evalF(gq,x(idx),x,coef,g0); y=Yall(idx,ii(j)); if any(~isfinite(F)|~isfinite(y)),J=1e12;return;end; e=[e;y-cGain*F]; %#ok<AGROW>
 end
 J=mean(e.^2);
end
function F=evalF(g,xq,xfull,coef,g0)
 c=zeros(numel(xq),3); for k=1:3,ok=isfinite(coef(:,k));c(:,k)=interp1(xfull(ok),coef(ok,k),xq,'linear','extrap');end
 F=c(:,1)+c(:,2)./g+c(:,3).*log(g./g0);
end
