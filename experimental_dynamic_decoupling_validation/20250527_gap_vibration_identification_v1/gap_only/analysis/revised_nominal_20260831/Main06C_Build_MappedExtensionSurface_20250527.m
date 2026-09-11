%% Main06C: build a reference-anchored mapped extension (diagnostic only)
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'outputs');
A=load(fullfile(outDir,'Step05_Response_Surface_20250527.mat'),'responseSurface'); R=A.responseSurface;
B=load(fullfile(outDir,'BladeRelativeState_20250527_Blade2Anchor_20260831.mat'),'BladeCondition'); B=B.BladeCondition;
x=R.xGrid(:); gaps=R.trueGapMm(:); nX=numel(x); nB=numel(R.bladeIds); nG=numel(gaps); x0=0;
gMap=nan(nX,nB,nG); Y=R.waveforms; direct=nan(nX,nB,nG);
for ib=1:nB
 st=B([B.bladeId]==R.bladeIds(ib));
 for ig=1:nG
  gMap(:,ib,ig)=gaps(ig)+st.deltaGapMm+st.muGapPerXMm*(x-x0);
  direct(:,ib,ig)=abs(gMap(:,ib,ig)-gaps(ig))<1e-12 & ib==find(R.bladeIds==2,1);
 end
end
coef=nan(nX,3); gMinX=nan(nX,1); gMaxX=nan(nX,1); nObs=zeros(nX,1); supportLevel=strings(nX,1);
for ix=1:nX
 gs=squeeze(gMap(ix,:,:)); ys=squeeze(Y(ix,:,:)); gs=gs(:); ys=ys(:);
 ok=isfinite(gs)&isfinite(ys)&gs>0; gs=gs(ok); ys=ys(ok); if isempty(gs),continue;end
 w=ones(size(gs)); % Blade-2 direct anchor has priority inside direct interval
 labels=repmat(R.bladeIds(:),nG,1); %#ok<NASGU>
 % reshape ordering is blade-major after squeeze(:)
 bladeRep=repelem(R.bladeIds(:),nG); w(bladeRep==2 & gs>=min(gaps)-1e-9 & gs<=max(gaps)+1e-9)=5;
 H=[ones(size(gs)),1./gs,log(gs./R.g0Mm)];
 % weighted least squares (written explicitly for numerical clarity)
 W=spdiags(w,0,numel(w),numel(w)); coef(ix,:)=(H'*W*H)\(H'*W*ys);
 gMinX(ix)=min(gs); gMaxX(ix)=max(gs); nObs(ix)=numel(gs);
 directGs=gs(gs>=min(gaps)-1e-9 & gs<=max(gaps)+1e-9);
 if isempty(directGs), supportLevel(ix)="mapped_extension_only"; elseif min(gs)<min(gaps)-1e-9 || max(gs)>max(gaps)+1e-9, supportLevel(ix)="direct_plus_mapped_extension"; else, supportLevel(ix)="direct"; end
end
ExtendedResponseSurface=struct('dataset',R.dataset,'referenceBladeId',2,'xGrid',x,'gReferenceMm',gaps,'gDirectDomainMm',[min(gaps) max(gaps)],'gMappedMinByX',gMinX,'gMappedMaxByX',gMaxX,'coeff',coef,'nObservationsByX',nObs,'supportLevelByX',supportLevel,'sourceDirectLibrary',fullfile(outDir,'Step05_Response_Surface_20250527.mat'),'sourceBladeState',fullfile(outDir,'BladeRelativeState_20250527_Blade2Anchor_20260831.mat'),'note','Diagnostic mapped extension only; Blade-2 direct nominal gap axis remains the physical anchor.');
save(fullfile(outDir,'MappedExtensionSurface_20250527_Blade2Anchor_20260831.mat'),'ExtendedResponseSurface','-v7');
T=table(x,gMinX,gMaxX,nObs,supportLevel,'VariableNames',{'xMm','mappedGapMinMm','mappedGapMaxMm','nObservations','supportLevel'}); writetable(T,fullfile(outDir,'MappedExtensionSurface_ByX_20250527_Blade2Anchor_20260831.csv'));
fprintf('Saved mapped extension diagnostic surface. g support %.4f..%.4f mm\n',min(gMinX,[],'omitnan'),max(gMaxX,[],'omitnan'));
