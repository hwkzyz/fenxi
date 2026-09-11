%% Build six donor surfaces on the frozen B2 nominal-gap axis
clear; clc;
thisDir=fileparts(mfilename('fullpath')); rootDir=fileparts(fileparts(thisDir)); outDir=fullfile(thisDir,'results');
U=load(fullfile(outDir,'StaticWaveformFamily_20251222.mat')); F=U.StaticWaveformFamily;
baseFile=fullfile(rootDir,'inputs','calibration','Step05I_OffsetTilt_Shared_Response_Surface_20251222.mat');
B=load(baseFile,'OffsetTiltResponseSurface'); base=B.OffsetTiltResponseSurface;
Y=F.waveformsB2FixedPeakAlignedMv; x=F.xB2OprAxisMm(:); gaps=F.nominalGapMm(:); blades=F.bladeIds(:).';
surfaceDir=fullfile(outDir,'per_blade_response_surfaces'); if exist(surfaceDir,'dir')~=7, mkdir(surfaceDir); end
rows=cell(numel(blades),1);
for ib=1:numel(blades)
    Yi=squeeze(Y(:,ib,:)); C=nan(numel(x),3); A=[ones(numel(gaps),1),1./gaps,log(gaps)];
    for ix=1:numel(x), y=Yi(ix,:).'; ok=isfinite(y); if nnz(ok)>=4, C(ix,:)=A(ok,:)\y(ok); end, end
    Yfit=C(:,1)+C(:,2)./gaps.'+C(:,3).*log(gaps.');
    dFdg=-C(:,2)./(gaps.'.^2)+C(:,3)./gaps.';
    dFdx=nan(size(Yfit)); for j=1:numel(gaps), dFdx(:,j)=gradient(Yfit(:,j),x); end
    OffsetTiltResponseSurface=base; %#ok<NASGU>
    OffsetTiltResponseSurface.method='per_blade_B2_fixed_coordinate_basis_surface';
    OffsetTiltResponseSurface.description=sprintf('Blade %d donor surface; B2 fixed spatial coordinate and frozen nominal gap axis.',blades(ib));
    OffsetTiltResponseSurface.sourceStep05File=F.sourceStep05File;
    OffsetTiltResponseSurface.xGrid=x;
    OffsetTiltResponseSurface.g0Mm=1;
    OffsetTiltResponseSurface.gTrainMm=gaps;
    OffsetTiltResponseSurface.trueGapMm=gaps;
    OffsetTiltResponseSurface.recordedGapMm=gaps;
    OffsetTiltResponseSurface.recordedGapTrainMm=gaps;
    OffsetTiltResponseSurface.coeff=C;
    OffsetTiltResponseSurface.Yref=Yi;
    OffsetTiltResponseSurface.Yfit=Yfit;
    OffsetTiltResponseSurface.dFdgGrid=dFdg;
    OffsetTiltResponseSurface.dFdxGrid=dFdx;
    OffsetTiltResponseSurface.bladeIds=blades(ib);
    OffsetTiltResponseSurface.referenceBladeId=2;
    OffsetTiltResponseSurface.referenceBladeIndex=find(blades==2,1);
    OffsetTiltResponseSurface.note=sprintf('Donor blade %d; nominal gaps are the frozen B2 settings.',blades(ib));
    file=fullfile(surfaceDir,sprintf('PerBlade_ResponseSurface_20251222_B%d.mat',blades(ib)));
    save(file,'OffsetTiltResponseSurface','-v7.3');
    e=Yfit-Yi; rows{ib}=table(blades(ib),sqrt(mean(e(isfinite(e)).^2)),string(file), ...
        'VariableNames',{'blade','surfaceFitRMSEmV','file'});
end
T=vertcat(rows{:}); writetable(T,fullfile(outDir,'PerBlade_ResponseSurface_Summary_20251222.csv')); disp(T);
