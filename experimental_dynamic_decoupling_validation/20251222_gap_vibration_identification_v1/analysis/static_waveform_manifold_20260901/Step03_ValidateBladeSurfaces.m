%% 03: Interleaved leave-gap-out validation for every blade surface
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir,'results');
familyFile = fullfile(outDir,'StaticWaveformFamily_20251222.mat');
S=load(familyFile,'StaticWaveformFamily'); F=S.StaticWaveformFamily;
xi=F.xiAxisMm(:); Y=F.waveformsPeakAlignedMv; blades=F.bladeIds(:).'; gaps=F.nominalGapMm(:).';
nX=numel(xi); nB=numel(blades); nG=numel(gaps); g0=1.0; rows=cell(nB*ceil(nG/2),1); k=0;
for ib=1:nB
    testIdx=2:2:nG; trainIdx=setdiff(1:nG,testIdx);
    Atrain=[ones(numel(trainIdx),1),1./gaps(trainIdx).',log(gaps(trainIdx).'/g0)];
    C=nan(nX,3);
    for ix=1:nX
        y=squeeze(Y(ix,ib,trainIdx)); valid=isfinite(y);
        if nnz(valid)>=3, C(ix,:)= (Atrain(valid,:) \ y(valid)).'; end
    end
    for ig=testIdx
        pred=C(:,1)+C(:,2)./gaps(ig)+C(:,3).*log(gaps(ig)/g0);
        obs=Y(:,ib,ig); valid=isfinite(obs)&isfinite(pred);
        if nnz(valid)<20, continue; end
        k=k+1; rmse=sqrt(mean((obs(valid)-pred(valid)).^2)); rel=rmse/max(abs(obs(valid)));
        rows{k}=table(blades(ib),gaps(ig),nnz(valid),rmse,rel,'VariableNames',{'bladeId','heldOutGapMm','pointCount','rmseMv','relativeRmse'});
    end
end
validationTable=vertcat(rows{1:k}); summary=groupsummary(validationTable,'bladeId',{'mean','std'},{'rmseMv','relativeRmse'});
writetable(validationTable,fullfile(outDir,'03_LeaveGapOutValidation_20251222.csv')); writetable(summary,fullfile(outDir,'03_LeaveGapOutValidationSummary_20251222.csv'));
fig=figure('Color','w','Position',[120 120 950 480]); boxchart(categorical(validationTable.bladeId),validationTable.rmseMv); grid on; xlabel('Blade ID'); ylabel('Held-out waveform RMSE (mV)'); title('Within-blade interleaved leave-gap-out validation');
saveas(fig,fullfile(outDir,'03_LeaveGapOutValidation_20251222.png')); disp(summary);
