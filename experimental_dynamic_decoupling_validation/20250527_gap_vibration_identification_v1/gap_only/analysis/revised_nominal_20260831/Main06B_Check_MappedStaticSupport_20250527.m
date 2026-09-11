%% Main06B: classify all static states after freezing blade-relative states
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'outputs');
R=load(fullfile(outDir,'Step05_Response_Surface_20250527.mat'),'responseSurface'); R=R.responseSurface;
B=load(fullfile(outDir,'BladeRelativeState_20250527_Blade2Anchor_20260831.mat'),'BladeCondition'); B=B.BladeCondition;
x=R.xGrid(:); gaps=R.trueGapMm(:); gLo=min(gaps); gHi=max(gaps); x0=0; nB=numel(R.bladeIds); nG=numel(gaps);
rows={};
for ib=1:nB
  b=R.bladeIds(ib); st=B([B.bladeId]==b); d=st.deltaGapMm; m=st.muGapPerXMm;
  for ig=1:nG
    gPath=gaps(ig)+d+m*(x-x0); inside=(gPath>=gLo-1e-9 & gPath<=gHi+1e-9);
    if all(inside), level="direct_or_mapped_inside"; elseif any(inside), level="partial_mapped_extension"; else, level="unsupported"; end
    rows{end+1,1}=table(b,gaps(ig),min(gPath),max(gPath),nnz(inside)/numel(inside),level,...
      'VariableNames',{'bladeId','gapSetpointMm','mappedGapMinMm','mappedGapMaxMm','directSupportFraction','supportLevel'}); %#ok<AGROW>
  end
end
T=vertcat(rows{:}); writetable(T,fullfile(outDir,'BladeMappedStaticSupport_20250527_Blade2Anchor_20260831.csv'));
S=groupsummary(T,'bladeId',{'mean','min','max'},{'directSupportFraction'}); writetable(S,fullfile(outDir,'BladeMappedStaticSupport_Summary_20250527_Blade2Anchor_20260831.csv'));
disp(T); disp(S);
