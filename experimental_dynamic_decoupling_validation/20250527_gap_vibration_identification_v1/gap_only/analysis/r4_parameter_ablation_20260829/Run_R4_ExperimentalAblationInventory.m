function T=Run_R4_ExperimentalAblationInventory()
% Inventory-only comparison of existing formal fixed/gap-aware fits.
% No formal result is modified and no identification is rerun.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
cfg=Config_20250527(); f=cfg.files.formalGapResult;
if exist(f,'file')~=2, error('Formal gap-aware result not found: %s',f); end
s=load(f); R=s.Result; W=R.WindowResult; rows={};
for i=1:numel(W)
 w=W(i); if ~isfield(w,'modelFits')||~isfield(w.modelFits,'fixed')||~isfield(w.modelFits,'gap_only'),continue;end
 a=w.modelFits.fixed; b=w.modelFits.gap_only;
 rows(end+1,:)={w.windowId,a.EO,b.EO,a.freqHz,b.freqHz,a.amplitudeMm,b.amplitudeMm,...
     a.deltaGapMm,b.deltaGapMm,a.plainRmseMv,b.plainRmseMv,...
     b.amplitudeMm-a.amplitudeMm,b.plainRmseMv-a.plainRmseMv}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'window','EO_fixed','EO_gapOnly','freq_fixed_Hz','freq_gapOnly_Hz',...
 'A_fixed_mm','A_gapOnly_mm','dg_fixed_mm','dg_gapOnly_mm','rmse_fixed_mV','rmse_gapOnly_mV',...
 'deltaA_gapMinusFixed_mm','deltaRmse_gapMinusFixed_mV'});
outDir=fileparts(mfilename('fullpath')); if exist(outDir,'dir')~=7,mkdir(outDir);end
writetable(T,fullfile(outDir,'formal_fixed_vs_gap_only.csv')); save(fullfile(outDir,'formal_fixed_vs_gap_only.mat'),'T','-v7.3');
fprintf('Inventory complete: %d comparable windows.\n',height(T));
end
