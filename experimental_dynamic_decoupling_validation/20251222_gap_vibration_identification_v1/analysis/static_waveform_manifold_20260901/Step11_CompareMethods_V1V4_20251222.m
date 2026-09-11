%% Compare method families under the V1--V4 evidence ledger
clear; clc; thisDir=fileparts(mfilename('fullpath')); rootDir=fileparts(fileparts(thisDir)); outDir=fullfile(thisDir,'results'); gapDir=fullfile(rootDir,'results','gap_aware'); cmpDir=fullfile(rootDir,'results','comparison');
rows={};
% Current method: V1--V3 are direct static leave-out tests. The available
% V4 pair is a free-b_s versus b_s=0 closure regression, not an independent
% high-speed truth set for the proposed route; keep that distinction explicit.
G=readtable(fullfile(outDir,'ValidationGate_20251222.csv')); passAll=all(G.pass); ev=strjoin(G.evidence,'; ');
rows(end+1,:)={'Method-A_B2-fixed_increment', 'B2 fixed coordinate; adjacent same-blade and cross-blade increment transfer', status(G.pass(1)),status(G.pass(2)),status(G.pass(3)),'provisional', [ev '; V4 is closure-equivalence only']}; %#ok<AGROW>
% Legacy fixed-gap/gap-aware comparison: useful V4 diagnostic, not V1--V3.
L=readtable(fullfile(cmpDir,'Summary_R1_R7_R4_20251222_V1.csv'));
for i=1:height(L)
 nWin = 18; % all three legacy regions contain 18 sliding windows
 rows(end+1,:)={['Legacy_fixed_gap_' char(string(L.Region(i)))], 'Original no-gap/fixed-gap baseline', 'not_run','not_run','not_run','reported',sprintf('EO %d/%d; RMSE %.3f mV',L.FixedCorrectCount(i),nWin,L.FixedMeanPlainRmseMv(i))}; %#ok<AGROW>
 rows(end+1,:)={['Legacy_gap_aware_' char(string(L.Region(i)))], 'Original static gap-aware increment method', 'not_run','not_run','not_run','reported',sprintf('EO %d/%d; RMSE %.3f mV; mean dG %.5f mm',L.GapCorrectCount(i),nWin,L.GapMeanPlainRmseMv(i),L.GapMeanDgMm(i))}; %#ok<AGROW>
end
% Existing experimental ablations: V4 metrics only, because their static
% validation was not run under the fixed B2-coordinate protocol.
 [rows,ok] = addSummary(rows,'Legacy_free_bs_B1','free intercept, B1',fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_unified_structured_voltage_r01_20260811.csv'));
 [rows,ok] = addSummary(rows,'No_bs_B1','b_s=0, B1',fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_nob_formal_B1_20260831.csv'));
 [rows,ok] = addSummary(rows,'Gain1_B1','a_s=1 ablation, B1',fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_gain1_ablation.csv'));
 [rows,ok] = addSummary(rows,'Legacy_free_bs_B5','free intercept, B5',fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B5_S123_unified_structured_voltage_r04_20260811.csv'));
 [rows,ok] = addSummary(rows,'No_bs_B5','b_s=0, B5',fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B5_S123_nob_formal_B5_20260830.csv'));
T=cell2table(rows,'VariableNames',{'method','description','V1','V2','V3','V4','evidence'}); writetable(T,fullfile(outDir,'MethodComparison_V1V4_20251222.csv')); disp(T);
% Keep static evidence (V1--V3) separate from dynamic closure evidence (V4).
writetable(T(:,{'method','description','V1','V2','V3','evidence'}),fullfile(outDir,'MethodComparison_Static_V1V3_20251222.csv'));
writetable(T(:,{'method','description','V4','evidence'}),fullfile(outDir,'MethodComparison_Dynamic_V4Closure_20251222.csv'));
% Also write a compact numeric table for the legacy and ablation results.
summaryIds = {'Legacy_free_bs_B1','No_bs_B1','Gain1_B1','Legacy_free_bs_B5','No_bs_B5'};
summaryPaths = {fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_unified_structured_voltage_r01_20260811.csv'), fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_nob_formal_B1_20260831.csv'), fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B1_S123_gain1_ablation.csv'), fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B5_S123_unified_structured_voltage_r04_20260811.csv'), fullfile(gapDir,'Main_GapAware_VPTopK_FullWave_Summary_20251222_B5_S123_nob_formal_B5_20260830.csv')};
hs = cell(0,6);
for k=1:numel(summaryIds)
    if isfile(summaryPaths{k})
        q=readtable(summaryPaths{k}); q=q(1,:); dg=NaN;
        if ismember('mean_delta_gap_mm',q.Properties.VariableNames), dg=q.mean_delta_gap_mm; end
        hs(end+1,:)={summaryIds{k},q.dominant_EO,q.mean_frequency_hz,q.mean_amplitude_mm,q.mean_rmse_mV,dg}; %#ok<AGROW>
    end
end
for i=1:height(L)
    hs(end+1,:)={['Legacy_fixed_gap_' char(string(L.Region(i)))],L.ExpectedEO(i),L.FixedMeanFreqCorrectHz(i),L.FixedMeanAmpAllMm(i),L.FixedMeanPlainRmseMv(i),NaN}; %#ok<AGROW>
    hs(end+1,:)={['Legacy_gap_aware_' char(string(L.Region(i)))],L.ExpectedEO(i),L.GapMeanFreqCorrectHz(i),L.GapMeanAmpAllMm(i),L.GapMeanPlainRmseMv(i),L.GapMeanDgMm(i)}; %#ok<AGROW>
end
HS=cell2table(hs,'VariableNames',{'method','expectedOrDominantEO','meanFrequencyHz','meanAmplitudeMm','meanRMSEmV','meanDeltaGapMm'}); writetable(HS,fullfile(outDir,'MethodComparison_HighSpeed_20251222.csv')); disp(HS);
function [rows,ok] = addSummary(rows,id,desc,p)
    ok=false;
    if ~isfile(p), return, end
    Q=readtable(p); q=Q(1,:);
    rows(end+1,1:7) = {id,desc,'not_run','not_run','not_run','reported', ...
        sprintf('EO %g; f %.6f Hz; A %.6f mm; RMSE %.6f mV',q.dominant_EO,q.mean_frequency_hz,q.mean_amplitude_mm,q.mean_rmse_mV)};
    ok=true;
end
function s=status(tf)
    if tf, s='pass'; else, s='fail'; end
end
