function [A,B] = Method_Main10_B2Donor_vs_Formal_Load_20251222(dataDir)
% A: B2-anchored per-blade surface candidate; B: formal shared no-b_s result.
aFile=fullfile(dataDir,'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_manifold_donorB2_nob_R4_full_20260901.csv');
bFile=fullfile(dataDir,'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_nob_formal_B5_20260830.csv');
if ~isfile(aFile)||~isfile(bFile), A=table; B=table; return; end
A=readtable(aFile); B=readtable(bFile);
end
