%% V4 engineering closure for the selected B2-anchored donor candidate
clear; clc;
thisDir=fileparts(mfilename('fullpath')); addpath(fullfile(thisDir,'methods')); addpath(fullfile(thisDir,'validation'));
rootDir=fileparts(fileparts(thisDir)); dataDir=fullfile(rootDir,'results','gap_aware'); outDir=fullfile(thisDir,'results');
[A,B]=Method_Main10_B2Donor_vs_Formal_Load_20251222(dataDir);
[S,T]=V4_HighSpeedClosedLoop_Framework_20251222(A,B,'B2_donor_no_b_vs_formal_shared_no_b');
writetable(S,fullfile(outDir,'V4_B2DonorCandidate_Summary_20251222.csv'));
writetable(T,fullfile(outDir,'V4_B2DonorCandidate_PerWindow_20251222.csv'));
fprintf('B2 donor V4 framework complete: %d windows, %d EO matches.\n',S.nWindows,S.nEOMatches);
