%% V4 wrapper: high-speed closure validation protocol only
clear; clc; thisDir=fileparts(mfilename('fullpath')); addpath(fullfile(thisDir,'methods')); addpath(fullfile(thisDir,'validation'));
rootDir=fileparts(fileparts(thisDir)); dataDir=fullfile(rootDir,'results','gap_aware');
[A,B]=Method_Main10Pair_Load_20251222(dataDir);
[S,T]=V4_HighSpeedClosedLoop_Framework_20251222(A,B,'free_b_vs_no_b_closure');
outDir=fullfile(thisDir,'results'); writetable(S,fullfile(outDir,'V4_ExperimentalClosure_Summary_20251222.csv'));
writetable(T,fullfile(outDir,'V4_ExperimentalClosure_PerWindow_20251222.csv')); fprintf('V4 framework complete: %d windows.\n',height(T));
