%% Select the best deployable candidate under an explicit operational V3 screen.
% V3 has no absolute mV gate until a repeatability reference is available.
% This screen reuses the fixed V1/V2 relative-RMSE criterion and requires
% that the transfer improve over the zero-increment baseline in >=95% cases.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); outDir=fullfile(thisDir,'results');
S=readtable(fullfile(outDir,'StrictV3_MethodComparison_Summary_20251222.csv'),'TextType','string');
S=S(~contains(S.method,"oracle") & ~startsWith(S.method,"R1_"),:);
S.V3OperationalPass=S.P90RelativeRMSEpct<=10 & S.fractionBetterThanZeroIncrement>=0.95;
S=sortrows(S,{'V3OperationalPass','P90TransferRMSEmV'},{'descend','ascend'});
writetable(S,fullfile(outDir,'V1V2V3_OperationalCandidateAssessment_20251222.csv'));
disp(S(:,{'method','V3OperationalPass','P90RelativeRMSEpct','P90TransferRMSEmV', ...
    'P90TransferToBaselineRatio','fractionBetterThanZeroIncrement'}));
fprintf('Operational V3 candidates: %d/%d\n',nnz(S.V3OperationalPass),height(S));
