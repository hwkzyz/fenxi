%% V3 wrapper: fixed validation framework + replaceable method interface
clear; clc; close all;
thisDir=fileparts(mfilename('fullpath')); addpath(fullfile(thisDir,'methods'));
T=V3_LeaveOneBladeTransfer_Framework_20251222(@V3_DirectIncrementBasis_Train_20251222,@V3_DirectIncrementBasis_Predict_20251222,'direct_increment_basis_v1');
outDir=fullfile(thisDir,'results'); writetable(T,fullfile(outDir,'V3_LeaveOneBlade_IncrementTransfer_20251222.csv'));
writetable(groupsummary(T,'heldOutBlade',{'mean','std'},T.Properties.VariableNames(4:7)),fullfile(outDir,'V3_LeaveOneBlade_IncrementTransfer_Summary_20251222.csv'));
fprintf('V3 framework complete: %d held-out transitions.\n',height(T));
