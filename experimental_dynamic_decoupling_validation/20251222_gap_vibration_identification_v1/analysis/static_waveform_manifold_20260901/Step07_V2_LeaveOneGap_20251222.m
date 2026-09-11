%% V2 wrapper: validation protocol only
clear; clc; close all; thisDir=fileparts(mfilename('fullpath')); addpath(fullfile(thisDir,'methods')); addpath(fullfile(thisDir,'validation'));
T=V2_LeaveOneGap_Framework_20251222(@Method_StaticSurface_Train_20251222,@Method_StaticSurface_Predict_20251222,'static_surface_basis_v1'); outDir=fullfile(thisDir,'results'); writetable(T,fullfile(outDir,'V2_LeaveOneGap_20251222.csv')); writetable(groupsummary(T,'blade',{'mean','std'},T.Properties.VariableNames(3:6)),fullfile(outDir,'V2_LeaveOneGap_Summary_20251222.csv'));
