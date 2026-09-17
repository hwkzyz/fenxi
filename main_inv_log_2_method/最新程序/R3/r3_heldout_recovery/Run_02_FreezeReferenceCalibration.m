function Result = Run_02_FreezeReferenceCalibration()
%RUN_02_FREEZEREFERENCECALIBRATION Freeze sigma_V and the main low-speed C0.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
outDir=fullfile(root,'output','02_reference_calibration');if ~exist(outDir,'dir'),mkdir(outDir);end
[sigmaV,noiseInfo]=r3_reference_noise_std(R3,P.referenceSnrDb);
[C0,low0]=r3_freeze_reference_calibration(R3,sigmaV,P.mainLowSeed);
Result=struct('P',P,'sigmaV',sigmaV,'noiseInfo',noiseInfo,'C0',C0,...
    'low0',low0,'calibration_gap_states_mm',P.gCalMm,...
    'held_out_operating_states_mm',P.gOperatingMm,'outputDir',outDir);
save(fullfile(outDir,'r3_reference_calibration_C0.mat'),'Result','-v7.3');
fprintf('Frozen C0: g0 truth %.3f mm, estimate %.6f mm, sigma %.6g V\n',...
    P.gReferenceMm,C0.pathCal.g0,sigmaV);
end
