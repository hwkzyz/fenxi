%% Summarize low-speed donor selection and completed R4 high-speed reruns
clear; clc;
thisDir=fileparts(mfilename('fullpath')); rootDir=fileparts(fileparts(thisDir)); outDir=fullfile(thisDir,'results');
calDir=fullfile(rootDir,'results','prepared','calibration'); gapDir=fullfile(rootDir,'results','gap_aware');

L=table();
for b=1:6
    f=fullfile(calDir,sprintf('Step06I_OffsetTiltShared_GapLibrary_Calibration_20251222_B5_S123_manifold_donorB%d_nob_20260901.csv',b));
    if ~isfile(f), continue, end
    q=readtable(f); q.donorBlade=repmat(b,height(q),1); L=[L;q]; %#ok<AGROW>
end
if ~isempty(L)
    L=movevars(L,'donorBlade','Before',1);
    writetable(L,fullfile(outDir,'DonorSurface_LowSpeedCalibration_20251222.csv'));
    S=groupsummary(L,'donorBlade',{'mean','max'},{'lowFitRmseMv'});
    writetable(S,fullfile(outDir,'DonorSurface_LowSpeedCalibration_Summary_20251222.csv')); disp(S);
end

ids={'formal_shared_no_b','per_blade_B1','per_blade_B2'};
files={ ...
    'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_nob_formal_B5_20260830.csv', ...
    'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_manifold_donorB1_nob_R4_full_20260901.csv', ...
    'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_manifold_donorB2_nob_R4_full_20260901.csv'};
rows=cell(0,9); ref=[];
for i=1:numel(files)
    f=fullfile(gapDir,files{i}); if ~isfile(f), continue, end
    q=readtable(f); if i==1, ref=q; end
    nEO18=sum(q.gap_EO==18); rows(end+1,:)={ids{i},height(q),nEO18,mean(q.gap_frequency_hz), ... %#ok<AGROW>
        std(q.gap_frequency_hz),mean(q.gap_amplitude_mm),mean(q.gap_rmse_mV), ...
        mean(q.gap_mean_delta_gap_mm),string(f)};
end
H=cell2table(rows,'VariableNames',{'method','nWindows','nEO18','meanFrequencyHz','stdFrequencyHz', ...
    'meanAmplitudeMm','meanRMSEmV','meanDeltaGapMm','file'});
writetable(H,fullfile(outDir,'DonorSurface_HighSpeedV4_Summary_20251222.csv')); disp(H(:,1:8));

if ~isempty(ref)
    rows=cell(0,7);
    for i=2:numel(files)
        f=fullfile(gapDir,files{i}); if ~isfile(f), continue, end
        q=readtable(f); n=min(height(q),height(ref)); q=q(1:n,:); r=ref(1:n,:);
        rows(end+1,:)={ids{i},sum(q.gap_EO==r.gap_EO),n,mean(abs(q.gap_frequency_hz-r.gap_frequency_hz)), ... %#ok<AGROW>
            mean(abs(q.gap_amplitude_mm-r.gap_amplitude_mm)),mean(abs(q.gap_rmse_mV-r.gap_rmse_mV)), ...
            mean(abs(q.gap_mean_delta_gap_mm-r.gap_mean_delta_gap_mm))};
    end
    D=cell2table(rows,'VariableNames',{'method','nEOMatches','nWindows','meanAbsFrequencyDiffHz', ...
        'meanAbsAmplitudeDiffMm','meanAbsRMSEDiffmV','meanAbsDeltaGapDiffMm'});
    writetable(D,fullfile(outDir,'DonorSurface_HighSpeedV4_VsFormal_20251222.csv')); disp(D);
end
