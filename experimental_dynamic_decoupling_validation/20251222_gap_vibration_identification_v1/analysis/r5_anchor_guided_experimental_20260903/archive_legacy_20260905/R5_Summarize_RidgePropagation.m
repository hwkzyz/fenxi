function T=R5_Summarize_RidgePropagation(manifestFile,resultsDir,outputFile)
% Summarize the one-window ridge-to-dynamics propagation smoke test.
if nargin<1||isempty(manifestFile), manifestFile=fullfile(fileparts(mfilename('fullpath')),'results','ridge_sidecars','R5_RidgeSidecar_Manifest.csv'); end
if nargin<2||isempty(resultsDir), resultsDir=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'results','gap_aware'); end
if nargin<3||isempty(outputFile), outputFile=fullfile(fileparts(mfilename('fullpath')),'results','ridge_sidecars','R5_RidgePropagation_OneWindow.csv'); end
M=readtable(manifestFile,'TextType','string'); rows=table();
for k=1:height(M)
    b=M.target_blade(k); p=M.ridge_index(k); f=fullfile(resultsDir,sprintf('Main_GapAware_VPTopK_FullWave_Summary_20251222_B%d_S123_r5_ridge_b%d_p%d_w1_20260904.csv',b,b,p));
    assert(isfile(f),'R5:RidgeResultMissing','Missing result %s.',f); s=readtable(f);
    rows=[rows; table(b,p,M.latent_coordinate(k),M.joint_rmse_mv(k),M.delta_j_mv(k),M.gap_s1_mm(k),M.gap_s2_mm(k),M.gap_s3_mm(k),s.dominant_EO(1),s.mean_frequency_hz(1),s.mean_rmse_mV(1),string(f), ...
        'VariableNames',{'target_blade','ridge_index','latent_coordinate','joint_rmse_mv','delta_j_mv','gap_s1_mm','gap_s2_mm','gap_s3_mm','dominant_EO','frequency_hz','rmse_mv','result_file'})]; %#ok<AGROW>
end
if ~isfolder(fileparts(outputFile)), mkdir(fileparts(outputFile)); end
writetable(rows,outputFile); T=rows;
end
