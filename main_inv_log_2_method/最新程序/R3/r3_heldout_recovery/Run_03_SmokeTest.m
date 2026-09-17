function Result = Run_03_SmokeTest(maxCells)
%RUN_03_SMOKETEST Sparse engineering validation; not paper evidence.

if nargin<1||isempty(maxCells),maxCells=Inf;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
outDir=fullfile(root,'output','03_smoke');if ~exist(outDir,'dir'),mkdir(outDir);end
[sigmaV,noiseInfo]=r3_reference_noise_std(R3,P.referenceSnrDb);
[C0,~]=r3_freeze_reference_calibration(R3,sigmaV,P.mainLowSeed);
thresholds=struct();rows=[];cellCount=0;
for g=P.smokeGapsMm
    for A=P.smokeAmplitudesMm
        cellCount=cellCount+1;if cellCount>maxCells,break;end
        for rep=P.smokeReplicateIndices
            q=r3_fit_rows(R3,C0,g,A,P.phaseValuesRad(rep),sigmaV,...
                P.mainHighSeeds(rep),thresholds);
            rows=[rows;q(:)]; %#ok<AGROW>
        end
    end
    if cellCount>=maxCells,break;end
end
Detail=struct2table(rows);writetable(Detail,fullfile(outDir,'smoke_detail.csv'));
Result=struct('P',P,'noiseInfo',noiseInfo,'Detail',Detail,'outputDir',outDir,...
    'warning',"Engineering smoke test only; detection thresholds are not applied.");
save(fullfile(outDir,'smoke_result.mat'),'Result','-v7.3');disp(Detail(:,...
    {'method','g_truth_mm','A_true_mm','g_est_mm','A_est_mm','eo_est','vib_success'}));
end
