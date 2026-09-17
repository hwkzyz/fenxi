function Result = Run_04_CalibrateDetectionThreshold(maxReplicates)
%RUN_04_CALIBRATEDETECTIONTHRESHOLD Conditional null thresholds using frozen C0.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
if nargin<1||isempty(maxReplicates),maxReplicates=P.detectionCalibrationCount;end
maxReplicates=min(maxReplicates,P.detectionCalibrationCount);
calFile=fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat');
if ~exist(calFile,'file'),error('r3:MissingC0','Run Run_02_FreezeReferenceCalibration first.');end
S=load(calFile,'Result');C0=S.Result.C0;sigmaV=S.Result.sigmaV;
outDir=fullfile(root,'output','04_detection_threshold');if ~exist(outDir,'dir'),mkdir(outDir);end
checkpoint=fullfile(outDir,'detection_detail.csv');
if exist(checkpoint,'file'),Detail=readtable(checkpoint,'TextType','string');else,Detail=table();end
for rep=1:maxReplicates
    if ~isempty(Detail)&&nnz(Detail.high_seed==P.detectionHighSeeds(rep))==3,continue;end
    q=r3_fit_rows(R3,C0,P.gReferenceMm,0,0,sigmaV,P.detectionHighSeeds(rep),struct());
    Detail=[Detail;struct2table(q)]; %#ok<AGROW>
    if mod(rep,10)==0||rep==maxReplicates,writetable(Detail,checkpoint);end
end
names=["fixed","adaptive","state_matched"];Thresholds=struct();
for i=1:3
    q=Detail.method==names(i);Thresholds.(char(names(i)))=quantile(Detail.A_est_mm(q),P.detectionQuantile);
end
formal=maxReplicates==P.detectionCalibrationCount&&...
    numel(unique(Detail.high_seed))>=P.detectionCalibrationCount;
Result=struct('P',P,'Thresholds',Thresholds,'Detail',Detail,'formal',formal,...
    'replicate_count',numel(unique(Detail.high_seed)),'outputDir',outDir);
save(fullfile(outDir,'detection_thresholds.mat'),'Result','-v7.3');disp(Thresholds);
end
