function out = run_R5_20251222(varargin)
%RUN_R5_20251222 Run the packaged 20251222 R5 dynamic identification.
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'adapters'),'-begin');
addpath(fullfile(rootDir,'core'),'-begin');
verify_R5_Installation('case','20251222');
sidecar = fullfile(rootDir,'inputs','20251222','frozen','R5_SensorConditionedSidecar_B1_S123.mat');
foundation = fullfile(rootDir,'inputs','20251222','foundation', ...
    'Result_Step05_NoEtaVPTopKDirectTemplate_B1_S123_20251222.mat');
outFile = fullfile(rootDir,'results','20251222','r5_sensor_conditioned_dynamic_B1_result.mat');
out = Run_R5_SensorConditionedDynamic_20251222(sidecar,foundation,outFile,1);
end
