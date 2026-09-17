function out = run_R5_20250527(varargin)
%RUN_R5_20250527 Run the packaged 20250527 R5 dynamic identification.
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'adapters'),'-begin');
addpath(fullfile(rootDir,'core'),'-begin');
verify_R5_Installation('case','20250527');
sidecar = fullfile(rootDir,'inputs','20250527','frozen','r5_sensor_conditioned_sidecar.mat');
foundation = fullfile(rootDir,'inputs','20250527','foundation', ...
    'Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat');
outFile = fullfile(rootDir,'results','20250527','r5_sensor_conditioned_dynamic_result.mat');
out = Run_R5_SensorConditionedDynamic_20250527(sidecar,foundation,outFile);
end
