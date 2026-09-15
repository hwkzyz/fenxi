function out=Run_R5_SensorConditionedDynamic_20250527(sidecarFile,foundationFile,outputFile)
rootDir=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(rootDir); addpath(fullfile(rootDir,'common'));
cfg=Config_20250527(); if nargin<1||isempty(sidecarFile),sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');end
if nargin<2||isempty(foundationFile),foundationFile=fullfile(cfg.paths.preparedInputs,'foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat');end
if nargin<3||isempty(outputFile),outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat');end
out=R5_Run_SensorConditionedWindowIdentification(cfg,@R5_Build_SensorConditionedDynamicModels_20250527,sidecarFile,cfg.files.lowSpeedTemplate,foundationFile,outputFile);
end
