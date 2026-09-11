function out=Run_R5_SensorConditionedDynamic_20241106(sidecarFile,foundationFile,outputFile)
rootDir=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(rootDir);
cfg=Config_20241106(); if nargin<1||isempty(sidecarFile),sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');end
if nargin<2||isempty(foundationFile),error('R5:MissingFoundation','Pass the 20241106 foundation WindowResult file.');end
if nargin<3||isempty(outputFile),outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat');end
templateFile=fullfile(cfg.paths.results,'r5_low_speed_template_bank_adapter.mat');
assert(isfile(templateFile),'R5:MissingTemplateAdapter','Run Run_R5_Preparation_20241106 first.');
out=R5_Run_SensorConditionedWindowIdentification(cfg,@R5_Build_SensorConditionedDynamicModels_20241106,sidecarFile,templateFile,foundationFile,outputFile);
end
