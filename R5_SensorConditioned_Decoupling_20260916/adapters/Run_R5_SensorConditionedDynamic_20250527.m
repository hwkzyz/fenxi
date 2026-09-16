function out=Run_R5_SensorConditionedDynamic_20250527(sidecarFile,foundationFile,outputFile)
%RUN_R5_SENSORCONDITIONEDDYNAMIC_20250527 Package-local formal entry point.
rootDir=fileparts(mfilename('fullpath'));
addpath(rootDir,'-begin');
addpath(fullfile(rootDir,'..','core'),'-begin');
cfg=Config_20250527();
if nargin<1||isempty(sidecarFile),sidecarFile=fullfile(rootDir,'..','inputs','20250527','frozen','r5_sensor_conditioned_sidecar.mat');end
if nargin<2||isempty(foundationFile),foundationFile=fullfile(rootDir,'..','inputs','20250527','foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat');end
if nargin<3||isempty(outputFile),outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat');end
assert(isfile(sidecarFile),'R5:MissingSidecar','Missing package sidecar: %s',sidecarFile);
assert(isfile(foundationFile),'R5:MissingFoundation','Missing package foundation: %s',foundationFile);
assert(isfile(cfg.files.lowSpeedTemplate),'R5:MissingTemplate','Missing package template: %s',cfg.files.lowSpeedTemplate);
out=R5_Run_SensorConditionedWindowIdentification(cfg,@R5_Build_SensorConditionedDynamicModels_20250527,sidecarFile,cfg.files.lowSpeedTemplate,foundationFile,outputFile);
end
