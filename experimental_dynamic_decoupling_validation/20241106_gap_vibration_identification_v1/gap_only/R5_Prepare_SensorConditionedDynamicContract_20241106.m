function C = R5_Prepare_SensorConditionedDynamicContract_20241106(sidecarFile, foundationResultFile, templateFile)
% Freeze the static-to-dynamic contract without running identification.
cfg=Config_20241106();
if nargin<1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin<2 || isempty(foundationResultFile), foundationResultFile=''; end
if nargin<3 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
R5_Validate_SensorConditionedSidecar_20241106(sidecarFile);
models=R5_Build_SensorConditionedDynamicModels_20241106(sidecarFile,templateFile,cfg.case.targetBlade);
assert(isequal(sort([models.sensorId]),sort(cfg.case.gapSensors)),'R5:SensorRoleMismatch');
windowCount=NaN;
if ~isempty(foundationResultFile)
    S=load(foundationResultFile,'Result'); assert(isfield(S.Result,'WindowResult'),'R5:MissingWindows');
    windowCount=numel(S.Result.WindowResult);
end
C=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_CONTRACT_V1','dataset',cfg.dataset,...
    'targetBlade',cfg.case.targetBlade,'sensorIds',[models.sensorId],'windowCount',windowCount,...
    'sidecarFile',sidecarFile,'templateFile',templateFile,...
    'foundationResultFile',foundationResultFile,'freeParameters',{{'EO','A','phase','dx','dg_s'}},...
    'forbiddenDynamicParameters',{{'deltaMu','deltaTau'}},'status','prepared');
save(fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_contract.mat'),'C','-v7.3');
end
