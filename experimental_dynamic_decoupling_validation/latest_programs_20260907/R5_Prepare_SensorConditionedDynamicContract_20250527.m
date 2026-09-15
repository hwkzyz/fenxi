function C = R5_Prepare_SensorConditionedDynamicContract_20250527(sidecarFile, foundationResultFile)
% Freeze the static-to-dynamic contract without running identification.
cfg=Setup_Paths_20250527();
if nargin<1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin<2 || isempty(foundationResultFile)
    foundationResultFile=fullfile(cfg.paths.preparedInputs,'foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat');
end
R5_Validate_SensorConditionedSidecar_20250527(sidecarFile);
models=R5_Build_SensorConditionedDynamicModels_20250527(sidecarFile,cfg.files.lowSpeedTemplate,cfg.case.targetBlade);
S=load(foundationResultFile,'Result'); R=S.Result;
assert(isfield(R,'WindowResult') && ~isempty(R.WindowResult),'R5:MissingWindows');
assert(all(isfinite([R.WindowResult.window_id])),'R5:WindowId');
assert(all(arrayfun(@(w)isfield(w,'CoreBundlePreview'),R.WindowResult)),'R5:MissingBundle');
ids=[models.sensorId]; assert(isequal(sort(ids),sort(cfg.case.gapSensors)),'R5:SensorRoleMismatch');
C=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_CONTRACT_V1','dataset',cfg.dataset,...
    'targetBlade',cfg.case.targetBlade,'sensorIds',ids,'windowCount',numel(R.WindowResult),...
    'sidecarFile',sidecarFile,'templateFile',cfg.files.lowSpeedTemplate,...
    'foundationResultFile',foundationResultFile,'freeParameters',{{'EO','A','phase','dx','dg_s'}},...
    'forbiddenDynamicParameters',{{'deltaMu','deltaTau'}},'status','prepared');
save(fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_contract.mat'),'C','-v7.3');
end
