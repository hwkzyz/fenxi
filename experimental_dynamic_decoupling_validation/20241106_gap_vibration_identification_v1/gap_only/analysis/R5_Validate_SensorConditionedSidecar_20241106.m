function report = R5_Validate_SensorConditionedSidecar_20241106(sidecarFile)
% Validate the explicit per-sensor dynamic forward-model contract.
cfg = Config_20241106();
if nargin < 1 || isempty(sidecarFile)
    sidecarFile = fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');
end
S = load(sidecarFile,'SensorConditionedLibrary');
L = S.SensorConditionedLibrary;
assert(strcmp(L.schema,'R5_SENSOR_CONDITIONED_SIDECAR_V1'),'R5:Schema');
assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
assert(~isempty(L.sensor),'R5:EmptySidecar');
rows = table();
for i = 1:numel(L.sensor)
    s = L.sensor{i};
    assert(~isempty(s.registration),'R5:MissingRegistration');
    assert(numel(s.state)==height(s.rows),'R5:StateRowMismatch');
    for j = 1:numel(s.state)
        B = s.state(j).B;
        assert(isequal(size(B),[numel(L.x_mm),3]),'R5:CoefficientShape');
        assert(all(isfinite(B(:))),'R5:NonfiniteCoefficient');
        rows = [rows; table(double(s.sensor_id),double(s.state(j).blade_id), ...
            double(s.state(j).z),double(s.state(j).gap_mm),double(s.state(j).rmse_mv), ...
            'VariableNames',{'sensor_id','blade_id','z','gap_mm','rmse_mv'})]; %#ok<AGROW>
    end
end
assert(all(isfinite(rows.z)) && all(isfinite(rows.gap_mm)),'R5:NonfiniteState');
report = struct('file',sidecarFile,'schema',L.schema,'latent_mode',L.latent_mode, ...
    'sensor_ids',unique(rows.sensor_id).','state_count',height(rows),'states',rows, ...
    'status','pass');
fprintf('Sensor-conditioned sidecar validated: %s (%d states, sensors %s)\n', ...
    sidecarFile,report.state_count,mat2str(report.sensor_ids));
end
