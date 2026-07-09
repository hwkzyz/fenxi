function eta = fixed_sensor_eta_by_position(sensorIds, cfg)
%FIXED_SENSOR_ETA_BY_POSITION Return eta in the same order as sensorIds.

sensorIds = sensorIds(:).';
nSensor = numel(sensorIds);
eta = zeros(nSensor, 1);

if nargin < 2 || isempty(cfg) || ~isstruct(cfg)
    return;
end

if isfield(cfg, 'fixedSensorEtaMm') && ~isempty(cfg.fixedSensorEtaMm)
    raw = double(cfg.fixedSensorEtaMm(:));
elseif isfield(cfg, 'fixed_sensor_eta_mm') && ~isempty(cfg.fixed_sensor_eta_mm)
    raw = double(cfg.fixed_sensor_eta_mm(:));
elseif isfield(cfg, 'staticSensorEtaMm') && ~isempty(cfg.staticSensorEtaMm)
    raw = double(cfg.staticSensorEtaMm(:));
else
    raw = [];
end

if isempty(raw)
    return;
end

n = min(numel(raw), nSensor);
eta(1:n) = raw(1:n);
eta(~isfinite(eta)) = 0;
if nSensor >= 1
    eta(1) = 0;
end
end
