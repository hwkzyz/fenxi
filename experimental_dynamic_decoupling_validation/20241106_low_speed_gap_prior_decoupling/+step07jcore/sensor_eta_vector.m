function etaVec = sensor_eta_vector(sensorKey, sensorIds, etaByPosition)
%SENSOR_ETA_VECTOR Expand per-sensor eta to one value per observation.

sensorKey = sensorKey(:);
sensorIds = sensorIds(:).';
etaByPosition = etaByPosition(:).';
etaVec = zeros(size(sensorKey));

if isempty(sensorKey) || isempty(sensorIds) || isempty(etaByPosition)
    return;
end

nSensor = min(numel(sensorIds), numel(etaByPosition));
sensorIds = sensorIds(1:nSensor);
etaByPosition = etaByPosition(1:nSensor);

finiteKey = sensorKey(isfinite(sensorKey));
useLocalIndex = ~isempty(finiteKey) && ...
    all(abs(finiteKey - round(finiteKey)) < 1e-9) && ...
    all(finiteKey >= 1 & finiteKey <= nSensor) && ...
    ~all(ismember(finiteKey(:).', sensorIds));

if useLocalIndex
    for is = 1:nSensor
        etaVec(sensorKey == is) = etaByPosition(is);
    end
else
    for is = 1:nSensor
        etaVec(sensorKey == sensorIds(is)) = etaByPosition(is);
    end
end

etaVec(~isfinite(etaVec)) = 0;
end
