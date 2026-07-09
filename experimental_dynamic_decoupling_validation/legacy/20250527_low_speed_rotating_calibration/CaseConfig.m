function C = CaseConfig()
%CASECONFIG Single-point case configuration for 20250527 rotating calibration.
%
% Edit only C.bladeId and C.sensorIds below when switching the target case.
% Optional temporary overrides:
%   BLADE_CASE_BLADE_ID
%   BLADE_CASE_SENSOR_IDS   e.g. '1 3 6'

C = struct();
C.dataset = '20250527';
C.bladeId = 1;
C.sensorIds = [1 3 6];

[overrideBladeId, hasBladeOverride] = parse_scalar_override_local('BLADE_CASE_BLADE_ID');
if hasBladeOverride
    C.bladeId = overrideBladeId;
end
[overrideSensorIds, hasSensorOverride] = parse_sensor_override_local('BLADE_CASE_SENSOR_IDS');
if hasSensorOverride
    C.sensorIds = overrideSensorIds;
end

C.sensorIds = reshape(C.sensorIds, 1, []);
C.sensorTag = ['S', sprintf('%d', C.sensorIds)];
C.caseTag = sprintf('B%d_%s', C.bladeId, C.sensorTag);
end

function [value, hasOverride] = parse_scalar_override_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    value = NaN;
    hasOverride = false;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 1 || abs(value - round(value)) > eps(value)
    error('%s must be a positive integer.', name);
end
value = round(value);
hasOverride = true;
end

function [values, hasOverride] = parse_sensor_override_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    values = [];
    hasOverride = false;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values) || any(values < 1)
    error('%s must contain positive integer sensor IDs separated by spaces.', name);
end
hasOverride = true;
end
