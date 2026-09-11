function G = R5_ImportLegacyP2Registration_20250527(csvFile, outputFile)
%R5_IMPORTLEGACYP2REGISTRATION_20250527 Freeze the accepted old P2 registration.
cfg = Config_20250527();
if nargin < 1 || isempty(csvFile)
    csvFile = fullfile(cfg.paths.calibration, 'NestedPlatformModels_20250527.csv');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results, 'r5_legacy_p2_registration.mat');
end
T = readtable(csvFile, 'TextType', 'string');
assert(all(ismember({'sensorId','model','tauMm','xScale','voltageGain'}, T.Properties.VariableNames)), ...
    'R5:LegacyP2Schema', 'Legacy registration CSV lacks required fields.');
T = T(strcmpi(string(T.model), 'P2_registration'), :);
assert(all(ismember(cfg.case.gapSensors, T.sensorId)), ...
    'R5:LegacyP2Sensors', 'Legacy P2 registration does not cover every gap sensor.');
registration = repmat(struct('sensor_id',NaN,'tau_mm',NaN,'x_scale',NaN, ...
    'voltage_gain',NaN,'voltage_offset_mv',0,'rmse_mv',NaN,'point_count',NaN, ...
    'b2_gap_mm',NaN), height(T), 1);
for k = 1:height(T)
    registration(k).sensor_id = T.sensorId(k);
    registration(k).tau_mm = T.tauMm(k);
    registration(k).x_scale = T.xScale(k);
    registration(k).voltage_gain = T.voltageGain(k);
    registration(k).rmse_mv = T.rmseMv(k);
end
G = struct('schema','R5_LEGACY_P2_FROZEN_REGISTRATION_V1', ...
    'source_csv',csvFile,'anchor_blade',NaN,'registration',registration, ...
    'summary',T,'note','P2 coordinate registration is frozen; no new B2 affine refit is performed.');
if ~exist(fileparts(outputFile), 'dir'), mkdir(fileparts(outputFile)); end
save(outputFile, 'G', '-v7.3');
end
