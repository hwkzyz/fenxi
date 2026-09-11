%% 06: calibrate the all-blade gap library for 20241106
% This is the calibration-first bank builder. It consumes the all-blade
% low-speed template bank from Step06A and fits gap-correction parameters
% only for capacitive sensors CH5/CH7. CH2/CH3 remain in the low-speed
% template bank as direct-only channels and are never added to this gap bank.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir, fullfile(thisDir, 'functions', 'gap_aware'), ...
    fullfile(thisDir, 'functions', 'utilities'));
packageCfg = Config_20241106();
outDir = packageCfg.paths.calibrationRuntime;
figDir = fullfile(outDir, 'figures_step06i_all_blade_gap_library_20241106');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

P = ProjectionFlow_Config_20241106();
responseFile = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_20241106.mat');
templateBankFile = fullfile(outDir, P.calibration.runtimeTemplateBankFile);
if ~isfile(responseFile)
    error('Run Step05I first. Missing file: %s', responseFile);
end
if ~isfile(templateBankFile)
    fprintf('Low-speed template bank is missing. Running program 05 first...\n');
    run(fullfile(thisDir, 'Main05_Build_LowSpeed_TemplateBank_20241106.m'));
end
if ~isfile(templateBankFile)
    error('Missing low-speed template bank: %s', templateBankFile);
end

Srf = load(responseFile, 'OffsetTiltResponseSurface');
responseSurface = Srf.OffsetTiltResponseSurface;
Sb = load(templateBankFile, 'LowSpeedTemplateBank');
LowSpeedTemplateBank = Sb.LowSpeedTemplateBank;

cfg = struct();
cfg.gBoundsMm = [min(responseSurface.gTrainMm(:)), max(responseSurface.gTrainMm(:))];
cfg.tauBoundsMm = [-0.80, 0.80];
cfg.xScaleBounds = [0.65, 1.45];
cfg.muBounds = [-0.35, 0.35];
cfg.lowFitMode = 'affine';
lowFitModeOverride = lower(strtrim(getenv('STEP06I_LOW_FIT_MODE')));
if ~isempty(lowFitModeOverride)
    if ~ismember(lowFitModeOverride, {'affine','slope_only'})
        error('STEP06I_LOW_FIT_MODE must be affine or slope_only.');
    end
    cfg.lowFitMode = lowFitModeOverride;
end
cfg.overshootPenaltyMv = 800;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 450, ...
    'MaxFunEvals', 1600, 'TolX', 1e-6, 'TolFun', 1e-6);

% Calibrate exactly the blades represented by the formal template bank.
% The current bundled common-window source is the validated B4 case; using
% P.calibration.buildBladeIds here would request absent, fabricated entries.
bladeIds = parse_int_env_local('STEP06I_BANK_BLADE_IDS', LowSpeedTemplateBank.bladeIds);
gapSensors = parse_int_env_local('STEP06I_BANK_GAP_SENSORS', P.calibration.gapSensorIds);
assert_capacitive_gap_sensors_local(gapSensors);

bankDir = fullfile(outDir, P.calibration.gapBankDirName);
if exist(bankDir, 'dir') ~= 7
    mkdir(bankDir);
end

fprintf('\n=== Step06I: all-blade gap calibration bank 20241106 ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Template bank    : %s\n', templateBankFile);
fprintf('Blades: %s | gap sensors: %s | low-speed fit: %s\n', mat2str(bladeIds), mat2str(gapSensors), cfg.lowFitMode);

entry = repmat(struct('bladeId', NaN, 'sensorId', NaN, 'sensor', struct(), ...
    'gapCalibFile', "", 'templateFile', "", 'sourceTemplateBank', string(templateBankFile)), ...
    numel(bladeIds) * numel(gapSensors), 1);
rows = cell(numel(entry), 1);
k = 0;

for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    sensorCorr = repmat(empty_sensor_corr_local(), numel(gapSensors), 1);
    tableRows = cell(numel(gapSensors), 1);
    Template = build_template_from_bank_local(LowSpeedTemplateBank, bladeId, gapSensors);

    for is = 1:numel(gapSensors)
        sid = gapSensors(is);
        Tsen = get_template_sensor_local(Template, sid);
        validate_template_sensor_fields_local(Tsen, sid);
        fit = fit_low_template_local(Tsen, responseSurface, cfg);
        sensorCorr(is) = fit;
        tableRows{is} = calibration_row_local(bladeId, fit);

        k = k + 1;
        gapCalibFile = fullfile(bankDir, sprintf('%s_20241106_B%d_S%d.mat', ...
            P.calibration.gapEntryPrefix, bladeId, sid));
        GapCalib = struct(); %#ok<NASGU>
        GapCalib.dataset = P.dataset;
        GapCalib.bladeId = bladeId;
        GapCalib.sensorId = sid;
        GapCalib.sensor = fit;
        GapCalib.responseFile = responseFile;
        GapCalib.templateEntry = find_template_bank_entry_local(LowSpeedTemplateBank, bladeId, sid);
        GapCalib.sourceTemplateBank = templateBankFile;
        GapCalib.method = 'all_blade_offset_tilt_shared_gap_calibration';
        save(gapCalibFile, 'GapCalib', '-v7.3');

        entry(k).bladeId = bladeId;
        entry(k).sensorId = sid;
        entry(k).sensor = fit;
        entry(k).gapCalibFile = string(gapCalibFile);
        entry(k).templateFile = string(GapCalib.templateEntry.templateFile);
        rows{k} = table(string(P.dataset), bladeId, sid, string(gapCalibFile), ...
            string(GapCalib.templateEntry.templateFile), fit.g0Mm, fit.tauMm, ...
            fit.xScale, fit.muGapPerXMm, fit.lowFitRmseMv, ...
            'VariableNames', {'dataset','bladeId','sensorId','gapCalibFile', ...
            'templateFile','g0Mm','tauMm','xScale','muGapPerXMm','lowFitRmseMv'});

        fprintf('  B%d CH%d: g0 %.4f, tau %.4f, k %.4f, mu %.4f, RMSE %.2f mV\n', ...
            bladeId, sid, fit.g0Mm, fit.tauMm, fit.xScale, fit.muGapPerXMm, fit.lowFitRmseMv);
    end

    calibrationTable = vertcat(tableRows{:}); %#ok<NASGU>
    CorrectedGapLibrary = build_corrected_gap_library_local(P, bladeId, gapSensors, ...
        Template, sensorCorr, responseSurface, responseFile, templateBankFile); %#ok<NASGU>
    sensorTag = ['S', sprintf('%d', gapSensors)];
    combinedFile = fullfile(outDir, sprintf( ...
        'Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', bladeId, sensorTag));
    save(combinedFile, 'CorrectedGapLibrary', 'Template', 'calibrationTable', ...
        'responseFile', 'templateBankFile', '-v7.3');
end

indexTable = vertcat(rows{:}); %#ok<NASGU>
GapCalibrationBank = struct(); %#ok<NASGU>
GapCalibrationBank.dataset = P.dataset;
GapCalibrationBank.method = 'all_blade_offset_tilt_shared_gap_calibration_bank';
GapCalibrationBank.description = ['Gap-correction bank for capacitive probes only. ' ...
    'Direct-only probes CH2/CH3 are intentionally excluded and remain in LowSpeedTemplateBank.'];
GapCalibrationBank.bladeIds = bladeIds(:).';
GapCalibrationBank.gapSensors = gapSensors(:).';
GapCalibrationBank.directOnlySensors = P.calibration.directOnlySensorIds(:).';
GapCalibrationBank.lowSpeedTemplateSensors = P.calibration.lowSpeedTemplateSensors(:).';
GapCalibrationBank.responseFile = responseFile;
GapCalibrationBank.templateBankFile = templateBankFile;
GapCalibrationBank.responseSurface = responseSurface;
GapCalibrationBank.entry = entry;
GapCalibrationBank.indexTable = indexTable;
GapCalibrationBank.cfg = P;

bankFile = fullfile(outDir, P.calibration.runtimeGapBankFile);
indexFile = fullfile(bankDir, 'GapCalibrationBank_Index_20241106.csv');
save(bankFile, 'GapCalibrationBank', 'indexTable', '-v7.3');
writetable(indexTable, indexFile);

fprintf('\nStep06I all-blade bank complete.\n');
fprintf('Bank : %s\n', bankFile);
fprintf('Index: %s\n', indexFile);

function item = empty_sensor_corr_local()
item = struct('sensorId', NaN, 'g0Mm', NaN, 'tauMm', NaN, ...
    'xScale', NaN, 'muGapPerXMm', NaN, 'tiltAngleDeg', NaN, ...
    'voltageGain', NaN, 'voltageOffsetMv', NaN, 'lowFitRmseMv', NaN, ...
    'overshootRmseMm', NaN, 'pointCount', NaN, 'x', [], 'vLowMv', [], ...
    'vFitMv', [], 'xLib', [], 'gEff', [], 'residualMv', [], ...
    'etaMedianMm', NaN, 'etaIqrMm', NaN, 'etaLimitMm', NaN);
end

function Template = build_template_from_bank_local(bank, bladeId, sensorIds)
Template = struct();
Template.dataset = bank.dataset;
Template.targetBlade = bladeId;
Template.analysisSensors = sensorIds(:).';
Template.SensorIDs = sensorIds(:).';
Template.SensorTag = ['S', sprintf('%d', sensorIds)];
Template.Sensor = struct([]);
for sid = sensorIds(:).'
    e = find_template_bank_entry_local(bank, bladeId, sid);
    sensor = e.sensorTemplate;
    Template.Sensor = append_sensor_local(Template.Sensor, sensor);
end
end

function e = find_template_bank_entry_local(bank, bladeId, sid)
idx = find([bank.entry.bladeId] == bladeId & [bank.entry.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('LowSpeedTemplateBank does not contain B%d CH%d.', bladeId, sid);
end
e = bank.entry(idx);
end

function sensorArray = append_sensor_local(sensorArray, sensor)
if isempty(sensorArray)
    sensorArray = sensor;
    return;
end
allFields = unique([fieldnames(sensorArray); fieldnames(sensor)], 'stable');
for i = 1:numel(sensorArray)
    sensorArray(i) = ensure_fields_local(sensorArray(i), allFields);
end
sensor = ensure_fields_local(sensor, allFields);
sensorArray(end + 1) = orderfields(sensor, fieldnames(sensorArray)); %#ok<AGROW>
end

function s = ensure_fields_local(s, fields)
for i = 1:numel(fields)
    if ~isfield(s, fields{i})
        s.(fields{i}) = [];
    end
end
end

function row = calibration_row_local(bladeId, fit)
row = table(bladeId, fit.sensorId, fit.g0Mm, fit.tauMm, fit.xScale, ...
    fit.muGapPerXMm, fit.tiltAngleDeg, fit.voltageGain, fit.voltageOffsetMv, ...
    fit.lowFitRmseMv, fit.overshootRmseMm, fit.pointCount, fit.etaMedianMm, ...
    fit.etaIqrMm, fit.etaLimitMm, ...
    'VariableNames', {'bladeId','sensorId','g0Mm','tauMm','xScale', ...
    'muGapPerXMm','tiltAngleDeg','voltageGain','voltageOffsetMv', ...
    'lowFitRmseMv','overshootRmseMm','pointCount','etaMedianMm', ...
    'etaIqrMm','etaLimitMm'});
end

function Lib = build_corrected_gap_library_local(P, bladeId, gapSensors, Template, sensorCorr, responseSurface, responseFile, templateBankFile)
Lib = struct();
Lib.dataset = P.dataset;
Lib.method = 'all_blade_bank_offset_tilt_shared_gap_library';
Lib.description = ['Per-blade view of the all-blade gap bank. ' ...
    'Only capacitive sensors are included.'];
Lib.responseFile = responseFile;
Lib.templateBankFile = templateBankFile;
Lib.templateFile = string(templateBankFile);
Lib.templateSourceMode = 'low_speed_template_bank';
Lib.templateSourceFile = string(templateBankFile);
Lib.targetBlade = bladeId;
Lib.analysisSensors = gapSensors(:).';
Lib.gapSensors = gapSensors(:).';
Lib.directOnlySensors = P.calibration.directOnlySensorIds(:).';
Lib.responseSurface = responseSurface;
Lib.lowSpeedTemplate = Template;
Lib.sensor = sensorCorr;
Lib.formula = 'T_low(x)+a*(F_raw(g+mu*(x-tau),k*(x-tau))-F_raw(g0+mu*(x-tau),k*(x-tau)))';
Lib.cfg = P;
end

function fit = fit_low_template_local(Tsen, responseSurface, cfg)
x = Tsen.x_grid(:);
vLow = Tsen.v_grid(:) - Tsen.baseline;
mask = isfinite(x) & isfinite(vLow);
x = x(mask);
v = vLow(mask);
if numel(x) < 50
    error('Low-speed template for sensor %d has too few points.', Tsen.sensor_id);
end
gInit = median(responseSurface.gTrainMm(:), 'omitnan');
p0 = [gInit, 0, 1, 0];
obj = @(p) low_template_objective_local(p, x, v, responseSurface, cfg);
p = fminsearch(obj, p0, cfg.fminOptions);
p = clip_params_local(p, cfg);
[model, over] = eval_raw_tilt_path_local(responseSurface, p(1), p(2), p(3), p(4), x);
if strcmpi(cfg.lowFitMode, 'slope_only')
    [gain, offset, vFit] = solve_gain_local(model, v);
else
    [gain, offset, vFit] = solve_gain_offset_local(model, v);
end
residual = v - vFit;
fit = empty_sensor_corr_local();
fit.sensorId = Tsen.sensor_id;
fit.g0Mm = p(1);
fit.tauMm = p(2);
fit.xScale = p(3);
fit.muGapPerXMm = p(4);
fit.tiltAngleDeg = atan(p(4)) * 180 / pi;
fit.voltageGain = gain;
fit.voltageOffsetMv = offset;
fit.lowFitRmseMv = sqrt(mean(residual.^2, 'omitnan'));
fit.overshootRmseMm = sqrt(mean(over.^2, 'omitnan'));
fit.pointCount = numel(x);
fit.x = x;
fit.vLowMv = v;
fit.vFitMv = vFit;
fit.xLib = p(3) .* (x - p(2));
fit.gEff = p(1) + p(4) .* (x - p(2));
fit.residualMv = residual;
fit.etaMedianMm = get_eta_scalar_local(Tsen, 'eta_median_mm');
fit.etaIqrMm = get_eta_scalar_local(Tsen, 'eta_iqr_mm');
fit.etaLimitMm = get_eta_scalar_local(Tsen, 'eta_limit_mm');
end

function value = low_template_objective_local(p, x, v, responseSurface, cfg)
if any(~isfinite(p))
    value = inf;
    return;
end
pClip = clip_params_local(p, cfg);
penalty = sum((p - pClip).^2) * cfg.overshootPenaltyMv;
[model, over] = eval_raw_tilt_path_local(responseSurface, pClip(1), pClip(2), pClip(3), pClip(4), x);
if strcmpi(cfg.lowFitMode, 'slope_only')
    [~, ~, vFit] = solve_gain_local(model, v);
else
    [~, ~, vFit] = solve_gain_offset_local(model, v);
end
residual = v - vFit;
value = sqrt(mean(residual.^2, 'omitnan')) + penalty + ...
    cfg.overshootPenaltyMv * sqrt(mean(over.^2, 'omitnan'));
end

function p = clip_params_local(p, cfg)
p(1) = min(max(p(1), cfg.gBoundsMm(1)), cfg.gBoundsMm(2));
p(2) = min(max(p(2), cfg.tauBoundsMm(1)), cfg.tauBoundsMm(2));
p(3) = min(max(p(3), cfg.xScaleBounds(1)), cfg.xScaleBounds(2));
p(4) = min(max(p(4), cfg.muBounds(1)), cfg.muBounds(2));
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLib = k .* (xOpr - tau);
gEff = g0 + mu .* (xOpr - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
overX = max(0, xMin - xLib) + max(0, xLib - xMax);
overG = max(0, gMin - gEff) + max(0, gEff - gMax);
overshoot = hypot(overX, overG);
xEval = min(max(xLib, xMin), xMax);
gEval = min(max(gEff, gMin), gMax);
model = eval_response_surface_local(responseSurface, gEval, xEval);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function [gain, offset, yFit] = solve_gain_offset_local(model, y)
mask = isfinite(model) & isfinite(y);
if nnz(mask) < 4
    gain = 1;
    offset = 0;
    yFit = NaN(size(y));
    return;
end
A = [model(mask), ones(nnz(mask), 1)];
ab = A \ y(mask);
gain = ab(1);
offset = ab(2);
yFit = gain .* model + offset;
end

function [gain, offset, yFit] = solve_gain_local(model, y)
mask = isfinite(model) & isfinite(y);
if nnz(mask) < 4
    gain = 1; offset = 0; yFit = NaN(size(y)); return;
end
gain = model(mask) \ y(mask);
offset = 0;
yFit = gain .* model;
end

function validate_template_sensor_fields_local(Tsen, sid)
need = {'x_grid','v_grid','baseline'};
for i = 1:numel(need)
    if ~isfield(Tsen, need{i}) || isempty(Tsen.(need{i}))
        error('Template CH%d is missing %s.', sid, need{i});
    end
end
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function value = get_eta_scalar_local(Tsen, fieldName)
value = NaN;
if ~isfield(Tsen, fieldName)
    return;
end
raw = Tsen.(fieldName);
if isempty(raw)
    return;
end
raw = raw(isfinite(raw));
if isempty(raw)
    return;
end
value = raw(1);
end

function values = parse_int_env_local(name, defaultValues)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValues;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    error('%s must contain integer values.', name);
end
end

function assert_capacitive_gap_sensors_local(sensorIds)
allowed = [5 7];
if ~all(ismember(sensorIds(:).', allowed))
    error('20241106 gap calibration bank may only include capacitive sensors [5 7]. Requested: %s', ...
        mat2str(sensorIds));
end
end
