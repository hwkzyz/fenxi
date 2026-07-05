%% Step06E_AuditIdentificationEquation_20241106
% Audit whether the Step06 identification equation can distinguish EO
% candidates for the current blade/sensor/window data.
%
% This script does not change the main identification result. It reads the
% saved Step06 result and reports:
%   1. nonlinear candidate RMSE gaps from Step06,
%   2. a local linearized template-equation RMSE for selected EO values,
%   3. the weighted basis similarity between EO candidates.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

%% Parameters to tune
S06E = struct();
S06E.analysisSensors = [2 5 7];
S06E.startTimeSec = 75.0;
S06E.outputLabel = 'B4_only';
S06E.resultFileOverride = '';     % optional: audit a specific IdentificationResult_*.mat
S06E.bladeId = 4;
S06E.windowIds = [];             % [] -> all saved windows
S06E.auditEOs = [4 5 12];
S06E.referenceEO = 12;
S06E.phaseModes = {'stored_theta', 'window_time_theta', 'absolute_time_theta'};
S06E.segmentGapSec = 2e-4;        % separate points into one waveform pulse
S06E.minSegmentPoints = 20;
S06E.plotEnable = true;
S06E.saveFigures = false;

P = apply_step06e_local_options_local(P, S06E);
if ~isempty(strtrim(S06E.resultFileOverride))
    P.files.identificationResult = char(S06E.resultFileOverride);
end
require_file_local(P.files.identificationResult, 'Step06 identification result');

loaded = load(P.files.identificationResult, 'IdentificationResult');
IdentificationResult = loaded.IdentificationResult;
windowResults = IdentificationResult.WindowResult;
windowResults = windowResults([windowResults.blade_id] == S06E.bladeId);
if isempty(windowResults)
    error('No Step06 window results found for B%d.', S06E.bladeId);
end
if ~isempty(S06E.windowIds)
    keep = ismember([windowResults.window_id], S06E.windowIds);
    windowResults = windowResults(keep);
end

fprintf('\n=== Step06E: equation audit ===\n');
fprintf('Source: %s\n', P.files.identificationResult);
fprintf('Blade: B%d, sensors: %s, EO audit set: %s\n', ...
    S06E.bladeId, mat2str(S06E.analysisSensors), mat2str(S06E.auditEOs));
drawnow;

auditRows = repmat(empty_audit_row_local(S06E.auditEOs), numel(windowResults), 1);
phaseRows = struct([]);
shiftRows = struct([]);
sensorRows = struct([]);
for iw = 1:numel(windowResults)
    WR = windowResults(iw);
    bundle = WR.Bundle;
    audit = audit_one_window_local(WR, bundle, S06E.auditEOs, S06E.referenceEO);
    auditRows(iw) = audit.WindowRow;
    sensorRows = [sensorRows; audit.SensorRows(:)]; %#ok<AGROW>
    phaseRows = [phaseRows; audit_phase_modes_one_window_local( ...
        WR, bundle, S06E.auditEOs, S06E.referenceEO, S06E.phaseModes)]; %#ok<AGROW>
    shiftRows = [shiftRows; audit_template_shift_points_one_window_local( ...
        WR, bundle, S06E.auditEOs, S06E.referenceEO, S06E.phaseModes, S06E)]; %#ok<AGROW>
end

EquationAudit = struct2table(auditRows);
SensorAudit = struct2table(sensorRows);
PhaseModeAudit = struct2table(phaseRows);
ShiftPointAudit = struct2table(shiftRows);

outDir = fullfile(P.outputDir, '06_identification', 'equation_audit', sensor_tag_local(S06E.analysisSensors));
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
timeLabel = time_label_local(S06E.startTimeSec);
suffix = output_suffix_local(S06E.outputLabel);
outCsv = fullfile(outDir, sprintf('Step06E_EquationAudit_%s_20241106%s.csv', timeLabel, suffix));
outSensorCsv = fullfile(outDir, sprintf('Step06E_SensorResidualAudit_%s_20241106%s.csv', timeLabel, suffix));
outPhaseCsv = fullfile(outDir, sprintf('Step06E_PhaseModeAudit_%s_20241106%s.csv', timeLabel, suffix));
outShiftCsv = fullfile(outDir, sprintf('Step06E_TemplateShiftPointAudit_%s_20241106%s.csv', timeLabel, suffix));
writetable(EquationAudit, outCsv);
writetable(SensorAudit, outSensorCsv);
writetable(PhaseModeAudit, outPhaseCsv);
writetable(ShiftPointAudit, outShiftCsv);

fprintf('\nWindow-level equation audit:\n');
disp(EquationAudit(:, {'WindowID','ChosenEO','CandidateBestEO','CandidateGapToRefV','LinearBestEO','LinearGapToRefV','AliasMaxVsRef'}));
fprintf('\nPhase-mode audit, lower alias means EO candidates are easier to separate:\n');
disp(PhaseModeAudit(:, {'WindowID','PhaseMode','LinearBestEO','LinearGapToRefV','AliasMaxVsRef','ThetaSpanCycle'}));
fprintf('\nTemplate-shift point audit, each waveform pulse is compressed to one displacement point:\n');
disp(ShiftPointAudit(:, {'WindowID','PhaseMode','ShiftPointCount','ShiftBestEO','ShiftGapToRefMM','ShiftAliasMaxVsRef'}));
fprintf('Saved: %s\n', outCsv);
fprintf('Saved: %s\n', outSensorCsv);
fprintf('Saved: %s\n', outPhaseCsv);
fprintf('Saved: %s\n', outShiftCsv);

if S06E.plotEnable
    plot_equation_audit_local(EquationAudit, S06E, outDir, timeLabel, suffix);
end

function P = apply_step06e_local_options_local(P, S06E)
P.sensors.analysis = S06E.analysisSensors(:).';
P.region.startTimeSec = S06E.startTimeSec;
sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
suffix = output_suffix_local(S06E.outputLabel);
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
end

function audit = audit_one_window_local(WR, bundle, auditEOs, referenceEO)
T = get_candidate_table_local(WR);
candidateBestEO = NaN;
candidateGapToRef = NaN;
candidateRefRMSE = NaN;
candidateBestRMSE = NaN;
if ~isempty(T)
    [candidateBestRMSE, idx] = min(T.weighted_voltage_rmse, [], 'omitnan');
    candidateBestEO = T.EO(idx);
    refIdx = find(T.EO == referenceEO, 1, 'first');
    if ~isempty(refIdx)
        candidateRefRMSE = T.weighted_voltage_rmse(refIdx);
        candidateGapToRef = candidateRefRMSE - candidateBestRMSE;
    end
end

linearRows = repmat(struct('EO', NaN, 'RMSE', NaN, 'DxCMM', NaN, 'AmplitudeMM', NaN), numel(auditEOs), 1);
sensorRows = struct([]);
for ie = 1:numel(auditEOs)
    eo = auditEOs(ie);
    fit = solve_linearized_template_eo_local(bundle, eo);
    linearRows(ie).EO = eo;
    linearRows(ie).RMSE = fit.rmse;
    linearRows(ie).DxCMM = fit.dx_c;
    linearRows(ie).AmplitudeMM = hypot(fit.a_sin, fit.b_cos);
    sensorRows = [sensorRows; build_sensor_residual_rows_local(WR, bundle, eo, fit)]; %#ok<AGROW>
end
[linearBestRMSE, idxBest] = min([linearRows.RMSE]);
linearBestEO = linearRows(idxBest).EO;
refLinearIdx = find([linearRows.EO] == referenceEO, 1, 'first');
linearGapToRef = NaN;
linearRefRMSE = NaN;
if ~isempty(refLinearIdx)
    linearRefRMSE = linearRows(refLinearIdx).RMSE;
    linearGapToRef = linearRefRMSE - linearBestRMSE;
end

aliasMaxVsRef = max_alias_to_reference_local(bundle, auditEOs, referenceEO);

row = empty_audit_row_local(auditEOs);
row.WindowID = WR.window_id;
row.LapStart = WR.lap_range(1);
row.LapEnd = WR.lap_range(end);
row.ChosenEO = WR.Result.EO_id;
row.ChosenRMSE = WR.Result.weighted_voltage_rmse;
row.CandidateBestEO = candidateBestEO;
row.CandidateBestRMSE = candidateBestRMSE;
row.CandidateRefRMSE = candidateRefRMSE;
row.CandidateGapToRefV = candidateGapToRef;
row.LinearBestEO = linearBestEO;
row.LinearBestRMSE = linearBestRMSE;
row.LinearRefRMSE = linearRefRMSE;
row.LinearGapToRefV = linearGapToRef;
row.AliasMaxVsRef = aliasMaxVsRef;
row.PointCount = bundle.point_count;
for ie = 1:numel(linearRows)
    eo = linearRows(ie).EO;
    row.(sprintf('LinearRMSE_EO%d', eo)) = linearRows(ie).RMSE;
    row.(sprintf('LinearAmpMM_EO%d', eo)) = linearRows(ie).AmplitudeMM;
end

audit = struct('WindowRow', row, 'SensorRows', sensorRows);
end

function rows = audit_phase_modes_one_window_local(WR, bundle, auditEOs, referenceEO, phaseModes)
rows = repmat(empty_phase_row_local(), numel(phaseModes), 1);
for im = 1:numel(phaseModes)
    modeName = char(phaseModes{im});
    theta = build_theta_for_phase_mode_local(bundle, modeName);
    linearRows = repmat(struct('EO', NaN, 'RMSE', NaN), numel(auditEOs), 1);
    for ie = 1:numel(auditEOs)
        eo = auditEOs(ie);
        fit = solve_linearized_template_eo_local(bundle, eo, theta);
        linearRows(ie).EO = eo;
        linearRows(ie).RMSE = fit.rmse;
    end
    [linearBestRMSE, idxBest] = min([linearRows.RMSE]);
    linearBestEO = linearRows(idxBest).EO;
    refLinearIdx = find([linearRows.EO] == referenceEO, 1, 'first');
    linearRefRMSE = NaN;
    linearGapToRef = NaN;
    if ~isempty(refLinearIdx)
        linearRefRMSE = linearRows(refLinearIdx).RMSE;
        linearGapToRef = linearRefRMSE - linearBestRMSE;
    end

    rows(im).WindowID = WR.window_id;
    rows(im).LapStart = WR.lap_range(1);
    rows(im).LapEnd = WR.lap_range(end);
    rows(im).PhaseMode = string(modeName);
    rows(im).LinearBestEO = linearBestEO;
    rows(im).LinearBestRMSE = linearBestRMSE;
    rows(im).LinearRefRMSE = linearRefRMSE;
    rows(im).LinearGapToRefV = linearGapToRef;
    rows(im).AliasMaxVsRef = max_alias_to_reference_local(bundle, auditEOs, referenceEO, theta);
    rows(im).ThetaSpanCycle = range(theta(isfinite(theta))) / (2 * pi);
end
end

function theta = build_theta_for_phase_mode_local(bundle, modeName)
switch lower(strtrim(modeName))
    case 'stored_theta'
        theta = bundle.Theta;
    case 'window_time_theta'
        theta = 2 * pi * bundle.rot_freq_mean_hz .* (bundle.T - min(bundle.T, [], 'omitnan'));
    case 'absolute_time_theta'
        theta = 2 * pi * bundle.rot_freq_mean_hz .* bundle.T;
    otherwise
        error('Unknown phase mode: %s', modeName);
end
theta = theta(:);
end

function rows = audit_template_shift_points_one_window_local(WR, bundle, auditEOs, referenceEO, phaseModes, S06E)
shiftPoints = estimate_template_shift_points_local(bundle, S06E);
rows = repmat(empty_shift_point_row_local(), numel(phaseModes), 1);
for im = 1:numel(phaseModes)
    modeName = char(phaseModes{im});
    thetaPoint = build_shift_point_theta_local(shiftPoints, bundle, modeName);
    fitRows = repmat(struct('EO', NaN, 'RMSE', NaN), numel(auditEOs), 1);
    for ie = 1:numel(auditEOs)
        eo = auditEOs(ie);
        fitRows(ie).EO = eo;
        fitRows(ie).RMSE = solve_shift_point_eo_rmse_local(shiftPoints, thetaPoint, eo);
    end
    [bestRmse, idxBest] = min([fitRows.RMSE]);
    bestEO = fitRows(idxBest).EO;
    refIdx = find([fitRows.EO] == referenceEO, 1, 'first');
    refRmse = NaN;
    gapToRef = NaN;
    if ~isempty(refIdx)
        refRmse = fitRows(refIdx).RMSE;
        gapToRef = refRmse - bestRmse;
    end

    rows(im).WindowID = WR.window_id;
    rows(im).LapStart = WR.lap_range(1);
    rows(im).LapEnd = WR.lap_range(end);
    rows(im).PhaseMode = string(modeName);
    rows(im).ShiftPointCount = numel(shiftPoints.d_mm);
    rows(im).ShiftBestEO = bestEO;
    rows(im).ShiftBestRMSEMM = bestRmse;
    rows(im).ShiftRefRMSEMM = refRmse;
    rows(im).ShiftGapToRefMM = gapToRef;
    rows(im).ShiftAliasMaxVsRef = max_shift_point_alias_to_reference_local(shiftPoints, thetaPoint, auditEOs, referenceEO);
end
end

function shiftPoints = estimate_template_shift_points_local(bundle, S06E)
shiftPoints = struct('sensor_id', [], 't_sec', [], 'theta_stored', [], 'd_mm', [], 'point_count', []);
sensorIds = bundle.sensor_ids(:).';
for is = 1:numel(sensorIds)
    keepSensor = bundle.sensor_index == is & isfinite(bundle.T) & isfinite(bundle.V) & ...
        isfinite(bundle.F0) & isfinite(bundle.Fx) & isfinite(bundle.W) & isfinite(bundle.Theta);
    idxAll = find(keepSensor);
    if isempty(idxAll)
        continue;
    end
    [~, order] = sort(bundle.T(idxAll));
    idxAll = idxAll(order);
    dt = [inf; diff(bundle.T(idxAll))];
    segmentStart = find(dt > S06E.segmentGapSec);
    segmentEnd = [segmentStart(2:end) - 1; numel(idxAll)];
    for k = 1:numel(segmentStart)
        idx = idxAll(segmentStart(k):segmentEnd(k));
        if numel(idx) < S06E.minSegmentPoints
            continue;
        end
        w = max(bundle.W(idx), 0);
        fx = bundle.Fx(idx);
        y = bundle.V(idx) - bundle.F0(idx);
        valid = isfinite(w) & isfinite(fx) & isfinite(y) & w > 0;
        if nnz(valid) < S06E.minSegmentPoints
            continue;
        end
        fx = fx(valid);
        y = y(valid);
        w = w(valid);
        denom = sum(w .* fx.^2);
        if denom <= eps
            continue;
        end
        d = -sum(w .* fx .* y) / denom;
        idxValid = idx(valid);
        shiftPoints.sensor_id(end + 1, 1) = sensorIds(is); %#ok<AGROW>
        shiftPoints.t_sec(end + 1, 1) = median(bundle.T(idxValid), 'omitnan'); %#ok<AGROW>
        shiftPoints.theta_stored(end + 1, 1) = median(bundle.Theta(idxValid), 'omitnan'); %#ok<AGROW>
        shiftPoints.d_mm(end + 1, 1) = d; %#ok<AGROW>
        shiftPoints.point_count(end + 1, 1) = nnz(valid); %#ok<AGROW>
    end
end
end

function thetaPoint = build_shift_point_theta_local(shiftPoints, bundle, modeName)
switch lower(strtrim(modeName))
    case 'stored_theta'
        thetaPoint = shiftPoints.theta_stored;
    case 'window_time_theta'
        thetaPoint = 2 * pi * bundle.rot_freq_mean_hz .* ...
            (shiftPoints.t_sec - min(bundle.T, [], 'omitnan'));
    case 'absolute_time_theta'
        thetaPoint = 2 * pi * bundle.rot_freq_mean_hz .* shiftPoints.t_sec;
    otherwise
        error('Unknown phase mode: %s', modeName);
end
thetaPoint = thetaPoint(:);
end

function rmse = solve_shift_point_eo_rmse_local(shiftPoints, thetaPoint, eo)
valid = isfinite(shiftPoints.d_mm) & isfinite(thetaPoint);
y = shiftPoints.d_mm(valid);
theta = thetaPoint(valid);
if numel(y) < 4
    rmse = NaN;
    return;
end
A = [ones(size(theta)), sin(eo .* theta), cos(eo .* theta)];
p = A \ y;
res = y - A * p;
rmse = sqrt(mean(res.^2, 'omitnan'));
end

function aliasMax = max_shift_point_alias_to_reference_local(shiftPoints, thetaPoint, auditEOs, referenceEO)
refBasis = shift_point_basis_local(shiftPoints, thetaPoint, referenceEO);
aliasMax = NaN;
for eo = auditEOs(:).'
    if eo == referenceEO
        continue;
    end
    basis = shift_point_basis_local(shiftPoints, thetaPoint, eo);
    if isempty(refBasis) || isempty(basis)
        continue;
    end
    s = svd(refBasis' * basis);
    aliasMax = max([aliasMax; s(1)], [], 'omitnan');
end
end

function Q = shift_point_basis_local(shiftPoints, thetaPoint, eo)
valid = isfinite(shiftPoints.d_mm) & isfinite(thetaPoint);
theta = thetaPoint(valid);
B = [sin(eo .* theta), cos(eo .* theta)];
if size(B, 1) < 4 || rank(B) < 2
    Q = [];
    return;
end
[Q, ~] = qr(B, 0);
end

function T = get_candidate_table_local(WR)
T = table();
if isfield(WR, 'Result') && isfield(WR.Result, 'CandidateTable') && istable(WR.Result.CandidateTable)
    T = WR.Result.CandidateTable;
end
end

function fit = solve_linearized_template_eo_local(bundle, eo, thetaOverride)
if nargin < 3 || isempty(thetaOverride)
    thetaOverride = bundle.Theta;
end
finite = isfinite(bundle.V) & isfinite(bundle.F0) & isfinite(bundle.Fx) & ...
    isfinite(bundle.W) & isfinite(thetaOverride);
y = bundle.V(finite) - bundle.F0(finite);
fx = bundle.Fx(finite);
theta = thetaOverride(finite);
w = max(bundle.W(finite), 0);
A = [-fx, -fx .* sin(eo .* theta), -fx .* cos(eo .* theta)];
valid = all(isfinite(A), 2) & isfinite(y) & w > 0;
A = A(valid, :);
y = y(valid);
w = w(valid);
if size(A, 1) < size(A, 2) + 2
    fit = struct('rmse', NaN, 'dx_c', NaN, 'a_sin', NaN, 'b_cos', NaN, ...
        'prediction', nan(size(bundle.V)), 'valid_mask', false(size(bundle.V)));
    return;
end
Aw = A .* sqrt(w);
yw = y .* sqrt(w);
p = Aw \ yw;
res = y - A * p;
rmse = sqrt(sum(w .* res.^2) / max(sum(w), eps));
prediction = nan(size(bundle.V));
validIndex = find(finite);
validIndex = validIndex(valid);
prediction(validIndex) = bundle.F0(validIndex) + A * p;
validMask = false(size(bundle.V));
validMask(validIndex) = true;
fit = struct('rmse', rmse, 'dx_c', p(1), 'a_sin', p(2), 'b_cos', p(3), ...
    'prediction', prediction, 'valid_mask', validMask);
end

function rows = build_sensor_residual_rows_local(WR, bundle, eo, fit)
sensorIds = bundle.sensor_ids(:).';
rows = repmat(struct( ...
    'WindowID', WR.window_id, ...
    'SensorID', NaN, ...
    'EO', eo, ...
    'PointCount', NaN, ...
    'SensorRMSE', NaN), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    keep = fit.valid_mask & bundle.sensor_index == is & isfinite(fit.prediction) & isfinite(bundle.V);
    w = max(bundle.W(keep), 0);
    res = bundle.V(keep) - fit.prediction(keep);
    rows(is).SensorID = sid;
    rows(is).PointCount = nnz(keep);
    rows(is).SensorRMSE = sqrt(sum(w .* res.^2) / max(sum(w), eps));
end
end

function aliasMax = max_alias_to_reference_local(bundle, auditEOs, referenceEO, thetaOverride)
if nargin < 4 || isempty(thetaOverride)
    thetaOverride = bundle.Theta;
end
refBasis = weighted_basis_local(bundle, referenceEO, thetaOverride);
aliasMax = NaN;
for eo = auditEOs(:).'
    if eo == referenceEO
        continue;
    end
    basis = weighted_basis_local(bundle, eo, thetaOverride);
    if isempty(refBasis) || isempty(basis)
        continue;
    end
    s = svd(refBasis' * basis);
    aliasMax = max([aliasMax; s(1)], [], 'omitnan');
end
end

function Q = weighted_basis_local(bundle, eo, thetaOverride)
if nargin < 3 || isempty(thetaOverride)
    thetaOverride = bundle.Theta;
end
finite = isfinite(bundle.Fx) & isfinite(bundle.W) & isfinite(thetaOverride);
fx = bundle.Fx(finite);
theta = thetaOverride(finite);
w = max(bundle.W(finite), 0);
B = [-fx .* sin(eo .* theta), -fx .* cos(eo .* theta)];
valid = all(isfinite(B), 2) & w > 0;
B = B(valid, :) .* sqrt(w(valid));
if size(B, 1) < 4 || rank(B) < 2
    Q = [];
    return;
end
[Q, ~] = qr(B, 0);
end

function row = empty_phase_row_local()
row = struct( ...
    'WindowID', NaN, ...
    'LapStart', NaN, ...
    'LapEnd', NaN, ...
    'PhaseMode', "", ...
    'LinearBestEO', NaN, ...
    'LinearBestRMSE', NaN, ...
    'LinearRefRMSE', NaN, ...
    'LinearGapToRefV', NaN, ...
    'AliasMaxVsRef', NaN, ...
    'ThetaSpanCycle', NaN);
end

function row = empty_shift_point_row_local()
row = struct( ...
    'WindowID', NaN, ...
    'LapStart', NaN, ...
    'LapEnd', NaN, ...
    'PhaseMode', "", ...
    'ShiftPointCount', NaN, ...
    'ShiftBestEO', NaN, ...
    'ShiftBestRMSEMM', NaN, ...
    'ShiftRefRMSEMM', NaN, ...
    'ShiftGapToRefMM', NaN, ...
    'ShiftAliasMaxVsRef', NaN);
end

function row = empty_audit_row_local(auditEOs)
row = struct( ...
    'WindowID', NaN, ...
    'LapStart', NaN, ...
    'LapEnd', NaN, ...
    'ChosenEO', NaN, ...
    'ChosenRMSE', NaN, ...
    'CandidateBestEO', NaN, ...
    'CandidateBestRMSE', NaN, ...
    'CandidateRefRMSE', NaN, ...
    'CandidateGapToRefV', NaN, ...
    'LinearBestEO', NaN, ...
    'LinearBestRMSE', NaN, ...
    'LinearRefRMSE', NaN, ...
    'LinearGapToRefV', NaN, ...
    'AliasMaxVsRef', NaN, ...
    'PointCount', NaN);
for eo = auditEOs(:).'
    row.(sprintf('LinearRMSE_EO%d', eo)) = NaN;
    row.(sprintf('LinearAmpMM_EO%d', eo)) = NaN;
end
end

function plot_equation_audit_local(EquationAudit, S06E, outDir, timeLabel, suffix)
fig = figure('Name', 'Step06E equation audit', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 13], 'NumberTitle', 'off');
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, EquationAudit.WindowID, EquationAudit.ChosenEO, 'o-', 'LineWidth', 1.1);
yline(ax1, S06E.referenceEO, '--', 'EO ref');
grid(ax1, 'on');
ylabel(ax1, 'Chosen EO');
title(ax1, sprintf('B%d %s equation audit', S06E.bladeId, timeLabel), 'FontWeight', 'normal');

ax2 = nexttile;
plot(ax2, EquationAudit.WindowID, EquationAudit.CandidateGapToRefV, 'o-', 'LineWidth', 1.1);
yline(ax2, 0, '--');
grid(ax2, 'on');
ylabel(ax2, '\Delta RMSE to ref (V)');

ax3 = nexttile;
plot(ax3, EquationAudit.WindowID, EquationAudit.AliasMaxVsRef, 'o-', 'LineWidth', 1.1);
ylim(ax3, [0, 1.05]);
grid(ax3, 'on');
xlabel(ax3, 'Window ID');
ylabel(ax3, 'max basis similarity');

if S06E.saveFigures
    exportgraphics(fig, fullfile(outDir, sprintf('Step06E_EquationAudit_%s_20241106%s.png', timeLabel, suffix)), ...
        'Resolution', 200);
end
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = time_label_local(tSec)
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

function suffix = output_suffix_local(labelText)
labelText = strtrim(char(labelText));
if isempty(labelText)
    suffix = '';
else
    suffix = ['_' labelText];
end
end

