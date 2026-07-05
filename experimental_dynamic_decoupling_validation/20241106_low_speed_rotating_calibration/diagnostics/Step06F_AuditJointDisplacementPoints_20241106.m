%% Step06F_AuditJointDisplacementPoints_20241106
% Diagnose whether the current Step06 waveform equation has an EO
% identifiability problem.
%
% This script reads a saved Step06 result, compresses each waveform pulse to
% one low-speed-template shift point, then compares EO candidates using two
% simple displacement-point models:
%
%   shared_time_phase:
%       d = blade/sensor offset + a*sin(EO*theta_time) + b*cos(EO*theta_time)
%
%   blade_free_phase:
%       d = blade/sensor offset + per-blade sin/cos terms
%
% If shared_time_phase points to EO12 but blade_free_phase is ambiguous, the
% current single-blade waveform fit is throwing away global time-phase
% information.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

%% Parameters to tune
S06F = struct();
S06F.resultFile = fullfile(P.outputDir, '06_identification', 'S2357', ...
    'IdentificationResult_T075s_20241106.mat');
S06F.auditEOs = 1:20;
S06F.referenceEO = 12;
S06F.segmentGapSec = 2e-4;
S06F.minSegmentPoints = 20;
S06F.fitModels = {'shared_time_phase', 'blade_free_phase'};
S06F.phaseModes = {'stored_theta', 'window_time_theta', 'absolute_time_theta'};
S06F.plotEnable = true;

require_file_local(S06F.resultFile, 'Step06 identification result');
loaded = load(S06F.resultFile, 'IdentificationResult');
IdentificationResult = loaded.IdentificationResult;

fprintf('\n=== Step06F: joint displacement-point EO audit ===\n');
fprintf('Source: %s\n', S06F.resultFile);
fprintf('EO candidates: %s, reference EO%d\n', mat2str(S06F.auditEOs), S06F.referenceEO);
drawnow;

pointRows = struct([]);
for iw = 1:numel(IdentificationResult.WindowResult)
    WR = IdentificationResult.WindowResult(iw);
    if isempty(WR.Bundle) || ~isfield(WR.Bundle, 'T')
        continue;
    end
    pts = estimate_template_shift_points_local(WR, S06F);
    pointRows = [pointRows; pts(:)]; %#ok<AGROW>
end

if isempty(pointRows)
    error('No template-shift points were produced from %s.', S06F.resultFile);
end

PointTable = struct2table(pointRows);
fprintf('Template-shift points: %d | blades=%s | sensors=%s | windows=%s\n', ...
    height(PointTable), mat2str(unique(PointTable.BladeID).'), ...
    mat2str(unique(PointTable.SensorID).'), mat2str(unique(PointTable.WindowID).'));

auditRows = struct([]);
windowIds = unique(PointTable.WindowID).';
for iw = windowIds
    keepWindow = PointTable.WindowID == iw;
    for ip = 1:numel(S06F.phaseModes)
        phaseMode = S06F.phaseModes{ip};
        Tw = attach_phase_to_points_local(PointTable(keepWindow, :), phaseMode);
        for im = 1:numel(S06F.fitModels)
            modelName = S06F.fitModels{im};
            rows = audit_one_window_model_local(Tw, S06F.auditEOs, ...
                S06F.referenceEO, modelName, phaseMode);
            auditRows = [auditRows; rows(:)]; %#ok<AGROW>
        end
    end
end

JointEOAudit = struct2table(auditRows);
outDir = fullfile(P.outputDir, '06_identification', 'equation_audit', 'joint_displacement_points');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outPointCsv = fullfile(outDir, 'Step06F_TemplateShiftPoints_20241106.csv');
outAuditCsv = fullfile(outDir, 'Step06F_JointEOAudit_20241106.csv');
writetable(PointTable, outPointCsv);
writetable(JointEOAudit, outAuditCsv);

fprintf('\nJoint EO audit summary:\n');
disp(JointEOAudit(:, {'WindowID','PhaseMode','ModelName','PointCount','BestEO','GapToRefMM','AliasMaxVsRef','BestRMSEMM','RefRMSEMM'}));
fprintf('Saved points: %s\n', outPointCsv);
fprintf('Saved audit:  %s\n', outAuditCsv);

if S06F.plotEnable
    plot_joint_eo_audit_local(JointEOAudit, S06F);
end

function rows = audit_one_window_model_local(T, auditEOs, referenceEO, modelName, phaseMode)
rows = repmat(empty_audit_row_local(), 1, 1);
rmse = nan(numel(auditEOs), 1);
aliasToRef = nan(numel(auditEOs), 1);
refBasis = model_basis_local(T, referenceEO, modelName);
for ie = 1:numel(auditEOs)
    eo = auditEOs(ie);
    rmse(ie) = fit_model_rmse_local(T, eo, modelName);
    if eo ~= referenceEO
        basis = model_basis_local(T, eo, modelName);
        if ~isempty(refBasis) && ~isempty(basis)
            s = svd(refBasis' * basis);
            aliasToRef(ie) = s(1);
        end
    end
end
[bestRmse, idxBest] = min(rmse, [], 'omitnan');
bestEO = auditEOs(idxBest);
refIdx = find(auditEOs == referenceEO, 1, 'first');
refRmse = NaN;
gapToRef = NaN;
if ~isempty(refIdx)
    refRmse = rmse(refIdx);
    gapToRef = refRmse - bestRmse;
end
rows.WindowID = T.WindowID(1);
rows.PhaseMode = string(phaseMode);
rows.ModelName = string(modelName);
rows.PointCount = height(T);
rows.BestEO = bestEO;
rows.BestRMSEMM = bestRmse;
rows.RefRMSEMM = refRmse;
rows.GapToRefMM = gapToRef;
rows.AliasMaxVsRef = max(aliasToRef, [], 'omitnan');
end

function rmse = fit_model_rmse_local(T, eo, modelName)
[A, y] = model_design_local(T, eo, modelName);
valid = all(isfinite(A), 2) & isfinite(y);
A = A(valid, :);
y = y(valid);
if size(A, 1) <= size(A, 2)
    rmse = NaN;
    return;
end
p = A \ y;
res = y - A * p;
rmse = sqrt(mean(res.^2, 'omitnan'));
end

function Q = model_basis_local(T, eo, modelName)
[A, ~, vibCols] = model_design_local(T, eo, modelName);
B = A(:, vibCols);
valid = all(isfinite(B), 2);
B = B(valid, :);
if size(B, 1) <= size(B, 2) || rank(B) < min(size(B))
    Q = [];
    return;
end
[Q, ~] = qr(B, 0);
end

function [A, y, vibCols] = model_design_local(T, eo, modelName)
y = T.ShiftMM;
theta = T.PhaseRad;
groups = findgroups(T.BladeID, T.SensorID);
nGroup = max(groups);
offsetA = zeros(height(T), nGroup);
for i = 1:height(T)
    offsetA(i, groups(i)) = 1;
end

switch lower(strtrim(modelName))
    case 'shared_time_phase'
        vibA = [sin(eo .* theta), cos(eo .* theta)];
    case 'blade_free_phase'
        bladeGroups = findgroups(T.BladeID);
        nBlade = max(bladeGroups);
        vibA = zeros(height(T), 2 * nBlade);
        for i = 1:height(T)
            b = bladeGroups(i);
            vibA(i, 2*b - 1) = sin(eo .* theta(i));
            vibA(i, 2*b) = cos(eo .* theta(i));
        end
    otherwise
        error('Unknown modelName: %s', modelName);
end

A = [offsetA, vibA];
vibCols = (size(offsetA, 2) + 1):size(A, 2);
end

function pts = estimate_template_shift_points_local(WR, S06F)
bundle = WR.Bundle;
pts = repmat(empty_point_row_local(), 0, 1);
sensorIds = bundle.sensor_ids(:).';
for is = 1:numel(sensorIds)
    keepSensor = bundle.sensor_index == is & isfinite(bundle.T) & isfinite(bundle.V) & ...
        isfinite(bundle.F0) & isfinite(bundle.Fx) & isfinite(bundle.W);
    idxAll = find(keepSensor);
    if isempty(idxAll)
        continue;
    end
    [~, order] = sort(bundle.T(idxAll));
    idxAll = idxAll(order);
    dt = [inf; diff(bundle.T(idxAll))];
    segmentStart = find(dt > S06F.segmentGapSec);
    segmentEnd = [segmentStart(2:end) - 1; numel(idxAll)];
    for k = 1:numel(segmentStart)
        idx = idxAll(segmentStart(k):segmentEnd(k));
        if numel(idx) < S06F.minSegmentPoints
            continue;
        end
        w = max(bundle.W(idx), 0);
        fx = bundle.Fx(idx);
        y = bundle.V(idx) - bundle.F0(idx);
        valid = isfinite(w) & isfinite(fx) & isfinite(y) & w > 0;
        if nnz(valid) < S06F.minSegmentPoints
            continue;
        end
        denom = sum(w(valid) .* fx(valid).^2);
        if denom <= eps
            continue;
        end
        d = -sum(w(valid) .* fx(valid) .* y(valid)) / denom;
        idxValid = idx(valid);
        row = empty_point_row_local();
        row.BladeID = WR.blade_id;
        row.WindowID = WR.window_id;
        row.LapStart = WR.lap_range(1);
        row.LapEnd = WR.lap_range(end);
        row.SensorID = sensorIds(is);
        row.TimeSec = median(bundle.T(idxValid), 'omitnan');
        row.TimeRelSec = median(bundle.T(idxValid), 'omitnan') - min(bundle.T, [], 'omitnan');
        row.ThetaStoredRad = median(bundle.Theta(idxValid), 'omitnan');
        row.RotFreqHz = bundle.rot_freq_mean_hz;
        row.ShiftMM = d;
        row.PointCount = nnz(valid);
        pts(end + 1, 1) = row; %#ok<AGROW>
    end
end
end

function row = empty_point_row_local()
row = struct('BladeID', NaN, 'WindowID', NaN, 'LapStart', NaN, 'LapEnd', NaN, ...
    'SensorID', NaN, 'TimeSec', NaN, 'TimeRelSec', NaN, 'ThetaStoredRad', NaN, ...
    'RotFreqHz', NaN, 'ShiftMM', NaN, 'PointCount', NaN);
end

function row = empty_audit_row_local()
row = struct('WindowID', NaN, 'PhaseMode', "", 'ModelName', "", 'PointCount', NaN, ...
    'BestEO', NaN, 'BestRMSEMM', NaN, 'RefRMSEMM', NaN, ...
    'GapToRefMM', NaN, 'AliasMaxVsRef', NaN);
end

function T = attach_phase_to_points_local(T, phaseMode)
switch lower(strtrim(phaseMode))
    case 'stored_theta'
        T.PhaseRad = T.ThetaStoredRad;
    case 'window_time_theta'
        T.PhaseRad = 2 * pi * T.RotFreqHz .* T.TimeRelSec;
    case 'absolute_time_theta'
        T.PhaseRad = 2 * pi * T.RotFreqHz .* T.TimeSec;
    otherwise
        error('Unknown phaseMode: %s', phaseMode);
end
end

function plot_joint_eo_audit_local(JointEOAudit, S06F)
fig = figure('Name', 'Step06F joint displacement-point EO audit', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 24, 13], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; grid on; box on;
models = unique(JointEOAudit.ModelName, 'stable');
for i = 1:numel(models)
    keep = JointEOAudit.ModelName == models(i) & JointEOAudit.PhaseMode == "stored_theta";
    plot(JointEOAudit.WindowID(keep), JointEOAudit.BestEO(keep), 'o-', ...
        'LineWidth', 1.1, 'DisplayName', char(models(i)));
end
yline(S06F.referenceEO, '--', sprintf('EO%d ref', S06F.referenceEO));
ylabel('Best EO');
legend('Location', 'best');
title('Joint displacement-point EO audit', 'FontWeight', 'normal');

nexttile;
hold on; grid on; box on;
for i = 1:numel(models)
    keep = JointEOAudit.ModelName == models(i) & JointEOAudit.PhaseMode == "stored_theta";
    plot(JointEOAudit.WindowID(keep), JointEOAudit.AliasMaxVsRef(keep), 'o-', ...
        'LineWidth', 1.1, 'DisplayName', char(models(i)));
end
ylim([0, 1.05]);
xlabel('Window ID');
ylabel('max alias vs ref');
legend('Location', 'best');
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

