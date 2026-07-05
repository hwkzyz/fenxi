%% Step07B_AuditBladeIdentityAndEquation_20241106
% Audit two questions without changing the Step06 solver:
%   1) is the selected high-speed window really better matched to this blade
%      than to the other low-speed blade templates,
%   2) does the full waveform equation actually improve the fit compared with
%      simpler template-mapping models.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local audit parameter block. Edit here directly.
A08B = struct();
A08B.analysisSensors = [5 7];
A08B.startTimeSec = 75.0;
A08B.outputLabel = 'B6_only';
A08B.referenceBladeId = 6;
A08B.windowIds = [];            % [] -> audit all windows in Step06 result
A08B.compareBladeIds = 1:6;     % audit all low-speed templates against the same high-speed window
A08B.viewEnable = true;

P = apply_local_options_local(P, A08B);

IdentificationResult = load_identification_result_local(P);
[LowSpeedTemplateLibrary, sourceInfo] = load_compatible_template_library_local(P);

windowResults = IdentificationResult.WindowResult;
windowResults = windowResults([windowResults.blade_id] == A08B.referenceBladeId);
if ~isempty(A08B.windowIds)
    keep = ismember([windowResults.window_id], A08B.windowIds);
    windowResults = windowResults(keep);
end
if isempty(windowResults)
    error('No Step06 window result found for B%d.', A08B.referenceBladeId);
end

BladeCompetitionTable = build_blade_competition_table_local( ...
    windowResults, LowSpeedTemplateLibrary, A08B.compareBladeIds);
EquationAuditTable = build_equation_audit_table_local(windowResults, A08B.referenceBladeId);

fprintf('\n=== Step07B: blade identity and equation audit ===\n');
fprintf('Reference blade in Step06 result: B%d\n', A08B.referenceBladeId);
fprintf('Sensors: %s\n', mat2str(P.sensors.analysis));
fprintf('Template source: %s\n', sourceInfo.template_file);
fprintf('Identification source: %s\n', P.files.identificationResult);
disp(EquationAuditTable);
disp(BladeCompetitionTable);

if A08B.viewEnable
    plot_blade_competition_heatmap_local(BladeCompetitionTable, A08B.referenceBladeId, A08B.compareBladeIds);
    plot_equation_audit_local(EquationAuditTable, A08B.referenceBladeId);
end

function P = apply_local_options_local(P, A08B)
P.sensors.analysis = A08B.analysisSensors(:).';
P.region.startTimeSec = A08B.startTimeSec;

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
suffix = output_suffix_local(A08B.outputLabel);
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');
end

function IdentificationResult = load_identification_result_local(P)
require_file_local(P.files.identificationResult, 'Step06 identification result');
loaded = load(P.files.identificationResult, 'IdentificationResult');
IdentificationResult = loaded.IdentificationResult;
end

function [LowSpeedTemplateLibrary, sourceInfo] = load_compatible_template_library_local(P)
sourceInfo = struct('template_file', "");
requested = unique(P.sensors.analysis(:).', 'stable');
if exist(P.files.lowSpeedTemplateLibrary, 'file') == 2
    loaded = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
    LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
    sourceInfo.template_file = string(P.files.lowSpeedTemplateLibrary);
    return;
end

candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'LowSpeedTemplateLibrary');
    if ~isfield(loaded, 'LowSpeedTemplateLibrary') || ~isfield(loaded.LowSpeedTemplateLibrary, 'sensor_ids')
        continue;
    end
    available = unique(loaded.LowSpeedTemplateLibrary.sensor_ids(:).', 'stable');
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing compatible low-speed template library for sensors %s.', mat2str(requested));
end

fileNow = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(fileNow, 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
sourceInfo.template_file = string(fileNow);
end

function T = build_blade_competition_table_local(windowResults, LowSpeedTemplateLibrary, compareBladeIds)
rows = [];
for iw = 1:numel(windowResults)
    WR = windowResults(iw);
    bundle = WR.Bundle;
    result = WR.Result;
    for bladeId = compareBladeIds(:).'
        [rmseWeighted, rmsePlain, sensorWise] = evaluate_blade_template_rmse_local(bundle, result, LowSpeedTemplateLibrary, bladeId);
        row = struct();
        row.WindowID = WR.window_id;
        row.ReferenceBladeID = WR.blade_id;
        row.CandidateBladeID = bladeId;
        row.ReferenceEO = result.EO_id;
        row.ReferenceFrequencyHz = result.fn_id;
        row.ReferenceAmplitudeMM = result.A_id;
        row.CandidateWeightedRMSE = rmseWeighted;
        row.CandidatePlainRMSE = rmsePlain;
        for is = 1:numel(bundle.sensor_ids)
            sid = bundle.sensor_ids(is);
            row.(sprintf('CH%d_RMSE', sid)) = sensorWise(is);
        end
        rows = [rows; row]; %#ok<AGROW>
    end
end
T = struct2table(rows);
end

function T = build_equation_audit_table_local(windowResults, referenceBladeId)
rows = [];
for iw = 1:numel(windowResults)
    WR = windowResults(iw);
    bundle = WR.Bundle;
    result = WR.Result;
    [rmseTemplateOnly, rmseStaticShift, rmseFull] = evaluate_equation_levels_local(bundle, result);
    rows = [rows; struct( ... %#ok<AGROW>
        'BladeID', referenceBladeId, ...
        'WindowID', WR.window_id, ...
        'EO', result.EO_id, ...
        'FrequencyHz', result.fn_id, ...
        'AmplitudeMM', result.A_id, ...
        'DxCMM', result.dx_c_id, ...
        'TemplateOnlyRMSE', rmseTemplateOnly, ...
        'StaticShiftRMSE', rmseStaticShift, ...
        'FullEquationRMSE', rmseFull, ...
        'StaticImprovementV', rmseTemplateOnly - rmseStaticShift, ...
        'DynamicImprovementV', rmseStaticShift - rmseFull, ...
        'ChosenResultRMSE', result.weighted_voltage_rmse)];
    end
T = struct2table(rows);
end

function [rmseWeighted, rmsePlain, sensorWise] = evaluate_blade_template_rmse_local(bundle, result, LowSpeedTemplateLibrary, bladeId)
sensorIds = bundle.sensor_ids(:).';
sensorWise = nan(1, numel(sensorIds));
vPred = nan(size(bundle.V));
weights = bundle.W(:);
uEst = result.u_est(:);
sensorEta = result.sensor_eta_id(:).';

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    rows = find(bundle.S == sid);
    Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid);
    xQuery = bundle.X(rows) - result.dx_c_id - sensorEta(min(is, numel(sensorEta))) - uEst(rows);
    vPred(rows) = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xQuery, 'pchip', 'extrap');
    sensorWise(is) = rmse_local(bundle.V(rows), vPred(rows), ones(size(rows)));
end

rmseWeighted = rmse_local(bundle.V(:), vPred(:), weights);
rmsePlain = rmse_local(bundle.V(:), vPred(:), ones(size(vPred)));
end

function [rmseTemplateOnly, rmseStaticShift, rmseFull] = evaluate_equation_levels_local(bundle, result)
vTemplateOnly = nan(size(bundle.V));
vStaticShift = nan(size(bundle.V));
vFull = nan(size(bundle.V));
weights = bundle.W(:);
sensorEta = result.sensor_eta_id(:).';

for is = 1:numel(bundle.sensor_ids)
    rows = find(bundle.S == bundle.sensor_ids(is));
    interpV = bundle.interp_v{is};
    xObs = bundle.X(rows);
    eta = sensorEta(min(is, numel(sensorEta)));
    vTemplateOnly(rows) = interpV(xObs);
    vStaticShift(rows) = interpV(xObs - result.dx_c_id - eta);
    vFull(rows) = interpV(xObs - result.dx_c_id - eta - result.u_est(rows));
end

rmseTemplateOnly = rmse_local(bundle.V(:), vTemplateOnly(:), weights);
rmseStaticShift = rmse_local(bundle.V(:), vStaticShift(:), weights);
rmseFull = rmse_local(bundle.V(:), vFull(:), weights);
end

function plot_blade_competition_heatmap_local(BladeCompetitionTable, referenceBladeId, compareBladeIds)
windowIds = unique(BladeCompetitionTable.WindowID(:).', 'stable');
rmseMat = nan(numel(compareBladeIds), numel(windowIds));
for ib = 1:numel(compareBladeIds)
    for iw = 1:numel(windowIds)
        row = BladeCompetitionTable( ...
            BladeCompetitionTable.CandidateBladeID == compareBladeIds(ib) & ...
            BladeCompetitionTable.WindowID == windowIds(iw), :);
        if height(row) == 1
            rmseMat(ib, iw) = row.CandidateWeightedRMSE;
        end
    end
end

fig = figure('Name', sprintf('Step07B blade competition B%d', referenceBladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 10], 'NumberTitle', 'off');
imagesc(windowIds, compareBladeIds, rmseMat);
set(gca, 'YDir', 'normal');
xlabel('Window ID');
ylabel('Candidate low-speed blade template');
title(sprintf('B%d high-speed windows tested against B1-B6 low-speed templates', referenceBladeId), ...
    'FontWeight', 'normal');
cb = colorbar;
cb.Label.String = 'Weighted RMSE (V)';
box on;
style_axes_local();
end

function plot_equation_audit_local(EquationAuditTable, referenceBladeId)
fig = figure('Name', sprintf('Step07B equation audit B%d', referenceBladeId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 14], 'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(EquationAuditTable.WindowID, EquationAuditTable.TemplateOnlyRMSE, 'o-', 'LineWidth', 1.1, ...
    'DisplayName', 'template only'); hold on;
plot(EquationAuditTable.WindowID, EquationAuditTable.StaticShiftRMSE, 's-', 'LineWidth', 1.1, ...
    'DisplayName', 'static shift'); 
plot(EquationAuditTable.WindowID, EquationAuditTable.FullEquationRMSE, 'd-', 'LineWidth', 1.2, ...
    'DisplayName', 'full equation');
xlabel('Window ID');
ylabel('Weighted RMSE (V)');
title(sprintf('B%d waveform equation level comparison', referenceBladeId), 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on; style_axes_local();

nexttile;
bar(EquationAuditTable.WindowID, [EquationAuditTable.StaticImprovementV, EquationAuditTable.DynamicImprovementV], ...
    'grouped');
xlabel('Window ID');
ylabel('RMSE reduction (V)');
title('How much each equation term improves the fit', 'FontWeight', 'normal');
legend({'template -> static shift', 'static shift -> full equation'}, 'Location', 'best', 'Box', 'off');
grid on; box on; style_axes_local();
end

function Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid)
entry = LowSpeedTemplateLibrary.entry(bladeId);
if ~entry.exists
    error('Missing low-speed template bundle for B%d.', bladeId);
end
loadedBundle = load(char(entry.template_file), 'TemplateBundle');
bundle = loadedBundle.TemplateBundle;
row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
if height(row) ~= 1
    error('Low-speed template bundle for B%d does not contain CH%d.', bladeId, sid);
end
loadedSensor = load(char(row.template_file(1)), 'sensor_template');
Tpl = loadedSensor.sensor_template.Sensor;
end

function y = rmse_local(vObs, vPred, weights)
keep = isfinite(vObs) & isfinite(vPred) & isfinite(weights);
if ~any(keep)
    y = NaN;
    return;
end
vObs = vObs(keep);
vPred = vPred(keep);
weights = weights(keep);
weights = max(weights, eps);
y = sqrt(sum(weights .* (vObs - vPred).^2) / sum(weights));
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

function suffix = output_suffix_local(label)
suffix = '';
if isempty(label)
    return;
end
label = regexprep(char(label), '[^\w\d-]', '_');
label = regexprep(label, '_+', '_');
label = strtrim(label);
if ~isempty(label)
    suffix = ['_', label];
end
end

