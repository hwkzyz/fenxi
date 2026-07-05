%% Step06H_AuditCH3TemplateNumberingAngle_20241106
% CH3 focused health check: template, blade numbering, and mounting angle.
%
% This script does not rerun numbering or identification. It only reads the
% current artifacts and answers a narrow question:
%   Does CH3 look abnormal before it enters the multi-sensor identification?

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

%% Parameters to tune
A = struct();
A.auditSensor = 3;
A.compareSensors = [2 3 5 7];
A.bladeId = 4;
A.startTimeSec = 75.0;
A.sensorTag = 'S2357';
A.templateSuffix = 'GradientXRange030_OPRCenterStd';
A.targetEO = 12;
A.eoRange = 1:30;
A.maxHighSpeedLapsToPlot = 6;
A.showFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(routeDir, 'output', 'new_flow');
timeLabel = time_label_local(A.startTimeSec);

sensorConfigFile = fullfile(routeDir, 'output', 'reference', '900_reference', ...
    'Sensor_Config_20241106.mat');
lowReferenceFile = fullfile(outputDir, '01_low_speed_reference', ...
    'LowSpeedReferenceFingerprint_20241106.mat');
lowNumberingFile = fullfile(outputDir, '02_low_speed_numbering', A.sensorTag, ...
    'LowSpeedNumbering_20241106.mat');
highNumberingFile = fullfile(outputDir, '04_high_speed_numbering', A.sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
waveformFile = fullfile(outputDir, '05_waveform_library', A.sensorTag, ...
    sprintf('WaveformLibrary_%s_20241106.mat', timeLabel));

require_file_local(sensorConfigFile, 'Sensor_Config');
require_file_local(lowReferenceFile, 'LowSpeedReference');
require_file_local(lowNumberingFile, 'LowSpeedNumbering');
require_file_local(highNumberingFile, 'HighSpeedNumbering');

loadedCfg = load(sensorConfigFile, 'Sensor_Config');
Sensor_Config = loadedCfg.Sensor_Config;
loadedRef = load(lowReferenceFile, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;
loadedLow = load(lowNumberingFile, 'LowSpeedNumbering');
LowSpeedNumbering = loadedLow.LowSpeedNumbering;
loadedHigh = load(highNumberingFile, 'HighSpeedNumbering');
HighSpeedNumbering = loadedHigh.HighSpeedNumbering;

%% Read artifact tables
TemplateTable = build_template_table_local(outputDir, A);
LowNumberingTable = summarize_low_numbering_local(LowSpeedNumbering, A.compareSensors);
HighNumberingTable = summarize_high_numbering_local(HighSpeedNumbering, A.compareSensors);
AngleTable = build_angle_table_local(Sensor_Config, LowSpeedReference, A);
TemplateTable = attach_angle_to_template_table_local(TemplateTable, AngleTable);
AliasTable = build_alias_table_local(AngleTable, A);
EtaSummaryTable = summarize_existing_eta_results_local(outputDir, A);
WaveformResidualTable = table();
EtaInitialTable = table();
if exist(waveformFile, 'file') == 2
    loadedWaveForResidual = load(waveformFile, 'WaveformLibrary');
    WaveformResidualTable = summarize_waveform_residual_local(TemplateTable, loadedWaveForResidual.WaveformLibrary, A);
    EtaInitialTable = summarize_eta_initial_guess_local(TemplateTable, loadedWaveForResidual.WaveformLibrary, A);
end

fprintf('\n=== Step06H: CH%d template / numbering / mounting-angle health check ===\n', A.auditSensor);
fprintf('Blade: B%d, sensor set artifact: %s, start time: %.3f s\n', ...
    A.bladeId, A.sensorTag, A.startTimeSec);
fprintf('This script reads existing artifacts only; it does not rerun Step04/Step06.\n\n');

fprintf('1) Low-speed template metrics\n');
disp(TemplateTable);
fprintf('2) Low-speed blade numbering summary\n');
disp(LowNumberingTable);
fprintf('3) High-speed blade numbering summary\n');
disp(HighNumberingTable);
fprintf('4) OPR-center mounting-angle table for B%d\n', A.bladeId);
disp(AngleTable);
fprintf('5) EO%d mounting-angle alias check\n', A.targetEO);
disp(AliasTable);
fprintf('6) Existing identification eta summary\n');
disp(EtaSummaryTable);
if ~isempty(WaveformResidualTable)
    fprintf('7) Low-speed template vs high-speed raw waveform residual, first %d laps\n', ...
        A.maxHighSpeedLapsToPlot);
    disp(WaveformResidualTable);
end
if ~isempty(EtaInitialTable)
    fprintf('8) Eta initial guess from raw waveform vs low-speed template inverse, first %d laps\n', ...
        A.maxHighSpeedLapsToPlot);
    disp(EtaInitialTable);
end

print_ch3_interpretation_local(TemplateTable, LowNumberingTable, HighNumberingTable, ...
    AngleTable, AliasTable, EtaSummaryTable, WaveformResidualTable, EtaInitialTable, A);

%% Figures
if A.showFigures
    plot_template_comparison_local(TemplateTable, A);
    if exist(waveformFile, 'file') == 2
        plot_low_high_overlay_local(TemplateTable, loadedWaveForResidual.WaveformLibrary, A);
    else
        fprintf('\nWaveform library not found, skip low/high overlay:\n  %s\n', waveformFile);
    end
    if ~isempty(EtaSummaryTable)
        plot_eta_summary_local(EtaSummaryTable, A);
    end
    plot_angle_signature_local(AngleTable, AliasTable, A);
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function label = time_label_local(tSec)
if abs(tSec - round(tSec)) < 1e-9
    label = sprintf('T%03ds', round(tSec));
else
    label = sprintf('T%07.3fs', tSec);
    label = strrep(label, '.', 'p');
end
end

function TemplateTable = build_template_table_local(outputDir, A)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'AngleDeg', NaN, ...
    'XcMM', NaN, ...
    'XLeftMM', NaN, ...
    'XRightMM', NaN, ...
    'WidthMM', NaN, ...
    'ThresholdV', NaN, ...
    'BaselineV', NaN, ...
    'PointCount', NaN, ...
    'WidePointCount', NaN, ...
    'LapCount', NaN, ...
    'VMin', NaN, ...
    'VMax', NaN, ...
    'VRange', NaN, ...
    'TemplateFile', ""), numel(A.compareSensors), 1);

for i = 1:numel(A.compareSensors)
    sid = A.compareSensors(i);
    templateFile = fullfile(outputDir, 'templates', A.sensorTag, sprintf( ...
        'Template_LowSpeedRotating_B%d_CH%d_%s_20241106.mat', ...
        A.bladeId, sid, A.templateSuffix));
    require_file_local(templateFile, sprintf('B%d CH%d low-speed template', A.bladeId, sid));
    loaded = load(templateFile, 'sensor_template');
    T = loaded.sensor_template.Sensor;

    rows(i).SensorID = sid;
    rows(i).XcMM = T.xc;
    rows(i).XLeftMM = T.x_domain(1);
    rows(i).XRightMM = T.x_domain(2);
    rows(i).WidthMM = diff(T.x_domain);
    rows(i).ThresholdV = T.threshold;
    rows(i).BaselineV = T.baseline;
    rows(i).PointCount = T.point_count;
    rows(i).WidePointCount = T.wide_point_count;
    rows(i).LapCount = T.lap_count;
    rows(i).VMin = min(T.v_grid, [], 'omitnan');
    rows(i).VMax = max(T.v_grid, [], 'omitnan');
    rows(i).VRange = rows(i).VMax - rows(i).VMin;
    rows(i).TemplateFile = string(templateFile);
end

TemplateTable = struct2table(rows);
end

function LowNumberingTable = summarize_low_numbering_local(LowSpeedNumbering, sensorIds)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ReferenceRevolutionID', NaN, ...
    'MeanBestScore', NaN, ...
    'MeanScoreMargin', NaN, ...
    'MinScoreMargin', NaN, ...
    'RevolutionCount', NaN), numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    keep = LowSpeedNumbering.SensorID == sid;
    rowsNow = LowSpeedNumbering(keep, :);
    if isempty(rowsNow)
        continue;
    end
    rows(i).SensorID = sid;
    rows(i).DominantShift = mode(rowsNow.DominantShift);
    rows(i).LocalSlotOfB1 = mode(rowsNow.LocalSlotOfB1);
    rows(i).ReferenceRevolutionID = mode(rowsNow.ReferenceRevolutionID);
    rows(i).MeanBestScore = mean(rowsNow.BestScore, 'omitnan');
    rows(i).MeanScoreMargin = mean(rowsNow.ScoreMargin, 'omitnan');
    rows(i).MinScoreMargin = min(rowsNow.ScoreMargin, [], 'omitnan');
    rows(i).RevolutionCount = height(rowsNow);
end

LowNumberingTable = struct2table(rows);
end

function HighNumberingTable = summarize_high_numbering_local(HighSpeedNumbering, sensorIds)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ValidRevolutionCount', NaN, ...
    'MeanBestScore', NaN, ...
    'MeanScoreMargin', NaN, ...
    'MinBestScore', NaN), numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    keep = HighSpeedNumbering.summary.SensorID == sid;
    if ~any(keep)
        continue;
    end
    rowNow = HighSpeedNumbering.summary(keep, :);
    rows(i).SensorID = sid;
    rows(i).DominantShift = rowNow.DominantShift(1);
    rows(i).LocalSlotOfB1 = rowNow.LocalSlotOfB1(1);
    rows(i).ValidRevolutionCount = rowNow.ValidRevolutionCount(1);
    rows(i).MeanBestScore = rowNow.MeanBestScore(1);
    rows(i).MeanScoreMargin = rowNow.MeanScoreMargin(1);
    rows(i).MinBestScore = rowNow.MinBestScore(1);
end

HighNumberingTable = struct2table(rows);
end

function AngleTable = build_angle_table_local(Sensor_Config, LowSpeedReference, A)
if isfield(LowSpeedReference, 'standard_angles_opr_center')
    angleMat = LowSpeedReference.standard_angles_opr_center;
elseif isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    angleMat = Sensor_Config.Standard_Relative_Angles_OPRCenter;
elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
    angleMat = Sensor_Config.Standard_Relative_Angles;
else
    error('No standard relative angle table found.');
end

sensorIds = A.compareSensors(:);
angleDeg = angleMat(sensorIds, A.bladeId);
relToAudit = wrap_to_180_local(angleDeg - angleDeg(sensorIds == A.auditSensor));
relToFirst = wrap_to_180_local(angleDeg - angleDeg(1));
AngleTable = table(sensorIds, angleDeg(:), relToAudit(:), relToFirst(:), ...
    'VariableNames', {'SensorID', 'B4AngleDeg', 'RelToCH3Deg', 'RelToFirstSensorDeg'});
end

function TemplateTable = attach_angle_to_template_table_local(TemplateTable, AngleTable)
[isIn, loc] = ismember(TemplateTable.SensorID, AngleTable.SensorID);
TemplateTable.AngleDeg(isIn) = AngleTable.B4AngleDeg(loc(isIn));
end

function AliasTable = build_alias_table_local(AngleTable, A)
combos = { ...
    [2 5 7], ...
    [3 5 7], ...
    [2 3 5 7], ...
    [2 3], ...
    [3 5], ...
    [3 7], ...
    [2 3 5], ...
    [2 3 7]};

rows = repmat(struct( ...
    'SensorCombo', "", ...
    'TargetEO', NaN, ...
    'NearestAliasEO', NaN, ...
    'NearestAliasDistanceDeg', NaN, ...
    'TargetSignatureDeg', "", ...
    'NearestAliasSignatureDeg', ""), numel(combos), 1);

for i = 1:numel(combos)
    combo = combos{i};
    [isIn, loc] = ismember(combo, AngleTable.SensorID);
    if ~all(isIn)
        continue;
    end
    angles = AngleTable.B4AngleDeg(loc);
    targetSig = phase_signature_deg_local(A.targetEO, angles);
    bestEO = NaN;
    bestDist = inf;
    bestSig = [];
    for eo = A.eoRange
        if eo == A.targetEO
            continue;
        end
        sig = phase_signature_deg_local(eo, angles);
        distNow = sqrt(mean(wrap_to_180_local(sig - targetSig) .^ 2, 'omitnan'));
        if distNow < bestDist
            bestDist = distNow;
            bestEO = eo;
            bestSig = sig;
        end
    end
    rows(i).SensorCombo = string(mat2str(combo));
    rows(i).TargetEO = A.targetEO;
    rows(i).NearestAliasEO = bestEO;
    rows(i).NearestAliasDistanceDeg = bestDist;
    rows(i).TargetSignatureDeg = string(vec_to_text_local(targetSig, '%.1f'));
    rows(i).NearestAliasSignatureDeg = string(vec_to_text_local(bestSig, '%.1f'));
end

AliasTable = struct2table(rows);
end

function EtaSummaryTable = summarize_existing_eta_results_local(outputDir, A)
combos = { ...
    [2 5 7], ...
    [3 5 7], ...
    [2 3 5 7]};
comboTags = strings(numel(combos), 1);
for i = 1:numel(combos)
    comboTags(i) = string(sensor_tag_local(combos{i}));
end

rows = repmat(struct( ...
    'SensorCombo', "", ...
    'SensorID', NaN, ...
    'WindowCount', NaN, ...
    'EO12Count', NaN, ...
    'ModeEO', NaN, ...
    'MeanEtaMM', NaN, ...
    'MedianEtaMM', NaN, ...
    'MaxAbsEtaMM', NaN), 0, 1);

for ic = 1:numel(combos)
    combo = combos{ic};
    tag = char(comboTags(ic));
    resultFile = fullfile(outputDir, '06_identification', tag, sprintf( ...
        'IdentificationResult_%s_20241106_B%d_%s_L03_eta200.mat', ...
        time_label_local(A.startTimeSec), A.bladeId, tag));
    if exist(resultFile, 'file') ~= 2
        continue;
    end
    loaded = load(resultFile, 'IdentificationResult');
    R = loaded.IdentificationResult;
    if ~isfield(R, 'WindowResult') || isempty(R.WindowResult)
        continue;
    end
    nWin = numel(R.WindowResult);
    eo = nan(nWin, 1);
    eta = nan(nWin, numel(combo));
    for iw = 1:nWin
        eo(iw) = R.WindowResult(iw).Result.EO_id;
        etaNow = R.WindowResult(iw).Result.sensor_eta_id;
        eta(iw, 1:min(numel(combo), numel(etaNow))) = etaNow(1:min(numel(combo), numel(etaNow)));
    end
    for is = 1:numel(combo)
        rows(end + 1, 1) = struct( ... %#ok<AGROW>
            'SensorCombo', string(tag), ...
            'SensorID', combo(is), ...
            'WindowCount', nWin, ...
            'EO12Count', sum(eo == A.targetEO), ...
            'ModeEO', mode(eo), ...
            'MeanEtaMM', mean(eta(:, is), 'omitnan'), ...
            'MedianEtaMM', median(eta(:, is), 'omitnan'), ...
            'MaxAbsEtaMM', max(abs(eta(:, is)), [], 'omitnan'));
    end
end

EtaSummaryTable = struct2table(rows);
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function WaveformResidualTable = summarize_waveform_residual_local(TemplateTable, WaveformLibrary, A)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'LapCountUsed', NaN, ...
    'PointCount', NaN, ...
    'RMSEV', NaN, ...
    'MedianAbsErrorV', NaN, ...
    'BiasV', NaN, ...
    'XMinMM', NaN, ...
    'XMaxMM', NaN), numel(A.compareSensors), 1);

B = WaveformLibrary.Blade(A.bladeId);
for i = 1:numel(A.compareSensors)
    sid = A.compareSensors(i);
    sensorIdx = find([B.Sensor.sensor_id] == sid, 1);
    tplRow = TemplateTable(TemplateTable.SensorID == sid, :);
    if isempty(sensorIdx) || isempty(tplRow)
        continue;
    end
    loaded = load(char(tplRow.TemplateFile), 'sensor_template');
    T = loaded.sensor_template.Sensor;
    Hs = B.Sensor(sensorIdx);
    lapIds = 1:min(A.maxHighSpeedLapsToPlot, numel(Hs.Lap));
    errAll = [];
    xAll = [];
    for il = 1:numel(lapIds)
        Lap = Hs.Lap(lapIds(il));
        x = Lap.x_rel(:);
        v = Lap.V(:);
        vTpl = interp1(T.x_grid(:), T.v_grid(:), x, 'linear', NaN);
        keep = isfinite(x) & isfinite(v) & isfinite(vTpl) & x >= T.x_domain(1) & x <= T.x_domain(2);
        errAll = [errAll; v(keep) - vTpl(keep)]; %#ok<AGROW>
        xAll = [xAll; x(keep)]; %#ok<AGROW>
    end
    rows(i).SensorID = sid;
    rows(i).LapCountUsed = numel(lapIds);
    rows(i).PointCount = numel(errAll);
    rows(i).RMSEV = sqrt(mean(errAll .^ 2, 'omitnan'));
    rows(i).MedianAbsErrorV = median(abs(errAll), 'omitnan');
    rows(i).BiasV = mean(errAll, 'omitnan');
    rows(i).XMinMM = min(xAll, [], 'omitnan');
    rows(i).XMaxMM = max(xAll, [], 'omitnan');
end

WaveformResidualTable = struct2table(rows);
end

function sigDeg = phase_signature_deg_local(eo, anglesDeg)
refDeg = anglesDeg(1);
sigDeg = wrap_to_180_local(eo * (anglesDeg(2:end) - refDeg));
end

function txt = vec_to_text_local(x, fmt)
if isempty(x)
    txt = "";
    return;
end
parts = compose(fmt, x(:).');
txt = strjoin(cellstr(parts), ', ');
end

function y = wrap_to_180_local(x)
y = mod(x + 180, 360) - 180;
end

function print_ch3_interpretation_local(TemplateTable, LowNumberingTable, HighNumberingTable, ...
        AngleTable, AliasTable, EtaSummaryTable, WaveformResidualTable, EtaInitialTable, A)
sid = A.auditSensor;
tpl = TemplateTable(TemplateTable.SensorID == sid, :);
low = LowNumberingTable(LowNumberingTable.SensorID == sid, :);
high = HighNumberingTable(HighNumberingTable.SensorID == sid, :);
fprintf('\n--- CH%d quick interpretation ---\n', sid);
fprintf('Template: width %.3f mm, xc %.4f mm, threshold %.3f V, V range %.3f V.\n', ...
    tpl.WidthMM, tpl.XcMM, tpl.ThresholdV, tpl.VRange);
fprintf('Numbering: low shift/localB1 = %d/%d; high shift/localB1 = %d/%d.\n', ...
    low.DominantShift, low.LocalSlotOfB1, high.DominantShift, high.LocalSlotOfB1);
if low.DominantShift == high.DominantShift && low.LocalSlotOfB1 == high.LocalSlotOfB1
    fprintf('Low/high numbering is internally consistent for CH%d.\n', sid);
else
    fprintf('WARNING: low/high numbering differs for CH%d.\n', sid);
end

ang = AngleTable(AngleTable.SensorID == sid, :);
fprintf('Mounting angle for B%d CH%d = %.3f deg under OPR-center reference.\n', ...
    A.bladeId, sid, ang.B4AngleDeg);

keepCH3Combos = contains(AliasTable.SensorCombo, sprintf('%d', sid));
badAlias = AliasTable(keepCH3Combos, :);
badAlias = sortrows(badAlias, 'NearestAliasDistanceDeg', 'ascend');
if ~isempty(badAlias)
    fprintf('Closest CH%d-including EO%d alias: combo %s -> EO%d, %.2f deg away.\n', ...
        sid, A.targetEO, badAlias.SensorCombo(1), ...
        badAlias.NearestAliasEO(1), badAlias.NearestAliasDistanceDeg(1));
end
if ~isempty(EtaSummaryTable)
    etaRows = EtaSummaryTable(EtaSummaryTable.SensorID == sid, :);
    if ~isempty(etaRows)
        fprintf('CH%d eta across existing combo results: maxAbs range %.4f to %.4f mm.\n', ...
            sid, min(etaRows.MaxAbsEtaMM, [], 'omitnan'), max(etaRows.MaxAbsEtaMM, [], 'omitnan'));
    end
end
if ~isempty(WaveformResidualTable)
    res = WaveformResidualTable(WaveformResidualTable.SensorID == sid, :);
    if ~isempty(res)
        fprintf('CH%d direct low/high residual: RMSE %.4f V, median |err| %.4f V, bias %.4f V.\n', ...
            sid, res.RMSEV, res.MedianAbsErrorV, res.BiasV);
    end
end
if ~isempty(EtaInitialTable)
    etaInit = EtaInitialTable(EtaInitialTable.SensorID == sid, :);
    if ~isempty(etaInit)
        fprintf(['CH%d eta initial guess from template-inverse alignment: median %.4f mm, ' ...
            'IQR [%.4f, %.4f] mm.\n'], ...
            sid, etaInit.MedianEtaGuessMM, etaInit.Q25EtaGuessMM, etaInit.Q75EtaGuessMM);
    end
end
fprintf(['Interpretation cue: if CH%d alone has good numbering scores but CH%d-combos lose EO%d, ' ...
    'the next suspect is waveform residual/eta compatibility, not low/high blade numbering.\n'], ...
    sid, sid, A.targetEO);
end

function plot_template_comparison_local(TemplateTable, A)
fig = figure('Name', sprintf('Step06H B%d CH%d template health', A.bladeId, A.auditSensor), ...
    'Color', 'w', 'NumberTitle', 'off', 'Units', 'centimeters', ...
    'Position', [1, 1, 22, 13]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
for i = 1:height(TemplateTable)
    loaded = load(char(TemplateTable.TemplateFile(i)), 'sensor_template');
    T = loaded.sensor_template.Sensor;
    if TemplateTable.SensorID(i) == A.auditSensor
        lw = 2.0;
    else
        lw = 1.0;
    end
    plot(T.x_grid, T.v_grid, 'LineWidth', lw, ...
        'DisplayName', sprintf('CH%d', TemplateTable.SensorID(i)));
end
xlabel('x_{rel} (mm)');
ylabel('Voltage (V)');
title(sprintf('B%d low-speed templates, raw voltage', A.bladeId), 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on;

nexttile;
bar(categorical(compose('CH%d', TemplateTable.SensorID)), TemplateTable.WidthMM);
ylabel('Template domain width (mm)');
title('Template width comparison', 'FontWeight', 'normal');
grid on; box on;
end

function plot_low_high_overlay_local(TemplateTable, WaveformLibrary, A)
sensorIdx = find([WaveformLibrary.Blade(A.bladeId).Sensor.sensor_id] == A.auditSensor, 1);
if isempty(sensorIdx)
    fprintf('CH%d not found in waveform library. Skip overlay.\n', A.auditSensor);
    return;
end
loaded = load(char(TemplateTable.TemplateFile(TemplateTable.SensorID == A.auditSensor)), 'sensor_template');
Tpl = loaded.sensor_template.Sensor;
Hs = WaveformLibrary.Blade(A.bladeId).Sensor(sensorIdx);
lapIds = 1:min(A.maxHighSpeedLapsToPlot, numel(Hs.Lap));

fig = figure('Name', sprintf('Step06H B%d CH%d low/high overlay', A.bladeId, A.auditSensor), ...
    'Color', 'w', 'NumberTitle', 'off', 'Units', 'centimeters', ...
    'Position', [2, 2, 23, 10]);
hold on;
plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.8, 'DisplayName', 'low-speed template');
cmap = lines(max(1, numel(lapIds)));
for i = 1:numel(lapIds)
    Lap = Hs.Lap(lapIds(i));
    [xRel, order] = sort(Lap.x_rel(:));
    plot(xRel, Lap.V(order), '-', 'Color', cmap(i, :), 'LineWidth', 0.9, ...
        'DisplayName', sprintf('high lap %d', lapIds(i)));
end
xline(Tpl.x_domain(1), 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xline(Tpl.x_domain(2), 'b-.', 'LineWidth', 1.0, 'DisplayName', 'template domain');
xline(0, 'k--', 'LineWidth', 0.8, 'DisplayName', 'template center');
xlabel('x_{rel} (mm)');
ylabel('Voltage (V)');
title(sprintf('B%d CH%d template vs high-speed waveform, thetaStd %.3f deg', ...
    A.bladeId, A.auditSensor, Hs.theta_std_deg), 'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
grid on; box on;
end

function EtaInitialTable = summarize_eta_initial_guess_local(TemplateTable, WaveformLibrary, A)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'LapCountUsed', NaN, ...
    'PointCount', NaN, ...
    'MeanEtaGuessMM', NaN, ...
    'MedianEtaGuessMM', NaN, ...
    'Q25EtaGuessMM', NaN, ...
    'Q75EtaGuessMM', NaN, ...
    'MaxAbsEtaGuessMM', NaN), numel(A.compareSensors), 1);

B = WaveformLibrary.Blade(A.bladeId);
for i = 1:numel(A.compareSensors)
    sid = A.compareSensors(i);
    sensorIdx = find([B.Sensor.sensor_id] == sid, 1);
    tplRow = TemplateTable(TemplateTable.SensorID == sid, :);
    if isempty(sensorIdx) || isempty(tplRow)
        continue;
    end
    loaded = load(char(tplRow.TemplateFile), 'sensor_template');
    T = loaded.sensor_template.Sensor;
    Hs = B.Sensor(sensorIdx);
    lapIds = 1:min(A.maxHighSpeedLapsToPlot, numel(Hs.Lap));
    etaGuessAll = [];
    for il = 1:numel(lapIds)
        Lap = Hs.Lap(lapIds(il));
        x = Lap.x_rel(:);
        v = Lap.V(:);
        xTpl = invert_template_voltage_local(T, v, x);
        etaGuess = x - xTpl;
        keep = isfinite(etaGuess) & isfinite(x) & isfinite(v) & ...
            x >= T.x_domain(1) & x <= T.x_domain(2);
        etaGuessAll = [etaGuessAll; etaGuess(keep)]; %#ok<AGROW>
    end
    rows(i).SensorID = sid;
    rows(i).LapCountUsed = numel(lapIds);
    rows(i).PointCount = numel(etaGuessAll);
    rows(i).MeanEtaGuessMM = mean(etaGuessAll, 'omitnan');
    rows(i).MedianEtaGuessMM = median(etaGuessAll, 'omitnan');
    rows(i).Q25EtaGuessMM = prctile(etaGuessAll, 25);
    rows(i).Q75EtaGuessMM = prctile(etaGuessAll, 75);
    rows(i).MaxAbsEtaGuessMM = max(abs(etaGuessAll), [], 'omitnan');
end

EtaInitialTable = struct2table(rows);
end

function xStat = invert_template_voltage_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diffV = vGrid - vv;
    crossingX = [];
    exactIdx = find(abs(diffV) <= 1e-10);
    if ~isempty(exactIdx)
        crossingX = xGrid(exactIdx);
    end
    for k = 1:numel(diffV) - 1
        if ~isfinite(diffV(k)) || ~isfinite(diffV(k + 1))
            continue;
        end
        if diffV(k) == 0 || diffV(k) * diffV(k + 1) > 0
            continue;
        end
        denom = vGrid(k + 1) - vGrid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - vGrid(k)) / denom;
        crossingX(end + 1, 1) = xGrid(k) + alpha * (xGrid(k + 1) - xGrid(k)); %#ok<AGROW>
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function plot_eta_summary_local(EtaSummaryTable, A)
fig = figure('Name', sprintf('Step06H B%d eta summary', A.bladeId), ...
    'Color', 'w', 'NumberTitle', 'off', 'Units', 'centimeters', ...
    'Position', [4, 4, 23, 10]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
groups = categorical(strcat(EtaSummaryTable.SensorCombo, "_CH", string(EtaSummaryTable.SensorID)));
bar(groups, EtaSummaryTable.MaxAbsEtaMM);
ylabel('max |eta| (mm)');
title('Static sensor offset used by existing eta-enabled runs', 'FontWeight', 'normal');
grid on; box on;

nexttile;
comboRows = unique(EtaSummaryTable(:, {'SensorCombo', 'EO12Count', 'WindowCount', 'ModeEO'}), 'rows');
bar(categorical(comboRows.SensorCombo), comboRows.EO12Count ./ comboRows.WindowCount);
ylabel(sprintf('EO%d window ratio', A.targetEO));
xlabel('Sensor combo');
ylim([0 1]);
title('Existing identification success ratio', 'FontWeight', 'normal');
grid on; box on;
end

function plot_angle_signature_local(AngleTable, AliasTable, A)
fig = figure('Name', sprintf('Step06H B%d angle and EO alias check', A.bladeId), ...
    'Color', 'w', 'NumberTitle', 'off', 'Units', 'centimeters', ...
    'Position', [3, 3, 23, 12]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
stem(AngleTable.SensorID, AngleTable.B4AngleDeg, 'filled', 'LineWidth', 1.2);
xlabel('Sensor ID');
ylabel('OPR-center angle (deg)');
title(sprintf('B%d sensor mounting angles', A.bladeId), 'FontWeight', 'normal');
grid on; box on;

nexttile;
AliasTable = sortrows(AliasTable, 'NearestAliasDistanceDeg', 'ascend');
bar(categorical(AliasTable.SensorCombo), AliasTable.NearestAliasDistanceDeg);
ylabel(sprintf('Nearest alias distance to EO%d (deg)', A.targetEO));
xlabel('Sensor combo');
title('Smaller distance means stronger EO alias risk from mounting geometry alone', ...
    'FontWeight', 'normal');
grid on; box on;
end

