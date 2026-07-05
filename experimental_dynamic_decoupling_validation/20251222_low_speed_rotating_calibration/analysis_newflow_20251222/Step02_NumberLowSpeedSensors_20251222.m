%% Step02_NumberLowSpeedSensors_20251222
% Build the low-speed physical blade numbering artifact for NewFlow.
%
% The 20251222 Sensor_Config already contains the accepted blade numbering
% from the mature preprocessing route. This step makes that numbering
% explicit: per sensor and per physical blade, it exports fingerprint
% values, OPR-center standard angles, and target pulse anchors.

clear; close all; clc;

%% Parameters to tune
P = NewFlow_Config_20251222();
analysisSensors = P.sensors.analysis;
bladeCount = P.machine.bladeCount;
viewEnable = true;
saveFigures = true;

%% Paths
lowSpeedReferenceFile = P.files.lowSpeedReference;
outFile = P.files.lowSpeedNumbering;
numberingCsv = strrep(outFile, '.mat', '_ByBlade.csv');
summaryCsv = strrep(outFile, '.mat', '_Summary.csv');
figureDir = fullfile(P.outputDir, 'figures', '02_low_speed_numbering');

ensure_parent_dir_local(outFile);
require_file_local(lowSpeedReferenceFile, 'Step01 low-speed reference fingerprint');

loaded = load(lowSpeedReferenceFile, 'LowSpeedReference');
LowSpeedReference = loaded.LowSpeedReference;
validate_low_speed_reference_local(LowSpeedReference, analysisSensors, bladeCount, lowSpeedReferenceFile);

[LowSpeedNumberingTable, LowSpeedNumberingSummary] = build_low_speed_numbering_local( ...
    LowSpeedReference, analysisSensors, bladeCount);

LowSpeedNumbering = struct();
LowSpeedNumbering.dataset = P.dataset;
LowSpeedNumbering.mode = 'explicit_physical_blade_numbering_from_sensor_config';
LowSpeedNumbering.reference_file = lowSpeedReferenceFile;
LowSpeedNumbering.analysis_sensors = analysisSensors;
LowSpeedNumbering.blade_count = bladeCount;
LowSpeedNumbering.numbering_table = LowSpeedNumberingTable;
LowSpeedNumbering.summary_table = LowSpeedNumberingSummary;
LowSpeedNumbering.standard_angles_opr_center = LowSpeedReference.standard_angles_opr_center;
LowSpeedNumbering.fingerprints = LowSpeedReference.fingerprints;
LowSpeedNumbering.target_indices = LowSpeedReference.target_indices;

save(outFile, 'LowSpeedNumbering', 'LowSpeedNumberingTable', 'LowSpeedNumberingSummary');
writetable(LowSpeedNumberingTable, numberingCsv);
writetable(LowSpeedNumberingSummary, summaryCsv);

fprintf('\n=== Step02: low-speed physical blade numbering ===\n');
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Source:  %s\n', lowSpeedReferenceFile);
fprintf('Saved:   %s\n', outFile);
fprintf('Table:   %s\n', numberingCsv);
disp(LowSpeedNumberingSummary);

if viewEnable
    visualize_low_speed_numbering_local( ...
        LowSpeedNumberingTable, LowSpeedNumberingSummary, bladeCount, saveFigures, figureDir);
end

function validate_low_speed_reference_local(LowSpeedReference, sensorIds, bladeCount, sourceFile)
requiredFields = {'fingerprints', 'target_indices', 'standard_angles_opr_center', ...
    'analysis_sensors', 'blade_count'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(LowSpeedReference, name)
        error('LowSpeedReference is missing field %s:\n  %s', name, sourceFile);
    end
end
if ~isa(LowSpeedReference.fingerprints, 'containers.Map')
    error('LowSpeedReference.fingerprints must be containers.Map.');
end
if ~isa(LowSpeedReference.target_indices, 'containers.Map')
    error('LowSpeedReference.target_indices must be containers.Map.');
end
for sid = sensorIds(:).'
    if ~isKey(LowSpeedReference.fingerprints, sid)
        error('LowSpeedReference.fingerprints does not contain CH%d.', sid);
    end
    fp = LowSpeedReference.fingerprints(sid);
    if numel(fp) < bladeCount
        error('CH%d fingerprint has %d values; expected at least %d.', ...
            sid, numel(fp), bladeCount);
    end
end
end

function [numberingTable, summaryTable] = build_low_speed_numbering_local( ...
        LowSpeedReference, sensorIds, bladeCount)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'PhysicalBladeID', NaN, ...
    'FingerprintPeak', NaN, ...
    'FingerprintNormalized', NaN, ...
    'OPRCenterAngleDeg', NaN, ...
    'TargetIndexB1', NaN), numel(sensorIds) * bladeCount, 1);
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'TargetIndexB1', NaN, ...
    'MaxFingerprintPeak', NaN, ...
    'MinFingerprintPeak', NaN, ...
    'MeanFingerprintPeak', NaN, ...
    'AngleSpanDeg', NaN, ...
    'NumberedBladeCount', NaN), numel(sensorIds), 1);

rowIdx = 0;
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    fp = force_row_vector_local(LowSpeedReference.fingerprints(sid));
    fp = fp(1:bladeCount);
    angles = force_row_vector_local(LowSpeedReference.standard_angles_opr_center(sid, 1:bladeCount));
    fpMax = max(fp, [], 'omitnan');
    if ~isfinite(fpMax) || fpMax == 0
        fpMax = 1;
    end
    targetIndexB1 = NaN;
    if isKey(LowSpeedReference.target_indices, sid)
        targetIndexB1 = LowSpeedReference.target_indices(sid);
    end

    for bladeId = 1:bladeCount
        rowIdx = rowIdx + 1;
        rows(rowIdx).SensorID = sid;
        rows(rowIdx).PhysicalBladeID = bladeId;
        rows(rowIdx).FingerprintPeak = fp(bladeId);
        rows(rowIdx).FingerprintNormalized = fp(bladeId) / fpMax;
        rows(rowIdx).OPRCenterAngleDeg = angles(bladeId);
        rows(rowIdx).TargetIndexB1 = targetIndexB1;
    end

    summaryRows(is).SensorID = sid;
    summaryRows(is).TargetIndexB1 = targetIndexB1;
    summaryRows(is).MaxFingerprintPeak = max(fp, [], 'omitnan');
    summaryRows(is).MinFingerprintPeak = min(fp, [], 'omitnan');
    summaryRows(is).MeanFingerprintPeak = mean(fp, 'omitnan');
    summaryRows(is).AngleSpanDeg = max(angles, [], 'omitnan') - min(angles, [], 'omitnan');
    summaryRows(is).NumberedBladeCount = nnz(isfinite(fp) & isfinite(angles));
end

numberingTable = struct2table(rows);
summaryTable = struct2table(summaryRows);
end

function visualize_low_speed_numbering_local(numberingTable, summaryTable, bladeCount, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end
sensorIds = summaryTable.SensorID(:).';
fpMat = nan(numel(sensorIds), bladeCount);
angleMat = nan(numel(sensorIds), bladeCount);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    for bladeId = 1:bladeCount
        idx = numberingTable.SensorID == sid & numberingTable.PhysicalBladeID == bladeId;
        fpMat(is, bladeId) = numberingTable.FingerprintNormalized(idx);
        angleMat(is, bladeId) = numberingTable.OPRCenterAngleDeg(idx);
    end
end

fig = figure('Name', 'Step02 low-speed numbering 20251222', 'Color', 'w', ...
    'Position', [100 100 1200 740], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(1:bladeCount, sensorIds, fpMat);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Normalized fingerprint by physical blade');
colorbar; box on;

nexttile;
plot(1:bladeCount, angleMat.', '-o', 'LineWidth', 1.1);
xlabel('Physical blade ID');
ylabel('OPR-center angle (deg)');
title('Per-sensor standard angles');
legend(compose('CH%d', sensorIds), 'Location', 'best');
grid on; box on;

nexttile;
bar(categorical(compose('CH%d', sensorIds)), summaryTable.NumberedBladeCount);
ylabel('Numbered blade count');
title('Numbering completeness');
grid on; box on;

nexttile;
axis off;
txt = evalc('disp(summaryTable)');
text(0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Consolas', 'FontSize', 8.5, 'Interpreter', 'none');
title('Low-speed numbering summary');

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step02_LowSpeedNumbering_20251222.png'), ...
        'Resolution', 300);
end
end

function x = force_row_vector_local(x)
x = x(:).';
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function ensure_parent_dir_local(pathText)
folder = fileparts(pathText);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
