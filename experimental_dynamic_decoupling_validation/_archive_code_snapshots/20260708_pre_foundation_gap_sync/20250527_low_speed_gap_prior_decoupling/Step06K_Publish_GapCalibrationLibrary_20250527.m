%% Step06K: Publish existing Step06I gap calibration as blade-sensor entries
% This script does not recalibrate and does not change Step07J inputs. It
% reads the existing combined Step06I CorrectedGapLibrary and publishes
% per blade-sensor GapCalib_*.mat files plus an index CSV, matching the
% calibration-first workflow used by the identification route.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
sensorTag = C0.sensorTag;

sourceFile = find_corrected_gap_library_file_local(outDir, C0.bladeId, sensorTag);
if ~isfile(sourceFile)
    error('Missing combined Step06I calibration file: %s', sourceFile);
end

S = load(sourceFile, 'CorrectedGapLibrary', 'calibrationTable', 'highCheckTable');
if ~isfield(S, 'CorrectedGapLibrary')
    error('File does not contain CorrectedGapLibrary: %s', sourceFile);
end
CorrectedGapLibrary = S.CorrectedGapLibrary;
if isfield(S, 'calibrationTable')
    calibrationTable = S.calibrationTable;
else
    calibrationTable = build_calibration_table_from_library_local(CorrectedGapLibrary);
end
if isfield(S, 'highCheckTable')
    highCheckTable = S.highCheckTable;
else
    highCheckTable = table();
end

[indexTable, indexFile] = publish_gap_calibration_entries_local( ...
    CorrectedGapLibrary, calibrationTable, highCheckTable, outDir, C0, sourceFile);

fprintf('\nStep06K complete.\n');
fprintf('Source: %s\n', sourceFile);
fprintf('Index : %s\n', indexFile);
disp(indexTable);

function correctedLibFile = find_corrected_gap_library_file_local(outDir, targetBlade, sensorTag)
patterns = {
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s_OPRCenterStd.mat', targetBlade, sensorTag)
    sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s*.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    files = dir(fullfile(outDir, patterns{i}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        correctedLibFile = fullfile(files(idx).folder, files(idx).name);
        return;
    end
end
correctedLibFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20250527_B%d_%s.mat', ...
    targetBlade, sensorTag));
end

function T = build_calibration_table_from_library_local(Lib)
rows = cell(numel(Lib.sensor), 1);
for i = 1:numel(Lib.sensor)
    s = Lib.sensor(i);
    rows{i} = table(s.sensorId, s.g0Mm, s.tauMm, s.xScale, s.muGapPerXMm, ...
        s.tiltAngleDeg, s.voltageGain, s.voltageOffsetMv, s.lowFitRmseMv, ...
        s.overshootRmseMm, s.pointCount, s.etaMedianMm, s.etaIqrMm, s.etaLimitMm, ...
        'VariableNames', {'sensorId','g0Mm','tauMm','xScale','muGapPerXMm', ...
        'tiltAngleDeg','voltageGain','voltageOffsetMv','lowFitRmseMv', ...
        'overshootRmseMm','pointCount','etaMedianMm','etaIqrMm','etaLimitMm'});
end
T = vertcat(rows{:});
end

function [indexTable, indexFile] = publish_gap_calibration_entries_local( ...
        CorrectedGapLibrary, calibrationTable, highCheckTable, outDir, C0, sourceFile)
if isfield(C0, 'flowConfig') && isfield(C0.flowConfig, 'calibration')
    libDirName = C0.flowConfig.calibration.libraryDirName;
    filePrefix = C0.flowConfig.calibration.filePrefix;
else
    libDirName = 'gap_calibration_library';
    filePrefix = 'GapCalib';
end
libDir = fullfile(outDir, libDirName);
if exist(libDir, 'dir') ~= 7
    mkdir(libDir);
end

rows = cell(numel(CorrectedGapLibrary.sensor), 1);
for i = 1:numel(CorrectedGapLibrary.sensor)
    sensorCal = CorrectedGapLibrary.sensor(i);
    sid = sensorCal.sensorId;
    sensorFile = fullfile(libDir, sprintf('%s_%s_B%d_S%d.mat', ...
        filePrefix, C0.dataset, CorrectedGapLibrary.targetBlade, sid));

    SingleCorrectedGapLibrary = struct(); %#ok<NASGU>
    SingleCorrectedGapLibrary.dataset = CorrectedGapLibrary.dataset;
    SingleCorrectedGapLibrary.method = CorrectedGapLibrary.method;
    SingleCorrectedGapLibrary.description = CorrectedGapLibrary.description;
    SingleCorrectedGapLibrary.sourceCombinedFile = sourceFile;
    SingleCorrectedGapLibrary.responseFile = CorrectedGapLibrary.responseFile;
    SingleCorrectedGapLibrary.templateFile = CorrectedGapLibrary.templateFile;
    SingleCorrectedGapLibrary.targetBlade = CorrectedGapLibrary.targetBlade;
    SingleCorrectedGapLibrary.analysisSensors = sid;
    SingleCorrectedGapLibrary.sensor = sensorCal;

    calibRow = calibrationTable(calibrationTable.sensorId == sid, :);
    if isempty(highCheckTable)
        highCheckRow = table();
    else
        highCheckRow = highCheckTable(highCheckTable.sensorId == sid, :);
    end

    GapCalib = struct(); %#ok<NASGU>
    GapCalib.dataset = C0.dataset;
    GapCalib.bladeId = CorrectedGapLibrary.targetBlade;
    GapCalib.sensorId = sid;
    GapCalib.method = CorrectedGapLibrary.method;
    GapCalib.sourceCombinedFile = sourceFile;
    GapCalib.responseFile = CorrectedGapLibrary.responseFile;
    GapCalib.templateFile = CorrectedGapLibrary.templateFile;
    GapCalib.correctedGapLibrary = SingleCorrectedGapLibrary;
    GapCalib.sensor = sensorCal;
    GapCalib.calibrationRow = calibRow;
    GapCalib.highCheckRow = highCheckRow;

    save(sensorFile, 'GapCalib', 'SingleCorrectedGapLibrary', 'calibRow', 'highCheckRow', '-v7.3');
    rows{i} = table(string(C0.dataset), CorrectedGapLibrary.targetBlade, sid, ...
        string(sensorFile), string(sourceFile), sensorCal.g0Mm, sensorCal.tauMm, ...
        sensorCal.xScale, sensorCal.muGapPerXMm, sensorCal.lowFitRmseMv, ...
        'VariableNames', {'dataset','bladeId','sensorId','calibrationFile', ...
        'sourceCombinedFile','g0Mm','tauMm','xScale','muGapPerXMm','lowFitRmseMv'});
end

indexTable = vertcat(rows{:});
indexFile = fullfile(libDir, sprintf('%s_Index_%s_%s.csv', filePrefix, C0.dataset, C0.caseTag));
writetable(indexTable, indexFile);
end
