%% Step06C_ScanSensorCombinations_20241106
% Diagnose how sensor combinations affect Step06 EO stability.
%
% This driver temporarily edits the Step06 local parameter block, runs the
% same identification method for several sensor combinations, and summarizes
% the EO sequence, RMSE, and fitted static sensor offsets.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

routeDir = fileparts(mfilename('fullpath'));
step06File = fullfile(routeDir, 'Step06_RunIdentificationByBlade_20241106.m');

%% Parameters to tune
S06C = struct();
S06C.sensorSets = {[5 7], [2 5 7], [3 5 7], [2 3 5 7]};
S06C.windowLaps = 3;
S06C.sensorEtaLimitListMM = [0 0.20];
S06C.sensorEtaRegWeightVPerMM = 0.02;
S06C.targetEO = 12;
S06C.startTimeSec = 75.0;
S06C.blades = 4;
S06C.viewEnable = false;

fprintf('\n=== Step06C: scan sensor combinations ===\n');
fprintf('windowLaps: %d\n', S06C.windowLaps);
fprintf('sensorEtaLimitMM: %s\n', mat2str(S06C.sensorEtaLimitListMM));
fprintf('Target EO: %d\n', S06C.targetEO);
drawnow;

originalText = fileread(step06File);
cleanupObj = onCleanup(@() restore_file_text_local(step06File, originalText)); %#ok<NASGU>

nCase = numel(S06C.sensorSets) * numel(S06C.sensorEtaLimitListMM);
summaryRows = repmat(empty_summary_row_local(), nCase, 1);
caseIndex = 0;

for iset = 1:numel(S06C.sensorSets)
    sensors = S06C.sensorSets{iset};
    for ie = 1:numel(S06C.sensorEtaLimitListMM)
        etaLimit = S06C.sensorEtaLimitListMM(ie);
        caseIndex = caseIndex + 1;
        sensorTag = sensor_tag_local(sensors);
        outputLabel = sprintf('B4_%s_L%02d_eta%03d', sensorTag, S06C.windowLaps, round(etaLimit * 1000));

        fprintf('\n--- Case %d/%d: sensors=%s, etaLimit=%.3f mm ---\n', ...
            caseIndex, nCase, mat2str(sensors), etaLimit);
        drawnow;

        write_step06_case_local(step06File, originalText, S06C, sensors, etaLimit, outputLabel);
        STEP06_SKIP_CLEAR_FOR_DRIVER = true; %#ok<NASGU>
        Step06_RunIdentificationByBlade_20241106;

        summaryRows(caseIndex) = summarize_step06_case_local(routeDir, S06C, sensors, outputLabel, etaLimit);
        fprintf('Result %s eta=%.3f: targetEO=%d/%d, dominantEO=%d, meanRMSE=%.5f V, etaMax=%.4f mm\n', ...
            sensorTag, etaLimit, summaryRows(caseIndex).TargetEOCount, ...
            summaryRows(caseIndex).WindowCount, summaryRows(caseIndex).DominantEO, ...
            summaryRows(caseIndex).MeanRMSE, summaryRows(caseIndex).EtaMaxAbsMM);
        fprintf('EO sequence: %s\n', summaryRows(caseIndex).EOSequence);
        drawnow;
    end
end

SensorComboSummary = struct2table(summaryRows);
outDir = fullfile(routeDir, 'output', 'new_flow', '06_identification', 'scan');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outFile = fullfile(outDir, 'Step06C_SensorCombinationScan_20241106.csv');
writetable(SensorComboSummary, outFile);

fprintf('\n=== Step06C summary ===\n');
disp(SensorComboSummary(:, {'SensorTag','SensorEtaLimitMM','WindowCount','TargetEOCount', ...
    'TargetEORatio','DominantEO','MeanRMSE','EtaMaxAbsMM','EtaHitLimit'}));
fprintf('Saved sensor-combination scan summary: %s\n', outFile);

function write_step06_case_local(step06File, originalText, S06C, sensors, etaLimit, outputLabel)
newText = originalText;
newText = regex_replace_once_local(newText, 'S06\.analysisSensors\s*=\s*\[[^\]]*\];', ...
    sprintf('S06.analysisSensors = [%s];', num2str(sensors)));
newText = regex_replace_once_local(newText, 'S06\.startTimeSec\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.startTimeSec = %.6g;', S06C.startTimeSec));
newText = regex_replace_once_local(newText, 'S06\.blades\s*=\s*[^;]+;', ...
    sprintf('S06.blades = %s;', mat2str(S06C.blades)));
newText = regex_replace_once_local(newText, 'S06\.windowLaps\s*=\s*\d+;', ...
    sprintf('S06.windowLaps = %d;', S06C.windowLaps));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaLimitMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaLimitMM = %.6g;', etaLimit));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaRegWeightVPerMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaRegWeightVPerMM = %.6g;', S06C.sensorEtaRegWeightVPerMM));
newText = regex_replace_once_local(newText, 'S06\.outputLabel\s*=\s*''[^'']*'';', ...
    sprintf('S06.outputLabel = ''%s'';', outputLabel));
newText = regex_replace_once_local(newText, 'S06\.viewEnable\s*=\s*(true|false);', ...
    sprintf('S06.viewEnable = %s;', logical_text_local(S06C.viewEnable)));

fid = fopen(step06File, 'w');
if fid < 0
    error('Could not open Step06 file for writing: %s', step06File);
end
fwrite(fid, newText, 'char');
fclose(fid);
clear('Step06_RunIdentificationByBlade_20241106');
end

function textOut = regex_replace_once_local(textIn, pattern, replacement)
matches = regexp(textIn, pattern, 'match');
if numel(matches) ~= 1
    error('Could not uniquely replace pattern: %s', pattern);
end
textOut = regexprep(textIn, pattern, replacement, 'once');
end

function restore_file_text_local(step06File, originalText)
fid = fopen(step06File, 'w');
if fid < 0
    warning('Could not restore Step06 file after scan: %s', step06File);
    return;
end
fwrite(fid, originalText, 'char');
fclose(fid);
clear('Step06_RunIdentificationByBlade_20241106');
end

function text = logical_text_local(value)
if value
    text = 'true';
else
    text = 'false';
end
end

function row = summarize_step06_case_local(routeDir, S06C, sensors, outputLabel, etaLimit)
sensorTag = sensor_tag_local(sensors);
timeLabel = time_label_local(S06C.startTimeSec);
resultFile = fullfile(routeDir, 'output', 'new_flow', '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106_%s.mat', timeLabel, outputLabel));
if exist(resultFile, 'file') ~= 2
    error('Missing Step06 result after scan: %s', resultFile);
end
loaded = load(resultFile, 'IdentificationResult');
R = loaded.IdentificationResult;
T = R.Trend;
eo = T.EO(:);
uniqueEO = unique(eo(:).');
counts = arrayfun(@(x) nnz(eo == x), uniqueEO);
[~, imax] = max(counts);
etaMax = 0;
etaMedianBySensor = "";
if isfield(R, 'WindowResult') && ~isempty(R.WindowResult)
    etaAll = [];
    for i = 1:numel(R.WindowResult)
        if isfield(R.WindowResult(i).Result, 'sensor_eta_id') && ~isempty(R.WindowResult(i).Result.sensor_eta_id)
            eta = R.WindowResult(i).Result.sensor_eta_id(:).';
            etaAll = [etaAll; eta]; %#ok<AGROW>
            etaMax = max(etaMax, max(abs(eta), [], 'omitnan'));
        end
    end
    if ~isempty(etaAll)
        etaMed = median(etaAll, 1, 'omitnan');
        etaMedianBySensor = string(mat2str(etaMed, 4));
    end
end

row = empty_summary_row_local();
row.SensorTag = string(sensorTag);
row.SensorSet = string(mat2str(sensors));
row.WindowLaps = S06C.windowLaps;
row.SensorEtaLimitMM = etaLimit;
row.WindowCount = height(T);
row.TargetEOCount = nnz(eo == S06C.targetEO);
row.TargetEORatio = row.TargetEOCount / max(row.WindowCount, 1);
row.DominantEO = uniqueEO(imax);
row.MeanRMSE = mean(T.WeightedVoltageRMSE, 'omitnan');
row.MedianRMSE = median(T.WeightedVoltageRMSE, 'omitnan');
row.BestRMSE = min(T.WeightedVoltageRMSE, [], 'omitnan');
row.WorstRMSE = max(T.WeightedVoltageRMSE, [], 'omitnan');
row.EtaMaxAbsMM = etaMax;
row.EtaHitLimit = etaLimit > 0 && etaMax >= 0.98 * etaLimit;
row.EtaMedianBySensorMM = etaMedianBySensor;
row.OutputLabel = string(outputLabel);
row.ResultFile = string(resultFile);
row.EOSequence = string(strjoin(string(eo(:).'), ','));
end

function row = empty_summary_row_local()
row = struct( ...
    'SensorTag', "", ...
    'SensorSet', "", ...
    'WindowLaps', NaN, ...
    'SensorEtaLimitMM', NaN, ...
    'WindowCount', NaN, ...
    'TargetEOCount', NaN, ...
    'TargetEORatio', NaN, ...
    'DominantEO', NaN, ...
    'MeanRMSE', NaN, ...
    'MedianRMSE', NaN, ...
    'BestRMSE', NaN, ...
    'WorstRMSE', NaN, ...
    'EtaMaxAbsMM', NaN, ...
    'EtaHitLimit', false, ...
    'EtaMedianBySensorMM', "", ...
    'OutputLabel', "", ...
    'ResultFile', "", ...
    'EOSequence', "");
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

