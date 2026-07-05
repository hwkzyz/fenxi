%% Step06D_ScanRequestedCombosWithEta_20241106
% Scan the requested sensor combinations with eta enabled and more laps.
%
% Focus:
%   S57 / S357 / S2357
%   windowLaps = 5, 7
%   sensorEtaLimitMM = 0.20

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

routeDir = fileparts(mfilename('fullpath'));
step06File = fullfile(routeDir, 'Step06_RunIdentificationByBlade_20241106.m');

%% Parameters to tune
S06D = struct();
S06D.sensorSets = {[5 7], [3 5 7], [2 3 5 7]};
S06D.windowLapsList = [5 7];
S06D.sensorEtaLimitMM = 0.20;
S06D.sensorEtaRegWeightVPerMM = 0.02;
S06D.targetEO = 12;
S06D.startTimeSec = 75.0;
S06D.blades = 4;
S06D.viewEnable = false;

fprintf('\n=== Step06D: scan requested combos with eta ===\n');
fprintf('windowLaps: %s\n', mat2str(S06D.windowLapsList));
fprintf('sensorSets: S57 / S357 / S2357\n');
fprintf('sensorEtaLimitMM: %.3f\n', S06D.sensorEtaLimitMM);
drawnow;

originalText = fileread(step06File);
cleanupObj = onCleanup(@() restore_file_text_local(step06File, originalText)); %#ok<NASGU>

nCase = numel(S06D.sensorSets) * numel(S06D.windowLapsList);
summaryRows = repmat(empty_summary_row_local(), nCase, 1);
caseIndex = 0;

for iset = 1:numel(S06D.sensorSets)
    sensors = S06D.sensorSets{iset};
    for il = 1:numel(S06D.windowLapsList)
        windowLaps = S06D.windowLapsList(il);
        caseIndex = caseIndex + 1;
        sensorTag = sensor_tag_local(sensors);
        outputLabel = sprintf('B4_%s_L%02d_eta200', sensorTag, windowLaps);

        fprintf('\n--- Case %d/%d: sensors=%s, windowLaps=%d ---\n', ...
            caseIndex, nCase, mat2str(sensors), windowLaps);
        drawnow;

        write_step06_case_local(step06File, originalText, S06D, sensors, windowLaps, outputLabel);
        STEP06_SKIP_CLEAR_FOR_DRIVER = true; %#ok<NASGU>
        Step06_RunIdentificationByBlade_20241106;

        summaryRows(caseIndex) = summarize_step06_case_local(routeDir, S06D, sensors, outputLabel, windowLaps);
        fprintf('Result %s L=%d: targetEO=%d/%d, dominantEO=%d, meanRMSE=%.5f V, etaMax=%.4f mm\n', ...
            sensorTag, windowLaps, summaryRows(caseIndex).TargetEOCount, ...
            summaryRows(caseIndex).WindowCount, summaryRows(caseIndex).DominantEO, ...
            summaryRows(caseIndex).MeanRMSE, summaryRows(caseIndex).EtaMaxAbsMM);
        fprintf('EO sequence: %s\n', summaryRows(caseIndex).EOSequence);
        drawnow;
    end
end

ScanSummary = struct2table(summaryRows);
outDir = fullfile(routeDir, 'output', 'new_flow', '06_identification', 'scan');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outFile = fullfile(outDir, 'Step06D_RequestedCombosWithEtaScan_20241106.csv');
writetable(ScanSummary, outFile);

fprintf('\n=== Step06D summary ===\n');
disp(ScanSummary(:, {'SensorTag','WindowLaps','WindowCount','TargetEOCount', ...
    'TargetEORatio','DominantEO','MeanRMSE','EtaMaxAbsMM','EtaHitLimit'}));
fprintf('Saved scan summary: %s\n', outFile);

function write_step06_case_local(step06File, originalText, S06D, sensors, windowLaps, outputLabel)
newText = originalText;
newText = regex_replace_once_local(newText, 'S06\.analysisSensors\s*=\s*\[[^\]]*\];', ...
    sprintf('S06.analysisSensors = [%s];', num2str(sensors)));
newText = regex_replace_once_local(newText, 'S06\.startTimeSec\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.startTimeSec = %.6g;', S06D.startTimeSec));
newText = regex_replace_once_local(newText, 'S06\.blades\s*=\s*[^;]+;', ...
    sprintf('S06.blades = %s;', mat2str(S06D.blades)));
newText = regex_replace_once_local(newText, 'S06\.windowLaps\s*=\s*\d+;', ...
    sprintf('S06.windowLaps = %d;', windowLaps));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaLimitMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaLimitMM = %.6g;', S06D.sensorEtaLimitMM));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaRegWeightVPerMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaRegWeightVPerMM = %.6g;', S06D.sensorEtaRegWeightVPerMM));
newText = regex_replace_once_local(newText, 'S06\.outputLabel\s*=\s*''[^'']*'';', ...
    sprintf('S06.outputLabel = ''%s'';', outputLabel));
newText = regex_replace_once_local(newText, 'S06\.viewEnable\s*=\s*(true|false);', ...
    sprintf('S06.viewEnable = %s;', logical_text_local(S06D.viewEnable)));

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

function row = summarize_step06_case_local(routeDir, S06D, sensors, outputLabel, windowLaps)
sensorTag = sensor_tag_local(sensors);
timeLabel = time_label_local(S06D.startTimeSec);
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
if isfield(R, 'WindowResult') && ~isempty(R.WindowResult)
    for i = 1:numel(R.WindowResult)
        if isfield(R.WindowResult(i).Result, 'sensor_eta_id') && ~isempty(R.WindowResult(i).Result.sensor_eta_id)
            eta = R.WindowResult(i).Result.sensor_eta_id(:).';
            etaMax = max(etaMax, max(abs(eta), [], 'omitnan'));
        end
    end
end

row = empty_summary_row_local();
row.SensorTag = string(sensorTag);
row.SensorSet = string(mat2str(sensors));
row.WindowLaps = windowLaps;
row.SensorEtaLimitMM = S06D.sensorEtaLimitMM;
row.WindowCount = height(T);
row.TargetEOCount = nnz(eo == S06D.targetEO);
row.TargetEORatio = row.TargetEOCount / max(row.WindowCount, 1);
row.DominantEO = uniqueEO(imax);
row.MeanRMSE = mean(T.WeightedVoltageRMSE, 'omitnan');
row.MedianRMSE = median(T.WeightedVoltageRMSE, 'omitnan');
row.BestRMSE = min(T.WeightedVoltageRMSE, [], 'omitnan');
row.WorstRMSE = max(T.WeightedVoltageRMSE, [], 'omitnan');
row.EtaMaxAbsMM = etaMax;
row.EtaHitLimit = etaMax >= 0.98 * S06D.sensorEtaLimitMM;
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

