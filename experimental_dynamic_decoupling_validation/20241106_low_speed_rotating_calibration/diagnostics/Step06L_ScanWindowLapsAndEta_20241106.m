%% Step06L_ScanWindowLapsAndEta_20241106
% Diagnose whether Step06 EO instability comes from too few laps or from
% missing static sensor-offset freedom.
%
% This is a diagnostic driver. It temporarily edits the local parameter
% block in Step06_RunIdentificationByBlade_20241106, runs Step06, summarizes
% the saved result, and restores Step06 at the end.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

routeDir = fileparts(mfilename('fullpath'));
step06File = fullfile(routeDir, 'Step06_RunIdentificationByBlade_20241106.m');

%% Parameters to tune
S06L = struct();
S06L.windowLapsList = [3 5 7];
S06L.sensorEtaLimitListMM = [0 0.20];
S06L.sensorEtaRegWeightVPerMM = 0.02;
S06L.targetEO = 12;
S06L.startTimeSec = 75.0;
S06L.analysisSensors = [2 5 7];
S06L.blades = 4;
S06L.viewEnable = false;

fprintf('\n=== Step06L: scan window laps and sensor eta ===\n');
fprintf('windowLaps: %s\n', mat2str(S06L.windowLapsList));
fprintf('sensorEtaLimitMM: %s\n', mat2str(S06L.sensorEtaLimitListMM));
fprintf('Target EO: %d\n', S06L.targetEO);
drawnow;

originalText = fileread(step06File);
cleanupObj = onCleanup(@() restore_file_text_local(step06File, originalText)); %#ok<NASGU>

nCase = numel(S06L.windowLapsList) * numel(S06L.sensorEtaLimitListMM);
summaryRows = repmat(empty_summary_row_local(), nCase, 1);
caseIndex = 0;

for il = 1:numel(S06L.windowLapsList)
    windowLaps = S06L.windowLapsList(il);
    for ie = 1:numel(S06L.sensorEtaLimitListMM)
        etaLimit = S06L.sensorEtaLimitListMM(ie);
        caseIndex = caseIndex + 1;
        outputLabel = sprintf('B4_L%02d_eta%03d', windowLaps, round(etaLimit * 1000));

        fprintf('\n--- Case %d/%d: windowLaps=%d, etaLimit=%.3f mm ---\n', ...
            caseIndex, nCase, windowLaps, etaLimit);
        drawnow;

        write_step06_case_local(step06File, originalText, S06L, windowLaps, etaLimit, outputLabel);
        STEP06_SKIP_CLEAR_FOR_DRIVER = true; %#ok<NASGU>
        Step06_RunIdentificationByBlade_20241106;

        summaryRows(caseIndex) = summarize_step06_case_local(routeDir, S06L, outputLabel, windowLaps, etaLimit);
        fprintf('Result: targetEO=%d/%d, dominantEO=%d, meanRMSE=%.5f V, etaMax=%.4f mm\n', ...
            summaryRows(caseIndex).TargetEOCount, summaryRows(caseIndex).WindowCount, ...
            summaryRows(caseIndex).DominantEO, summaryRows(caseIndex).MeanRMSE, ...
            summaryRows(caseIndex).EtaMaxAbsMM);
        fprintf('EO sequence: %s\n', summaryRows(caseIndex).EOSequence);
        drawnow;
    end
end

ScanSummary = struct2table(summaryRows);
outDir = fullfile(routeDir, 'output', 'new_flow', '06_identification', 'scan');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outFile = fullfile(outDir, 'Step06L_WindowLapsEtaScan_20241106.csv');
writetable(ScanSummary, outFile);

fprintf('\n=== Step06L summary ===\n');
disp(ScanSummary(:, {'WindowLaps','SensorEtaLimitMM','WindowCount','TargetEOCount', ...
    'TargetEORatio','DominantEO','MeanRMSE','EtaMaxAbsMM','EtaHitLimit'}));
fprintf('Saved scan summary: %s\n', outFile);

function write_step06_case_local(step06File, originalText, S06L, windowLaps, etaLimit, outputLabel)
newText = originalText;
newText = regex_replace_once_local(newText, 'S06\.analysisSensors\s*=\s*\[[^\]]*\];', ...
    sprintf('S06.analysisSensors = [%s];', num2str(S06L.analysisSensors)));
newText = regex_replace_once_local(newText, 'S06\.startTimeSec\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.startTimeSec = %.6g;', S06L.startTimeSec));
newText = regex_replace_once_local(newText, 'S06\.blades\s*=\s*[^;]+;', ...
    sprintf('S06.blades = %s;', mat2str(S06L.blades)));
newText = regex_replace_once_local(newText, 'S06\.windowLaps\s*=\s*\d+;', ...
    sprintf('S06.windowLaps = %d;', windowLaps));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaLimitMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaLimitMM = %.6g;', etaLimit));
newText = regex_replace_once_local(newText, 'S06\.sensorEtaRegWeightVPerMM\s*=\s*[-+]?\d+(\.\d+)?;', ...
    sprintf('S06.sensorEtaRegWeightVPerMM = %.6g;', S06L.sensorEtaRegWeightVPerMM));
newText = regex_replace_once_local(newText, 'S06\.outputLabel\s*=\s*''[^'']*'';', ...
    sprintf('S06.outputLabel = ''%s'';', outputLabel));
newText = regex_replace_once_local(newText, 'S06\.viewEnable\s*=\s*(true|false);', ...
    sprintf('S06.viewEnable = %s;', logical_text_local(S06L.viewEnable)));

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
count = numel(matches);
if count ~= 1
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

function row = summarize_step06_case_local(routeDir, S06L, outputLabel, windowLaps, etaLimit)
sensorTag = ['S', sprintf('%d', S06L.analysisSensors)];
timeLabel = time_label_local(S06L.startTimeSec);
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
for i = 1:numel(R.WindowResult)
    if isfield(R.WindowResult(i).Result, 'sensor_eta_id') && ~isempty(R.WindowResult(i).Result.sensor_eta_id)
        etaMax = max(etaMax, max(abs(R.WindowResult(i).Result.sensor_eta_id), [], 'omitnan'));
    end
end

row = empty_summary_row_local();
row.WindowLaps = windowLaps;
row.SensorEtaLimitMM = etaLimit;
row.WindowCount = height(T);
row.TargetEOCount = nnz(eo == S06L.targetEO);
row.TargetEORatio = row.TargetEOCount / max(row.WindowCount, 1);
row.DominantEO = uniqueEO(imax);
row.MeanRMSE = mean(T.WeightedVoltageRMSE, 'omitnan');
row.MedianRMSE = median(T.WeightedVoltageRMSE, 'omitnan');
row.BestRMSE = min(T.WeightedVoltageRMSE, [], 'omitnan');
row.WorstRMSE = max(T.WeightedVoltageRMSE, [], 'omitnan');
row.EtaMaxAbsMM = etaMax;
row.EtaHitLimit = etaLimit > 0 && etaMax >= 0.98 * etaLimit;
row.OutputLabel = string(outputLabel);
row.ResultFile = string(resultFile);
row.EOSequence = string(strjoin(string(eo(:).'), ','));
end

function row = empty_summary_row_local()
row = struct( ...
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

function label = time_label_local(tSec)
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

