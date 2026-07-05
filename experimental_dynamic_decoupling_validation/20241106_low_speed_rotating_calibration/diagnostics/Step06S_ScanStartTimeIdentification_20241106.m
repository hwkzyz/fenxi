%% Step06S_ScanStartTimeIdentification_20241106
% Scan several high-speed start times with the current Step06 settings.
% Step06 remains the single implementation of the identification method.
% This script only edits the Step06 local start-time line, runs Step06, and
% summarizes the saved results.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

routeDir = fileparts(mfilename('fullpath'));
step06File = fullfile(routeDir, 'Step06_RunIdentificationByBlade_20241106.m');

%% Parameters to tune
startTimeListSec = 75.0:0.1:75.8;
targetEO = 12;
showPerWindowSequence = true;

fprintf('\n=== Step06S: scan Step06 startTimeSec ===\n');
fprintf('Start times: %s\n', mat2str(startTimeListSec));
fprintf('Target EO: %d\n', targetEO);
drawnow;

originalText = fileread(step06File);
cleanupObj = onCleanup(@() restore_file_text_local(step06File, originalText)); %#ok<NASGU>

summaryRows = repmat(struct( ...
    'StartTimeSec', NaN, ...
    'TimeLabel', "", ...
    'WindowCount', NaN, ...
    'TargetEOCount', NaN, ...
    'EO4Count', NaN, ...
    'EO5Count', NaN, ...
    'TargetEORatio', NaN, ...
    'DominantEO', NaN, ...
    'MeanRMSE', NaN, ...
    'MedianRMSE', NaN, ...
    'BestRMSE', NaN, ...
    'WorstRMSE', NaN, ...
    'EOSequence', ""), numel(startTimeListSec), 1);

for i = 1:numel(startTimeListSec)
    t0 = startTimeListSec(i);
    fprintf('\n--- Step06S scan %d/%d: startTimeSec = %.3f s ---\n', ...
        i, numel(startTimeListSec), t0);
    drawnow;

    write_step06_start_time_local(step06File, originalText, t0);
    Step06_RunIdentificationByBlade_20241106;

    row = summarize_latest_step06_result_local(routeDir, t0, targetEO);
    summaryRows(i) = row;
    fprintf('Summary %.3f s: targetEO=%d/%d, dominantEO=%d, meanRMSE=%.5f V\n', ...
        t0, row.TargetEOCount, row.WindowCount, row.DominantEO, row.MeanRMSE);
    if showPerWindowSequence
        fprintf('EO sequence: %s\n', row.EOSequence);
    end
    drawnow;
end

ScanSummary = struct2table(summaryRows);
outDir = fullfile(routeDir, 'output', 'new_flow', '06_identification', 'scan');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outFile = fullfile(outDir, 'Step06S_StartTimeScan_B4_20241106.csv');
writetable(ScanSummary, outFile);

fprintf('\n=== Step06S scan summary ===\n');
disp(ScanSummary(:, {'StartTimeSec','WindowCount','TargetEOCount','EO4Count','EO5Count','TargetEORatio','DominantEO','MeanRMSE','BestRMSE','WorstRMSE'}));
fprintf('Saved scan summary: %s\n', outFile);

function write_step06_start_time_local(step06File, originalText, startTimeSec)
pattern = 'S06\.startTimeSec\s*=\s*[-+]?\d+(\.\d+)?;';
replacement = sprintf('S06.startTimeSec = %.6g;', startTimeSec);
[newText, count] = regexprep(originalText, pattern, replacement, 'once');
if count ~= 1
    error('Could not uniquely update S06.startTimeSec in %s.', step06File);
end
fid = fopen(step06File, 'w');
if fid < 0
    error('Could not open Step06 file for writing: %s', step06File);
end
fwrite(fid, newText, 'char');
fclose(fid);
clear('Step06_RunIdentificationByBlade_20241106');
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

function row = summarize_latest_step06_result_local(routeDir, startTimeSec, targetEO)
sensorTag = 'S257';
timeLabel = time_label_local(startTimeSec);
summaryFile = fullfile(routeDir, 'output', 'new_flow', '06_identification', sensorTag, ...
    sprintf('IdentificationSummary_%s_20241106_B4_only.csv', timeLabel));
if exist(summaryFile, 'file') ~= 2
    error('Missing Step06 summary after scan: %s', summaryFile);
end
T = readtable(summaryFile);
eo = T.EO(:);
uniqueEO = unique(eo(:).');
counts = arrayfun(@(x) nnz(eo == x), uniqueEO);
[~, imax] = max(counts);
row = struct();
row.StartTimeSec = startTimeSec;
row.TimeLabel = string(timeLabel);
row.WindowCount = height(T);
row.TargetEOCount = nnz(eo == targetEO);
row.EO4Count = nnz(eo == 4);
row.EO5Count = nnz(eo == 5);
row.TargetEORatio = row.TargetEOCount / max(row.WindowCount, 1);
row.DominantEO = uniqueEO(imax);
row.MeanRMSE = mean(T.WeightedVoltageRMSE, 'omitnan');
row.MedianRMSE = median(T.WeightedVoltageRMSE, 'omitnan');
row.BestRMSE = min(T.WeightedVoltageRMSE, [], 'omitnan');
row.WorstRMSE = max(T.WeightedVoltageRMSE, [], 'omitnan');
row.EOSequence = string(strjoin(string(eo(:).'), ','));
end

function label = time_label_local(tSec)
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

