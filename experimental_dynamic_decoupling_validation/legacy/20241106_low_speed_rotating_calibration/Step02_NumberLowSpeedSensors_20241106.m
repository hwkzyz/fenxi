%% Step02_NumberLowSpeedSensors_20241106
% Export low-speed numbering in an OPR-defined revolution view.
% Each valid revolution contributes one six-peak row in physical blade order.

clear; close all; clc;

%% Parameters to tune
analysisSensors = [2 3 5 7];
bladeCount = 6;
lowSpeedCase = '900';

viewEnable = true;
saveFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
sensorTag = ['S', sprintf('%d', analysisSensors)];

lowSpeedFingerprintFile = fullfile(routeDir, 'output', 'new_flow', '01_low_speed_reference', ...
    'LowSpeedReferenceFingerprint_20241106.mat');
outDir = fullfile(routeDir, 'output', 'new_flow', '02_low_speed_numbering', sensorTag);
outFile = fullfile(outDir, 'LowSpeedNumbering_20241106.mat');
summaryFile = fullfile(outDir, 'LowSpeedNumbering_Summary_20241106.csv');
figureDir = fullfile(routeDir, 'output', 'new_flow', 'figures', '02_low_speed_numbering');

if exist(lowSpeedFingerprintFile, 'file') ~= 2
    error(['Missing Step01 artifact:\n  %s\n\n' ...
        'Run Step01_BuildLowSpeedReferenceFingerprint_20241106 first.'], ...
        lowSpeedFingerprintFile);
end
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

loaded = load(lowSpeedFingerprintFile, 'LowSpeedReference');
if ~isfield(loaded, 'LowSpeedReference')
    error('File does not contain variable LowSpeedReference:\n  %s', lowSpeedFingerprintFile);
end
LowSpeedReference = loaded.LowSpeedReference;

validate_low_speed_reference_local(LowSpeedReference, analysisSensors, bladeCount, lowSpeedFingerprintFile);
[LowSpeedNumbering, LowSpeedSummary] = build_low_speed_numbering_tables_local( ...
    LowSpeedReference, analysisSensors, bladeCount);

save(outFile, 'LowSpeedNumbering', 'LowSpeedSummary');
writetable(LowSpeedNumbering, strrep(outFile, '.mat', '.csv'));
writetable(LowSpeedSummary, summaryFile);

fprintf('\n=== Step02: low-speed sensor numbering ===\n');
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Low-speed case: %s\n', lowSpeedCase);
fprintf('Source: %s\n', lowSpeedFingerprintFile);
fprintf('Saved table: %s\n', outFile);
fprintf('Saved summary: %s\n', summaryFile);
disp(LowSpeedSummary);

if viewEnable
    visualize_low_speed_numbering_local( ...
        LowSpeedReference, LowSpeedNumbering, LowSpeedSummary, bladeCount, saveFigures, figureDir);
end

function validate_low_speed_reference_local(LowSpeedReference, sensorIds, bladeCount, sourceFile)
requiredFields = {'fingerprints', 'revolution_ids', 'revolution_physical_peaks', ...
    'revolution_best_shifts', 'revolution_best_scores', 'revolution_score_margins', ...
    'sensor_dominant_shifts', 'sensor_physical_to_local'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(LowSpeedReference, name)
        error('LowSpeedReference is missing field %s:\n  %s', name, sourceFile);
    end
end

mustBeMap = {'fingerprints', 'revolution_ids', 'revolution_physical_peaks', ...
    'revolution_best_shifts', 'revolution_best_scores', 'revolution_score_margins', ...
    'sensor_dominant_shifts', 'sensor_physical_to_local'};
for i = 1:numel(mustBeMap)
    name = mustBeMap{i};
    if ~isa(LowSpeedReference.(name), 'containers.Map')
        error('LowSpeedReference.%s must be containers.Map.', name);
    end
end

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    requiredKeys = {'fingerprints', 'revolution_ids', 'revolution_physical_peaks', ...
        'revolution_best_shifts', 'revolution_best_scores', 'revolution_score_margins', ...
        'sensor_dominant_shifts', 'sensor_physical_to_local'};
    for j = 1:numel(requiredKeys)
        keyName = requiredKeys{j};
        if ~isKey(LowSpeedReference.(keyName), sid)
            error('LowSpeedReference.%s does not contain CH%d.', keyName, sid);
        end
    end
    fp = LowSpeedReference.fingerprints(sid);
    if numel(fp) < bladeCount
        error('Fingerprint for CH%d has only %d values, but bladeCount=%d.', ...
            sid, numel(fp), bladeCount);
    end
end
end

function [numberingTable, summaryTable] = build_low_speed_numbering_tables_local( ...
        LowSpeedReference, sensorIds, bladeCount)
numberingRows = struct([]);
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'ReferenceRevolutionID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ValidRevolutionCount', NaN, ...
    'MeanBestScore', NaN, ...
    'MinBestScore', NaN, ...
    'MeanScoreMargin', NaN, ...
    'FingerprintB1', NaN, ...
    'FingerprintB2', NaN, ...
    'FingerprintB3', NaN, ...
    'FingerprintB4', NaN, ...
    'FingerprintB5', NaN, ...
    'FingerprintB6', NaN), numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    revIds = force_row_vector_local(LowSpeedReference.revolution_ids(sid));
    peakMat = LowSpeedReference.revolution_physical_peaks(sid);
    bestShift = force_row_vector_local(LowSpeedReference.revolution_best_shifts(sid));
    bestScore = force_row_vector_local(LowSpeedReference.revolution_best_scores(sid));
    scoreMargin = force_row_vector_local(LowSpeedReference.revolution_score_margins(sid));
    physicalToLocal = force_row_vector_local(LowSpeedReference.sensor_physical_to_local(sid));
    dominantShift = LowSpeedReference.sensor_dominant_shifts(sid);
    fingerprint = force_row_vector_local(LowSpeedReference.fingerprints(sid));

    nRows = size(peakMat, 1);
    if numel(revIds) ~= nRows || numel(bestShift) ~= nRows || ...
            numel(bestScore) ~= nRows || numel(scoreMargin) ~= nRows
        error('CH%d low-speed revolution arrays are size-inconsistent.', sid);
    end

    for k = 1:nRows
        row = struct();
        row.SensorID = sid;
        row.RevolutionID = revIds(k);
        row.ReferenceRevolutionID = LowSpeedReference.reference_revolution_id;
        row.DominantShift = dominantShift;
        row.BestShift = bestShift(k);
        row.BestScore = bestScore(k);
        row.ScoreMargin = scoreMargin(k);
        row.LocalSlotOfB1 = physicalToLocal(1);
        for b = 1:bladeCount
            row.(sprintf('B%d', b)) = peakMat(k, b);
        end
        numberingRows = [numberingRows; row]; %#ok<AGROW>
    end

    summaryRows(i).SensorID = sid;
    summaryRows(i).ReferenceRevolutionID = LowSpeedReference.reference_revolution_id;
    summaryRows(i).DominantShift = dominantShift;
    summaryRows(i).LocalSlotOfB1 = physicalToLocal(1);
    summaryRows(i).ValidRevolutionCount = nRows;
    summaryRows(i).MeanBestScore = mean(bestScore, 'omitnan');
    summaryRows(i).MinBestScore = min(bestScore, [], 'omitnan');
    summaryRows(i).MeanScoreMargin = mean(scoreMargin, 'omitnan');
    for b = 1:bladeCount
        summaryRows(i).(sprintf('FingerprintB%d', b)) = fingerprint(b);
    end
end

if isempty(numberingRows)
    numberingTable = table();
else
    numberingTable = struct2table(numberingRows);
end
summaryTable = struct2table(summaryRows);
end

function visualize_low_speed_numbering_local( ...
        LowSpeedReference, LowSpeedNumbering, LowSpeedSummary, bladeCount, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

sensorIds = LowSpeedSummary.SensorID;
fpMat = nan(numel(sensorIds), bladeCount);
for b = 1:bladeCount
    fpMat(:, b) = LowSpeedSummary.(sprintf('FingerprintB%d', b));
end
rowMax = max(fpMat, [], 2, 'omitnan');
rowMax(~isfinite(rowMax) | rowMax == 0) = 1;
fpNorm = fpMat ./ rowMax;

fig = figure('Name', 'Step02 low-speed sensor numbering', 'Color', 'w', ...
    'Position', [100, 100, 1200, 760], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(1:bladeCount, sensorIds, fpNorm);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Mean six-peak fingerprint by physical blade');
colorbar;
box on;

nexttile;
yyaxis left;
bar(categorical(compose('CH%d', sensorIds)), LowSpeedSummary.LocalSlotOfB1);
ylabel('Local pulse slot assigned to B1');
yyaxis right;
plot(categorical(compose('CH%d', sensorIds)), LowSpeedSummary.ValidRevolutionCount, ...
    'o-k', 'LineWidth', 1.0);
ylabel('Valid OPR-defined revolutions');
title('B1 slot and valid revolution count');
grid on; box on;

nexttile;
scoreMat = nan(numel(sensorIds), bladeCount);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    revScores = LowSpeedReference.revolution_best_scores(sid);
    revMargins = LowSpeedReference.revolution_score_margins(sid);
    scoreMat(i, 1) = mean(revScores, 'omitnan');
    scoreMat(i, 2) = min(revScores, [], 'omitnan');
    scoreMat(i, 3) = mean(revMargins, 'omitnan');
    scoreMat(i, 4) = LowSpeedReference.sensor_dominant_shifts(sid);
end
bar(scoreMat);
xticklabels(compose('CH%d', sensorIds));
legend({'Mean score', 'Min score', 'Mean margin', 'Dominant shift'}, 'Location', 'best');
title('Matching quality summary');
grid on; box on;

nexttile;
axis off;
summaryText = evalc('disp(LowSpeedSummary)');
text(0, 1, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Consolas', 'FontSize', 8.5, 'Interpreter', 'none');
title(sprintf('Per-sensor summary and %d-row revolution table', height(LowSpeedNumbering)));

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step02_LowSpeedSensorNumbering_20241106.png'), ...
        'Resolution', 300);
end
end

function x = force_row_vector_local(x)
x = x(:).';
end
