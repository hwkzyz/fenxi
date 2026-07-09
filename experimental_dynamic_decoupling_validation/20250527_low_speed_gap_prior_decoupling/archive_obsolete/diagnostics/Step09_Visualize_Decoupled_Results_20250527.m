%% Step09: Visualize main clearance-vibration decoupling results
% Loads the Step07 main result and exports a compact sliding-window figure.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step09');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

targetBlade = 1;
analysisSensors = [1, 3, 6];
sensorOverride = strtrim(getenv('STEP09_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP09_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', analysisSensors)];

resultFile = fullfile(outDir, sprintf( ...
    'Step07_Main_Decoupled_Identification_20250527_B%d_%s.mat', ...
    targetBlade, sensorTag));
if ~isfile(resultFile)
    error('Run Step07 first. Missing file: %s', resultFile);
end

S = load(resultFile, 'result', 'summaryTable', 'sensorTable');
result = S.result;
trendTable = result.Trend;
sensorTable = S.sensorTable;
hasGapPrior = isfield(result, 'gapPrior') && isfield(result.gapPrior, 'sensorPrior');

fig = figure('Name', '20250527 Step09 decoupled result visualization', ...
    'Color', 'w', 'Position', [80, 80, 1450, 900], 'NumberTitle', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot(trendTable.window_id, trendTable.frequency_hz, '-o', ...
    'LineWidth', 1.3, 'MarkerSize', 5);
plot(trendTable.window_id, trendTable.rot_freq_hz * result.referenceOrder, '--', ...
    'LineWidth', 1.1, 'Color', [0.35, 0.35, 0.35]);
xlabel('Sliding window');
ylabel('Frequency (Hz)');
title('Identified frequency');
legend({'identified', sprintf('EO%d local reference', result.referenceOrder)}, 'Location', 'best');

nexttile; hold on; grid on; box on;
plot(trendTable.window_id, trendTable.amplitude_mm, '-o', ...
    'LineWidth', 1.3, 'MarkerSize', 5);
xlabel('Sliding window');
ylabel('Amplitude (mm)');
title('Identified vibration amplitude');

nexttile; hold on; grid on; box on;
plot(trendTable.window_id, trendTable.mean_gap_mm, '-o', ...
    'LineWidth', 1.3, 'MarkerSize', 5);
xlabel('Sliding window');
ylabel('Mean effective gap (mm)');
title('Mean effective clearance');

nexttile; hold on; grid on; box on;
xNum = 1:numel(sensorTable.sensorId);
xLabels = "CH" + string(sensorTable.sensorId);
bar(xNum, sensorTable.gHatMm, 'DisplayName', 'dynamic estimate');
if hasGapPrior
    [gCenter, gLower, gUpper] = get_gap_prior_vectors_local(result.gapPrior, sensorTable.sensorId(:).');
    errorbar(xNum, gCenter, gCenter - gLower, gUpper - gCenter, ...
        'k.', 'LineWidth', 1.2, 'CapSize', 12, 'DisplayName', 'low-speed prior');
    plot(xNum, gCenter, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'DisplayName', 'prior center');
end
set(gca, 'XTick', xNum, 'XTickLabel', cellstr(xLabels));
xlabel('Sensor');
ylabel('Effective gap (mm)');
title(sprintf('Best window %d sensor gaps', result.BestWindowIndex));
legend('Location', 'best');

figFile = fullfile(figDir, sprintf( ...
    'Step09_Decoupled_Result_20250527_B%d_%s.png', targetBlade, sensorTag));
saveas(fig, figFile);

fprintf('\nStep09 complete.\n');
fprintf('Loaded result:\n  %s\n', resultFile);
fprintf('Figure saved to:\n  %s\n', figFile);
disp(S.summaryTable);
disp(sensorTable);

function [gCenter, gLower, gUpper] = get_gap_prior_vectors_local(GapPrior, sensorIds)
gCenter = nan(numel(sensorIds), 1);
gLower = nan(numel(sensorIds), 1);
gUpper = nan(numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    idx = find([GapPrior.sensorPrior.sensorId] == sid, 1, 'first');
    if isempty(idx)
        continue;
    end
    gCenter(i) = GapPrior.sensorPrior(idx).gCenterMm;
    gLower(i) = GapPrior.sensorPrior(idx).gLowerMm;
    gUpper(i) = GapPrior.sensorPrior(idx).gUpperMm;
end
end
