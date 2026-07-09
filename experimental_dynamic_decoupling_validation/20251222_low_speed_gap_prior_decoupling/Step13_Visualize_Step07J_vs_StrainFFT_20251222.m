%% Visualize Step07J identification results against local strain FFT evidence
% This script reads the current phase-safe Step07J result and the Step04
% local strain FFT table.  The comparison is frequency-based: the strain
% FFT amplitude is still in voltage units, so it is not used as a direct
% displacement-amplitude reference.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_strain_fft_comparison');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

datasetLabel = sprintf('%s %s', C0.dataset, C0.caseTag);
step07jFile = resolve_step07j_result_local(outDir, C0);
step04File = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.mat');
step04FftCsv = fullfile(outDir, 'Step04_LocalStrainFFT_20251222.csv');

if ~isfile(step04File) || ~isfile(step04FftCsv)
    error('Run Step04 first. Missing Step04 local FFT outputs.');
end

S7 = load(step07jFile, 'Result');
R7 = S7.Result;
T7 = R7.Trend;
W7 = R7.WindowResult;
fprintf('Using Step07J result: %s\n', step07jFile);

S4 = load(step04File, 'regionTable');
regionTable = S4.regionTable;
localFftTable = readtable(step04FftCsv, 'TextType', 'string');

timeStart = arrayfun(@(s) s.timeWindow(1), W7(:));
timeEnd = arrayfun(@(s) s.timeWindow(2), W7(:));
timeMid = 0.5 * (timeStart + timeEnd);
windowEnvelope = [min(timeStart), max(timeEnd)];

[methodNames, methodLabels, eoMat, freqMat, rmseMat] = collect_available_methods_local(T7);
if isempty(methodNames)
    error('The selected Step07J Trend table has no recognizable method columns.');
end
mainMethodIdx = find(methodNames == "gap_tilt", 1, 'first');
if isempty(mainMethodIdx)
    mainMethodIdx = 1;
end
mainEoValues = eoMat(:, mainMethodIdx);
if ~any(isfinite(mainEoValues))
    mainEoValues = eoMat(isfinite(eoMat));
end
if isempty(mainEoValues)
    error('No finite EO values found in Step07J Trend.');
end
mainEo = mode(round(mainEoValues(isfinite(mainEoValues))));

selectedRegionIdx = select_strain_fft_reference_region_local( ...
    regionTable, localFftTable, windowEnvelope, mainEo);
strainRefFreqHz = localFftTable.fftPeakFreqHz(selectedRegionIdx);
strainSyncFreqHz = localFftTable.syncFreqAtRegionOrderHz(selectedRegionIdx);
strainRegionStart = regionTable.regionStart(selectedRegionIdx);
strainRegionEnd = regionTable.regionEnd(selectedRegionIdx);
strainRegionOrder = regionTable.regionPeakOrder(selectedRegionIdx);

freqErrMat = freqMat - strainRefFreqHz;

comparisonTable = table(T7.window_id, timeStart(:), timeEnd(:), timeMid(:), ...
    repmat(selectedRegionIdx, height(T7), 1), ...
    repmat(strainRegionOrder, height(T7), 1), ...
    repmat(strainRefFreqHz, height(T7), 1), ...
    repmat(strainSyncFreqHz, height(T7), 1), ...
    'VariableNames', {'window_id', 'time_start_s', 'time_end_s', 'time_mid_s', ...
    'strain_region_index', 'strain_region_order', 'strain_fft_peak_hz', ...
    'strain_sync_freq_hz'});
for m = 1:numel(methodNames)
    prefix = char(methodNames(m));
    comparisonTable.([prefix, '_EO']) = eoMat(:, m);
    comparisonTable.([prefix, '_freq_hz']) = freqMat(:, m);
    comparisonTable.([prefix, '_freq_error_hz']) = freqErrMat(:, m);
    comparisonTable.([prefix, '_rmse_mV']) = rmseMat(:, m);
end

csvFile = fullfile(outDir, ...
    sprintf('Step07J_vs_StrainFFT_FrequencyComparison_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(comparisonTable, csvFile);

fig = figure('Name', '20251222 Step07J vs strain FFT', 'Color', 'w', ...
    'Position', [80, 80, 1320, 880], 'NumberTitle', 'off');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; grid on; box on;
for i = 1:height(regionTable)
    if regionTable.regionPeakOrder(i) == mainEo
        faceColor = [1.00, 0.86, 0.35];
        faceAlpha = 0.22;
    else
        faceColor = [0.70, 0.70, 0.70];
        faceAlpha = 0.10;
    end
    patch([regionTable.regionStart(i), regionTable.regionEnd(i), ...
        regionTable.regionEnd(i), regionTable.regionStart(i)], ...
        [0, 0, 1, 1], faceColor, 'FaceAlpha', faceAlpha, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for i = 1:numel(timeStart)
    plot([timeStart(i), timeEnd(i)], [0.55, 0.55], '-', ...
        'Color', [0.10, 0.35, 0.75], 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
end
plot(timeMid, 0.55 * ones(size(timeMid)), 'o', 'Color', [0.10, 0.35, 0.75], ...
    'MarkerFaceColor', [0.10, 0.35, 0.75], 'DisplayName', 'Step07J windows');
xline(strainRegionStart, 'k:', 'DisplayName', 'selected strain FFT region');
xline(strainRegionEnd, 'k:', 'HandleVisibility', 'off');
ylim([0, 1]);
yticks([]);
xlabel('BTT-aligned time (s)');
title(sprintf('%s: Step07J windows inside Step04 strain region %d (%dX)', ...
    datasetLabel, selectedRegionIdx, strainRegionOrder));
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
colors = method_colors_local(numel(methodNames));
for m = 1:numel(methodNames)
    plot(T7.window_id, freqMat(:, m), '-o', 'LineWidth', 1.2, ...
        'Color', colors(m, :), 'MarkerSize', 4, 'DisplayName', methodLabels(m));
end
yline(strainRefFreqHz, 'r--', sprintf('strain local FFT %.3f Hz', strainRefFreqHz), ...
    'LineWidth', 1.2, 'DisplayName', 'strain FFT peak');
yline(strainSyncFreqHz, 'k:', sprintf('%dX mean RPM %.3f Hz', ...
    strainRegionOrder, strainSyncFreqHz), 'LineWidth', 1.0, ...
    'DisplayName', 'strain region order line');
xlabel('Step07J window');
ylabel('Frequency (Hz)');
title('Identified frequency compared with local strain FFT');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
for m = 1:numel(methodNames)
    plot(T7.window_id, freqErrMat(:, m), '-o', 'LineWidth', 1.2, ...
        'Color', colors(m, :), 'MarkerSize', 4, 'DisplayName', methodLabels(m));
end
yline(0, 'k-', 'HandleVisibility', 'off');
xlabel('Step07J window');
ylabel('Frequency error (Hz)');
title('Frequency error relative to strain local FFT peak');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
for m = 1:numel(methodNames)
    plot(T7.window_id, rmseMat(:, m), '-o', 'LineWidth', 1.2, ...
        'Color', colors(m, :), 'MarkerSize', 4, 'DisplayName', methodLabels(m));
end
xlabel('Step07J window');
ylabel('Waveform RMSE (mV)');
title('Waveform objective from Step07J');
legend('Location', 'best');

figFile = fullfile(figDir, ...
    sprintf('Step07J_vs_StrainFFT_FrequencyComparison_%s_%s.png', C0.dataset, C0.caseTag));
saveas(fig, figFile);

fprintf('Step07J vs strain FFT comparison complete.\n');
fprintf('Selected Step04 region %d: %.3f-%.3f s, %dX, FFT %.6f Hz, order line %.6f Hz.\n', ...
    selectedRegionIdx, strainRegionStart, strainRegionEnd, ...
    strainRegionOrder, strainRefFreqHz, strainSyncFreqHz);
summaryVars = {'window_id'};
for m = 1:numel(methodNames)
    prefix = char(methodNames(m));
    summaryVars{end+1} = [prefix, '_freq_error_hz']; %#ok<SAGROW>
    summaryVars{end+1} = [prefix, '_rmse_mV']; %#ok<SAGROW>
end
disp(comparisonTable(:, summaryVars));
fprintf('Saved:\n  %s\n  %s\n', csvFile, figFile);

function step07jFile = resolve_step07j_result_local(outDir, C0)
regionTag = '';
if isfield(C0, 'flowConfig') && isfield(C0.flowConfig, 'identification') && ...
        isfield(C0.flowConfig.identification, 'resonanceRegionShortTag')
    regionTag = C0.flowConfig.identification.resonanceRegionShortTag;
end

candidates = {
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', C0.dataset, C0.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_%s_main_gaptilt.mat', C0.dataset, C0.bladeId, C0.sensorTag, regionTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaptilt.mat', C0.dataset, C0.bladeId, C0.sensorTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_method_compare.mat', C0.dataset, C0.bladeId, C0.sensorTag)
    };
candidates = candidates(~cellfun(@(s) contains(s, '__'), candidates));
for i = 1:numel(candidates)
    candidate = fullfile(outDir, candidates{i});
    if isfile(candidate)
        step07jFile = candidate;
        return;
    end
end

files = dir(fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s*.mat', C0.dataset, C0.caseTag)));
if isempty(files)
    files = dir(fullfile(outDir, sprintf( ...
        'Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s*.mat', ...
        C0.dataset, C0.bladeId, C0.sensorTag)));
end
if ~isempty(files)
    [~, idx] = max([files.datenum]);
    step07jFile = fullfile(files(idx).folder, files(idx).name);
    warning('Using latest matching Step07J result because canonical files were not found:\n  %s', step07jFile);
    return;
end

error(['Missing Step07J result for %s. Expected a main_gaptilt or method_compare file under:\n' ...
    '  %s\nRun Step07J/Step00 first.'], C0.caseTag, outDir);
end

function [methodNames, methodLabels, eoMat, freqMat, rmseMat] = collect_available_methods_local(T)
spec = {
    "direct", "direct/fixed", "direct_EO", "direct_frequency_hz", "direct_rmse_mV"
    "fixed", "fixed", "fixed_EO", "fixed_frequency_hz", "fixed_rmse_mV"
    "gap_only", "gap only", "gap_EO", "gap_frequency_hz", "gap_rmse_mV"
    "gap_tilt", "gap + tilt", "tilt_EO", "tilt_frequency_hz", "tilt_rmse_mV"
    "gap_tilt_shift", "gap + tilt + shift", "shift_EO", "shift_frequency_hz", "shift_rmse_mV"
    };
vars = string(T.Properties.VariableNames);
methodNames = strings(1, 0);
methodLabels = strings(1, 0);
eoMat = [];
freqMat = [];
rmseMat = [];
for i = 1:size(spec, 1)
    if all(ismember(string(spec(i, 3:5)), vars))
        methodNames(end+1) = spec{i, 1}; %#ok<AGROW>
        methodLabels(end+1) = spec{i, 2}; %#ok<AGROW>
        eoMat(:, end+1) = T.(spec{i, 3}); %#ok<AGROW>
        freqMat(:, end+1) = T.(spec{i, 4}); %#ok<AGROW>
        rmseMat(:, end+1) = T.(spec{i, 5}); %#ok<AGROW>
    end
end
end

function colors = method_colors_local(n)
base = [
    0.35, 0.35, 0.35
    0.10, 0.45, 0.72
    0.10, 0.55, 0.25
    0.70, 0.30, 0.10
    0.45, 0.25, 0.65
    ];
if n <= size(base, 1)
    colors = base(1:n, :);
else
    colors = lines(n);
end
end

function selectedIdx = select_strain_fft_reference_region_local( ...
    regionTable, localFftTable, windowEnvelope, mainEo)
n = height(regionTable);
score = -inf(n, 1);
for i = 1:n
    overlapStart = max(regionTable.regionStart(i), windowEnvelope(1));
    overlapEnd = min(regionTable.regionEnd(i), windowEnvelope(2));
    overlap = max(0, overlapEnd - overlapStart);
    orderMatch = double(regionTable.regionPeakOrder(i) == mainEo);
    fftOrderMatch = double(abs(localFftTable.fftPeakOrderApprox(i) - mainEo) <= 0.25);
    score(i) = 100 * orderMatch + 50 * fftOrderMatch + overlap;
end
[~, selectedIdx] = max(score);
end
