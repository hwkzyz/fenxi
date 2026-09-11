%% Diagnose the voltage baseline used by the 20251222 static calibration library.
clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(thisDir);
rootDir = fileparts(projectDir);
libraryDir = fullfile(rootDir, 'static_gap_waveform_library_20260524_2000Hz');
calibrationFile = fullfile(projectDir, 'inputs', 'calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B5_S123.mat');
outDir = fullfile(thisDir, 'static_calibration_baseline_20260830');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

L = load(calibrationFile, 'CorrectedGapLibrary');
R = L.CorrectedGapLibrary.responseSurface;
B = R.staticBaselineByGap;

% Use the middle calibration condition as a representative raw record.
row = find(abs(B.trueReferenceGapMm - 1.3) < 1e-9, 1, 'first');
sourceName = char(B.sourceName(row));
binFile = fullfile(libraryDir, 'raw_binary_source', [sourceName, '.bin']);
headerFile = fullfile(libraryDir, 'raw_binary_source', [sourceName, '_header.txt']);
[t, yRawMv, fs] = read_acts1000_binary_channel_local(binFile, headerFile, 1);
baselineMv = B.baselineMv(row);

gapFolder = sprintf('recorded_%.2fmm_true_%.2fmm', ...
    B.recordedGapMm(row), B.trueReferenceGapMm(row));
gapFolder = strrep(gapFolder, '.', 'p');
W = readtable(fullfile(libraryDir, 'by_recorded_gap', gapFolder, 'waveforms.csv'), ...
    'TextType', 'string');
bladeId = R.referenceBladeId;
eventNo = W.eventNo(find(W.bladeId == bladeId, 1));
mask = W.bladeId == bladeId & W.eventNo == eventNo;
tPulse = W.timeRelative_s(mask) * 1000;
yPulseRawMv = W.voltageMv(mask);

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 15]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl); hold(ax1, 'on');
stride = max(1, floor(numel(t) / 30000));
plot(ax1, t(1:stride:end), yRawMv(1:stride:end), 'Color', [0.18 0.18 0.18], 'LineWidth', 0.55);
yline(ax1, baselineMv, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('Baseline %.1f mV', baselineMv));
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Raw voltage (mV)');
title(ax1, sprintf('(a) Continuous calibration record, true gap %.2f mm', B.trueReferenceGapMm(row)));
legend(ax1, 'Measured signal', sprintf('Baseline %.1f mV', baselineMv), 'Location', 'best');

ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, tPulse, yPulseRawMv, 'Color', [0.18 0.18 0.18], 'LineWidth', 0.9);
yline(ax2, baselineMv, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
plot(ax2, tPulse, yPulseRawMv - baselineMv, 'Color', [0 0.45 0.74], 'LineWidth', 1.0);
xlabel(ax2, 'Relative time (ms)'); ylabel(ax2, 'Voltage (mV)');
title(ax2, sprintf('(b) Extracted blade-%d waveform before and after baseline subtraction', bladeId));
legend(ax2, 'Raw extracted waveform', 'Estimated baseline', 'Baseline-subtracted waveform', ...
    'Location', 'best');

ax3 = nexttile(tl); hold(ax3, 'on');
colors = turbo(numel(R.gTrainMm));
for ig = 1:numel(R.gTrainMm)
    plot(ax3, R.xGrid, R.Yref(:, ig), 'Color', colors(ig, :), 'LineWidth', 0.65);
end
yline(ax3, 0, ':', 'Color', [0.35 0.35 0.35]);
xline(ax3, min(R.xGrid(R.effectiveWindow)), '--', 'Color', [0.5 0.5 0.5]);
xline(ax3, max(R.xGrid(R.effectiveWindow)), '--', 'Color', [0.5 0.5 0.5]);
xlabel(ax3, 'Spatial coordinate, x (mm)'); ylabel(ax3, 'Baseline-subtracted voltage (mV)');
title(ax3, '(c) Waveforms stored in the static response library');
cb = colorbar(ax3); cb.Label.String = 'True gap (mm)';
colormap(ax3, colors); caxis(ax3, [min(R.gTrainMm), max(R.gTrainMm)]);

axesList = [ax1, ax2, ax3];
set(axesList, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'Box', 'on', 'TickDir', 'in', 'LineWidth', 0.8);
set(findall(fig, 'Type', 'Legend'), 'FontName', 'Times New Roman', 'FontSize', 8);
set(findall(fig, 'Type', 'Text'), 'FontName', 'Times New Roman');

pngFile = fullfile(outDir, 'StaticCalibration_Baseline_RawAndCentered_20251222.png');
pdfFile = fullfile(outDir, 'StaticCalibration_Baseline_RawAndCentered_20251222.pdf');
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

summary = B;
writetable(summary, fullfile(outDir, 'StaticCalibration_Baseline_Summary_20251222.csv'));

fprintf('Representative raw file: %s\n', binFile);
fprintf('Raw range: %.3f to %.3f mV; estimated baseline: %.3f mV\n', ...
    min(yRawMv), max(yRawMv), baselineMv);
fprintf('Library Yref range: %.3f to %.3f mV\n', min(R.Yref(:)), max(R.Yref(:)));
fprintf('Saved: %s\n', pngFile);

function [t, y, fs] = read_acts1000_binary_channel_local(binPath, headerPath, selectedChannel)
headerText = fileread(headerPath);
chanCount = max(1, round(read_header_number_local(headerText, 'chanEnableCount', 1)));
fs = read_header_number_local(headerText, 'SampleRate', 2000);
resolution = round(read_header_number_local(headerText, 'resolution', 12));
rangeMax = read_header_number_local(headerText, sprintf('%drangeMaxValue', selectedChannel - 1), 5000);
rangeMin = read_header_number_local(headerText, sprintf('%drangeMinValue', selectedChannel - 1), -5000);
fid = fopen(binPath, 'rb');
cleanup = onCleanup(@() fclose(fid));
raw = fread(fid, inf, 'uint16=>uint16');
validLen = floor(numel(raw) / chanCount) * chanCount;
raw = reshape(raw(1:validLen), chanCount, []);
code = double(bitand(raw(selectedChannel, :), uint16(2^resolution - 1)));
y = ((rangeMax - rangeMin) / 2^resolution .* code + rangeMin).';
t = (0:numel(y)-1).' / fs;
end

function value = read_header_number_local(headerText, key, defaultValue)
expr = [regexptranslate('escape', key), '\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'];
tok = regexp(headerText, expr, 'tokens', 'once');
if isempty(tok), value = defaultValue; else, value = str2double(tok{1}); end
end
