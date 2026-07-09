%% Step07_AuditNumberingAndCorrespondence_20241106
% Audit low-speed numbering, high-speed numbering, and cross-sensor physical
% blade correspondence for the selected region.

clear; close all; clc;

P = NewFlow_Config_20241106();
require_file_local(P.files.lowSpeedNumbering, 'Step02 low-speed numbering');
require_file_local(P.files.highSpeedNumbering, 'Step04 high-speed numbering');
require_file_local(P.files.waveformLibrary, 'Step06 waveform library');
require_file_local(P.files.filteredWaveformLibrary, 'Step06 filtered waveform library');

loadedLow = load(P.files.lowSpeedNumbering, 'LowSpeedNumbering');
LowSpeedNumbering = loadedLow.LowSpeedNumbering;
loadedHigh = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
HighSpeedNumbering = loadedHigh.HighSpeedNumbering;
loadedWave = load(P.files.waveformLibrary, 'WaveformLibrary');
WaveformLibrary = loadedWave.WaveformLibrary;
loadedFiltered = load(P.files.filteredWaveformLibrary, 'FilteredWaveformLibrary');
FilteredWaveformLibrary = loadedFiltered.FilteredWaveformLibrary;

outFile = P.files.auditReport;
outDir = fileparts(outFile);
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

rows = [];
for bladeId = 1:P.machine.bladeCount
    for is = 1:numel(P.sensors.analysis)
        sid = P.sensors.analysis(is);
        high = HighSpeedNumbering.sensor([HighSpeedNumbering.sensor.sensor_id] == sid);
        wave = WaveformLibrary.Blade(bladeId).Sensor([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid);
        filteredCounts = filtered_point_counts_local(FilteredWaveformLibrary, bladeId, sid);
        rows = [rows; struct( ... %#ok<AGROW>
            'BladeID', bladeId, ...
            'SensorID', sid, ...
            'LowNumberingPresent', any(LowSpeedNumbering.SensorID == sid), ...
            'HighNumberingFiniteRows', all(isfinite(high.selected_rows_by_physical(:, bladeId))), ...
            'LapCount', numel(wave.Lap), ...
            'MinPointCount', min(arrayfun(@(x) numel(x.t), wave.Lap)), ...
            'MedianPointCount', median(arrayfun(@(x) numel(x.t), wave.Lap)), ...
            'MedianFilteredPointCount', median(filteredCounts, 'omitnan'), ...
            'MatchScore', high.best_score, ...
            'MatchMargin', high.score_margin)];
    end
end

Audit = struct2table(rows);
writetable(Audit, outFile);

fprintf('\n=== Step07: audit numbering and correspondence ===\n');
fprintf('Saved: %s\n', outFile);
disp(Audit);
if P.view.enable
    visualize_audit_local(Audit, P);
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function counts = filtered_point_counts_local(FilteredWaveformLibrary, bladeId, sid)
counts = nan(numel(FilteredWaveformLibrary.Blade(bladeId).Window), 1);
for w = 1:numel(FilteredWaveformLibrary.Blade(bladeId).Window)
    Bdl = FilteredWaveformLibrary.Blade(bladeId).Window(w).Bundle;
    if isempty(Bdl)
        continue;
    end
    counts(w) = nnz(Bdl.S == sid);
end
end

function visualize_audit_local(Audit, P)
figDir = fullfile(P.view.figureDir, '07_audit');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end
bladeIds = unique(Audit.BladeID).';
sensorIds = unique(Audit.SensorID).';
passMat = nan(numel(sensorIds), numel(bladeIds));
scoreMat = nan(numel(sensorIds), numel(bladeIds));
marginMat = nan(numel(sensorIds), numel(bladeIds));
for is = 1:numel(sensorIds)
    for ib = 1:numel(bladeIds)
        row = Audit(Audit.SensorID == sensorIds(is) & Audit.BladeID == bladeIds(ib), :);
        if height(row) == 1
            passMat(is, ib) = double(row.LowNumberingPresent && row.HighNumberingFiniteRows);
            scoreMat(is, ib) = row.MatchScore;
            marginMat(is, ib) = row.MatchMargin;
        end
    end
end

fig = figure('Name', 'Step07 numbering audit', 'Color', 'w');
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(bladeIds, sensorIds, passMat);
set(gca, 'YDir', 'normal');
clim([0 1]);
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Correspondence pass');
colorbar;
box on;

nexttile;
imagesc(bladeIds, sensorIds, scoreMat);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Match score');
colorbar;
box on;

nexttile;
imagesc(bladeIds, sensorIds, marginMat);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Score margin');
colorbar;
box on;

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step07_NumberingAudit_20241106.png'), 'Resolution', 300);
end
end
