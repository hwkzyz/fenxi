%% Step02: Calculate full-time BTT displacement and strain synchronization
% This step reuses the verified legacy Step3 displacement calculation.
% It should be run once after Step01. Later decoupling steps should read
% its saved result instead of recalculating the full-time displacement.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20250527();
legacyDir = packageCfg.paths.preparation;
addpath(legacyDir);

cfg = Get_20250527_BTT_Config();
caseName = packageCfg.case.dynamicCase;
showPlots = packageCfg.run.showPlots;
forceRebuild = packageCfg.run.forceRebuild;

caseOutputDir = fullfile(cfg.output_root, caseName);
resultFile = fullfile(caseOutputDir, 'Step3_Result_20250527.mat');

if isfile(resultFile) && ~forceRebuild
    fprintf('Full-time displacement result already exists. Loading:\n  %s\n', resultFile);
    load(resultFile, 'result');
else
    fprintf('Running legacy full-time displacement calculation for %s...\n', caseName);
    result = Step3_Displacement_Calculation_20250527(caseName, ...
        cfg.step3_default_alignment_offset_sec, showPlots);
end

outDir = packageCfg.paths.results;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
save(fullfile(outDir, 'Step02_FullTime_BTT_Displacement_20250527.mat'), ...
    'result', '-v7.3');

%% Visualization: BTT displacement, strain waveform and RPM
sensorIds = result.loaded_sensor_ids;
if isempty(sensorIds)
    sensorIds = cfg.step3_default_btt_sensor_ids;
end

figure('Name', '20250527 Step02 full-time BTT displacement', ...
    'Color', 'w', 'Position', [80, 80, 1320, 860], 'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; grid on; box on;
for sid = sensorIds
    vibFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d_vib_final.mat', sid));
    if ~isfile(vibFile)
        continue;
    end
    S = load(vibFile, 'jilublade');
    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 6
        continue;
    end
    tBtt = S.jilublade(:, 3);
    yBtt = S.jilublade(:, 6);
    plot(tBtt, yBtt, '.', 'MarkerSize', 3, ...
        'DisplayName', sprintf('CH%d', sid));
end
if isfield(result, 'btt_full_xlim') && numel(result.btt_full_xlim) == 2
    xlim(result.btt_full_xlim);
end
ylabel('BTT disp. (mm)');
title(sprintf('Full-time BTT displacement: %s, strain shift = %.3f s', ...
    caseName, result.btt_time_shift_sec), 'Interpreter', 'none');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
if isfield(result, 'strain_file') && isfile(result.strain_file)
    load(result.strain_file, 'Datas');
    tStrain = Datas(:, 1) + result.strain_time_offset_total_sec;
    vStrain = detrend(Datas(:, 2));
    stepPlot = max(1, floor(numel(tStrain) / 80000));
    idxPlot = 1:stepPlot:numel(tStrain);
    plot(tStrain(idxPlot), vStrain(idxPlot), 'Color', [0.15, 0.35, 0.80], 'LineWidth', 0.7);
    clear Datas;
end
if isfield(result, 'sync_detail_xlim') && numel(result.sync_detail_xlim) == 2
    xline(result.sync_detail_xlim(1), '--', 'Sync detail start');
    xline(result.sync_detail_xlim(2), '--', 'Sync detail end');
end
if isfield(result, 'fft_window') && numel(result.fft_window) == 2
    xline(result.fft_window(1), ':', 'FFT start');
    xline(result.fft_window(2), ':', 'FFT end');
end
ylabel('Strain signal');
title('Aligned strain waveform');

nexttile;
plot(result.rpm_time_s, result.rpm_values, 'k-', 'LineWidth', 1.1);
grid on; box on;
xlabel('Time (s)');
ylabel('RPM');
title('OPR-derived rotor speed');

figFile = fullfile(outDir, 'Step02_FullTime_BTT_Displacement_20250527.png');
saveas(gcf, figFile);

fprintf('\nStep02 complete. Saved validation copy:\n  %s\n', ...
    fullfile(outDir, 'Step02_FullTime_BTT_Displacement_20250527.mat'));
fprintf('Saved visualization:\n  %s\n', figFile);
