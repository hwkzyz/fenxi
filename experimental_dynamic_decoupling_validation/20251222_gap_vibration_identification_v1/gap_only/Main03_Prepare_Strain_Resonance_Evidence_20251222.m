%% Step03: Locate vibration windows from strain spectrum and OPR/RPM
% This step reuses the verified legacy Step3Align logic.
% It should be run after Step02 and produces the strain-side frequency
% reference used by later decoupling validation.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20251222();
legacyDir = packageCfg.paths.preparation;
if exist(legacyDir, 'dir') ~= 7
    error('Legacy 20251222 helper folder not found: %s', legacyDir);
end
addpath(legacyDir);

cfg = Get_20251222_BTT_Config();
caseName = cfg.dynamic_cases{1};
showPlots = true;
forceRebuild = false;
strainPreprocessMethod = 'linear';

caseOutputDir = fullfile(cfg.output_root, caseName);
resultFile = fullfile(caseOutputDir, 'Step3_Spectrum_RPM_Alignment_20251222.mat');

needRebuild = forceRebuild;
if isfile(resultFile) && ~forceRebuild
    fprintf('Strain-OPR alignment result already exists. Loading:\n  %s\n', resultFile);
    load(resultFile, 'result');
    oldMethod = '';
    if isfield(result, 'strain_preprocess_method')
        oldMethod = result.strain_preprocess_method;
    end
    if ~strcmpi(string(oldMethod), string(strainPreprocessMethod))
        fprintf('Existing Step03 result used strain preprocess "%s"; rebuild as "%s".\n', ...
            string_or_empty_local(oldMethod), strainPreprocessMethod);
        needRebuild = true;
    end
else
    needRebuild = true;
end

if needRebuild
    fprintf('Running legacy strain-OPR alignment for %s...\n', caseName);
    result = Step3_Align_StrainSpectrum_With_RPM_20251222(caseName, ...
        cfg.step3_default_alignment_offset_sec, cfg.step3_align_order_candidates, ...
        cfg.step3_align_freq_plot_hz, showPlots);
end

outDir = packageCfg.paths.prepared;
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
save(fullfile(outDir, 'Step03_Strain_OPR_Vibration_Window_20251222.mat'), ...
    'result', '-v7.3');

%% Visualization: strain spectrum, tracked order and RPM
freqMask = result.stft_freq_hz >= cfg.step3_align_freq_plot_hz(1) & ...
    result.stft_freq_hz <= cfg.step3_align_freq_plot_hz(2);
ampDb = 20 * log10(result.stft_amp(freqMask, :) + eps);

figure('Name', '20251222 Step03 strain-OPR vibration location', ...
    'Color', 'w', 'Position', [80, 80, 1320, 900], 'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(result.stft_time_s, result.stft_freq_hz(freqMask), ampDb);
axis xy; grid on; box on;
colormap(gca, turbo);
colorbar;
hold on;
plot(result.stft_time_s, result.black_line_hz, 'w-', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('strain tracked %.0fX', result.best_order));
plot(result.shifted_rpm_time_s, result.shifted_btt_black_hz, 'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'OPR order reference');
if isfield(result, 'focus_info') && isfield(result.focus_info, 'time_window_sec')
    xline(result.focus_info.time_window_sec(1), 'c:', 'Focus start');
    xline(result.focus_info.time_window_sec(2), 'c:', 'Focus end');
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Strain STFT and OPR alignment: %s', caseName), 'Interpreter', 'none');
legend('Location', 'best');

nexttile;
plot(result.rpm_time_s, result.rpm_values, 'k-', 'LineWidth', 1.1);
grid on; box on;
ylabel('RPM');
title(sprintf('Rotor speed, best time shift = %.3f s', result.best_tau_sec));

nexttile;
bar(result.order_candidates, result.order_energy_scores, 0.75, 'FaceColor', [0.20, 0.45, 0.75]);
grid on; box on;
xlabel('Order candidate');
ylabel('Energy score');
title(sprintf('Order selection, best order = %.0f', result.best_order));

figFile = fullfile(outDir, 'Step03_Strain_OPR_Vibration_Window_20251222.png');
saveas(gcf, figFile);

fprintf('\nStep03 complete. Saved validation copy:\n  %s\n', ...
    fullfile(outDir, 'Step03_Strain_OPR_Vibration_Window_20251222.mat'));
fprintf('Saved visualization:\n  %s\n', figFile);

function txt = string_or_empty_local(x)
if isempty(x)
    txt = '(empty)';
else
    txt = char(string(x));
end
end

