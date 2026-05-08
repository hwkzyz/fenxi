%% Step 6B: SNR sweep using full-time-domain waveform SNR
% SNR definition:
%   SNR_all = 20*log10(rms(V_clean - baseline) / rms(noise))
% This script converts target SNR_all values to cfg.noiseRatio, then runs
% the low-template-guided blind high-speed identification pipeline.
clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

targetSNRdB = [50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
numRepeat = 3;
baseSeed = 20260430;
overridePath = fullfile(outDir, '__blind_template_config_override.mat');

oldFigureVisible = get(0, 'DefaultFigureVisible');
cleanupObj = onCleanup(@() restore_sweep_state(oldFigureVisible, overridePath));
set(0, 'DefaultFigureVisible', 'off');

% Build a noise-free reference once to calibrate full-waveform SNR.
overrideCfg = struct();
overrideCfg.makeFigures = false;
overrideCfg.noiseRatio = 0;
save(overridePath, 'overrideCfg');
rng(baseSeed, 'twister');
run(fullfile(scriptDir, 'Step_0_Config_BlindTemplate_HighSpeed.m'));
run(fullfile(scriptDir, 'Step_1_Simulate_BlindTemplate_Data.m'));
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'Data_High', 'truth');

signalAll = Data_High.V_cap_clean(:) - Data_High.baseline;
signalRmsAll = sqrt(mean(signalAll.^2));
templateProbe = interp1(truth.templateX, truth.templateY, ...
    linspace(truth.domain(1), truth.domain(2), 400)', 'pchip');
templateRange = range(templateProbe);
noiseRatioForTarget = signalRmsAll ./ (templateRange .* 10.^(targetSNRdB(:) / 20));

numCases = numel(targetSNRdB) * numRepeat;
records = repmat(make_empty_record(), numCases, 1);
idxCase = 0;

for iSnr = 1:numel(targetSNRdB)
    for iRep = 1:numRepeat
        idxCase = idxCase + 1;
        overrideCfg = struct();
        overrideCfg.makeFigures = false;
        overrideCfg.noiseRatio = noiseRatioForTarget(iSnr);
        save(overridePath, 'overrideCfg');

        seed = baseSeed + 1000 * iSnr + iRep;
        rng(seed, 'twister');
        caseTimer = tic;
        run(fullfile(scriptDir, 'Step_0_Config_BlindTemplate_HighSpeed.m'));
        run(fullfile(scriptDir, 'Step_1_Simulate_BlindTemplate_Data.m'));
        run(fullfile(scriptDir, 'Step_2B_Build_HighSpeed_PointCloud.m'));
        run(fullfile(scriptDir, 'Step_3B_Waveform_Profile_Identify.m'));
        elapsed = toc(caseTimer);

        load(fullfile(outDir, 'stage1_simulated_data.mat'), 'Data_High');
        load(fullfile(outDir, 'stage2B_point_cloud.mat'), 'cloud');
        load(fullfile(outDir, 'stage3B_waveform_profile_result.mat'), 'summary');

        noise = Data_High.V_cap(:) - Data_High.V_cap_clean(:);
        measuredSnr = 20 * log10(signalRmsAll / sqrt(mean(noise.^2)));

        records(idxCase).target_snr_db = targetSNRdB(iSnr);
        records(idxCase).measured_snr_db = measuredSnr;
        records(idxCase).noise_ratio = noiseRatioForTarget(iSnr);
        records(idxCase).repeat_id = iRep;
        records(idxCase).seed = seed;
        records(idxCase).mean_freq_error_Hz = summary.mean_freq_error_Hz;
        records(idxCase).mean_amp_error_mm = summary.mean_amp_error_mm;
        records(idxCase).profile_rmse = summary.profile_rmse;
        records(idxCase).profile_elapsed_s = summary.profile_elapsed_s;
        records(idxCase).case_elapsed_s = elapsed;
        records(idxCase).num_point_samples = numel(cloud.V);

        fprintf('[Step 6B] target SNR %.1f dB, repeat %d/%d: freq err %.4g Hz, amp err %.4g mm, RMSE %.4g, elapsed %.3f s\n', ...
            targetSNRdB(iSnr), iRep, numRepeat, summary.mean_freq_error_Hz, ...
            summary.mean_amp_error_mm, summary.profile_rmse, elapsed);
        close all;
    end
end

resultTable = struct2table(records);
summaryTable = summarize_by_snr(resultTable);
save(fullfile(outDir, 'stage6B_snr_sweep_fullwaveform.mat'), ...
    'resultTable', 'summaryTable', 'targetSNRdB', 'numRepeat', ...
    'signalRmsAll', 'templateRange', 'noiseRatioForTarget');
writetable(resultTable, fullfile(outDir, 'stage6B_snr_sweep_fullwaveform_trials.csv'));
writetable(summaryTable, fullfile(outDir, 'stage6B_snr_sweep_fullwaveform_summary.csv'));

set(0, 'DefaultFigureVisible', oldFigureVisible);
figure('Name', 'Step 6B - Full-Waveform SNR Sweep', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 14]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
errorbar(summaryTable.target_snr_db, summaryTable.mean_freq_error_Hz_mean, ...
    summaryTable.mean_freq_error_Hz_std, 'o-', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Full-waveform SNR (dB)');
ylabel('Mean frequency error (Hz)');
title('Frequency identification');
grid on;

nexttile;
errorbar(summaryTable.target_snr_db, summaryTable.mean_amp_error_mm_mean, ...
    summaryTable.mean_amp_error_mm_std, 'o-', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Full-waveform SNR (dB)');
ylabel('Mean amplitude error (mm)');
title('Amplitude identification');
grid on;

nexttile;
errorbar(summaryTable.target_snr_db, summaryTable.profile_rmse_mean, ...
    summaryTable.profile_rmse_std, 'o-', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Full-waveform SNR (dB)');
ylabel('Profile RMSE');
title('Template-profile fitting');
grid on;

nexttile;
errorbar(summaryTable.target_snr_db, summaryTable.profile_elapsed_s_mean, ...
    summaryTable.profile_elapsed_s_std, 'o-', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Full-waveform SNR (dB)');
ylabel('Step 3B elapsed (s)');
title('Runtime');
grid on;

exportgraphics(gcf, fullfile(outDir, 'stage6B_snr_sweep_fullwaveform.png'), 'Resolution', 300);

fprintf('\n[Step 6B] Full-waveform SNR sweep summary:\n');
disp(summaryTable);
fprintf('[Step 6B] Saved results to: %s\n', outDir);

function record = make_empty_record()
record = struct();
record.target_snr_db = NaN;
record.measured_snr_db = NaN;
record.noise_ratio = NaN;
record.repeat_id = NaN;
record.seed = NaN;
record.mean_freq_error_Hz = NaN;
record.mean_amp_error_mm = NaN;
record.profile_rmse = NaN;
record.profile_elapsed_s = NaN;
record.case_elapsed_s = NaN;
record.num_point_samples = NaN;
end

function summaryTable = summarize_by_snr(resultTable)
snrList = unique(resultTable.target_snr_db, 'stable');
summary = repmat(struct( ...
    'target_snr_db', NaN, ...
    'measured_snr_db_mean', NaN, ...
    'noise_ratio', NaN, ...
    'mean_freq_error_Hz_mean', NaN, ...
    'mean_freq_error_Hz_std', NaN, ...
    'mean_amp_error_mm_mean', NaN, ...
    'mean_amp_error_mm_std', NaN, ...
    'profile_rmse_mean', NaN, ...
    'profile_rmse_std', NaN, ...
    'profile_elapsed_s_mean', NaN, ...
    'profile_elapsed_s_std', NaN, ...
    'num_trials', NaN), numel(snrList), 1);
for i = 1:numel(snrList)
    idx = resultTable.target_snr_db == snrList(i);
    summary(i).target_snr_db = snrList(i);
    summary(i).measured_snr_db_mean = mean(resultTable.measured_snr_db(idx));
    summary(i).noise_ratio = mean(resultTable.noise_ratio(idx));
    summary(i).mean_freq_error_Hz_mean = mean(resultTable.mean_freq_error_Hz(idx));
    summary(i).mean_freq_error_Hz_std = std(resultTable.mean_freq_error_Hz(idx));
    summary(i).mean_amp_error_mm_mean = mean(resultTable.mean_amp_error_mm(idx));
    summary(i).mean_amp_error_mm_std = std(resultTable.mean_amp_error_mm(idx));
    summary(i).profile_rmse_mean = mean(resultTable.profile_rmse(idx));
    summary(i).profile_rmse_std = std(resultTable.profile_rmse(idx));
    summary(i).profile_elapsed_s_mean = mean(resultTable.profile_elapsed_s(idx));
    summary(i).profile_elapsed_s_std = std(resultTable.profile_elapsed_s(idx));
    summary(i).num_trials = nnz(idx);
end
summaryTable = struct2table(summary);
end

function restore_sweep_state(oldFigureVisible, overridePath)
set(0, 'DefaultFigureVisible', oldFigureVisible);
if isfile(overridePath)
    delete(overridePath);
end
end
