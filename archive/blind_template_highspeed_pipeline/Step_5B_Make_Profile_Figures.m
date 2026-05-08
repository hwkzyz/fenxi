%% Step 5B: visualize waveform-level blind profile identification
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage2B_point_cloud.mat'), 'cloud');
load(fullfile(outDir, 'stage3B_waveform_profile_result.mat'), 'profileResult', 'summary');

figure('Name', 'Blind-Template Step 5B - Waveform Profile Summary', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 25, 14]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(truth.templateX, truth.templateY, 'k-', 'LineWidth', 1.1, 'DisplayName', 'truth hidden template');
plot(profileResult.xGrid, profileResult.T, 'b-', 'LineWidth', 1.2, 'DisplayName', 'blind profile template');
xlabel('Corrected coordinate \xi (mm)');
ylabel('Voltage');
title('Shared template');
legend('Location', 'best', 'Box', 'off');

nexttile;
plot(profileResult.costHist, 'o-', 'LineWidth', 1.2);
xlabel('Outer iteration');
ylabel('Profile RMSE');
title('Convergence');

nexttile; hold on;
idx = 1:min(3000, numel(cloud.V));
scatter(cloud.x(idx), cloud.V(idx), 5, [0.70 0.70 0.70], 'filled');
scatter(profileResult.xi(idx), cloud.V(idx), 5, profileResult.residual(idx), 'filled');
xlabel('Coordinate (mm)');
ylabel('Voltage');
title('Nominal and dynamically corrected cloud');
colorbar;

nexttile;
bar(categorical({'f1', 'f2'}), [summary.f_est(:), summary.f_true(:)]);
ylabel('Frequency (Hz)');
title('Frequency identification');
legend({'estimated', 'truth'}, 'Location', 'best', 'Box', 'off');

nexttile;
bar(categorical({'A1', 'A2'}), [summary.A_est(:), summary.A_true(:)]);
ylabel('Amplitude (mm)');
title('Amplitude identification');
legend({'estimated', 'truth'}, 'Location', 'best', 'Box', 'off');

nexttile;
text(0.02, 0.82, sprintf('Mean frequency error: %.4f Hz', summary.mean_freq_error_Hz), 'FontSize', 11);
text(0.02, 0.62, sprintf('Mean amplitude error: %.4f mm', summary.mean_amp_error_mm), 'FontSize', 11);
text(0.02, 0.42, sprintf('Profile RMSE: %.4e', summary.profile_rmse), 'FontSize', 11);
text(0.02, 0.22, sprintf('Samples used: %d', numel(cloud.V)), 'FontSize', 11);
axis off;
title('Summary metrics');

fprintf('[Step 5B] Waveform-level blind profile final summary:\n');
disp(summary);
