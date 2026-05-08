%% Step 5: visualize blind-template reconstruction and vibration fit
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage2_segments.mat'), 'segments');
load(fullfile(outDir, 'stage3_blind_template.mat'), 'blindResult');
load(fullfile(outDir, 'stage4_vibration_fit.mat'), 'vibFit', 'summary');

tVec = [segments.tCenter]';
uTruth = zeros(size(tVec));
for i = 1:numel(truth.A)
    uTruth = uTruth + truth.A(i) * sin(2*pi*truth.f(i)*tVec + truth.phi(i));
end

figure('Name', 'Blind-Template Step 5 - Summary', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 24, 13]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(truth.templateX, truth.templateY, 'k-', 'LineWidth', 1.2, 'DisplayName', 'truth hidden template');
plot(blindResult.xRef, blindResult.T, 'b-', 'LineWidth', 1.2, 'DisplayName', 'blind reconstructed template');
xlabel('Position x (mm)');
ylabel('Voltage / capacitance');
title('Template comparison');
legend('Location', 'best', 'Box', 'off');

nexttile;
plot(blindResult.costHist, 'o-', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Mean registration RMSE');
title('Blind-template convergence');

nexttile; hold on;
scatter(tVec * 1000, blindResult.d, 12, 'b', 'filled', 'DisplayName', 'recovered d_q');
plot(tVec * 1000, vibFit.dHat, 'k-', 'LineWidth', 1.0, 'DisplayName', 'dual-frequency fit');
xlabel('Time (ms)');
ylabel('Shift (mm)');
title('Recovered shifts and fitted vibration');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
scatter(uTruth, blindResult.d, 14, 'b', 'filled');
xlabel('Truth u(t_q) (mm)');
ylabel('Recovered shift d_q (mm)');
title('Recovered shift vs truth');
grid on;

nexttile;
bar(categorical({'f1', 'f2'}), [summary.f_est(:), summary.f_true(:)]);
ylabel('Frequency (Hz)');
title('Estimated vs truth frequencies');
legend({'estimated', 'truth'}, 'Location', 'best', 'Box', 'off');

nexttile;
text(0.02, 0.82, sprintf('Mean frequency error: %.4f Hz', summary.mean_freq_error_Hz), 'FontSize', 11);
text(0.02, 0.62, sprintf('Mean amplitude error: %.4f mm', summary.mean_amp_error_mm), 'FontSize', 11);
text(0.02, 0.42, sprintf('Shift-fit RMSE: %.4e mm', summary.shift_rmse), 'FontSize', 11);
text(0.02, 0.22, sprintf('Segments used: %d', numel(segments)), 'FontSize', 11);
axis off;
title('Summary metrics');

fprintf('[Step 5] Blind-template branch final summary:\n');
disp(summary);
fprintf('[Step 5] Identifiability note:\n%s\n', summary.identifiability_note);
