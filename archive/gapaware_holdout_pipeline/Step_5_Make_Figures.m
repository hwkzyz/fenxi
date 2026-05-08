%% Step 5: final visual check of the gap-aware holdout result
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'gapList', 'xCell', 'yCell');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage3_gap_estimation.mat'), 'highMap', 'staticState');
load(fullfile(outDir, 'stage4_identification_results.mat'), ...
    'fixedResult', 'rawScanResult', 'gapAwareResult', 'oracleResult', 'summaryTable');

figure('Name', 'Step 5 - Final Visual Summary', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 14]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for i = 1:numel(gapList)
    if abs(gapList(i) - truth.g_high) < 1e-12
        plot(xCell{i}, yCell{i}, 'r--', 'LineWidth', 1.6, 'DisplayName', 'hidden high-speed gap');
    elseif abs(gapList(i) - truth.g_low) < 1e-12
        plot(xCell{i}, yCell{i}, 'b-', 'LineWidth', 1.4, 'DisplayName', 'low-speed reference gap');
    else
        plot(xCell{i}, yCell{i}, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    end
end
xlabel('Position x (mm)');
ylabel('Capacitance');
title('Holdout design');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(staticState.gapQueryGrid, staticState.J_static, 'b-', 'DisplayName', 'static cost');
xline(staticState.gHat, 'g-', 'LineWidth', 1.2, 'DisplayName', sprintf('gHat = %.3f', staticState.gHat));
xline(truth.g_high, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('truth = %.3f', truth.g_high));
xlabel('Candidate gap g (mm)');
ylabel('Cost');
title(sprintf('Estimated static gap before VARPRO (%s)', staticState.supportType));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
idxPlot = 1:min(1200, numel(highMap.t_v));
plot(highMap.t_v(idxPlot) * 1000, highMap.V_a(idxPlot), 'k.', 'MarkerSize', 4, 'DisplayName', 'observed');
plot(highMap.t_v(idxPlot) * 1000, fixedResult.VFit(idxPlot), 'r-', 'DisplayName', 'fixed low gap');
plot(highMap.t_v(idxPlot) * 1000, rawScanResult.VFit(idxPlot), 'c--', 'DisplayName', 'raw-scan initial');
plot(highMap.t_v(idxPlot) * 1000, gapAwareResult.VFit(idxPlot), 'b-', 'DisplayName', 'gap-aware (joint refined)');
plot(highMap.t_v(idxPlot) * 1000, oracleResult.VFit(idxPlot), 'g--', 'DisplayName', 'oracle');
xlabel('Time (ms)');
ylabel('Voltage');
title('Final waveform fit');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
tCommon = gapAwareResult.t;
uTruth = cfg.A_true(1) * sin(2*pi*cfg.f_true(1)*tCommon + cfg.phi_true(1)) + ...
    cfg.A_true(2) * sin(2*pi*cfg.f_true(2)*tCommon + cfg.phi_true(2));
plot(tCommon(idxPlot) * 1000, uTruth(idxPlot), 'k-', 'DisplayName', 'truth');
plot(tCommon(idxPlot) * 1000, fixedResult.uHat(idxPlot), 'r-', 'DisplayName', 'fixed low gap');
plot(tCommon(idxPlot) * 1000, rawScanResult.uHat(idxPlot), 'c--', 'DisplayName', 'raw-scan initial');
plot(tCommon(idxPlot) * 1000, gapAwareResult.uHat(idxPlot), 'b-', 'DisplayName', 'gap-aware (joint refined)');
xlabel('Time (ms)');
ylabel('u(t) (mm)');
title('Vibration reconstruction');
legend('Location', 'best', 'Box', 'off');

fprintf('[Step 5] Final summary:\n');
disp(summaryTable);
fprintf('Hidden gap truth = %.3f mm, initial gHat = %.3f mm, final gap-aware g = %.3f mm\n', ...
    truth.g_high, staticState.gHat, gapAwareResult.g_used);
