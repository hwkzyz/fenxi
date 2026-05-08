%% Step 03: displacement reconstruction with the fixed trusted window.
%
% This post-processing step converts the final identified vibration
% parameters into a displacement time history and compares it with the
% simulation truth.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gTrue', 'var') || isempty(gTrue)
    gTrue = 0.2;
end

D = compute_inv_log_2_fixed_demo(ctx, gTrue);
t = D.highMap.t_v(:);
uIdentified = fit_u(D.result.fit.p, t);
uTruth = interp1(D.dataHigh.t(:), D.dataHigh.u_truth(:), t, ...
    'linear', 'extrap');
uError = uIdentified - uTruth;

tDense = linspace(min(t), max(t), 4000)';
uTruthDense = interp1(D.dataHigh.t(:), D.dataHigh.u_truth(:), tDense, ...
    'linear', 'extrap');
uIdentifiedDense = fit_u(D.result.fit.p, tDense);
uErrorDense = uIdentifiedDense - uTruthDense;
VFit = D.result.VFit(:);
voltageResidual = D.highMap.V_a(:) - VFit;
previewMask = tDense <= min(tDense) + 0.018;

[fTrue, trueOrder] = sort(D.cfgCase.f_true(:));
ATrue = D.cfgCase.A_true(trueOrder).';
fIdentified = D.result.f_id(:);
AIdentified = D.result.A_id(:);

displacementTable = table(t, uTruth, uIdentified, uError, ...
    D.highMap.rev_v(:), D.highMap.S_v(:), ...
    'VariableNames', {'t_s','u_truth_mm','u_identified_mm', ...
    'u_error_mm','rev_index','sensor_index'});

summary = table(D.gTrue, D.result.g_used, abs(D.result.g_used - D.gTrue), ...
    D.result.f_id(1), D.result.f_id(2), D.result.A_id(1), D.result.A_id(2), ...
    rms(uError), max(abs(uError)), D.result.rmse, ...
    'VariableNames', {'g_true_mm','g_identified_mm','gap_error_mm', ...
    'f1_Hz','f2_Hz','A1_mm','A2_mm','u_rmse_mm','u_max_abs_error_mm', ...
    'voltage_rmse_V'});

outCsv = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_displacement.csv');
outSummary = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_displacement_summary.csv');
outMat = fullfile(ctx.outDir, 'inv_log_2_fixed_trust_displacement.mat');
writetable(displacementTable, outCsv);
writetable(summary, outSummary);
save(outMat, 'displacementTable', 'summary', 'D', '-v7.3');

disp(summary);
fprintf('Saved displacement reconstruction:\n%s\n', outCsv);

fig = figure('Name', 'fixed trust displacement reconstruction', 'Color', 'w');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, tDense(previewMask), uTruthDense(previewMask), 'k--', ...
    'LineWidth', 1.1, 'DisplayName', 'truth');
hold(ax1, 'on');
plot(ax1, tDense(previewMask), uIdentifiedDense(previewMask), '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.2, 'DisplayName', 'identified');
ylabel(ax1, 'Displacement (mm)', 'FontSize', 9);
title(ax1, '(a) Continuous displacement reconstruction', ...
    'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
plot(ax2, tDense(previewMask), uErrorDense(previewMask), '-', ...
    'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(ax2, 0, 'k--', 'LineWidth', 0.8);
xlabel(ax2, 'Time t (s)', 'FontSize', 9);
ylabel(ax2, 'Error (mm)', 'FontSize', 9);
title(ax2, '(b) Continuous displacement error', ...
    'FontWeight', 'normal', 'FontSize', 9);

ax3 = nexttile;
plot(ax3, uTruth, uIdentified, 'k.', 'MarkerSize', 5);
hold(ax3, 'on');
lims = [min([uTruth; uIdentified]), max([uTruth; uIdentified])];
plot(ax3, lims, lims, 'k--', 'LineWidth', 0.8);
xlabel(ax3, 'Truth displacement (mm)', 'FontSize', 9);
ylabel(ax3, 'Identified displacement (mm)', 'FontSize', 9);
title(ax3, '(c) Sample-wise displacement agreement', ...
    'FontWeight', 'normal', 'FontSize', 9);

ax4 = nexttile;
edges = linspace(min(uError), max(uError), 28);
histogram(ax4, uError, edges, 'FaceColor', [0.00 0.45 0.74], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.85);
xline(ax4, 0, 'k--', 'LineWidth', 0.8);
xlabel(ax4, 'Displacement error (mm)', 'FontSize', 9);
ylabel(ax4, 'Count', 'FontSize', 9);
title(ax4, '(d) Error distribution at fitting samples', ...
    'FontWeight', 'normal', 'FontSize', 9);

ax5 = nexttile;
bar(ax5, categorical({'Mode 1','Mode 2'}), [fTrue(:), fIdentified(:)]);
ylabel(ax5, 'Frequency (Hz)', 'FontSize', 9);
title(ax5, '(e) Identified frequencies', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax5, {'truth','identified'}, 'Location', 'best');

ax6 = nexttile;
bar(ax6, categorical({'Mode 1','Mode 2'}), [ATrue(:), AIdentified(:)]);
ylabel(ax6, 'Amplitude (mm)', 'FontSize', 9);
title(ax6, '(f) Identified amplitudes', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax6, {'truth','identified'}, 'Location', 'best');

fig2 = figure('Name', 'fixed trust voltage residual after displacement reconstruction', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax7 = nexttile;
[tSort, sortIdx] = sort(t, 'ascend');
plot(ax7, tSort, D.highMap.V_a(sortIdx), 'k.', 'MarkerSize', 4, ...
    'DisplayName', 'high-speed samples');
hold(ax7, 'on');
plot(ax7, tSort, VFit(sortIdx), '-', 'Color', [0.47 0.67 0.19], ...
    'LineWidth', 1.1, 'DisplayName', 'voltage fit');
ylabel(ax7, 'Voltage (V)', 'FontSize', 9);
title(ax7, '(a) Voltage fit after displacement reconstruction', ...
    'FontWeight', 'normal', 'FontSize', 9);
legend(ax7, 'Location', 'best');

ax8 = nexttile;
plot(ax8, tSort, voltageResidual(sortIdx), '.', 'Color', [0.85 0.33 0.10], ...
    'MarkerSize', 4);
yline(ax8, 0, 'k--', 'LineWidth', 0.8);
xlabel(ax8, 'Time t (s)', 'FontSize', 9);
ylabel(ax8, 'Residual (V)', 'FontSize', 9);
title(ax8, '(b) Voltage residual', 'FontWeight', 'normal', 'FontSize', 9);

apply_inv_log_2_figure_style(ax1, 180, 150);
apply_inv_log_2_figure_style(ax7, 170, 110);
