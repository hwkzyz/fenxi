%% Step 5: Validation — compare blind ID results against ground truth
%  The algorithm never sees: template library, low-speed reference, true gap.
%  Validation checks:
%    1. Template reconstruction accuracy (vs true F_true)
%    2. Vibration displacement estimation accuracy (d_est vs d_true)
%    3. Vibration velocity estimation accuracy (dot_d_est vs dot_d_true)
%    4. Modal parameter accuracy (f, A, φ vs truth)
%    5. Displacement-velocity consistency (harmonic constraint check)

clc;
fprintf('=== Step 5: Validation against ground truth ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');
load(fullfile(outDir, 'step1_data.mat'), 'F_true');
load(fullfile(outDir, 'step2_passages.mat'), 'Passages');
load(fullfile(outDir, 'step3_blind_result.mat'), 'BlindResult');
load(fullfile(outDir, 'step4_vib_result.mat'), 'VibResult');

d_true     = Passages.d_true;
alpha_true = Passages.alpha_true;
dot_d_true = Passages.v0 * (alpha_true - 1);
t_passage  = Passages.t_centre;
xGrid      = Passages.x;

d_est      = BlindResult.d_est;
alpha_est  = BlindResult.alpha_est;
dot_d_est  = BlindResult.dot_d_est;
T_est      = BlindResult.T;

%% ---- 1. Template reconstruction assessment -----------------------
% Compare reconstructed template with true template
% Normalize both to [0, 1] range (absolute scale depends on unknown gap)
T_true_raw = F_true(xGrid);
T_est_norm = (T_est - min(T_est)) / (max(T_est) - min(T_est));
T_true_norm = (T_true_raw - min(T_true_raw)) / (max(T_true_raw) - min(T_true_raw));

% Shape similarity: Pearson correlation
T_corr = corr(T_est_norm, T_true_norm);
T_rmse_shape = rms(T_est_norm - T_true_norm);

fprintf('\n--- Template reconstruction ---\n');
fprintf('  Shape correlation with true:  %.4f\n', T_corr);
fprintf('  Normalized shape RMSE:        %.4f\n', T_rmse_shape);

%% ---- 2. Displacement estimation ----------------------------------
d_error = d_est - d_true;
d_rmse = rms(d_error);
d_mae = mean(abs(d_error));
d_corr = corr(d_est, d_true);
d_bias = mean(d_error);

fprintf('\n--- Displacement d_k estimation ---\n');
fprintf('  RMSE:      %.4f mm\n', d_rmse);
fprintf('  MAE:       %.4f mm\n', d_mae);
fprintf('  Bias:      %.4f mm (gauge absorbs static offset)\n', d_bias);
fprintf('  Correlation: %.4f\n', d_corr);

%% ---- 3. Velocity estimation --------------------------------------
dot_d_error = dot_d_est - dot_d_true;
dot_d_rmse = rms(dot_d_error);
dot_d_corr = corr(dot_d_est, dot_d_true);

fprintf('\n--- Velocity ḋ_k estimation ---\n');
fprintf('  RMSE:      %.1f mm/s\n', dot_d_rmse);
fprintf('  Correlation: %.4f\n', dot_d_corr);

%% ---- 4. Modal parameter accuracy ---------------------------------
% Match identified frequencies to true frequencies (nearest neighbour)
f_true = cfg.f_true;
A_true = cfg.A_true;
phi_true = cfg.phi_true;

f_est = VibResult.f_est;
A_est = VibResult.A_est;
phi_est = VibResult.phi_est;

% Match by frequency
[f_true_sorted, trueOrder] = sort(f_true);
A_true_sorted = A_true(trueOrder);
phi_true_sorted = phi_true(trueOrder);

f_err = zeros(size(f_true_sorted));
A_err = zeros(size(f_true_sorted));
phi_err = zeros(size(f_true_sorted));
matched = false(size(f_est));

for i = 1:numel(f_true_sorted)
    [~, j] = min(abs(f_est - f_true_sorted(i)));
    if ~matched(j)
        f_err(i) = abs(f_est(j) - f_true_sorted(i));
        A_err(i) = abs(A_est(j) - A_true_sorted(i));
        phi_err(i) = angular_diff(phi_est(j), phi_true_sorted(i));
        matched(j) = true;
    end
end

fprintf('\n--- Modal parameter accuracy ---\n');
for i = 1:numel(f_true_sorted)
    fprintf('  Mode %d:  Δf = %6.2f Hz,  ΔA = %.4f mm,  Δφ = %.4f rad\n', ...
        i, f_err(i), A_err(i), phi_err(i));
end
fprintf('  Mean Δf = %.2f Hz,  Mean ΔA = %.4f mm,  Mean Δφ = %.4f rad\n', ...
    mean(f_err), mean(A_err), mean(phi_err));

%% ---- 5. Self-consistency: displacement-velocity relation ---------
% From harmonic fit, compute ḋ from d parameters and compare with measured ḋ
% This checks internal consistency of the affine model
dot_d_from_harmonic = VibResult.dot_d_fit;
consistency_corr = corr(dot_d_est, dot_d_from_harmonic);
fprintf('\n--- Internal consistency ---\n');
fprintf('  corr(ḋ_measured, ḋ_harmonic): %.4f\n', consistency_corr);

%% ---- Package validation ------------------------------------------
Validation = struct();
Validation.d_rmse = d_rmse;
Validation.d_mae = d_mae;
Validation.d_corr = d_corr;
Validation.d_bias = d_bias;
Validation.dot_d_rmse = dot_d_rmse;
Validation.dot_d_corr = dot_d_corr;
Validation.f_err_Hz = f_err;
Validation.A_err_mm = A_err;
Validation.phi_err_rad = phi_err;
Validation.mean_f_err = mean(f_err);
Validation.mean_A_err = mean(A_err);
Validation.mean_phi_err = mean(phi_err);
Validation.T_shape_corr = T_corr;
Validation.T_shape_rmse = T_rmse_shape;
Validation.consistency_corr = consistency_corr;
Validation.nPassages = Passages.nPassages;
Validation.nRevs = cfg.NumRevs;

% Summary table
summaryVars = {'nPassages', 'nRevs', 'd_rmse', 'd_corr', 'dot_d_rmse', ...
    'dot_d_corr', 'mean_f_err', 'mean_A_err', 'mean_phi_err', ...
    'T_shape_corr', 'consistency_corr'};
summaryData = zeros(1, numel(summaryVars));
for i = 1:numel(summaryVars)
    val = Validation.(summaryVars{i});
    if numel(val) > 1, val = mean(val); end
    summaryData(i) = val;
end
SummaryTable = array2table(summaryData, 'VariableNames', summaryVars);
disp(SummaryTable);

save(fullfile(outDir, 'step5_validation.mat'), 'Validation', 'SummaryTable', '-v7.3');
fprintf('[Step 5] Validation saved.\n');

%% ---- Diagnostic plots --------------------------------------------
figure('Name', 'Step 5 - Validation', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 14]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) Template shape comparison
nexttile; hold on;
plot(xGrid, T_est_norm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'blind reconstructed');
plot(xGrid, T_true_norm, 'r--', 'LineWidth', 1.2, 'DisplayName', 'true (hidden gap)');
xlabel('x (mm)'); ylabel('Normalized capacitance');
title(sprintf('(a) Template shape:  corr = %.4f', T_corr));
legend('Location', 'best', 'Box', 'off');

% (b) d_k: estimated vs true
nexttile; hold on;
plot(d_true, d_est, 'b.', 'MarkerSize', 4);
xPlot = [min(d_true), max(d_true)];
plot(xPlot, xPlot, 'r--', 'LineWidth', 1.0);
xlabel('d_k true (mm)'); ylabel('d_k estimated (mm)');
title(sprintf('(b) Displacement: RMSE=%.4f, corr=%.4f', d_rmse, d_corr));

% (c) ḋ_k: estimated vs true
nexttile; hold on;
plot(dot_d_true, dot_d_est, 'b.', 'MarkerSize', 4);
xPlot = [min(dot_d_true), max(dot_d_true)];
plot(xPlot, xPlot, 'r--', 'LineWidth', 1.0);
xlabel('ḋ_k true (mm/s)'); ylabel('ḋ_k estimated (mm/s)');
title(sprintf('(c) Velocity: RMSE=%.1f, corr=%.4f', dot_d_rmse, dot_d_corr));

% (d) Displacement time series (first 200 passages)
nexttile; hold on;
nShow = min(200, numel(d_true));
tShow = (1:nShow)';
plot(tShow, d_true(1:nShow), 'b-', 'LineWidth', 1.0, 'DisplayName', 'true');
plot(tShow, d_est(1:nShow), 'r--', 'LineWidth', 1.0, 'DisplayName', 'estimated');
xlabel('Passage index'); ylabel('d_k (mm)');
title('(d) Displacement sequence');
legend('Location', 'best', 'Box', 'off');

% (e) Frequency identification
nexttile; hold on;
stem(f_true, A_true, 'b', 'LineWidth', 1.5, 'DisplayName', 'true');
stem(f_est, A_est, 'r--', 'LineWidth', 1.5, 'DisplayName', 'estimated');
xlabel('Frequency (Hz)'); ylabel('Amplitude (mm)');
title(sprintf('(e) Modal ID:  mean Δf=%.1f Hz, mean ΔA=%.3f mm', ...
    mean(f_err), mean(A_err)));
legend('Location', 'best', 'Box', 'off');

% (f) Error histograms
nexttile; hold on;
histogram(d_error, 30, 'FaceColor', [0.2 0.6 0.6], 'FaceAlpha', 0.6, ...
    'DisplayName', 'd error');
histogram(dot_d_error, 30, 'FaceColor', [0.85 0.35 0.10], 'FaceAlpha', 0.6, ...
    'DisplayName', 'ḋ error');
xline(0, 'k-', 'LineWidth', 0.8);
xlabel('Error');
ylabel('Count');
title('(f) Error distributions');
legend('Location', 'best', 'Box', 'off');

saveas(gcf, fullfile(outDir, 'fig_step5_validation.png'));

%% ============ Local functions ======================================
function delta = angular_diff(phi1, phi2)
    delta = mod(phi1 - phi2 + pi, 2*pi) - pi;
end
