%% Debug: analyse per-passage registration quality
%  Compares estimated (alpha, d) vs true for each sensor separately
%  to understand why d_k estimation is poor.

clc;
fprintf('=== Registration diagnostic ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');
load(fullfile(outDir, 'step1_data.mat'), 'F_true');
load(fullfile(outDir, 'step2_passages.mat'), 'Passages');
load(fullfile(outDir, 'step3_blind_result.mat'), 'BlindResult');

xGrid      = Passages.x;
d_true     = Passages.d_true;
alpha_true = Passages.alpha_true;
t_passage  = Passages.t_centre;
sensor_id  = Passages.sensor_id;
rev_id     = Passages.rev_id;
v0         = Passages.v0;

d_est     = BlindResult.d_est;
alpha_est = BlindResult.alpha_est;
T_est     = BlindResult.T;
T_xGrid   = BlindResult.T_xGrid;

nPassages = numel(d_true);

%% 1. Overall stats
fprintf('\n--- Overall ---\n');
fprintf('True d:   mean=%.4f, std=%.4f, range=[%.3f, %.3f] mm\n', ...
    mean(d_true), std(d_true), min(d_true), max(d_true));
fprintf('Est  d:   mean=%.4f, std=%.4f, range=[%.3f, %.3f] mm\n', ...
    mean(d_est), std(d_est), min(d_est), max(d_est));
fprintf('True alpha: std=%.4f\n', std(alpha_true));
fprintf('Est  alpha: std=%.4f\n', std(alpha_est));

fprintf('\nd correlation: %.4f\n', corr(d_est, d_true));
fprintf('alpha correlation: %.4f\n', corr(alpha_est, alpha_true));

%% 2. Per-sensor analysis
% Each sensor sees different vibration phases
for s = 1:3
    idx = sensor_id == s;
    d_t_s = d_true(idx);
    d_e_s = d_est(idx);
    a_t_s = alpha_true(idx);
    a_e_s = alpha_est(idx);

    fprintf('\n--- Sensor %d (%d passages) ---\n', s, sum(idx));
    fprintf('  d:  true std=%.4f,  est std=%.4f,  corr=%.4f\n', ...
        std(d_t_s), std(d_e_s), corr(d_e_s, d_t_s));
    fprintf('  alpha: true std=%.4f, est std=%.4f, corr=%.4f\n', ...
        std(a_t_s), std(a_e_s), corr(a_e_s, a_t_s));
end

%% 3. Registration cost landscape for a few sample passages
% Check if the cost function has a clear minimum
fprintf('\n--- Cost landscape analysis (5 random passages) ---\n');
rng(1);
sampleIdx = randperm(nPassages, min(5, nPassages));

for si = 1:numel(sampleIdx)
    k = sampleIdx(si);
    yk = Passages.Y(k, :)';
    activeMask = (yk - min(yk)) > 0.1 * range(yk);

    % Sweep alpha and d around the true values
    a_test = linspace(alpha_true(k) - 0.1, alpha_true(k) + 0.1, 41);
    d_test = linspace(d_true(k) - 0.5, d_true(k) + 0.5, 41);
    J_sweep = zeros(numel(a_test), numel(d_test));

    for ia = 1:numel(a_test)
        for id = 1:numel(d_test)
            xTrans = a_test(ia) * xGrid(activeMask) - d_test(id);
            T_pred = interp1(T_xGrid, T_est, xTrans, 'pchip', NaN);
            valid = isfinite(T_pred);
            J_sweep(ia, id) = mean((yk(activeMask(valid)) - T_pred(valid)).^2);
        end
    end

    [~, ia_min] = min(min(J_sweep, [], 2));
    [~, id_min] = min(min(J_sweep, [], 1));
    alpha_found = a_test(ia_min);
    d_found = d_test(id_min);

    fprintf('  Passage %d (sensor %d):\n', k, sensor_id(k));
    fprintf('    True:     alpha=%.4f, d=%.4f\n', alpha_true(k), d_true(k));
    fprintf('    Grid opt: alpha=%.4f, d=%.4f\n', alpha_found, d_found);
    fprintf('    Est:      alpha=%.4f, d=%.4f\n', alpha_est(k), d_est(k));

    % Check curvature
    J_center = mean(J_sweep(ia_min, :));
    J_plus_d = J_sweep(ia_min, min(id_min + 5, numel(d_test)));
    curv_d = (J_plus_d - J_center) / (d_test(min(id_min + 5, numel(d_test))) - d_test(id_min))^2;
    fprintf('    d-direction curvature: ~%.2e\n', curv_d);
end

%% 4. Check if alpha estimates are just absorbing residual d error
% Theoretically, if α = 1 is used (no velocity), d absorbs all waveform variation
% Check: what fraction of waveform variation is captured by d alone vs d+α?
fprintf('\n--- Model fit quality: d-only vs full affine ---\n');
for si = 1:min(3, nPassages)
    k = sampleIdx(si);
    yk = Passages.Y(k, :)';

    % d-only fit (α = 1)
    x_d = xGrid - d_true(k);
    y_fit_d = interp1(T_xGrid, T_est, x_d, 'pchip', NaN);
    err_d = rms(yk - y_fit_d, 'omitnan');

    % Full affine fit
    x_a = alpha_true(k) * xGrid - d_true(k);
    y_fit_a = interp1(T_xGrid, T_est, x_a, 'pchip', NaN);
    err_a = rms(yk - y_fit_a, 'omitnan');

    % Noise level
    err_noise = rms(yk - F_true(xGrid - d_true(k)), 'omitnan');

    fprintf('  Passage %d: d-only RMS=%.4e, affine RMS=%.4e, noise=%.4e\n', ...
        k, err_d, err_a, err_noise);
end
