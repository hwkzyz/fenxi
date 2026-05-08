%% Step 4: Vibration parameter extraction from {d_k} and {dot_d_k}
%  - Spectral analysis of displacement and velocity sequences
%  - Multi-harmonic model fitting
%  - Consistency check between displacement and velocity channels
%
%  The two sequences are NOT independent — they come from the same d(t).
%  This redundancy is exploited for robust modal parameter estimation.

clc;
fprintf('=== Step 4: Vibration parameter extraction ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');
load(fullfile(outDir, 'step2_passages.mat'), 'Passages');
load(fullfile(outDir, 'step3_blind_result.mat'), 'BlindResult');

d_est     = BlindResult.d_est;       % estimated displacement (mm)
dot_d_est = BlindResult.dot_d_est;    % estimated velocity (mm/s)
t_passage = BlindResult.t_centre;     % passage times (s)
sensor_id = BlindResult.sensor_id;
v0        = BlindResult.v0;

nPassages = numel(d_est);
nSensors  = numel(unique(sensor_id));

fprintf('Passages: %d,  sensors: %d\n', nPassages, nSensors);

%% ---- 1. Spectral analysis ----------------------------------------
% Use Lomb-Scargle periodogram (handles non-uniform sampling)
freqAxis = linspace(cfg.vib.freqSearchHz(1), cfg.vib.freqSearchHz(2), 4000)';

% Displacement spectrum
[Pxx_d, ~] = plomb_custom(d_est - mean(d_est), t_passage, freqAxis);
Pxx_d = Pxx_d / max(Pxx_d);

% Velocity spectrum
[Pxx_v, ~] = plomb_custom(dot_d_est - mean(dot_d_est), t_passage, freqAxis);
Pxx_v = Pxx_v / max(Pxx_v);

% Combined spectrum (product emphasizes shared frequencies)
Pxx_combined = sqrt(Pxx_d .* Pxx_v);
Pxx_combined = Pxx_combined / max(Pxx_combined);

%% ---- 2. Peak detection → mode frequencies ------------------------
[peakFreqs, peakAmps] = detect_spectral_peaks(freqAxis, Pxx_combined, cfg.vib.maxModes);
M = numel(peakFreqs);
fprintf('Detected %d spectral peaks:\n', M);
for m = 1:M
    fprintf('  f_%d = %.1f Hz  (rel amp = %.3f)\n', m, peakFreqs(m), peakAmps(m));
end

%% ---- 3. Multi-harmonic fitting -----------------------------------
% Model:  d(t) = sum_m A_m sin(2π f_m t + φ_m)
%          ḋ(t) = sum_m 2π f_m A_m cos(2π f_m t + φ_m)
%
% Joint fit: use BOTH displacement AND velocity data to constrain
% the shared parameters {A_m, f_m, φ_m}.

% Step 3a: nonlinear refinement of frequencies (starting from spectral peaks)
f_init = peakFreqs(:)';
f_opt = f_init;
[pAll_opt, f_opt, cost_trace] = refine_frequencies_joint(...
    d_est, dot_d_est, t_passage, f_init, v0);

% Step 3b: final linear fit of amplitudes and phases at refined frequencies
[A_opt, phi_opt, d_fit, dot_d_fit] = fit_amplitudes_phases(...
    d_est, dot_d_est, t_passage, f_opt);

% Sort by frequency
[f_opt, order] = sort(f_opt);
A_opt = A_opt(order);
phi_opt = phi_opt(order);

%% ---- 4. Quality metrics ------------------------------------------
% Residual RMS
res_d = d_est - d_fit;
res_dot = dot_d_est - dot_d_fit;
rms_d = rms(res_d);
rms_dot = rms(res_dot);

% R² (coefficient of determination)
R2_d = 1 - sum(res_d.^2) / sum((d_est - mean(d_est)).^2);
R2_dot = 1 - sum(res_dot.^2) / sum((dot_d_est - mean(dot_d_est)).^2);

fprintf('\nHarmonic fit quality:\n');
fprintf('  d(t) fit:   RMS = %.4f mm,  R^2 = %.4f\n', rms_d, R2_d);
fprintf('  ḋ(t) fit:   RMS = %.1f mm/s, R^2 = %.4f\n', rms_dot, R2_dot);
fprintf('Identified parameters:\n');
for m = 1:numel(f_opt)
    fprintf('  Mode %d:  f = %8.2f Hz,  A = %.4f mm,  φ = %+.3f rad\n', ...
        m, f_opt(m), A_opt(m), phi_opt(m));
end

%% ---- Package results ---------------------------------------------
VibResult = struct();
VibResult.f_est     = f_opt;
VibResult.A_est     = A_opt;
VibResult.phi_est   = phi_opt;
VibResult.M         = numel(f_opt);
VibResult.d_fit     = d_fit;
VibResult.dot_d_fit = dot_d_fit;
VibResult.rms_d     = rms_d;
VibResult.rms_dot   = rms_dot;
VibResult.R2_d      = R2_d;
VibResult.R2_dot    = R2_dot;
VibResult.freqAxis  = freqAxis;
VibResult.Pxx_d     = Pxx_d;
VibResult.Pxx_v     = Pxx_v;
VibResult.Pxx_combined = Pxx_combined;
VibResult.peakFreqs_raw = peakFreqs;
VibResult.d_est     = d_est;
VibResult.dot_d_est = dot_d_est;
VibResult.t_passage = t_passage;

save(fullfile(outDir, 'step4_vib_result.mat'), 'VibResult', '-v7.3');
fprintf('[Step 4] Vibration result saved.\n');

%% ---- Diagnostic plots --------------------------------------------
figure('Name', 'Step 4 - Vibration Analysis', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 14]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) Spectra
nexttile; hold on;
plot(freqAxis, Pxx_d, 'b-', 'DisplayName', 'displacement', 'LineWidth', 1.0);
plot(freqAxis, Pxx_v, 'r-', 'DisplayName', 'velocity', 'LineWidth', 1.0);
plot(freqAxis, Pxx_combined, 'k-', 'DisplayName', 'combined', 'LineWidth', 1.5);
for m = 1:numel(f_opt)
    xline(f_opt(m), 'g--', 'LineWidth', 1.0);
end
xlabel('Frequency (Hz)'); ylabel('Normalized power');
title('(a) Lomb-Scargle spectra');
legend('Location', 'best', 'Box', 'off');

% (b) d_k time series and fit (first 500 passages)
nexttile; hold on;
nShow = min(500, nPassages);
plot(1:nShow, d_est(1:nShow), 'b.', 'MarkerSize', 4, 'DisplayName', 'estimated');
plot(1:nShow, d_fit(1:nShow), 'r-', 'LineWidth', 1.0, 'DisplayName', 'harmonic fit');
xlabel('Passage index'); ylabel('d_k (mm)');
title(sprintf('(b) Displacement: RMS=%.3f mm, R^2=%.3f', rms_d, R2_d));
legend('Location', 'best', 'Box', 'off');

% (c) ḋ_k time series and fit
nexttile; hold on;
plot(1:nShow, dot_d_est(1:nShow), 'b.', 'MarkerSize', 4, 'DisplayName', 'estimated');
plot(1:nShow, dot_d_fit(1:nShow), 'r-', 'LineWidth', 1.0, 'DisplayName', 'harmonic fit');
xlabel('Passage index'); ylabel('ḋ_k (mm/s)');
title(sprintf('(c) Velocity: RMS=%.1f mm/s, R^2=%.3f', rms_dot, R2_dot));
legend('Location', 'best', 'Box', 'off');

% (d) Displacement residual
nexttile;
plot(1:nShow, res_d(1:nShow), 'k.', 'MarkerSize', 3);
xlabel('Passage index'); ylabel('Residual (mm)');
title('(d) Displacement residual');

% (e) Phase portrait: d vs ḋ
nexttile; hold on;
plot(d_est, dot_d_est, 'b.', 'MarkerSize', 3, 'DisplayName', 'estimated');
plot(d_fit, dot_d_fit, 'r.', 'MarkerSize', 2, 'DisplayName', 'harmonic fit');
xlabel('d (mm)'); ylabel('ḋ (mm/s)');
title('(e) Phase portrait');
legend('Location', 'best', 'Box', 'off');

% (f) Lissajous between two dominant modes
nexttile; hold on;
if M >= 2
    % Extract individual mode contributions
    d1 = A_opt(1) * sin(2*pi*f_opt(1)*t_passage + phi_opt(1));
    d2 = A_opt(2) * sin(2*pi*f_opt(2)*t_passage + phi_opt(2));
    plot(d1(1:nShow), d2(1:nShow), 'k.', 'MarkerSize', 3);
    xlabel(sprintf('Mode 1: %.0f Hz (mm)', f_opt(1)));
    ylabel(sprintf('Mode 2: %.0f Hz (mm)', f_opt(2)));
    title('(f) Mode 1 vs Mode 2');
end

saveas(gcf, fullfile(outDir, 'fig_step4_vibration.png'));

%% ============ Local functions ======================================
function [Pxx, freq] = plomb_custom(x, t, freqAxis)
    % Simplified Lomb-Scargle periodogram
    % x: signal, t: observation times (non-uniform), freqAxis: frequencies to evaluate
    x = x(:); t = t(:);
    x = x - mean(x);
    omega = 2*pi * freqAxis(:)';
    tau = 0.5 * atan2(sum(sin(2*omega.*t), 1), sum(cos(2*omega.*t), 1)) ./ omega;
    tau(isnan(tau)) = 0;
    cosTerm = cos(omega .* (t - tau));
    sinTerm = sin(omega .* (t - tau));
    P_num = (sum(x .* cosTerm, 1).^2 ./ sum(cosTerm.^2, 1) + ...
             sum(x .* sinTerm, 1).^2 ./ sum(sinTerm.^2, 1));
    P_den = sum(x.^2);
    Pxx = (P_num ./ max(P_den, eps))';
    freq = freqAxis(:);
end

function [peakFreqs, peakAmps] = detect_spectral_peaks(freqAxis, Pxx, maxModes)
    % Detect dominant peaks in spectrum with minimum separation
    minSep = 40;  % Hz, minimum separation between peaks
    [pks, locs] = findpeaks(Pxx, 'MinPeakHeight', 0.05, ...
        'MinPeakDistance', max(1, minSep / (freqAxis(2) - freqAxis(1))));
    [~, order] = sort(pks, 'descend');
    locs = locs(order);
    pks = pks(order);
    nKeep = min(maxModes, numel(pks));
    peakFreqs = freqAxis(locs(1:nKeep));
    peakAmps = pks(1:nKeep);
end

function [pAll_opt, f_opt, cost_trace] = refine_frequencies_joint(...
        d_obs, dot_d_obs, t, f_init, v0)
    % Joint nonlinear refinement of frequencies using displacement + velocity
    % p = [A1, phi1, f1, A2, phi2, f2, ...]
    M = numel(f_init);
    p0 = zeros(3*M, 1);

    % Initial amplitudes from RMS of filtered signal around each frequency
    for m = 1:M
        % Crude amplitude estimate
        f0 = f_init(m);
        % Bandpass filter around f0 (simple)
        bw = 20; % Hz
        % Use rough amplitude estimate
        A_guess = sqrt(2) * rms(d_obs) / sqrt(M);
        p0(3*m-2) = A_guess;
        p0(3*m-1) = 0;  % phase (will be refined)
        p0(3*m)   = f0;
    end

    % Bounded optimization
    lb = zeros(3*M, 1);
    ub = zeros(3*M, 1);
    for m = 1:M
        lb(3*m-2) = 0.001;  ub(3*m-2) = 1.0;     % amplitude
        lb(3*m-1) = -pi;    ub(3*m-1) = pi;       % phase
        lb(3*m)   = f_init(m) - 30;  ub(3*m) = f_init(m) + 30;  % frequency
    end

    % Cost function: weighted joint fit
    w_d = 1 / max(var(d_obs), 1e-6);
    w_v = v0^2 / max(var(dot_d_obs), 1e-6);  % scale velocity to displacement units
    w_v = w_v / 10;  % downweight velocity (noisier)

    costFun = @(p) joint_harmonic_cost(p, d_obs, dot_d_obs, t, w_d, w_v);

    % Optimize
    p0_clamped = min(max(p0, lb), ub);
    opts = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000, ...
        'TolX', 1e-8, 'TolFun', 1e-10);
    [pAll_opt, ~] = fminsearch(costFun, p0_clamped, opts);

    % Extract frequencies
    f_opt = zeros(M, 1);
    for m = 1:M
        f_opt(m) = pAll_opt(3*m);
    end
    cost_trace = NaN;  % simplified
end

function J = joint_harmonic_cost(p, d_obs, dot_d_obs, t, w_d, w_v)
    M = numel(p) / 3;
    d_model = zeros(size(d_obs));
    dot_d_model = zeros(size(dot_d_obs));
    for m = 1:M
        A = p(3*m-2);
        phi = p(3*m-1);
        f = p(3*m);
        d_model = d_model + A * sin(2*pi*f*t + phi);
        dot_d_model = dot_d_model + 2*pi*f*A * cos(2*pi*f*t + phi);
    end
    J = w_d * mean((d_obs - d_model).^2) + w_v * mean((dot_d_obs - dot_d_model).^2);
end

function [A_opt, phi_opt, d_fit, dot_d_fit] = fit_amplitudes_phases(...
        d_obs, dot_d_obs, t, f_opt)
    % Given frequencies, fit amplitudes and phases via linear least squares
    M = numel(f_opt);
    n = numel(d_obs);

    % Design matrix for displacement: [sin(2π f_m t), cos(2π f_m t)]
    X = zeros(n, 2*M);
    for m = 1:M
        X(:, 2*m-1) = sin(2*pi*f_opt(m)*t);
        X(:, 2*m)   = cos(2*pi*f_opt(m)*t);
    end

    % Solve for displacement channel
    coeffs_d = X \ d_obs;

    % Convert [sin_coeff, cos_coeff] → (A, φ)
    A_opt = zeros(M, 1);
    phi_opt = zeros(M, 1);
    for m = 1:M
        s = coeffs_d(2*m-1);
        c = coeffs_d(2*m);
        A_opt(m) = sqrt(s^2 + c^2);
        phi_opt(m) = atan2(c, s);  % A sin(ωt + φ) = A sin ωt cos φ + A cos ωt sin φ
                                   % = (A cos φ) sin ωt + (A sin φ) cos ωt
                                   % so s = A cos φ, c = A sin φ → φ = atan2(c, s)
    end

    % Reconstruct fitted signals
    d_fit = X * coeffs_d;
    dot_d_fit = zeros(size(dot_d_obs));
    for m = 1:M
        dot_d_fit = dot_d_fit + 2*pi*f_opt(m)*A_opt(m) * cos(2*pi*f_opt(m)*t + phi_opt(m));
    end
end
