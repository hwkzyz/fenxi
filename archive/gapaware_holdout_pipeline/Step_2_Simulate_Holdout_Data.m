%% Step 2: simulate low-speed reference and high-speed hidden-gap data
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), ...
    'gapList', 'xCell', 'yCell', 'templateLib');

rng(7);

idxLow = find(abs(gapList - cfg.g_low) < 1e-12, 1);
idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
F_low_true = griddedInterpolant(xCell{idxLow}, yCell{idxLow}, 'pchip', 'nearest');
F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

Data_Low = simulate_rotating_waveform_from_template(F_low_true, cfg.RPM_low, ...
    cfg.NumRevs_low, cfg.fs, cfg.R_tip, cfg.alpha_k, cfg.noiseRatio, ...
    templateLib.domain, [0, 0], [0, 0], [0, 0]);

Data_High = simulate_rotating_waveform_from_template(F_high_true, cfg.RPM_high, ...
    cfg.NumRevs_high, cfg.fs, cfg.R_tip, cfg.alpha_k, cfg.noiseRatio, ...
    templateLib.domain, cfg.A_true, cfg.f_true, cfg.phi_true);

truth = struct();
truth.g_high = cfg.g_holdout;
truth.g_low = cfg.g_low;
truth.A = cfg.A_true;
truth.f = cfg.f_true;
truth.phi = cfg.phi_true;
truth.u_t = Data_High.u_t;

save(fullfile(outDir, 'stage2_simulated_data.mat'), ...
    'Data_Low', 'Data_High', 'truth', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Step 2 - Simulated Data', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 20, 12]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    idx = Data_Low.t <= min(Data_Low.t) + 2 * 60 / cfg.RPM_low;
    plot(Data_Low.t(idx) * 1000, Data_Low.V_cap(idx), 'Color', [0.55 0.55 0.55]);
    plot(Data_Low.t(idx) * 1000, Data_Low.V_cap_clean(idx), 'b-', 'LineWidth', 1.2);
    plot(Data_Low.t(idx) * 1000, Data_Low.V_OPR(idx) / max(Data_Low.V_OPR(idx)) * max(Data_Low.V_cap_clean(idx)), 'r-');
    ylabel('Voltage');
    title(sprintf('Low-speed no-vibration reference, g = %.1f mm', cfg.g_low));
    legend({'noisy', 'clean', 'OPR scaled'}, 'Location', 'best', 'Box', 'off');

    nexttile; hold on;
    idx = Data_High.t <= min(Data_High.t) + 4 * 60 / cfg.RPM_high;
    plot(Data_High.t(idx) * 1000, Data_High.V_cap(idx), 'Color', [0.55 0.55 0.55]);
    plot(Data_High.t(idx) * 1000, Data_High.V_cap_clean(idx), 'b-', 'LineWidth', 1.2);
    plot(Data_High.t(idx) * 1000, Data_High.V_OPR(idx) / max(Data_High.V_OPR(idx)) * max(Data_High.V_cap_clean(idx)), 'r-');
    ylabel('Voltage');
    title(sprintf('High-speed vibrating waveform from hidden gap, g = %.1f mm', cfg.g_holdout));
    legend({'noisy', 'clean', 'OPR scaled'}, 'Location', 'best', 'Box', 'off');

    nexttile;
    plot(Data_High.t(idx) * 1000, Data_High.u_t(idx), 'k-');
    xlabel('Time (ms)');
    ylabel('u(t) (mm)');
    title('True two-frequency vibration');
end

fprintf('[Step 2] Simulated hidden-gap high-speed data. Truth gap = %.3f mm\n', truth.g_high);

function Data = simulate_rotating_waveform_from_template(Fx, RPM, NumRevs, fs, R_tip, alpha_k, noiseRatio, domain, A, f, phi)
Omega = RPM * 2*pi / 60;
T = 2*pi / Omega;
dt = 1 / fs;
tEnd = (NumRevs + 1) * T;
t = (0:dt:tEnd)';
oprDelay = 0.08 * T;
T_opr = oprDelay + (0:NumRevs)' * T;
xHalf = 0.52 * diff(domain);
probe = linspace(domain(1), domain(2), 400)';
yProbe = Fx(probe);
baseline = min(yProbe);
V_clean = baseline * ones(size(t));
u_t = zeros(size(t));
for im = 1:numel(A)
    u_t = u_t + A(im) * sin(2*pi*f(im)*t + phi(im));
end
for m = 1:NumRevs
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = abs(t - Te) <= xHalf / (Omega * R_tip);
        xNom = Omega * R_tip * (t(idx) - Te);
        xPhys = xNom - u_t(idx);
        inDomain = xPhys >= domain(1) & xPhys <= domain(2);
        idxAll = find(idx);
        idxUse = idxAll(inDomain);
        V_clean(idxUse) = max(V_clean(idxUse), Fx(xPhys(inDomain)));
    end
end
noiseSigma = noiseRatio * max(range(yProbe), eps);
V_cap = V_clean + noiseSigma * randn(size(V_clean));
V_OPR = zeros(size(t));
oprWidth = 0.006 * T;
for m = 1:numel(T_opr)
    V_OPR(abs(t - T_opr(m)) <= oprWidth / 2) = 5.0;
end
V_OPR = V_OPR + 0.02 * randn(size(V_OPR));
Data = struct('t', t, 'V_cap', V_cap, 'V_cap_clean', V_clean, ...
    'V_OPR', V_OPR, 'RPM', RPM, 'T_opr_truth', T_opr(:), ...
    'u_t', u_t, 'baseline', baseline);
end

