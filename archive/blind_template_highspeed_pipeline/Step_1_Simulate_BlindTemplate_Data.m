%% Step 1: simulate high-speed data from the hidden static gap
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
if isempty(idxHoldout)
    error('Hidden high-speed gap %.3f mm does not exist in the stacked static curves.', cfg.g_holdout);
end

xHold = xCell{idxHoldout}(:);
yHold = yCell{idxHoldout}(:);
xMin = min(xHold);
xMax = max(xHold);
xPad = cfg.domainPadRatio * (xMax - xMin);
domain = [xMin + xPad, xMax - xPad];

F_high_true = griddedInterpolant(xHold, yHold, 'pchip', 'nearest');
Data_High = simulate_rotating_waveform_from_template(F_high_true, cfg.RPM_high, ...
    cfg.NumRevs_high, cfg.fs, cfg.R_tip, cfg.alpha_k, cfg.noiseRatio, ...
    domain, cfg.A_true, cfg.f_true, cfg.phi_true);

truth = struct();
truth.g_high = cfg.g_holdout;
truth.A = cfg.A_true(:).';
truth.f = cfg.f_true(:).';
truth.phi = cfg.phi_true(:).';
truth.domain = domain;
truth.templateX = xHold;
truth.templateY = yHold;
truth.u_t = Data_High.u_t;

save(fullfile(outDir, 'stage1_simulated_data.mat'), ...
    'gapList', 'xCell', 'yCell', 'Data_High', 'truth', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Blind-Template Step 1 - Simulated Data', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 20, 12]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    idxPlot = Data_High.t <= min(Data_High.t) + 4 * 60 / cfg.RPM_high;
    nexttile; hold on;
    plot(Data_High.t(idxPlot) * 1000, Data_High.V_cap(idxPlot), 'Color', [0.55 0.55 0.55]);
    plot(Data_High.t(idxPlot) * 1000, Data_High.V_cap_clean(idxPlot), 'b-', 'LineWidth', 1.2);
    plot(Data_High.t(idxPlot) * 1000, Data_High.V_OPR(idxPlot) / max(Data_High.V_OPR(idxPlot)) * max(Data_High.V_cap_clean(idxPlot)), 'r-');
    ylabel('Voltage');
    title(sprintf('High-speed vibrating waveform from hidden gap, g = %.1f mm', cfg.g_holdout));
    legend({'noisy', 'clean', 'OPR scaled'}, 'Location', 'best', 'Box', 'off');

    nexttile; hold on;
    plot(xHold, yHold, 'k-', 'LineWidth', 1.2);
    xline(domain(1), 'r--', 'LineWidth', 1.0);
    xline(domain(2), 'r--', 'LineWidth', 1.0);
    xlabel('Position x (mm)');
    ylabel('Capacitance');
    title('Hidden static template used to generate high-speed data');

    nexttile;
    plot(Data_High.t(idxPlot) * 1000, Data_High.u_t(idxPlot), 'k-');
    xlabel('Time (ms)');
    ylabel('u(t) (mm)');
    title('True two-frequency vibration');
end

fprintf('[Step 1] Simulated hidden-gap high-speed data. Truth gap = %.3f mm\n', truth.g_high);

function [gapList, xCell, yCell] = load_stacked_gap_curves(filePath)
data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
x = data(:, 1);
y = data(:, 2);
breakIdx = [find(diff(x) < 0); numel(x)];
startIdx = [1; breakIdx(1:end-1) + 1];
nCurve = numel(breakIdx);
if nCurve == 7
    gapList = (0.2:0.2:1.4)';
else
    gapList = (1:nCurve)';
end
xCell = cell(nCurve, 1);
yCell = cell(nCurve, 1);
for i = 1:nCurve
    idx = startIdx(i):breakIdx(i);
    [xUnique, ia] = unique(x(idx), 'stable');
    yUnique = y(idx);
    xCell{i} = xUnique(:);
    yCell{i} = yUnique(ia);
end
end

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
        if isempty(idxUse), continue; end
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
