%% Step 1: Generate multi-revolution high-speed BTT waveform data
%  Uses COMSOL-simulated static gap curves as the physical template.
%  Simulates blade-tip-capacitance waveforms under dual-frequency vibration
%  at a HIDDEN gap — the algorithm never sees the template library.
%
%  Output: Data_High  — continuous-time waveform, OPR pulses, vibration truth
%          F_true     — the true static template (for validation only)

clc;
fprintf('=== Step 1: Generate high-speed multi-revolution data ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');

%% Load static gap curves from COMSOL data
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
fprintf('Loaded %d gap curves from: %s\n', numel(gapList), cfg.dataFile);
fprintf('Available gaps: %s mm\n', mat2str(gapList(:)', 3));

% Select the hidden-gap curve as the true physical template
idxHidden = find(abs(gapList - cfg.g_hidden) < 1e-12, 1);
if isempty(idxHidden)
    error('Hidden gap %.1f mm not found in available gaps.', cfg.g_hidden);
end
F_true = griddedInterpolant(xCell{idxHidden}, yCell{idxHidden}, 'pchip', 'nearest');
fprintf('True static template: gap = %.1f mm (HIDDEN from algorithm)\n', cfg.g_hidden);

%% Simulate multi-revolution high-speed waveform
rng(42);  % reproducible
Data_High = simulate_rotating_waveform(F_true, cfg, cfg.A_true, cfg.f_true, cfg.phi_true);

fprintf('Generated %.0f revolutions, %d samples, fs = %.0f Hz\n', ...
    cfg.NumRevs, numel(Data_High.t), cfg.fs);

save(fullfile(outDir, 'step1_data.mat'), 'Data_High', 'F_true', 'gapList', 'xCell', 'yCell', '-v7.3');
fprintf('[Step 1] Data saved to step1_data.mat\n');

%% Quick diagnostic plot
figure('Name', 'Step 1 - Generated Data', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 10]);
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Show a few revolutions
tShow = min(Data_High.t) + 4 * 60 / cfg.RPM_high + (0:0.5:2) * 60 / cfg.RPM_high;
idx = build_time_index(Data_High.t, tShow);

nexttile; hold on;
plot(Data_High.t(idx)*1000, Data_High.V_cap(idx), 'Color', [0.55 0.55 0.55]);
plot(Data_High.t(idx)*1000, Data_High.V_clean(idx), 'b-', 'LineWidth', 1.0);
ylabel('Capacitance');
legend({'noisy', 'clean'}, 'Location', 'best', 'Box', 'off');
title('High-speed waveform (segment)');

nexttile; hold on;
plot(Data_High.t(idx)*1000, Data_High.V_OPR(idx), 'r-');
ylabel('OPR');
title('Once-per-revolution pulses');

nexttile; hold on;
plot(Data_High.t(idx)*1000, Data_High.u_t(idx), 'k-');
xlabel('Time (ms)');
ylabel('u(t) (mm)');
title('True vibration displacement');

saveas(gcf, fullfile(outDir, 'fig_step1_data.png'));

%% ============ Local functions ======================================
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
        xCell{i} = xUnique(:);
        yCell{i} = y(idx(ia));
    end
end

function Data = simulate_rotating_waveform(Fx, cfg, A, f, phi)
    Omega = cfg.RPM_high * 2*pi / 60;
    T_rev = 2*pi / Omega;
    dt = 1 / cfg.fs;
    tEnd = (cfg.NumRevs + 1) * T_rev;
    t = (0:dt:tEnd)';

    % OPR timing (once per revolution, with a fixed delay offset)
    oprDelay = 0.08 * T_rev;
    T_opr = oprDelay + (0:cfg.NumRevs)' * T_rev;

    % Probe the static template to get its domain and baseline
    xProbe = linspace(Fx.GridVectors{1}(1), Fx.GridVectors{1}(end), 400)';
    yProbe = Fx(xProbe);
    domain = [min(xProbe), max(xProbe)];
    baseline = min(yProbe);
    xHalf = 0.52 * diff(domain);

    V_clean = baseline * ones(size(t));
    u_t = zeros(size(t));
    for im = 1:numel(A)
        u_t = u_t + A(im) * sin(2*pi*f(im)*t + phi(im));
    end

    % Generate waveform: for each revolution × sensor, place the template
    % at the physical position x_phys = x_nom - u(t)
    for m = 1:cfg.NumRevs
        for k = 1:numel(cfg.alpha_k)
            Te = T_opr(m) + cfg.alpha_k(k) / Omega;   % nominal arrival time
            idx = abs(t - Te) <= xHalf / (Omega * cfg.R_tip);
            xNom = Omega * cfg.R_tip * (t(idx) - Te);
            xPhys = xNom - u_t(idx);
            inDom = xPhys >= domain(1) & xPhys <= domain(2);
            idxAll = find(idx);
            idxUse = idxAll(inDom);
            V_clean(idxUse) = max(V_clean(idxUse), Fx(xPhys(inDom)));
        end
    end

    % Add noise
    noiseSigma = cfg.noiseRatio * range(yProbe);
    V_cap = V_clean + noiseSigma * randn(size(V_clean));

    % OPR signal
    V_OPR = zeros(size(t));
    oprWidth = 0.006 * T_rev;
    for m = 1:numel(T_opr)
        V_OPR(abs(t - T_opr(m)) <= oprWidth/2) = 5.0;
    end
    V_OPR = V_OPR + 0.02 * randn(size(V_OPR));

    Data = struct('t', t, 'V_cap', V_cap, 'V_clean', V_clean, ...
        'V_OPR', V_OPR, 'RPM', cfg.RPM_high, 'T_opr_truth', T_opr(:), ...
        'u_t', u_t, 'baseline', baseline, 'fs', cfg.fs);
end

function idx = build_time_index(tVec, tTarget)
    idx = false(size(tVec));
    for i = 1:2:numel(tTarget)-1
        idx = idx | (tVec >= tTarget(i) & tVec <= tTarget(i+1));
    end
end
