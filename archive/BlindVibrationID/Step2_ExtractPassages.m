%% Step 2: Extract individual blade-passage waveforms from continuous data
%  Detects OPR pulses → computes passage centre times → extracts spatial
%  windows → resamples all passages to a common x-grid.
%
%  Output: Passages —  K waveforms on a common grid, with ground-truth
%                       alpha_k and d_k for validation.

clc;
fprintf('=== Step 2: Extract blade passages ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');
load(fullfile(outDir, 'step1_data.mat'), 'Data_High', 'F_true');

%% Detect OPR pulses
[T_opr_est] = detect_opr_pulses(Data_High.t, Data_High.V_OPR);
fprintf('Detected %d OPR pulses (expected %d)\n', numel(T_opr_est), cfg.NumRevs + 1);

% Use at most cfg.NumRevs complete revolutions
nRev = min(cfg.NumRevs, numel(T_opr_est) - 1);
T_opr_est = T_opr_est(1:nRev+1);

%% Compute mean rotation speed (allows slow drift)
Omega_rev = zeros(nRev, 1);
for m = 1:nRev
    Omega_rev(m) = 2*pi / (T_opr_est(m+1) - T_opr_est(m));
end
Omega_mean = mean(Omega_rev);
v0 = Omega_mean * cfg.R_tip;   % nominal tip speed (mm/s)
fprintf('Mean tip speed v0 = %.1f mm/s\n', v0);

%% Determine spatial window
% Probe the true template to find where it's "active"
xProbe = linspace(F_true.GridVectors{1}(1), F_true.GridVectors{1}(end), 400)';
yProbe = F_true(xProbe);
baseline = min(yProbe);
span = range(yProbe);
activeThreshold = baseline + cfg.activeLevel * span;
activeMask = yProbe > activeThreshold;
activeIdx = find(activeMask);
xActive = xProbe([activeIdx(1), activeIdx(end)]);
% Add margin
xMargin = 0.15 * diff(xActive);
xWin = [xActive(1) - xMargin, xActive(2) + xMargin];
fprintf('Spatial extraction window: [%.2f, %.2f] mm\n', xWin(1), xWin(2));

%% Common spatial grid
xGrid = linspace(xWin(1), xWin(2), cfg.xGridN)';
dx = xGrid(2) - xGrid(1);

%% Extract all passages:  (sensor, revolution) pairs
nSensors = numel(cfg.alpha_k);
nTotal = nRev * nSensors;

Y_raw     = zeros(nTotal, cfg.xGridN);
d_true    = zeros(nTotal, 1);
alpha_true = zeros(nTotal, 1);
t_centre  = zeros(nTotal, 1);
sensor_id = zeros(nTotal, 1);
rev_id    = zeros(nTotal, 1);

row = 0;
for m = 1:nRev
    Omega = Omega_rev(m);
    for k = 1:nSensors
        row = row + 1;
        Te = T_opr_est(m) + cfg.alpha_k(k) / Omega;  % nominal arrival time
        t_centre(row) = Te;

        % True vibration at passage centre
        u_true = eval_vibration(cfg.A_true, cfg.f_true, cfg.phi_true, Te);
        du_true = eval_vibration_deriv(cfg.A_true, cfg.f_true, cfg.phi_true, Te);
        d_true(row) = u_true;
        alpha_true(row) = 1 - du_true / v0;

        % Extract spatial waveform
        tLocal = Te + xGrid / v0;  % time → space mapping
        VLocal = interp1(Data_High.t, Data_High.V_cap, tLocal, 'pchip', NaN);
        Y_raw(row, :) = VLocal(:)';

        sensor_id(row) = k;
        rev_id(row) = m;
    end
end

%% Remove bad passages (NaN due to boundary)
goodRows = all(isfinite(Y_raw), 2);
fprintf('Removed %d passages with NaN (boundary clipping)\n', sum(~goodRows));
Y_raw     = Y_raw(goodRows, :);
d_true    = d_true(goodRows);
alpha_true = alpha_true(goodRows);
t_centre  = t_centre(goodRows);
sensor_id = sensor_id(goodRows);
rev_id    = rev_id(goodRows);
nPassages = sum(goodRows);

fprintf('Total usable passages: %d  (%d revs × %d sensors, minus clipped)\n', ...
    nPassages, nRev, nSensors);

%% Package output
Passages = struct();
Passages.x       = xGrid;
Passages.Y       = Y_raw;
Passages.d_true  = d_true;
Passages.alpha_true = alpha_true;
Passages.t_centre = t_centre;
Passages.sensor_id = sensor_id;
Passages.rev_id  = rev_id;
Passages.nPassages = nPassages;
Passages.v0      = v0;
Passages.Omega_mean = Omega_mean;

save(fullfile(outDir, 'step2_passages.mat'), 'Passages', '-v7.3');
fprintf('[Step 2] %d passages saved to step2_passages.mat\n', nPassages);

%% Diagnostic
figure('Name', 'Step 2 - Extracted Passages', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Show raw passages (first 30 overlaid)
nexttile; hold on;
nShow = min(30, nPassages);
cmap = lines(nShow);
for i = 1:nShow
    plot(xGrid, Y_raw(i, :), 'Color', [cmap(i, :), 0.5], 'LineWidth', 0.5);
end
xlabel('x (mm)'); ylabel('Capacitance');
title(sprintf('Raw passages (first %d)', nShow));

% Passage envelope (min, mean, max at each x)
nexttile; hold on;
Y_min = min(Y_raw, [], 1);
Y_max = max(Y_raw, [], 1);
Y_mean = mean(Y_raw, 1);
hFill = fill([xGrid; flipud(xGrid)], [Y_min(:); flipud(Y_max(:))], ...
    [0.7 0.8 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(xGrid, Y_mean, 'b-', 'LineWidth', 1.5);
xlabel('x (mm)'); ylabel('Capacitance');
title('Passage envelope (min–max + mean)');
legend([hFill], {'min–max range'}, 'Location', 'best', 'Box', 'off');

% Distribution of true alpha and d
nexttile;
histogram(alpha_true, 25, 'FaceColor', [0.2 0.6 0.6]);
xlabel('\alpha_k true'); ylabel('Count');
title(sprintf('True \\alpha_k: mean=%.4f, std=%.4f', mean(alpha_true), std(alpha_true)));

nexttile;
histogram(d_true, 25, 'FaceColor', [0.85 0.35 0.10]);
xlabel('d_k true (mm)'); ylabel('Count');
title(sprintf('True d_k: mean=%.4f, std=%.4f mm', mean(d_true), std(d_true)));

saveas(gcf, fullfile(outDir, 'fig_step2_passages.png'));

%% ============ Local functions ======================================
function T_opr = detect_opr_pulses(t, V_OPR)
    % Detect rising edges of OPR signal
    threshold = 2.5;
    above = V_OPR > threshold;
    edges = diff(above);
    riseIdx = find(edges == 1);
    % Refine: find peak near each rising edge
    T_opr = zeros(numel(riseIdx), 1);
    winHalf = 20;  % samples
    for i = 1:numel(riseIdx)
        lo = max(1, riseIdx(i) - winHalf);
        hi = min(numel(V_OPR), riseIdx(i) + winHalf);
        [~, imax] = max(V_OPR(lo:hi));
        T_opr(i) = t(lo + imax - 1);
    end
    T_opr = unique(T_opr);  % remove duplicates
end

function u = eval_vibration(A, f, phi, t)
    u = zeros(size(t));
    for m = 1:numel(A)
        u = u + A(m) * sin(2*pi*f(m)*t + phi(m));
    end
end

function du = eval_vibration_deriv(A, f, phi, t)
    du = zeros(size(t));
    for m = 1:numel(A)
        du = du + 2*pi*f(m)*A(m) * cos(2*pi*f(m)*t + phi(m));
    end
end
