%% Step 3: Blind template reconstruction + affine parameter estimation
%  Core algorithm: iterative affine registration and template averaging.
%
%  Model:  y_k(x) = T(alpha_k * x + d_k)   for each passage k
%
%  Alternating estimation:
%    1. Given T, solve (alpha_k, d_k) for each passage (affine registration)
%    2. Apply gauge conditions: mean(alpha)=1, mean(d)=0
%    3. Given (alpha_k, d_k), update T by inverse-transform averaging
%
%  NO template library, NO low-speed reference, NO gap calibration.

clc;
fprintf('=== Step 3: Blind affine template reconstruction ===\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
outDir = fullfile(scriptDir, 'results');
load(fullfile(outDir, 'cfg.mat'), 'cfg');
load(fullfile(outDir, 'step2_passages.mat'), 'Passages');

xGrid   = Passages.x(:);
Y       = Passages.Y;
nPassages = Passages.nPassages;
Nx      = numel(xGrid);

fprintf('Passages: %d,  grid points: %d\n', nPassages, Nx);

%% ---- 0. Initialization -------------------------------------------
% Baseline subtraction and coarse peak alignment
Y_base = min(Y, [], 2);
Y_cent = Y - Y_base;  % centred waveforms

% Find peak positions via 3-point quadratic interpolation
[~, peakIdx] = max(Y_cent, [], 2);
peakPos = zeros(nPassages, 1);
for k = 1:nPassages
    i0 = peakIdx(k);
    if i0 > 1 && i0 < Nx
        yVals = Y_cent(k, i0-1:i0+1);
        denom = 2*yVals(2) - yVals(1) - yVals(3);
        if abs(denom) > 1e-12
            frac = (yVals(3) - yVals(1)) / (2 * denom);
        else
            frac = 0;
        end
        peakPos(k) = xGrid(i0) + frac * (xGrid(2) - xGrid(1));
    else
        peakPos(k) = xGrid(i0);
    end
end

% Align by peak: each waveform is shifted so peak is at median peak position
peakRef = median(peakPos);
d_init = peakPos - peakRef;
alpha_init = ones(nPassages, 1);

% Build initial template: average peak-aligned waveforms
T_init = build_template_by_alignment(Y_cent, xGrid, alpha_init, d_init, cfg.blind.smoothLambda);

fprintf('Initial template built. Peak positions: mean=%.3f, std=%.3f mm\n', ...
    mean(peakPos), std(peakPos));

%% ---- Iterative reconstruction ------------------------------------
alpha_cur = alpha_init;
d_cur     = d_init;
T_cur     = T_init;
T_xGrid   = xGrid;  % template lives on the same grid

nIter = cfg.blind.maxIter;
history = struct('iter', {}, 'dAlpha', {}, 'dTemplate', {}, ...
    'meanAlpha', {}, 'stdAlpha', {}, 'meanD', {}, 'stdD', {});

fprintf('\nStarting blind reconstruction (max %d iterations)...\n', nIter);

for iter = 1:nIter
    %% (a) Affine registration: given T, estimate (alpha_k, d_k) for each passage
    alpha_new = zeros(nPassages, 1);
    d_new     = zeros(nPassages, 1);
    cost_reg  = zeros(nPassages, 1);

    % Use active region for registration (where waveform > 10% of span)
    activeMask = Y_cent > 0.1 * max(Y_cent, [], 2);

    % Process passages (can be parallelized with parfor)
    for k = 1:nPassages
        yk = Y(k, :)';
        [alpha_new(k), d_new(k), cost_reg(k)] = ...
            affine_register_passage(yk, xGrid, T_cur, T_xGrid, ...
            alpha_cur(k), d_cur(k), activeMask(k, :), cfg.blind, Nx);
    end

    % Remove outliers (relaxed in early iterations when template is poor)
    costThresh = median(cost_reg) + 5 * mad(cost_reg, 1);
    goodMask = cost_reg <= costThresh;
    nBad = sum(~goodMask);
    if nBad > 0 && nBad < nPassages * 0.4  % never remove more than 40%
        fprintf('  iter %2d: removed %d outlier passages (cost > %.3e)\n', iter, nBad, costThresh);
        alpha_new(~goodMask) = alpha_cur(~goodMask);
        d_new(~goodMask) = d_cur(~goodMask);
    elseif nBad >= nPassages * 0.4
        fprintf('  iter %2d: %d outliers (>40%%) — suppressing removal, template likely poor\n', iter, nBad);
    end

    %% (b) Gauge fixing: enforce mean(alpha)=1, mean(d)=0
    aGauge = mean(alpha_new);
    bGauge = mean(d_new);

    % Transform template:  T_gauge(x) = T_old(a * x - b)
    T_gauge = interp_waveform(T_xGrid, T_cur, aGauge * T_xGrid - bGauge);

    % Transform parameters:  alpha_new = alpha / a,  d_new = (d - b) / a
    alpha_gauge = alpha_new / aGauge;
    d_gauge     = (d_new - bGauge) / aGauge;

    %% (c) Template update: average inverse-transformed passages
    T_new = build_template_by_alignment(Y, xGrid, alpha_gauge, d_gauge, cfg.blind.smoothLambda);

    %% (d) Convergence check
    dAlpha = mean(abs(alpha_gauge - alpha_cur));
    dTemplate = max(abs(T_new - T_cur)) / (max(abs(T_cur)) + eps);

    history(iter).iter = iter;
    history(iter).dAlpha = dAlpha;
    history(iter).dTemplate = dTemplate;
    history(iter).meanAlpha = mean(alpha_gauge);
    history(iter).stdAlpha = std(alpha_gauge);
    history(iter).meanD = mean(d_gauge);
    history(iter).stdD = std(d_gauge);

    fprintf('  iter %2d:  dAlpha=%.2e  dTemplate=%.2e  std(alpha)=%.4f  std(d)=%.4f mm\n', ...
        iter, dAlpha, dTemplate, std(alpha_gauge), std(d_gauge));

    % Update
    alpha_cur = alpha_gauge;
    d_cur     = d_gauge;
    T_cur     = T_new;

    if dAlpha < cfg.blind.tolParam && dTemplate < cfg.blind.tolTemplate
        fprintf('  Converged at iteration %d.\n', iter);
        break;
    end
end

if iter == nIter
    fprintf('  Reached max iterations (%d).\n', nIter);
end

%% ---- Package results ---------------------------------------------
BlindResult = struct();
BlindResult.T         = T_cur;
BlindResult.T_xGrid   = T_xGrid;
BlindResult.alpha_est = alpha_cur;
BlindResult.d_est     = d_cur;
BlindResult.t_centre  = Passages.t_centre;
BlindResult.sensor_id = Passages.sensor_id;
BlindResult.rev_id    = Passages.rev_id;
BlindResult.v0        = Passages.v0;
BlindResult.xGrid     = xGrid;
BlindResult.nIter     = iter;
BlindResult.history   = history;
BlindResult.converged = (iter < nIter);

% Reconstructed vibration velocity
BlindResult.dot_d_est = Passages.v0 * (1 - BlindResult.alpha_est);  % alpha=1-dot_d/v0

save(fullfile(outDir, 'step3_blind_result.mat'), 'BlindResult', '-v7.3');
fprintf('[Step 3] Blind reconstruction saved.  Converged = %d\n', BlindResult.converged);

%% ---- Diagnostic plots --------------------------------------------
figure('Name', 'Step 3 - Blind Reconstruction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 14]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% (a) Template evolution
nexttile; hold on;
plot(xGrid, T_init, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, ...
    'DisplayName', 'initial (peak-aligned)');
plot(xGrid, T_cur, 'b-', 'LineWidth', 1.5, 'DisplayName', 'final (blind)');
xlabel('x (mm)'); ylabel('Capacitance');
title('(a) Reconstructed template');
legend('Location', 'best', 'Box', 'off');

% (b) Convergence metrics
nexttile; hold on;
yyaxis left;
plot(1:iter, [history.dAlpha], 'b-o', 'LineWidth', 1.2);
ylabel('\Delta\alpha');
yyaxis right;
plot(1:iter, [history.dTemplate], 'r-s', 'LineWidth', 1.2);
ylabel('\Delta Template');
xlabel('Iteration');
title('(b) Convergence');

% (c) Estimated alpha_k distribution
nexttile;
histogram(alpha_cur, 30, 'FaceColor', [0.2 0.6 0.6]);
xline(1, 'r--', 'LineWidth', 1.2);
xlabel('\alpha_k');
ylabel('Count');
title(sprintf('(c) Est \\alpha_k: std=%.4f', std(alpha_cur)));

% (d) Estimated d_k distribution
nexttile;
histogram(d_cur, 30, 'FaceColor', [0.85 0.35 0.10]);
xline(0, 'r--', 'LineWidth', 1.2);
xlabel('d_k (mm)');
ylabel('Count');
title(sprintf('(d) Est d_k: std=%.4f mm', std(d_cur)));

% (e) A few fitted passages vs raw
nexttile; hold on;
showK = [1, round(nPassages/3), round(2*nPassages/3), nPassages];
cmatAll = lines(numel(showK));
for i = 1:numel(showK)
    k = showK(i);
    yFit = interp_waveform(T_xGrid, T_cur, alpha_cur(k) * xGrid - d_cur(k));
    plot(xGrid, Y(k, :), '.', 'MarkerSize', 3, ...
        'Color', cmatAll(i, :), 'DisplayName', sprintf('k=%d raw', k));
    plot(xGrid, yFit, '-', 'LineWidth', 1.5, ...
        'Color', cmatAll(i, :)*0.6, 'DisplayName', sprintf('k=%d fit', k));
end
xlabel('x (mm)'); ylabel('Capacitance');
title('(e) Raw vs fitted passages');
legend('Location', 'best', 'Box', 'off');

% (f) Template overlay with all passages (inverse-transformed)
nexttile; hold on;
nShow = min(50, nPassages);
cmap = lines(nShow);
for k = 1:nShow
    xInv = (xGrid + d_cur(k)) / alpha_cur(k);
    yInv = interp1(xGrid, Y(k, :), xInv, 'pchip', NaN);
    plot(xGrid, yInv, 'Color', [cmap(k, :), 0.3], 'LineWidth', 0.3);
end
plot(xGrid, T_cur, 'k-', 'LineWidth', 2.0);
xlabel('x (mm)'); ylabel('Capacitance');
title('(f) Inverse-transformed passages + template');

saveas(gcf, fullfile(outDir, 'fig_step3_blind_recon.png'));

%% ============ Core functions ======================================
function T = build_template_by_alignment(Y, xGrid, alpha, d, smoothLambda)
    % Build template by averaging inverse-transformed passages
    % T(x) = mean_k [ y_k( (x - d_k) / alpha_k ) ]
    [K, Nx] = size(Y);
    T_acc = zeros(1, Nx);
    count = zeros(1, Nx);

    for k = 1:K
        xInv = (xGrid(:)' + d(k)) / alpha(k);
        % Interpolate y_k at xInv positions
        yInv = interp1(xGrid, Y(k, :), xInv, 'pchip', NaN);
        valid = isfinite(yInv);
        T_acc(valid) = T_acc(valid) + yInv(valid);
        count(valid) = count(valid) + 1;
    end

    T_raw = T_acc ./ max(count, 1);
    T_raw(count < K/4) = NaN;  % mark poorly supported regions
    T_raw = fillmissing(T_raw, 'linear', 'EndValues', 'nearest');

    % Light smoothing via penalized spline (ridge on 2nd difference)
    if smoothLambda > 0
        D2 = diff(eye(Nx), 2);
        penalty = smoothLambda * (D2' * D2);
        T = (eye(Nx) + penalty) \ T_raw(:);
    else
        T = T_raw(:);
    end
end

function [alpha_opt, d_opt, cost_opt] = affine_register_passage(...
        yk, xGrid, T_vals, T_xGrid, alpha0, d0, activeMask, blindCfg, Nx)
    % Affine registration: find (alpha, d) that best map template to y_k
    % Uses coarse grid search + fminsearch refinement

    % Active region (where waveform is above noise)
    activeIdx = find(activeMask);
    if numel(activeIdx) < 10
        % Fall back to using all points if active region too small
        activeIdx = (1:Nx)';
    end
    xAct = xGrid(activeIdx);
    yAct = yk(activeIdx);

    % Cost function for optimization
    costFun = @(p) compute_registration_cost(p(1), p(2), xAct, yAct, T_vals, T_xGrid);

    % --- Coarse grid search ---
    aGrid = linspace(blindCfg.alphaSearch(1), blindCfg.alphaSearch(2), blindCfg.alphaGridN);
    dGrid = linspace(blindCfg.dSearch(1), blindCfg.dSearch(2), blindCfg.dGridN);

    bestCost = inf;
    bestAlpha = alpha0;
    bestD = d0;

    for ia = 1:numel(aGrid)
        a = aGrid(ia);
        for id = 1:numel(dGrid)
            dVal = dGrid(id);
            c = compute_registration_cost(a, dVal, xAct, yAct, T_vals, T_xGrid);
            if c < bestCost
                bestCost = c;
                bestAlpha = a;
                bestD = dVal;
            end
        end
    end

    % --- Local refinement ---
    p0 = [bestAlpha, bestD];
    [pOpt, cost_opt] = fminsearch(costFun, p0, ...
        optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10, 'MaxIter', 200));

    alpha_opt = pOpt(1);
    d_opt = pOpt(2);

    % Sanity check: if alpha went crazy, fall back to grid result
    if alpha_opt < blindCfg.alphaSearch(1)*0.8 || alpha_opt > blindCfg.alphaSearch(2)*1.2
        alpha_opt = bestAlpha;
        d_opt = bestD;
        cost_opt = bestCost;
    end
end

function J = compute_registration_cost(alpha, d, x, y, T_vals, T_xGrid)
    % Compute || y(x) - T(alpha*x + d) ||^2
    xTrans = alpha * x - d;
    T_pred = interp_waveform(T_xGrid, T_vals, xTrans);
    valid = isfinite(T_pred);
    if sum(valid) < 5
        J = 1e10;
        return;
    end
    res = y(valid) - T_pred(valid);
    J = mean(res.^2);
end

function yq = interp_waveform(xGrid, yVals, xq)
    % Interpolate a waveform defined on xGrid at query points xq
    % Uses pchip, falls back to linear for speed if many queries
    yq = interp1(xGrid, yVals, xq, 'pchip', NaN);
end
