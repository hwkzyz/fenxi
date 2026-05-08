%% Step 3B: experimental gap estimation using g + [b,a,c] joint search
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'templateLib');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'Data_High', 'truth');

highMap = map_highspeed_to_space(Data_High, cfg.alpha_k, cfg.R_tip, ...
    templateLib.domain, cfg.fitActiveLevel);
staticStateABC = estimate_highspeed_static_gap_abcjoint(highMap, templateLib, cfg);

save(fullfile(outDir, 'stage3b_abcjoint_gap_estimation.mat'), 'highMap', 'staticStateABC', '-v7.3');

figure('Name', 'Step 3B - ABC Joint Gap Search', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 10]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
sensorIds = unique(highMap.S_v(:)');
colors = lines(numel(sensorIds));
for ii = 1:numel(sensorIds)
    sid = sensorIds(ii);
    idx = highMap.S_v == sid;
    plot(highMap.x_v(idx), highMap.V_a(idx), '.', 'MarkerSize', 4, ...
        'Color', colors(ii, :), 'DisplayName', sprintf('sensor %d', sid));
end
xlabel('Mapped coordinate x (mm)');
ylabel('High-speed voltage');
title('Mapped samples');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(staticStateABC.gapQueryGrid, staticStateABC.J_g, 'b-', 'DisplayName', 'ABC joint cost');
xline(staticStateABC.gHat, 'g-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('gHat = %.3f', staticStateABC.gHat));
xline(truth.g_high, 'r--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('truth = %.3f', truth.g_high));
xlabel('Candidate gap g (mm)');
ylabel('Joint-search cost');
title(sprintf('Gap search with [b,a,c] (%s)', staticStateABC.supportType));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(staticStateABC.gapQueryGrid, staticStateABC.dxHatByGap, 'k-', 'DisplayName', 'dx(g)');
yline(staticStateABC.dx0, 'g--', 'LineWidth', 1.2, 'DisplayName', sprintf('dx0 = %.4f', staticStateABC.dx0));
xlabel('Candidate gap g (mm)');
ylabel('Recovered dx (mm)');
title('Recovered displacement nuisance term');
legend('Location', 'best', 'Box', 'off');

fprintf('[Step 3B] Mapped samples: %d\n', numel(highMap.t_v));
fprintf('[Step 3B] gHat = %.5f mm, truth = %.5f mm, error = %.5f mm\n', ...
    staticStateABC.gHat, truth.g_high, abs(staticStateABC.gHat - truth.g_high));
fprintf('[Step 3B] dx0 = %.6f mm, aHat = %.6f, bHat = %.6f, cHat = %.6f\n', ...
    staticStateABC.dx0, staticStateABC.aHat, staticStateABC.bHat, staticStateABC.cHat);

function highMap = map_highspeed_to_space(Data_High, alpha_k, R_tip, domain, activeLevel)
T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
xHalf = 0.48 * diff(domain);
base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + activeLevel * span;
tAll = [];
vAll = [];
xAll = [];
revAll = [];
sAll = [];
for m = 1:numel(T_opr)-1
    Omega = 2*pi / (T_opr(m+1) - T_opr(m));
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = find(abs(t - Te) <= xHalf / (Omega * R_tip));
        if isempty(idx), continue; end
        xLocal = Omega * R_tip * (t(idx) - Te);
        keep = xLocal >= domain(1) & xLocal <= domain(2) & V(idx) > activeThreshold;
        idx = idx(keep);
        xLocal = xLocal(keep);
        if isempty(idx), continue; end
        tAll = [tAll; t(idx)]; %#ok<AGROW>
        vAll = [vAll; V(idx)]; %#ok<AGROW>
        xAll = [xAll; xLocal]; %#ok<AGROW>
        revAll = [revAll; m * ones(numel(idx), 1)]; %#ok<AGROW>
        sAll = [sAll; k * ones(numel(idx), 1)]; %#ok<AGROW>
    end
end
highMap = struct('t_v', tAll, 'V_a', vAll, 'x_v', xAll, ...
    'rev_v', revAll, 'S_v', sAll, 'activeThreshold', activeThreshold);
end

function staticState = estimate_highspeed_static_gap_abcjoint(highMap, templateLib, cfg)
gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gapQueryGrid = linspace(gapMin, gapMax, cfg.gapQueryN);
y = highMap.V_a(:);
x = highMap.x_v(:);
J_g = zeros(numel(gapQueryGrid), 1);
thetaMat = zeros(numel(gapQueryGrid), 3);
dxHatByGap = zeros(numel(gapQueryGrid), 1);
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    [F0, Fx0] = evaluate_gap_curve_samples(templateLib, g, x);
    valid = isfinite(F0) & isfinite(Fx0);
    X = [ones(nnz(valid), 1), F0(valid), -Fx0(valid)];
    yv = y(valid);
    theta = X \ yv;
    r = yv - X * theta;
    J_g(ig) = mean(r.^2);
    thetaMat(ig, :) = theta(:)';
    if abs(theta(2)) > 1e-9
        dxHatByGap(ig) = theta(3) / theta(2);
    else
        dxHatByGap(ig) = 0;
    end
end
[~, bestIdx] = min(J_g);
gHat = gapQueryGrid(bestIdx);
bHat = thetaMat(bestIdx, 1);
aHat = thetaMat(bestIdx, 2);
cHat = thetaMat(bestIdx, 3);
dx0 = dxHatByGap(bestIdx);
staticState = struct('gHat', gHat, 'dx0', dx0, ...
    'bHat', bHat, 'aHat', aHat, 'cHat', cHat, ...
    'gapQueryGrid', gapQueryGrid(:), 'J_g', J_g(:), ...
    'thetaMat', thetaMat, 'dxHatByGap', dxHatByGap(:), ...
    'supportType', classify_gap_support(gHat, templateLib.gapTrain));
end

function supportType = classify_gap_support(g, gapTrain)
gMin = min(gapTrain);
gMax = max(gapTrain);
edgeBand = 0.08 * (gMax - gMin);
if g < gMin || g > gMax
    supportType = 'extrapolated';
elseif g <= gMin + edgeBand || g >= gMax - edgeBand
    supportType = 'boundary-supported';
else
    supportType = 'interior-supported';
end
end

function [curveVal, dVal] = evaluate_gap_curve_samples(templateLib, g, x)
[curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g);
curveVal = interp1(templateLib.xGrid, curveGrid(:), x, 'pchip', 'extrap');
dVal = interp1(templateLib.xGrid, dGrid(:), x, 'pchip', 'extrap');
end

function [curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g)
zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
SSort = templateLib.S(order, :);
FxSort = templateLib.FxMat(order, :);
zq = 1 / max(g, 1e-6);
curveGrid = interp1(zSort, SSort, zq, 'linear', 'extrap');
dGrid = interp1(zSort, FxSort, zq, 'linear', 'extrap');
end

function p = local_percentile(x, pct)
x = sort(x(isfinite(x)));
if isempty(x), p = NaN; return; end
q = 1 + (numel(x) - 1) * pct / 100;
lo = floor(q);
hi = ceil(q);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (q - lo) * (x(hi) - x(lo));
end
end
