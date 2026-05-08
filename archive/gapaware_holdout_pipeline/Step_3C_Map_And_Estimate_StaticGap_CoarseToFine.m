%% Step 3C: experimental coarse-to-fine raw-scan gap estimation
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
staticStateCF = estimate_highspeed_static_gap_coarse_to_fine(highMap, templateLib, cfg);

save(fullfile(outDir, 'stage3c_coarsetofine_gap_estimation.mat'), ...
    'highMap', 'staticStateCF', '-v7.3');

figure('Name', 'Step 3C - Coarse-to-Fine Gap Estimation', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 10]);
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
title('High-speed samples after spatial mapping');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(staticStateCF.gapQueryGrid, staticStateCF.J_static, 'b-', 'DisplayName', 'coarse cost');
xline(staticStateCF.gHat, 'g-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('refined gHat = %.3f', staticStateCF.gHat));
xline(staticStateCF.gCoarse, 'c--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('coarse g = %.3f', staticStateCF.gCoarse));
xline(truth.g_high, 'r--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('truth = %.3f', truth.g_high));
xlabel('Candidate gap g (mm)');
ylabel('Static matching cost');
title(sprintf('Coarse-to-fine gap search (%s)', staticStateCF.supportType));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
imagesc(staticStateCF.dxGrid, staticStateCF.gapQueryGrid, staticStateCF.J2D);
set(gca, 'YDir', 'normal');
plot(staticStateCF.dxCoarse, staticStateCF.gCoarse, 'wo', 'MarkerSize', 7, ...
    'LineWidth', 1.2, 'DisplayName', 'coarse best');
plot(staticStateCF.dx0, staticStateCF.gHat, 'gx', 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'DisplayName', 'refined best');
xlabel('Candidate dx (mm)');
ylabel('Candidate gap g (mm)');
title('Coarse joint-search surface');
colorbar;

fprintf('[Step 3C] Mapped samples: %d\n', numel(highMap.t_v));
fprintf('[Step 3C] coarse  g = %.5f mm, dx = %.5f mm\n', ...
    staticStateCF.gCoarse, staticStateCF.dxCoarse);
fprintf('[Step 3C] refined gHat = %.5f mm, truth = %.5f mm, error = %.5f mm\n', ...
    staticStateCF.gHat, truth.g_high, abs(staticStateCF.gHat - truth.g_high));
fprintf('[Step 3C] refined dx0 = %.5f mm, support = %s\n', ...
    staticStateCF.dx0, staticStateCF.supportType);

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

function staticState = estimate_highspeed_static_gap_coarse_to_fine(highMap, templateLib, cfg)
gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gapQueryGrid = linspace(gapMin, gapMax, cfg.gapQueryN);
dxGrid = cfg.dxStaticGrid(:)';
y = highMap.V_a(:);
x = highMap.x_v(:);
J2D = zeros(numel(gapQueryGrid), numel(dxGrid));
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    for id = 1:numel(dxGrid)
        dx = dxGrid(id);
        J2D(ig, id) = evaluate_joint_cost(templateLib, x, y, g, dx);
    end
end
[J_static, dxIdx] = min(J2D, [], 2);
[coarseCost, bestIdx] = min(J_static);
gCoarse = gapQueryGrid(bestIdx);
dxCoarse = dxGrid(dxIdx(bestIdx));

lb = [gapMin, min(dxGrid)];
ub = [gapMax, max(dxGrid)];
theta0 = inv_sigmoid_transform([gCoarse, dxCoarse], lb, ub);
obj = @(theta) objective_in_unconstrained(theta, lb, ub, templateLib, x, y);
thetaOpt = fminsearch(obj, theta0, optimset('Display', 'off', ...
    'MaxIter', 1200, 'MaxFunEvals', 2000, 'TolX', 1e-9, 'TolFun', 1e-10));
[paramsOpt, ~] = forward_sigmoid_transform(thetaOpt, lb, ub);
gHat = paramsOpt(1);
dx0 = paramsOpt(2);
refinedCost = evaluate_joint_cost(templateLib, x, y, gHat, dx0);

staticState = struct('gHat', gHat, 'dx0', dx0, ...
    'gCoarse', gCoarse, 'dxCoarse', dxCoarse, ...
    'coarseCost', coarseCost, 'refinedCost', refinedCost, ...
    'gapQueryGrid', gapQueryGrid(:), 'dxGrid', dxGrid(:), ...
    'J2D', J2D, 'J_static', J_static(:), ...
    'supportType', classify_gap_support(gHat, templateLib.gapTrain));
end

function J = objective_in_unconstrained(theta, lb, ub, templateLib, x, y)
[params, ~] = forward_sigmoid_transform(theta, lb, ub);
J = evaluate_joint_cost(templateLib, x, y, params(1), params(2));
end

function [params, placeholder] = forward_sigmoid_transform(theta, lb, ub)
span = ub - lb;
params = lb + span ./ (1 + exp(-theta));
placeholder = NaN;
end

function theta = inv_sigmoid_transform(params, lb, ub)
span = max(ub - lb, 1e-12);
ratio = (params - lb) ./ span;
ratio = min(max(ratio, 1e-9), 1 - 1e-9);
theta = log(ratio ./ (1 - ratio));
end

function J = evaluate_joint_cost(templateLib, x, y, g, dx)
yPred = evaluate_template_surface(templateLib, g, x - dx);
valid = isfinite(yPred);
res = y(valid) - yPred(valid);
J = mean(res.^2);
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

function yPred = evaluate_template_surface(templateLib, g, x)
[curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
yPred = interp1(templateLib.xGrid, curveGrid(:), x, 'pchip', 'extrap');
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
