%% Step 3: map high-speed data to space and estimate hidden static gap
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'templateLib');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'Data_High', 'truth');

stepTimer = tic;
highMap = map_highspeed_to_space(Data_High, cfg.alpha_k, cfg.R_tip, ...
    templateLib.domain, cfg.fitActiveLevel);
staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
staticState.elapsedTimeSeconds = toc(stepTimer);

save(fullfile(outDir, 'stage3_gap_estimation.mat'), 'highMap', 'staticState', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Step 3 - Mapping and Gap Estimation', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 20, 10]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

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
    plot(staticState.gapQueryGrid, staticState.J_static, 'b-', 'DisplayName', 'static cost');
    xline(staticState.gHat, 'g-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('gHat = %.3f', staticState.gHat));
    xline(truth.g_high, 'r--', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('truth = %.3f', truth.g_high));
    xlabel('Candidate gap g (mm)');
    ylabel('Static matching cost');
    title(sprintf('Blind static-gap estimation (%s)', staticState.supportType));
    legend('Location', 'best', 'Box', 'off');
end

fprintf('[Step 3] Mapped samples: %d\n', numel(highMap.t_v));
fprintf('[Step 3] gHat = %.5f mm, truth = %.5f mm, error = %.5f mm, support = %s\n', ...
    staticState.gHat, truth.g_high, abs(staticState.gHat - truth.g_high), staticState.supportType);
fprintf('[Step 3] Mapping + raw gap scan elapsed %.3f s\n', staticState.elapsedTimeSeconds);

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

function staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg)
gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gapQueryGrid = linspace(gapMin, gapMax, cfg.gapQueryN);
dxGrid = cfg.dxStaticGrid(:)';
y = highMap.V_a(:);
x = highMap.x_v(:);
J2D = zeros(numel(gapQueryGrid), numel(dxGrid));
xShiftGrid = x - dxGrid;
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    [curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
    yPredMat = interp1(templateLib.xGrid, curveGrid(:), xShiftGrid, 'pchip', 'extrap');
    resMat = y - yPredMat;
    J2D(ig, :) = mean(resMat.^2, 1, 'omitnan');
end
[J_static, dxIdx] = min(J2D, [], 2);
[~, bestIdx] = min(J_static);
gHat = gapQueryGrid(bestIdx);
staticState = struct('gHat', gHat, ...
    'dx0', dxGrid(dxIdx(bestIdx)), 'gapQueryGrid', gapQueryGrid(:), ...
    'dxGrid', dxGrid(:), 'J2D', J2D, 'J_static', J_static(:), ...
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
