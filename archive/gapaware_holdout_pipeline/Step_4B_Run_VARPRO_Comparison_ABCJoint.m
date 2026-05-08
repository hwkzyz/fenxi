%% Step 4B: compare original raw-scan and experimental ABC-joint gap search
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'templateLib', 'templateLibFull');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage3_gap_estimation.mat'), 'highMap', 'staticState');
load(fullfile(outDir, 'stage3b_abcjoint_gap_estimation.mat'), 'staticStateABC');

fixedResult = run_varpro_for_gap(highMap, templateLibFull, cfg.g_low, cfg, 'fixed-low-gap');
rawGapAwareResult = run_varpro_for_gap(highMap, templateLib, staticState.gHat, cfg, 'raw-scan-gap');
abcJointResult = run_varpro_for_gap(highMap, templateLib, staticStateABC.gHat, cfg, 'abc-joint-gap');
oracleResult = run_varpro_for_gap(highMap, templateLibFull, cfg.g_holdout, cfg, 'oracle-holdout');

summaryTable = make_summary_table(truth, fixedResult, rawGapAwareResult, abcJointResult, oracleResult);
save(fullfile(outDir, 'stage4b_abcjoint_identification_results.mat'), ...
    'summaryTable', 'fixedResult', 'rawGapAwareResult', 'abcJointResult', 'oracleResult', ...
    'staticState', 'staticStateABC', '-v7.3');

figure('Name', 'Step 4B - ABC Joint Comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

labels = categorical(summaryTable.method);
nexttile;
bar(labels, summaryTable.gap_error_mm);
ylabel('Gap error (mm)');
title('Static-gap estimation error');

nexttile;
bar(labels, summaryTable.rmse);
ylabel('RMSE');
title('Waveform fitting error');

nexttile;
bar(labels, summaryTable.mean_freq_error_Hz);
ylabel('Mean frequency error (Hz)');
title('Frequency identification error');

nexttile; hold on;
idxPlot = 1:min(1200, numel(highMap.t_v));
plot(highMap.t_v(idxPlot) * 1000, highMap.V_a(idxPlot), 'k.', 'MarkerSize', 4, 'DisplayName', 'observed');
plot(highMap.t_v(idxPlot) * 1000, rawGapAwareResult.VFit(idxPlot), 'b-', 'DisplayName', sprintf('raw scan %.3f', staticState.gHat));
plot(highMap.t_v(idxPlot) * 1000, abcJointResult.VFit(idxPlot), 'm-', 'DisplayName', sprintf('ABC joint %.3f', staticStateABC.gHat));
plot(highMap.t_v(idxPlot) * 1000, oracleResult.VFit(idxPlot), 'g--', 'DisplayName', 'oracle');
xlabel('Time (ms)');
ylabel('Voltage');
title('Waveform-fit comparison');
legend('Location', 'best', 'Box', 'off');

fprintf('[Step 4B] Experimental ABC-joint comparison:\n');
disp(summaryTable);

function result = run_varpro_for_gap(highMap, templateLib, gFit, cfg, label)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
V0 = evaluate_template_surface(templateLib, gFit, x);
K1 = evaluate_template_derivative(templateLib, gFit, x);
DV = V - V0;
[pCandidates, coarseCosts] = coarse_varpro_2freq_candidates(t, DV, K1, ...
    cfg.f1Grid, cfg.f2Grid, cfg.numVarproCandidates);
lb = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid)];
ub = [0.800,  pi, max(cfg.f1Grid), 0.800,  pi, max(cfg.f2Grid)];
resFun = @(p) residual_calc_6d_gapaware(p, t, V, x, templateLib, gFit);
bestCost = inf;
pOpt = pCandidates(1, :);
for ic = 1:size(pCandidates, 1)
    pTry = refine_bounded_least_squares(resFun, pCandidates(ic, :), lb, ub);
    costTry = sum(resFun(pTry).^2);
    if costTry < bestCost
        bestCost = costTry;
        pOpt = pTry;
    end
end
r = resFun(pOpt);
uHat = displacement_2freq(pOpt, t);
VFit = V - r;
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
result = struct('label', label, 'g_used', gFit, 'p', pOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)), 'coarseCosts', coarseCosts, ...
    't', t, 'VFit', VFit, 'uHat', uHat);
end

function summaryTable = make_summary_table(truth, fixedResult, rawResult, abcResult, oracleResult)
labels = ["fixed_low_gap"; "raw_scan_gap"; "abc_joint_gap"; "oracle_holdout"];
methods = {fixedResult, rawResult, abcResult, oracleResult};
[fTrue, trueOrder] = sort(truth.f);
ATrue = truth.A(trueOrder);
n = numel(methods);
gUsed = zeros(n, 1);
gapError = zeros(n, 1);
fErrMean = zeros(n, 1);
AErrMean = zeros(n, 1);
rmse = zeros(n, 1);
for i = 1:n
    m = methods{i};
    gUsed(i) = m.g_used;
    gapError(i) = abs(m.g_used - truth.g_high);
    [fSort, order] = sort(m.f_id);
    ASort = m.A_id(order);
    fErrMean(i) = mean(abs(fSort - fTrue));
    AErrMean(i) = mean(abs(ASort - ATrue));
    rmse(i) = m.rmse;
end
summaryTable = table(labels, gUsed, gapError, fErrMean, AErrMean, rmse, ...
    'VariableNames', {'method', 'gap_used_mm', 'gap_error_mm', ...
    'mean_freq_error_Hz', 'mean_amp_error_mm', 'rmse'});
end

function [pCandidates, bestCosts] = coarse_varpro_2freq_candidates(t, DV, K1, f1Grid, f2Grid, numCandidates)
lambda = 1e-7 * max(numel(t), 1);
maxKeep = max(numCandidates * 4, numCandidates);
candCost = inf(maxKeep, 1);
candCoef = zeros(maxKeep, 4);
candFreq = zeros(maxKeep, 2);
for f1 = f1Grid
    s1 = sin(2*pi*f1*t);
    c1 = cos(2*pi*f1*t);
    for f2 = f2Grid
        if abs(f1 - f2) < 20, continue; end
        H = [-K1 .* s1, -K1 .* c1, ...
             -K1 .* sin(2*pi*f2*t), -K1 .* cos(2*pi*f2*t)];
        coef = (H' * H + lambda * eye(4)) \ (H' * DV);
        res = DV - H * coef;
        cost = sum(res.^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            candCost(worstIdx) = cost;
            candCoef(worstIdx, :) = coef(:)';
            candFreq(worstIdx, :) = [f1, f2];
        end
    end
end
[candCost, order] = sort(candCost, 'ascend');
candCoef = candCoef(order, :);
candFreq = candFreq(order, :);
keep = isfinite(candCost);
candCost = candCost(keep);
candCoef = candCoef(keep, :);
candFreq = candFreq(keep, :);
uniqueRows = true(size(candCost));
for i = 2:numel(candCost)
    prev = candFreq(1:i-1, :);
    if any(abs(prev(:, 1) - candFreq(i, 1)) < 8 & abs(prev(:, 2) - candFreq(i, 2)) < 8)
        uniqueRows(i) = false;
    end
end
candCost = candCost(uniqueRows);
candCoef = candCoef(uniqueRows, :);
candFreq = candFreq(uniqueRows, :);
nOut = min(numCandidates, numel(candCost));
pCandidates = zeros(nOut, 6);
for i = 1:nOut
    c = candCoef(i, :);
    pCandidates(i, :) = [hypot(c(1), c(2)), atan2(c(2), c(1)), candFreq(i, 1), ...
        hypot(c(3), c(4)), atan2(c(4), c(3)), candFreq(i, 2)];
end
bestCosts = candCost(1:nOut);
end

function pOpt = refine_bounded_least_squares(resFun, p0, lb, ub)
p0 = min(max(p0, lb), ub);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, 'MaxIterations', 400);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1000, 'TolX', 1e-8));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function r = residual_calc_6d_gapaware(p, t, V, x, templateLib, gFit)
u = displacement_2freq(p, t);
VSim = evaluate_template_surface(templateLib, gFit, x - u);
r = V - VSim;
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function yPred = evaluate_template_surface(templateLib, g, x)
[curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
yPred = interp1(templateLib.xGrid, curveGrid(:), x, 'pchip', 'extrap');
end

function dPred = evaluate_template_derivative(templateLib, g, x)
[~, dGrid] = evaluate_gap_curve_grid(templateLib, g);
dPred = interp1(templateLib.xGrid, dGrid(:), x, 'pchip', 'extrap');
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

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end
