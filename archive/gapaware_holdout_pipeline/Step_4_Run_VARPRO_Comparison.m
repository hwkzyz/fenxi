%% Step 4: run fixed-template, gap-aware, and oracle VARPRO comparison
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'templateLib', 'templateLibFull');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage3_gap_estimation.mat'), 'highMap', 'staticState');

stepTimer = tic;
rawScanResult = run_varpro_for_gap(highMap, templateLib, staticState.gHat, cfg, 'raw-scan-initial');
gapAwareResult = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfg, cfg.jointRefine);
if cfg.computeStep4Baselines
    fixedResult = run_varpro_for_gap(highMap, templateLibFull, cfg.g_low, cfg, 'fixed-low-gap');
    oracleResult = run_varpro_for_gap(highMap, templateLibFull, cfg.g_holdout, cfg, 'oracle-holdout');
else
    fixedResult = [];
    oracleResult = [];
end
stepElapsed = toc(stepTimer);

Result = struct();
Result.truth = truth;
Result.staticState = staticState;
Result.fixed = fixedResult;
Result.rawScanInitial = rawScanResult;
Result.gapAware = gapAwareResult;
Result.oracle = oracleResult;
Result.config = cfg;
Result.step4ElapsedSeconds = stepElapsed;

summaryTable = make_summary_table(Result);
save(fullfile(outDir, 'stage4_identification_results.mat'), ...
    'Result', 'fixedResult', 'rawScanResult', 'gapAwareResult', 'oracleResult', 'summaryTable', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Step 4 - VARPRO Comparison', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 20, 11]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    labels = categorical(summaryTable.method);
    nexttile;
    bar(labels, summaryTable.rmse);
    ylabel('RMSE');
    title('Waveform fitting error');

    nexttile;
    bar(labels, summaryTable.mean_freq_error_Hz);
    ylabel('Mean frequency error (Hz)');
    title('Frequency identification error');

    nexttile;
    bar(labels, summaryTable.mean_amp_error_mm);
    ylabel('Mean amplitude error (mm)');
    title('Amplitude identification error');

    nexttile; hold on;
    idxPlot = 1:min(1200, numel(highMap.t_v));
    plot(highMap.t_v(idxPlot) * 1000, highMap.V_a(idxPlot), 'k.', 'MarkerSize', 4, 'DisplayName', 'observed');
    plot(highMap.t_v(idxPlot) * 1000, rawScanResult.VFit(idxPlot), 'c--', 'DisplayName', 'raw-scan initial');
    plot(highMap.t_v(idxPlot) * 1000, gapAwareResult.VFit(idxPlot), 'b-', 'DisplayName', 'gap-aware (joint refined)');
    if cfg.computeStep4Baselines
        plot(highMap.t_v(idxPlot) * 1000, fixedResult.VFit(idxPlot), 'r-', 'DisplayName', 'fixed low gap');
        plot(highMap.t_v(idxPlot) * 1000, oracleResult.VFit(idxPlot), 'g--', 'DisplayName', 'oracle');
    end
    xlabel('Time (ms)');
    ylabel('Voltage');
    title('Waveform fit comparison');
    legend('Location', 'best', 'Box', 'off');
end

fprintf('[Step 4] Identification summary:\n');
disp(summaryTable);
fprintf('[Step 4] VARPRO + joint refinement elapsed %.3f s\n', stepElapsed);

function result = run_varpro_for_gap(highMap, templateLib, gFit, cfg, label)
methodTimer = tic;
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
fixedEval = make_fixed_gap_evaluator(templateLib, gFit);
V0 = fixedEval.F(x);
K1 = fixedEval.dF(x);
DV = V - V0;
[pCandidates, coarseCosts] = coarse_varpro_2freq_candidates(t, DV, K1, ...
    cfg.f1Grid, cfg.f2Grid, cfg.numVarproCandidates);
lb = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid)];
ub = [0.800,  pi, max(cfg.f1Grid), 0.800,  pi, max(cfg.f2Grid)];
resFun = @(p) residual_calc_6d_fixed_gap(p, t, V, x, fixedEval);
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
    't', t, 'x', x, 'V', V, 'V0', V0, 'K1', K1, ...
    'DV', DV, 'VFit', VFit, 'uHat', uHat, ...
    'elapsedTimeSeconds', toc(methodTimer));
end

function result = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfg, jointCfg)
methodTimer = tic;
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
p0 = rawScanResult.p(:).';
theta0 = [staticState.gHat, staticState.dx0, p0];

gLo = max(0.05, staticState.gHat - jointCfg.gHalfWidth);
gHi = staticState.gHat + jointCfg.gHalfWidth;
dxLo = staticState.dx0 - jointCfg.dxHalfWidth;
dxHi = staticState.dx0 + jointCfg.dxHalfWidth;
a1 = max(0.001, jointCfg.ampScaleLo * abs(p0(1)));
a2 = max(0.001, jointCfg.ampScaleLo * abs(p0(4)));
b1 = min(0.8, max(a1 + 1e-6, jointCfg.ampScaleHi * abs(p0(1))));
b2 = min(0.8, max(a2 + 1e-6, jointCfg.ampScaleHi * abs(p0(4))));
f1Lo = max(min(cfg.f1Grid), p0(3) - jointCfg.fHalfWidth);
f1Hi = min(max(cfg.f1Grid), p0(3) + jointCfg.fHalfWidth);
f2Lo = max(min(cfg.f2Grid), p0(6) - jointCfg.fHalfWidth);
f2Hi = min(max(cfg.f2Grid), p0(6) + jointCfg.fHalfWidth);

lb = [gLo, dxLo, a1, -pi, f1Lo, a2, -pi, f2Lo];
ub = [gHi, dxHi, b1,  pi, f1Hi, b2,  pi, f2Hi];
templateEval = make_template_evaluator(templateLib);
resFun = @(theta) residual_calc_joint(theta, t, V, x, templateEval);
thetaOpt = refine_bounded_least_squares_joint(resFun, theta0, lb, ub);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));

r = resFun(thetaOpt);
uHat = displacement_2freq(thetaOpt(3:end), t);
VFit = V - r;
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
phi_raw = [thetaOpt(4), thetaOpt(7)];

result = struct('label', 'gap-aware-joint-refined', ...
    'g_used', thetaOpt(1), 'dx_used', thetaOpt(2), ...
    'p', thetaOpt(3:end), 'theta', thetaOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)), 't', t, 'x', x, 'V', V, 'VFit', VFit, 'uHat', uHat, ...
    'elapsedTimeSeconds', toc(methodTimer));
end

function [pCandidates, bestCosts] = coarse_varpro_2freq_candidates(t, DV, K1, f1Grid, f2Grid, numCandidates)
lambda = 1e-7 * max(numel(t), 1);
maxKeep = max(numCandidates * 4, numCandidates);
candCost = inf(maxKeep, 1);
candCoef = zeros(maxKeep, 4);
candFreq = zeros(maxKeep, 2);
f1Grid = f1Grid(:)';
f2Grid = f2Grid(:)';
S1 = sin(2*pi*t*f1Grid);
C1 = cos(2*pi*t*f1Grid);
S2 = sin(2*pi*t*f2Grid);
C2 = cos(2*pi*t*f2Grid);
KS1 = -K1 .* S1;
KC1 = -K1 .* C1;
KS2 = -K1 .* S2;
KC2 = -K1 .* C2;
for i1 = 1:numel(f1Grid)
    f1 = f1Grid(i1);
    h1s = KS1(:, i1);
    h1c = KC1(:, i1);
    for i2 = 1:numel(f2Grid)
        f2 = f2Grid(i2);
        if abs(f1 - f2) < 20, continue; end
        H = [h1s, h1c, KS2(:, i2), KC2(:, i2)];
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

function thetaOpt = refine_bounded_least_squares_joint(resFun, theta0, lb, ub)
theta0 = min(max(theta0, lb), ub);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, 'MaxIterations', 500);
    thetaOpt = lsqnonlin(resFun, theta0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((theta0 - lb + 1e-9) ./ max(ub - theta0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1500, 'TolX', 1e-9));
    thetaOpt = lb + width ./ (1 + exp(-zOpt));
end
end

function r = residual_calc_6d_fixed_gap(p, t, V, x, fixedEval)
u = displacement_2freq(p, t);
VSim = fixedEval.F(x - u);
r = V - VSim;
end

function r = residual_calc_joint(theta, t, V, x, templateEval)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = displacement_2freq(p, t);
VSim = templateEval.F(gFit, x - dx0 - u);
r = V - VSim;
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function summaryTable = make_summary_table(Result)
labels = strings(0, 1);
methods = {};
if ~isempty(Result.fixed)
    labels(end+1, 1) = "fixed_low_gap"; %#ok<AGROW>
    methods{end+1, 1} = Result.fixed; %#ok<AGROW>
end
labels(end+1, 1) = "gap_aware";
methods{end+1, 1} = Result.gapAware;
if ~isempty(Result.oracle)
    labels(end+1, 1) = "oracle_holdout"; %#ok<AGROW>
    methods{end+1, 1} = Result.oracle; %#ok<AGROW>
end
[fTrue, trueOrder] = sort(Result.truth.f);
ATrue = Result.truth.A(trueOrder);
n = numel(methods);
gUsed = zeros(n, 1); gapError = zeros(n, 1); f1 = zeros(n, 1); f2 = zeros(n, 1);
A1 = zeros(n, 1); A2 = zeros(n, 1); fErrMean = zeros(n, 1); AErrMean = zeros(n, 1); rmse = zeros(n, 1);
for i = 1:n
    m = methods{i};
    gUsed(i) = m.g_used;
    gapError(i) = abs(m.g_used - Result.truth.g_high);
    [fSort, methodOrder] = sort(m.f_id);
    ASort = m.A_id(methodOrder);
    f1(i) = fSort(1); f2(i) = fSort(2);
    A1(i) = ASort(1); A2(i) = ASort(2);
    fErrMean(i) = mean(abs(fSort - fTrue));
    AErrMean(i) = mean(abs(ASort - ATrue));
    rmse(i) = m.rmse;
end
summaryTable = table(labels, gUsed, gapError, f1, f2, A1, A2, ...
    fErrMean, AErrMean, rmse, ...
    'VariableNames', {'method', 'gap_used_mm', 'gap_error_mm', ...
    'f1_Hz', 'f2_Hz', 'A1_mm', 'A2_mm', ...
    'mean_freq_error_Hz', 'mean_amp_error_mm', 'rmse'});
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

function fixedEval = make_fixed_gap_evaluator(templateLib, g)
[curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g);
xGrid = templateLib.xGrid(:);
F = griddedInterpolant(xGrid, curveGrid(:), 'pchip', 'linear');
dF = griddedInterpolant(xGrid, dGrid(:), 'pchip', 'linear');
fixedEval = struct();
fixedEval.F = @(xq) F(xq);
fixedEval.dF = @(xq) dF(xq);
end

function templateEval = make_template_evaluator(templateLib)
zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
templateEval = struct();
templateEval.xGrid = templateLib.xGrid(:);
templateEval.zSort = zSort(:);
templateEval.SSort = templateLib.S(order, :);
templateEval.FxSort = templateLib.FxMat(order, :);
templateEval.F = @(g, xq) evaluate_template_from_cache(templateEval, g, xq);
end

function yPred = evaluate_template_from_cache(templateEval, g, xq)
zq = 1 / max(g, 1e-6);
curveGrid = interp1(templateEval.zSort, templateEval.SSort, zq, 'linear', 'extrap');
F = griddedInterpolant(templateEval.xGrid, curveGrid(:), 'pchip', 'linear');
yPred = F(xq);
end

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end
