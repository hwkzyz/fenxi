%% Step 6G targeted recheck for previously failed joint_wide_free scenarios
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
cfgAna = cfg;
cfgAna.fs = 1e5;
cfgAna.NumRevs_high = 6;
cfgAna.gapQueryN = 121;
cfgAna.dxStaticGrid = linspace(-0.6, 0.6, 81);
cfgAna.f1Grid = 440:10:640;
cfgAna.f2Grid = 1160:10:1400;
cfgAna.numVarproCandidates = 6;

jointWide = make_setting("joint_wide_free", true, true, 0.05, 0.05, 0.0, 0.0);

gHoldout = cfg.g_holdout;
gLow = cfg.g_low;
templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, cfgAna.xGridN);
templateLibFull = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);
idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

cases = {
    struct('scenario',"dual_default",'noise',0,'A',[0.25,0.15],'f',[500,1300],'phi',[pi/4,-pi/3]), ...
    struct('scenario',"dual_shifted_freq",'noise',0.02,'A',[0.25,0.15],'f',[460,1180],'phi',[pi/7,-pi/4])
};

rows = struct('scenario', {}, 'noise_ratio', {}, 'method', {}, ...
    'g_used_mm', {}, 'gap_error_mm', {}, ...
    'f1_error_Hz', {}, 'f2_error_Hz', {}, ...
    'A1_error_mm', {}, 'A2_error_mm', {}, 'rmse', {});

for ic = 1:numel(cases)
    c = cases{ic};
    cfgCase = cfgAna;
    cfgCase.noiseRatio = c.noise;
    cfgCase.A_true = c.A;
    cfgCase.f_true = c.f;
    cfgCase.phi_true = c.phi;

    rng(7);
    Data_High = simulate_rotating_waveform_from_template(F_high_true, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true);

    highMap = map_highspeed_to_space(Data_High, cfgCase.alpha_k, cfgCase.R_tip, ...
        templateLib.domain, cfgCase.fitActiveLevel);
    staticRaw = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
    staticCF = estimate_highspeed_static_gap_coarse_to_fine(highMap, templateLib, cfgCase);

    fixedResult = run_varpro_for_gap(highMap, templateLibFull, gLow, cfgCase, 'fixed');
    rawResult = run_varpro_for_gap(highMap, templateLib, staticRaw.gHat, cfgCase, 'raw');
    coarseFineResult = run_varpro_for_gap(highMap, templateLib, staticCF.gHat, cfgCase, 'cf');
    jointWideResult = run_joint_refinement_with_setting(highMap, templateLib, staticRaw, rawResult, cfgCase, jointWide);
    oracleResult = run_varpro_for_gap(highMap, templateLibFull, gHoldout, cfgCase, 'oracle');

    resultMap = {
        "fixed_low_gap", fixedResult;
        "raw_scan_only", rawResult;
        "coarse_to_fine_gap", coarseFineResult;
        "joint_wide_free", jointWideResult;
        "oracle_holdout", oracleResult
    };

    for im = 1:size(resultMap, 1)
        methodName = resultMap{im, 1};
        res = resultMap{im, 2};
        [~, ~, ~, ~, fErrVec, aErrVec] = calc_param_errors(res, cfgCase.f_true, cfgCase.A_true);
        rows(end+1).scenario = c.scenario; %#ok<SAGROW>
        rows(end).noise_ratio = c.noise;
        rows(end).method = methodName;
        rows(end).g_used_mm = res.g_used;
        rows(end).gap_error_mm = abs(res.g_used - gHoldout);
        rows(end).f1_error_Hz = fErrVec(1);
        rows(end).f2_error_Hz = fErrVec(2);
        rows(end).A1_error_mm = aErrVec(1);
        rows(end).A2_error_mm = aErrVec(2);
        rows(end).rmse = res.rmse;
    end
end

summaryTable = struct2table(rows);
save(fullfile(outDir, 'stage6g_targeted_recheck_failedcases.mat'), 'summaryTable', '-v7.3');

fprintf('[Step 6G targeted] Rechecked failed cases after fixing joint refinement bug:\n');
disp(summaryTable);

function sc = make_setting(name, useJoint, refineDx, gHalfWidth, dxHalfWidth, wG, wDx)
sc = struct();
sc.name = string(name);
sc.useJoint = useJoint;
sc.refineDx = refineDx;
sc.gHalfWidth = gHalfWidth;
sc.dxHalfWidth = dxHalfWidth;
sc.fHalfWidth = 20;
sc.ampScaleLo = 0.50;
sc.ampScaleHi = 1.50;
sc.weightG = wG;
sc.weightDx = wDx;
end

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

function templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN)
if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
trainIds = find(trainIdx);
G_train = gapList(trainIds);
xMin = -inf;
xMax = inf;
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(G_train), numel(xGrid));
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    S(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ii = 1:size(S, 1)
    FxMat(ii, :) = gradient(S(ii, :), dx);
end
templateLib = struct('gapTrain', G_train(:), 'xGrid', xGrid(:), ...
    'S', S, 'FxMat', FxMat, 'domain', [min(xGrid), max(xGrid)]);
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
        V_clean(idxUse) = Fx(xPhys(inDomain));
    end
end
signalRange = max(V_clean) - min(V_clean);
noiseStd = noiseRatio * signalRange;
Data = struct();
Data.t = t;
Data.V_cap = V_clean + noiseStd * randn(size(V_clean));
Data.V_clean = V_clean;
Data.T_opr_truth = T_opr;
Data.u_truth = u_t;
end

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
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    for id = 1:numel(dxGrid)
        dx = dxGrid(id);
        yPred = evaluate_template_surface(templateLib, g, x - dx);
        valid = isfinite(yPred);
        res = y(valid) - yPred(valid);
        J2D(ig, id) = mean(res.^2);
    end
end
[J_static, dxIdx] = min(J2D, [], 2);
[~, bestIdx] = min(J_static);
gHat = gapQueryGrid(bestIdx);
staticState = struct('gHat', gHat, 'dx0', dxGrid(dxIdx(bestIdx)));
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
        yPred = evaluate_template_surface(templateLib, g, x - dx);
        valid = isfinite(yPred);
        res = y(valid) - yPred(valid);
        J2D(ig, id) = mean(res.^2);
    end
end
[J_static, dxIdx] = min(J2D, [], 2);
[~, bestIdx] = min(J_static);
gCoarse = gapQueryGrid(bestIdx);
dxCoarse = dxGrid(dxIdx(bestIdx));
lb = [gapMin, min(dxGrid)];
ub = [gapMax, max(dxGrid)];
theta0 = inv_sigmoid_transform([gCoarse, dxCoarse], lb, ub);
obj = @(theta) objective_in_unconstrained(theta, lb, ub, templateLib, x, y);
thetaOpt = fminsearch(obj, theta0, optimset('Display', 'off', ...
    'MaxIter', 1200, 'MaxFunEvals', 2000, 'TolX', 1e-9, 'TolFun', 1e-10));
paramsOpt = forward_sigmoid_transform(thetaOpt, lb, ub);
staticState = struct('gHat', paramsOpt(1), 'dx0', paramsOpt(2));
end

function J = objective_in_unconstrained(theta, lb, ub, templateLib, x, y)
params = forward_sigmoid_transform(theta, lb, ub);
J = evaluate_joint_cost(templateLib, x, y, params(1), params(2));
end

function params = forward_sigmoid_transform(theta, lb, ub)
span = ub - lb;
params = lb + span ./ (1 + exp(-theta));
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

function result = run_varpro_for_gap(highMap, templateLib, gFit, cfg, label)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
V0 = evaluate_template_surface(templateLib, gFit, x);
K1 = evaluate_template_derivative(templateLib, gFit, x);
DV = V - V0;
[pCandidates, ~] = coarse_varpro_2freq_candidates(t, DV, K1, ...
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
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
result = struct('label', label, 'g_used', gFit, 'p', pOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)));
end

function result = run_joint_refinement_with_setting(highMap, templateLib, staticState, rawResult, cfg, setting)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
p0 = rawResult.p(:).';
theta0 = [staticState.gHat, staticState.dx0, p0];
gLo = max(0.05, staticState.gHat - setting.gHalfWidth);
gHi = staticState.gHat + setting.gHalfWidth;
dxLo = staticState.dx0 - setting.dxHalfWidth;
dxHi = staticState.dx0 + setting.dxHalfWidth;
a1 = max(0.001, setting.ampScaleLo * abs(p0(1)));
a2 = max(0.001, setting.ampScaleLo * abs(p0(4)));
b1 = min(0.8, max(a1 + 1e-6, setting.ampScaleHi * abs(p0(1))));
b2 = min(0.8, max(a2 + 1e-6, setting.ampScaleHi * abs(p0(4))));
f1Lo = max(min(cfg.f1Grid), p0(3) - setting.fHalfWidth);
f1Hi = min(max(cfg.f1Grid), p0(3) + setting.fHalfWidth);
f2Lo = max(min(cfg.f2Grid), p0(6) - setting.fHalfWidth);
f2Hi = min(max(cfg.f2Grid), p0(6) + setting.fHalfWidth);
lb = [gLo, dxLo, a1, -pi, f1Lo, a2, -pi, f2Lo];
ub = [gHi, dxHi, b1,  pi, f1Hi, b2,  pi, f2Hi];
resFun = @(theta) residual_calc_joint(theta, t, V, x, templateLib);
thetaOpt = refine_bounded_least_squares(resFun, theta0, lb, ub);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));
r = resFun(thetaOpt);
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
phi_raw = [thetaOpt(4), thetaOpt(7)];
result = struct('label', 'joint-refined', 'g_used', thetaOpt(1), 'dx_used', thetaOpt(2), ...
    'p', thetaOpt(3:end), 'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)));
end

function r = residual_calc_joint(theta, t, V, x, templateLib)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = displacement_2freq(p, t);
VSim = evaluate_template_surface(templateLib, gFit, x - dx0 - u);
r = V - VSim;
end

function [freqErr, ampErr, fSort, ASort, fErrVec, aErrVec] = calc_param_errors(result, fTrueIn, ATrueIn)
[fTrue, orderTrue] = sort(fTrueIn(:));
ATrue = ATrueIn(orderTrue);
[fSort, orderFit] = sort(result.f_id(:));
ASort = result.A_id(orderFit);
fErrVec = abs(fSort - fTrue);
aErrVec = abs(ASort - ATrue);
freqErr = mean(fErrVec);
ampErr = mean(aErrVec);
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
