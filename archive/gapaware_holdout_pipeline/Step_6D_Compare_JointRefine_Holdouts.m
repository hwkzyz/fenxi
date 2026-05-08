%% Step 6D: compare raw-scan and final local joint refinement on key holdouts
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
cfgAna.f1Grid = 440:20:560;
cfgAna.f2Grid = 1180:20:1380;
cfgAna.numVarproCandidates = 4;
holdoutSweepList = unique([min(gapList), cfg.g_holdout, max(gapList)]).';

jointCfg = struct('gHalfWidth', 0.05, 'dxHalfWidth', 0.05, ...
    'fHalfWidth', 20, 'ampScaleLo', 0.50, 'ampScaleHi', 1.50);

rows = struct('holdout_gap_mm', {}, 'raw_gHat_mm', {}, 'raw_gap_error_mm', {}, ...
    'joint_gHat_mm', {}, 'joint_gap_error_mm', {}, ...
    'raw_rmse', {}, 'joint_rmse', {}, ...
    'raw_freq_error_Hz', {}, 'joint_freq_error_Hz', {}, ...
    'raw_amp_error_mm', {}, 'joint_amp_error_mm', {});

for ih = 1:numel(holdoutSweepList)
    gHoldout = holdoutSweepList(ih);
    gLow = pick_low_gap(gapList, gHoldout, cfg.g_low);
    cfgCase = cfgAna;
    cfgCase.g_holdout = gHoldout;
    cfgCase.g_low = gLow;

    templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, cfgCase.xGridN);
    idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
    F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

    rng(7);
    Data_High = simulate_rotating_waveform_from_template(F_high_true, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true);

    highMap = map_highspeed_to_space(Data_High, cfgCase.alpha_k, cfgCase.R_tip, ...
        templateLib.domain, cfgCase.fitActiveLevel);
    staticRaw = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
    rawResult = run_varpro_for_gap(highMap, templateLib, staticRaw.gHat, cfgCase, 'raw-scan');
    jointResult = run_joint_refinement(highMap, templateLib, staticRaw, rawResult, cfgCase, jointCfg);

    [rawFreqErr, rawAmpErr] = calc_param_errors(rawResult, cfgCase.f_true, cfgCase.A_true);
    [jointFreqErr, jointAmpErr] = calc_param_errors(jointResult, cfgCase.f_true, cfgCase.A_true);

    rows(end+1).holdout_gap_mm = gHoldout; %#ok<SAGROW>
    rows(end).raw_gHat_mm = staticRaw.gHat;
    rows(end).raw_gap_error_mm = abs(staticRaw.gHat - gHoldout);
    rows(end).joint_gHat_mm = jointResult.g_used;
    rows(end).joint_gap_error_mm = abs(jointResult.g_used - gHoldout);
    rows(end).raw_rmse = rawResult.rmse;
    rows(end).joint_rmse = jointResult.rmse;
    rows(end).raw_freq_error_Hz = rawFreqErr;
    rows(end).joint_freq_error_Hz = jointFreqErr;
    rows(end).raw_amp_error_mm = rawAmpErr;
    rows(end).joint_amp_error_mm = jointAmpErr;
end

summaryTable = struct2table(rows);
save(fullfile(outDir, 'stage6d_jointrefine_holdouts.mat'), 'summaryTable', '-v7.3');

figure('Name', 'Step 6D - Raw vs Joint-Refined Holdouts', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 10]);
tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = categorical(string(summaryTable.holdout_gap_mm));

nexttile;
bar(labels, [summaryTable.raw_gap_error_mm, summaryTable.joint_gap_error_mm], 'grouped');
ylabel('|gHat - gTruth| (mm)');
title('Gap estimation error');
legend({'raw-scan', 'joint-refined'}, 'Location', 'best', 'Box', 'off');

nexttile;
bar(labels, [summaryTable.raw_rmse, summaryTable.joint_rmse], 'grouped');
ylabel('RMSE');
title('Waveform fitting error');

nexttile;
bar(labels, [summaryTable.raw_freq_error_Hz, summaryTable.joint_freq_error_Hz], 'grouped');
ylabel('Mean frequency error (Hz)');
title('Frequency identification error');

nexttile;
bar(labels, [summaryTable.raw_amp_error_mm, summaryTable.joint_amp_error_mm], 'grouped');
ylabel('Mean amplitude error (mm)');
title('Amplitude identification error');

fprintf('[Step 6D] Raw vs joint-refined summary:\n');
disp(summaryTable);

function gLow = pick_low_gap(gapList, gHoldout, defaultLow)
available = gapList(abs(gapList - gHoldout) > 1e-12);
if any(abs(available - defaultLow) < 1e-12)
    gLow = defaultLow;
else
    [~, idx] = min(abs(available - defaultLow));
    gLow = available(idx);
end
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
trainIdx = abs(gapList - gHoldout) > 1e-12;
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

function result = run_joint_refinement(highMap, templateLib, staticState, gapAwareResult, cfg, jointCfg)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
p0 = gapAwareResult.p(:).';
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
resFun = @(theta) residual_calc_joint(theta, t, V, x, templateLib);
thetaOpt = refine_bounded_least_squares(resFun, theta0, lb, ub);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));
r = resFun(thetaOpt);
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
phi_raw = [thetaOpt(4), thetaOpt(7)];
result = struct('label', 'joint-refined', 'g_used', thetaOpt(1), ...
    'dx_used', thetaOpt(2), 'p', thetaOpt(3:end), 'theta', thetaOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)));
end

function r = residual_calc_joint(theta, t, V, x, templateLib)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
VSim = evaluate_template_surface(templateLib, gFit, x - dx0 - u);
r = V - VSim;
end

function [freqErr, ampErr] = calc_param_errors(result, fTrueIn, ATrueIn)
[fTrue, orderTrue] = sort(fTrueIn(:));
ATrue = ATrueIn(orderTrue);
[fFit, orderFit] = sort(result.f_id(:));
AFit = result.A_id(orderFit);
freqErr = mean(abs(fFit - fTrue));
ampErr = mean(abs(AFit - ATrue));
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
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, 'MaxIterations', 500);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1500, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
end

function r = residual_calc_6d_gapaware(p, t, V, x, templateLib, gFit)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
VSim = evaluate_template_surface(templateLib, gFit, x - u);
r = V - VSim;
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
