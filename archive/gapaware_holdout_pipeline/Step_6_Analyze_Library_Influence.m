%% Step 6: analyze holdout position and library sparsity for the joint-refined main pipeline
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
maxCombosPerSize = 6;
holdoutSweepList = unique([min(gapList), cfg.g_holdout, max(gapList)]).';

%% Part A: holdout position sweep
holdoutRows = struct('holdout_gap_mm', {}, 'low_gap_mm', {}, 'is_edge_holdout', {}, ...
    'gHat_raw_mm', {}, 'g_used_final_mm', {}, 'gap_error_mm', {}, 'supportType', {}, ...
    'fixed_freq_error_Hz', {}, 'gapaware_freq_error_Hz', {}, ...
    'fixed_amp_error_mm', {}, 'gapaware_amp_error_mm', {}, ...
    'fixed_rmse', {}, 'gapaware_rmse', {});

for ih = 1:numel(holdoutSweepList)
    gHoldout = holdoutSweepList(ih);
    gLow = pick_low_gap(gapList, gHoldout, cfg.g_low);
    cfgCase = cfgAna;
    cfgCase.g_holdout = gHoldout;
    cfgCase.g_low = gLow;

    templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, cfgCase.xGridN);
    templateLibFull = build_gap_template_library(gapList, xCell, yCell, NaN, cfgCase.xGridN);

    idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
    F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

    rng(7);
    Data_High = simulate_rotating_waveform_from_template(F_high_true, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true);

    highMap = map_highspeed_to_space(Data_High, cfgCase.alpha_k, cfgCase.R_tip, ...
        templateLib.domain, cfgCase.fitActiveLevel);
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);

    fixedResult = run_varpro_for_gap(highMap, templateLibFull, gLow, cfgCase, 'fixed-low-gap');
    rawScanResult = run_varpro_for_gap(highMap, templateLib, staticState.gHat, cfgCase, 'raw-scan-initial');
    gapAwareResult = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfgCase, cfgCase.jointRefine);

    [freqErrFixed, ampErrFixed] = calc_param_errors(fixedResult, cfgCase.f_true, cfgCase.A_true);
    [freqErrGap, ampErrGap] = calc_param_errors(gapAwareResult, cfgCase.f_true, cfgCase.A_true);

    holdoutRows(end+1).holdout_gap_mm = gHoldout; %#ok<SAGROW>
    holdoutRows(end).low_gap_mm = gLow;
    holdoutRows(end).is_edge_holdout = gHoldout == min(gapList) || gHoldout == max(gapList);
    holdoutRows(end).gHat_raw_mm = staticState.gHat;
    holdoutRows(end).g_used_final_mm = gapAwareResult.g_used;
    holdoutRows(end).gap_error_mm = abs(gapAwareResult.g_used - gHoldout);
    holdoutRows(end).supportType = string(staticState.supportType);
    holdoutRows(end).fixed_freq_error_Hz = freqErrFixed;
    holdoutRows(end).gapaware_freq_error_Hz = freqErrGap;
    holdoutRows(end).fixed_amp_error_mm = ampErrFixed;
    holdoutRows(end).gapaware_amp_error_mm = ampErrGap;
    holdoutRows(end).fixed_rmse = fixedResult.rmse;
    holdoutRows(end).gapaware_rmse = gapAwareResult.rmse;
end

holdoutTable = struct2table(holdoutRows);

%% Part B: library sparsity sweep around the default interior holdout
gHoldout = cfg.g_holdout;
availableTrain = gapList(abs(gapList - gHoldout) > 1e-12);
idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

rng(7);
templateLibDense = build_gap_template_library(gapList, xCell, yCell, gHoldout, cfgAna.xGridN);
Data_High_dense = simulate_rotating_waveform_from_template(F_high_true, cfgAna.RPM_high, ...
    cfgAna.NumRevs_high, cfgAna.fs, cfgAna.R_tip, cfgAna.alpha_k, cfgAna.noiseRatio, ...
    templateLibDense.domain, cfgAna.A_true, cfgAna.f_true, cfgAna.phi_true);

sparseRows = struct('subset_size', {}, 'subset_gaps', {}, 'brackets_holdout', {}, ...
    'gHat_raw_mm', {}, 'g_used_final_mm', {}, 'gap_error_mm', {}, 'supportType', {}, ...
    'freq_error_Hz', {}, 'amp_error_mm', {}, 'rmse', {});

subsetSizesToTest = intersect([2, 3, 4, numel(availableTrain)], 2:numel(availableTrain));
for k = subsetSizesToTest
    combos = nchoosek(1:numel(availableTrain), k);
    if size(combos, 1) > maxCombosPerSize
        keepIdx = unique(round(linspace(1, size(combos, 1), maxCombosPerSize)));
        combos = combos(keepIdx, :);
    end
    for ic = 1:size(combos, 1)
        subsetGaps = availableTrain(combos(ic, :));
        trainMask = ismember(gapList, subsetGaps);
        templateLibSparse = build_gap_template_library_ids(gapList, xCell, yCell, trainMask, cfgAna.xGridN);
        highMap = map_highspeed_to_space(Data_High_dense, cfgAna.alpha_k, cfgAna.R_tip, ...
            templateLibSparse.domain, cfgAna.fitActiveLevel);
        staticState = estimate_highspeed_static_gap_raw(highMap, templateLibSparse, cfgAna);
        rawScanResult = run_varpro_for_gap(highMap, templateLibSparse, staticState.gHat, cfgAna, 'raw-scan-initial');
        gapAwareResult = run_joint_refinement(highMap, templateLibSparse, staticState, rawScanResult, cfgAna, cfgAna.jointRefine);
        [freqErrGap, ampErrGap] = calc_param_errors(gapAwareResult, cfgAna.f_true, cfgAna.A_true);

        sparseRows(end+1).subset_size = k; %#ok<SAGROW>
        sparseRows(end).subset_gaps = join(string(subsetGaps.'), ',');
        sparseRows(end).brackets_holdout = any(subsetGaps < gHoldout) && any(subsetGaps > gHoldout);
        sparseRows(end).gHat_raw_mm = staticState.gHat;
        sparseRows(end).g_used_final_mm = gapAwareResult.g_used;
        sparseRows(end).gap_error_mm = abs(gapAwareResult.g_used - gHoldout);
        sparseRows(end).supportType = string(staticState.supportType);
        sparseRows(end).freq_error_Hz = freqErrGap;
        sparseRows(end).amp_error_mm = ampErrGap;
        sparseRows(end).rmse = gapAwareResult.rmse;
    end
end

sparseTable = struct2table(sparseRows);

save(fullfile(outDir, 'stage6_library_influence.mat'), 'holdoutTable', 'sparseTable', '-v7.3');

%% Visualization
figure('Name', 'Step 6 - Holdout Position Influence', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

edgeMask = holdoutTable.is_edge_holdout;

nexttile; hold on;
plot(holdoutTable.holdout_gap_mm(~edgeMask), holdoutTable.gap_error_mm(~edgeMask), 'bo-', ...
    'DisplayName', 'interior holdout');
plot(holdoutTable.holdout_gap_mm(edgeMask), holdoutTable.gap_error_mm(edgeMask), 'rs', ...
    'MarkerFaceColor', 'r', 'DisplayName', 'edge holdout');
xlabel('Held-out gap g (mm)');
ylabel('|gHat - gTruth| (mm)');
title('Static-gap error vs holdout position');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(holdoutTable.holdout_gap_mm, holdoutTable.fixed_freq_error_Hz, 'r--o', ...
    'DisplayName', 'fixed low-gap');
plot(holdoutTable.holdout_gap_mm, holdoutTable.gapaware_freq_error_Hz, 'b-o', ...
    'DisplayName', 'gap-aware');
xlabel('Held-out gap g (mm)');
ylabel('Mean frequency error (Hz)');
title('Frequency error vs holdout position');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(holdoutTable.holdout_gap_mm, holdoutTable.fixed_rmse, 'r--o', ...
    'DisplayName', 'fixed low-gap');
plot(holdoutTable.holdout_gap_mm, holdoutTable.gapaware_rmse, 'b-o', ...
    'DisplayName', 'gap-aware');
xlabel('Held-out gap g (mm)');
ylabel('RMSE');
title('Waveform fitting error vs holdout position');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
ratio = holdoutTable.fixed_rmse ./ max(holdoutTable.gapaware_rmse, eps);
bar(categorical(string(holdoutTable.holdout_gap_mm)), ratio, 'FaceColor', [0.25 0.55 0.85]);
xlabel('Held-out gap g (mm)');
ylabel('RMSE(fixed) / RMSE(gap-aware)');
title('Gap-aware improvement factor');

figure('Name', 'Step 6 - Library Sparsity Influence', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot_subset_scatter(sparseTable.subset_size, sparseTable.gap_error_mm, sparseTable.brackets_holdout);
xlabel('Training subset size');
ylabel('|gHat - gTruth| (mm)');
title('Gap error under library sparsity');

nexttile; hold on;
plot_subset_scatter(sparseTable.subset_size, sparseTable.freq_error_Hz, sparseTable.brackets_holdout);
xlabel('Training subset size');
ylabel('Mean frequency error (Hz)');
title('Frequency error under library sparsity');

nexttile; hold on;
plot_subset_scatter(sparseTable.subset_size, sparseTable.rmse, sparseTable.brackets_holdout);
xlabel('Training subset size');
ylabel('RMSE');
title('Waveform fitting error under library sparsity');

nexttile; hold on;
subsetSizes = unique(sparseTable.subset_size);
medVals = zeros(size(subsetSizes));
minVals = zeros(size(subsetSizes));
maxVals = zeros(size(subsetSizes));
for i = 1:numel(subsetSizes)
    idx = sparseTable.subset_size == subsetSizes(i);
    vals = sparseTable.freq_error_Hz(idx);
    medVals(i) = median(vals);
    minVals(i) = min(vals);
    maxVals(i) = max(vals);
end
errorbar(subsetSizes, medVals, medVals-minVals, maxVals-medVals, 'o-', 'Color', [0.10 0.35 0.75], ...
    'MarkerFaceColor', [0.10 0.35 0.75]);
xlabel('Training subset size');
ylabel('Frequency error range (Hz)');
title('Median / min / max frequency error');

%% Console summary
edgeSummary = summarize_group(holdoutTable(edgeMask, :), "edge");
interiorSummary = summarize_group(holdoutTable(~edgeMask, :), "interior");
subsetSummary = summarize_subset_sizes(sparseTable);
bracketSummary = summarize_bracketing(sparseTable);

fprintf('\n[Step 6] Holdout position summary:\n');
disp(holdoutTable);
fprintf('\n[Step 6] Edge vs interior summary:\n');
disp(edgeSummary);
disp(interiorSummary);
fprintf('\n[Step 6] Library sparsity summary by subset size:\n');
disp(subsetSummary);
fprintf('\n[Step 6] Bracketing effect summary:\n');
disp(bracketSummary);

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
if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
templateLib = build_gap_template_library_ids(gapList, xCell, yCell, trainIdx, xGridN);
end

function templateLib = build_gap_template_library_ids(gapList, xCell, yCell, trainMask, xGridN)
trainIds = find(trainMask);
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
templateLib = struct();
templateLib.gapTrain = G_train(:);
templateLib.xGrid = xGrid(:);
templateLib.S = S;
templateLib.FxMat = FxMat;
templateLib.domain = [min(xGrid), max(xGrid)];
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
        V_clean(idxUse) = max(V_clean(idxUse), Fx(xPhys(inDomain)));
    end
end
noiseSigma = noiseRatio * max(range(yProbe), eps);
V_cap = V_clean + noiseSigma * randn(size(V_clean));
V_OPR = zeros(size(t));
oprWidth = 0.006 * T;
for m = 1:numel(T_opr)
    V_OPR(abs(t - T_opr(m)) <= oprWidth / 2) = 5.0;
end
V_OPR = V_OPR + 0.02 * randn(size(V_OPR));
Data = struct('t', t, 'V_cap', V_cap, 'V_cap_clean', V_clean, ...
    'V_OPR', V_OPR, 'RPM', RPM, 'T_opr_truth', T_opr(:), ...
    'u_t', u_t, 'baseline', baseline);
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
[xProf, yProf, nProf] = build_raw_profile(x, y, templateLib.domain, cfg);
J2D = zeros(numel(gapQueryGrid), numel(dxGrid));
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    for id = 1:numel(dxGrid)
        dx = dxGrid(id);
        yPred = evaluate_template_surface(templateLib, g, xProf - dx);
        dPred = evaluate_template_derivative(templateLib, g, xProf - dx);
        J2D(ig, id) = robust_profile_cost(yProf, yPred, dPred, nProf, cfg);
    end
end
[J_static, dxIdx] = min(J2D, [], 2);
[~, bestIdx] = min(J_static);
gHat = gapQueryGrid(bestIdx);
staticState = struct('gHat', gHat, ...
    'dx0', dxGrid(dxIdx(bestIdx)), 'gapQueryGrid', gapQueryGrid(:), ...
    'dxGrid', dxGrid(:), 'J2D', J2D, 'J_static', J_static(:), ...
    'supportType', classify_gap_support(gHat, templateLib.gapTrain), ...
    'profileX', xProf, 'profileY', yProf, 'profileCount', nProf);
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
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
result = struct('label', label, 'g_used', gFit, 'p', pOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, ...
    'rmse', sqrt(mean(r.^2)), 'coarseCosts', coarseCosts);
end

function result = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfg, jointCfg)
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
resFun = @(theta) residual_calc_joint(theta, t, V, x, templateLib);
thetaOpt = refine_bounded_least_squares_joint(resFun, theta0, lb, ub);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));

r = resFun(thetaOpt);
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
result = struct('label', 'gap-aware-joint-refined', ...
    'g_used', thetaOpt(1), 'dx_used', thetaOpt(2), ...
    'p', thetaOpt(3:end), 'theta', thetaOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'rmse', sqrt(mean(r.^2)));
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

function r = residual_calc_6d_gapaware(p, t, V, x, templateLib, gFit)
u = displacement_2freq(p, t);
VSim = evaluate_template_surface(templateLib, gFit, x - u);
r = V - VSim;
end

function r = residual_calc_joint(theta, t, V, x, templateLib)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = displacement_2freq(p, t);
VSim = evaluate_template_surface(templateLib, gFit, x - dx0 - u);
r = V - VSim;
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + p(4) * sin(2*pi*p(6)*t + p(5));
end

function [freqErr, ampErr] = calc_param_errors(result, fTrue, ATrue)
[fSort, order] = sort(result.f_id);
ASort = result.A_id(order);
[fTrueSort, trueOrder] = sort(fTrue);
ATrueSort = ATrue(trueOrder);
freqErr = mean(abs(fSort - fTrueSort));
ampErr = mean(abs(ASort - ATrueSort));
end

function yPred = evaluate_template_surface(templateLib, g, x)
[curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
yPred = interp1(templateLib.xGrid, curveGrid(:), x, 'pchip', 'extrap');
end

function [xProf, yProf, nProf] = build_raw_profile(x, y, domain, cfg)
xProf = linspace(domain(1), domain(2), cfg.rawProfileGridN)';
yRaw = nan(size(xProf));
nProf = zeros(size(xProf));
edges = [-inf; 0.5 * (xProf(1:end-1) + xProf(2:end)); inf];
bin = discretize(x, edges);
for ib = 1:numel(xProf)
    vals = y(bin == ib);
    nProf(ib) = numel(vals);
    if numel(vals) >= cfg.rawProfileMinCount
        yRaw(ib) = local_percentile(vals, cfg.rawProfilePercentile);
    elseif ~isempty(vals)
        yRaw(ib) = median(vals);
    end
end
valid = isfinite(yRaw);
yProf = interp1(xProf(valid), yRaw(valid), xProf, 'pchip', 'extrap');
yProf = smoothdata(yProf, 'movmean', cfg.rawProfileSmoothWindow);
nProf = max(nProf, 1);
end

function cost = robust_profile_cost(yObs, yPred, dPred, nCount, cfg)
valid = isfinite(yObs) & isfinite(yPred) & isfinite(dPred);
yObs = yObs(valid);
yPred = yPred(valid);
dPred = dPred(valid);
nCount = nCount(valid);
wSlope = abs(dPred);
wSlope = wSlope / max(max(wSlope), eps);
wBase = sqrt(nCount(:));
wBase = wBase / max(max(wBase), eps);
w = wBase .* (cfg.rawWeightFloor + (1 - cfg.rawWeightFloor) * wSlope);
X = [yPred(:), ones(numel(yPred), 1)];
theta = (X' * (w .* X)) \ (X' * (w .* yObs(:)));
for iter = 1:4
    r = yObs(:) - X * theta;
    sigma = 1.4826 * median(abs(r - median(r))) + 1e-9;
    u = abs(r) / (cfg.rawHuberDelta * sigma);
    wHuber = ones(size(u));
    idx = u > 1;
    wHuber(idx) = 1 ./ u(idx);
    wAll = w .* wHuber;
    theta = (X' * (wAll .* X)) \ (X' * (wAll .* yObs(:)));
end
r = yObs(:) - X * theta;
sigma = 1.4826 * median(abs(r - median(r))) + 1e-9;
u = abs(r) / (cfg.rawHuberDelta * sigma);
rho = 0.5 * r.^2;
idx = u > 1;
rho(idx) = cfg.rawHuberDelta * sigma * (abs(r(idx)) - 0.5 * cfg.rawHuberDelta * sigma);
cost = sum(w .* rho) / sum(w);
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

function plot_subset_scatter(x, y, bracketsHoldout)
jitter = 0.08 * (rand(size(x)) - 0.5);
idxB = logical(bracketsHoldout);
plot(x(idxB) + jitter(idxB), y(idxB), 'bo', 'MarkerFaceColor', 'b', ...
    'DisplayName', 'subset brackets holdout');
plot(x(~idxB) + jitter(~idxB), y(~idxB), 'rs', 'MarkerFaceColor', 'r', ...
    'DisplayName', 'subset does not bracket');
legend('Location', 'best', 'Box', 'off');
end

function summary = summarize_group(T, groupName)
summary = table(string(groupName), mean(T.gap_error_mm), mean(T.gapaware_freq_error_Hz), ...
    mean(T.gapaware_amp_error_mm), mean(T.gapaware_rmse), ...
    'VariableNames', {'group', 'mean_gap_error_mm', 'mean_freq_error_Hz', ...
    'mean_amp_error_mm', 'mean_rmse'});
end

function summary = summarize_subset_sizes(T)
subsetSizes = unique(T.subset_size);
rows = repmat(struct('subset_size', 0, 'median_gap_error_mm', 0, ...
    'median_freq_error_Hz', 0, 'median_amp_error_mm', 0, 'median_rmse', 0), numel(subsetSizes), 1);
for i = 1:numel(subsetSizes)
    idx = T.subset_size == subsetSizes(i);
    rows(i).subset_size = subsetSizes(i);
    rows(i).median_gap_error_mm = median(T.gap_error_mm(idx));
    rows(i).median_freq_error_Hz = median(T.freq_error_Hz(idx));
    rows(i).median_amp_error_mm = median(T.amp_error_mm(idx));
    rows(i).median_rmse = median(T.rmse(idx));
end
summary = struct2table(rows);
end

function summary = summarize_bracketing(T)
rows = repmat(struct('brackets_holdout', false, 'mean_gap_error_mm', 0, ...
    'mean_freq_error_Hz', 0, 'mean_amp_error_mm', 0, 'mean_rmse', 0), 2, 1);
for i = 1:2
    flag = i == 1;
    idx = T.brackets_holdout == flag;
    rows(i).brackets_holdout = flag;
    rows(i).mean_gap_error_mm = mean(T.gap_error_mm(idx));
    rows(i).mean_freq_error_Hz = mean(T.freq_error_Hz(idx));
    rows(i).mean_amp_error_mm = mean(T.amp_error_mm(idx));
    rows(i).mean_rmse = mean(T.rmse(idx));
end
summary = struct2table(rows);
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
