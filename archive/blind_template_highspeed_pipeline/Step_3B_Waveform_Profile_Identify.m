%% Step 3B: waveform-level blind profile optimization
% Unknown static clearance effects are absorbed into a shared template T_H(x).
% The dual-frequency vibration is identified by making all dynamically
% corrected high-speed samples collapse onto that shared template.
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'truth', 'gapList', 'xCell', 'yCell');
load(fullfile(outDir, 'stage2B_point_cloud.mat'), 'cloud');

lowTemplate = build_low_template_from_library(gapList, xCell, yCell, cfg.profile.lowTemplateGap);
profileTimer = tic;
profileResult = identify_waveform_profile(cloud, cfg.profile, lowTemplate);
summary = make_profile_summary(profileResult, truth);
profileElapsed = toc(profileTimer);
profileResult.elapsedTimeSeconds = profileElapsed;
summary.profile_elapsed_s = profileElapsed;

save(fullfile(outDir, 'stage3B_waveform_profile_result.mat'), ...
    'profileResult', 'summary', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Blind-Template Step 3B - Waveform Profile Identification', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 24, 12]);
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    plot(truth.templateX, truth.templateY, 'k-', 'LineWidth', 1.1, 'DisplayName', 'truth hidden template');
    plot(profileResult.xGrid, profileResult.T, 'b-', 'LineWidth', 1.2, 'DisplayName', 'blind profile template');
    xlabel('Corrected coordinate \xi (mm)');
    ylabel('Voltage');
    title('Recovered shared template');
    legend('Location', 'best', 'Box', 'off');

    nexttile;
    plot(profileResult.costHist, 'o-', 'LineWidth', 1.2);
    xlabel('Outer iteration');
    ylabel('Profile RMSE');
    title('Profile optimization history');

    nexttile; hold on;
    idxShow = 1:min(2500, numel(cloud.V));
    scatter(cloud.x(idxShow), cloud.V(idxShow), 5, [0.65 0.65 0.65], 'filled', 'DisplayName', 'nominal');
    scatter(profileResult.xi(idxShow), cloud.V(idxShow), 5, profileResult.residual(idxShow), 'filled', 'DisplayName', 'corrected');
    xlabel('Coordinate (mm)');
    ylabel('Voltage');
    title('Point cloud before/after dynamic correction');
    legend('Location', 'best', 'Box', 'off');
    colorbar;

    nexttile;
    bar(categorical({'f1', 'f2'}), [summary.f_est(:), summary.f_true(:)]);
    ylabel('Frequency (Hz)');
    title('Estimated vs truth frequencies');
    legend({'estimated', 'truth'}, 'Location', 'best', 'Box', 'off');

    nexttile;
    bar(categorical({'A1', 'A2'}), [summary.A_est(:), summary.A_true(:)]);
    ylabel('Amplitude (mm)');
    title('Estimated vs truth amplitudes');
    legend({'estimated', 'truth'}, 'Location', 'best', 'Box', 'off');

    nexttile;
    text(0.02, 0.82, sprintf('Mean frequency error: %.4f Hz', summary.mean_freq_error_Hz), 'FontSize', 11);
    text(0.02, 0.62, sprintf('Mean amplitude error: %.4f mm', summary.mean_amp_error_mm), 'FontSize', 11);
    text(0.02, 0.42, sprintf('Profile RMSE: %.4e', summary.profile_rmse), 'FontSize', 11);
    text(0.02, 0.22, sprintf('Samples used: %d', numel(cloud.V)), 'FontSize', 11);
    axis off;
    title('Summary metrics');
end

fprintf('[Step 3B] Waveform-level blind profile summary:\n');
disp(summary);
fprintf('[Step 3B] Profile identification elapsed %.3f s\n', profileElapsed);

function result = identify_waveform_profile(cloud, profileCfg, lowTemplate)
xGrid = cloud.xGrid(:);
T0 = estimate_template_from_cloud(cloud.x, cloud.V, xGrid, profileCfg);
templateAnchor = [];
if get_profile_field(profileCfg, 'useLowTemplateInit', false)
    templateAnchor = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
    [pInitList, initInfo] = initial_candidates_by_low_template(cloud, lowTemplate, xGrid, profileCfg);
    pInitList = expand_low_template_starts(pInitList, profileCfg);
else
    [pInitList, initInfo] = initial_candidates_by_linearized_varpro(cloud, T0, xGrid, profileCfg);
    pInitList = screen_initial_candidates_by_profile(cloud, xGrid, pInitList, profileCfg);
    if profileCfg.useGlobalProfileSearch && exist('particleswarm', 'file') == 2
        [pGlobalList, globalInfo] = global_profile_search(cloud, profileCfg, pInitList);
        pInitList = [pGlobalList; pInitList]; %#ok<AGROW>
        initInfo.globalInfo = globalInfo;
    end
end

candidateList = cell(size(pInitList, 1), 1);
baseScore = inf(size(pInitList, 1), 1);
for ic = 1:size(pInitList, 1)
    candidate = run_profile_outer_loop(cloud, xGrid, T0, pInitList(ic, :), profileCfg, templateAnchor);
    candidate.cvRmse = NaN;
    candidate.selectionScore = candidate.rmse * frequency_penalty_multiplier(candidate.p, profileCfg);
    candidateList{ic} = candidate;
    baseScore(ic) = candidate.selectionScore;
end
[~, candidateOrder] = sort(baseScore, 'ascend');
nCv = min(get_profile_field(profileCfg, 'cvCandidateCount', numel(candidateOrder)), numel(candidateOrder));
for ii = 1:nCv
    ic = candidateOrder(ii);
    candidate = candidateList{ic};
    candidate.cvRmse = calc_template_cv_rmse(cloud, templateAnchor, xGrid, candidate.p, profileCfg);
    candidate.selectionScore = blend_selection_rmse(candidate.rmse, candidate.cvRmse, profileCfg) * ...
        frequency_penalty_multiplier(candidate.p, profileCfg);
    candidateList{ic} = candidate;
end
selectionScore = cellfun(@(c) c.selectionScore, candidateList);
[~, bestIdx] = min(selectionScore);

result = candidateList{bestIdx};
result.T0 = T0(:);
result.initInfo = initInfo;
result.pInitList = pInitList;
end

function lowTemplate = build_low_template_from_library(gapList, xCell, yCell, lowGap)
idxLow = find(abs(gapList(:) - lowGap) < 1e-12, 1);
if isempty(idxLow)
    [~, idxLow] = min(abs(gapList(:) - lowGap));
end
lowTemplate = struct();
lowTemplate.g = gapList(idxLow);
lowTemplate.x = xCell{idxLow}(:);
lowTemplate.y = yCell{idxLow}(:);
end

function [pList, initInfo] = initial_candidates_by_low_template(cloud, lowTemplate, xGrid, profileCfg)
Tlow = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
TAtX = interp1(xGrid, Tlow, cloud.x, 'pchip', 'extrap');
dT = gradient(Tlow(:), mean(diff(xGrid)));
dTAtX = interp1(xGrid, dT, cloud.x, 'pchip', 'extrap');
valid = isfinite(TAtX) & isfinite(dTAtX) & isfinite(cloud.V);
derivScale = max(abs(dTAtX(valid)));
valid = valid & abs(dTAtX) > profileCfg.derivativeFloorFrac * max(derivScale, eps);

if get_profile_field(profileCfg, 'lowInitUseGainBias', true)
    theta = [TAtX(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
    gain0 = theta(1);
    bias0 = theta(2);
else
    gain0 = 1;
    bias0 = 0;
end

f1Grid = profileCfg.f1Grid(:)';
f2Grid = profileCfg.f2Grid(:)';
keepCount = get_profile_field(profileCfg, 'lowInitKeep', 8);
maxKeep = max(keepCount * 6, keepCount);
candCost = inf(maxKeep, 1);
candP = zeros(maxKeep, 6);
lambda = 1e-9 * max(nnz(valid), 1);

t = cloud.t(valid);
y = cloud.V(valid) - (gain0 * TAtX(valid) + bias0);
kx = gain0 * dTAtX(valid);
for f1 = f1Grid
    s1 = sin(2*pi*f1*t);
    c1 = cos(2*pi*f1*t);
    for f2 = f2Grid
        if abs(f1 - f2) < profileCfg.minFreqSeparation
            continue;
        end
        s2 = sin(2*pi*f2*t);
        c2 = cos(2*pi*f2*t);
        H = [-kx .* s1, -kx .* c1, -kx .* s2, -kx .* c2];
        beta = (H' * H + lambda * eye(4)) \ (H' * y);
        res = y - H * beta;
        cost = sum(res.^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            candCost(worstIdx) = cost;
            candP(worstIdx, :) = [hypot(beta(1), beta(2)), atan2(beta(2), beta(1)), f1, ...
                hypot(beta(3), beta(4)), atan2(beta(4), beta(3)), f2];
        end
    end
end

[candCost, order] = sort(candCost, 'ascend');
candP = candP(order, :);
keep = isfinite(candCost);
candCost = candCost(keep);
candP = candP(keep, :);
if isempty(candP)
    error('Low-template initialization did not produce any valid frequency candidate.');
end

lb = [0, -pi, min(f1Grid), 0, -pi, min(f2Grid)];
ub = [profileCfg.ampUpperBound, pi, max(f1Grid), profileCfg.ampUpperBound, pi, max(f2Grid)];
refineCount = get_profile_field(profileCfg, 'lowInitRefineCount', maxKeep);
numRaw = min(refineCount, size(candP, 1));
refinedP = zeros(numRaw, 6);
refinedCost = inf(numRaw, 1);
for i = 1:numRaw
    p0 = min(max(candP(i, :), lb), ub);
    pOpt = refine_p_given_low_template(cloud, Tlow, xGrid, p0, lb, ub, profileCfg);
    refinedP(i, :) = pOpt;
    refinedCost(i) = calc_low_template_rmse(cloud, Tlow, xGrid, pOpt, profileCfg);
end

[refinedCost, order] = sort(refinedCost, 'ascend');
refinedP = refinedP(order, :);
refinedScore = refinedCost .* frequency_penalty_multiplier(refinedP, profileCfg);
[~, scoreOrder] = sort(refinedScore, 'ascend');
refinedP = refinedP(scoreOrder, :);
refinedCost = refinedCost(scoreOrder);
key = round([refinedP(:, 1), refinedP(:, 3), refinedP(:, 4), refinedP(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
refinedP = refinedP(ia, :);
refinedCost = refinedCost(ia);
nOut = min(keepCount, size(refinedP, 1));
pList = refinedP(1:nOut, :);
initInfo = struct('method', 'low_template_guided', ...
    'lowTemplateGap', lowTemplate.g, ...
    'lowTemplateGain0', gain0, ...
    'lowTemplateBias0', bias0, ...
    'lowTemplateRmse', refinedCost(1:nOut), ...
    'linearizedCost', candCost(1:min(numRaw, numel(candCost))), ...
    'numValidLinearizedSamples', nnz(valid));
end

function pStart = expand_low_template_starts(pInitList, profileCfg)
ampScales = get_profile_field(profileCfg, 'lowInitAmpScales', 1);
freqJitters = get_profile_field(profileCfg, 'lowInitFreqJitters', 0);
numStarts = get_profile_field(profileCfg, 'numProfileStarts', size(pInitList, 1));
starts = [];
baseCount = min(2, size(pInitList, 1));
for ib = 1:baseCount
    p = pInitList(ib, :);
    for ia = 1:numel(ampScales)
        for jf1 = 1:numel(freqJitters)
            for jf2 = 1:numel(freqJitters)
                ps = p;
                ps(1) = min(profileCfg.ampUpperBound, max(0, ampScales(ia) * p(1)));
                ps(4) = min(profileCfg.ampUpperBound, max(0, ampScales(ia) * p(4)));
                ps(3) = min(max(p(3) + freqJitters(jf1), min(profileCfg.f1Grid)), max(profileCfg.f1Grid));
                ps(6) = min(max(p(6) + freqJitters(jf2), min(profileCfg.f2Grid)), max(profileCfg.f2Grid));
                starts = [starts; ps]; %#ok<AGROW>
            end
        end
    end
end
starts = [pInitList; starts]; %#ok<AGROW>
key = round([starts(:, 1), starts(:, 2), starts(:, 3), starts(:, 4), starts(:, 5), starts(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
starts = starts(ia, :);
pStart = starts(1:min(numStarts, size(starts, 1)), :);
end

function pOpt = refine_p_given_low_template(cloud, Tlow, xGrid, p0, lb, ub, profileCfg)
p0 = min(max(p0, lb), ub);
resFun = @(p) low_template_residual(p, cloud, Tlow, xGrid, profileCfg);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'MaxIterations', profileCfg.maxLsqIter);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 5 * profileCfg.maxLsqIter, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function rmse = calc_low_template_rmse(cloud, Tlow, xGrid, p, profileCfg)
r = low_template_residual(p, cloud, Tlow, xGrid, profileCfg);
rmse = sqrt(mean(r.^2));
end

function r = low_template_residual(p, cloud, Tlow, xGrid, profileCfg)
xi = cloud.x - displacement_2freq(p, cloud.t);
Tfit = interp1(xGrid, Tlow, xi, 'pchip', NaN);
valid = isfinite(Tfit) & isfinite(cloud.V);
if nnz(valid) < 0.75 * numel(cloud.V)
    templateRange = max(Tlow) - min(Tlow);
    r = 5 * max(templateRange, eps) * ones(size(cloud.V));
    return;
end
if get_profile_field(profileCfg, 'lowInitUseGainBias', true)
    theta = [Tfit(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
    Vpred = theta(1) * Tfit + theta(2);
else
    Vpred = Tfit;
end
r = cloud.V - Vpred;
invalid = ~isfinite(r);
if any(invalid)
    templateRange = max(Tlow) - min(Tlow);
    r(invalid) = 5 * max(templateRange, eps);
end
end

function [pGlobalList, globalInfo] = global_profile_search(cloud, profileCfg, pSeedList)
[lb, ub] = make_profile_refine_bounds(p0, profileCfg);
xGridGlobal = linspace(cloud.domain(1), cloud.domain(2), profileCfg.globalTemplateGridN)';
obj = @(p) profile_objective_for_global(p, cloud, xGridGlobal, profileCfg);
opts = optimoptions('particleswarm', 'Display', 'none', ...
    'SwarmSize', profileCfg.globalSwarmSize, ...
    'MaxIterations', profileCfg.globalMaxIter, ...
    'FunctionTolerance', 1e-9, ...
    'UseParallel', false);

numRuns = get_profile_field(profileCfg, 'globalNumRuns', 1);
keepCount = get_profile_field(profileCfg, 'globalKeepCount', numRuns);
usePattern = get_profile_field(profileCfg, 'globalUsePatternSearch', false) && exist('patternsearch', 'file') == 2;
patMaxIter = get_profile_field(profileCfg, 'globalPatternMaxIter', 60);

allP = zeros(0, 6);
allF = zeros(0, 1);
exitflag = zeros(numRuns, 1);
iter = zeros(numRuns, 1);
funccount = zeros(numRuns, 1);
for ir = 1:numRuns
    rng(7300 + ir, 'twister');
    [pBest, fBest, exitflag(ir), output] = particleswarm(obj, 6, lb, ub, opts);
    iter(ir) = output.iterations;
    funccount(ir) = output.funccount;
    if usePattern
        patOpts = optimoptions('patternsearch', 'Display', 'none', ...
            'MaxIterations', patMaxIter, ...
            'FunctionTolerance', 1e-10, ...
            'StepTolerance', 1e-8, ...
            'UseCompletePoll', true, ...
            'UseCompleteSearch', false);
        [pPat, fPat] = patternsearch(obj, pBest, [], [], [], [], lb, ub, [], patOpts);
        if fPat < fBest
            pBest = pPat;
            fBest = fPat;
        end
    end
    pBest(2) = wrap_to_pi_local(pBest(2));
    pBest(5) = wrap_to_pi_local(pBest(5));
    allP(end+1, :) = pBest(:).'; %#ok<AGROW>
    allF(end+1, 1) = fBest; %#ok<AGROW>
end

useSurrogate = get_profile_field(profileCfg, 'globalUseSurrogateSearch', false) && exist('surrogateopt', 'file') == 2;
surrogateInfo = [];
if useSurrogate
    pSeed = sanitize_seed_points(pSeedList, lb, ub);
    if ~isempty(allP)
        pSeed = [allP; pSeed]; %#ok<AGROW>
    end
    maxEvals = get_profile_field(profileCfg, 'globalSurrogateMaxEvals', 1000);
    surOpts = optimoptions('surrogateopt', 'Display', 'none', ...
        'MaxFunctionEvaluations', maxEvals, ...
        'UseParallel', false);
    if ~isempty(pSeed)
        surOpts.InitialPoints = pSeed;
    end
    rng(9100, 'twister');
    [pSur, fSur, surExit, surOutput] = surrogateopt(obj, lb, ub, surOpts);
    pSur(2) = wrap_to_pi_local(pSur(2));
    pSur(5) = wrap_to_pi_local(pSur(5));
    allP(end+1, :) = pSur(:).'; %#ok<AGROW>
    allF(end+1, 1) = fSur; %#ok<AGROW>
    surrogateInfo = struct('objective', fSur, 'exitflag', surExit, ...
        'funccount', surOutput.funccount);
end

[allF, order] = sort(allF, 'ascend');
allP = allP(order, :);
key = round([allP(:, 1), allP(:, 3), allP(:, 4), allP(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
allP = allP(ia, :);
allF = allF(ia);
nKeep = min(keepCount, size(allP, 1));
pGlobalList = allP(1:nKeep, :);
globalInfo = struct('objective', allF(1:nKeep), 'allObjective', allF, ...
    'exitflag', exitflag, 'iterations', iter, 'funccount', funccount, ...
    'usedPatternSearch', usePattern, 'usedSurrogateSearch', useSurrogate, ...
    'surrogateInfo', surrogateInfo);
end

function pSeed = sanitize_seed_points(pSeedList, lb, ub)
if isempty(pSeedList)
    pSeed = zeros(0, numel(lb));
    return;
end
pSeed = pSeedList(:, 1:numel(lb));
pSeed(:, 2) = arrayfun(@wrap_to_pi_local, pSeed(:, 2));
pSeed(:, 5) = arrayfun(@wrap_to_pi_local, pSeed(:, 5));
pSeed = min(max(pSeed, lb), ub);
key = round([pSeed(:, 1), pSeed(:, 3), pSeed(:, 4), pSeed(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
pSeed = pSeed(ia, :);
end

function val = profile_objective_for_global(p, cloud, xGrid, profileCfg)
if abs(p(3) - p(6)) < profileCfg.minFreqSeparation
    val = 1e6;
    return;
end
xi = cloud.x - displacement_2freq(p, cloud.t);
[rmse, validFrac] = smooth_profile_rmse(xi, cloud.V, xGrid, profileCfg);
if ~isfinite(rmse) || validFrac < 0.75
    val = 1e3 + (1 - validFrac);
    return;
end
val = rmse;
end

function [rmse, validFrac] = smooth_profile_rmse(xi, V, xGrid, profileCfg)
xi = xi(:);
V = V(:);
xGrid = xGrid(:);
[S, valid] = linear_interp_matrix(xi, xGrid);
valid = valid & isfinite(V);
validFrac = nnz(valid) / max(numel(V), 1);
if nnz(valid) < 0.75 * numel(V)
    rmse = inf;
    return;
end

S = S(valid, :);
v = V(valid);
nGrid = numel(xGrid);
lambda = get_profile_field(profileCfg, 'globalSmoothPenalty', 0);
ridge = 1e-10 * speye(nGrid);
if nGrid > 2 && lambda > 0
    e = ones(nGrid, 1);
    D2 = spdiags([e, -2*e, e], 0:2, nGrid-2, nGrid);
    lhs = S' * S + lambda * (D2' * D2) + ridge;
else
    lhs = S' * S + ridge;
end
T = lhs \ (S' * v);
r = v - S * T;
rmse = sqrt(mean(r.^2));
end

function [S, valid] = linear_interp_matrix(xi, xGrid)
n = numel(xi);
nGrid = numel(xGrid);
dx = mean(diff(xGrid));
q = (xi - xGrid(1)) / dx + 1;
j = floor(q);
a = q - j;
valid = isfinite(q) & j >= 1 & j < nGrid;
row = find(valid);
jj = j(valid);
aa = a(valid);
S = sparse([row; row], [jj; jj + 1], [1 - aa; aa], n, nGrid);
end

function pStart = screen_initial_candidates_by_profile(cloud, xGrid, pInitList, profileCfg)
aug = [];
sourceIdx = [];
scaleVal = [];
for i = 1:size(pInitList, 1)
    p = pInitList(i, :);
    phaseFlips = [0, 0; pi, 0; 0, pi; pi, pi];
    for scale = profileCfg.initAmpScales
        for iph = 1:size(phaseFlips, 1)
            ps = p;
            ps(1) = min(profileCfg.ampUpperBound, max(0, scale * p(1)));
            ps(4) = min(profileCfg.ampUpperBound, max(0, scale * p(4)));
            ps(2) = wrap_to_pi_local(ps(2) + phaseFlips(iph, 1));
            ps(5) = wrap_to_pi_local(ps(5) + phaseFlips(iph, 2));
            aug = [aug; ps]; %#ok<AGROW>
            sourceIdx = [sourceIdx; i]; %#ok<AGROW>
            scaleVal = [scaleVal; scale]; %#ok<AGROW>
        end
    end
end
if isempty(aug)
    pStart = pInitList;
    return;
end
key = round([aug(:, 1), aug(:, 2), aug(:, 3), aug(:, 4), aug(:, 5), aug(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
aug = aug(ia, :);
sourceIdx = sourceIdx(ia);
scaleVal = scaleVal(ia);

score = inf(size(aug, 1), 1);
for i = 1:size(aug, 1)
    xi = cloud.x - displacement_2freq(aug(i, :), cloud.t);
    T = estimate_template_from_cloud(xi, cloud.V, xGrid, profileCfg);
    VFit = interp1(xGrid, T, xi, 'pchip', NaN);
    r = cloud.V - VFit;
    score(i) = sqrt(mean(r(isfinite(r)).^2));
end
[~, order] = sort(score, 'ascend');
profileKeep = order(1:min(profileCfg.numProfileStarts, numel(order)));

fullScaleCount = get_profile_field(profileCfg, 'numLinearFullScaleStarts', 0);
fullKeep = find(sourceIdx <= min(fullScaleCount, size(pInitList, 1)));

robustCount = get_profile_field(profileCfg, 'numLinearRobustStarts', 0);
robustScales = get_profile_field(profileCfg, 'linearRobustAmpScales', []);
if isempty(robustScales)
    robustKeep = [];
else
    isRobustScale = false(size(scaleVal));
    for is = 1:numel(robustScales)
        isRobustScale = isRobustScale | abs(scaleVal - robustScales(is)) < 1e-12;
    end
    robustKeep = find(sourceIdx <= min(robustCount, size(pInitList, 1)) & isRobustScale);
end

diverseCount = get_profile_field(profileCfg, 'numLinearDiverseStarts', 0);
diverseKeep = [];
for src = 1:min(diverseCount, size(pInitList, 1))
    idx = find(sourceIdx == src);
    if isempty(idx)
        continue;
    end
    [~, bestLocal] = min(score(idx));
    diverseKeep = [diverseKeep; idx(bestLocal)]; %#ok<AGROW>
end

keep = unique([fullKeep(:); robustKeep(:); diverseKeep(:); profileKeep(:)], 'stable');
if numel(keep) > profileCfg.numProfileStarts
    keep = keep(1:profileCfg.numProfileStarts);
end
pStart = aug(keep, :);
end

function candidate = run_profile_outer_loop(cloud, xGrid, T0, p0, profileCfg, templateAnchor)
p = p0(:).';
xi0 = cloud.x - displacement_2freq(p, cloud.t);
if get_profile_field(profileCfg, 'useLowToHighTemplateTransfer', false) && ~isempty(templateAnchor)
    T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, profileCfg);
else
    T = estimate_template_from_cloud(xi0, cloud.V, xGrid, profileCfg);
end
if any(~isfinite(T))
    T = T0(:);
end
costHist = zeros(profileCfg.maxOuterIter, 1);
paramHist = zeros(profileCfg.maxOuterIter, numel(p));
templateChangeHist = zeros(profileCfg.maxOuterIter, 1);
bestLoop = struct('rmse', inf, 'p', p, 'T', T, 'iter', 0);

lb = [0, -pi, min(profileCfg.f1Grid), 0, -pi, min(profileCfg.f2Grid)];
ub = [profileCfg.ampUpperBound, pi, max(profileCfg.f1Grid), ...
    profileCfg.ampUpperBound, pi, max(profileCfg.f2Grid)];

for iter = 1:profileCfg.maxOuterIter
    pPrev = p;
    Tprev = T;
    p = refine_p_given_template(cloud, T, xGrid, p, lb, ub, profileCfg);
    xi = cloud.x - displacement_2freq(p, cloud.t);
    if get_profile_field(profileCfg, 'useLowToHighTemplateTransfer', false) && ~isempty(templateAnchor)
        T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, profileCfg);
    else
        T = estimate_template_from_cloud(xi, cloud.V, xGrid, profileCfg);
    end
    [rmse, ~, ~] = calc_profile_rmse(cloud, T, xGrid, p);
    costHist(iter) = rmse;
    paramHist(iter, :) = p;
    templateChangeHist(iter) = norm(T - Tprev) / max(norm(Tprev), eps);
    if rmse < bestLoop.rmse
        bestLoop.rmse = rmse;
        bestLoop.p = p;
        bestLoop.T = T;
        bestLoop.iter = iter;
    end
    relParam = norm(p - pPrev) / max(norm(pPrev), eps);
    if iter > 1 && relParam < profileCfg.tolParam && templateChangeHist(iter) < profileCfg.tolTemplate
        costHist = costHist(1:iter);
        paramHist = paramHist(1:iter, :);
        templateChangeHist = templateChangeHist(1:iter);
        break;
    end
end

p = bestLoop.p;
T = bestLoop.T;
if p(3) > p(6)
    p = [p(4), p(5), p(6), p(1), p(2), p(3)];
end
[rmse, residual, VFit, xi] = calc_profile_rmse(cloud, T, xGrid, p);
candidate = struct();
candidate.p = p;
candidate.A = [p(1), p(4)];
candidate.phi = [p(2), p(5)];
candidate.f = [p(3), p(6)];
candidate.T = T(:);
candidate.xGrid = xGrid(:);
candidate.xi = xi(:);
candidate.VFit = VFit(:);
candidate.residual = residual(:);
candidate.rmse = rmse;
candidate.costHist = costHist(:);
candidate.paramHist = paramHist;
candidate.templateChangeHist = templateChangeHist(:);
candidate.bestIter = bestLoop.iter;
end

function T = transfer_template_from_anchor(cloud, Tanchor, xGrid, p, profileCfg)
xi = cloud.x - displacement_2freq(p, cloud.t);
TanchorAtXi = interp1(xGrid, Tanchor, xi, 'pchip', NaN);
valid = isfinite(TanchorAtXi) & isfinite(cloud.V);
if nnz(valid) < 0.75 * numel(cloud.V)
    T = estimate_template_from_cloud(xi, cloud.V, xGrid, profileCfg);
    return;
end

if get_profile_field(profileCfg, 'transferUseGainBias', true)
    theta = [TanchorAtXi(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
    gain = theta(1);
    bias = theta(2);
else
    gain = 1;
    bias = 0;
end

residual = cloud.V - (gain * TanchorAtXi + bias);
residualGrid = estimate_template_from_cloud(xi(valid), residual(valid), xGrid, profileCfg);
smoothWin = get_profile_field(profileCfg, 'transferResidualSmoothWindow', profileCfg.smoothWindow);
if smoothWin > 2
    residualGrid = smoothdata(residualGrid, 'sgolay', smoothWin);
end
clipLimit = get_profile_field(profileCfg, 'transferResidualClipFrac', inf) * max(range(Tanchor), eps);
if isfinite(clipLimit)
    residualGrid = min(max(residualGrid, -clipLimit), clipLimit);
end
T = gain * Tanchor(:) + bias + residualGrid(:);
T = fill_template_gaps(T, xGrid, isfinite(T));
T = smoothdata(T, 'sgolay', profileCfg.smoothWindow);
end

function [pList, initInfo] = initial_candidates_by_linearized_varpro(cloud, T, xGrid, profileCfg)
TAtX = interp1(xGrid, T, cloud.x, 'pchip', 'extrap');
dT = gradient(T(:), mean(diff(xGrid)));
dTAtX = interp1(xGrid, dT, cloud.x, 'pchip', 'extrap');
DV = cloud.V - TAtX;
valid = isfinite(DV) & isfinite(dTAtX);
derivScale = max(abs(dTAtX(valid)));
valid = valid & abs(dTAtX) > profileCfg.derivativeFloorFrac * max(derivScale, eps);

f1Grid = profileCfg.f1Grid(:)';
f2Grid = profileCfg.f2Grid(:)';
maxKeep = max(profileCfg.numInitCandidates * 4, profileCfg.numInitCandidates);
candCost = inf(maxKeep, 1);
candP = zeros(maxKeep, 6);
lambda = 1e-9 * max(nnz(valid), 1);

t = cloud.t(valid);
y = DV(valid);
kx = dTAtX(valid);
for f1 = f1Grid
    s1 = sin(2*pi*f1*t);
    c1 = cos(2*pi*f1*t);
    for f2 = f2Grid
        if abs(f1 - f2) < profileCfg.minFreqSeparation
            continue;
        end
        s2 = sin(2*pi*f2*t);
        c2 = cos(2*pi*f2*t);
        H = [-kx .* s1, -kx .* c1, -kx .* s2, -kx .* c2];
        beta = (H' * H + lambda * eye(4)) \ (H' * y);
        res = y - H * beta;
        cost = sum(res.^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            candCost(worstIdx) = cost;
            candP(worstIdx, :) = [hypot(beta(1), beta(2)), atan2(beta(2), beta(1)), f1, ...
                hypot(beta(3), beta(4)), atan2(beta(4), beta(3)), f2];
        end
    end
end

[candCost, order] = sort(candCost, 'ascend');
candP = candP(order, :);
keep = isfinite(candCost);
candCost = candCost(keep);
candP = candP(keep, :);
nOut = min(profileCfg.numInitCandidates, size(candP, 1));
pList = candP(1:nOut, :);
if isempty(pList)
    error('Could not generate waveform-profile initial candidates.');
end
initInfo = struct('cost', candCost(1:nOut), 'numValidLinearizedSamples', nnz(valid));
end

function pOpt = refine_p_given_template(cloud, T, xGrid, p0, lb, ub, profileCfg)
p0 = min(max(p0, lb), ub);
templateRange = max(T) - min(T);
templateEval = make_profile_template_eval(T, xGrid, profileCfg);
resFun = @(p) fixed_template_residual(p, cloud, templateRange, templateEval);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'MaxIterations', profileCfg.maxLsqIter);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 5 * profileCfg.maxLsqIter, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function [lb, ub] = make_profile_refine_bounds(p0, profileCfg)
lbWide = [0, -pi, min(profileCfg.f1Grid), 0, -pi, min(profileCfg.f2Grid)];
ubWide = [profileCfg.ampUpperBound, pi, max(profileCfg.f1Grid), ...
    profileCfg.ampUpperBound, pi, max(profileCfg.f2Grid)];
if ~get_profile_field(profileCfg, 'useLocalRefineBounds', false)
    lb = lbWide;
    ub = ubWide;
    return;
end
freqHalfWidth = get_profile_field(profileCfg, 'localFreqHalfWidth', 20);
ampScaleLo = get_profile_field(profileCfg, 'localAmpScaleLo', 0.5);
ampScaleHi = get_profile_field(profileCfg, 'localAmpScaleHi', 1.5);
ampFloor = 0.005;
lb = lbWide;
ub = ubWide;
lb(1) = max(lbWide(1), min(p0(1) * ampScaleLo, max(p0(1) - ampFloor, 0)));
ub(1) = min(ubWide(1), max(p0(1) * ampScaleHi, p0(1) + ampFloor));
lb(4) = max(lbWide(4), min(p0(4) * ampScaleLo, max(p0(4) - ampFloor, 0)));
ub(4) = min(ubWide(4), max(p0(4) * ampScaleHi, p0(4) + ampFloor));
lb(3) = max(lbWide(3), p0(3) - freqHalfWidth);
ub(3) = min(ubWide(3), p0(3) + freqHalfWidth);
lb(6) = max(lbWide(6), p0(6) - freqHalfWidth);
ub(6) = min(ubWide(6), p0(6) + freqHalfWidth);
end

function T = estimate_template_from_cloud(xi, V, xGrid, profileCfg)
xi = xi(:);
V = V(:);
xGrid = xGrid(:);
edges = [-inf; 0.5 * (xGrid(1:end-1) + xGrid(2:end)); inf];
bin = discretize(xi, edges);
validBin = isfinite(V) & isfinite(bin);
T = accumarray(bin(validBin), V(validBin), [numel(xGrid), 1], @mean, NaN);
valid = isfinite(T);
if nnz(valid) < 5
    error('Too few samples to estimate the shared template.');
end
T = fill_template_gaps(T, xGrid, valid);
T = smoothdata(T, 'sgolay', profileCfg.smoothWindow);
end

function [rmse, residual, VFit, xi] = calc_profile_rmse(cloud, T, xGrid, p)
xi = cloud.x - displacement_2freq(p, cloud.t);
VFit = interp1(xGrid, T, xi, 'pchip', NaN);
residual = cloud.V - VFit;
valid = isfinite(residual);
rmse = sqrt(mean(residual(valid).^2));
end

function cvRmse = calc_template_cv_rmse(cloud, templateAnchor, xGrid, p, profileCfg)
if isempty(templateAnchor)
    cvRmse = NaN;
    return;
end
rev = cloud.rev(:);
splitA = mod(rev, 2) == 0;
splitB = ~splitA;
if nnz(splitA) < 0.25 * numel(rev) || nnz(splitB) < 0.25 * numel(rev)
    cvRmse = NaN;
    return;
end
cloudA = subset_cloud(cloud, splitA);
cloudB = subset_cloud(cloud, splitB);
TA = transfer_template_from_anchor(cloudA, templateAnchor, xGrid, p, profileCfg);
TB = transfer_template_from_anchor(cloudB, templateAnchor, xGrid, p, profileCfg);
rmseAB = calc_profile_rmse(cloudB, TA, xGrid, p);
rmseBA = calc_profile_rmse(cloudA, TB, xGrid, p);
cvRmse = mean([rmseAB, rmseBA], 'omitnan');
end

function cloudOut = subset_cloud(cloud, idx)
cloudOut = cloud;
cloudOut.t = cloud.t(idx);
cloudOut.x = cloud.x(idx);
cloudOut.V = cloud.V(idx);
cloudOut.rev = cloud.rev(idx);
cloudOut.sensor = cloud.sensor(idx);
end

function scoreRmse = blend_selection_rmse(trainRmse, cvRmse, profileCfg)
if ~isfinite(cvRmse)
    scoreRmse = trainRmse;
    return;
end
w = get_profile_field(profileCfg, 'selectionCvWeight', 0);
scoreRmse = (1 - w) * trainRmse + w * cvRmse;
end

function r = fixed_template_residual(p, cloud, templateRange, templateEval)
xi = cloud.x - displacement_2freq(p, cloud.t);
VFit = templateEval.F(xi);
w = informative_profile_weights(xi, templateEval);
r = sqrt(w) .* (cloud.V - VFit);
invalid = ~isfinite(r);
if any(invalid)
    r(invalid) = 5 * max(templateRange, eps);
end
end

function templateEval = make_profile_template_eval(T, xGrid, profileCfg)
xGrid = xGrid(:);
T = T(:);
dT = gradient(T, mean(diff(xGrid)));
templateEval = struct();
templateEval.F = griddedInterpolant(xGrid, T, 'pchip', 'none');
templateEval.dF = griddedInterpolant(xGrid, dT, 'pchip', 'nearest');
templateEval.weightFloor = get_profile_field(profileCfg, 'fitWeightFloor', 0);
templateEval.weightPower = get_profile_field(profileCfg, 'fitWeightPower', 1);
end

function w = informative_profile_weights(xi, templateEval)
dTAtXi = templateEval.dF(xi);
slopeScale = max(abs(dTAtXi));
if slopeScale <= eps
    w = ones(size(xi));
    return;
end
wShape = (abs(dTAtXi) ./ slopeScale) .^ templateEval.weightPower;
w = templateEval.weightFloor + (1 - templateEval.weightFloor) * wShape;
w(~isfinite(w)) = templateEval.weightFloor;
end

function mult = frequency_penalty_multiplier(p, profileCfg)
if isempty(p)
    mult = [];
    return;
end
p = reshape(p, [], 6);
f1 = p(:, 3);
f2 = p(:, 6);
edgeGuardHz = get_profile_field(profileCfg, 'edgeGuardHz', 0);
edgePenaltyWeight = get_profile_field(profileCfg, 'edgePenaltyWeight', 0);
hit1 = f1 <= min(profileCfg.f1Grid) + edgeGuardHz | ...
    f1 >= max(profileCfg.f1Grid) - edgeGuardHz;
hit2 = f2 <= min(profileCfg.f2Grid) + edgeGuardHz | ...
    f2 >= max(profileCfg.f2Grid) - edgeGuardHz;
mult = 1 + edgePenaltyWeight * (double(hit1) + double(hit2));
end

function T = fill_template_gaps(T, xGrid, valid)
idx = (1:numel(T))';
valid = valid(:) & isfinite(T(:));
if nnz(valid) < 2
    T(~valid) = 0;
    return;
end
first = find(valid, 1, 'first');
last = find(valid, 1, 'last');
left = idx < first;
right = idx > last;
midMissing = ~valid & idx >= first & idx <= last;
T(left) = T(first);
T(right) = T(last);
T(midMissing) = interp1(xGrid(valid), T(valid), xGrid(midMissing), 'pchip');
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end

function value = get_profile_field(profileCfg, name, defaultValue)
if isfield(profileCfg, name) && ~isempty(profileCfg.(name))
    value = profileCfg.(name);
else
    value = defaultValue;
end
end

function summary = make_profile_summary(profileResult, truth)
[fTrue, trueOrder] = sort(truth.f(:).');
ATrue = truth.A(trueOrder);
[fEst, estOrder] = sort(profileResult.f(:).');
AEst = profileResult.A(estOrder);
summary = struct();
summary.f_true = fTrue;
summary.f_est = fEst;
summary.A_true = ATrue;
summary.A_est = AEst;
summary.mean_freq_error_Hz = mean(abs(fEst - fTrue));
summary.mean_amp_error_mm = mean(abs(AEst - ATrue));
summary.profile_rmse = profileResult.rmse;
summary.p = profileResult.p;
end
