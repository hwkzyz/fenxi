function R = run_vp_main(highMap, templateLib, cfg, staticState)
%run_vp_main  Joint variable-projection estimate over gap and frequencies.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);

gapGrid = linspace(max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin), ...
    max(templateLib.gapTrain) + cfg.rawGapSearchMargin, cfg.gapQueryN);
if isfield(staticState, 'gHat') && ~isempty(staticState.gHat)
    gCenter = staticState.gHat;
else
    gCenter = staticState.g0;
    if isempty(gCenter)
        gCenter = mean(templateLib.gapTrain);
    end
end
if ~isempty(gCenter) && isfinite(gCenter)
    gapGrid = linspace(max(0.05, gCenter - 0.06), gCenter + 0.06, 61);
    gapGrid = gapGrid(gapGrid >= min(templateLib.gapTrain) - cfg.rawGapSearchMargin & ...
        gapGrid <= max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
end
precomp = build_vp_precompute(t, V, x, templateLib, gapGrid, cfg);

bestFit = [];
bestG = gapGrid(1);
bestSSE = inf;
coarseCount = min(21, numel(gapGrid));
coarseIdx = round(linspace(1, numel(gapGrid), coarseCount));
coarseIdx = unique(coarseIdx(:));
coarseRows = zeros(numel(coarseIdx), 3);

for ic = 1:numel(coarseIdx)
    ig = coarseIdx(ic);
    gTry = gapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    coarseRows(ic, :) = [gTry, fitTry.J1, fitTry.rhoJ];
end

[~, coarseOrder] = sort(coarseRows(:, 2), 'ascend');
bestCoarseJ = coarseRows(coarseOrder(1), 2);
keepCoarse = min(5, numel(coarseOrder));
relTol = 0.02;
extraKeep = coarseOrder(coarseRows(coarseOrder, 2) <= bestCoarseJ * (1 + relTol));
keepList = unique([coarseOrder(1:keepCoarse); extraKeep(:)], 'stable');
keepList = keepList(1:min(numel(keepList), 7));
fineMask = false(numel(gapGrid), 1);
for ic = 1:numel(keepList)
    ig0 = coarseIdx(keepList(ic));
    lo = max(1, ig0 - 5);
    hi = min(numel(gapGrid), ig0 + 5);
    fineMask(lo:hi) = true;
end
fineIdx = find(fineMask);

scoreRows = zeros(numel(fineIdx), 3);
fitBank = cell(numel(fineIdx), 1);
for jf = 1:numel(fineIdx)
    ig = fineIdx(jf);
    gTry = gapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    scoreRows(jf, :) = [gTry, fitTry.J1, fitTry.rhoJ];
    fitBank{jf} = fitTry;
end

[~, order] = sort(scoreRows(:, 2), 'ascend');
numRefine = min(5, numel(order));
for ir = 1:numRefine
    jf = order(ir);
    gTry = scoreRows(jf, 1);
    fitTry = fitBank{jf};
    fitTry = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    resTry = V - eval_gap_template(templateLib, gTry, x - fitTry.delta_sample);
    sseTry = dot(resTry, resTry);
    if sseTry < bestSSE
        bestSSE = sseTry;
        bestFit = fitTry;
        bestG = gTry;
    end
end

if isempty(bestFit)
    jf = order(1);
    gTry = scoreRows(jf, 1);
    fitTry = fitBank{jf};
    fitTry = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    bestFit = fitTry;
    bestG = gTry;
end

VFit = eval_gap_template(templateLib, bestG, x - bestFit.delta_sample);
R.method = "vp_main";
R.g_used = bestG;
R.dx_used = bestFit.dx;
R.fit = bestFit;
R.f_id = bestFit.f;
R.A_id = bestFit.A;
R.phi_id = bestFit.phi;
R.VFit = VFit;
R.rmse = sqrt(mean((V - VFit).^2));
R.eta_g = eta_gap(templateLib, bestG);
R.deltaJ = bestFit.deltaJ;
R.rhoJ = bestFit.rhoJ;
R.J = scoreRows(:, 2);
R.gap_grid = scoreRows(:, 1);
R.score_table = array2table(scoreRows, 'VariableNames', {'gap_mm','J1','rhoJ'});
R.static_dx0 = staticState.dx0;
R.refined_gap_count = numRefine;
R.coarse_gap_count = numel(coarseIdx);
R.fine_gap_count = numel(fineIdx);
end
