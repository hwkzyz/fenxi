function baseState = estimate_projected_gap_only(highMap, templateLib, cfg, staticState)
%estimate_projected_gap_only  Hierarchical projected-gap center without vibration fit.

x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

dxHat = staticState.dx0;
gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gAnchor = staticState.gHat;
if isempty(gAnchor) || ~isfinite(gAnchor)
    gAnchor = 0.5 * (gapMin + gapMax);
end

coarseHalfWidth = get_cfg_field(cfg, 'gapProjHalfWidth', 0.08);
coarseN = get_cfg_field(cfg, 'gapProjCoarseN', 61);
fineHalfWidth = get_cfg_field(cfg, 'gapProjFineHalfWidth', 0.02);
fineN = get_cfg_field(cfg, 'gapProjFineN', 31);
projLambda = get_cfg_field(cfg, 'gapProjectionLambda', 1);

coarseGapGrid = linspace(max(gapMin, gAnchor - coarseHalfWidth), ...
    min(gapMax, gAnchor + coarseHalfWidth), coarseN);
sampleIdx = pick_evenly_spaced_indices(numel(x), get_cfg_field(cfg, 'gapProjTargetSamples', 5000));
xCoarse = x(sampleIdx);
VCoarse = V(sampleIdx);
[~, ~, gidxCoarse] = unique([highMap.rev_v(sampleIdx), highMap.S_v(sampleIdx)], 'rows');
ngCoarse = max(gidxCoarse);
shapeCoarse = build_gap_shape_precompute(xCoarse - dxHat, templateLib, coarseGapGrid);
JCoarse = zeros(numel(coarseGapGrid), 1);
for ig = 1:numel(coarseGapGrid)
    JCoarse(ig) = projected_shape_cost(VCoarse, xCoarse - dxHat, coarseGapGrid(ig), ...
        templateLib, gidxCoarse, ngCoarse, shapeCoarse, ig, projLambda);
end
[~, iCoarseBest] = min(JCoarse);
gCoarse = coarseGapGrid(iCoarseBest);

coarseStep = coarseGapGrid(min(end, 2)) - coarseGapGrid(1);
fineGapGrid = linspace(max(gapMin, gCoarse - max(fineHalfWidth, coarseStep)), ...
    min(gapMax, gCoarse + max(fineHalfWidth, coarseStep)), fineN);
shapeFine = build_gap_shape_precompute(x - dxHat, templateLib, fineGapGrid);
JFine = zeros(numel(fineGapGrid), 1);
for ig = 1:numel(fineGapGrid)
    JFine(ig) = projected_shape_cost(V, x - dxHat, fineGapGrid(ig), ...
        templateLib, gidx, ng, shapeFine, ig, projLambda);
end
[~, igBest] = min(JFine);

baseState = struct();
baseState.g_used = fineGapGrid(igBest);
baseState.dx_used = dxHat;
baseState.gap_grid = fineGapGrid(:);
baseState.J = JFine;
baseState.coarse_gap_grid = coarseGapGrid(:);
baseState.coarse_J = JCoarse;
baseState.coarse_best_g = gCoarse;
baseState.projection_lambda = projLambda;
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
