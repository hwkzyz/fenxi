function gapState = estimate_gap_init_smooth_projected(highMap, templateLib, cfg, staticState)
%estimate_gap_init_smooth_projected  Gap correction with smooth unknown passage shifts.
%
% This estimator does not require vibration-frequency priors. For each
% candidate gap, it models the first-order vibration contamination as an
% unknown passage-wise translation sequence and regularizes that sequence
% to vary smoothly across passages.

if nargin < 4 || isempty(staticState)
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

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
    JCoarse(ig) = smooth_projected_shape_cost(VCoarse, xCoarse - dxHat, ...
        coarseGapGrid(ig), templateLib, gidxCoarse, ngCoarse, shapeCoarse, ig, cfg);
end
[~, iCoarseBest] = min(JCoarse);
gCoarse = coarseGapGrid(iCoarseBest);

coarseStep = coarseGapGrid(min(end, 2)) - coarseGapGrid(1);
fineGapGrid = linspace(max(gapMin, gCoarse - max(fineHalfWidth, coarseStep)), ...
    min(gapMax, gCoarse + max(fineHalfWidth, coarseStep)), fineN);
shapeFine = build_gap_shape_precompute(x - dxHat, templateLib, fineGapGrid);
JFine = zeros(numel(fineGapGrid), 1);
deltaBank = cell(numel(fineGapGrid), 1);
for ig = 1:numel(fineGapGrid)
    [JFine(ig), deltaBank{ig}] = smooth_projected_shape_cost(V, x - dxHat, ...
        fineGapGrid(ig), templateLib, gidx, ng, shapeFine, ig, cfg);
end
[bestCost, iBest] = min(JFine);

gapState = struct();
gapState.gHat = fineGapGrid(iBest);
gapState.g_used = fineGapGrid(iBest);
gapState.dx0 = dxHat;
gapState.dx_used = dxHat;
gapState.g_static = staticState.gHat;
gapState.dx_static = staticState.dx0;
gapState.gap_grid = fineGapGrid(:);
gapState.J = JFine;
gapState.coarse_gap_grid = coarseGapGrid(:);
gapState.coarse_J = JCoarse;
gapState.coarse_best_g = gCoarse;
gapState.best_cost = bestCost;
gapState.delta_passage = deltaBank{iBest};
gapState.staticState = staticState;
gapState.smoothness = get_cfg_field(cfg, 'gapSmoothProjSmoothness', 1);
gapState.ridge = get_cfg_field(cfg, 'gapSmoothProjRidge', 0.01);
end

function [cost, delta] = smooth_projected_shape_cost(V, xq, gap, templateLib, gidx, ng, precomp, gapIdx, cfg)
if nargin >= 8 && ~isempty(precomp) && ~isempty(gapIdx)
    FAll = precomp.FMat(:, gapIdx);
    FxAll = precomp.FxMat(:, gapIdx);
else
    FAll = eval_gap_template(templateLib, gap, xq);
    FxAll = eval_gap_derivative(templateLib, gap, xq);
end

a = zeros(ng, 1);
b = zeros(ng, 1);
c = zeros(ng, 1);
count = 0;
for k = 1:ng
    idx = gidx == k;
    r = V(idx) - FAll(idx);
    Fx = FxAll(idx);
    valid = isfinite(r) & isfinite(Fx);
    r = r(valid);
    Fx = Fx(valid);
    a(k) = dot(Fx, Fx);
    b(k) = dot(Fx, r);
    c(k) = dot(r, r);
    count = count + numel(r);
end

scale = median(a(a > 0));
if isempty(scale) || ~isfinite(scale)
    scale = max(mean(a), eps);
end
ridge = get_cfg_field(cfg, 'gapSmoothProjRidge', 0.01) * scale;
smoothness = get_cfg_field(cfg, 'gapSmoothProjSmoothness', 1) * scale;

H = diag(a + ridge);
if ng > 1 && smoothness > 0
    e = ones(ng, 1);
    D = spdiags([-e, e], [0, 1], ng - 1, ng);
    H = H + smoothness * full(D' * D);
end
delta = H \ b;

residualCost = sum(c - 2 * delta .* b + (delta.^2) .* a);
ridgeCost = ridge * sum(delta.^2);
smoothCost = 0;
if ng > 1 && smoothness > 0
    smoothCost = smoothness * sum(diff(delta).^2);
end
cost = (residualCost + ridgeCost + smoothCost) / max(count, 1);
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
