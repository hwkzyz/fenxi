function staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg)
%estimate_highspeed_static_gap_raw  Hierarchical static gap and shift search.

gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
y = highMap.V_a(:);
x = highMap.x_v(:);
dxGridFull = cfg.dxStaticGrid(:)';

coarseGapGrid = linspace(gapMin, gapMax, get_cfg_field(cfg, 'gapStaticCoarseGapN', 41));
coarseDxGrid = linspace(min(dxGridFull), max(dxGridFull), get_cfg_field(cfg, 'gapStaticCoarseDxN', 41));
sampleIdx = pick_evenly_spaced_indices(numel(x), get_cfg_field(cfg, 'gapStaticTargetSamples', 4000));
xCoarse = x(sampleIdx);
yCoarse = y(sampleIdx);
JCoarse = eval_static_gap_cost_grid(xCoarse, yCoarse, templateLib, coarseGapGrid, coarseDxGrid);
[gCoarse, dxCoarse] = pick_best_gap_dx(coarseGapGrid, coarseDxGrid, JCoarse);

fineGapHalfWidth = get_cfg_field(cfg, 'gapStaticFineGapHalfWidth', 0.03);
fineDxHalfWidth = get_cfg_field(cfg, 'gapStaticFineDxHalfWidth', 0.06);
fineGapN = get_cfg_field(cfg, 'gapStaticFineGapN', 31);
fineDxN = get_cfg_field(cfg, 'gapStaticFineDxN', 31);

gapStepCoarse = coarseGapGrid(min(end, 2)) - coarseGapGrid(1);
dxStepCoarse = coarseDxGrid(min(end, 2)) - coarseDxGrid(1);
fineGapGrid = linspace(max(gapMin, gCoarse - max(fineGapHalfWidth, gapStepCoarse)), ...
    min(gapMax, gCoarse + max(fineGapHalfWidth, gapStepCoarse)), fineGapN);
fineDxGrid = linspace(max(min(dxGridFull), dxCoarse - max(fineDxHalfWidth, dxStepCoarse)), ...
    min(max(dxGridFull), dxCoarse + max(fineDxHalfWidth, dxStepCoarse)), fineDxN);
JFine = eval_static_gap_cost_grid(x, y, templateLib, fineGapGrid, fineDxGrid);
[gBest, dxBest, bestCost] = pick_best_gap_dx(fineGapGrid, fineDxGrid, JFine);

staticState = struct( ...
    'gHat', gBest, ...
    'dx0', dxBest, ...
    'coarse_gap_grid', coarseGapGrid(:), ...
    'coarse_dx_grid', coarseDxGrid(:), ...
    'fine_gap_grid', fineGapGrid(:), ...
    'fine_dx_grid', fineDxGrid(:), ...
    'coarse_best_g', gCoarse, ...
    'coarse_best_dx', dxCoarse, ...
    'best_cost', bestCost);
end

function J2D = eval_static_gap_cost_grid(x, y, templateLib, gapGrid, dxGrid)
J2D = zeros(numel(gapGrid), numel(dxGrid));
for ig = 1:numel(gapGrid)
    g = gapGrid(ig);
    for id = 1:numel(dxGrid)
        dx = dxGrid(id);
        yPred = eval_gap_template(templateLib, g, x - dx);
        valid = isfinite(yPred);
        res = y(valid) - yPred(valid);
        J2D(ig, id) = mean(res.^2);
    end
end
end

function [gBest, dxBest, bestCost] = pick_best_gap_dx(gapGrid, dxGrid, J2D)
[JStatic, dxIdx] = min(J2D, [], 2);
[bestCost, bestIdx] = min(JStatic);
gBest = gapGrid(bestIdx);
dxBest = dxGrid(dxIdx(bestIdx));
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
