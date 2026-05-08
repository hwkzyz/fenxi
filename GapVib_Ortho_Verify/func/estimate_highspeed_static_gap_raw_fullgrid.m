function staticState = estimate_highspeed_static_gap_raw_fullgrid(highMap, templateLib, cfg)
%estimate_highspeed_static_gap_raw_fullgrid  Original exhaustive static gap and shift grid search.

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
        yPred = eval_gap_template(templateLib, g, x - dx);
        valid = isfinite(yPred);
        res = y(valid) - yPred(valid);
        J2D(ig, id) = mean(res.^2);
    end
end

[JStatic, dxIdx] = min(J2D, [], 2);
[bestCost, bestIdx] = min(JStatic);
staticState = struct( ...
    'gHat', gapQueryGrid(bestIdx), ...
    'dx0', dxGrid(dxIdx(bestIdx)), ...
    'best_cost', bestCost);
end
