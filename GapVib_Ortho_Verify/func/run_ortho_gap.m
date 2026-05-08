function R = run_ortho_gap(highMap, templateLib, cfg, staticState)
%run_ortho_gap  Estimate gap from translation-orthogonal shape residuals.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

gapGrid = linspace(max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin), ...
    max(templateLib.gapTrain) + cfg.rawGapSearchMargin, cfg.gapQueryN);
dxHat = staticState.dx0;
J = zeros(numel(gapGrid), 1);

for ig = 1:numel(gapGrid)
    J(ig) = projected_shape_cost(V, x - dxHat, gapGrid(ig), templateLib, gidx, ng);
end

[~, igBest] = min(J);
gHat = gapGrid(igBest);
dxFitAnchor = 0;
fit = fit_wave_varpro(t, V, x, templateLib, gHat, dxFitAnchor, cfg);
fit = refine_wave_fit(fit, t, V, x, templateLib, gHat, dxFitAnchor, cfg);

uSample = fit.delta_sample;
VFit = eval_gap_template(templateLib, gHat, x - dxFitAnchor - uSample);

R.method = "ortho_gap";
R.g_used = gHat;
R.dx_used = dxFitAnchor + fit.dx;
R.gap_grid = gapGrid(:);
R.J = J;
R.fit = fit;
R.f_id = fit.f;
R.A_id = fit.A;
R.phi_id = fit.phi;
R.VFit = VFit;
R.rmse = sqrt(mean((V - VFit).^2));
R.eta_g = eta_gap(templateLib, gHat);
end
