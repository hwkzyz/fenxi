function gapState = estimate_gap_init_static_projected(highMap, templateLib, cfg)
%estimate_gap_init_static_projected  Static registration followed by projected gap correction.
%
% Stage 1 estimates a complete geometric initialization (g, dx) using the
% hierarchical static template registration. Stage 2 keeps dx fixed and
% refines g with a projected residual that suppresses the first-order
% translation-like vibration component.

staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
projState = estimate_projected_gap_only(highMap, templateLib, cfg, staticState);

gapState = struct();
gapState.gHat = projState.g_used;
gapState.dx0 = staticState.dx0;
gapState.g_static = staticState.gHat;
gapState.g_projected = projState.g_used;
gapState.dx_static = staticState.dx0;
gapState.staticState = staticState;
gapState.projectedState = projState;
gapState.gap_grid = projState.gap_grid;
gapState.J = projState.J;
gapState.coarse_gap_grid = projState.coarse_gap_grid;
gapState.coarse_J = projState.coarse_J;
gapState.best_cost = min(projState.J);
end
