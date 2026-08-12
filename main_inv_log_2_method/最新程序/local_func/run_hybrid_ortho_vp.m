function R = run_hybrid_ortho_vp(highMap, templateLib, cfg, staticState)
%run_hybrid_ortho_vp  Orthogonal gap narrowing followed by exhaustive VP search.

timerBase = tic;
base = estimate_projected_gap_only(highMap, templateLib, cfg, staticState);
timeBase = toc(timerBase);
t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);

localGrid = linspace(max(0.05, base.g_used - 0.04), base.g_used + 0.04, 81);
localGrid = localGrid(localGrid >= min(templateLib.gapTrain) - cfg.rawGapSearchMargin & ...
    localGrid <= max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
timerPrecomp = tic;
precomp = build_vp_precompute(t, V, x, templateLib, localGrid, cfg);
timePrecomp = toc(timerPrecomp);

bestFit = [];
bestG = localGrid(1);
bestSSE = inf;
scoreRows = zeros(numel(localGrid), 3);
timeFitVp = 0;
timeRefine = 0;
timeEval = 0;
refineCandidateKeep = 1;

for ig = 1:numel(localGrid)
    gTry = localGrid(ig);
    timerFitVp = tic;
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    timeFitVp = timeFitVp + toc(timerFitVp);
    if isfield(fitTry, 'candidates') && numel(fitTry.candidates) > refineCandidateKeep
        fitTry.candidates = fitTry.candidates(1:refineCandidateKeep);
    end
    timerRefine = tic;
    fitTry = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    timeRefine = timeRefine + toc(timerRefine);
    timerEval = tic;
    resTry = V - eval_gap_template(templateLib, gTry, x - fitTry.delta_sample);
    timeEval = timeEval + toc(timerEval);
    sseTry = dot(resTry, resTry);
    scoreRows(ig, :) = [gTry, fitTry.J1, fitTry.rhoJ];
    if sseTry < bestSSE
        bestSSE = sseTry;
        bestFit = fitTry;
        bestG = gTry;
    end
end

VFit = eval_gap_template(templateLib, bestG, x - bestFit.delta_sample);
R.method = "hybrid_ortho_vp";
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
R.gap_grid = localGrid(:);
R.score_table = array2table(scoreRows, 'VariableNames', {'gap_mm','J1','rhoJ'});
R.base_result = base;
R.coarse_gap_count = numel(localGrid);
R.fine_gap_count = numel(localGrid);
R.refined_gap_count = numel(localGrid);
R.refine_candidate_keep = refineCandidateKeep;
R.time_base_s = timeBase;
R.time_precomp_s = timePrecomp;
R.time_fit_vp_s = timeFitVp;
R.time_refine_s = timeRefine;
R.time_eval_s = timeEval;
end
