function [orthoResult, orthoIterResult, summaryTable] = RunOrthoMethodsOnMap(highMap, templateLib, cfg, truth, staticState)
%RunOrthoMethodsOnMap  Compute orthogonal gap-vibration methods on one map.

baseTimer = tic;
orthoResult = run_ortho_gap(highMap, templateLib, cfg, staticState);
orthoResult.elapsed_s = toc(baseTimer);

iterTimer = tic;
orthoIterResult = run_ortho_iter_gap(highMap, templateLib, cfg, staticState, orthoResult);
orthoIterResult.elapsed_s = toc(iterTimer);

summaryTable = [make_ortho_summary(orthoResult, truth); ...
    make_ortho_summary(orthoIterResult, truth)];
end
