function R = run_confidence_gap_vp_joint(highMap, templateLib, cfg, staticState)
%run_confidence_gap_vp_joint
% Single static gap initialization -> local VP screening ->
% fixed-gap nonlinear refinement -> final refine_projected_joint.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

if nargin < 4 || isempty(staticState) || ~isfield(staticState, 'gHat')
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

timerGap = tic;
gStatic = staticState.gHat;
gCenter = gStatic;
gHalfWidth = cfg.confGapLocalHalfWidth;
gMin = gCenter - gHalfWidth;
gMax = gCenter + gHalfWidth;
gMin = max(0.05, gMin);
gMin = max(gMin, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gMax = min(gMax, max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
coarseGapGrid = linspace(gMin, gMax, cfg.confGapCoarseN);
timeGapInit = toc(timerGap);

timerPrecomp = tic;
precompCoarse = build_vp_precompute(t, V, x, templateLib, coarseGapGrid, cfg);
timePrecomp = toc(timerPrecomp);

fitBankCoarse = cell(numel(coarseGapGrid), 1);
scoreRowsCoarse = zeros(numel(coarseGapGrid), 5);
timerVp = tic;
for ig = 1:numel(coarseGapGrid)
    gTry = coarseGapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precompCoarse, ig);
    fitBankCoarse{ig} = fitTry;
    scoreRowsCoarse(ig, :) = [gTry, fitTry.J1, fitTry.J2, fitTry.deltaJ, fitTry.rhoJ];
end

[~, orderJ] = sort(scoreRowsCoarse(:, 2), 'ascend');
keepN = min(cfg.confTopGapKeep, numel(orderJ));
keepIdxCoarse = orderJ(1:keepN);

gRefMin = max(gMin, min(coarseGapGrid(keepIdxCoarse)) - 0.5 * (coarseGapGrid(2) - coarseGapGrid(1)));
gRefMax = min(gMax, max(coarseGapGrid(keepIdxCoarse)) + 0.5 * (coarseGapGrid(2) - coarseGapGrid(1)));
localGapGrid = linspace(gRefMin, gRefMax, cfg.confGapConsensusN);
precompFine = build_vp_precompute(t, V, x, templateLib, localGapGrid, cfg);

fitBank = cell(numel(localGapGrid), 1);
scoreRows = zeros(numel(localGapGrid), 5);
for ig = 1:numel(localGapGrid)
    gTry = localGapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precompFine, ig);
    fitBank{ig} = fitTry;
    scoreRows(ig, :) = [gTry, fitTry.J1, fitTry.J2, fitTry.deltaJ, fitTry.rhoJ];
end
timeVp = toc(timerVp);

refineMask = false(numel(localGapGrid), 1);
halfWindow = get_cfg_field(cfg, 'confRefineGapHalfWindow', 1);
for ik = 1:numel(keepIdxCoarse)
    g0 = coarseGapGrid(keepIdxCoarse(ik));
    [~, iNear] = min(abs(localGapGrid - g0));
    lo = max(1, iNear - halfWindow);
    hi = min(numel(localGapGrid), iNear + halfWindow);
    refineMask(lo:hi) = true;
end
refineIdx = find(refineMask);
if isempty(refineIdx)
    refineIdx = 1:numel(localGapGrid);
end

timerRefineWave = tic;
bestFit = [];
bestGap = localGapGrid(refineIdx(1));
bestSSE = inf;
refineRows = zeros(numel(refineIdx), 3);
for ir = 1:numel(refineIdx)
    ig = refineIdx(ir);
    gTry = localGapGrid(ig);
    fitTry = fitBank{ig};
    fitTry = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    resTry = V - eval_gap_template(templateLib, gTry, x - fitTry.delta_sample);
    sseTry = dot(resTry, resTry);
    refineRows(ir, :) = [gTry, sseTry, scoreRows(ig, 5)];
    if sseTry < bestSSE
        bestSSE = sseTry;
        bestFit = fitTry;
        bestGap = gTry;
    end
end
timeRefineWave = toc(timerRefineWave);

muList = cfg.confUseProjectedJointMuList;
bestJoint = [];
bestRmse = inf;
timerJoint = tic;
for imu = 1:numel(muList)
    jointTry = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
        bestGap, bestFit.dx, bestFit.p, muList(imu));
    VTry = eval_gap_template(templateLib, jointTry.g, ...
        x - jointTry.dx - fit_u(jointTry.p, t));
    rmseTry = sqrt(mean((V - VTry).^2));
    jointTry.mu = muList(imu);
    jointTry.wave_rmse = rmseTry;
    if rmseTry < bestRmse
        bestRmse = rmseTry;
        bestJoint = jointTry;
    end
end
timeJoint = toc(timerJoint);

VFit = eval_gap_template(templateLib, bestJoint.g, x - bestJoint.dx - fit_u(bestJoint.p, t));

R.method = "conf_gap_vp_joint";
R.g_used = bestJoint.g;
R.dx_used = bestJoint.dx;
R.fit = bestJoint;
R.f_id = bestJoint.f;
R.A_id = bestJoint.A;
R.phi_id = bestJoint.phi;
R.VFit = VFit;
R.rmse = bestRmse;
R.eta_g = eta_gap(templateLib, bestJoint.g);
R.deltaJ = bestFit.deltaJ;
R.rhoJ = bestFit.rhoJ;
R.local_gap_grid = localGapGrid(:);
R.coarse_gap_grid = coarseGapGrid(:);
R.coarse_score_table = array2table(scoreRowsCoarse, ...
    'VariableNames', {'gap_mm','J1','J2','deltaJ','rhoJ'});
R.score_table = array2table(scoreRows, ...
    'VariableNames', {'gap_mm','J1','J2','deltaJ','rhoJ'});
R.consensus_keep_idx = keepIdxCoarse;
R.refine_gap_idx = refineIdx;
R.refine_score_table = array2table(refineRows, ...
    'VariableNames', {'gap_mm','refined_sse','rhoJ'});
R.g_static = gStatic;
R.time_gap_init_s = timeGapInit;
R.time_precomp_s = timePrecomp;
R.time_vp_screen_s = timeVp;
R.time_refine_wave_s = timeRefineWave;
R.time_joint_s = timeJoint;
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
