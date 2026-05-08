function R = run_main_vp_joint(highMap, templateLib, cfg, staticState)
%run_main_vp_joint
% Fast static gap init -> VP-style local gap screening ->
% fixed-gap nonlinear refinement -> final projected joint refinement.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

if nargin < 4 || isempty(staticState) || ~isfield(staticState, 'gHat')
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

gCenter = staticState.gHat;
gapHalfWidth = get_cfg_field(cfg, 'mainVpGapHalfWidth', 0.06);
gapGrid = linspace(max(0.05, gCenter - gapHalfWidth), gCenter + gapHalfWidth, ...
    get_cfg_field(cfg, 'mainVpGapN', 61));
gapGrid = gapGrid(gapGrid >= min(templateLib.gapTrain) - cfg.rawGapSearchMargin & ...
    gapGrid <= max(templateLib.gapTrain) + cfg.rawGapSearchMargin);

timerPrecomp = tic;
precomp = build_vp_precompute(t, V, x, templateLib, gapGrid, cfg);
timePrecomp = toc(timerPrecomp);

bestFit = [];
bestG = gapGrid(1);
bestSSE = inf;
coarseCount = min(get_cfg_field(cfg, 'mainVpCoarseCount', 21), numel(gapGrid));
coarseIdx = round(linspace(1, numel(gapGrid), coarseCount));
coarseIdx = unique(coarseIdx(:));
coarseRows = zeros(numel(coarseIdx), 3);

timerScreen = tic;
for ic = 1:numel(coarseIdx)
    ig = coarseIdx(ic);
    gTry = gapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    coarseRows(ic, :) = [gTry, fitTry.J1, fitTry.rhoJ];
end

[~, coarseOrder] = sort(coarseRows(:, 2), 'ascend');
bestCoarseJ = coarseRows(coarseOrder(1), 2);
keepCoarse = min(get_cfg_field(cfg, 'mainVpKeepCoarse', 5), numel(coarseOrder));
relTol = get_cfg_field(cfg, 'mainVpRelTol', 0.02);
extraKeep = coarseOrder(coarseRows(coarseOrder, 2) <= bestCoarseJ * (1 + relTol));
keepList = unique([coarseOrder(1:keepCoarse); extraKeep(:)], 'stable');
keepList = keepList(1:min(numel(keepList), get_cfg_field(cfg, 'mainVpMaxKeep', 7)));
fineMask = false(numel(gapGrid), 1);
halfWindow = get_cfg_field(cfg, 'mainVpFineHalfWindow', 5);
for ic = 1:numel(keepList)
    ig0 = coarseIdx(keepList(ic));
    lo = max(1, ig0 - halfWindow);
    hi = min(numel(gapGrid), ig0 + halfWindow);
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
timeScreen = toc(timerScreen);

[~, order] = sort(scoreRows(:, 2), 'ascend');
numRefine = min(get_cfg_field(cfg, 'mainVpRefineCount', 5), numel(order));
timerRefineWave = tic;
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
timeRefineWave = toc(timerRefineWave);

if isempty(bestFit)
    jf = order(1);
    gTry = scoreRows(jf, 1);
    fitTry = fitBank{jf};
    bestFit = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    bestG = gTry;
end

muList = get_cfg_field(cfg, 'mainUseProjectedJointMuList', [0, 0.15, 0.75]);
bestJoint = [];
bestRmse = inf;
timerJoint = tic;
for imu = 1:numel(muList)
    jointTry = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
        bestG, bestFit.dx, bestFit.p, muList(imu));
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

R.method = "main_vp_joint";
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
R.J = scoreRows(:, 2);
R.gap_grid = scoreRows(:, 1);
R.score_table = array2table(scoreRows, 'VariableNames', {'gap_mm','J1','rhoJ'});
R.static_dx0 = staticState.dx0;
R.refined_gap_count = numRefine;
R.coarse_gap_count = numel(coarseIdx);
R.fine_gap_count = numel(fineIdx);
R.time_precomp_s = timePrecomp;
R.time_screen_s = timeScreen;
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
