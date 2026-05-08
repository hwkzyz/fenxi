function R = run_ortho_iter_gap(highMap, templateLib, cfg, ~, initResult)
%run_ortho_iter_gap  Alternate projected gap update and vibration fitting.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

gHat = initResult.g_used;
fit = initResult.fit;
history = table();
timeGapSearch = 0;
timeVarpro = 0;
timeRefineWave = 0;

for iter = 1:3
    uSample = fit.delta_sample;
    localGrid = linspace(max(0.05, gHat - 0.04), gHat + 0.04, 81);
    localGrid = localGrid(localGrid >= min(templateLib.gapTrain) - cfg.rawGapSearchMargin & ...
        localGrid <= max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
    shapePrecomp = build_gap_shape_precompute(x - uSample, templateLib, localGrid);
    J = zeros(numel(localGrid), 1);
    timerGap = tic;
    for ig = 1:numel(localGrid)
        J(ig) = projected_shape_cost(V, x - uSample, localGrid(ig), templateLib, gidx, ng, shapePrecomp, ig);
    end
    timeGapSearch = timeGapSearch + toc(timerGap);
    [~, bestIdx] = min(J);
    gHat = localGrid(bestIdx);

    timerVarpro = tic;
    wavePrecomp = build_wave_varpro_precompute(t, V, x, templateLib, gHat, 0, cfg);
    fit = fit_wave_varpro(t, V, x, templateLib, gHat, 0, cfg, wavePrecomp);
    timeVarpro = timeVarpro + toc(timerVarpro);
    timerRefineWave = tic;
    fit = refine_wave_fit(fit, t, V, x, templateLib, gHat, 0, cfg);
    timeRefineWave = timeRefineWave + toc(timerRefineWave);
    history = [history; table(iter, gHat, fit.f(1), fit.f(2), fit.A(1), fit.A(2), ...
        'VariableNames', {'iter','gap_mm','f1_Hz','f2_Hz','A1_mm','A2_mm'})]; %#ok<AGROW>
end

muList = [0, 0.05, 0.15, 0.35, 0.75];
bestJoint = [];
bestRmse = inf;
timerJoint = tic;
for imu = 1:numel(muList)
    jointTry = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
        gHat, fit.dx, fit.p, muList(imu));
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

joint = bestJoint;
VFit = eval_gap_template(templateLib, joint.g, x - joint.dx - fit_u(joint.p, t));

R.method = "ortho_iter_gap";
R.g_used = joint.g;
R.dx_used = joint.dx;
R.fit = joint;
R.f_id = joint.f;
R.A_id = joint.A;
R.phi_id = joint.phi;
R.VFit = VFit;
R.rmse = bestRmse;
R.eta_g = eta_gap(templateLib, joint.g);
R.history = history;
R.time_gap_search_s = timeGapSearch;
R.time_varpro_s = timeVarpro;
R.time_refine_wave_s = timeRefineWave;
R.time_joint_s = timeJoint;
end
