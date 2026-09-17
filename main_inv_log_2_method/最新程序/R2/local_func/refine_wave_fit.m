function fit = refine_wave_fit(fit, t, V, x, templateLib, gHat, dxHat, cfg)
%refine_wave_fit  Nonlinear vibration refinement after fixed-gap VARPRO.

if isfield(fit, 'candidates') && numel(fit.candidates) > 1
    candList = fit.candidates;
    bestFit = [];
    bestCost = inf;
    for ic = 1:numel(candList)
        fitTry = refine_one_wave_fit(candList(ic), t, V, x, templateLib, gHat, dxHat, cfg);
        resTry = V - eval_gap_template(templateLib, gHat, ...
            x - dxHat - fitTry.delta_sample);
        costTry = dot(resTry, resTry);
        fitTry.refined_sse = costTry;
        if costTry < bestCost
            bestCost = costTry;
            bestFit = fitTry;
        end
    end
    fit = bestFit;
    fit.candidate_count = numel(candList);
    return;
end

fit = refine_one_wave_fit(fit, t, V, x, templateLib, gHat, dxHat, cfg);
end

function fit = refine_one_wave_fit(fit, t, V, x, templateLib, gHat, dxHat, cfg)
theta0 = [fit.dx, fit.p];
lb = [fit.dx - 0.12, 0.001, -pi, max(min(cfg.f1Grid), fit.f(1) - 20), ...
    0.001, -pi, max(min(cfg.f2Grid), fit.f(2) - 20)];
ub = [fit.dx + 0.12, 0.800, pi, min(max(cfg.f1Grid), fit.f(1) + 20), ...
    0.800, pi, min(max(cfg.f2Grid), fit.f(2) + 20)];

theta = solve_lsq_bounded(@(th) V - eval_gap_template(templateLib, gHat, ...
    x - dxHat - th(1) - fit_u(th(2:end), t)), theta0, lb, ub, 300);
fit = pack_nonlinear_fit(theta(1), theta(2:end), t);
end
