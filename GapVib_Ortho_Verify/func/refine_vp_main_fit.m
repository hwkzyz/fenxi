function fit = refine_vp_main_fit(fit, t, V, x, templateLib, gHat, cfg)
%refine_vp_main_fit  Nonlinear refinement for VP-main candidates.

meta.J1 = getfield_with_default(fit, 'J1', NaN); %#ok<GFLD>
meta.J2 = getfield_with_default(fit, 'J2', NaN); %#ok<GFLD>
meta.deltaJ = getfield_with_default(fit, 'deltaJ', NaN); %#ok<GFLD>
meta.rhoJ = getfield_with_default(fit, 'rhoJ', NaN); %#ok<GFLD>

if isfield(fit, 'candidates') && numel(fit.candidates) > 1
    candList = fit.candidates;
    bestFit = [];
    bestCost = inf;
    for ic = 1:numel(candList)
        fitTry = refine_one_vp_fit(candList(ic), t, V, x, templateLib, gHat, cfg);
        resTry = V - eval_gap_template(templateLib, gHat, x - fitTry.delta_sample);
        costTry = dot(resTry, resTry);
        fitTry.refined_sse = costTry;
        if costTry < bestCost
            bestCost = costTry;
            bestFit = fitTry;
        end
    end
    fit = bestFit;
    fit.candidate_count = numel(candList);
    fit.J1 = meta.J1;
    fit.J2 = meta.J2;
    fit.deltaJ = meta.deltaJ;
    fit.rhoJ = meta.rhoJ;
    return;
end

fit = refine_one_vp_fit(fit, t, V, x, templateLib, gHat, cfg);
fit.J1 = meta.J1;
fit.J2 = meta.J2;
fit.deltaJ = meta.deltaJ;
fit.rhoJ = meta.rhoJ;
end

function fit = refine_one_vp_fit(fit, t, V, x, templateLib, gHat, cfg)
theta0 = [fit.dx, fit.p];
lb = [fit.dx - 0.12, 0.001, -pi, max(min(cfg.f1Grid), fit.f(1) - 20), ...
    0.001, -pi, max(min(cfg.f2Grid), fit.f(2) - 20)];
ub = [fit.dx + 0.12, 0.800, pi, min(max(cfg.f1Grid), fit.f(1) + 20), ...
    0.800, pi, min(max(cfg.f2Grid), fit.f(2) + 20)];

theta = solve_lsq_bounded(@(th) V - eval_gap_template(templateLib, gHat, ...
    x - th(1) - fit_u(th(2:end), t)), theta0, lb, ub, 120, 1e-6, 1e-6);
fit = pack_nonlinear_fit(theta(1), theta(2:end), t);
end

function val = getfield_with_default(S, name, defaultVal)
if isfield(S, name)
    val = S.(name);
else
    val = defaultVal;
end
end
