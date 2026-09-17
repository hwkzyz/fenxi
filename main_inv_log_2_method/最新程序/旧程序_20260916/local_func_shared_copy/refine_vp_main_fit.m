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
        resTry(~isfinite(resTry)) = 10 * max(std(V), 1e-3);
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
dxHalfWidth = getfield_with_default(cfg, 'vpDxHalfWidth', 0.12);
ampUpperBound = getfield_with_default(cfg, 'nonlinearAmpUpperBound', 0.8); %#ok<GFLD>
freqHalfWidth = getfield_with_default(cfg, 'nonlinearFreqHalfWidth', 20); %#ok<GFLD>
ampSeed = getfield_with_default(cfg, 'nonlinearAmpSeedMm', 0.12 * ampUpperBound); %#ok<GFLD>
pSeed = fit.p;
pSeed([1,4]) = max(pSeed([1,4]), ampSeed);
theta0 = [fit.dx, pSeed];
lb = [fit.dx - dxHalfWidth, 0.001, -pi, max(min(cfg.f1Grid), fit.f(1) - freqHalfWidth), ...
    0.001, -pi, max(min(cfg.f2Grid), fit.f(2) - freqHalfWidth)];
ub = [fit.dx + dxHalfWidth, ampUpperBound, pi, min(max(cfg.f1Grid), fit.f(1) + freqHalfWidth), ...
    ampUpperBound, pi, min(max(cfg.f2Grid), fit.f(2) + freqHalfWidth)];
% The VP amplitudes are derivatives-based linearization coefficients.  They
% can be much smaller than the displacement amplitude when the gap response
% is nonlinear, so they must not define the nonlinear search interval.
ampMinimumUpper = getfield_with_default(cfg, 'nonlinearAmpMinimumUpperBound', 0); %#ok<GFLD>
ub(2) = max(ub(2), ampMinimumUpper);
ub(5) = max(ub(5), ampMinimumUpper);

[theta, solveInfo] = solve_lsq_bounded(@(th) vp_refine_residual(th, t, V, x, ...
    templateLib, gHat), theta0, lb, ub, 120, 1e-6, 1e-6);
fit = pack_nonlinear_fit(theta(1), theta(2:end), t);
fit.solve_info = solveInfo;
end

function r = vp_refine_residual(th, t, V, x, templateLib, gHat)
y = eval_gap_template(templateLib, gHat, x - th(1) - fit_u(th(2:end), t));
r = V - y;
r(~isfinite(r)) = 10 * max(std(V), 1e-3);
end

function val = getfield_with_default(S, name, defaultVal)
if isfield(S, name)
    val = S.(name);
else
    val = defaultVal;
end
end
