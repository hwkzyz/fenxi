function joint = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, g0, dx0, p0, mu)
%refine_projected_joint  Joint nonlinear update with projected residuals.

theta0 = [g0, dx0, p0];
dxHalfWidth = get_cfg_field(cfg, 'projectedJointDxHalfWidth', 0.05);
gapHalfWidth = get_cfg_field(cfg, 'projectedJointGapHalfWidth', 0.05);
ampScaleLower = get_cfg_field(cfg, 'nonlinearAmpScaleLower', 0.5);
ampScaleUpper = get_cfg_field(cfg, 'nonlinearAmpScaleUpper', 1.5);
ampUpperBound = get_cfg_field(cfg, 'nonlinearAmpUpperBound', 0.8);
ampMinimumUpper = get_cfg_field(cfg, 'nonlinearAmpMinimumUpperBound', 0);
freqHalfWidth = get_cfg_field(cfg, 'nonlinearFreqHalfWidth', 20);
lb = [max(0.05, g0 - gapHalfWidth), dx0 - dxHalfWidth, ...
    max(0.001, ampScaleLower * p0(1)), -pi, max(min(cfg.f1Grid), p0(3) - freqHalfWidth), ...
    max(0.001, ampScaleLower * p0(4)), -pi, max(min(cfg.f2Grid), p0(6) - freqHalfWidth)];
ub = [g0 + gapHalfWidth, dx0 + dxHalfWidth, ...
    min(ampUpperBound, max(ampScaleUpper * p0(1), ampMinimumUpper)), pi, min(max(cfg.f1Grid), p0(3) + freqHalfWidth), ...
    min(ampUpperBound, max(ampScaleUpper * p0(4), ampMinimumUpper)), pi, min(max(cfg.f2Grid), p0(6) + freqHalfWidth)];
if ~(isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement)
    lb(1) = max(lb(1), min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
    ub(1) = min(ub(1), max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
end

[theta, solveInfo] = solve_lsq_bounded(@(th) projected_joint_residual(th, t, V, x, ...
    templateLib, gidx, ng, mu), theta0, lb, ub, ...
    get_cfg_field(cfg, 'projectedJointMaxIter', 500));
joint = pack_nonlinear_fit(theta(2), theta(3:end), t);
joint.g = theta(1);
joint.dx = theta(2);
joint.theta = theta;
joint.solve_info = solveInfo;
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
