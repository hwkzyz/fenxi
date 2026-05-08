function joint = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, g0, dx0, p0, mu)
%refine_projected_joint  Joint nonlinear update with projected residuals.

theta0 = [g0, dx0, p0];
dxHalfWidth = get_cfg_field(cfg, 'projectedJointDxHalfWidth', 0.05);
lb = [max(0.05, g0 - 0.05), dx0 - dxHalfWidth, ...
    max(0.001, 0.5*p0(1)), -pi, max(min(cfg.f1Grid), p0(3) - 20), ...
    max(0.001, 0.5*p0(4)), -pi, max(min(cfg.f2Grid), p0(6) - 20)];
ub = [g0 + 0.05, dx0 + dxHalfWidth, ...
    min(0.8, 1.5*p0(1)), pi, min(max(cfg.f1Grid), p0(3) + 20), ...
    min(0.8, 1.5*p0(4)), pi, min(max(cfg.f2Grid), p0(6) + 20)];
lb(1) = max(lb(1), min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
ub(1) = min(ub(1), max(templateLib.gapTrain) + cfg.rawGapSearchMargin);

theta = solve_lsq_bounded(@(th) projected_joint_residual(th, t, V, x, ...
    templateLib, gidx, ng, mu), theta0, lb, ub, 500);
joint = pack_nonlinear_fit(theta(2), theta(3:end), t);
joint.g = theta(1);
joint.dx = theta(2);
joint.theta = theta;
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
