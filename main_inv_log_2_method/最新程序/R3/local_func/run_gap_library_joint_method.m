function R = run_gap_library_joint_method(highMap, templateLib, cfg, staticState)
%run_gap_library_joint_method  Static gap scan + fixed-gap VARPRO + joint refine.

t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);

if nargin < 4 || isempty(staticState) || ~isfield(staticState, 'gHat')
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

rawTimer = tic;
rawFit = fit_wave_varpro(t, V, x, templateLib, staticState.gHat, 0, cfg);
rawFit = refine_wave_fit(rawFit, t, V, x, templateLib, staticState.gHat, 0, cfg);
timeRaw = toc(rawTimer);

[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);
jointTimer = tic;
joint = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
    staticState.gHat, staticState.dx0, rawFit.p, 0);
timeJoint = toc(jointTimer);

VFit = eval_gap_template(templateLib, joint.g, x - joint.dx - fit_u(joint.p, t));

R.method = "gap_library_joint";
R.g_used = joint.g;
R.dx_used = joint.dx;
R.fit = joint;
R.raw_fit = rawFit;
R.staticState = staticState;
R.f_id = joint.f;
R.A_id = joint.A;
R.phi_id = joint.phi;
R.VFit = VFit;
R.rmse = sqrt(mean((V - VFit).^2));
R.eta_g = eta_gap(templateLib, joint.g);
R.elapsed_s = timeRaw + timeJoint;
R.time_raw_fit_s = timeRaw;
R.time_joint_s = timeJoint;
end
