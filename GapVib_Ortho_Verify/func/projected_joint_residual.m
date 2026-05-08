function rAll = projected_joint_residual(theta, t, V, x, templateLib, gidx, ng, mu)
%projected_joint_residual  Full residual plus weighted shape residual.

gap = theta(1);
dx = theta(2);
p = theta(3:end);
xq = x - dx - fit_u(p, t);
rFull = V - eval_gap_template(templateLib, gap, xq);
rPerp = cell(ng, 1);
for k = 1:ng
    idx = gidx == k;
    Fx = eval_gap_derivative(templateLib, gap, xq(idx));
    rPerp{k} = remove_translation_component(rFull(idx), Fx);
end
rAll = [rFull; sqrt(mu) * vertcat(rPerp{:})];
end
