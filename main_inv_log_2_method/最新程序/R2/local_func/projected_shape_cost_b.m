function cost = projected_shape_cost_b(V, xq, gap, templateLib, gidx, ng)
%projected_shape_cost_b  Cost after removing per-passage bias and translation.

cost = 0;
count = 0;
for k = 1:ng
    idx = gidx == k;
    F = eval_gap_template(templateLib, gap, xq(idx));
    Fx = eval_gap_derivative(templateLib, gap, xq(idx));
    B = [ones(nnz(idx), 1), Fx];
    rPerp = remove_component_basis(V(idx) - F, B);
    cost = cost + dot(rPerp, rPerp);
    count = count + nnz(idx);
end
cost = cost / max(count, 1);
end
