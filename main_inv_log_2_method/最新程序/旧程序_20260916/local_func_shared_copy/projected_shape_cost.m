function cost = projected_shape_cost(V, xq, gap, templateLib, gidx, ng, precomp, gapIdx, lambda)
%projected_shape_cost  Cost after removing each passage's translation mode.

if nargin < 9 || isempty(lambda)
    lambda = 1;
end

if nargin >= 8 && ~isempty(precomp) && ~isempty(gapIdx)
    FAll = precomp.FMat(:, gapIdx);
    FxAll = precomp.FxMat(:, gapIdx);
else
    FAll = eval_gap_template(templateLib, gap, xq);
    FxAll = eval_gap_derivative(templateLib, gap, xq);
end

cost = 0;
count = 0;
for k = 1:ng
    idx = gidx == k;
    rPerp = remove_translation_component(V(idx) - FAll(idx), FxAll(idx), lambda);
    cost = cost + dot(rPerp, rPerp);
    count = count + nnz(idx);
end
cost = cost / max(count, 1);
end
