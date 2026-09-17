function precomp = build_gap_shape_precompute(xq, templateLib, gapGrid)
%build_gap_shape_precompute  Precompute template and derivative on a gap grid.

xq = xq(:);
gapGrid = gapGrid(:);
nSample = numel(xq);
nGap = numel(gapGrid);

precomp.gapGrid = gapGrid;
precomp.FMat = zeros(nSample, nGap);
precomp.FxMat = zeros(nSample, nGap);

for ig = 1:nGap
    precomp.FMat(:, ig) = eval_gap_template(templateLib, gapGrid(ig), xq);
    precomp.FxMat(:, ig) = eval_gap_derivative(templateLib, gapGrid(ig), xq);
end
end
