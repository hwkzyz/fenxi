function precomp = build_vp_precompute(t, V, x, templateLib, gapGrid, cfg)
%build_vp_precompute  Precompute reusable VP terms for a gap grid.

t = t(:);
V = V(:);
x = x(:);
gapGrid = gapGrid(:);

precomp.gapGrid = gapGrid;
precomp.S1 = sin(2*pi*t*cfg.f1Grid(:).');
precomp.C1 = cos(2*pi*t*cfg.f1Grid(:).');
precomp.S2 = sin(2*pi*t*cfg.f2Grid(:).');
precomp.C2 = cos(2*pi*t*cfg.f2Grid(:).');

nSample = numel(x);
nGap = numel(gapGrid);
precomp.F0Mat = zeros(nSample, nGap);
precomp.FxMat = zeros(nSample, nGap);
precomp.YMat = zeros(nSample, nGap);

for ig = 1:nGap
    F0 = eval_gap_template(templateLib, gapGrid(ig), x);
    Fx = eval_gap_derivative(templateLib, gapGrid(ig), x);
    precomp.F0Mat(:, ig) = F0;
    precomp.FxMat(:, ig) = Fx;
    precomp.YMat(:, ig) = V - F0;
end
end
