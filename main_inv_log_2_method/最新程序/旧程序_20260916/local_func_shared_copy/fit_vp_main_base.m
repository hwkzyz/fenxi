function fit = fit_vp_main_base(t, V, x, templateLib, gHat, cfg, precomp, gapIdx)
%fit_vp_main  Variable projection over frequencies under fixed gap.

numCandidates = max(cfg.numVarproCandidates, 10);
maxKeep = max(numCandidates * 4, numCandidates);
candSSE = inf(maxKeep, 1);
candBeta = zeros(maxKeep, 5);
candF = zeros(maxKeep, 2);
if nargin >= 8 && ~isempty(precomp) && ~isempty(gapIdx)
    F0 = precomp.F0Mat(:, gapIdx);
    Fx = precomp.FxMat(:, gapIdx);
    y = precomp.YMat(:, gapIdx);
    S1 = precomp.S1;
    C1 = precomp.C1;
    S2 = precomp.S2;
    C2 = precomp.C2;
else
    F0 = eval_gap_template(templateLib, gHat, x);
    Fx = eval_gap_derivative(templateLib, gHat, x);
    y = V - F0;
    S1 = sin(2*pi*t*cfg.f1Grid(:).');
    C1 = cos(2*pi*t*cfg.f1Grid(:).');
    S2 = sin(2*pi*t*cfg.f2Grid(:).');
    C2 = cos(2*pi*t*cfg.f2Grid(:).');
end

w = Fx.^2;
fy = -Fx .* y;
yNorm2 = dot(y, y);

g11 = sum(w);
g1s1 = w.' * S1;
g1c1 = w.' * C1;
g1s2 = w.' * S2;
g1c2 = w.' * C2;

gs1s1 = sum(w .* (S1.^2), 1);
gs1c1 = sum(w .* (S1 .* C1), 1);
gc1c1 = sum(w .* (C1.^2), 1);

gs2s2 = sum(w .* (S2.^2), 1);
gs2c2 = sum(w .* (S2 .* C2), 1);
gc2c2 = sum(w .* (C2.^2), 1);

gs1s2 = S1.' * (w .* S2);
gs1c2 = S1.' * (w .* C2);
gc1s2 = C1.' * (w .* S2);
gc1c2 = C1.' * (w .* C2);

r1 = sum(fy);
rs1 = fy.' * S1;
rc1 = fy.' * C1;
rs2 = fy.' * S2;
rc2 = fy.' * C2;

for i1 = 1:numel(cfg.f1Grid)
    for i2 = 1:numel(cfg.f2Grid)
        f1 = cfg.f1Grid(i1);
        f2 = cfg.f2Grid(i2);
        % The two sinusoidal components are exchangeable when the two grids
        % are identical.  Scan one triangular half only when explicitly
        % enabled; the default preserves legacy behavior for asymmetric grids.
        if isfield(cfg, 'vpUniqueUnorderedPairs') && cfg.vpUniqueUnorderedPairs && ...
                numel(cfg.f1Grid) == numel(cfg.f2Grid) && isequal(cfg.f1Grid(:), cfg.f2Grid(:)) && f2 < f1
            continue;
        end
        if abs(f2 - f1) < 20
            continue;
        end
        G = [g11,           g1s1(i1),        g1c1(i1),        g1s2(i2),        g1c2(i2); ...
             g1s1(i1),      gs1s1(i1),       gs1c1(i1),       gs1s2(i1, i2),   gs1c2(i1, i2); ...
             g1c1(i1),      gs1c1(i1),       gc1c1(i1),       gc1s2(i1, i2),   gc1c2(i1, i2); ...
             g1s2(i2),      gs1s2(i1, i2),   gc1s2(i1, i2),   gs2s2(i2),       gs2c2(i2); ...
             g1c2(i2),      gs1c2(i1, i2),   gc1c2(i1, i2),   gs2c2(i2),       gc2c2(i2)];
        rhs = [r1; rs1(i1); rc1(i1); rs2(i2); rc2(i2)];
        beta = G \ rhs;
        sse = max(yNorm2 - dot(rhs, beta), 0);
        [worstSSE, worstIdx] = max(candSSE);
        if sse < worstSSE
            candSSE(worstIdx) = sse;
            candBeta(worstIdx, :) = beta(:).';
            candF(worstIdx, :) = [f1, f2];
        end
    end
end

[candSSE, order] = sort(candSSE, 'ascend');
candBeta = candBeta(order, :);
candF = candF(order, :);
keep = isfinite(candSSE);
candSSE = candSSE(keep);
candBeta = candBeta(keep, :);
candF = candF(keep, :);

uniqueRows = true(size(candSSE));
diversityHz = 8;
if isfield(cfg,'vpCandidateDiversityHz') && ~isempty(cfg.vpCandidateDiversityHz), diversityHz=cfg.vpCandidateDiversityHz; end
for i = 2:numel(candSSE)
    prevF = candF(1:i-1, :);
    if any(abs(prevF(:, 1) - candF(i, 1)) < diversityHz & abs(prevF(:, 2) - candF(i, 2)) < diversityHz)
        uniqueRows(i) = false;
    end
end
candSSE = candSSE(uniqueRows);
candBeta = candBeta(uniqueRows, :);
candF = candF(uniqueRows, :);

nOut = min(numCandidates, numel(candSSE));
if nOut == 0
    error('No valid VP-main candidates were found.');
end

for i = 1:nOut
    fitOne = pack_varpro_fit(candBeta(i, :).', candF(i, :), t);
    fitOne.linear_sse = candSSE(i);
    if i == 1
        candidates = fitOne;
    else
        candidates(i) = fitOne;
    end
end

fit = candidates(1);
fit.candidates = candidates;
fit.J1 = candSSE(1);
if numel(candSSE) >= 2
    fit.J2 = candSSE(2);
else
    fit.J2 = inf;
end
fit.deltaJ = fit.J2 - fit.J1;
fit.rhoJ = fit.deltaJ / max(fit.J1, eps);
fit.y = y;
fit.F0 = F0;
fit.Fx = Fx;
end
