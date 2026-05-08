function fit = fit_wave_varpro(t, V, x, templateLib, gHat, dxHat, cfg, precomp)
%fit_wave_varpro  Frequency-grid variable projection under fixed gap.

numCandidates = max(cfg.numVarproCandidates, 10);
maxKeep = max(numCandidates * 4, numCandidates);
candSSE = inf(maxKeep, 1);
candBeta = zeros(maxKeep, 5);
candF = zeros(maxKeep, 2);
if nargin >= 8 && ~isempty(precomp)
    Fx = precomp.Fx;
    DV = precomp.DV;
    S1 = precomp.S1;
    C1 = precomp.C1;
    S2 = precomp.S2;
    C2 = precomp.C2;
else
    F0 = eval_gap_template(templateLib, gHat, x - dxHat);
    Fx = eval_gap_derivative(templateLib, gHat, x - dxHat);
    DV = V - F0;
    S1 = sin(2*pi*t*cfg.f1Grid(:).');
    C1 = cos(2*pi*t*cfg.f1Grid(:).');
    S2 = sin(2*pi*t*cfg.f2Grid(:).');
    C2 = cos(2*pi*t*cfg.f2Grid(:).');
end

w = Fx.^2;
fy = -Fx .* DV;
dvNorm2 = dot(DV, DV);

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
        sse = max(dvNorm2 - dot(rhs, beta), 0);
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
for i = 2:numel(candSSE)
    prevF = candF(1:i-1, :);
    if any(abs(prevF(:, 1) - candF(i, 1)) < 8 & abs(prevF(:, 2) - candF(i, 2)) < 8)
        uniqueRows(i) = false;
    end
end
candSSE = candSSE(uniqueRows);
candBeta = candBeta(uniqueRows, :);
candF = candF(uniqueRows, :);

nOut = min(numCandidates, numel(candSSE));
if nOut == 0
    error('No valid frequency candidates were found.');
end
for i = 1:nOut
    fitOne = pack_linear_fit(candBeta(i, :).', candF(i, :), t);
    fitOne.linear_sse = candSSE(i);
    if i == 1
        candidates = fitOne;
    else
        candidates(i) = fitOne;
    end
end

fit = candidates(1);
fit.candidates = candidates;
end
