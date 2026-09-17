function fit = fit_vp_main(t, V, x, templateLib, gHat, cfg, precomp, gapIdx)
%fit_vp_main  Variable projection over frequencies under fixed gap.
if isfield(cfg,'vpUseStagedGrid') && cfg.vpUseStagedGrid
    fit = fit_vp_main_staged(t,V,x,templateLib,gHat,cfg,precomp,gapIdx);
    return;
end

numCandidates = max(cfg.numVarproCandidates, 10);
maxKeep = max(numCandidates * 4, numCandidates);
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

[candSSE,candBeta,candF,searchPairCount] = solve_frequency_pairs_batched(cfg,maxKeep,yNorm2, ...
    g11,g1s1,g1c1,g1s2,g1c2,gs1s1,gs1c1,gc1c1, ...
    gs2s2,gs2c2,gc2c2,gs1s2,gs1c2,gc1s2,gc1c2, ...
    r1,rs1,rc1,rs2,rc2);

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
fit.search_pair_count = searchPairCount;
end

function [bestSSE,bestBeta,bestF,pairCount] = solve_frequency_pairs_batched(cfg,maxKeep,yNorm2, ...
    g11,g1s1,g1c1,g1s2,g1c2,gs1s1,gs1c1,gc1c1, ...
    gs2s2,gs2c2,gc2c2,gs1s2,gs1c2,gc1s2,gc1c2, ...
    r1,rs1,rc1,rs2,rc2)
% Solve the same 5-by-5 normal equations as the scalar loop, but in pages.
% Every trusted waveform sample has already contributed to the sufficient
% statistics above; batching only removes MATLAB loop/solver overhead.
f1=cfg.f1Grid(:); f2=cfg.f2Grid(:);
requested=[];
if isfield(cfg,'vpCandidateFrequencyPairsHz')&&~isempty(cfg.vpCandidateFrequencyPairsHz)
    requested=sort(cfg.vpCandidateFrequencyPairsHz,2);
end
if ~isempty(requested)
    [ok1,i1]=ismember(requested(:,1),f1);
    [ok2,i2]=ismember(requested(:,2),f2);
    ok=ok1&ok2&abs(f1(i1)-f2(i2))>=20;
    i1=i1(ok);i2=i2(ok);
else
    [i1,i2]=ndgrid(1:numel(f1),1:numel(f2));
    valid=abs(f1(i1)-f2(i2))>=20;
    sameGrid=numel(f1)==numel(f2)&&isequal(f1,f2);
    if isfield(cfg,'vpUniqueUnorderedPairs')&&cfg.vpUniqueUnorderedPairs&&sameGrid
        valid=valid&(f2(i2)>=f1(i1));
    end
    i1=i1(valid);i2=i2(valid);
end
pairCount=numel(i1);
bestSSE=inf(maxKeep,1);bestBeta=zeros(maxKeep,5);bestF=zeros(maxKeep,2);
chunkSize=4096;
if isfield(cfg,'vpPairBatchSize')&&~isempty(cfg.vpPairBatchSize)
    chunkSize=max(128,round(cfg.vpPairBatchSize));
end
n2=numel(f2);
for first=1:chunkSize:numel(i1)
    q=first:min(first+chunkSize-1,numel(i1));a=i1(q);b=i2(q);m=numel(q);
    lin=sub2ind([numel(f1),n2],a,b);
    G=zeros(5,5,m);
    G(1,1,:)=g11;
    G(1,2,:)=g1s1(a); G(2,1,:)=g1s1(a);
    G(1,3,:)=g1c1(a); G(3,1,:)=g1c1(a);
    G(1,4,:)=g1s2(b); G(4,1,:)=g1s2(b);
    G(1,5,:)=g1c2(b); G(5,1,:)=g1c2(b);
    G(2,2,:)=gs1s1(a); G(2,3,:)=gs1c1(a); G(3,2,:)=gs1c1(a);
    G(3,3,:)=gc1c1(a);
    G(4,4,:)=gs2s2(b); G(4,5,:)=gs2c2(b); G(5,4,:)=gs2c2(b);
    G(5,5,:)=gc2c2(b);
    G(2,4,:)=gs1s2(lin); G(4,2,:)=gs1s2(lin);
    G(2,5,:)=gs1c2(lin); G(5,2,:)=gs1c2(lin);
    G(3,4,:)=gc1s2(lin); G(4,3,:)=gc1s2(lin);
    G(3,5,:)=gc1c2(lin); G(5,3,:)=gc1c2(lin);
    rhs=[repmat(r1,1,m);reshape(rs1(a),1,m);reshape(rc1(a),1,m); ...
        reshape(rs2(b),1,m);reshape(rc2(b),1,m)];
    beta=reshape(pagemldivide(G,reshape(rhs,5,1,m)),5,m);
    sse=max(yNorm2-sum(rhs.*beta,1),0).';
    finite=isfinite(sse)&all(isfinite(beta),1).';
    if ~any(finite),continue;end
    local=[sse(finite),beta(:,finite).',f1(a(finite)),f2(b(finite))];
    combined=[bestSSE,bestBeta,bestF;local];
    [~,ord]=sort(combined(:,1),'ascend');
    combined=combined(ord(1:min(maxKeep,size(combined,1))),:);
    if size(combined,1)<maxKeep
        combined(end+1:maxKeep,:)=repmat([Inf,zeros(1,7)],maxKeep-size(combined,1),1);
    end
    bestSSE=combined(:,1);bestBeta=combined(:,2:6);bestF=combined(:,7:8);
end
end
