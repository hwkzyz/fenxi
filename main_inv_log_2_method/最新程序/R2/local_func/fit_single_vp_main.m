function fit = fit_single_vp_main(highMap, templateLib, cfg, staticState)
%FIT_SINGLE_VP_MAIN Independent one-frequency voltage-domain VP scan.
%
% For a fixed gap and frequency, the first-order voltage model is linear in
% [dx, a, b]. No nuisance second frequency is introduced.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
gCenter = staticState.gHat;
if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    gapHalf = get_field(cfg, 'route30GapHalfWidthMm', 0.70);
    gapGrid = linspace(gCenter-gapHalf, gCenter+gapHalf, ...
        get_field(cfg, 'route30GapCount', 61));
else
    gapMin = max(0.05, min(templateLib.gapTrain)-cfg.rawGapSearchMargin);
    gapMax = max(templateLib.gapTrain)+cfg.rawGapSearchMargin;
    gapGrid = linspace(max(gapMin,gCenter-get_field(cfg,'mainVpGapHalfWidth',.45)), ...
        min(gapMax,gCenter+get_field(cfg,'mainVpGapHalfWidth',.45)), ...
        get_field(cfg,'mainVpGapN',61));
end
freqGrid = cfg.route30SingleGrid(:).';
SAll = sin(2*pi*t*freqGrid);
CAll = cos(2*pi*t*freqGrid);
scores = inf(numel(gapGrid),1);
best = repmat(struct('g',NaN,'dx',NaN,'a',NaN,'b',NaN,'f',NaN, ...
    'sse',Inf), numel(gapGrid), 1);
candidateBank=repmat(best(1),0,1);
keepPerGap=get_field(cfg,'singleVpKeepPerGap',5);

for ig = 1:numel(gapGrid)
    g = gapGrid(ig);
    F0 = eval_gap_template(templateLib, g, x);
    Fx = eval_gap_derivative(templateLib, g, x);
    y = V-F0;
    valid = isfinite(Fx) & isfinite(y) & abs(Fx)>eps;
    if nnz(valid) < 8, continue; end
    q = -Fx(valid);
    yv = y(valid); S=SAll(valid,:); C=CAll(valid,:);
    W = abs(Fx(valid)); qW=q.*W; yW=yv.*W;
    d=qW.^2; rhsSignal=qW.*yW;
    nf=numel(freqGrid);G=zeros(3,3,nf);
    G(1,1,:)=sum(d);
    G(1,2,:)=d.'*S;G(2,1,:)=G(1,2,:);
    G(1,3,:)=d.'*C;G(3,1,:)=G(1,3,:);
    G(2,2,:)=sum(d.*S.^2,1);
    G(2,3,:)=sum(d.*S.*C,1);G(3,2,:)=G(2,3,:);
    G(3,3,:)=sum(d.*C.^2,1);
    rhs=[repmat(sum(rhsSignal),1,nf);rhsSignal.'*S;rhsSignal.'*C];
    beta=reshape(pagemldivide(G,reshape(rhs,3,1,nf)),3,nf);
    sse=max(sum(yW.^2)-sum(rhs.*beta,1),0);
    sse(~isfinite(sse)|~all(isfinite(beta),1))=Inf;
    [bestSse,jf]=min(sse);
    bestRow=struct('g',g,'dx',beta(1,jf),'a',beta(2,jf), ...
        'b',beta(3,jf),'f',freqGrid(jf),'sse',bestSse);
    [~,localOrder]=sort(sse,'ascend');
    localOrder=localOrder(1:min(keepPerGap,numel(localOrder)));
    for q=localOrder
        candidateBank(end+1,1)=struct('g',g,'dx',beta(1,q),'a',beta(2,q),...
            'b',beta(3,q),'f',freqGrid(q),'sse',sse(q)); %#ok<AGROW>
    end
    scores(ig) = bestSse;
    best(ig) = bestRow;
end

[J1, ig] = min(scores);
if ~isfinite(J1)
    error('invlog2:NoSingleVPCandidate','No valid single-frequency VP candidate.');
end
row = best(ig);
fit = struct();
fit.g_used = row.g;
fit.dx = row.dx;
fit.dx_used = row.dx;
fit.f_id = row.f;
fit.A_id = hypot(row.a,row.b);
fit.phi_id = atan2(row.b,row.a);
fit.p = [fit.A_id, fit.phi_id, fit.f_id];
fit.J1 = J1;
fit.J2 = second_best(scores, J1);
fit.deltaJ = fit.J2-fit.J1;
fit.rhoJ = fit.deltaJ/max(fit.J1,eps);
if isempty(candidateBank)
    fit.candidates=fit;
else
    [~,order]=sort([candidateBank.sse],'ascend');
    maxCandidates=min(get_field(cfg,'route30SingleCandidateCount',50),numel(order));
    selected=repmat(pack_candidate(candidateBank(order(1))),0,1);
    diversity=get_field(cfg,'route30SingleCandidateDiversityHz',5);
    for ii=order(:).'
        c=pack_candidate(candidateBank(ii));
        if isempty(selected)||all(abs([selected.f_id]-c.f_id)>=diversity | ...
                abs([selected.g_used]-c.g_used)>.02)
            selected(end+1,1)=c; %#ok<AGROW>
            if numel(selected)>=maxCandidates,break;end
        end
    end
    fit.candidates=selected;
end
fit.gap_grid = gapGrid(:);
fit.score = scores;
end

function c=pack_candidate(row)
A=hypot(row.a,row.b);phi=atan2(row.b,row.a);
c=struct('g_used',row.g,'dx_used',row.dx,'f_id',row.f,'A_id',A,...
    'phi_id',phi,'p',[A,phi,row.f],'J1',row.sse,'linear_sse',row.sse);
end

function value = second_best(scores, first)
scores = scores(isfinite(scores) & scores > first*(1+1e-12));
if isempty(scores), value = Inf; else, value = min(scores); end
end

function value = get_field(s,name,defaultValue)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), value=s.(name); else, value=defaultValue; end
end
