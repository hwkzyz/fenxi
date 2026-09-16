function [p,c] = R5_Evaluate_CanonicalForward(z,f,B,M)
%R5_EVALUATECANONICALFORWARD Unified low-template differential-gap operator.
% z = [A, phase, dx, dg_1 ... dg_S].  Static registration, template,
% baseline-gap path and gain live in M; only z is dynamic.  This function
% is the shared forward backend for experimental R5 and synthetic checks.
eo=double(B.EO); if ~isfinite(eo), eo=f/B.rotFreqMeanHz; end
phase=eo*B.Theta+z(2)+2*pi*(f-eo*B.rotFreqMeanHz)*(B.T-mean(B.T));
u=z(1)*sin(phase); p=nan(size(B.V));
c=struct('foundation',nan(size(B.V)),'noGapCurrent',nan(size(B.V)), ...
    'noGapBase',nan(size(B.V)),'gapIncrement',nan(size(B.V)), ...
    'xCurrent',B.X-z(3)-u,'xBase',nan(size(B.V)),'uCurrent',u);
for i=1:numel(M)
    q=B.sensorIndex==i;
    [p(q),d]=M(i).evaluate(z(3+i),c.xCurrent(q));
    c.noGapCurrent(q)=d.noGapMv; c.gapIncrement(q)=d.gapIncrementMv;
end
if isfield(B,'useNestedFoundationAnchor') && B.useNestedFoundationAnchor
    assert(all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx])), ...
        'R5:IncompleteFoundationAnchor');
    z0=[B.anchorA B.anchorPhi B.anchorDx zeros(1,numel(M))];
    eo0=B.anchorEO; f0=B.anchorF;
    u0=z0(1)*sin(eo0*B.Theta+z0(2)+2*pi*(f0-eo0*B.rotFreqMeanHz)*(B.T-mean(B.T)));
    p0=nan(size(B.V));
    for i=1:numel(M)
        q=B.sensorIndex==i;
        [p0(q),~]=M(i).evaluate(0,B.X(q)-z0(3)-u0(q));
    end
    q=isfinite(B.anchorV)&isfinite(p)&isfinite(p0);
    p(q)=B.anchorV(q)+p(q)-p0(q);
    c.foundation=B.anchorV; c.noGapBase=p0;
else
    z0=z; z0(4:end)=0; u0=z0(1)*sin(phase-z(2));
    for i=1:numel(M)
        q=B.sensorIndex==i;
        [c.noGapBase(q),~]=M(i).evaluate(0,B.X(q)-z0(3)-u0(q));
    end
end
end
