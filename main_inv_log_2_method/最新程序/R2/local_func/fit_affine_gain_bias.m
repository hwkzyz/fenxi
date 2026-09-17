function [cHat, bHat, sse] = fit_affine_gain_bias(F, V, cBounds)
%fit_affine_gain_bias  Fit affine gain/bias under optional gain bounds.

if nargin < 3 || isempty(cBounds)
    cBounds = [-inf, inf];
end

F = F(:);
V = V(:);
H = [F, ones(size(F))];
beta = H \ V;
cHat = beta(1);
bHat = beta(2);

if cHat < cBounds(1) || cHat > cBounds(2)
    cHat = min(max(cHat, cBounds(1)), cBounds(2));
    bHat = mean(V - cHat * F);
end

res = V - (cHat * F + bHat);
sse = dot(res, res);
end
