function p = local_percentile(x, pct)
%local_percentile  Small percentile helper that avoids toolbox assumptions.

x = sort(x(isfinite(x)));
if isempty(x)
    p = NaN;
    return;
end

q = 1 + (numel(x) - 1) * pct / 100;
lo = floor(q);
hi = ceil(q);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (q - lo) * (x(hi) - x(lo));
end
end
