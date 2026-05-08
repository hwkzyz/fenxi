function idx = pick_evenly_spaced_indices(nTotal, targetCount)
%pick_evenly_spaced_indices  Evenly sample indices with endpoints included.

if nargin < 2 || isempty(targetCount) || targetCount <= 0 || targetCount >= nTotal
    idx = (1:nTotal).';
    return;
end

idx = round(linspace(1, nTotal, targetCount));
idx = unique(idx(:), 'stable');
if idx(1) ~= 1
    idx = [1; idx];
end
if idx(end) ~= nTotal
    idx = [idx; nTotal];
end
end
