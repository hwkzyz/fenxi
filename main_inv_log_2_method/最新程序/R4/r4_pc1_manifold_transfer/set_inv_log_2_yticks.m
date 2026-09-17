function set_inv_log_2_yticks(ax, nTicks)
%set_inv_log_2_yticks  Keep y-axis major ticks sparse for paper figures.

if nargin < 2 || isempty(nTicks)
    nTicks = 5;
end

yl = ylim(ax);
if any(~isfinite(yl)) || yl(1) == yl(2)
    return;
end

if yl(1) >= 0 && yl(2) > 0
    yMax = yl(2);
    rawStep = yMax / max(nTicks - 1, 1);
    step = nice_step_local(rawStep);
    yMaxNice = ceil(yMax / step) * step;
    ticks = 0:step:yMaxNice;
    if numel(ticks) < 4
        step = step / 2;
        ticks = 0:step:yMaxNice;
    elseif numel(ticks) > 6
        step = step * 2;
        yMaxNice = ceil(yMax / step) * step;
        ticks = 0:step:yMaxNice;
    end
    ylim(ax, [0, yMaxNice]);
    yticks(ax, ticks);
    return;
end

rawStep = (yl(2) - yl(1)) / max(nTicks - 1, 1);
step = nice_step_local(rawStep);
yMinNice = floor(yl(1) / step) * step;
yMaxNice = ceil(yl(2) / step) * step;
ticks = yMinNice:step:yMaxNice;
if numel(ticks) < 4
    step = step / 2;
    ticks = yMinNice:step:yMaxNice;
elseif numel(ticks) > 6
    step = step * 2;
    yMinNice = floor(yl(1) / step) * step;
    yMaxNice = ceil(yl(2) / step) * step;
    ticks = yMinNice:step:yMaxNice;
end
ylim(ax, [yMinNice, yMaxNice]);
yticks(ax, ticks);
end

function step = nice_step_local(rawStep)
if ~isfinite(rawStep) || rawStep <= 0
    step = 1;
    return;
end
exponent = floor(log10(rawStep));
base = 10 ^ exponent;
candidates = base * [1, 2, 2.5, 5, 10];
step = candidates(find(candidates >= rawStep, 1, 'first'));
if isempty(step)
    step = candidates(end);
end
end
