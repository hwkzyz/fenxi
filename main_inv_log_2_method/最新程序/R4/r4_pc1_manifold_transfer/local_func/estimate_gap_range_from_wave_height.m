function gapRange = estimate_gap_range_from_wave_height(highMap, pathModel, cfg)
%ESTIMATE_GAP_RANGE_FROM_WAVE_HEIGHT  Coarse gap interval from high-speed height.
% Height is evaluated for each sensor-turn group.  A bounded horizontal shift
% is treated as a nuisance variable so that vibration-induced sampling shifts
% are not directly interpreted as clearance changes.

if nargin < 3 || isempty(cfg), cfg = struct(); end
xAll = highMap.x_v(:);
vAll = highMap.V_a(:);
revAll = highMap.rev_v(:);
sensorAll = highMap.S_v(:);
valid = isfinite(xAll) & isfinite(vAll) & isfinite(revAll) & isfinite(sensorAll);
xAll = xAll(valid);
vAll = vAll(valid);
revAll = revAll(valid);
sensorAll = sensorAll(valid);

[~, ~, groupId] = unique([revAll, sensorAll], 'rows', 'stable');
nGroup = max(groupId);
xGroup = cell(nGroup, 1);
vGroup = cell(nGroup, 1);
obsHeight = nan(nGroup, 1);
for k = 1:nGroup
    idx = groupId == k;
    xGroup{k} = xAll(idx);
    vGroup{k} = vAll(idx);
    if nnz(idx) >= 10
        obsHeight(k) = robust_wave_height(vGroup{k});
    end
end
keepGroup = isfinite(obsHeight);
xGroup = xGroup(keepGroup);
vGroup = vGroup(keepGroup); %#ok<NASGU>
obsHeight = obsHeight(keepGroup);
if numel(obsHeight) < 2
    error('Too few valid sensor-turn groups for waveform-height localization.');
end

heightSigma = robust_scale(obsHeight);
heightSigma = max(heightSigma, ...
    get_cfg_field(cfg, 'heightNoiseFloorFraction', 0.02) * median(obsHeight));
heightSigma = max(heightSigma, eps);

xTemplate = pathModel.xTemplate(:);
lowHeight = robust_wave_height(pathModel.templateLow(:));
g0 = pathModel.pathCal.g0;
dMin = get_cfg_field(cfg, 'heightGapMin', min(pathModel.templateLib.gapTrain) - g0);
dMax = get_cfg_field(cfg, 'heightGapMax', max(pathModel.templateLib.gapTrain) - g0);
deltaGrid = linspace(dMin, dMax, get_cfg_field(cfg, 'heightGapGridN', 81));
shiftHalfWidth = get_cfg_field(cfg, 'heightShiftHalfWidth', 0.50);
shiftGrid = linspace(-shiftHalfWidth, shiftHalfWidth, ...
    get_cfg_field(cfg, 'heightShiftGridN', 21));

groupCost = nan(numel(deltaGrid), numel(obsHeight));
bestShift = nan(size(groupCost));
predMedianHeight = nan(size(deltaGrid));
for ig = 1:numel(deltaGrid)
    predTemplate = eval_path_increment_template(pathModel, deltaGrid(ig), xTemplate);
    predHeights = nan(numel(obsHeight), numel(shiftGrid));
    for k = 1:numel(obsHeight)
        xk = xGroup{k};
        for is = 1:numel(shiftGrid)
            vk = interp1(xTemplate, predTemplate, xk - shiftGrid(is), ...
                'pchip', NaN);
            predHeights(k, is) = robust_wave_height(vk);
        end
        mismatch = abs(predHeights(k, :) - obsHeight(k)) ./ heightSigma;
        [groupCost(ig, k), iShift] = min(mismatch);
        bestShift(ig, k) = shiftGrid(iShift);
    end
    predMedianHeight(ig) = median(predHeights(:), 'omitnan');
end

% A median across sensor-turn groups prevents one noisy turn from dominating
% the coarse interval.  This cost is used only for localization.
cost = median(groupCost, 2, 'omitnan');
[minCost, iBest] = min(cost);
costTol = get_cfg_field(cfg, 'heightCostTolerance', 2.5);
inside = cost <= minCost + costTol;
if any(inside)
    dLo = min(deltaGrid(inside));
    dHi = max(deltaGrid(inside));
else
    dLo = deltaGrid(iBest);
    dHi = deltaGrid(iBest);
end
dCenter = deltaGrid(iBest);
gridStep = mean(diff(deltaGrid));
minHalf = get_cfg_field(cfg, 'heightMinHalfWidth', 2 * gridStep);
halfWidth = max([dCenter - dLo, dHi - dCenter, minHalf]);

gapRange = struct();
gapRange.delta_center = dCenter;
gapRange.delta_half_width = halfWidth;
gapRange.delta_lower = max(dMin, dCenter - halfWidth);
gapRange.delta_upper = min(dMax, dCenter + halfWidth);
gapRange.g_center = g0 + dCenter;
gapRange.g_lower = g0 + gapRange.delta_lower;
gapRange.g_upper = g0 + gapRange.delta_upper;
gapRange.observed_height_median = median(obsHeight);
gapRange.observed_height_std = heightSigma;
gapRange.low_height = lowHeight;
gapRange.delta_grid = deltaGrid(:);
gapRange.predicted_height_median = predMedianHeight(:);
gapRange.cost = cost(:);
gapRange.group_cost = groupCost;
gapRange.best_shift = bestShift;
gapRange.shift_grid = shiftGrid(:);
gapRange.num_groups = numel(obsHeight);
gapRange.cost_tolerance = costTol;
gapRange.method = "sensor_turn_height_with_bounded_shift";
end

function h = robust_wave_height(v)
v = v(isfinite(v));
if numel(v) < 5
    h = NaN;
    return;
end
q = prctile(v, [5, 95]);
h = q(2) - q(1);
end

function s = robust_scale(v)
v = v(isfinite(v));
med = median(v);
s = 1.4826 * median(abs(v - med));
end

function value = get_cfg_field(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
