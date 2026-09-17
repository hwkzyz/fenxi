function gapRange = estimate_gap_range_from_wave_features(highMap, pathModel, cfg)
%ESTIMATE_GAP_RANGE_FROM_WAVE_FEATURES  Coarse gap range from three waveform features.
% Features: peak height, positive area, and full width at half maximum.
% The result is a search interval only; final identification uses the full
% waveform objective.

if nargin < 3 || isempty(cfg), cfg = struct(); end
xAll = highMap.x_v(:);
vAll = highMap.V_a(:);
revAll = highMap.rev_v(:);
sensorAll = highMap.S_v(:);
valid = isfinite(xAll) & isfinite(vAll) & isfinite(revAll) & isfinite(sensorAll);
xAll = xAll(valid); vAll = vAll(valid);
revAll = revAll(valid); sensorAll = sensorAll(valid);
[~, ~, groupId] = unique([revAll, sensorAll], 'rows', 'stable');
nGroup = max(groupId);
xGroup = cell(nGroup, 1);
obsFeature = nan(nGroup, 3);
for k = 1:nGroup
    idx = groupId == k;
    xGroup{k} = xAll(idx);
    obsFeature(k, :) = waveform_features(xAll(idx), vAll(idx));
end
keep = all(isfinite(obsFeature), 2);
xGroup = xGroup(keep);
obsFeature = obsFeature(keep, :);
if size(obsFeature, 1) < 2
    error('Too few valid sensor-turn groups for multi-feature localization.');
end

featureScale = robust_feature_scale(obsFeature);
featureFloor = [get_cfg_field(cfg, 'featureHeightFloor', 0.02), ...
    get_cfg_field(cfg, 'featureAreaFloor', 0.02), ...
    get_cfg_field(cfg, 'featureWidthFloor', 0.02)];
featureScale = max(featureScale, featureFloor);

xTemplate = pathModel.xTemplate(:);
g0 = pathModel.pathCal.g0;
dMin = get_cfg_field(cfg, 'featureGapMin', min(pathModel.templateLib.gapTrain) - g0);
dMax = get_cfg_field(cfg, 'featureGapMax', max(pathModel.templateLib.gapTrain) - g0);
deltaGrid = linspace(dMin, dMax, get_cfg_field(cfg, 'featureGapGridN', 41));
shiftHalfWidth = get_cfg_field(cfg, 'featureShiftHalfWidth', 0.50);
shiftGrid = linspace(-shiftHalfWidth, shiftHalfWidth, ...
    get_cfg_field(cfg, 'featureShiftGridN', 11));

groupCost = nan(numel(deltaGrid), size(obsFeature, 1));
bestShift = nan(size(groupCost));
for ig = 1:numel(deltaGrid)
    predTemplate = eval_path_increment_template(pathModel, deltaGrid(ig), xTemplate);
    for k = 1:size(obsFeature, 1)
        costShift = nan(size(shiftGrid));
        for is = 1:numel(shiftGrid)
            vk = interp1(xTemplate, predTemplate, xGroup{k} - shiftGrid(is), ...
                'pchip', NaN);
            predFeature = waveform_features(xGroup{k}, vk);
            costShift(is) = norm((predFeature - obsFeature(k, :)) ./ featureScale);
        end
        [groupCost(ig, k), iShift] = min(costShift);
        bestShift(ig, k) = shiftGrid(iShift);
    end
end

cost = median(groupCost, 2, 'omitnan');
[minCost, iBest] = min(cost);
tol = get_cfg_field(cfg, 'featureCostTolerance', 2.5);
inside = cost <= minCost + tol;
if any(inside)
    dLo = min(deltaGrid(inside));
    dHi = max(deltaGrid(inside));
else
    dLo = deltaGrid(iBest);
    dHi = deltaGrid(iBest);
end
dCenter = deltaGrid(iBest);
minHalf = get_cfg_field(cfg, 'featureMinHalfWidth', 0.10);
halfWidth = max([dCenter - dLo, dHi - dCenter, minHalf]);

gapRange = struct();
gapRange.delta_center = dCenter;
gapRange.delta_half_width = halfWidth;
gapRange.delta_lower = max(dMin, dCenter - halfWidth);
gapRange.delta_upper = min(dMax, dCenter + halfWidth);
gapRange.g_center = g0 + dCenter;
gapRange.g_lower = g0 + gapRange.delta_lower;
gapRange.g_upper = g0 + gapRange.delta_upper;
gapRange.observed_features = obsFeature;
gapRange.feature_scale = featureScale;
gapRange.delta_grid = deltaGrid(:);
gapRange.cost = cost(:);
gapRange.group_cost = groupCost;
gapRange.best_shift = bestShift;
gapRange.method = "sensor_turn_peak_area_width";
end

function f = waveform_features(x, v)
x = x(:); v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid); v = v(valid);
if numel(x) < 10
    f = [NaN, NaN, NaN];
    return;
end
[x, order] = sort(x);
v = v(order);
[x, keep] = unique(x, 'stable');
v = v(keep);
if numel(x) < 10
    f = [NaN, NaN, NaN];
    return;
end
baseline = prctile(v, 5);
vp = max(v - baseline, 0);
[peak, iPeak] = max(vp);
if peak <= eps
    f = [0, 0, NaN];
    return;
end
area = trapz(x, vp);
half = 0.5 * peak;
left = find(vp(1:iPeak) <= half, 1, 'last');
rightRel = find(vp(iPeak:end) <= half, 1, 'first');
if isempty(left) || isempty(rightRel)
    width = x(end) - x(1);
else
    iRight = iPeak + rightRel - 1;
    xLeft = interpolate_crossing(x, vp, left, left + 1, half);
    xRight = interpolate_crossing(x, vp, iRight - 1, iRight, half);
    width = max(xRight - xLeft, eps);
end
f = [peak, area, width];
end

function xc = interpolate_crossing(x, y, i1, i2, level)
if i1 < 1 || i2 > numel(x) || y(i2) == y(i1)
    xc = x(max(1, min(numel(x), i1)));
else
    xc = x(i1) + (level - y(i1)) * (x(i2) - x(i1)) / (y(i2) - y(i1));
end
end

function s = robust_feature_scale(F)
s = 1.4826 * median(abs(F - median(F, 1)), 1, 'omitnan');
end

function value = get_cfg_field(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
