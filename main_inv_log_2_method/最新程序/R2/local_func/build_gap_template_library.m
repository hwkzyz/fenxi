function templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN)
%build_gap_template_library  Build an interpolated static-gap template set.

if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
trainIds = find(trainIdx);
G_train = gapList(trainIds);

xMin = -inf;
xMax = inf;
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end

xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(G_train), numel(xGrid));
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    S(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end

dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ii = 1:size(S, 1)
    FxMat(ii, :) = gradient(S(ii, :), dx);
end

templateLib = struct('gapTrain', G_train(:), 'xGrid', xGrid(:), ...
    'S', S, 'FxMat', FxMat, 'domain', [min(xGrid), max(xGrid)]);
end
