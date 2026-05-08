function trustInfo = select_gap_trust_domain(templateLib, gapList, xCell, yCell, ...
    gHoldout, xGridN, opts)
%select_gap_trust_domain  Pick a fitting x-window from static-library validation.
%
% The rule follows the same spirit as the weighted SG program: use the full
% static library to find a pulse, but fit only in the x-range where the
% static response model is stable under leave-one-gap reconstruction.

domainFull = templateLib.domain;
maxHalf = min(abs(domainFull));
if isempty(maxHalf) || ~isfinite(maxHalf) || maxHalf <= 0
    maxHalf = 0.5 * diff(domainFull);
end
if isfield(opts, 'maxHalfWidth') && ~isempty(opts.maxHalfWidth)
    maxHalf = min(maxHalf, opts.maxHalfWidth);
end
effectiveDomain = estimate_effective_domain(templateLib, opts);
maxHalf = min(maxHalf, min(abs(effectiveDomain)));

minHalf = min(get_opt(opts, 'minHalfWidth', 2.5), maxHalf);
stepHalf = get_opt(opts, 'stepHalfWidth', 0.25);
candidateHalf = minHalf:stepHalf:maxHalf;
if isempty(candidateHalf) || candidateHalf(end) < maxHalf - 1e-9
    candidateHalf = [candidateHalf, maxHalf];
end
candidateHalf = unique(candidateHalf(:).', 'stable');

templateRmse = nan(numel(candidateHalf), 1);
derivativeRmse = nan(numel(candidateHalf), 1);
gapSensitivity = nan(numel(candidateHalf), 1);
score = nan(numel(candidateHalf), 1);
validGapCount = zeros(numel(candidateHalf), 1);

for ih = 1:numel(candidateHalf)
    halfWidth = candidateHalf(ih);
    domainNow = [max(domainFull(1), -halfWidth), min(domainFull(2), halfWidth)];
    [templateRmse(ih), derivativeRmse(ih), validGapCount(ih)] = ...
        evaluate_leave_one_error(templateLib, gapList, xCell, yCell, ...
        gHoldout, xGridN, domainNow, opts);
    gapSensitivity(ih) = evaluate_gap_sensitivity(templateLib, domainNow);
end

valid = isfinite(templateRmse) & isfinite(derivativeRmse) & ...
    isfinite(gapSensitivity) & gapSensitivity > 0 & validGapCount > 0;
if ~any(valid)
    selectedHalf = maxHalf;
else
    tmplNorm = templateRmse ./ max(min(templateRmse(valid)), eps);
    derivNorm = derivativeRmse ./ max(min(derivativeRmse(valid)), eps);
    sensNorm = gapSensitivity ./ max(max(gapSensitivity(valid)), eps);
    narrowPenalty = maxHalf ./ max(candidateHalf(:), eps) - 1;
    score = tmplNorm + get_opt(opts, 'derivativeWeight', 0.20) .* derivNorm + ...
        get_opt(opts, 'sensitivityWeight', 0.30) .* (1 ./ max(sensNorm, eps) - 1) + ...
        get_opt(opts, 'narrowPenaltyWeight', 0.03) .* narrowPenalty;
    bestScore = min(score(valid));
    eligible = valid & score <= bestScore * (1 + get_opt(opts, 'scoreTolerance', 0.08));
    eligibleIdx = find(eligible);
    [~, widestIdx] = max(candidateHalf(eligibleIdx));
    selectedHalf = candidateHalf(eligibleIdx(widestIdx));
end

selectedDomain = [max(domainFull(1), -selectedHalf), min(domainFull(2), selectedHalf)];
trustInfo = struct();
trustInfo.domain = selectedDomain;
trustInfo.selectedHalfWidth = selectedHalf;
trustInfo.fullDomain = domainFull;
trustInfo.effectiveDomain = effectiveDomain;
trustInfo.selectionMode = get_opt(opts, 'selectionMode', "leave_one_gap_static_validation");
trustInfo.gHoldout = gHoldout;
trustInfo.metrics = table(candidateHalf(:), templateRmse, derivativeRmse, score, ...
    gapSensitivity, validGapCount, 'VariableNames', {'half_width_mm','template_rmse', ...
    'derivative_rmse','score','gap_sensitivity','valid_gap_count'});
end

function [templateErr, derivativeErr, validCount] = evaluate_leave_one_error( ...
    refLib, gapList, xCell, yCell, ~, xGridN, domainNow, opts)
templateErrList = [];
derivativeErrList = [];
validCount = 0;
minGridPoints = get_opt(opts, 'minGridPoints', 40);

for ig = 1:numel(gapList)
    gTest = gapList(ig);
    testHoldout = gTest;
    testLib = build_gap_template_library(gapList, xCell, yCell, testHoldout, xGridN);
    testLib = copy_interp_settings(refLib, testLib);

    xEval = testLib.xGrid(:);
    xEval = xEval(xEval >= domainNow(1) & xEval <= domainNow(2));
    if numel(xEval) < minGridPoints
        continue;
    end
    yTruth = interp1(xCell{ig}(:), yCell{ig}(:), xEval, 'pchip');
    yPred = eval_gap_template(testLib, gTest, xEval);
    dyTruth = gradient(yTruth, mean(diff(xEval)));
    dyPred = eval_gap_derivative(testLib, gTest, xEval);

    yScale = max(range(yTruth), eps);
    dyScale = max(range(dyTruth), eps);
    templateErrList(end + 1, 1) = sqrt(mean((yPred(:) - yTruth(:)).^2)) / yScale; %#ok<AGROW>
    derivativeErrList(end + 1, 1) = sqrt(mean((dyPred(:) - dyTruth(:)).^2)) / dyScale; %#ok<AGROW>
    validCount = validCount + 1;
end

templateErr = mean(templateErrList, 'omitnan');
derivativeErr = mean(derivativeErrList, 'omitnan');
end

function effectiveDomain = estimate_effective_domain(templateLib, opts)
x = templateLib.xGrid(:);
S = templateLib.S;
meanResponse = mean(S, 1, 'omitnan')';
responseActivity = meanResponse - min(meanResponse);
responseActivity = responseActivity ./ max(max(responseActivity), eps);
gapSensitivity = std(S, 0, 1, 'omitnan')';
gapSensitivity = gapSensitivity ./ max(max(gapSensitivity), eps);

mask = responseActivity >= get_opt(opts, 'minResponseFraction', 0.12) & ...
    gapSensitivity >= get_opt(opts, 'minSensitivityFraction', 0.12);
if ~any(mask)
    effectiveDomain = templateLib.domain;
    return;
end

effectiveDomain = [min(x(mask)), max(x(mask))];
end

function sens = evaluate_gap_sensitivity(templateLib, domainNow)
mask = templateLib.xGrid(:) >= domainNow(1) & templateLib.xGrid(:) <= domainNow(2);
if nnz(mask) < 4
    sens = NaN;
    return;
end
SNow = templateLib.S(:, mask);
pointSensitivity = std(SNow, 0, 1, 'omitnan');
sens = mean(pointSensitivity, 'omitnan');
end

function dst = copy_interp_settings(src, dst)
fields = {'gapInterpMode','gapInterpParameter','gapInterpBasisOrder','gapInterpBasisName'};
for ii = 1:numel(fields)
    name = fields{ii};
    if isfield(src, name)
        dst.(name) = src.(name);
    end
end
end

function val = get_opt(opts, name, defaultVal)
if isfield(opts, name) && ~isempty(opts.(name))
    val = opts.(name);
else
    val = defaultVal;
end
end
