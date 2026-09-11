%% Production calibration gate: reference anchor -> frozen blade maps -> one joint update.
% This analysis is isolated from the production package and never overwrites
% the existing Step05/Step05I results.
clear; clc;

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
sourceFile = fullfile(caseDir, 'analysis', 'revised_nominal_20260831', ...
    'outputs', 'Step05_Response_Surface_20250527.mat');
outDir = fullfile(thisDir, 'outputs');
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end

S = load(sourceFile, 'responseSurface');
R = S.responseSurface;
x = R.xGrid(:);
gaps = R.trueGapMm(:);
Y = R.waveforms;
bladeIds = R.bladeIds(:).';
refId = 2;
refIdx = find(bladeIds == refId, 1);
assert(~isempty(refIdx), 'Internal reference blade %d was not found.', refId);
assert(R.referenceBladeId == refId, 'Step05 reference blade does not match the frozen reference.');

cfg = struct();
cfg.referenceBladeId = refId;
cfg.gDomainMm = [min(gaps), max(gaps)];
cfg.anchorGapsMm = [1.10, 1.30, 1.50, 1.70];
cfg.holdoutGapsMm = [1.20, 1.40, 1.60, 1.80];
cfg.deltaBoundsMm = [-0.10, 0.55];
cfg.muBounds = [-0.020, 0.020];
cfg.deltaGridStepMm = 0.01;
cfg.muGridStep = 0.001;
cfg.profileRelativeTolerance = 0.05;
cfg.maxProfileSpanFraction = 0.60;
cfg.maxHoldoutRatio = 1.50;
cfg.maxHoldoutIncreaseMv = 20;
cfg.maxResidualCorrelation = 0.50;
cfg.otherBladeTotalWeight = 0.50;
cfg.maxRefLooRatio = 1.05;
cfg.maxSensitivityRelativeChange = 0.25;

fitX = R.effectiveWindow(:) & all(isfinite(R.coeff), 2);
xFit = x(fitX);
xCenter = 0.5 * (min(xFit) + max(xFit));
anchorMask = ismembertol(gaps, cfg.anchorGapsMm, 1e-9);
holdoutMask = ismembertol(gaps, cfg.holdoutGapsMm, 1e-9);
assert(nnz(anchorMask) == numel(cfg.anchorGapsMm), 'Anchor gaps do not match the direct library.');
assert(nnz(holdoutMask) == numel(cfg.holdoutGapsMm), 'Holdout gaps do not match the direct library.');

fprintf('Reference blade: internal bladeId=%d (maximum-amplitude reference)\n', refId);
fprintf('Direct gap domain: %.2f..%.2f mm\n', cfg.gDomainMm);
fprintf('Frozen anchors: %s mm\n', mat2str(cfg.anchorGapsMm, 3));
fprintf('Frozen holdouts: %s mm\n', mat2str(cfg.holdoutGapsMm, 3));

oldOffset = read_old_offset_initialization_local(caseDir, bladeIds);
nB = numel(bladeIds);
mappingRows = cell(nB, 1);
profileRows = cell(nB, 1);
maps = repmat(struct('bladeId', NaN, 'deltaGapMm', NaN, 'muGapPerXMm', NaN, ...
    'mappingValid', false, 'modelLimited', false), nB, 1);

for ib = 1:nB
    b = bladeIds(ib);
    maps(ib).bladeId = b;
    if b == refId
        maps(ib).deltaGapMm = 0;
        maps(ib).muGapPerXMm = 0;
        maps(ib).mappingValid = true;
        mappingRows{ib} = table(b, 0, 0, 0, 0, 0, 0, 0, 0, true, false, ...
            'VariableNames', mapping_variable_names_local());
        profileRows{ib} = table();
        continue;
    end

    Yb = squeeze(Y(:, ib, :));
    objective = @(p) mapping_mse_local(p, Yb, x, gaps, R, anchorMask, fitX, xCenter, cfg.gDomainMm);
    deltaGrid = cfg.deltaBoundsMm(1):cfg.deltaGridStepMm:cfg.deltaBoundsMm(2);
    muGrid = cfg.muBounds(1):cfg.muGridStep:cfg.muBounds(2);
    J = nan(numel(deltaGrid), numel(muGrid));
    for id = 1:numel(deltaGrid)
        for im = 1:numel(muGrid)
            J(id, im) = objective([deltaGrid(id), muGrid(im)]);
        end
    end
    [~, idxBest] = min(J(:));
    [idBest, imBest] = ind2sub(size(J), idxBest);
    pGrid = [deltaGrid(idBest), muGrid(imBest)];
    starts = unique([pGrid; oldOffset(ib), 0; 0, 0], 'rows', 'stable');
    pBest = pGrid;
    jBest = objective(pBest);
    for is = 1:size(starts, 1)
        p0 = min(max(starts(is, :), [cfg.deltaBoundsMm(1), cfg.muBounds(1)]), ...
            [cfg.deltaBoundsMm(2), cfg.muBounds(2)]);
        if exist('fmincon', 'file') == 2
            opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                'MaxIterations', 500, 'MaxFunctionEvaluations', 3000, ...
                'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-9);
            p = fmincon(objective, p0, [], [], [], [], ...
                [cfg.deltaBoundsMm(1), cfg.muBounds(1)], ...
                [cfg.deltaBoundsMm(2), cfg.muBounds(2)], [], opts);
        else
            z0 = inverse_logistic_local(p0, [cfg.deltaBoundsMm(1), cfg.muBounds(1)], ...
                [cfg.deltaBoundsMm(2), cfg.muBounds(2)]);
            funZ = @(z) objective(logistic_bounds_local(z, ...
                [cfg.deltaBoundsMm(1), cfg.muBounds(1)], ...
                [cfg.deltaBoundsMm(2), cfg.muBounds(2)]));
            z = fminsearch(funZ, z0, optimset('Display', 'off', 'MaxIter', 1000, ...
                'MaxFunEvals', 3000, 'TolX', 1e-8, 'TolFun', 1e-6));
            p = logistic_bounds_local(z, [cfg.deltaBoundsMm(1), cfg.muBounds(1)], ...
                [cfg.deltaBoundsMm(2), cfg.muBounds(2)]);
        end
        j = objective(p);
        if j < jBest, pBest = p; jBest = j; end
    end

    fitRmse = sqrt(jBest);
    holdMse = mapping_mse_local(pBest, Yb, x, gaps, R, holdoutMask, fitX, xCenter, cfg.gDomainMm);
    holdRmse = sqrt(holdMse);
    [residual, xResidual, gResidual] = mapping_residual_local( ...
        pBest, Yb, x, gaps, R, anchorMask | holdoutMask, fitX, xCenter, cfg.gDomainMm);
    corrX = safe_corr_local(residual, xResidual);
    corrG = safe_corr_local(residual, gResidual);

    % Profile width is measured on the coarse profile itself. Using the
    % refined minimum here can leave the discrete grid with no selected point.
    jProfileMin = min(J(:), [], 'omitnan');
    jThreshold = jProfileMin + cfg.profileRelativeTolerance * max(jProfileMin, eps);
    near = J <= jThreshold;
    if any(near(:))
        [nearD, nearM] = find(near);
        deltaSpanFraction = (max(deltaGrid(nearD)) - min(deltaGrid(nearD))) / diff(cfg.deltaBoundsMm);
        muSpanFraction = (max(muGrid(nearM)) - min(muGrid(nearM))) / diff(cfg.muBounds);
    else
        deltaSpanFraction = 1;
        muSpanFraction = 1;
    end
    boundaryTol = 0.01;
    atBoundary = (pBest(1) - cfg.deltaBoundsMm(1)) / diff(cfg.deltaBoundsMm) < boundaryTol || ...
        (cfg.deltaBoundsMm(2) - pBest(1)) / diff(cfg.deltaBoundsMm) < boundaryTol || ...
        (pBest(2) - cfg.muBounds(1)) / diff(cfg.muBounds) < boundaryTol || ...
        (cfg.muBounds(2) - pBest(2)) / diff(cfg.muBounds) < boundaryTol;
    profileValid = ~atBoundary && deltaSpanFraction <= cfg.maxProfileSpanFraction && ...
        muSpanFraction <= cfg.maxProfileSpanFraction;
    holdoutValid = isfinite(holdRmse) && holdRmse <= max( ...
        cfg.maxHoldoutRatio * fitRmse, fitRmse + cfg.maxHoldoutIncreaseMv);
    residualValid = abs(corrX) <= cfg.maxResidualCorrelation && ...
        abs(corrG) <= cfg.maxResidualCorrelation;
    mappingValid = profileValid && holdoutValid && residualValid;

    maps(ib).deltaGapMm = pBest(1);
    maps(ib).muGapPerXMm = pBest(2);
    maps(ib).mappingValid = mappingValid;
    maps(ib).modelLimited = ~mappingValid;
    mappingRows{ib} = table(b, pBest(1), pBest(2), fitRmse, holdRmse, ...
        deltaSpanFraction, muSpanFraction, corrX, corrG, mappingValid, ~mappingValid, ...
        'VariableNames', mapping_variable_names_local());

    [DG, MG] = ndgrid(deltaGrid, muGrid);
    profileRows{ib} = table(repmat(b, numel(J), 1), DG(:), MG(:), sqrt(J(:)), ...
        'VariableNames', {'bladeId','deltaGapMm','muGapPerXMm','rmseMv'});
    fprintf('B%d: delta=%+.5f mm, mu=%+.6f, fit=%.2f, hold=%.2f mV, corr=[%+.2f,%+.2f], valid=%d\n', ...
        b, pBest(1), pBest(2), fitRmse, holdRmse, corrX, corrG, mappingValid);
end

MappingGate = vertcat(mappingRows{:});
Profile = vertcat(profileRows{:});
writetable(MappingGate, fullfile(outDir, 'ProductionMappingGate_20250527.csv'));
writetable(Profile, fullfile(outDir, 'ProductionMappingProfiles_20250527.csv'));

% Freeze mappings and classify complete paths. Only complete effective-window
% paths inside the direct support can update the joint response surface.
supportRows = {};
eligible = false(nB, numel(gaps));
gPath = nan(numel(x), nB, numel(gaps));
for ib = 1:nB
    for ig = 1:numel(gaps)
        gp = gaps(ig) + maps(ib).deltaGapMm + maps(ib).muGapPerXMm * (x - xCenter);
        gPath(:, ib, ig) = gp;
        fullInside = all(gp(fitX) >= cfg.gDomainMm(1) - 1e-10 & ...
            gp(fitX) <= cfg.gDomainMm(2) + 1e-10);
        eligible(ib, ig) = (ib == refIdx) || (maps(ib).mappingValid && fullInside);
        supportRows{end+1,1} = table(bladeIds(ib), gaps(ig), min(gp(fitX)), ...
            max(gp(fitX)), fullInside, eligible(ib, ig), ...
            'VariableNames', {'bladeId','nominalGapMm','mappedGapMinMm', ...
            'mappedGapMaxMm','fullPathInDirectSupport','eligibleForJoint'}); %#ok<AGROW>
    end
end
Support = vertcat(supportRows{:});
writetable(Support, fullfile(outDir, 'ProductionMappedSupport_20250527.csv'));

coeffJoint = fit_joint_coeff_local(Y, gPath, bladeIds, refIdx, maps, eligible, ...
    fitX, R.coeff, cfg.otherBladeTotalWeight, []);

% Reference-blade leave-one-gap-out comparison with frozen mappings.
looRows = cell(numel(gaps), 1);
for ig = 1:numel(gaps)
    coeff0Loo = fit_reference_coeff_local(squeeze(Y(:, refIdx, :)), gaps, fitX, R.coeff, ig);
    coeffJointLoo = fit_joint_coeff_local(Y, gPath, bladeIds, refIdx, maps, eligible, ...
        fitX, R.coeff, cfg.otherBladeTotalWeight, ig);
    y = squeeze(Y(:, refIdx, ig));
    p0 = eval_surface_local(coeff0Loo, R.g0Mm, gaps(ig) * ones(size(x)), fitX);
    pj = eval_surface_local(coeffJointLoo, R.g0Mm, gaps(ig) * ones(size(x)), fitX);
    looRows{ig} = table(gaps(ig), rms_local(y(fitX) - p0(fitX)), ...
        rms_local(y(fitX) - pj(fitX)), ...
        'VariableNames', {'gapMm','referenceOnlyLooRmseMv','jointLooRmseMv'});
end
ReferenceLoo = vertcat(looRows{:});
writetable(ReferenceLoo, fullfile(outDir, 'ProductionReferenceLOO_20250527.csv'));

% Mapped holdout comparison. The reported holdout set was never used to
% estimate delta/mu; this is the required transfer gate.
mappedRows = {};
for ib = 1:nB
    if ib == refIdx || ~maps(ib).mappingValid, continue; end
    for ig = find(holdoutMask(:)).'
        if ~eligible(ib, ig), continue; end
        y = squeeze(Y(:, ib, ig));
        p0 = eval_surface_local(R.coeff, R.g0Mm, gPath(:, ib, ig), fitX);
        % Exclude this nominal state from every blade while evaluating it,
        % so the joint number is a prediction rather than a training error.
        coeffJointCv = fit_joint_coeff_local(Y, gPath, bladeIds, refIdx, maps, eligible, ...
            fitX, R.coeff, cfg.otherBladeTotalWeight, ig);
        pj = eval_surface_local(coeffJointCv, R.g0Mm, gPath(:, ib, ig), fitX);
        mappedRows{end+1,1} = table(bladeIds(ib), gaps(ig), ...
            rms_local(y(fitX) - p0(fitX)), rms_local(y(fitX) - pj(fitX)), ...
            'VariableNames', {'bladeId','nominalGapMm','referenceOnlyRmseMv','jointRmseMv'}); %#ok<AGROW>
    end
end
if isempty(mappedRows)
    MappedHoldout = table();
else
    MappedHoldout = vertcat(mappedRows{:});
end
writetable(MappedHoldout, fullfile(outDir, 'ProductionMappedHoldout_20250527.csv'));

[fgRel, fxRel, sensitivityTable] = compare_sensitivities_local( ...
    R.coeff, coeffJoint, R.g0Mm, gaps, x, fitX);
writetable(sensitivityTable, fullfile(outDir, 'ProductionSensitivityComparison_20250527.csv'));

refLoo0 = mean(ReferenceLoo.referenceOnlyLooRmseMv, 'omitnan');
refLooJoint = mean(ReferenceLoo.jointLooRmseMv, 'omitnan');
if isempty(MappedHoldout)
    mapped0 = NaN; mappedJoint = NaN; mappedPass = false;
else
    mapped0 = mean(MappedHoldout.referenceOnlyRmseMv, 'omitnan');
    mappedJoint = mean(MappedHoldout.jointRmseMv, 'omitnan');
    mappedPass = mappedJoint <= mapped0;
end
refPass = refLooJoint <= cfg.maxRefLooRatio * refLoo0;
sensitivityPass = fgRel <= cfg.maxSensitivityRelativeChange && ...
    fxRel <= cfg.maxSensitivityRelativeChange;
jointAccepted = refPass && mappedPass && sensitivityPass && any(MappingGate.mappingValid & MappingGate.bladeId ~= refId);

if jointAccepted
    productionCoeff = coeffJoint;
    productionModel = "joint";
else
    productionCoeff = R.coeff;
    productionModel = "reference_only";
end

ProductionResponseSurface = R;
ProductionResponseSurface.coeff = productionCoeff;
ProductionResponseSurface.gTrainMm = gaps;
ProductionResponseSurface.directGapDomainMm = cfg.gDomainMm;
ProductionResponseSurface.referenceBladeId = refId;
ProductionResponseSurface.referenceDefinition = 'internal maximum-amplitude blade';
ProductionResponseSurface.mapping = maps;
ProductionResponseSurface.gPathMm = gPath;
ProductionResponseSurface.mappedStateEligible = eligible;
ProductionResponseSurface.productionModel = char(productionModel);
ProductionResponseSurface.aBlade = 1;
ProductionResponseSurface.additionalBladeOffset = 0;
ProductionResponseSurface.sourceReferenceSurface = sourceFile;
ProductionResponseSurface.note = ['One-way production calibration gate. Mappings are frozen after F0; ' ...
    'the joint surface is updated once and never feeds back into delta/mu.'];

GateSummary = table(refLoo0, refLooJoint, mapped0, mappedJoint, fgRel, fxRel, ...
    refPass, mappedPass, sensitivityPass, jointAccepted, productionModel, ...
    'VariableNames', {'referenceOnlyLooRmseMv','jointLooRmseMv', ...
    'mappedReferenceOnlyRmseMv','mappedJointRmseMv','relativeFgChange', ...
    'relativeFxChange','referenceLooPass','mappedHoldoutPass', ...
    'sensitivityPass','jointAccepted','productionModel'});
writetable(GateSummary, fullfile(outDir, 'ProductionCalibrationGateSummary_20250527.csv'));
save(fullfile(outDir, 'ProductionResponseSurface_20250527.mat'), ...
    'ProductionResponseSurface', 'MappingGate', 'Support', 'ReferenceLoo', ...
    'MappedHoldout', 'GateSummary', 'cfg', '-v7.3');

fprintf('\nProduction calibration gate:\n');
disp(GateSummary);
fprintf('Selected production model: %s\n', productionModel);
fprintf('Outputs: %s\n', outDir);

function names = mapping_variable_names_local()
names = {'bladeId','deltaGapMm','muGapPerXMm','fitRmseMv','holdoutRmseMv', ...
    'profileDeltaSpanFraction','profileMuSpanFraction','residualCorrX', ...
    'residualCorrGap','mappingValid','modelLimited'};
end

function old = read_old_offset_initialization_local(caseDir, bladeIds)
f = fullfile(fileparts(caseDir), 'reference_blade_gap_analysis', 'results', ...
    'analysis_04_blade_offset_summary_matched.csv');
if ~isfile(f)
    f = fullfile(fileparts(fileparts(caseDir)), 'reference_blade_gap_analysis', ...
        'results', 'analysis_04_blade_offset_summary_matched.csv');
end
old = zeros(numel(bladeIds), 1);
if ~isfile(f), return; end
T = readtable(f);
for ib = 1:numel(bladeIds)
    row = T.bladeId == bladeIds(ib);
    if any(row), old(ib) = T.meanOffsetZeroedMm(find(row, 1)); end
end
end

function mse = mapping_mse_local(p, Yb, x, gaps, R, gapMask, fitX, xCenter, domain)
[residual, ~, ~] = mapping_residual_local(p, Yb, x, gaps, R, gapMask, fitX, xCenter, domain);
if isempty(residual), mse = 1e12; else, mse = mean(residual.^2, 'omitnan'); end
end

function [residual, xResidual, gResidual] = mapping_residual_local(p, Yb, x, gaps, R, gapMask, fitX, xCenter, domain)
residual = []; xResidual = []; gResidual = [];
for ig = find(gapMask(:)).'
    gp = gaps(ig) + p(1) + p(2) * (x - xCenter);
    if any(gp(fitX) < domain(1) | gp(fitX) > domain(2))
        residual = []; xResidual = []; gResidual = []; return;
    end
    pred = eval_surface_local(R.coeff, R.g0Mm, gp, fitX);
    y = Yb(:, ig);
    ok = fitX & isfinite(y) & isfinite(pred);
    residual = [residual; y(ok) - pred(ok)]; %#ok<AGROW>
    xResidual = [xResidual; x(ok)]; %#ok<AGROW>
    gResidual = [gResidual; gp(ok)]; %#ok<AGROW>
end
end

function coeff = fit_reference_coeff_local(Yref, gaps, fitX, coeffFallback, excludeGap)
coeff = coeffFallback;
use = true(numel(gaps), 1); use(excludeGap) = false;
X = [ones(nnz(use),1), 1 ./ gaps(use), log(gaps(use))];
for ix = find(fitX(:)).'
    y = Yref(ix, use).';
    ok = isfinite(y) & all(isfinite(X), 2);
    if nnz(ok) >= 3, coeff(ix,:) = (X(ok,:) \ y(ok)).'; end
end
end

function coeff = fit_joint_coeff_local(Y, gPath, bladeIds, refIdx, maps, eligible, fitX, coeffFallback, eta, excludeNominalGap)
coeff = coeffFallback;
nB = numel(bladeIds); nG = size(Y,3);
for ix = find(fitX(:)).'
    design = []; value = []; weight = [];
    validOther = find([maps.mappingValid]);
    validOther(validOther == refIdx) = [];
    validOther = validOther(any(eligible(validOther,:),2));
    for ib = 1:nB
        useG = find(eligible(ib,:));
        if ~isempty(excludeNominalGap), useG(useG == excludeNominalGap) = []; end
        if isempty(useG), continue; end
        g = squeeze(gPath(ix,ib,useG)); g = g(:);
        y = squeeze(Y(ix,ib,useG)); y = y(:);
        X = [ones(numel(g),1), 1./g, log(g)];
        ok = isfinite(y) & all(isfinite(X),2);
        X = X(ok,:); y = y(ok);
        if isempty(y), continue; end
        if ib == refIdx
            w = ones(numel(y),1) / numel(y);
        elseif ismember(ib, validOther)
            w = eta * ones(numel(y),1) / max(numel(validOther),1) / numel(y);
        else
            continue;
        end
        design = [design; X]; value = [value; y]; weight = [weight; w]; %#ok<AGROW>
    end
    if size(design,1) >= 3
        sw = sqrt(weight);
        coeff(ix,:) = ((design .* sw) \ (value .* sw)).';
    end
end
end

function pred = eval_surface_local(coeff, g0, g, fitX)
pred = nan(size(g));
ok = fitX & isfinite(g) & g > 0 & all(isfinite(coeff),2);
pred(ok) = coeff(ok,1) + coeff(ok,2)./g(ok) + coeff(ok,3).*log(g(ok)./g0);
end

function [fgRel, fxRel, T] = compare_sensitivities_local(c0, cj, g0, gaps, x, fitX)
fg0 = []; fgj = []; fx0 = []; fxj = [];
ixFit = find(fitX);
xUse = x(ixFit);
for ig = 1:numel(gaps)
    g = gaps(ig);
    f0 = c0(:,1) + c0(:,2)./g + c0(:,3).*log(g./g0);
    fj = cj(:,1) + cj(:,2)./g + cj(:,3).*log(g./g0);
    dg0 = -c0(:,2)./(g.^2) + c0(:,3)./g;
    dgj = -cj(:,2)./(g.^2) + cj(:,3)./g;
    dfx0 = gradient(f0(ixFit), xUse); dfxj = gradient(fj(ixFit), xUse);
    fg0 = [fg0; dg0(fitX)]; fgj = [fgj; dgj(fitX)]; %#ok<AGROW>
    fx0 = [fx0; dfx0]; fxj = [fxj; dfxj]; %#ok<AGROW>
end
fgRel = norm(fgj-fg0) / max(norm(fg0), eps);
fxRel = norm(fxj-fx0) / max(norm(fx0), eps);
T = table(fgRel, fxRel, rms_local(fg0), rms_local(fgj), rms_local(fx0), rms_local(fxj), ...
    'VariableNames', {'relativeFgChange','relativeFxChange','referenceFgRms', ...
    'jointFgRms','referenceFxRms','jointFxRms'});
end

function c = safe_corr_local(a, b)
ok = isfinite(a) & isfinite(b);
if nnz(ok) < 3 || std(a(ok)) == 0 || std(b(ok)) == 0, c = 0; else, c = corr(a(ok), b(ok)); end
end

function value = rms_local(v)
value = sqrt(mean(v(:).^2, 'omitnan'));
end

function p = logistic_bounds_local(z, lb, ub)
p = lb + (ub-lb) ./ (1 + exp(-z));
end

function z = inverse_logistic_local(p, lb, ub)
q = min(max((p-lb)./(ub-lb), 1e-8), 1-1e-8);
z = log(q./(1-q));
end
