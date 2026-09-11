%% Step06: Estimate frozen static gap paths from the B2-anchored initial surface
% B2 defines the physical gap coordinate: delta_B2 = q_B2 = 0.  Each other
% blade has one shared (delta_b, q_b) across all nominal-gap measurements.
% This program does not modify the Step05 reference-anchor output.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20251222();
outDir = packageCfg.paths.calibrationWork;
step05File = fullfile(outDir, 'Step05_Response_Surface_20251222_RefBladeAnchor.mat');
if ~isfile(step05File)
    error('Missing B2-anchored Step05 output: %s', step05File);
end

cfg = struct();
cfg.referenceBladeId = 2;
cfg.deltaBoundsMm = [-0.25, 0.25];
cfg.qBoundsMmPerMm = [-0.06, 0.06];
cfg.minValidGapStates = 8;
cfg.minValidPoints = 300;
cfg.maxIter = 500;
cfg.maxFunEvals = 1000;
cfg.saveFigures = true;

S = load(step05File, 'responseSurface');
base = S.responseSurface;
xGrid = base.xGrid(:);
waveforms = base.waveforms;
bladeIds = base.bladeIds(:).';
trueGaps = base.trueGapMm(:).';
referenceIndex = find(bladeIds == cfg.referenceBladeId, 1);
if isempty(referenceIndex)
    error('Reference blade B%d is absent.', cfg.referenceBladeId);
end

gMin = min(base.gTrainMm);
gMax = max(base.gTrainMm);
fitMask = base.effectiveWindow(:) & all(isfinite(base.coeff), 2);
xCenter = 0;

deltaBlade = nan(size(bladeIds));
qBlade = nan(size(bladeIds));
deltaBlade(referenceIndex) = 0;
qBlade(referenceIndex) = 0;
objectiveRmseMv = nan(size(bladeIds));
validStateCount = zeros(size(bladeIds));

for ib = 1:numel(bladeIds)
    if ib == referenceIndex
        [objectiveRmseMv(ib), validStateCount(ib)] = path_objective_local( ...
            [0, 0], waveforms(:, ib, :), xGrid, trueGaps, base, fitMask, ...
            xCenter, gMin, gMax, cfg);
        continue;
    end

    obj = @(p) path_objective_local(p, waveforms(:, ib, :), xGrid, trueGaps, ...
        base, fitMask, xCenter, gMin, gMax, cfg);
    options = optimset('Display', 'off', 'MaxIter', cfg.maxIter, ...
        'MaxFunEvals', cfg.maxFunEvals, 'TolX', 1e-6, 'TolFun', 1e-4);
    p0 = [0, 0];
    pOpt = fminsearch(@(p) bounded_objective_local(p, obj, cfg), p0, options);
    pOpt(1) = min(max(pOpt(1), cfg.deltaBoundsMm(1)), cfg.deltaBoundsMm(2));
    pOpt(2) = min(max(pOpt(2), cfg.qBoundsMmPerMm(1)), cfg.qBoundsMmPerMm(2));
    [objectiveRmseMv(ib), validStateCount(ib)] = obj(pOpt);
    deltaBlade(ib) = pOpt(1);
    qBlade(ib) = pOpt(2);
end

gapPathMm = nan(numel(xGrid), numel(bladeIds), numel(trueGaps));
pathValid = false(size(gapPathMm));
for ib = 1:numel(bladeIds)
    for ig = 1:numel(trueGaps)
        gPath = trueGaps(ig) + deltaBlade(ib) + qBlade(ib) .* (xGrid - xCenter);
        gapPathMm(:, ib, ig) = gPath;
        pathValid(:, ib, ig) = fitMask & gPath >= gMin & gPath <= gMax;
    end
end

pathTable = table(bladeIds(:), deltaBlade(:), qBlade(:), ...
    atan(qBlade(:)) * 180 / pi, objectiveRmseMv(:), validStateCount(:), ...
    'VariableNames', {'bladeId','deltaGapMm','qGapPerXMm','tiltAngleDeg', ...
    'initialSurfaceRmseMv','validGapStateCount'});

BladePathModel = struct();
BladePathModel.method = 'B2_anchor_then_frozen_multi_path_geometry';
BladePathModel.description = ['B2 defines delta=0 and q=0. Each other blade is fitted ' ...
    'once across all nominal-gap waveforms against the B2 initial surface.'];
BladePathModel.sourceStep05File = step05File;
BladePathModel.cfg = cfg;
BladePathModel.gMinMm = gMin;
BladePathModel.gMaxMm = gMax;
BladePathModel.xCenterMm = xCenter;
BladePathModel.xGrid = xGrid;
BladePathModel.fitMask = fitMask;
BladePathModel.bladeIds = bladeIds;
BladePathModel.trueGapMm = trueGaps;
BladePathModel.pathTable = pathTable;
BladePathModel.gapPathMm = gapPathMm;
BladePathModel.pathValid = pathValid;

matFile = fullfile(outDir, 'Step06_BladePathModel_20251222_RefBladeAnchor.mat');
csvFile = fullfile(outDir, 'Step06_BladePathSummary_20251222_RefBladeAnchor.csv');
save(matFile, 'BladePathModel', '-v7.3');
writetable(pathTable, csvFile);

if cfg.saveFigures
    fig = figure('Color', 'w', 'Position', [80, 80, 1300, 700]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on; grid on; box on;
    colors = lines(numel(bladeIds));
    midGapIdx = max(1, round(numel(trueGaps)/2));
    for ib = 1:numel(bladeIds)
        plot(xGrid, gapPathMm(:, ib, midGapIdx), 'Color', colors(ib,:), ...
            'LineWidth', 1.3, 'DisplayName', sprintf('B%d', bladeIds(ib)));
    end
    yline(gMin, 'k--', 'HandleVisibility', 'off');
    yline(gMax, 'k--', 'HandleVisibility', 'off');
    xlabel('B2-anchored x (mm)'); ylabel('Path gap (mm)');
    title(sprintf('Frozen paths at nominal gap %.2f mm', trueGaps(midGapIdx)));
    legend('Location', 'best');
    nexttile; hold on; grid on; box on;
    yyaxis left; bar(bladeIds, deltaBlade, 0.55); ylabel('\delta_b (mm)');
    yyaxis right; plot(bladeIds, qBlade, 'ko-', 'LineWidth', 1.2); ylabel('q_b (mm/mm)');
    xlabel('Blade ID'); title('Estimated static path parameters');
    saveas(fig, fullfile(outDir, 'Step06_BladePaths_20251222_RefBladeAnchor.png'));
end

disp(pathTable);
fprintf('Step06 blade paths saved to:\n  %s\n', matFile);

function value = bounded_objective_local(p, obj, cfg)
    if p(1) < cfg.deltaBoundsMm(1) || p(1) > cfg.deltaBoundsMm(2) || ...
            p(2) < cfg.qBoundsMmPerMm(1) || p(2) > cfg.qBoundsMmPerMm(2)
        value = 1e9 + 1e8 * sum(max([cfg.deltaBoundsMm(1)-p(1), p(1)-cfg.deltaBoundsMm(2), ...
            cfg.qBoundsMmPerMm(1)-p(2), p(2)-cfg.qBoundsMmPerMm(2), 0]).^2);
        return;
    end
    value = obj(p);
end

function [rmseMv, validStateCount] = path_objective_local(p, Yraw, xGrid, trueGaps, base, fitMask, xCenter, gMin, gMax, cfg)
    Y = squeeze(Yraw);
    residual = [];
    validStateCount = 0;
    for ig = 1:numel(trueGaps)
        g = trueGaps(ig) + p(1) + p(2) .* (xGrid - xCenter);
        valid = fitMask & g >= gMin & g <= gMax & isfinite(Y(:, ig));
        if nnz(valid) < cfg.minValidPoints / numel(trueGaps)
            continue;
        end
        prediction = eval_surface_local(base.coeff, base.g0Mm, g);
        residual = [residual; Y(valid, ig) - prediction(valid)]; %#ok<AGROW>
        validStateCount = validStateCount + 1;
    end
    if validStateCount < cfg.minValidGapStates || numel(residual) < cfg.minValidPoints
        rmseMv = 1e8 + 1e6 * (cfg.minValidGapStates - validStateCount);
        return;
    end
    rmseMv = sqrt(mean(residual.^2));
end

function value = eval_surface_local(coeff, g0Mm, g)
    value = nan(size(g));
    valid = isfinite(g) & all(isfinite(coeff), 2);
    value(valid) = coeff(valid,1) + coeff(valid,2) ./ g(valid) + ...
        coeff(valid,3) .* log(g(valid) ./ g0Mm);
end
