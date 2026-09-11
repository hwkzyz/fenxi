%% Step06B: Build the production response surface from frozen multi-blade paths
% This is the final static calibration object. Blade ID enters only through
% the frozen path g_bj(x), never as a response-surface dimension.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20251222();
outDir = packageCfg.paths.calibrationWork;
step05File = fullfile(outDir, 'Step05_Response_Surface_20251222_RefBladeAnchor.mat');
pathFile = fullfile(outDir, 'Step06_BladePathModel_20251222_RefBladeAnchor.mat');
if ~isfile(step05File) || ~isfile(pathFile)
    error('Run Main05_RefBladeAnchor_20251222 and Main06_EstimateBladePaths_20251222 first.');
end

S = load(step05File, 'responseSurface');
P = load(pathFile, 'BladePathModel');
base = S.responseSurface;
paths = P.BladePathModel;
xGrid = base.xGrid(:);
waveforms = base.waveforms;
gPath = paths.gapPathMm;
pathValid = paths.pathValid;
bladeIds = base.bladeIds(:).';
trueGaps = base.trueGapMm(:).';
fitMask = paths.fitMask(:);
g0Mm = base.g0Mm;
minSamples = 12;

coeff = nan(numel(xGrid), 3);
Yfit = nan(size(waveforms));
sampleCount = zeros(numel(xGrid), 1);
conditionNumber = nan(numel(xGrid), 1);
for ix = 1:numel(xGrid)
    if ~fitMask(ix)
        continue;
    end
    g = squeeze(gPath(ix, :, :));
    y = squeeze(waveforms(ix, :, :));
    valid = squeeze(pathValid(ix, :, :)) & isfinite(g) & isfinite(y);
    g = g(valid);
    y = y(valid);
    if numel(g) < minSamples
        continue;
    end
    A = [ones(numel(g), 1), 1 ./ g, log(g ./ g0Mm)];
    if rank(A) < 3
        continue;
    end
    coeff(ix, :) = (A \ y).';
    conditionNumber(ix) = cond(A);
    sampleCount(ix) = numel(g);
    for ib = 1:numel(bladeIds)
        for ig = 1:numel(trueGaps)
            if pathValid(ix, ib, ig)
                gg = gPath(ix, ib, ig);
                Yfit(ix, ib, ig) = coeff(ix,1) + coeff(ix,2) / gg + ...
                    coeff(ix,3) * log(gg / g0Mm);
            end
        end
    end
end

fitRows = cell(numel(bladeIds) * numel(trueGaps), 1);
k = 0;
for ib = 1:numel(bladeIds)
    for ig = 1:numel(trueGaps)
        k = k + 1;
        valid = fitMask & pathValid(:,ib,ig) & isfinite(waveforms(:,ib,ig)) & isfinite(Yfit(:,ib,ig));
        residual = waveforms(valid,ib,ig) - Yfit(valid,ib,ig);
        gNow = gPath(valid,ib,ig);
        if isempty(residual)
            rmseMv = NaN;
            relRmse = NaN;
            gMinNow = NaN;
            gMaxNow = NaN;
        else
            rmseMv = sqrt(mean(residual.^2, 'omitnan'));
            relRmse = rmseMv / max(abs(waveforms(valid,ib,ig)), [], 'omitnan');
            gMinNow = min(gNow);
            gMaxNow = max(gNow);
        end
        fitRows{k} = table(bladeIds(ib), trueGaps(ig), nnz(valid), gMinNow, gMaxNow, rmseMv, relRmse, ...
            'VariableNames', {'bladeId','nominalGapMm','pointCount','pathGapMinMm','pathGapMaxMm','rmseMv','relativeRmse'});
    end
end
fitTable = vertcat(fitRows{:});
bladeSummary = groupsummary(fitTable, 'bladeId', {'mean','std'}, {'rmseMv','relativeRmse'});

dFdgGrid = nan(numel(xGrid), numel(trueGaps));
dFdxGrid = nan(numel(xGrid), numel(trueGaps));
for ig = 1:numel(trueGaps)
    g = trueGaps(ig);
    dFdgGrid(:,ig) = -coeff(:,2) ./ g.^2 + coeff(:,3) ./ g;
    yHorizontal = coeff(:,1) + coeff(:,2) ./ g + coeff(:,3) .* log(g ./ g0Mm);
    dFdxGrid(:,ig) = gradient(yHorizontal, xGrid);
end

JointResponseSurface = base;
JointResponseSurface.method = 'B2_anchor_frozen_paths_joint_multi_blade_response_surface';
JointResponseSurface.description = ['B2 defines the gap coordinate. Frozen blade paths map all ' ...
    'static waveforms onto one shared F(g,x); blade ID is not a model dimension.'];
JointResponseSurface.sourceStep05File = step05File;
JointResponseSurface.sourceBladePathFile = pathFile;
JointResponseSurface.pathTable = paths.pathTable;
JointResponseSurface.gapPathMm = gPath;
JointResponseSurface.pathValid = pathValid;
JointResponseSurface.coeff = coeff;
JointResponseSurface.YfitAllBlades = Yfit;
JointResponseSurface.fitTable = fitTable;
JointResponseSurface.bladeFitSummary = bladeSummary;
JointResponseSurface.sampleCountByX = sampleCount;
JointResponseSurface.conditionNumberByX = conditionNumber;
JointResponseSurface.dFdgGrid = dFdgGrid;
JointResponseSurface.dFdxGrid = dFdxGrid;
JointResponseSurface.waveformAggregationLabel = 'all valid frozen multi-blade paths';

matFile = fullfile(outDir, 'Step06B_Joint_Response_Surface_20251222_RefBladeAnchor.mat');
fitCsv = fullfile(outDir, 'Step06B_Joint_Response_Surface_Fit_20251222_RefBladeAnchor.csv');
summaryCsv = fullfile(outDir, 'Step06B_Joint_Response_Surface_BladeSummary_20251222_RefBladeAnchor.csv');
save(matFile, 'JointResponseSurface', '-v7.3');
writetable(fitTable, fitCsv);
writetable(bladeSummary, summaryCsv);

fig = figure('Color', 'w', 'Position', [80, 80, 1300, 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; imagesc(trueGaps, xGrid, coeff(:,1) + coeff(:,2) ./ trueGaps(round(end/2)) + ...
    coeff(:,3) .* log(trueGaps(round(end/2)) ./ g0Mm)); axis xy; colorbar;
xlabel('Nominal gap (mm)'); ylabel('x (mm)'); title('Representative horizontal reconstruction');
nexttile; hold on; grid on; box on;
for ib = 1:numel(bladeIds)
    row = fitTable.bladeId == bladeIds(ib);
    plot(fitTable.nominalGapMm(row), fitTable.rmseMv(row), 'o-', 'DisplayName', sprintf('B%d',bladeIds(ib)));
end
xlabel('Nominal gap (mm)'); ylabel('Joint-fit RMSE (mV)'); legend('Location','best');
nexttile; plot(xGrid, sampleCount, 'LineWidth', 1.2); grid on; box on;
xlabel('x (mm)'); ylabel('Valid samples'); title('Samples constraining each coefficient triplet');
nexttile; semilogy(xGrid, conditionNumber, 'LineWidth', 1.2); grid on; box on;
xlabel('x (mm)'); ylabel('cond(A)'); title('Local regression conditioning');
saveas(fig, fullfile(outDir, 'Step06B_Joint_Response_Surface_Diagnostics_20251222_RefBladeAnchor.png'));

fprintf('Step06B joint response surface saved to:\n  %s\n', matFile);
