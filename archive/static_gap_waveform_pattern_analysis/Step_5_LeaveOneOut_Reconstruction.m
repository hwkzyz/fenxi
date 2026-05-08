%% Step 5: validate whether a few static gap samples can reconstruct a held-out gap waveform
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S');

nGap = numel(gapList);
nX = numel(xGrid);
xiGrid = linspace(-3, 3, nX)';

looNrmse = zeros(nGap, 1);
sparseNrmse = zeros(nGap, 1);
S_pred_loo = nan(size(S));
S_pred_sparse = nan(size(S));

for iHold = 1:nGap
    trainIdx = setdiff(1:nGap, iHold);
    [S_pred_loo(iHold, :), looNrmse(iHold)] = reconstruct_holdout( ...
        gapList(trainIdx), S(trainIdx, :), gapList(iHold), S(iHold, :), xGrid, xiGrid);

    sparseTrainIdx = choose_local_training_subset(gapList, iHold, 3);
    [S_pred_sparse(iHold, :), sparseNrmse(iHold)] = reconstruct_holdout( ...
        gapList(sparseTrainIdx), S(sparseTrainIdx, :), gapList(iHold), S(iHold, :), xGrid, xiGrid);
end

isEdge = false(nGap, 1);
isEdge([1, nGap]) = true;
isInterior = ~isEdge;

looSummary = table( ...
    mean(looNrmse), mean(looNrmse(isInterior)), mean(looNrmse(isEdge)), max(looNrmse), ...
    'VariableNames', {'mean_nrmse_all', 'mean_nrmse_interior', 'mean_nrmse_edge', 'max_nrmse'});
sparseSummary = table( ...
    mean(sparseNrmse), mean(sparseNrmse(isInterior)), mean(sparseNrmse(isEdge)), max(sparseNrmse), ...
    'VariableNames', {'mean_nrmse_all', 'mean_nrmse_interior', 'mean_nrmse_edge', 'max_nrmse'});

holdoutTable = table(gapList, looNrmse, sparseNrmse, isEdge, ...
    'VariableNames', {'gap_mm', 'loo_nrmse', 'sparse3_nrmse', 'is_edge'});

save(fullfile(outDir, 'stage5_leaveoneout_reconstruction.mat'), ...
    'holdoutTable', 'looSummary', 'sparseSummary', ...
    'S_pred_loo', 'S_pred_sparse', '-v7.3');

figure('Name', 'Static Gap Step 5 - Holdout Reconstruction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 24, 13]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
bar(categorical(string(gapList)), [looNrmse, sparseNrmse], 'grouped');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Holdout reconstruction error');
legend({'leave-one-out (6 train)', 'sparse local (3 train)'}, 'Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(gapList, looNrmse, 'o-', 'LineWidth', 1.2, 'DisplayName', 'leave-one-out');
plot(gapList, sparseNrmse, 's-', 'LineWidth', 1.2, 'DisplayName', 'sparse local');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Error trend from edge to interior');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iMid = ceil(nGap / 2);
plot(xGrid, S(iMid, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iMid)));
plot(xGrid, S_pred_loo(iMid, :), '--', 'LineWidth', 1.3, 'DisplayName', 'leave-one-out pred');
plot(xGrid, S_pred_sparse(iMid, :), ':', 'LineWidth', 1.6, 'DisplayName', 'sparse local pred');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Interior holdout reconstruction (g = %.1f mm)', gapList(iMid)));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iEdge = 1;
plot(xGrid, S(iEdge, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iEdge)));
plot(xGrid, S_pred_loo(iEdge, :), '--', 'LineWidth', 1.3, 'DisplayName', 'leave-one-out pred');
plot(xGrid, S_pred_sparse(iEdge, :), ':', 'LineWidth', 1.6, 'DisplayName', 'sparse local pred');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Edge holdout reconstruction (g = %.1f mm)', gapList(iEdge)));
legend('Location', 'best', 'Box', 'off');

fprintf('[Static Step 5] Leave-one-out reconstruction summary:\n');
disp(looSummary);
fprintf('[Static Step 5] Sparse local reconstruction summary:\n');
disp(sparseSummary);

function [sPred, nrmse] = reconstruct_holdout(gTrain, sTrain, gQuery, sTrue, xGrid, xiGrid)
    [base, amp, xc, fwhm, hShape] = fit_standard_shape(gTrain, sTrain, xGrid, xiGrid);

    baseHat = predict_linear(gTrain, base, gQuery);
    ampHat = predict_linear(1 ./ gTrain, amp, 1 / gQuery);
    xcHat = mean(xc, 'omitnan');
    fwhmHat = predict_linear(1 ./ gTrain, fwhm, 1 / gQuery);

    if ~isfinite(fwhmHat) || fwhmHat <= 0
        fwhmHat = max(range(xGrid) / 6, eps);
    end

    xiQuery = (xGrid - xcHat) / max(fwhmHat, eps);
    hOnX = interp1(xiGrid, hShape, xiQuery, 'pchip', 'extrap');
    hOnX = hOnX(:).';
    sPred = baseHat + ampHat * hOnX;
    nrmse = sqrt(mean((sTrue(:).' - sPred).^2)) / max(range(sTrue), eps);
end

function [base, amp, xc, fwhm, hShape] = fit_standard_shape(gTrain, sTrain, xGrid, xiGrid)
    nTrain = numel(gTrain);
    nX = numel(xGrid);
    base = min(sTrain, [], 2);
    amp = max(sTrain - base, [], 2);
    xc = zeros(nTrain, 1);
    fwhm = zeros(nTrain, 1);
    sNorm = nan(nTrain, nX);

    for ii = 1:nTrain
        y = sTrain(ii, :).';
        [~, idxPeak] = max(y);
        xc(ii) = xGrid(idxPeak);
        fwhm(ii) = calc_fwhm(xGrid, y);
        if ~isfinite(fwhm(ii)) || fwhm(ii) <= 0
            fwhm(ii) = max(range(xGrid) / 6, eps);
        end
        yNorm = (y - base(ii)) / max(amp(ii), eps);
        xi = (xGrid - xc(ii)) / max(fwhm(ii), eps);
        sNorm(ii, :) = interp1(xi, yNorm, xiGrid, 'pchip', NaN).';
    end

    hShape = mean(fillmissing(sNorm, 'linear', 2, 'EndValues', 'nearest'), 1, 'omitnan').';
    hShape = hShape - min(hShape);
    hShape = hShape / max(max(hShape), eps);
end

function width = calc_fwhm(x, y)
    halfLevel = min(y) + 0.5 * range(y);
    above = y >= halfLevel;
    idx = find(diff(above) ~= 0);
    if numel(idx) < 2
        width = NaN;
        return;
    end
    xCross = zeros(numel(idx), 1);
    for ii = 1:numel(idx)
        jj = idx(ii);
        xCross(ii) = interp1(y(jj:jj+1), x(jj:jj+1), halfLevel, 'linear', 'extrap');
    end
    width = xCross(end) - xCross(1);
end

function yHat = predict_linear(xTrain, yTrain, xQuery)
    xTrain = xTrain(:);
    yTrain = yTrain(:);
    beta = [ones(size(xTrain)), xTrain] \ yTrain;
    yHat = [1, xQuery] * beta;
end

function subsetIdx = choose_local_training_subset(gapList, holdIdx, nKeep)
    gapQuery = gapList(holdIdx);
    candidateIdx = setdiff(1:numel(gapList), holdIdx);
    [~, order] = sort(abs(gapList(candidateIdx) - gapQuery), 'ascend');
    subsetIdx = candidateIdx(order(1:min(nKeep, numel(candidateIdx))));
end
