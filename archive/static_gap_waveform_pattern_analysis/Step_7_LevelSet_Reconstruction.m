%% Step 7: reconstruct holdout waveforms through equal-level left/right branch geometry
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S', 'cfg');
load(fullfile(outDir, 'stage5_leaveoneout_reconstruction.mat'), 'holdoutTable', 'S_pred_sparse');

nGap = numel(gapList);
qLevels = cfg.levelList(:);
levelsetNrmse = zeros(nGap, 1);
S_pred_levelset = nan(size(S));

for iHold = 1:nGap
    trainIdx = choose_local_training_subset(gapList, iHold, 4);
    [S_pred_levelset(iHold, :), levelsetNrmse(iHold)] = reconstruct_by_levelset( ...
        gapList(trainIdx), S(trainIdx, :), gapList(iHold), S(iHold, :), xGrid, qLevels);
end

isEdge = false(nGap, 1);
isEdge([1, nGap]) = true;
isInterior = ~isEdge;

levelsetSummary = table( ...
    mean(levelsetNrmse), mean(levelsetNrmse(isInterior)), mean(levelsetNrmse(isEdge)), max(levelsetNrmse), ...
    'VariableNames', {'mean_nrmse_all', 'mean_nrmse_interior', 'mean_nrmse_edge', 'max_nrmse'});

levelsetTable = table(gapList, holdoutTable.sparse3_nrmse, levelsetNrmse, ...
    'VariableNames', {'gap_mm', 'standard_shape_sparse3_nrmse', 'levelset_local4_nrmse'});

save(fullfile(outDir, 'stage7_levelset_reconstruction.mat'), ...
    'levelsetTable', 'levelsetSummary', 'S_pred_levelset', '-v7.3');

figure('Name', 'Static Gap Step 7 - Level-Set Reconstruction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 24, 13]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
bar(categorical(string(gapList)), [holdoutTable.sparse3_nrmse, levelsetNrmse], 'grouped');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Holdout error: standard-shape vs level-set geometry');
legend({'standard-shape sparse3', 'level-set local4'}, 'Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(gapList, holdoutTable.sparse3_nrmse, 'o-', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(gapList, levelsetNrmse, 'd-', 'LineWidth', 1.2, 'DisplayName', 'level-set local4');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Error trend from edge to interior');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iMid = ceil(nGap / 2);
plot(xGrid, S(iMid, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iMid)));
plot(xGrid, S_pred_sparse(iMid, :), '--', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(xGrid, S_pred_levelset(iMid, :), '-.', 'LineWidth', 1.5, 'DisplayName', 'level-set local4');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Interior holdout reconstruction (g = %.1f mm)', gapList(iMid)));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iEdge = 1;
plot(xGrid, S(iEdge, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iEdge)));
plot(xGrid, S_pred_sparse(iEdge, :), '--', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(xGrid, S_pred_levelset(iEdge, :), '-.', 'LineWidth', 1.5, 'DisplayName', 'level-set local4');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Edge holdout reconstruction (g = %.1f mm)', gapList(iEdge)));
legend('Location', 'best', 'Box', 'off');

fprintf('[Static Step 7] Level-set reconstruction summary:\n');
disp(levelsetSummary);

function [sPred, nrmse] = reconstruct_by_levelset(gTrain, sTrain, gQuery, sTrue, xGrid, qLevels)
    nTrain = numel(gTrain);
    nQ = numel(qLevels);
    base = min(sTrain, [], 2);
    peak = max(sTrain, [], 2);
    amp = peak - base;
    xPeak = zeros(nTrain, 1);
    xLeft = nan(nTrain, nQ);
    xRight = nan(nTrain, nQ);

    for ii = 1:nTrain
        y = sTrain(ii, :).';
        yNorm = (y - base(ii)) / max(amp(ii), eps);
        [~, idxPeak] = max(yNorm);
        xPeak(ii) = xGrid(idxPeak);
        [xLeft(ii, :), xRight(ii, :)] = extract_level_positions(xGrid, yNorm, qLevels, idxPeak);
    end

    baseHat = predict_linear(gTrain, base, gQuery);
    ampHat = predict_linear(1 ./ gTrain, amp, 1 / gQuery);
    xPeakHat = mean(xPeak, 'omitnan');

    xLeftHat = zeros(1, nQ);
    xRightHat = zeros(1, nQ);
    for iq = 1:nQ
        xLeftHat(iq) = predict_linear(1 ./ gTrain, xLeft(:, iq), 1 / gQuery);
        xRightHat(iq) = predict_linear(1 ./ gTrain, xRight(:, iq), 1 / gQuery);
    end

    xLeftHat = max(xLeftHat, xGrid(1));
    xRightHat = min(xRightHat, xGrid(end));
    xLeftHat = cummax(xLeftHat);
    xRightHat = -cummax(-xRightHat);
    xLeftHat(end) = xPeakHat;
    xRightHat(end) = xPeakHat;

    qFull = [0; qLevels; 1];
    yFull = baseHat + ampHat * qFull;
    xLeftFull = [xGrid(1); xLeftHat(:); xPeakHat];
    xRightFullDesc = [xGrid(end); xRightHat(:); xPeakHat];
    xRightFullAsc = flipud(xRightFullDesc);
    yRightFullAsc = flipud(yFull);

    xPts = [xLeftFull; xRightFullAsc(2:end)];
    yPts = [yFull; yRightFullAsc(2:end)];

    [xPtsUnique, ia] = unique(xPts, 'stable');
    yPtsUnique = yPts(ia);
    sPred = interp1(xPtsUnique, yPtsUnique, xGrid, 'pchip', 'extrap');
    sPred = sPred(:).';
    nrmse = sqrt(mean((sTrue(:).' - sPred).^2)) / max(range(sTrue), eps);
end

function [xLeft, xRight] = extract_level_positions(xGrid, yNorm, qLevels, idxPeak)
    xLeft = nan(1, numel(qLevels));
    xRight = nan(1, numel(qLevels));

    xL = xGrid(1:idxPeak);
    yL = yNorm(1:idxPeak);
    [yLUniq, iaL] = unique(yL, 'stable');
    xLUniq = xL(iaL);

    xR = xGrid(idxPeak:end);
    yR = yNorm(idxPeak:end);
    [yRUniq, iaR] = unique(flipud(yR), 'stable');
    xRUniq = flipud(xR);
    xRUniq = xRUniq(iaR);

    for iq = 1:numel(qLevels)
        q = qLevels(iq);
        xLeft(iq) = interp1(yLUniq, xLUniq, q, 'linear', 'extrap');
        xRight(iq) = interp1(yRUniq, xRUniq, q, 'linear', 'extrap');
    end
end

function subsetIdx = choose_local_training_subset(gapList, holdIdx, nKeep)
    gapQuery = gapList(holdIdx);
    candidateIdx = setdiff(1:numel(gapList), holdIdx);
    [~, order] = sort(abs(gapList(candidateIdx) - gapQuery), 'ascend');
    subsetIdx = candidateIdx(order(1:min(nKeep, numel(candidateIdx))));
end

function yHat = predict_linear(xTrain, yTrain, xQuery)
    xTrain = xTrain(:);
    yTrain = yTrain(:);
    beta = [ones(size(xTrain)), xTrain] \ yTrain;
    yHat = [1, xQuery] * beta;
end
