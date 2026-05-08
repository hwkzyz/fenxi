%% Step 6: compare holdout reconstruction models and see which law is actually usable
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S');
load(fullfile(outDir, 'stage5_leaveoneout_reconstruction.mat'), 'holdoutTable', 'S_pred_sparse');

nGap = numel(gapList);
nX = numel(xGrid);

S_pred_quad_all = nan(nGap, nX);
S_pred_quad_local = nan(nGap, nX);
quadAllNrmse = zeros(nGap, 1);
quadLocalNrmse = zeros(nGap, 1);

for iHold = 1:nGap
    trainIdxAll = setdiff(1:nGap, iHold);
    trainIdxLocal = choose_local_training_subset(gapList, iHold, 4);

    S_pred_quad_all(iHold, :) = predict_pointwise_quadratic( ...
        gapList(trainIdxAll), S(trainIdxAll, :), gapList(iHold));
    S_pred_quad_local(iHold, :) = predict_pointwise_quadratic( ...
        gapList(trainIdxLocal), S(trainIdxLocal, :), gapList(iHold));

    quadAllNrmse(iHold) = calc_nrmse(S(iHold, :), S_pred_quad_all(iHold, :));
    quadLocalNrmse(iHold) = calc_nrmse(S(iHold, :), S_pred_quad_local(iHold, :));
end

compareTable = table( ...
    gapList, holdoutTable.sparse3_nrmse, quadAllNrmse, quadLocalNrmse, ...
    'VariableNames', {'gap_mm', 'standard_shape_sparse3_nrmse', 'pointwise_quad_all_nrmse', 'pointwise_quad_local4_nrmse'});

compareSummary = table( ...
    mean(holdoutTable.sparse3_nrmse), mean(quadAllNrmse), mean(quadLocalNrmse), ...
    mean(holdoutTable.sparse3_nrmse([1 end])), mean(quadAllNrmse([1 end])), mean(quadLocalNrmse([1 end])), ...
    max(holdoutTable.sparse3_nrmse), max(quadAllNrmse), max(quadLocalNrmse), ...
    'VariableNames', {'mean_shape_sparse3', 'mean_quad_all', 'mean_quad_local4', ...
    'edge_mean_shape_sparse3', 'edge_mean_quad_all', 'edge_mean_quad_local4', ...
    'max_shape_sparse3', 'max_quad_all', 'max_quad_local4'});

save(fullfile(outDir, 'stage6_reconstruction_model_compare.mat'), ...
    'compareTable', 'compareSummary', 'S_pred_quad_all', 'S_pred_quad_local', '-v7.3');

figure('Name', 'Static Gap Step 6 - Reconstruction Model Compare', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 24, 13]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
bar(categorical(string(gapList)), [holdoutTable.sparse3_nrmse, quadAllNrmse, quadLocalNrmse], 'grouped');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Holdout error by reconstruction model');
legend({'standard-shape sparse3', 'pointwise quad all', 'pointwise quad local4'}, ...
    'Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(gapList, holdoutTable.sparse3_nrmse, 'o-', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(gapList, quadAllNrmse, 's-', 'LineWidth', 1.2, 'DisplayName', 'pointwise quad all');
plot(gapList, quadLocalNrmse, 'd-', 'LineWidth', 1.2, 'DisplayName', 'pointwise quad local4');
xlabel('Held-out gap g (mm)');
ylabel('NRMSE');
title('Error trend from edge to interior');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iMid = ceil(nGap / 2);
plot(xGrid, S(iMid, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iMid)));
plot(xGrid, S_pred_sparse(iMid, :), '--', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(xGrid, S_pred_quad_all(iMid, :), ':', 'LineWidth', 1.6, 'DisplayName', 'pointwise quad all');
plot(xGrid, S_pred_quad_local(iMid, :), '-.', 'LineWidth', 1.3, 'DisplayName', 'pointwise quad local4');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Interior holdout reconstruction (g = %.1f mm)', gapList(iMid)));
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
iEdge = 1;
plot(xGrid, S(iEdge, :), 'k-', 'LineWidth', 1.4, 'DisplayName', sprintf('true g = %.1f', gapList(iEdge)));
plot(xGrid, S_pred_sparse(iEdge, :), '--', 'LineWidth', 1.2, 'DisplayName', 'standard-shape sparse3');
plot(xGrid, S_pred_quad_all(iEdge, :), ':', 'LineWidth', 1.6, 'DisplayName', 'pointwise quad all');
plot(xGrid, S_pred_quad_local(iEdge, :), '-.', 'LineWidth', 1.3, 'DisplayName', 'pointwise quad local4');
xlabel('Position x (mm)');
ylabel('Capacitance');
title(sprintf('Edge holdout reconstruction (g = %.1f mm)', gapList(iEdge)));
legend('Location', 'best', 'Box', 'off');

fprintf('[Static Step 6] Reconstruction model comparison summary:\n');
disp(compareSummary);

function sPred = predict_pointwise_quadratic(gTrain, sTrain, gQuery)
    gTrain = gTrain(:);
    X = [ones(size(gTrain)), gTrain, gTrain.^2];
    beta = X \ sTrain;
    sPred = [1, gQuery, gQuery^2] * beta;
end

function subsetIdx = choose_local_training_subset(gapList, holdIdx, nKeep)
    gapQuery = gapList(holdIdx);
    candidateIdx = setdiff(1:numel(gapList), holdIdx);
    [~, order] = sort(abs(gapList(candidateIdx) - gapQuery), 'ascend');
    subsetIdx = candidateIdx(order(1:min(nKeep, numel(candidateIdx))));
end

function value = calc_nrmse(yTrue, yPred)
    value = sqrt(mean((yTrue(:) - yPred(:)).^2)) / max(range(yTrue), eps);
end
