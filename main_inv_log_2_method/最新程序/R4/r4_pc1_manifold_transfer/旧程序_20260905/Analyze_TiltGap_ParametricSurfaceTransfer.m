%% Parametric transfer of tilt-gap waveform surfaces
% Independent analysis only. Raw source data and existing programs are untouched.

thisDir = fileparts(mfilename('fullpath'));
dataPath = fullfile(thisDir, 'output', 'waveform_families', ...
    'TiltGap_WaveformSurfaceData.mat');
outDir = fullfile(thisDir, 'output', 'parametric_surface_transfer');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

raw = load(dataPath, 'Data');
Data = raw.Data;
x = Data.x_mm(:);
gap = Data.gap_mm(:);
tilt = Data.tilt_deg(:);
V = Data.response_V;
[nX, nGap, nTilt] = size(V);
assert(isequal(size(V), [numel(x), numel(gap), numel(tilt)]), ...
    'Data dimensions do not match x, gap, and tilt coordinates.');
assert(all(diff(tilt) > 0), 'Tilt angles must be strictly increasing.');
assert(exist('sgolayfilt', 'file') == 2, ...
    'Signal Processing Toolbox function sgolayfilt is required for SG peak metrics.');

% Each tilt surface is one vector. The model uses a low-rank global basis and
% cubic coefficient curves versus tilt; both basis and curves are refit per holdout.
nSurface = nX * nGap;
Y = reshape(V, nSurface, nTilt);
modeCount = 2;
methodNames = ["Parametric rank-2", "Nearest surface", "Linear angle interpolation"];
nMethod = numel(methodNames);
Yhat = nan(nSurface, nTilt, nMethod);

for iHold = 1:nTilt
    train = setdiff(1:nTilt, iHold);
    Ytrain = Y(:, train);
    mu = mean(Ytrain, 2);
    [U, ~, ~] = svd(Ytrain - mu, 'econ');
    basis = U(:, 1:modeCount);
    score = basis' * (Ytrain - mu);
    alphaTrain = tilt(train);
    alphaTarget = tilt(iHold);
    scoreTarget = zeros(modeCount, 1);
    for iMode = 1:modeCount
        % Cubic is the lowest non-linear continuous model that can express
        % systematic curvature over the sampled tilt range.
        degree = min(3, numel(train) - 1);
        p = polyfit(alphaTrain, score(iMode, :), degree);
        scoreTarget(iMode) = polyval(p, alphaTarget);
    end
    Yhat(:, iHold, 1) = mu + basis * scoreTarget;

    [~, nearestLocal] = min(abs(alphaTrain - alphaTarget));
    Yhat(:, iHold, 2) = Ytrain(:, nearestLocal);

    left = find(alphaTrain < alphaTarget, 1, 'last');
    right = find(alphaTrain > alphaTarget, 1, 'first');
    if isempty(left) || isempty(right)
        % Endpoints are extrapolation cases: use the two closest retained angles.
        [~, order] = sort(abs(alphaTrain - alphaTarget), 'ascend');
        i1 = order(1);
        i2 = order(2);
    else
        i1 = left;
        i2 = right;
    end
    w = (alphaTarget - alphaTrain(i1)) / (alphaTrain(i2) - alphaTrain(i1));
    Yhat(:, iHold, 3) = (1 - w) * Ytrain(:, i1) + w * Ytrain(:, i2);
end

Yhat3 = reshape(Yhat, nX, nGap, nTilt, nMethod);
resultRows = nTilt * nMethod;
angleDeg = zeros(resultRows, 1);
method = strings(resultRows, 1);
rmseV = zeros(resultRows, 1);
nrmseRangePct = zeros(resultRows, 1);
nrmseStdPct = zeros(resultRows, 1);
corrR = zeros(resultRows, 1);
peakXMaeMm = zeros(resultRows, 1);
peakAmpMaeV = zeros(resultRows, 1);
peakXBiasMm = zeros(resultRows, 1);
peakAmpBiasV = zeros(resultRows, 1);

row = 0;
for iTilt = 1:nTilt
    truth = V(:, :, iTilt);
    truePeakX = nan(nGap, 1);
    truePeakAmp = nan(nGap, 1);
    for iGap = 1:nGap
        [truePeakAmp(iGap), idx] = sg_peak(truth(:, iGap));
        truePeakX(iGap) = x(idx);
    end
    for iMethod = 1:nMethod
        row = row + 1;
        pred = Yhat3(:, :, iTilt, iMethod);
        errorSurface = pred - truth;
        predPeakX = nan(nGap, 1);
        predPeakAmp = nan(nGap, 1);
        for iGap = 1:nGap
            [predPeakAmp(iGap), idx] = sg_peak(pred(:, iGap));
            predPeakX(iGap) = x(idx);
        end
        angleDeg(row) = tilt(iTilt);
        method(row) = methodNames(iMethod);
        rmseV(row) = sqrt(mean(errorSurface(:).^2));
        nrmseRangePct(row) = 100 * rmseV(row) / range(truth(:));
        nrmseStdPct(row) = 100 * rmseV(row) / std(truth(:));
        corrR(row) = corr(truth(:), pred(:));
        peakXMaeMm(row) = mean(abs(predPeakX - truePeakX));
        peakAmpMaeV(row) = mean(abs(predPeakAmp - truePeakAmp));
        peakXBiasMm(row) = mean(predPeakX - truePeakX);
        peakAmpBiasV(row) = mean(predPeakAmp - truePeakAmp);
    end
end

perAngle = table(angleDeg, method, rmseV, nrmseRangePct, nrmseStdPct, corrR, ...
    peakXMaeMm, peakAmpMaeV, peakXBiasMm, peakAmpBiasV);
writetable(perAngle, fullfile(outDir, 'LOAO_per_angle_metrics.csv'));

summary = table(methodNames', zeros(nMethod, 1), zeros(nMethod, 1), ...
    zeros(nMethod, 1), zeros(nMethod, 1), zeros(nMethod, 1), zeros(nMethod, 1), ...
    'VariableNames', {'method', 'mean_RMSE_V', 'mean_NRMSE_range_pct', ...
    'mean_correlation_r', 'mean_SG_peak_x_MAE_mm', 'mean_SG_peak_amp_MAE_V', ...
    'worst_NRMSE_range_pct'});
for iMethod = 1:nMethod
    rows = perAngle.method == methodNames(iMethod);
    summary.mean_RMSE_V(iMethod) = mean(perAngle.rmseV(rows));
    summary.mean_NRMSE_range_pct(iMethod) = mean(perAngle.nrmseRangePct(rows));
    summary.mean_correlation_r(iMethod) = mean(perAngle.corrR(rows));
    summary.mean_SG_peak_x_MAE_mm(iMethod) = mean(perAngle.peakXMaeMm(rows));
    summary.mean_SG_peak_amp_MAE_V(iMethod) = mean(perAngle.peakAmpMaeV(rows));
    summary.worst_NRMSE_range_pct(iMethod) = max(perAngle.nrmseRangePct(rows));
end
writetable(summary, fullfile(outDir, 'LOAO_summary_metrics.csv'));

% Rank evidence comes from the complete observed surface family only; it is
% descriptive rather than a held-out score.
Ycenter = Y - mean(Y, 2);
[~, Sfull, ~] = svd(Ycenter, 'econ');
energy = diag(Sfull).^2;
rankTable = table((1:numel(energy))', energy / sum(energy), cumsum(energy) / sum(energy), ...
    'VariableNames', {'mode', 'explained_variance_fraction', 'cumulative_fraction'});
writetable(rankTable, fullfile(outDir, 'surface_family_rank_evidence.csv'));

% Figure 1: LOAO comparison across every held-out tilt.
fig1 = figure('Name', 'LOAO surface-transfer comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 10]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
metrics = {'nrmseRangePct', 'corrR', 'peakXMaeMm', 'peakAmpMaeV'};
yLabels = {'NRMSE / range (%)', 'Surface correlation r', 'SG peak x MAE (mm)', 'SG peak amplitude MAE (V)'};
for iMetric = 1:numel(metrics)
    nexttile;
    hold on;
    for iMethod = 1:nMethod
        rows = perAngle.method == methodNames(iMethod);
        plot(perAngle.angleDeg(rows), perAngle.(metrics{iMetric})(rows), '-o', ...
            'LineWidth', 1.1, 'MarkerSize', 4, 'DisplayName', methodNames(iMethod));
    end
    xlabel('Held-out tilt (deg)', 'Interpreter', 'tex');
    ylabel(yLabels{iMetric}, 'Interpreter', 'tex');
    box on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 8, 'TickDir', 'in');
    if iMetric == 1
        legend('Location', 'best', 'FontSize', 7, 'Box', 'off');
    end
end
exportgraphics(fig1, fullfile(outDir, 'Fig_LOAO_metric_comparison.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'Fig_LOAO_metric_comparison.pdf'), 'ContentType', 'vector');

% Figure 2: strict interior interpolation example at 2 deg.
iExample = find(tilt == 2, 1);
if isempty(iExample), iExample = ceil(nTilt / 2); end
fig2 = figure('Name', 'Interior-angle surface prediction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [3 3 17 8]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
surfaces = cat(3, V(:, :, iExample), Yhat3(:, :, iExample, 1), ...
    Yhat3(:, :, iExample, 2), Yhat3(:, :, iExample, 3));
titles = {sprintf('Observed: %.1f deg', tilt(iExample)), 'Rank-2 parametric', ...
    'Nearest surface', 'Linear interpolation'};
cl = [min(surfaces(:)), max(surfaces(:))];
for k = 1:4
    nexttile;
    imagesc(gap, x, surfaces(:, :, k)); axis xy tight;
    clim(cl); title(titles{k}, 'FontWeight', 'normal');
    xlabel('Gap (mm)');
    if k == 1, ylabel('x (mm)'); end
    box on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 8, 'TickDir', 'in');
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'V (V)';
exportgraphics(fig2, fullfile(outDir, 'Fig_interior_angle_surface_prediction.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'Fig_interior_angle_surface_prediction.pdf'), 'ContentType', 'vector');

save(fullfile(outDir, 'TiltGap_ParametricSurfaceTransferResults.mat'), ...
    'perAngle', 'summary', 'rankTable', 'Yhat3', 'methodNames', 'modeCount', ...
    'x', 'gap', 'tilt', 'dataPath');

fprintf('\nLOAO mean metrics:\n');
disp(summary);
fprintf('Rank-1 explained variance: %.5f; rank-2 cumulative: %.5f\n', ...
    rankTable.explained_variance_fraction(1), rankTable.cumulative_fraction(2));
fprintf('Outputs written to:\n%s\n', outDir);

function [amplitude, index] = sg_peak(curve)
% Peak location and amplitude after a fixed, traceable SG smoothing rule.
smoothed = sgolayfilt(curve(:), 3, 11);
[amplitude, index] = max(smoothed);
end
