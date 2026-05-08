%% Step 4: fit a shared standard-shape model with gap-dependent parameters
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S');

nGap = numel(gapList);
nX = numel(xGrid);

base = min(S, [], 2);
amp = max(S - base, [], 2);
xc = zeros(nGap, 1);
fwhm = zeros(nGap, 1);
xiGrid = linspace(-3, 3, nX)';
S_norm = nan(nGap, nX);

for i = 1:nGap
    y = S(i, :).';
    [~, idxPeak] = max(y);
    xc(i) = xGrid(idxPeak);
    fwhm(i) = calc_fwhm(xGrid, y);
    if ~isfinite(fwhm(i)) || fwhm(i) <= 0
        fwhm(i) = max(range(xGrid) / 6, eps);
    end
    yNorm = (y - base(i)) / max(amp(i), eps);
    xi = (xGrid - xc(i)) / max(fwhm(i), eps);
    S_norm(i, :) = interp1(xi, yNorm, xiGrid, 'pchip', NaN).';
end

H = mean(fillmissing(S_norm, 'linear', 2, 'EndValues', 'nearest'), 1, 'omitnan').';
H = H - min(H);
H = H / max(max(H), eps);

S_model = nan(size(S));
fitRmse = zeros(nGap, 1);
fitNrmse = zeros(nGap, 1);
for i = 1:nGap
    xiNow = (xGrid - xc(i)) / max(fwhm(i), eps);
    H_on_x = interp1(xiGrid, H, xiNow, 'pchip', 'extrap');
    H_on_x = H_on_x(:).';
    S_model(i, :) = base(i) + amp(i) * H_on_x;
    fitRmse(i) = sqrt(mean((S(i, :) - S_model(i, :)).^2));
    fitNrmse(i) = fitRmse(i) / max(range(S(i, :)), eps);
end

invGap = 1 ./ gapList;
r2Base_g = fit_r2(base, [ones(size(gapList)), gapList]);
r2Base_inv = fit_r2(base, [ones(size(gapList)), invGap]);
r2Amp_g = fit_r2(amp, [ones(size(gapList)), gapList]);
r2Amp_inv = fit_r2(amp, [ones(size(gapList)), invGap]);
r2Xc_g = fit_r2(xc, [ones(size(gapList)), gapList]);
r2Xc_inv = fit_r2(xc, [ones(size(gapList)), invGap]);
r2Fwhm_g = fit_r2(fwhm, [ones(size(gapList)), gapList]);
r2Fwhm_inv = fit_r2(fwhm, [ones(size(gapList)), invGap]);

modelParamTable = table(gapList, base, amp, xc, fwhm, fitRmse, fitNrmse);
modelFitSummary = table( ...
    mean(fitRmse), mean(fitNrmse), max(fitNrmse), ...
    r2Base_g, r2Base_inv, r2Amp_g, r2Amp_inv, ...
    r2Xc_g, r2Xc_inv, r2Fwhm_g, r2Fwhm_inv, ...
    'VariableNames', {'mean_rmse', 'mean_nrmse', 'max_nrmse', ...
    'R2_base_g', 'R2_base_inv_g', 'R2_amp_g', 'R2_amp_inv_g', ...
    'R2_xc_g', 'R2_xc_inv_g', 'R2_fwhm_g', 'R2_fwhm_inv_g'});

save(fullfile(outDir, 'stage4_standard_shape_model.mat'), ...
    'xiGrid', 'H', 'base', 'amp', 'xc', 'fwhm', ...
    'S_norm', 'S_model', 'modelParamTable', 'modelFitSummary', '-v7.3');

figure('Name', 'Static Gap Step 4 - Standard Shape Model', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for i = 1:nGap
    plot(xiGrid, S_norm(i, :), 'HandleVisibility', 'off');
end
plot(xiGrid, H, 'k-', 'LineWidth', 2.0, 'DisplayName', 'shared H(\xi)');
xlabel('\xi = (x - x_c) / FWHM');
ylabel('Normalized response');
title('Shared standard shape H(\xi)');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
repIdx = unique([1, ceil(nGap/2), nGap]);
colors = lines(numel(repIdx));
for k = 1:numel(repIdx)
    i = repIdx(k);
    plot(xGrid, S(i, :), '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('true g = %.1f', gapList(i)));
    plot(xGrid, S_model(i, :), '--', 'Color', colors(k, :), ...
        'HandleVisibility', 'off');
end
xlabel('Position x (mm)');
ylabel('Capacitance');
title('True waveform vs standard-shape reconstruction');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(gapList, base, 'o-', 'DisplayName', 'baseline');
plot(gapList, amp / max(amp) * max(base), 's-', 'DisplayName', 'amplitude (scaled)');
plot(gapList, fwhm / max(fwhm) * max(base), 'd-', 'DisplayName', 'FWHM (scaled)');
xlabel('Gap g (mm)');
ylabel('Parameter / scaled');
title('Gap-dependent parameters');
legend('Location', 'best', 'Box', 'off');

nexttile;
bar(categorical(string(gapList)), fitNrmse);
xlabel('Gap g (mm)');
ylabel('NRMSE');
title('Reconstruction error by gap');

fprintf('[Static Step 4] Standard-shape model summary:\n');
disp(modelFitSummary);

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
        i = idx(ii);
        xCross(ii) = interp1(y(i:i+1), x(i:i+1), halfLevel, 'linear', 'extrap');
    end
    width = xCross(end) - xCross(1);
end

function r2 = fit_r2(y, X)
    beta = X \ y;
    yHat = X * beta;
    ssRes = sum((y - yHat).^2);
    ssTot = sum((y - mean(y)).^2);
    r2 = 1 - ssRes / max(ssTot, eps);
end
