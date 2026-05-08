%% Step 3: analyze geometric feature laws versus gap
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'cfg', 'gapList', 'xGrid', 'S');

nGap = numel(gapList);
peakVal = zeros(nGap, 1);
peakX = zeros(nGap, 1);
areaVal = zeros(nGap, 1);
fwhmVal = zeros(nGap, 1);
widthSig = zeros(nGap, numel(cfg.levelList));

for i = 1:nGap
    y = S(i, :).';
    [peakVal(i), idxPeak] = max(y);
    peakX(i) = xGrid(idxPeak);
    areaVal(i) = trapz(xGrid, y - min(y));
    fwhmVal(i) = calc_fwhm(xGrid, y);
    widthSig(i, :) = calc_width_signature(xGrid, y, cfg.levelList);
end

r2Peak_g = fit_r2(peakVal, [ones(size(gapList)), gapList]);
r2Peak_inv = fit_r2(peakVal, [ones(size(gapList)), 1 ./ gapList]);
r2Fwhm_g = fit_r2(fwhmVal, [ones(size(gapList)), gapList]);
r2Fwhm_inv = fit_r2(fwhmVal, [ones(size(gapList)), 1 ./ gapList]);
r2Area_g = fit_r2(areaVal, [ones(size(gapList)), gapList]);
r2Area_inv = fit_r2(areaVal, [ones(size(gapList)), 1 ./ gapList]);

widthR2_g = zeros(numel(cfg.levelList), 1);
widthR2_inv = zeros(numel(cfg.levelList), 1);
for iq = 1:numel(cfg.levelList)
    widthR2_g(iq) = fit_r2(widthSig(:, iq), [ones(size(gapList)), gapList]);
    widthR2_inv(iq) = fit_r2(widthSig(:, iq), [ones(size(gapList)), 1 ./ gapList]);
end

featureTable = table(gapList, peakVal, peakX, areaVal, fwhmVal);
featureLawSummary = table( ...
    r2Peak_g, r2Peak_inv, r2Area_g, r2Area_inv, r2Fwhm_g, r2Fwhm_inv, ...
    mean(widthR2_g), mean(widthR2_inv), ...
    'VariableNames', {'R2_peak_linear_g', 'R2_peak_linear_inv_g', ...
    'R2_area_linear_g', 'R2_area_linear_inv_g', ...
    'R2_fwhm_linear_g', 'R2_fwhm_linear_inv_g', ...
    'mean_R2_width_linear_g', 'mean_R2_width_linear_inv_g'});

save(fullfile(outDir, 'stage3_feature_laws.mat'), ...
    'featureTable', 'featureLawSummary', 'widthSig', 'widthR2_g', 'widthR2_inv', '-v7.3');

figure('Name', 'Static Gap Step 3 - Feature Laws', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(gapList, peakVal, 'o-', 'DisplayName', 'peak');
plot(gapList, areaVal / max(areaVal) * max(peakVal), 's-', 'DisplayName', 'area (scaled)');
plot(gapList, fwhmVal / max(fwhmVal) * max(peakVal), 'd-', 'DisplayName', 'FWHM (scaled)');
xlabel('Gap g (mm)');
ylabel('Feature value / scaled');
title('Main geometric features vs gap');
legend('Location', 'best', 'Box', 'off');

nexttile;
imagesc(cfg.levelList, gapList, widthSig);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Normalized level q');
ylabel('Gap g (mm)');
title('Multi-level width signature w(q,g)');

nexttile; hold on;
plot(cfg.levelList, widthR2_g, 'o-', 'DisplayName', 'linear in g');
plot(cfg.levelList, widthR2_inv, 's-', 'DisplayName', 'linear in 1/g');
xlabel('Normalized level q');
ylabel('R^2');
title('Width-signature law quality');
legend('Location', 'best', 'Box', 'off');

nexttile;
bar(categorical({'peak:g', 'peak:1/g', 'area:g', 'area:1/g', 'FWHM:g', 'FWHM:1/g'}), ...
    [r2Peak_g, r2Peak_inv, r2Area_g, r2Area_inv, r2Fwhm_g, r2Fwhm_inv]);
ylabel('R^2');
title('Feature-law comparison');

fprintf('[Static Step 3] Feature-law summary:\n');
disp(featureLawSummary);

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

function widthSig = calc_width_signature(x, y, levelList)
    yNorm = (y - min(y)) / max(range(y), eps);
    widthSig = nan(numel(levelList), 1);
    for iq = 1:numel(levelList)
        lv = levelList(iq);
        above = yNorm >= lv;
        idx = find(diff(above) ~= 0);
        if numel(idx) < 2
            continue;
        end
        xCross = zeros(numel(idx), 1);
        for ii = 1:numel(idx)
            i = idx(ii);
            xCross(ii) = interp1(yNorm(i:i+1), x(i:i+1), lv, 'linear', 'extrap');
        end
        widthSig(iq) = xCross(end) - xCross(1);
    end
end

function r2 = fit_r2(y, X)
    beta = X \ y;
    yHat = X * beta;
    ssRes = sum((y - yHat).^2);
    ssTot = sum((y - mean(y)).^2);
    r2 = 1 - ssRes / max(ssTot, eps);
end

