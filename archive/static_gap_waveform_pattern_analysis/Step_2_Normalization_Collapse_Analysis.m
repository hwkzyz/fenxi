%% Step 2: analyze whether the waveform family collapses after normalization
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S');

nGap = numel(gapList);
baseline = min(S, [], 2);
amp = max(S - baseline, [], 2);

Xpeak = zeros(nGap, 1);
Fwhm = zeros(nGap, 1);
S_peak = zeros(size(S));
S_peak_amp = zeros(size(S));
xiGrid = linspace(-3, 3, numel(xGrid))';
S_norm = nan(size(S));

for i = 1:nGap
    y = S(i, :).';
    [~, idxPeak] = max(y);
    Xpeak(i) = xGrid(idxPeak);
    Fwhm(i) = calc_fwhm(xGrid, y);
    xPeakCentered = xGrid - Xpeak(i);
    S_peak(i, :) = interp1(xPeakCentered, y, xGrid, 'pchip', NaN);
    yAmp = (y - baseline(i)) / max(amp(i), eps);
    S_peak_amp(i, :) = interp1(xPeakCentered, yAmp, xGrid, 'pchip', NaN);
    xi = (xGrid - Xpeak(i)) / max(Fwhm(i), eps);
    S_norm(i, :) = interp1(xi, yAmp, xiGrid, 'pchip', NaN);
end

stats = struct();
stats.raw_rank1_energy = rank1_energy(S);
stats.peak_rank1_energy = rank1_energy(fillmissing(S_peak, 'linear', 2, 'EndValues', 'nearest'));
stats.peak_amp_rank1_energy = rank1_energy(fillmissing(S_peak_amp, 'linear', 2, 'EndValues', 'nearest'));
stats.norm_rank1_energy = rank1_energy(fillmissing(S_norm, 'linear', 2, 'EndValues', 'nearest'));
stats.raw_mean_pair_rmse = mean_pair_rmse(S);
stats.peak_mean_pair_rmse = mean_pair_rmse(S_peak);
stats.peak_amp_mean_pair_rmse = mean_pair_rmse(S_peak_amp);
stats.norm_mean_pair_rmse = mean_pair_rmse(S_norm);

save(fullfile(outDir, 'stage2_normalization_collapse.mat'), ...
    'Xpeak', 'Fwhm', 'S_peak', 'S_peak_amp', 'xiGrid', 'S_norm', 'stats', '-v7.3');

figure('Name', 'Static Gap Step 2 - Normalization Collapse', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for i = 1:nGap
    plot(xGrid, S(i, :), 'DisplayName', sprintf('g = %.1f', gapList(i)));
end
xlabel('Raw x (mm)');
ylabel('Capacitance');
title('Raw waveform family');
legend('Location', 'eastoutside', 'Box', 'off');

nexttile; hold on;
for i = 1:nGap
    plot(xGrid, S_peak_amp(i, :));
end
xlabel('Peak-centered x (mm)');
ylabel('Amplitude-normalized');
title('Peak-centered + amplitude-normalized');

nexttile; hold on;
for i = 1:nGap
    plot(xiGrid, S_norm(i, :));
end
xlabel('\xi = (x - x_c) / FWHM');
ylabel('Normalized response');
title('Width-normalized collapse');

nexttile;
bar(categorical({'raw', 'peak', 'peak+amp', 'peak+amp+width'}), ...
    [stats.raw_rank1_energy, stats.peak_rank1_energy, ...
     stats.peak_amp_rank1_energy, stats.norm_rank1_energy]);
ylabel('Rank-1 energy fraction');
title('Low-rank collapse after normalization');

fprintf('[Static Step 2] Collapse statistics:\n');
disp(struct2table(stats));

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

function frac = rank1_energy(M)
    M = fillmissing(M, 'linear', 2, 'EndValues', 'nearest');
    M = M - mean(M, 2);
    [~, Ssvd, ~] = svd(M, 'econ');
    s = diag(Ssvd).^2;
    frac = s(1) / max(sum(s), eps);
end

function val = mean_pair_rmse(M)
    M = fillmissing(M, 'linear', 2, 'EndValues', 'nearest');
    c = 0;
    acc = 0;
    for i = 1:size(M, 1)-1
        for j = i+1:size(M, 1)
            valid = isfinite(M(i, :)) & isfinite(M(j, :));
            acc = acc + sqrt(mean((M(i, valid) - M(j, valid)).^2));
            c = c + 1;
        end
    end
    val = acc / max(c, 1);
end

