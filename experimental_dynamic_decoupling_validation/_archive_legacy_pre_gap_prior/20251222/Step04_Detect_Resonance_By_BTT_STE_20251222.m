%% Step04: Detect candidate resonance regions using BTT displacement STE
% This step uses blade vibration displacement as the main response signal.
% It keeps two BTT energy channels:
%   1. transient energy after moving-median bias removal;
%   2. level energy after only global median removal.
% The second channel is important for nearly steady vibration during a
% constant-speed segment, where moving-median removal can suppress the
% response that we actually want to detect.

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
legacyDir = fullfile(thisDir, 'legacy');
addpath(legacyDir);

cfg = Get_20251222_BTT_Config();
caseName = '1000_2500_3500';
caseOutputDir = fullfile(cfg.output_root, caseName);
outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sensorIds = cfg.step3_default_btt_sensor_ids;
selectedBlades = 1:cfg.blades_num;
steWindowSec = 1.20;
gridDtSec = 0.05;
baselineTimeSec = 3.0;
smoothTimeSec = 1.00;
minRegionDurationSec = 1.00;
mergeGapSec = 0.50;
transientWeight = 0.40;
levelWeight = 0.35;
strainWeight = 0.25;
syncFreqBandHz = cfg.step3_align_freq_plot_hz;
orderCandidates = cfg.step3_align_order_candidates(:).';
platformMinDurationSec = 3.0;
platformSlopeLimitRpmPerSec = 15;
platformStrainScoreThreshold = 0.62;
platformBttScoreThreshold = 0.00;

%% 2. Load BTT displacement and remove slow bias
allTime = [];
allRawDisp = [];
allVibDisp = [];
allLevelDisp = [];
allSensor = [];
allBlade = [];

for sid = sensorIds
    vibFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d_vib_final.mat', sid));
    if ~isfile(vibFile)
        warning('Missing BTT displacement file: %s', vibFile);
        continue;
    end

    S = load(vibFile, 'jilublade');
    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 6
        warning('Invalid BTT displacement file: %s', vibFile);
        continue;
    end

    for bid = selectedBlades
        bladeMask = S.jilublade(:, 4) == bid;
        t = S.jilublade(bladeMask, 3);
        y = S.jilublade(bladeMask, 6);
        valid = isfinite(t) & isfinite(y);
        t = t(valid);
        y = y(valid);
        if numel(t) < 10
            continue;
        end

        [t, orderIdx] = sort(t);
        y = y(orderIdx);
        dtBlade = median(diff(t), 'omitnan');
        if ~isfinite(dtBlade) || dtBlade <= 0
            dtBlade = steWindowSec;
        end
        baselinePts = max(5, round(baselineTimeSec / dtBlade));
        if mod(baselinePts, 2) == 0
            baselinePts = baselinePts + 1;
        end
        yBias = movmedian(y, baselinePts, 'omitnan');
        yVib = y - yBias;
        yLevel = y - median(y, 'omitnan');

        allTime = [allTime; t(:)]; %#ok<AGROW>
        allRawDisp = [allRawDisp; y(:)]; %#ok<AGROW>
        allVibDisp = [allVibDisp; yVib(:)]; %#ok<AGROW>
        allLevelDisp = [allLevelDisp; yLevel(:)]; %#ok<AGROW>
        allSensor = [allSensor; sid * ones(numel(t), 1)]; %#ok<AGROW>
        allBlade = [allBlade; bid * ones(numel(t), 1)]; %#ok<AGROW>
    end
end

if isempty(allTime)
    error('No valid BTT displacement points were loaded.');
end

%% 3. Compute BTT short-time energy on a regular time grid
tMin = min(allTime);
tMax = max(allTime);
gridTime = (tMin:gridDtSec:tMax).';
nGrid = numel(gridTime);
bttSTE = nan(nGrid, 1);
bttRMS = nan(nGrid, 1);
bttLevelSTE = nan(nGrid, 1);
bttLevelRMS = nan(nGrid, 1);
bttP95Abs = nan(nGrid, 1);
pointCount = zeros(nGrid, 1);

halfWin = 0.5 * steWindowSec;
for i = 1:nGrid
    mask = abs(allTime - gridTime(i)) <= halfWin;
    pointCount(i) = nnz(mask);
    if pointCount(i) < 8
        continue;
    end

    yWin = allVibDisp(mask);
    yWin = yWin(isfinite(yWin));
    yLevelWin = allLevelDisp(mask);
    yLevelWin = yLevelWin(isfinite(yLevelWin));
    if numel(yWin) < 8 || numel(yLevelWin) < 8
        continue;
    end

    bttSTE(i) = mean(yWin.^2);
    bttRMS(i) = sqrt(bttSTE(i));
    bttLevelSTE(i) = mean(yLevelWin.^2);
    bttLevelRMS(i) = sqrt(bttLevelSTE(i));
    sortedAbs = sort(abs(yWin));
    idx95 = max(1, min(numel(sortedAbs), round(0.95 * numel(sortedAbs))));
    bttP95Abs(i) = sortedAbs(idx95);
end

smoothPts = max(3, round(smoothTimeSec / gridDtSec));
bttSTEsm = movmean(bttSTE, smoothPts, 'omitnan');
bttRMSsm = sqrt(bttSTEsm);
bttLevelSTEsm = movmean(bttLevelSTE, smoothPts, 'omitnan');
bttLevelRMSsm = sqrt(bttLevelSTEsm);

finiteBtt = isfinite(bttSTEsm);
finiteLevelBtt = isfinite(bttLevelSTEsm);
if nnz(finiteBtt) < 5 || nnz(finiteLevelBtt) < 5
    error('Too few valid BTT STE grid points.');
end
bttLog = log1p(bttSTEsm);
bttLow = prctile(bttLog(finiteBtt), 5);
bttHigh = prctile(bttLog(finiteBtt), 95);
bttScore = (bttLog - bttLow) ./ max(bttHigh - bttLow, eps);
bttScore = min(max(bttScore, 0), 1);

bttLevelLog = log1p(bttLevelSTEsm);
levelLow = prctile(bttLevelLog(finiteLevelBtt), 5);
levelHigh = prctile(bttLevelLog(finiteLevelBtt), 95);
bttLevelScore = (bttLevelLog - levelLow) ./ max(levelHigh - levelLow, eps);
bttLevelScore = min(max(bttLevelScore, 0), 1);

%% 4. Project multi-order strain energy onto the BTT time grid
alignFile = fullfile(caseOutputDir, 'Step3_Spectrum_RPM_Alignment_20251222.mat');
if ~isfile(alignFile)
    error('Run Step03 first. Missing file: %s', alignFile);
end
A = load(alignFile, 'result');
alignResult = A.result;

strainOrderAmp = nan(numel(alignResult.stft_time_s), numel(orderCandidates));
for j = 1:numel(orderCandidates)
    orderNow = orderCandidates(j);
    orderLineHz = orderNow * alignResult.frot_strain_hz(:);
    for i = 1:numel(alignResult.stft_time_s)
        fLine = orderLineHz(i);
        if ~isfinite(fLine) || fLine < syncFreqBandHz(1) || fLine > syncFreqBandHz(2)
            continue;
        end
        [~, fIdx] = min(abs(alignResult.stft_freq_hz - fLine));
        strainOrderAmp(i, j) = alignResult.stft_amp(fIdx, i);
    end
end

strainOrderDb = 20 * log10(strainOrderAmp + eps);
validStrain = isfinite(strainOrderDb);
strainScoreByOrder = zeros(nGrid, numel(orderCandidates));
if nnz(validStrain) >= 5
    svAll = strainOrderDb(validStrain);
    sMin = prctile(svAll, 5);
    sMax = prctile(svAll, 95);
    for j = 1:numel(orderCandidates)
        validOrder = isfinite(strainOrderDb(:, j));
        if nnz(validOrder) < 5
            continue;
        end
        st = alignResult.stft_time_s(validOrder);
        svNorm = (strainOrderDb(validOrder, j) - sMin) ./ max(sMax - sMin, eps);
        svNorm = min(max(svNorm, 0), 1);
        strainScoreByOrder(:, j) = interp1(st - alignResult.best_tau_sec, svNorm, ...
            gridTime, 'linear', 0);
    end
end

rpmGrid = interp1(alignResult.rpm_time_s, alignResult.rpm_values, gridTime, 'linear', 'extrap');
syncFreqByOrder = (rpmGrid / 60) * orderCandidates;
freqGateByOrder = syncFreqByOrder >= syncFreqBandHz(1) & syncFreqByOrder <= syncFreqBandHz(2);
strainScoreByOrder(~freqGateByOrder) = 0;
[strainScore, selectedOrderIdx] = max(strainScoreByOrder, [], 2);
selectedOrderGrid = orderCandidates(selectedOrderIdx).';
syncFreqGrid = nan(nGrid, 1);
for i = 1:nGrid
    syncFreqGrid(i) = syncFreqByOrder(i, selectedOrderIdx(i));
end
strainEvidenceGate = strainScore > 0.05;
freqGate = any(freqGateByOrder, 2) & strainEvidenceGate;

%% 5. Build combined resonance score and candidate regions
combinedScore = transientWeight * bttScore + levelWeight * bttLevelScore + strainWeight * strainScore;
combinedScore(~isfinite(combinedScore)) = 0;
validScore = combinedScore(pointCount >= 8 & freqGate);
validScore = sort(validScore(isfinite(validScore)));
if isempty(validScore)
    error('No valid combined score points.');
end
q65 = validScore(max(1, min(numel(validScore), round(0.65 * numel(validScore)))));
scoreThreshold = max(0.30, min(0.70, q65));
strainSustainedThreshold = 0.72;
bttModerateThreshold = 0.22;
scoreMask = combinedScore >= scoreThreshold;
steadyMask = strainScore >= strainSustainedThreshold & ...
    max(bttScore, bttLevelScore) >= bttModerateThreshold;

rpmSmooth = movmedian(rpmGrid, max(3, round(1.0 / gridDtSec)), 'omitnan');
rpmSlope = gradient(rpmSmooth, gridDtSec);
platformRawMask = abs(rpmSlope) <= platformSlopeLimitRpmPerSec & pointCount >= 8 & freqGate;
platformAcceptMask = false(nGrid, 1);

platformEdge = diff([false; platformRawMask; false]);
platformStartIdx = find(platformEdge == 1);
platformEndIdx = find(platformEdge == -1) - 1;
for i = 1:numel(platformStartIdx)
    idx = platformStartIdx(i):platformEndIdx(i);
    platformDuration = gridTime(platformEndIdx(i)) - gridTime(platformStartIdx(i));
    if platformDuration < platformMinDurationSec
        continue;
    end
    platformStrain = mean(strainScore(idx), 'omitnan');
    platformBtt = max(mean(bttScore(idx), 'omitnan'), mean(bttLevelScore(idx), 'omitnan'));
    if platformStrain >= platformStrainScoreThreshold && platformBtt >= platformBttScoreThreshold
        platformAcceptMask(idx) = true;
    end
end

transientMask = (scoreMask | steadyMask) & ~platformRawMask;
activeMask = (transientMask | platformAcceptMask) & pointCount >= 8 & freqGate;

edge = diff([false; activeMask; false]);
startIdx = find(edge == 1);
endIdx = find(edge == -1) - 1;

if numel(startIdx) >= 2
    keepStart = startIdx(1);
    newStart = [];
    newEnd = [];
    for i = 1:numel(startIdx)
        if i == numel(startIdx)
            newStart = [newStart; keepStart]; %#ok<AGROW>
            newEnd = [newEnd; endIdx(i)]; %#ok<AGROW>
        else
            gapSec = gridTime(startIdx(i + 1)) - gridTime(endIdx(i));
            if gapSec > mergeGapSec
                newStart = [newStart; keepStart]; %#ok<AGROW>
                newEnd = [newEnd; endIdx(i)]; %#ok<AGROW>
                keepStart = startIdx(i + 1);
            end
        end
    end
    startIdx = newStart;
    endIdx = newEnd;
end

regionStart = [];
regionEnd = [];
regionDuration = [];
regionPeakTime = [];
regionPeakScore = [];
regionMeanScore = [];
regionPeakTransientRMS = [];
regionPeakLevelRMS = [];
regionPeakRPM = [];
regionPeakFreqHz = [];
regionPeakOrder = [];
regionEvidence = strings(0, 1);

for i = 1:numel(startIdx)
    durationSec = gridTime(endIdx(i)) - gridTime(startIdx(i));
    if durationSec < minRegionDurationSec
        continue;
    end

    localIdx = startIdx(i):endIdx(i);
    peakSearchIdx = localIdx(freqGate(localIdx));
    if isempty(peakSearchIdx)
        peakSearchIdx = localIdx;
    end
    [peakScore, relPeakIdx] = max(combinedScore(peakSearchIdx));
    peakIdx = peakSearchIdx(relPeakIdx);
    rpmPeak = rpmGrid(peakIdx);
    freqPeak = syncFreqGrid(peakIdx);
    orderPeak = selectedOrderGrid(peakIdx);

    regionStart = [regionStart; gridTime(startIdx(i))]; %#ok<AGROW>
    regionEnd = [regionEnd; gridTime(endIdx(i))]; %#ok<AGROW>
    regionDuration = [regionDuration; durationSec]; %#ok<AGROW>
    regionPeakTime = [regionPeakTime; gridTime(peakIdx)]; %#ok<AGROW>
    regionPeakScore = [regionPeakScore; peakScore]; %#ok<AGROW>
    regionMeanScore = [regionMeanScore; mean(combinedScore(localIdx), 'omitnan')]; %#ok<AGROW>
    regionPeakTransientRMS = [regionPeakTransientRMS; bttRMSsm(peakIdx)]; %#ok<AGROW>
    regionPeakLevelRMS = [regionPeakLevelRMS; bttLevelRMSsm(peakIdx)]; %#ok<AGROW>
    regionPeakRPM = [regionPeakRPM; rpmPeak]; %#ok<AGROW>
    regionPeakFreqHz = [regionPeakFreqHz; freqPeak]; %#ok<AGROW>
    regionPeakOrder = [regionPeakOrder; orderPeak]; %#ok<AGROW>
    if any(platformAcceptMask(localIdx))
        if mean(max(bttScore(localIdx), bttLevelScore(localIdx)), 'omitnan') >= 0.12
            regionEvidence = [regionEvidence; "strain_btt_platform"]; %#ok<AGROW>
        else
            regionEvidence = [regionEvidence; "strain_platform_btt_subwindow"]; %#ok<AGROW>
        end
    elseif any(steadyMask(localIdx))
        regionEvidence = [regionEvidence; "steady_order"]; %#ok<AGROW>
    else
        regionEvidence = [regionEvidence; "transient_peak"]; %#ok<AGROW>
    end
end

regionTable = table(regionStart, regionEnd, regionDuration, regionPeakTime, ...
    regionPeakScore, regionMeanScore, regionPeakTransientRMS, ...
    regionPeakLevelRMS, regionPeakRPM, regionPeakOrder, regionPeakFreqHz, ...
    regionEvidence);

if ~isempty(regionTable)
    regionTable = sortrows(regionTable, 'regionPeakScore', 'descend');
end

%% 6. Save outputs
scoreTable = table(gridTime, bttSTEsm, bttRMSsm, bttScore, ...
    bttLevelSTEsm, bttLevelRMSsm, bttLevelScore, strainScore, ...
    selectedOrderGrid, combinedScore, pointCount, rpmGrid, syncFreqGrid, ...
    freqGate, rpmSlope, platformRawMask, platformAcceptMask, steadyMask);

matFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.mat');
csvFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.csv');
save(matFile, 'scoreTable', 'regionTable', 'allTime', 'allRawDisp', ...
    'allVibDisp', 'allLevelDisp', 'allSensor', 'allBlade', 'alignResult', ...
    'orderCandidates', 'strainScoreByOrder', 'syncFreqByOrder', '-v7.3');
writetable(regionTable, csvFile);

%% 7. Visualization
figure('Name', '20251222 Step04 BTT STE resonance regions', ...
    'Color', 'w', 'Position', [70, 70, 1360, 940], 'NumberTitle', 'off');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on; grid on; box on;
for sid = sensorIds
    mask = allSensor == sid;
    plot(allTime(mask), allVibDisp(mask), '.', 'MarkerSize', 3, ...
        'DisplayName', sprintf('CH%d', sid));
end
yl = ylim;
for i = 1:height(regionTable)
    patch([regionTable.regionStart(i), regionTable.regionEnd(i), ...
        regionTable.regionEnd(i), regionTable.regionStart(i)], ...
        [yl(1), yl(1), yl(2), yl(2)], [1.00, 0.85, 0.40], ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
uistack(findobj(gca, 'Type', 'line'), 'top');
ylabel('BTT residual disp. (mm)');
title('Blade vibration displacement after slow-bias removal');
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(gridTime, bttRMSsm, 'Color', [0.10, 0.35, 0.75], 'LineWidth', 1.2, ...
    'DisplayName', 'transient RMS');
plot(gridTime, bttLevelRMSsm, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0, ...
    'DisplayName', 'level RMS');
for i = 1:height(regionTable)
    xline(regionTable.regionPeakTime(i), '--', sprintf('%.2f s', regionTable.regionPeakTime(i)), ...
        'HandleVisibility', 'off');
end
ylabel('BTT RMS (mm)');
title(sprintf('Short-time RMS, window = %.2f s', steWindowSec));
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(gridTime, bttScore, 'LineWidth', 1.1, 'DisplayName', 'transient BTT score');
plot(gridTime, bttLevelScore, 'LineWidth', 1.1, 'DisplayName', 'level BTT score');
plot(gridTime, strainScore, 'LineWidth', 1.1, 'DisplayName', 'strain order score');
plot(gridTime, combinedScore, 'k-', 'LineWidth', 1.5, 'DisplayName', 'combined score');
yline(scoreThreshold, 'r--', 'DisplayName', 'threshold');
ylabel('Score');
title(sprintf('Score: %.0f%% transient + %.0f%% level + %.0f%% strain, gated %.0f-%.0f Hz', ...
    100*transientWeight, 100*levelWeight, 100*strainWeight, ...
    syncFreqBandHz(1), syncFreqBandHz(2)));
legend('Location', 'best');

nexttile;
hold on; grid on; box on;
plot(alignResult.rpm_time_s, alignResult.rpm_values, 'k-', 'LineWidth', 1.0);
for i = 1:height(regionTable)
    xline(regionTable.regionStart(i), ':', 'Color', [0.64, 0.08, 0.18], 'HandleVisibility', 'off');
    xline(regionTable.regionEnd(i), ':', 'Color', [0.64, 0.08, 0.18], 'HandleVisibility', 'off');
end
xlabel('BTT time (s)');
ylabel('RPM');
title(sprintf('Candidate resonance regions, strain orders = %d-%dX, %.0f-%.0f Hz gate', ...
    min(orderCandidates), max(orderCandidates), syncFreqBandHz(1), syncFreqBandHz(2)));

figFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.png');
saveas(gcf, figFile);

%% 8. Validation visualization: aligned strain STFT
tStrainBtt = alignResult.stft_time_s(:) - alignResult.best_tau_sec;
freqMaskPlot = alignResult.stft_freq_hz >= syncFreqBandHz(1) & ...
    alignResult.stft_freq_hz <= syncFreqBandHz(2);
strainAmpDbPlot = 20 * log10(alignResult.stft_amp(freqMaskPlot, :) + eps);

figure('Name', '20251222 Step04 aligned strain STFT validation', ...
    'Color', 'w', 'Position', [90, 90, 1360, 900], 'NumberTitle', 'off');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(tStrainBtt, alignResult.stft_freq_hz(freqMaskPlot), strainAmpDbPlot);
axis xy; grid on; box on;
colormap(gca, turbo);
colorbar;
hold on;
yl = ylim;
for i = 1:height(regionTable)
    patch([regionTable.regionStart(i), regionTable.regionEnd(i), ...
        regionTable.regionEnd(i), regionTable.regionStart(i)], ...
        [yl(1), yl(1), yl(2), yl(2)], [1.00, 0.85, 0.30], ...
        'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for i = 1:height(regionTable)
    tLine = linspace(regionTable.regionStart(i), regionTable.regionEnd(i), 160);
    rpmLine = interp1(alignResult.rpm_time_s, alignResult.rpm_values, ...
        tLine, 'linear', 'extrap');
    fLine = regionTable.regionPeakOrder(i) * rpmLine / 60;
    validLine = fLine >= syncFreqBandHz(1) & fLine <= syncFreqBandHz(2);
    if any(validLine)
        plot(tLine(validLine), fLine(validLine), 'w-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%dX region %d', regionTable.regionPeakOrder(i), i));
        plot(regionTable.regionPeakTime(i), regionTable.regionPeakFreqHz(i), ...
            'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
    end
end
xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
title('Aligned strain STFT with Step04 resonance regions and selected order tracks');
legend('Location', 'eastoutside');

nexttile;
hold on; grid on; box on;
plot(gridTime, strainScore, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0, ...
    'DisplayName', 'multi-order strain score');
plot(gridTime, combinedScore, 'k-', 'LineWidth', 1.2, 'DisplayName', 'combined Step04 score');
yline(scoreThreshold, 'r--', 'DisplayName', 'threshold');
for i = 1:height(regionTable)
    xline(regionTable.regionStart(i), ':', 'Color', [0.45, 0.45, 0.45], 'HandleVisibility', 'off');
    xline(regionTable.regionEnd(i), ':', 'Color', [0.45, 0.45, 0.45], 'HandleVisibility', 'off');
end
ylabel('Score');
title('Strain-side confirmation score after time alignment');
legend('Location', 'best');

nexttile;
yyaxis left;
plot(gridTime, selectedOrderGrid, 'Color', [0.10, 0.35, 0.75], 'LineWidth', 1.0);
ylabel('Selected order');
yyaxis right;
plot(gridTime, syncFreqGrid, 'Color', [0.20, 0.55, 0.25], 'LineWidth', 1.0);
ylabel('Selected frequency (Hz)');
grid on; box on;
xlabel('BTT-aligned time (s)');
title('Order selected from strain STFT and its synchronous frequency');

stftFigFile = fullfile(outDir, 'Step04_Aligned_Strain_STFT_Validation_20251222.png');
saveas(gcf, stftFigFile);

fprintf('\nStep04 complete. Candidate resonance regions:\n');
disp(regionTable);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', matFile, csvFile, figFile, stftFigFile);
