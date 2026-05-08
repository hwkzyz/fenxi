%% Analyze how many distinct no-vibration gap states are needed
% We treat the available 7-gap dataset as the full static truth library.
% Then we subsample K gap states as "available calibration data" and test
% how well the remaining unseen gap states can be predicted.
%
% The pipeline is:
%   1) register all no-vibration waveforms to the refined common domain;
%   2) for each subset size K, choose all combinations of K gap states;
%   3) interpolate amplitude and normalized shape across gap using only
%      the chosen subset;
%   4) evaluate prediction error on the held-out gap states.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = resolve_data_file(scriptDir);
outDir = fullfile(scriptDir, 'gap_count_requirement_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Figure defaults
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.15, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% Refined registration settings
levelMap = linspace(0.005, 0.995, 199)';
nU = 1401;
interpMethod = 'linear';

%% Load and register all truth waveforms
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);
gRef = gapList(refIdx);

curve = preprocess_curves(gapList, xCell, yCell, levelMap);
uGrid = make_reference_grid(curve(refIdx), levelMap, nU);
Cmat = register_all_curves(curve, refIdx, uGrid, interpMethod);
Atrue = max(Cmat, [], 2);
Hmat = Cmat ./ max(Atrue, eps);

%% Enumerate subset calibration experiments
allRows = [];
for k = 2:(nGap - 1)
    combos = nchoosek(1:nGap, k);
    for ic = 1:size(combos, 1)
        subset = combos(ic, :);
        holdout = setdiff(1:nGap, subset);
        if isempty(holdout)
            continue;
        end

        subsetGaps = gapList(subset);
        Acal = Atrue(subset);
        Hcal = Hmat(subset, :);

        rmseList = nan(numel(holdout), 1);
        shapeRmseList = nan(numel(holdout), 1);
        ampRelErrList = nan(numel(holdout), 1);

        for ih = 1:numel(holdout)
            j = holdout(ih);
            gTarget = gapList(j);
            Ahat = interp1(subsetGaps, Acal, gTarget, 'pchip', 'extrap');
            Hhat = interp1(subsetGaps, Hcal, gTarget, 'pchip', 'extrap');
            Chat = Ahat .* Hhat;

            Ctrue = Cmat(j, :).';
            Htrue = Hmat(j, :).';
            rmseList(ih) = sqrt(mean((Chat(:) - Ctrue).^2, 'omitnan'));
            shapeRmseList(ih) = sqrt(mean((Hhat(:) - Htrue).^2, 'omitnan'));
            ampRelErrList(ih) = abs(Ahat - Atrue(j)) / max(Atrue(j), eps);
        end

        row.subset_size = k;
        row.subset_idx = join(string(subset), '-');
        row.subset_gaps_mm = join(compose('%.1f', subsetGaps), ',');
        row.n_holdout = numel(holdout);
        row.include_min_gap = double(any(subset == 1));
        row.include_ref_gap = double(any(subset == refIdx));
        row.include_both_endpoints = double(any(subset == 1) && any(subset == refIdx));
        row.max_gap_spacing_mm = max(diff(subsetGaps));
        row.mean_holdout_rmse_pF = mean(rmseList);
        row.max_holdout_rmse_pF = max(rmseList);
        row.mean_holdout_shape_rmse = mean(shapeRmseList);
        row.max_holdout_shape_rmse = max(shapeRmseList);
        row.mean_holdout_amp_rel_err = mean(ampRelErrList);
        row.max_holdout_amp_rel_err = max(ampRelErrList);
        allRows = [allRows; row]; %#ok<AGROW>
    end
end

detailTable = struct2table(allRows);
writetable(detailTable, fullfile(outDir, 'gap_subset_detail.csv'));

%% Build subset-size summaries
summaryRows = [];
groupDefs = { ...
    struct('name', "all_subsets", 'mask', true(height(detailTable), 1)), ...
    struct('name', "include_ref_gap", 'mask', logical(detailTable.include_ref_gap)), ...
    struct('name', "include_both_endpoints", 'mask', logical(detailTable.include_both_endpoints))};

for ig = 1:numel(groupDefs)
    group = groupDefs{ig};
    for k = unique(detailTable.subset_size).'
        mask = group.mask & detailTable.subset_size == k;
        if ~any(mask)
            continue;
        end
        T = detailTable(mask, :);
        [~, bestIdx] = min(T.mean_holdout_rmse_pF);
        rmseVals = sort(T.mean_holdout_rmse_pF);
        medVal = median(rmseVals);

        row.group_name = group.name;
        row.subset_size = k;
        row.n_combinations = height(T);
        row.best_mean_holdout_rmse_pF = T.mean_holdout_rmse_pF(bestIdx);
        row.best_max_holdout_rmse_pF = T.max_holdout_rmse_pF(bestIdx);
        row.best_subset_gaps_mm = string(T.subset_gaps_mm(bestIdx));
        row.mean_of_mean_holdout_rmse_pF = mean(T.mean_holdout_rmse_pF);
        row.median_of_mean_holdout_rmse_pF = medVal;
        row.worst_mean_holdout_rmse_pF = max(T.mean_holdout_rmse_pF);
        row.mean_max_gap_spacing_mm = mean(T.max_gap_spacing_mm);
        summaryRows = [summaryRows; row]; %#ok<AGROW>
    end
end

summaryTable = struct2table(summaryRows);
writetable(summaryTable, fullfile(outDir, 'gap_subset_summary.csv'));

%% Recommended minimum count under simple practical thresholds
thresholds = [1e-3; 5e-4; 1e-4; 5e-5];
thrGroup = strings(numel(thresholds), 1);
thrValue = nan(numel(thresholds), 1);
thrMinK = nan(numel(thresholds), 1);
thrSubset = strings(numel(thresholds), 1);
targetGroup = "include_both_endpoints";
for it = 1:numel(thresholds)
    thr = thresholds(it);
    mask = summaryTable.group_name == targetGroup & ...
           summaryTable.best_mean_holdout_rmse_pF <= thr;
    thrGroup(it) = targetGroup;
    thrValue(it) = thr;
    if any(mask)
        kMin = min(summaryTable.subset_size(mask));
        thrMinK(it) = kMin;
        thrSubset(it) = string(summaryTable.best_subset_gaps_mm( ...
            find(mask & summaryTable.subset_size == kMin, 1, 'first')));
    else
        thrMinK(it) = NaN;
        thrSubset(it) = "";
    end
end
recommendationTable = table(thrGroup, thrValue, thrMinK, thrSubset, ...
    'VariableNames', {'group_name', 'threshold_rmse_pF', ...
    'min_subset_size', 'best_subset_gaps_mm'});
writetable(recommendationTable, fullfile(outDir, 'gap_count_recommendation.csv'));

%% Figures
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 9]);
tl1 = tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
groupsToPlot = ["all_subsets", "include_ref_gap", "include_both_endpoints"];
colors = lines(numel(groupsToPlot));
for ig = 1:numel(groupsToPlot)
    mask = summaryTable.group_name == groupsToPlot(ig);
    plot(summaryTable.subset_size(mask), summaryTable.best_mean_holdout_rmse_pF(mask), '-o', ...
        'Color', colors(ig, :), 'MarkerFaceColor', colors(ig, :), ...
        'DisplayName', strrep(groupsToPlot(ig), '_', ' '));
end
xlabel('Number of available gap states K');
ylabel('Best mean holdout RMSE (pF)');
title('(a) Best achievable prediction error');
legend('Location', 'northeast', 'FontSize', 7.0);

nexttile; hold on;
mask = summaryTable.group_name == "include_both_endpoints";
plot(summaryTable.subset_size(mask), summaryTable.mean_of_mean_holdout_rmse_pF(mask), '-o', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], ...
    'DisplayName', 'Mean over endpoint-covered subsets');
plot(summaryTable.subset_size(mask), summaryTable.worst_mean_holdout_rmse_pF(mask), '-s', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], ...
    'DisplayName', 'Worst over endpoint-covered subsets');
xlabel('Number of available gap states K');
ylabel('Mean holdout RMSE (pF)');
title('(b) Practical robustness with endpoint coverage');
legend('Location', 'northeast', 'FontSize', 7.0);

title(tl1, sprintf('Required number of no-vibration gap states, reference state %.1f mm', gRef), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_gap_count_summary.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_gap_count_summary.pdf'), 'ContentType', 'vector');

fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

bestEndpointRows = summaryTable(summaryTable.group_name == "include_both_endpoints", :);
nexttile;
bar(categorical(string(bestEndpointRows.subset_size)), bestEndpointRows.best_mean_holdout_rmse_pF, ...
    0.65, 'FaceColor', [0.74 0.86 0.66], 'EdgeColor', 'k', 'LineWidth', 0.5);
xlabel('Subset size K');
ylabel('Best mean holdout RMSE (pF)');
title('(a) Best endpoint-covered subsets');

nexttile;
text(0.02, 0.95, format_recommendations(recommendationTable), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 8);
axis off;
title('(b) Threshold-based recommendations');

title(tl2, 'Gap-count recommendation summary', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_gap_count_recommendation.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_gap_count_recommendation.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Required number of no-vibration gap states ===\n');
fprintf('Data file: %s\n', dataFile);
fprintf('Reference state gap = %.1f mm\n\n', gRef);
disp(summaryTable(:, {'group_name', 'subset_size', 'best_mean_holdout_rmse_pF', ...
    'mean_of_mean_holdout_rmse_pF', 'best_subset_gaps_mm'}));
fprintf('\nThreshold-based recommendations:\n');
disp(recommendationTable);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function dataFile = resolve_data_file(scriptDir)
    candidateList = { ...
        fullfile(scriptDir, '直叶片2mm_不同间隙.txt'), ...
        fullfile(scriptDir, '鐩村彾鐗?mm_涓嶅悓闂撮殭.txt')};
    for i = 1:numel(candidateList)
        if isfile(candidateList{i})
            dataFile = candidateList{i};
            return;
        end
    end
    error('No supported stacked-curve data file was found in: %s', scriptDir);
end

function [gapList, xCell, yCell] = load_stacked_curves(filePath)
    data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
    data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
    x = data(:, 1);
    y = data(:, 2);
    breakIdx = [find(diff(x) < 0); numel(x)];
    startIdx = [1; breakIdx(1:end-1) + 1];
    nCurve = numel(breakIdx);
    if nCurve == 7
        gapList = (0.2:0.2:1.4)';
    else
        gapList = (1:nCurve)';
    end
    xCell = cell(nCurve, 1);
    yCell = cell(nCurve, 1);
    for i = 1:nCurve
        idx = startIdx(i):breakIdx(i);
        xCell{i} = x(idx);
        yCell{i} = y(idx);
    end
end

function curve = preprocess_curves(gapList, xCell, yCell, levelMap)
    nGap = numel(gapList);
    curve = repmat(struct( ...
        'gap', 0, 'x', [], 'yRaw', [], 'yBc', [], 'base', 0, ...
        'amp', 0, 'xPeak', 0, 'xCenter', [], 'yNorm', [], ...
        'xLeftMap', [], 'xRightMap', []), nGap, 1);
    for ig = 1:nGap
        x = xCell{ig}(:);
        y = yCell{ig}(:);
        [base, amp, yNorm] = normalize_waveform(y);
        yBc = y - base;
        [~, iPeak] = max(yBc);
        xPeak = x(iPeak);
        xCenter = x - xPeak;
        [xLeftMap, xRightMap] = level_crossings_normalized(xCenter, yNorm, levelMap);

        curve(ig).gap = gapList(ig);
        curve(ig).x = x;
        curve(ig).yRaw = y;
        curve(ig).yBc = yBc;
        curve(ig).base = base;
        curve(ig).amp = amp;
        curve(ig).xPeak = xPeak;
        curve(ig).xCenter = xCenter;
        curve(ig).yNorm = yNorm;
        curve(ig).xLeftMap = xLeftMap;
        curve(ig).xRightMap = xRightMap;
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function [xLeft, xRight] = level_crossings_normalized(x, yNorm, levels)
    x = x(:);
    yNorm = yNorm(:);
    levels = levels(:);
    [~, iPeak] = max(yNorm);
    xRise = x(1:iPeak);
    yRise = yNorm(1:iPeak);
    xFall = x(iPeak:end);
    yFall = yNorm(iPeak:end);
    [yRiseU, idxRiseU] = unique(yRise, 'stable');
    xRiseU = xRise(idxRiseU);
    [yFallU, idxFallU] = unique(yFall, 'stable');
    xFallU = xFall(idxFallU);
    xLeft = nan(size(levels));
    xRight = nan(size(levels));
    for k = 1:numel(levels)
        lv = levels(k);
        if lv >= min(yRiseU) && lv <= max(yRiseU) && lv >= min(yFallU) && lv <= max(yFallU)
            xLeft(k) = interp1(yRiseU, xRiseU, lv, 'linear');
            xRight(k) = interp1(flipud(yFallU), flipud(xFallU), lv, 'linear');
        end
    end
end

function uGrid = make_reference_grid(refCurve, levelMap, nU)
    valid = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
    uMin = refCurve.xLeftMap(find(valid, 1, 'first'));
    uMax = refCurve.xRightMap(find(valid, 1, 'first'));
    uGrid = linspace(uMin, uMax, nU)';
end

function Cmat = register_all_curves(curve, refIdx, uGrid, interpMethod)
    nGap = numel(curve);
    Cmat = nan(nGap, numel(uGrid));
    refCurve = curve(refIdx);
    for ig = 1:nGap
        [~, yOnU] = map_waveform_to_reference_domain(curve(ig), refCurve, uGrid, interpMethod);
        Cmat(ig, :) = yOnU(:)';
    end
end

function [uOfX, yOnU] = map_waveform_to_reference_domain(curve, refCurve, uGrid, interpMethod)
    valid = isfinite(curve.xLeftMap) & isfinite(curve.xRightMap) & ...
            isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
    xLeftAnchor = [curve.xLeftMap(valid); 0];
    uLeftAnchor = [refCurve.xLeftMap(valid); 0];
    [xLeftAnchor, leftIdx] = unique(xLeftAnchor, 'stable');
    uLeftAnchor = uLeftAnchor(leftIdx);

    xRightAnchor = [0; flipud(curve.xRightMap(valid))];
    uRightAnchor = [0; flipud(refCurve.xRightMap(valid))];
    [xRightAnchor, rightIdx] = unique(xRightAnchor, 'stable');
    uRightAnchor = uRightAnchor(rightIdx);

    xCenter = curve.xCenter(:);
    yBc = curve.yBc(:);
    [~, iPeak] = max(yBc);

    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), interpMethod, 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), interpMethod, 'extrap');

    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];
    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    [uSorted, uniqIdx] = unique(uSorted, 'stable');
    ySorted = ySorted(uniqIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, interpMethod, 'extrap');
end

function txt = format_recommendations(T)
    lines = strings(height(T), 1);
    for i = 1:height(T)
        if isnan(T.min_subset_size(i))
            lines(i) = sprintf('RMSE <= %.1e pF : not reached', T.threshold_rmse_pF(i));
        else
            lines(i) = sprintf('RMSE <= %.1e pF : K >= %d, subset %s', ...
                T.threshold_rmse_pF(i), T.min_subset_size(i), T.best_subset_gaps_mm(i));
        end
    end
    txt = strjoin(cellstr(lines), newline);
end
