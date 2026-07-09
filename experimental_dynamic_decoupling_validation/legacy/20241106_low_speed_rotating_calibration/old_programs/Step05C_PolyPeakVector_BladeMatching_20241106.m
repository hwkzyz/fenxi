%% Step05C_PolyPeakVector_BladeMatching_20241106
% Match blade waveforms across sensors using polynomial-fitted peak vectors.
%
% Each pulse is fitted by a local polynomial. The fitted peak values over
% the selected revolutions form a vector for each sensor/blade candidate.
% Cross-sensor blade correspondence is then checked by peak-vector
% correlation, which is more discriminative than a single-pulse shape match.
%
% Matching uses the fluctuation pattern after per-vector normalization,
% rather than the raw peak level, because sensor static gain offsets are
% large while the informative part is the lap-to-lap modulation.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
legacyDir = fullfile(validationRoot, '20241106_low_speed_gap_prior_decoupling', 'legacy');
addpath(legacyDir);

%% Parameters
P.targetBlade = 6;
P.analysisSensors = [2 3 5 7];
P.analysisStartTime = 75.0;
P.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.resultSuffix = sprintf('Direct_%s', make_time_label_local(P.analysisStartTime));
P.polyDegree = 6;
P.maxSelectedLaps = inf;
P.edgeGuardRatio = 0.10;
P.usePeakX = true;
P.featureWeights = struct('peak_v', 0.7, 'peak_x', 0.3);
P.saveFigures = true;

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
resultFile = fullfile(routeDir, 'output', 'identification', ...
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.resultSuffix));

legacyCfg = Get_20241106_BTT_Config();
legacyCfg.dynamic_data_dir = fullfile(legacyCfg.dataset_root, legacyCfg.dynamic_cases{1});
dynamicCaseOutputDir = fullfile(legacyCfg.output_root, legacyCfg.dynamic_cases{1});
sensorConfigFile = fullfile(legacyCfg.reference_output_dir, 'Sensor_Config_20241106.mat');
figureDir = fullfile(routeDir, 'output', 'figures');
if exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

fprintf('\n=== Step05C: polynomial peak-vector blade matching ===\n');
fprintf('Target blade: B%d, sensors: %s, polynomial degree: %d\n', ...
    P.targetBlade, mat2str(P.analysisSensors), P.polyDegree);

loadedResult = load(resultFile, 'Result');
Result = loadedResult.Result;
loadedDynamic = load(Result.DynamicMapFile, 'DynamicMap');
DynamicMap = loadedDynamic.DynamicMap;
loadedSensorConfig = load(sensorConfigFile, 'Sensor_Config');
Sensor_Config = loadedSensorConfig.Sensor_Config;

oprTimes = load_opr_center_times_local(dynamicCaseOutputDir);
oprReference = build_opr_reference_from_case_local(dynamicCaseOutputDir, legacyCfg);
F_omega_deg = build_phase_speed_local_from_case_local(dynamicCaseOutputDir, legacyCfg);
probeCache = load_probe_cache_local(dynamicCaseOutputDir, P.analysisSensors, ...
    oprTimes, legacyCfg.opr_pulses_per_rev);
selectedRevIds = get_anchor_revolutions_local(DynamicMap, probeCache, P.maxSelectedLaps);
rawStream = load_case_raw_streams_window_local(legacyCfg, P.analysisSensors, DynamicMap.GlobalTimeWindow);

PeakData = build_poly_peak_vectors_local(rawStream, probeCache, selectedRevIds, ...
    Sensor_Config, oprReference, F_omega_deg, legacyCfg, P);
[matchTables, matrixStruct] = build_anchor_matching_local(PeakData, DynamicMap.SelectionInfo.anchor_sensor_id, P);

for i = 1:numel(matchTables)
    fprintf('\nReference CH%d -> CH%d\n', matchTables(i).referenceSensor, matchTables(i).testSensor);
    disp(matchTables(i).Table);
end

plot_peak_vectors_local(PeakData, P, figureDir);
plot_anchor_match_matrices_local(matrixStruct, P, figureDir);

allMatchRows = vertcat(matchTables.Table);
writetable(allMatchRows, fullfile(figureDir, ...
    sprintf('Step05C_PolyPeakVector_Matching_B%d_%s_20241106.csv', P.targetBlade, sensorTag)));
save(fullfile(figureDir, sprintf('Step05C_PolyPeakVector_Data_B%d_%s_20241106.mat', ...
    P.targetBlade, sensorTag)), 'PeakData', 'matrixStruct', 'allMatchRows', 'P');

function PeakData = build_poly_peak_vectors_local(rawStream, probeCache, selectedRevIds, ...
        Sensor_Config, oprReference, F_omega_deg, cfg, P)
bladeCount = cfg.blades_num;
sensorIds = P.analysisSensors(:).';
PeakData = repmat(struct( ...
    'sensor_id', [], ...
    'blade_id', [], ...
    'revolution_ids', [], ...
        'peak_v', [], ...
        'peak_x_mm', [], ...
        'peak_v_norm', [], ...
        'peak_x_norm', [], ...
        'fit_ok', []), numel(sensorIds), bladeCount);

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    probe = probeCache([probeCache.sensor_id] == sid);
    for bladeId = 1:bladeCount
        thetaStd = read_opr_center_standard_angle_local(Sensor_Config, sid, bladeId, oprReference);
        selectedRows = project_rows_by_revolution_local(probe, selectedRevIds, thetaStd);
        n = numel(selectedRows);
        peakV = nan(n, 1);
        peakX = nan(n, 1);
        fitOk = false(n, 1);

        for k = 1:n
            rowId = selectedRows(k);
            tPeak = probe.jilublade(rowId, 3);
            [tStart, tEnd] = build_dynamic_segment_window_local(probe.jilublade, rowId, tPeak);
            keep = rawStream(sid).T >= tStart & rawStream(sid).T <= tEnd;
            tSeg = rawStream(sid).T(keep);
            vSeg = rawStream(sid).V(keep);
            if numel(tSeg) < max(8, P.polyDegree + 2)
                continue;
            end
            idxPrev = find(probe.opr_times < tPeak, 1, 'last');
            if isempty(idxPrev)
                continue;
            end
            xAbs = map_segment_to_relative_xabs_local(probe.opr_times(idxPrev), ...
                tSeg, F_omega_deg, thetaStd, cfg.r_tip_mm);
            [peakV(k), peakX(k), fitOk(k)] = fit_polynomial_peak_local(xAbs, vSeg, P.polyDegree);
        end

        PeakData(is, bladeId).sensor_id = sid;
        PeakData(is, bladeId).blade_id = bladeId;
        PeakData(is, bladeId).revolution_ids = selectedRevIds(1:n);
        PeakData(is, bladeId).peak_v = peakV;
        PeakData(is, bladeId).peak_x_mm = peakX;
        PeakData(is, bladeId).peak_v_norm = normalize_vector_local(peakV, fitOk);
        PeakData(is, bladeId).peak_x_norm = normalize_vector_local(peakX, fitOk);
        PeakData(is, bladeId).fit_ok = fitOk;
    end
end
end

function [matchTables, matrixStruct] = build_anchor_matching_local(PeakData, anchorSensor, P)
sensorIds = unique([PeakData(:, 1).sensor_id], 'stable');
anchorIdx = find(sensorIds == anchorSensor, 1, 'first');
if isempty(anchorIdx)
    anchorIdx = 1;
    anchorSensor = sensorIds(anchorIdx);
end
bladeCount = size(PeakData, 2);
targetBlade = P.targetBlade;
otherSensors = sensorIds(sensorIds ~= anchorSensor);
matchTables = repmat(struct('referenceSensor', [], 'testSensor', [], 'Table', []), numel(otherSensors), 1);
matrixStruct = repmat(struct('referenceSensor', [], 'testSensor', [], 'Corr', [], 'bladeIds', []), numel(otherSensors), 1);

for io = 1:numel(otherSensors)
    testSensor = otherSensors(io);
    testIdx = find(sensorIds == testSensor, 1, 'first');
    Corr = nan(bladeCount, bladeCount);
    for ba = 1:bladeCount
        for bb = 1:bladeCount
            corrV = compute_single_feature_corr_local( ...
                PeakData(anchorIdx, ba), PeakData(testIdx, bb), 'peak_v_norm');
            corrX = compute_single_feature_corr_local( ...
                PeakData(anchorIdx, ba), PeakData(testIdx, bb), 'peak_x_norm');
            Corr(ba, bb) = combine_feature_score_local(corrV, corrX, P.featureWeights, P.usePeakX);
        end
    end

    rows = repmat(struct( ...
        'referenceSensor', anchorSensor, ...
        'testSensor', testSensor, ...
        'referenceBlade', nan, ...
        'bestMatchedBlade', nan, ...
            'bestScore', nan, ...
            'diagonalScore', nan, ...
            'bestCorrPeakV', nan, ...
            'diagonalCorrPeakV', nan, ...
            'bestCorrPeakX', nan, ...
            'diagonalCorrPeakX', nan, ...
            'marginToSecond', nan, ...
            'isDiagonalBest', false, ...
            'isTargetBlade', false), bladeCount, 1);
    for ba = 1:bladeCount
        rowCorr = Corr(ba, :);
        [sortedCorr, order] = sort(rowCorr, 'descend', 'MissingPlacement', 'last');
        bestBlade = order(1);
        rows(ba).referenceBlade = ba;
        rows(ba).bestMatchedBlade = bestBlade;
        rows(ba).bestScore = sortedCorr(1);
        rows(ba).diagonalScore = rowCorr(ba);
        rows(ba).bestCorrPeakV = compute_single_feature_corr_local( ...
            PeakData(anchorIdx, ba), PeakData(testIdx, bestBlade), 'peak_v_norm');
        rows(ba).diagonalCorrPeakV = compute_single_feature_corr_local( ...
            PeakData(anchorIdx, ba), PeakData(testIdx, ba), 'peak_v_norm');
        rows(ba).bestCorrPeakX = compute_single_feature_corr_local( ...
            PeakData(anchorIdx, ba), PeakData(testIdx, bestBlade), 'peak_x_norm');
        rows(ba).diagonalCorrPeakX = compute_single_feature_corr_local( ...
            PeakData(anchorIdx, ba), PeakData(testIdx, ba), 'peak_x_norm');
        rows(ba).marginToSecond = sortedCorr(1) - sortedCorr(min(2, numel(sortedCorr)));
        rows(ba).isDiagonalBest = (bestBlade == ba);
        rows(ba).isTargetBlade = (ba == targetBlade);
    end
    matchTables(io).referenceSensor = anchorSensor;
    matchTables(io).testSensor = testSensor;
    matchTables(io).Table = struct2table(rows);
    matrixStruct(io).referenceSensor = anchorSensor;
    matrixStruct(io).testSensor = testSensor;
    matrixStruct(io).Corr = Corr;
    matrixStruct(io).bladeIds = 1:bladeCount;
end
end

function plot_peak_vectors_local(PeakData, P, figureDir)
sensorIds = unique([PeakData(:, 1).sensor_id], 'stable');
bladeCount = size(PeakData, 2);
fig = figure('Name', sprintf('20241106 polynomial peak vectors B%d', P.targetBlade), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 22, 16]);
tiledlayout(fig, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    nexttile;
    hold on;
    cmap = lines(bladeCount);
    for bladeId = 1:bladeCount
        v = PeakData(is, bladeId).peak_v_norm;
        ok = PeakData(is, bladeId).fit_ok;
        lap = 1:numel(v);
        if bladeId == P.targetBlade
            plot(lap(ok), v(ok), '-o', 'Color', cmap(bladeId, :), 'LineWidth', 1.6, ...
                'MarkerSize', 3.5, 'DisplayName', sprintf('B%d target', bladeId));
        else
            plot(lap(ok), v(ok), '-', 'Color', cmap(bladeId, :), 'LineWidth', 0.9, ...
                'DisplayName', sprintf('B%d', bladeId));
        end
    end
    ylabel(sprintf('CH%d norm. peak', sensorIds(is)));
    if is == 1
        legend('Location', 'eastoutside', 'Box', 'off');
    end
    if is == numel(sensorIds)
        xlabel('Selected revolution index');
    end
    title(sprintf('Normalized polynomial peak vector, CH%d', sensorIds(is)), 'FontWeight', 'normal');
    style_axes_local();
end

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05C_PolyPeakVectors_B%d_S%s_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors))), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05C_PolyPeakVectors_B%d_S%s_20241106.pdf', ...
        P.targetBlade, sprintf('%d', P.analysisSensors))), 'ContentType', 'vector');
end
end

function plot_anchor_match_matrices_local(matrixStruct, P, figureDir)
fig = figure('Name', sprintf('20241106 polynomial peak-vector matching B%d', P.targetBlade), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 5.5 * max(1, numel(matrixStruct))]);
tiledlayout(fig, 1, numel(matrixStruct), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(matrixStruct)
    nexttile;
    imagesc(matrixStruct(i).Corr, [-1, 1]);
    axis image;
    colormap(gca, parula);
    colorbar;
    xlabel(sprintf('CH%d blade', matrixStruct(i).testSensor));
    ylabel(sprintf('CH%d blade', matrixStruct(i).referenceSensor));
    title(sprintf('CH%d -> CH%d', matrixStruct(i).referenceSensor, matrixStruct(i).testSensor), ...
        'FontWeight', 'normal');
    xticks(1:numel(matrixStruct(i).bladeIds));
    yticks(1:numel(matrixStruct(i).bladeIds));
    hold on;
    for b = 1:numel(matrixStruct(i).bladeIds)
        plot(b, b, 'wo', 'MarkerSize', 7, 'LineWidth', 1.2);
    end
    plot(P.targetBlade, P.targetBlade, 'rs', 'MarkerSize', 9, 'LineWidth', 1.4);
    style_axes_local();
end

if P.saveFigures
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05C_PolyPeakMatchMatrix_B%d_S%s_20241106.png', ...
        P.targetBlade, sprintf('%d', P.analysisSensors))), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, ...
        sprintf('Step05C_PolyPeakMatchMatrix_B%d_S%s_20241106.pdf', ...
        P.targetBlade, sprintf('%d', P.analysisSensors))), 'ContentType', 'vector');
end
end

function [peakV, peakX, ok] = fit_polynomial_peak_local(x, v, degree)
peakV = nan;
peakX = nan;
ok = false;
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
if numel(x) < degree + 2 || range(x) <= eps
    return;
end
[x, order] = sort(x);
v = v(order);
xCenter = mean(x);
xScale = max(std(x), eps);
xn = (x - xCenter) / xScale;
coef = polyfit(xn, v, degree);
xDense = linspace(min(x), max(x), 500).';
vFit = polyval(coef, (xDense - xCenter) / xScale);
[peakV, idx] = max(vFit);
peakX = xDense(idx);
rangeX = max(xDense) - min(xDense);
leftGuard = min(xDense) + 0.10 * rangeX;
rightGuard = max(xDense) - 0.10 * rangeX;
ok = isfinite(peakV) && isfinite(peakX) && peakX >= leftGuard && peakX <= rightGuard;
if ~ok
    peakV = nan;
    peakX = nan;
end
end

function rawStream = load_case_raw_streams_window_local(cfg, sensorIds, timeWindow)
fileRanges = build_dynamic_file_ranges_local(cfg);
selectMask = [fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2);
fileRanges = fileRanges(selectMask);
maxSid = max(sensorIds);
rawStream(maxSid) = struct('T', [], 'V', []);
for i = 1:numel(fileRanges)
    fileId = fileRanges(i).file_id;
    offset = fileRanges(i).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_case_channel_local(cfg.dynamic_data_dir, sid, fileId, cfg.pinlv);
        tGlobal = tLocal(:) + offset;
        keep = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        rawStream(sid).T = [rawStream(sid).T; tGlobal(keep)]; %#ok<AGROW>
        rawStream(sid).V = [rawStream(sid).V; vLocal(keep)]; %#ok<AGROW>
    end
end
end

function fileRanges = build_dynamic_file_ranges_local(cfg)
fileIds = list_case_file_ids_local(cfg.dynamic_data_dir, cfg.opr_id);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = [];
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_case_channel_local(cfg.dynamic_data_dir, cfg.opr_id, fileIds(i), cfg.pinlv);
    if isempty(lastEnd)
        offset = 0;
    else
        offset = lastEnd + 1 / cfg.pinlv - tOpr(1);
    end
    lastEnd = tOpr(end) + offset;
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
end
end

function [tSec, v] = load_raw_case_channel_local(caseDir, sid, fileId, pinlv)
filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', sid, fileId));
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function fileIds = list_case_file_ids_local(caseDir, channelId)
d = dir(fullfile(caseDir, sprintf('4-%d-*.mat', channelId)));
fileIds = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channelId) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        fileIds(i) = str2double(tok{1});
    end
end
fileIds = sort(unique(fileIds(~isnan(fileIds))));
end

function probeCache = load_probe_cache_local(caseOutputDir, sensorIds, oprTimes, oprPulsesPerRev)
probeCache = repmat(struct( ...
    'sensor_id', NaN, ...
    'jilublade', [], ...
    'rel_angles_deg', [], ...
    'revolution_index', [], ...
    'opr_times', []), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    loadedProbe = load(fullfile(caseOutputDir, sprintf('jilublade_probe%d.mat', sid)), 'jilublade');
    jilublade = loadedProbe.jilublade;
    [relAnglesDeg, revolutionIndex] = compute_pulse_reference_geometry_local( ...
        jilublade(:, 3), oprTimes, oprPulsesPerRev);
    probeCache(is).sensor_id = sid;
    probeCache(is).jilublade = jilublade;
    probeCache(is).rel_angles_deg = relAnglesDeg;
    probeCache(is).revolution_index = revolutionIndex;
    probeCache(is).opr_times = oprTimes;
end
end

function selectedRevIds = get_anchor_revolutions_local(DynamicMap, probeCache, maxSelectedLaps)
anchorSid = DynamicMap.SelectionInfo.anchor_sensor_id;
anchorRows = DynamicMap.SelectionInfo.anchor_selected_rows(:);
if isfinite(maxSelectedLaps)
    anchorRows = anchorRows(1:min(numel(anchorRows), maxSelectedLaps));
end
anchorProbe = probeCache([probeCache.sensor_id] == anchorSid);
selectedRevIds = anchorProbe.revolution_index(anchorRows);
selectedRevIds = selectedRevIds(isfinite(selectedRevIds));
end

function selectedRows = project_rows_by_revolution_local(probe, selectedRevIds, thetaStd)
selectedRows = nan(numel(selectedRevIds), 1);
for k = 1:numel(selectedRevIds)
    candidates = find(probe.revolution_index == selectedRevIds(k) & isfinite(probe.rel_angles_deg));
    if isempty(candidates)
        continue;
    end
    candidateError = wrap_to_signed_period_local(probe.rel_angles_deg(candidates) - thetaStd, 360);
    [~, bestPos] = min(abs(candidateError));
    selectedRows(k) = candidates(bestPos);
end
selectedRows = selectedRows(isfinite(selectedRows));
end

function [relAnglesDeg, revolutionIndex] = compute_pulse_reference_geometry_local(arrivalTimes, oprTimes, oprPulsesPerRev)
relAnglesDeg = nan(size(arrivalTimes));
revolutionIndex = nan(size(arrivalTimes));
for i = 1:numel(arrivalTimes)
    idx = find(oprTimes < arrivalTimes(i), 1, 'last');
    if isempty(idx) || idx + oprPulsesPerRev > numel(oprTimes)
        continue;
    end
    dtRev = oprTimes(idx + oprPulsesPerRev) - oprTimes(idx);
    if ~isfinite(dtRev) || dtRev <= eps
        continue;
    end
    relAnglesDeg(i) = 360 * (arrivalTimes(i) - oprTimes(idx)) / dtRev;
    revolutionIndex(i) = idx;
end
end

function [tStart, tEnd] = build_dynamic_segment_window_local(jilublade, rowId, tPeak)
syncCfg = build_single_sync_experiment_config_20241106();
if strcmpi(syncCfg.dynamic_window_mode, 'legacy_row_bounds')
    tStart = jilublade(rowId, 1) - syncCfg.pulse_pad_sec;
    tEnd = jilublade(rowId, 2) + syncCfg.pulse_pad_sec;
else
    tStart = tPeak - syncCfg.pulse_window_sec;
    tEnd = tPeak + syncCfg.pulse_window_sec;
end
end

function xAbs = map_segment_to_relative_xabs_local(tRef, tSeg, F_omega_deg, thetaStd, rTipMm)
thetaPoints = map_segment_to_relative_angle_local(tRef, tSeg, F_omega_deg);
thetaDiff = wrap_to_signed_period_local(thetaPoints - thetaStd, 360);
xAbs = deg2rad(thetaDiff) * rTipMm;
end

function thetaPoints = map_segment_to_relative_angle_local(tRef, tSeg, F_omega_deg)
dtFirst = linspace(tRef, tSeg(1), 10);
thetaBase = trapz(dtFirst, F_omega_deg(dtFirst));
wSeg = F_omega_deg(tSeg);
thetaRel = cumtrapz(tSeg, wSeg);
thetaPoints = thetaBase + thetaRel;
end

function F = build_phase_speed_local_from_case_local(caseOutputDir, cfg)
loadedOmega = load(fullfile(caseOutputDir, 'omega.mat'));
if isfield(loadedOmega, 'omega') && ~isempty(loadedOmega.omega)
    omega = loadedOmega.omega;
    F = griddedInterpolant(omega(:, 1), omega(:, 2), 'linear', 'nearest');
    return;
end
oprTimes = load_opr_center_times_local(caseOutputDir);
spdT = oprTimes(1:(end - cfg.opr_pulses_per_rev));
spdV = 360 ./ max(oprTimes((cfg.opr_pulses_per_rev + 1):end) - oprTimes(1:(end - cfg.opr_pulses_per_rev)), eps);
F = griddedInterpolant(spdT, spdV, 'linear', 'nearest');
end

function oprReference = build_opr_reference_from_case_local(caseOutputDir, cfg)
loaded = load(fullfile(caseOutputDir, 'jiluOPR.mat'), 'jiluOPR');
jiluOPR = loaded.jiluOPR;
oprReference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0);
if size(jiluOPR, 2) < 2
    return;
end
centerTime = jiluOPR(:, 1);
startTime = jiluOPR(:, 2);
n = min(numel(centerTime) - cfg.opr_pulses_per_rev, numel(startTime));
if n < 1
    return;
end
dtCenter = centerTime(1:n) - startTime(1:n);
dtRev = centerTime((1:n) + cfg.opr_pulses_per_rev) - centerTime(1:n);
valid = isfinite(dtCenter) & isfinite(dtRev) & dtRev > eps;
shiftDeg = 360 * dtCenter(valid) ./ dtRev(valid);
oprReference.phase_shift_deg = median(shiftDeg, 'omitnan');
oprReference.phase_shift_mm = oprReference.phase_shift_deg * (pi / 180) * cfg.r_tip_mm;
end

function oprTimes = load_opr_center_times_local(caseOutputDir)
loaded = load(fullfile(caseOutputDir, 'jiluOPR.mat'), 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function thetaStd = read_opr_center_standard_angle_local(Sensor_Config, sid, bladeId, oprReference)
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    thetaStd = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, bladeId);
    return;
end
thetaStd = Sensor_Config.Standard_Relative_Angles(sid, bladeId);
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
thetaStd = thetaStd - oprReference.phase_shift_deg;
end

function value = local_corrcoef_scalar(a, b)
if numel(a) < 3 || numel(b) < 3 || std(a, 'omitnan') <= eps || std(b, 'omitnan') <= eps
    value = nan;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function corrValue = compute_single_feature_corr_local(A, B, fieldName)
va = A.(fieldName);
vb = B.(fieldName);
oka = A.fit_ok;
okb = B.fit_ok;
n = min(numel(va), numel(vb));
valid = oka(1:n) & okb(1:n) & isfinite(va(1:n)) & isfinite(vb(1:n));
corrValue = local_corrcoef_scalar(va(valid), vb(valid));
end

function score = combine_feature_score_local(corrV, corrX, weights, usePeakX)
parts = [];
partWeights = [];
if isfinite(corrV)
    parts(end + 1, 1) = corrV; %#ok<AGROW>
    partWeights(end + 1, 1) = weights.peak_v; %#ok<AGROW>
end
if usePeakX && isfinite(corrX)
    parts(end + 1, 1) = corrX; %#ok<AGROW>
    partWeights(end + 1, 1) = weights.peak_x; %#ok<AGROW>
end
if isempty(parts)
    score = nan;
    return;
end
score = sum(parts .* partWeights) / sum(partWeights);
end

function vNorm = normalize_vector_local(v, fitOk)
v = v(:);
fitOk = fitOk(:);
vNorm = nan(size(v));
valid = fitOk & isfinite(v);
if nnz(valid) < 3
    return;
end
mu = mean(v(valid), 'omitnan');
sigma = std(v(valid), 'omitnan');
if ~isfinite(sigma) || sigma <= eps
    vNorm(valid) = 0;
    return;
end
vNorm(valid) = (v(valid) - mu) / sigma;
end

function angle = wrap_to_signed_period_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function style_axes_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function label = make_time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end
