%% Step05E_Visualize_SixBlade_LowHigh_Sequence_20241106
% Compare one continuous 6-blade low-speed sequence against one continuous
% high-speed revolution, using polynomial-fitted peak values for each pulse.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
legacyDir = fullfile(validationRoot, '20241106_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacyDir, 'dir') ~= 7
    error('Legacy helper folder not found: %s', legacyDir);
end
addpath(legacyDir);

%% Parameters
P.targetBlade = 6;
P.analysisSensors = [2 3 5 7];
P.analysisStartTime = 75.0;
P.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', make_time_label_local(P.analysisStartTime), P.templateSuffix);
P.windowID = 1;
P.revolutionInWindow = 1;
P.polyDegree = 6;
P.saveFigures = true;
P.savePdf = false;

if use_env_overrides_local()
    targetBladeEnv = str2double(strtrim(getenv('STEP20241106_TARGET_BLADE')));
    if isfinite(targetBladeEnv) && targetBladeEnv >= 1
        P.targetBlade = round(targetBladeEnv);
    end
    analysisSensorsEnv = sscanf(strtrim(getenv('STEP20241106_ANALYSIS_SENSORS')), '%d').';
    if ~isempty(analysisSensorsEnv)
        P.analysisSensors = unique(analysisSensorsEnv, 'stable');
    end
    analysisStartEnv = str2double(strtrim(getenv('STEP20241106_ANALYSIS_START_TIME')));
    if isfinite(analysisStartEnv)
        P.analysisStartTime = analysisStartEnv;
        P.dynamicSuffix = sprintf('Main20L_W3S1_%s_%s', make_time_label_local(P.analysisStartTime), P.templateSuffix);
    end
    windowIdEnv = str2double(strtrim(getenv('STEP20241106_WINDOW_ID')));
    if isfinite(windowIdEnv) && windowIdEnv >= 1
        P.windowID = round(windowIdEnv);
    end
    revInWindowEnv = str2double(strtrim(getenv('STEP20241106_REVOLUTION_IN_WINDOW')));
    if isfinite(revInWindowEnv) && revInWindowEnv >= 1
        P.revolutionInWindow = round(revInWindowEnv);
    end
end

sensorTag = ['S', sprintf('%d', P.analysisSensors)];
dynamicMapFile = fullfile(routeDir, 'output', 'dynamic_maps', ...
    sprintf('DynamicMap_B%d_%s_SlidingWindows_%s_20241106.mat', ...
    P.targetBlade, sensorTag, P.dynamicSuffix));
numberingFile = fullfile(routeDir, 'output', 'blade_numbering', ...
    sprintf('BladeNumbering_%s_%s_20241106.mat', sensorTag, P.dynamicSuffix));
figureDir = fullfile(routeDir, 'output', 'figures');
if exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

fprintf('\n=== Step05E: six-blade low/high sequence audit ===\n');
fprintf('Target blade: B%d, sensors: %s\n', P.targetBlade, mat2str(P.analysisSensors));
fprintf('Dynamic map: %s\n', dynamicMapFile);
fprintf('Blade numbering: %s\n', numberingFile);

cfg = Get_20241106_BTT_Config();
DynamicMap = load_required_variable_local(dynamicMapFile, 'DynamicMap');
BladeNumbering = load_required_variable_local(numberingFile, 'BladeNumbering');
loadedSensorConfig = load(fullfile(cfg.reference_output_dir, 'Sensor_Config_20241106.mat'), 'Sensor_Config');
Sensor_Config = loadedSensorConfig.Sensor_Config;
lowNumberingSummaryFile = fullfile(cfg.reference_output_dir, 'Sensor_Config_Summary_20241106.csv');
if exist(lowNumberingSummaryFile, 'file') == 2
    lowNumberingSummary = readtable(lowNumberingSummaryFile);
else
    lowNumberingSummary = table();
end

caseOutputDir = fullfile(cfg.output_root, cfg.dynamic_cases{1});
oprTimesHigh = load_opr_center_times_local(caseOutputDir);
probeCacheHigh = load_probe_cache_local(caseOutputDir, P.analysisSensors, oprTimesHigh, cfg.opr_pulses_per_rev);
selectedRevIds = get_window_anchor_revolutions_local(DynamicMap, probeCacheHigh, P.windowID);
if isempty(selectedRevIds)
    error('No anchor revolutions found for window %d.', P.windowID);
end
selectedRevId = selectedRevIds(P.revolutionInWindow);

highCaseDir = fullfile(cfg.dataset_root, cfg.dynamic_cases{1});
lowCaseDir = fullfile(cfg.dataset_root, cfg.low_speed_case);
requiredLowPulseCount = max(arrayfun(@(sid) Sensor_Config.Target_Indices(sid), P.analysisSensors)) + cfg.blades_num;
lowCase = extract_reference_prefix_local(lowCaseDir, cfg, P.analysisSensors, requiredLowPulseCount);
highFileRanges = build_case_file_ranges_local(highCaseDir, cfg.pinlv);

bladeColors = lines(cfg.blades_num);
sequenceData = repmat(struct( ...
    'sensor_id', NaN, ...
    'low', struct(), ...
    'high', struct(), ...
    'local_to_physical', [], ...
    'physical_to_local', [], ...
    'corr_norm_peak_phys', NaN, ...
    'low_start_index', NaN, ...
    'low_fingerprint_corr', NaN), numel(P.analysisSensors), 1);

for i = 1:numel(P.analysisSensors)
    sid = P.analysisSensors(i);
    localToPhysical = BladeNumbering.local_to_physical_blade_ids(:, i).';
    physicalToLocal = BladeNumbering.physical_to_local_blade_ids(:, i).';
    lowSeq = build_low_speed_sequence_local(lowCase, lowCaseDir, sid, Sensor_Config, ...
        lowNumberingSummary, cfg, P.polyDegree);
    highSeq = build_high_speed_sequence_local(highCaseDir, highFileRanges, probeCacheHigh, sid, selectedRevId, ...
        localToPhysical, cfg, P.polyDegree);

    lowPeakNorm = lowSeq.peak_values_phys / max(lowSeq.peak_values_phys);
    highPeakNorm = highSeq.peak_values_phys / max(highSeq.peak_values_phys);
    corrNormPeakPhys = local_corrcoef_scalar_local(lowPeakNorm(:), highPeakNorm(:));

    sequenceData(i).sensor_id = sid;
    sequenceData(i).low = lowSeq;
    sequenceData(i).high = highSeq;
    sequenceData(i).local_to_physical = localToPhysical;
    sequenceData(i).physical_to_local = physicalToLocal;
    sequenceData(i).corr_norm_peak_phys = corrNormPeakPhys;
    sequenceData(i).low_start_index = lowSeq.start_index;
    sequenceData(i).low_fingerprint_corr = lowSeq.fingerprint_corr;
end

summaryTable = build_sequence_summary_table_local(sequenceData);
disp(summaryTable);
writetable(summaryTable, fullfile(figureDir, ...
    sprintf('Step05E_SixBladeSequenceSummary_B%d_%s_W%02d_R%02d_20241106.csv', ...
    P.targetBlade, sensorTag, P.windowID, P.revolutionInWindow)));

plot_six_blade_low_high_sequence_local(sequenceData, bladeColors, P, selectedRevId, figureDir);

function S = load_required_variable_local(filepath, varName)
if exist(filepath, 'file') ~= 2
    error('Required file not found: %s', filepath);
end
loaded = load(filepath, varName);
if ~isfield(loaded, varName)
    error('Variable %s is missing in %s.', varName, filepath);
end
S = loaded.(varName);
end

function tf = use_env_overrides_local()
text = lower(strtrim(getenv('STEP20241106_USE_ENV_OVERRIDES')));
tf = ismember(text, {'1', 'true', 'yes', 'y', 'on'});
end

function lowSeq = build_low_speed_sequence_local(caseData, caseDir, sid, Sensor_Config, lowNumberingSummary, cfg, polyDegree)
startIdx = Sensor_Config.Target_Indices(sid);
pulseIdx = startIdx + (0:(cfg.blades_num - 1));
ch = caseData.channels(sid);
if pulseIdx(end) > numel(ch.arrival_times)
    error('CH%d low-speed pulse indices exceed extracted pulse count.', sid);
end

pad = 2e-4;
timeWindow = [ch.start_times(pulseIdx(1)) - pad, ch.end_times(pulseIdx(end)) + pad];
raw = load_case_raw_stream_for_window_local(caseDir, caseData.file_ranges, sid, cfg.pinlv, timeWindow);

segments = repmat(struct( ...
    'physical_blade', NaN, ...
    't_seg', [], ...
    'v_seg', [], ...
    't_fit', [], ...
    'v_fit', [], ...
    'peak_t', NaN, ...
    'peak_v', NaN), cfg.blades_num, 1);
peakValues = nan(1, cfg.blades_num);

for b = 1:cfg.blades_num
    idx = pulseIdx(b);
    seg = build_polynomial_segment_local(raw, ch.start_times(idx), ch.end_times(idx), polyDegree);
    seg.physical_blade = b;
    segments(b) = seg;
    peakValues(b) = seg.peak_v;
end

lowSeq = struct();
lowSeq.time_window = timeWindow;
lowSeq.raw = raw;
lowSeq.segments = segments;
lowSeq.peak_values_phys = peakValues;
lowSeq.start_index = startIdx;
lowSeq.fingerprint_corr = lookup_low_fingerprint_corr_local(lowNumberingSummary, sid);
end

function caseData = extract_reference_prefix_local(caseDir, cfg, sensorIds, requiredPulseCount)
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR files found in %s.', caseDir);
end

maxChannelId = max([sensorIds(:).', cfg.opr_id]);
channels(maxChannelId) = struct('arrival_times', [], 'peak_features', [], 'start_times', [], 'end_times', []);
tailData(maxChannelId) = struct('t', [], 'v', []);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEndTime = [];
usedFileCount = 0;

for iFile = 1:numel(fileIds)
    fileId = fileIds(iFile);
    [tOprLocal, ~] = load_raw_channel_local(caseDir, cfg.opr_id, fileId, cfg.pinlv);
    if isempty(lastEndTime)
        offset = 0;
    else
        offset = lastEndTime + 1 / cfg.pinlv - tOprLocal(1);
    end
    tOprGlobal = tOprLocal + offset;
    lastEndTime = tOprGlobal(end);

    if iFile == 1 && cfg.initial_trim_points > 0 && numel(tOprGlobal) > cfg.initial_trim_points
        tOprGlobal = tOprGlobal((cfg.initial_trim_points + 1):end);
    end

    usedFileCount = usedFileCount + 1;
    fileRanges(usedFileCount).file_id = fileId;
    fileRanges(usedFileCount).offset = offset;
    fileRanges(usedFileCount).t_start = tOprGlobal(1);
    fileRanges(usedFileCount).t_end = tOprGlobal(end);

    for sid = sensorIds
        filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', sid, fileId));
        if ~isfile(filepath)
            continue;
        end

        [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, cfg.pinlv);
        tGlobal = tLocal + offset;

        if iFile == 1 && cfg.initial_trim_points > 0 && numel(tGlobal) > cfg.initial_trim_points
            keepIdx = (cfg.initial_trim_points + 1):numel(tGlobal);
            tGlobal = tGlobal(keepIdx);
            vLocal = vLocal(keepIdx);
        end

        [pulseInfo, nextTail] = extract_sensor_pulses_local( ...
            tGlobal, vLocal, cfg.sensor_thresholds(sid), cfg.gap_points, tailData(sid));
        tailData(sid) = nextTail;
        if isempty(pulseInfo.arrival_times)
            continue;
        end

        channels(sid).arrival_times = [channels(sid).arrival_times; pulseInfo.arrival_times]; %#ok<AGROW>
        channels(sid).peak_features = [channels(sid).peak_features; pulseInfo.peak_features]; %#ok<AGROW>
        channels(sid).start_times = [channels(sid).start_times; pulseInfo.start_times]; %#ok<AGROW>
        channels(sid).end_times = [channels(sid).end_times; pulseInfo.end_times]; %#ok<AGROW>
    end

    enoughMask = true(numel(sensorIds), 1);
    for iSensor = 1:numel(sensorIds)
        sid = sensorIds(iSensor);
        enoughMask(iSensor) = numel(channels(sid).arrival_times) >= requiredPulseCount;
    end
    if all(enoughMask)
        break;
    end
end

caseData = struct();
caseData.case_dir = caseDir;
caseData.file_ranges = fileRanges(1:usedFileCount);
caseData.channels = channels;
end

function highSeq = build_high_speed_sequence_local(caseDir, fileRanges, probeCache, sid, revId, localToPhysical, cfg, polyDegree)
probe = probeCache([probeCache.sensor_id] == sid);
candidateRows = find(probe.revolution_index == revId);
if numel(candidateRows) ~= cfg.blades_num
    error('CH%d revolution %d has %d pulses, expected %d.', sid, revId, numel(candidateRows), cfg.blades_num);
end
[~, order] = sort(probe.jilublade(candidateRows, 3), 'ascend');
candidateRows = candidateRows(order);

pad = 2e-5;
timeWindow = [probe.jilublade(candidateRows(1), 1) - pad, probe.jilublade(candidateRows(end), 2) + pad];
raw = load_case_raw_stream_for_window_local(caseDir, fileRanges, sid, cfg.pinlv, timeWindow);

segments = repmat(struct( ...
    'local_slot', NaN, ...
    'physical_blade', NaN, ...
    'row_id', NaN, ...
    't_seg', [], ...
    'v_seg', [], ...
    't_fit', [], ...
    'v_fit', [], ...
    'peak_t', NaN, ...
    'peak_v', NaN), cfg.blades_num, 1);
peakValuesPhys = nan(1, cfg.blades_num);

for k = 1:cfg.blades_num
    rowId = candidateRows(k);
    physicalBlade = localToPhysical(k);
    seg = build_polynomial_segment_local(raw, probe.jilublade(rowId, 1), probe.jilublade(rowId, 2), polyDegree);
    seg.local_slot = k;
    seg.physical_blade = physicalBlade;
    seg.row_id = rowId;
    segments(k) = seg;
    peakValuesPhys(physicalBlade) = seg.peak_v;
end

highSeq = struct();
highSeq.time_window = timeWindow;
highSeq.raw = raw;
highSeq.segments = segments;
highSeq.peak_values_phys = peakValuesPhys;
end

function raw = load_case_raw_stream_for_window_local(caseDir, fileRanges, sid, pinlv, timeWindow)
if nargin < 2 || isempty(fileRanges)
    fileRanges = build_case_file_ranges_local(caseDir, pinlv);
end
selectMask = [fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2);
fileRanges = fileRanges(selectMask);
raw = struct('T', [], 'V', []);
for i = 1:numel(fileRanges)
    fileId = fileRanges(i).file_id;
    offset = fileRanges(i).offset;
    [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, pinlv);
    tGlobal = tLocal(:) + offset;
    keepMask = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
    raw.T = [raw.T; tGlobal(keepMask)]; %#ok<AGROW>
    raw.V = [raw.V; vLocal(keepMask)]; %#ok<AGROW>
end
end

function fileRanges = build_case_file_ranges_local(caseDir, pinlv)
fileIds = list_case_file_ids_local(caseDir, 1);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = [];
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_channel_local(caseDir, 1, fileIds(i), pinlv);
    if isempty(lastEnd)
        offset = 0;
    else
        offset = lastEnd + 1 / pinlv - tOpr(1);
    end
    lastEnd = tOpr(end) + offset;
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
end
end

function seg = build_polynomial_segment_local(raw, tStart, tEnd, polyDegree)
keep = raw.T >= tStart & raw.T <= tEnd;
tSeg = raw.T(keep);
vSeg = raw.V(keep);
if numel(tSeg) < max(8, polyDegree + 2)
    error('Too few points between %.9f and %.9f for polynomial fitting.', tStart, tEnd);
end
[tSeg, order] = sort(tSeg(:));
vSeg = vSeg(order);
[peakT, peakV, tFit, vFit] = fit_polynomial_peak_local(tSeg, vSeg, polyDegree);
seg = struct( ...
    't_seg', tSeg, ...
    'v_seg', vSeg, ...
    't_fit', tFit, ...
    'v_fit', vFit, ...
    'peak_t', peakT, ...
    'peak_v', peakV);
end

function [peakT, peakV, tFit, vFit] = fit_polynomial_peak_local(tSeg, vSeg, degree)
realDegree = min(degree, numel(tSeg) - 1);
[p, ~, mu] = polyfit(tSeg, vSeg, realDegree);
tFit = linspace(min(tSeg), max(tSeg), 200).';
vFit = polyval(p, tFit, [], mu);
[peakV, idx] = max(vFit);
peakT = tFit(idx);
end

function summaryTable = build_sequence_summary_table_local(sequenceData)
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'low_start_index', NaN, ...
    'low_fingerprint_corr', NaN, ...
    'corr_norm_peak_phys', NaN, ...
    'low_peak_phys', "", ...
    'high_peak_phys', "", ...
    'local_to_physical', ""), numel(sequenceData), 1);
for i = 1:numel(sequenceData)
    rows(i).sensor_id = sequenceData(i).sensor_id;
    rows(i).low_start_index = sequenceData(i).low_start_index;
    rows(i).low_fingerprint_corr = sequenceData(i).low_fingerprint_corr;
    rows(i).corr_norm_peak_phys = sequenceData(i).corr_norm_peak_phys;
    rows(i).low_peak_phys = string(mat2str(sequenceData(i).low.peak_values_phys, 4));
    rows(i).high_peak_phys = string(mat2str(sequenceData(i).high.peak_values_phys, 4));
    rows(i).local_to_physical = string(mat2str(sequenceData(i).local_to_physical));
end
summaryTable = struct2table(rows);
end

function plot_six_blade_low_high_sequence_local(sequenceData, bladeColors, P, selectedRevId, figureDir)
fig = figure('Name', sprintf('20241106 six-blade low/high sequence B%d W%02d R%02d', ...
    P.targetBlade, P.windowID, P.revolutionInWindow), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 30, 24]);
tiledlayout(fig, numel(sequenceData), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sequenceData)
    S = sequenceData(i);
    sid = S.sensor_id;

    nexttile;
    hold on;
    plot_sequence_panel_local(S.low.raw, S.low.segments, bladeColors, 'low');
    xlabel('Low-speed time in selected 6-pulse block (ms)');
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('CH%d low-speed Step1 fingerprint-numbered sequence: B1-B6', sid), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 8.5);
    subtitle(sprintf('start idx=%d, Step1 fp corr=%.3f', ...
        S.low_start_index, S.low_fingerprint_corr), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 6.8);
    style_axes_local();

    nexttile;
    hold on;
    plot_sequence_panel_local(S.high.raw, S.high.segments, bladeColors, 'high');
    xlabel(sprintf('High-speed time in revolution %d (ms)', selectedRevId));
    ylabel(sprintf('CH%d V (V)', sid));
    title(sprintf('CH%d high-speed local sequence: L1-L6, corr=%.3f', ...
        sid, S.corr_norm_peak_phys), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 8.5);
    subtitle(sprintf('P->L %s | L->P %s', ...
        mat2str(S.physical_to_local), mat2str(S.local_to_physical)), ...
        'FontWeight', 'normal', 'FontName', 'Times New Roman', 'FontSize', 6.8);
    style_axes_local();
end

sgtitle(sprintf(['Six-blade low/high waveform audit with polynomial-fitted peaks, ' ...
    'B%d, W%02d lap %d'], P.targetBlade, P.windowID, P.revolutionInWindow), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

if P.saveFigures
    pngFile = fullfile(figureDir, sprintf( ...
        'Step05E_SixBladeLowHighSequence_B%d_%s_W%02d_R%02d_20241106.png', ...
        P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID, P.revolutionInWindow));
    exportgraphics(fig, pngFile, 'Resolution', 300);
    if isfield(P, 'savePdf') && P.savePdf
        pdfFile = fullfile(figureDir, sprintf( ...
            'Step05E_SixBladeLowHighSequence_B%d_%s_W%02d_R%02d_20241106.pdf', ...
            P.targetBlade, ['S', sprintf('%d', P.analysisSensors)], P.windowID, P.revolutionInWindow));
        exportgraphics(fig, pdfFile, 'ContentType', 'vector');
    end
end
end

function corrValue = lookup_low_fingerprint_corr_local(summaryTable, sid)
corrValue = NaN;
if isempty(summaryTable) || ~ismember('sensor_id', summaryTable.Properties.VariableNames)
    return;
end
row = summaryTable(summaryTable.sensor_id == sid, :);
if height(row) ~= 1 || ~ismember('best_corr', row.Properties.VariableNames)
    return;
end
corrValue = row.best_corr(1);
end

function plot_sequence_panel_local(raw, segments, bladeColors, modeText)
tRef = min(raw.T);
plot(1000 * (raw.T - tRef), raw.V, 'Color', [0.65, 0.65, 0.65], 'LineWidth', 0.7, ...
    'DisplayName', 'raw waveform');
peakY = nan(numel(segments), 1);
for i = 1:numel(segments)
    seg = segments(i);
    c = bladeColors(seg.physical_blade, :);
    plot(1000 * (seg.t_fit - tRef), seg.v_fit, '-', 'Color', c, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    scatter(1000 * (seg.peak_t - tRef), seg.peak_v, 20, c, 'filled', 'HandleVisibility', 'off');
    peakY(i) = seg.peak_v;

    if strcmpi(modeText, 'low')
        labelText = sprintf('B%d', seg.physical_blade);
    else
        labelText = sprintf('L%d/P%d', seg.local_slot, seg.physical_blade);
    end
    text(1000 * (seg.peak_t - tRef), seg.peak_v, ['  ' labelText], ...
        'Color', c, 'FontName', 'Times New Roman', 'FontSize', 7.2, ...
        'VerticalAlignment', 'bottom', 'Clipping', 'on');
end
pad_sequence_axis_for_labels_local(gca, raw.V, peakY);
end

function pad_sequence_axis_for_labels_local(ax, rawV, peakY)
valid = isfinite(rawV);
if ~any(valid)
    return;
end
yMin = min(rawV(valid));
yMax = max([rawV(valid); peakY(isfinite(peakY))]);
yRange = max(yMax - yMin, eps);
ylim(ax, [yMin - 0.08 * yRange, yMax + 0.30 * yRange]);
end

function [tSec, v] = load_raw_channel_local(caseDir, channelId, fileId, pinlv)
filepath = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function [pulseInfo, nextTail] = extract_sensor_pulses_local(t, v, threshold, gapPoints, tail)
if nargin < 5 || isempty(tail)
    tail.t = [];
    tail.v = [];
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

smoothV = smooth_signal_local(v);
idxAbove = find(smoothV > threshold);
pulseInfo = struct('arrival_times', [], 'peak_features', [], 'start_times', [], 'end_times', []);
nextTail = struct('t', [], 'v', []);
if isempty(idxAbove)
    return;
end

[segStart, segEnd] = split_segments_local(idxAbove, gapPoints);
if numel(v) - idxAbove(end) < gapPoints
    carryStart = max(1, segStart(end) - gapPoints);
    nextTail.t = t(carryStart:end);
    nextTail.v = v(carryStart:end);
    segStart(end) = [];
    segEnd(end) = [];
end
if isempty(segStart)
    return;
end

nSeg = numel(segStart);
arrivalTimes = nan(nSeg, 1);
peakFeatures = nan(nSeg, 1);
startTimes = nan(nSeg, 1);
endTimes = nan(nSeg, 1);
for i = 1:nSeg
    idxRange = segStart(i):segEnd(i);
    tSeg = t(idxRange);
    vSeg = v(idxRange);
    smoothSeg = smoothV(idxRange);
    [arrivalTimes(i), peakFeatures(i)] = fit_arrival_time_local(tSeg, vSeg, smoothSeg, threshold);
    startTimes(i) = tSeg(1);
    endTimes(i) = tSeg(end);
end

pulseInfo.arrival_times = arrivalTimes;
pulseInfo.peak_features = peakFeatures;
pulseInfo.start_times = startTimes;
pulseInfo.end_times = endTimes;
end

function smoothV = smooth_signal_local(v)
if numel(v) >= 21
    smoothV = sgolayfilt(v, 3, 21);
elseif numel(v) >= 5
    span = max(3, 2 * floor(numel(v) / 4) + 1);
    smoothV = smoothdata(v, 'movmean', span);
else
    smoothV = v;
end
end

function [segStart, segEnd] = split_segments_local(idxAbove, gapPoints)
jumps = find(diff(idxAbove) > gapPoints);
segStart = [idxAbove(1); idxAbove(jumps + 1)];
segEnd = [idxAbove(jumps); idxAbove(end)];
end

function [arrivalTime, peakFeature] = fit_arrival_time_local(tSeg, vSeg, smoothSeg, threshold)
validMask = smoothSeg > threshold;
if sum(validMask) < 4
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
    return;
end

tFit = tSeg(validMask);
vFit = vSeg(validMask);
mu = mean(tFit);
uniqueCount = numel(unique(tFit));
order = min(3, uniqueCount - 1);
if order < 1
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
    return;
end

try
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    p = polyfit(tFit - mu, vFit, order);
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    dt = linspace(min(tFit) - mu, max(tFit) - mu, 100);
    vf = polyval(p, dt);
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;
    if sum(weights) > 0
        arrivalTime = sum(tf .* weights) / sum(weights);
        peakFeature = max(vf);
    else
        [peakFeature, idxMax] = max(vf);
        arrivalTime = tf(idxMax);
    end
catch
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    [peakFeature, idxMax] = max(vSeg);
    arrivalTime = tSeg(idxMax);
end
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
    'prev_opr_index', [], ...
    'revolution_index', [], ...
    'opr_times', []), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    loadedProbe = load(fullfile(caseOutputDir, sprintf('jilublade_probe%d.mat', sid)), 'jilublade');
    jilublade = loadedProbe.jilublade;
    [relAnglesDeg, prevOprIndex, revolutionIndex] = ...
        compute_pulse_reference_geometry_local(jilublade(:, 3), oprTimes, oprPulsesPerRev);
    probeCache(is).sensor_id = sid;
    probeCache(is).jilublade = jilublade;
    probeCache(is).rel_angles_deg = relAnglesDeg;
    probeCache(is).prev_opr_index = prevOprIndex;
    probeCache(is).revolution_index = revolutionIndex;
    probeCache(is).opr_times = oprTimes;
end
end

function oprTimes = load_opr_center_times_local(caseOutputDir)
loaded = load(fullfile(caseOutputDir, 'jiluOPR.mat'), 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function selectedRevIds = get_window_anchor_revolutions_local(DynamicMap, probeCache, windowID)
anchorSid = DynamicMap.SelectionInfo.anchor_sensor_id;
anchorRows = DynamicMap.SelectionInfo.anchor_selected_rows(:);
lapRange = DynamicMap.Window(windowID).lap_range;
anchorProbe = probeCache([probeCache.sensor_id] == anchorSid);
selectedRevIds = anchorProbe.revolution_index(anchorRows(lapRange));
selectedRevIds = selectedRevIds(isfinite(selectedRevIds));
end

function [relAnglesDeg, prevOprIndex, revolutionIndex] = compute_pulse_reference_geometry_local(arrivalTimes, oprTimes, oprPulsesPerRev)
relAnglesDeg = nan(size(arrivalTimes));
prevOprIndex = nan(size(arrivalTimes));
revolutionIndex = nan(size(arrivalTimes));
if isempty(arrivalTimes) || numel(oprTimes) <= oprPulsesPerRev
    return;
end
for i = 1:numel(arrivalTimes)
    tMeas = arrivalTimes(i);
    idx = find(oprTimes < tMeas, 1, 'last');
    if isempty(idx) || idx + oprPulsesPerRev > numel(oprTimes)
        continue;
    end
    t0 = oprTimes(idx);
    t1 = oprTimes(idx + oprPulsesPerRev);
    dtRev = t1 - t0;
    if ~isfinite(dtRev) || dtRev <= eps
        continue;
    end
    relAnglesDeg(i) = 360 * (tMeas - t0) / dtRev;
    prevOprIndex(i) = idx;
    revolutionIndex(i) = idx;
end
end

function value = local_corrcoef_scalar_local(a, b)
if numel(a) < 3 || numel(b) < 3
    value = nan;
    return;
end
C = corrcoef(a(:), b(:));
if numel(C) < 4
    value = nan;
else
    value = C(1, 2);
end
end

function style_axes_local()
grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'LineWidth', 0.8);
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
