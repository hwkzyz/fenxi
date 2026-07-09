%% Visualize OPR anchor definitions for 20250527
% This script opens diagnostic figures only. It does not save image files.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
legacyDir = fullfile(thisDir, 'legacy');
addpath(legacyDir);

cfg = Get_20250527_BTT_Config();
caseName = cfg.dynamic_cases{1};
caseDir = fullfile(cfg.dataset_root, caseName);
maxPulses = parse_positive_integer_env_local('OPR_VIS_MAX_PULSES', 300);

fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR raw files found in %s.', caseDir);
end

[pulseTable, rawPreview] = extract_opr_anchor_table_local(caseDir, fileIds, cfg, maxPulses);
if isempty(pulseTable)
    error('No OPR pulses detected for %s.', caseName);
end

fprintf('\n=== OPR anchor definition diagnostic ===\n');
fprintf('Case: %s\nOPR channel: CH%d, fs = %.3g Hz, r_tip = %.3f mm\n', ...
    caseName, cfg.opr_id, cfg.pinlv, cfg.r_tip_mm);
fprintf('Legacy jiluOPR(:,1) uses segment start sample / fs, i.e. threshold rising-edge start.\n');
fprintf('Detected pulses shown: %d\n\n', height(pulseTable));

disp(pulseTable(1:min(12, height(pulseTable)), ...
    {'pulse_id','legacy_start_s','rise_interp_s','fall_interp_s','width_center_s','multi_center_s', ...
    'center_minus_legacy_us','multi_minus_legacy_us','center_minus_legacy_mm','multi_minus_legacy_mm'}));

summaryRows = table( ...
    ["width center - legacy"; "multi-threshold center - legacy"], ...
    [median(pulseTable.center_minus_legacy_us, 'omitnan'); median(pulseTable.multi_minus_legacy_us, 'omitnan')], ...
    [std(pulseTable.center_minus_legacy_us, 'omitnan'); std(pulseTable.multi_minus_legacy_us, 'omitnan')], ...
    [median(pulseTable.center_minus_legacy_mm, 'omitnan'); median(pulseTable.multi_minus_legacy_mm, 'omitnan')], ...
    [std(pulseTable.center_minus_legacy_mm, 'omitnan'); std(pulseTable.multi_minus_legacy_mm, 'omitnan')], ...
    'VariableNames', {'definition','median_dt_us','std_dt_us','median_dx_mm','std_dx_mm'});
disp(summaryRows);

plot_opr_anchor_preview_local(rawPreview, pulseTable, cfg);
plot_opr_anchor_offsets_local(pulseTable);

fprintf('\nInterpretation:\n');
fprintf('If the center-minus-legacy offset is not negligible, OPRAnchored currently carries a systematic x-origin shift.\n');
fprintf('A constant OPR timing shift maps to x shift by dx = dt * 2*pi*r_tip*f_rot.\n');
fprintf('No figures were saved. Close the figure windows when finished.\n');

%% Local functions
function [pulseTable, rawPreview] = extract_opr_anchor_table_local(caseDir, fileIds, cfg, maxPulses)
tail = [];
gapPoints = 1e4;
rows = repmat(struct( ...
    'pulse_id', NaN, 'legacy_start_s', NaN, 'rise_interp_s', NaN, ...
    'fall_interp_s', NaN, 'width_center_s', NaN, 'multi_center_s', NaN, ...
    'center_minus_legacy_us', NaN, 'multi_minus_legacy_us', NaN, ...
    'center_minus_legacy_mm', NaN, 'multi_minus_legacy_mm', NaN, ...
    'rot_freq_hz', NaN, 'peak_v', NaN), maxPulses, 1);
rowCount = 0;
rawPreview = struct('t', [], 'v', [], 'threshold', cfg.opr_threshold);

for iFile = 1:numel(fileIds)
    raw = load_raw_case_channel_local(caseDir, cfg.opr_id, fileIds(iFile));
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    if isempty(rawPreview.t)
        nPreview = min(size(raw, 1), round(0.08 * cfg.pinlv));
        rawPreview.t = raw(1:nPreview, 1) / cfg.pinlv;
        rawPreview.v = raw(1:nPreview, 2);
    end

    [segments, tail, gapPoints] = segment_signal_like_step2_local(raw, cfg.opr_threshold, gapPoints);
    for iSeg = 1:numel(segments.start_idx)
        if rowCount >= maxPulses
            break;
        end
        rowCount = rowCount + 1;
        seg = characterize_opr_segment_local(raw, segments.start_idx(iSeg), segments.end_idx(iSeg), cfg);
        rows(rowCount).pulse_id = rowCount;
        rows(rowCount).legacy_start_s = segments.start_sample(iSeg) / cfg.pinlv;
        rows(rowCount).rise_interp_s = seg.rise_s;
        rows(rowCount).fall_interp_s = seg.fall_s;
        rows(rowCount).width_center_s = seg.center_s;
        rows(rowCount).multi_center_s = seg.multi_center_s;
        rows(rowCount).center_minus_legacy_us = 1e6 * (seg.center_s - rows(rowCount).legacy_start_s);
        rows(rowCount).multi_minus_legacy_us = 1e6 * (seg.multi_center_s - rows(rowCount).legacy_start_s);
        rows(rowCount).peak_v = seg.peak_v;
    end
    if rowCount >= maxPulses
        break;
    end
end

rows = rows(1:rowCount);
legacy = [rows.legacy_start_s].';
rotFreq = estimate_rot_freq_by_pulse_local(legacy, cfg.blades_num);
for i = 1:rowCount
    rows(i).rot_freq_hz = rotFreq(i);
    scale = 2 * pi * cfg.r_tip_mm * rotFreq(i);
    rows(i).center_minus_legacy_mm = (rows(i).width_center_s - rows(i).legacy_start_s) * scale;
    rows(i).multi_minus_legacy_mm = (rows(i).multi_center_s - rows(i).legacy_start_s) * scale;
end
pulseTable = struct2table(rows);
end

function seg = characterize_opr_segment_local(raw, a, b, cfg)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(size(raw, 1), b + pad);
t = raw(a0:b0, 1) / cfg.pinlv;
v = raw(a0:b0, 2);
try
    vs = smooth(v, 16);
catch
    vs = smoothdata(v, 'movmean', 16);
end
[peakV, iPeak] = max(vs);
baseV = median(vs(1:max(3, min(20, numel(vs)))), 'omitnan');
rise = threshold_crossing_time_local(t(1:iPeak), vs(1:iPeak), cfg.opr_threshold, 'rising');
fall = threshold_crossing_time_local(t(iPeak:end), vs(iPeak:end), cfg.opr_threshold, 'falling');
levels = baseV + [0.30 0.40 0.50 0.60 0.70] .* max(peakV - baseV, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = threshold_crossing_time_local(t(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = threshold_crossing_time_local(t(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
seg = struct();
seg.rise_s = rise;
seg.fall_s = fall;
seg.center_s = 0.5 * (rise + fall);
seg.multi_center_s = median(centers, 'omitnan');
seg.peak_v = peakV;
end

function tc = threshold_crossing_time_local(t, v, threshold, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < threshold & v(2:end) >= threshold, 1, 'first');
else
    idx = find(v(1:end-1) >= threshold & v(2:end) < threshold, 1, 'last');
end
if isempty(idx)
    [~, idx] = min(abs(v - threshold));
    tc = t(idx);
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (threshold - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end

function rotFreq = estimate_rot_freq_by_pulse_local(oprTimes, pulsesPerRev)
rotFreq = nan(size(oprTimes));
for i = 1:numel(oprTimes)
    j0 = i - pulsesPerRev;
    j1 = i + pulsesPerRev;
    vals = [];
    if j0 >= 1
        vals(end+1) = 1 / max(oprTimes(i) - oprTimes(j0), eps); %#ok<AGROW>
    end
    if j1 <= numel(oprTimes)
        vals(end+1) = 1 / max(oprTimes(j1) - oprTimes(i), eps); %#ok<AGROW>
    end
    rotFreq(i) = median(vals, 'omitnan');
end
rotFreq = fillmissing(rotFreq, 'nearest');
end

function plot_opr_anchor_preview_local(rawPreview, pulseTable, cfg)
fig = figure('Name', 'OPR waveform anchor definitions', 'Color', 'w');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot(rawPreview.t, rawPreview.v, 'k-', 'LineWidth', 1.0, 'DisplayName', 'OPR raw');
yline(cfg.opr_threshold, '--', 'DisplayName', 'legacy threshold');
xlim([min(rawPreview.t), max(rawPreview.t)]);
xlabel('Time (s)'); ylabel('V');
title('Raw OPR waveform preview');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on; grid on; box on;
nShow = min(8, height(pulseTable));
for i = 1:nShow
    xline(pulseTable.legacy_start_s(i), '-', 'Color', [0.85 0.20 0.15], 'LineWidth', 0.9);
    xline(pulseTable.width_center_s(i), '--', 'Color', [0.10 0.35 0.80], 'LineWidth', 0.9);
    xline(pulseTable.multi_center_s(i), ':', 'Color', [0.10 0.55 0.20], 'LineWidth', 1.0);
end
plot(rawPreview.t, rawPreview.v, 'k-', 'LineWidth', 0.9);
yline(cfg.opr_threshold, '--', 'Color', [0.4 0.4 0.4]);
xlim([min(pulseTable.legacy_start_s(1), rawPreview.t(1)), pulseTable.legacy_start_s(nShow) + 0.02]);
xlabel('Time (s)'); ylabel('V');
title('Red: legacy start | Blue: width center | Green: multi-threshold center');
drawnow;
end

function plot_opr_anchor_offsets_local(pulseTable)
fig = figure('Name', 'OPR anchor time and spatial offsets', 'Color', 'w');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot(pulseTable.pulse_id, pulseTable.center_minus_legacy_us, '-', 'LineWidth', 1.0, 'DisplayName', 'width center');
plot(pulseTable.pulse_id, pulseTable.multi_minus_legacy_us, '-', 'LineWidth', 1.0, 'DisplayName', 'multi center');
xlabel('OPR pulse'); ylabel('\Delta t vs legacy (\mus)', 'Interpreter', 'tex');
title('Timing shift from legacy rising-edge start');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on; grid on; box on;
histogram(pulseTable.center_minus_legacy_us, 30, 'DisplayName', 'width center');
histogram(pulseTable.multi_minus_legacy_us, 30, 'DisplayName', 'multi center');
xlabel('\Delta t vs legacy (\mus)', 'Interpreter', 'tex'); ylabel('Count');
title('Timing offset distribution');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on; grid on; box on;
plot(pulseTable.pulse_id, pulseTable.center_minus_legacy_mm, '-', 'LineWidth', 1.0, 'DisplayName', 'width center');
plot(pulseTable.pulse_id, pulseTable.multi_minus_legacy_mm, '-', 'LineWidth', 1.0, 'DisplayName', 'multi center');
xlabel('OPR pulse'); ylabel('\Delta x at blade tip (mm)', 'Interpreter', 'tex');
title('Equivalent spatial anchor shift');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on; grid on; box on;
histogram(pulseTable.center_minus_legacy_mm, 30, 'DisplayName', 'width center');
histogram(pulseTable.multi_minus_legacy_mm, 30, 'DisplayName', 'multi center');
xlabel('\Delta x at blade tip (mm)', 'Interpreter', 'tex'); ylabel('Count');
title('Spatial offset distribution');
legend('Location', 'best', 'Box', 'off');
drawnow;
end

function raw = load_raw_case_channel_local(caseDir, channelId, fileId)
file = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
if ~isfile(file)
    raw = [];
    return;
end
loaded = load(file);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
end

function [segments, tail, nextGap] = segment_signal_like_step2_local(raw, threshold, gapPoints)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];
nextGap = gapPoints;
if isempty(raw)
    return;
end
sig = raw(:, 2);
try
    sigSmooth = smooth(sig, 16);
catch
    sigSmooth = smoothdata(sig, 'movmean', 16);
end
chase = find(sigSmooth > threshold);
if isempty(chase)
    return;
end
sampleOrder = raw(chase, 1);
segStarts = 1;
segEnds = [];
for ii = 1:(numel(sampleOrder) - 1)
    c = sampleOrder(ii + 1) - sampleOrder(ii);
    if c > nextGap
        segEnds(end + 1, 1) = ii; %#ok<AGROW>
        segStarts(end + 1, 1) = ii + 1; %#ok<AGROW>
        nextGap = 0.6 * c;
    end
end
segEnds(end + 1, 1) = numel(sampleOrder);
tailPoint = chase(segStarts(end)) - floor(nextGap / 2);
if tailPoint > 0 && tailPoint < size(raw, 1)
    tail = raw(tailPoint:end, :);
end
if numel(segStarts) < 2
    return;
end
compStarts = segStarts(1:end-1);
compEnds = segEnds(1:end-1);
segments.start_idx = chase(compStarts);
segments.end_idx = chase(compEnds);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end

function fileIds = list_case_file_ids_local(caseDir, channelId)
d = dir(fullfile(caseDir, sprintf('4-%d-*.mat', channelId)));
fileIds = [];
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channelId) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        fileIds(end + 1) = str2double(tok{1}); %#ok<AGROW>
    end
end
fileIds = sort(fileIds);
end

function n = parse_positive_integer_env_local(name, defaultValue)
txt = strtrim(getenv(name));
if isempty(txt)
    n = defaultValue;
    return;
end
v = str2double(txt);
if ~isfinite(v) || v < 1
    error('%s must be a positive integer.', name);
end
n = floor(v);
end
