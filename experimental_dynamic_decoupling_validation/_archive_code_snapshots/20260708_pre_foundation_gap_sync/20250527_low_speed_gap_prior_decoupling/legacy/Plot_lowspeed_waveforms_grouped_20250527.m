%% Plot low-speed waveforms for the 20250527 dataset
%
% This script is intended for direct visual inspection.
% It does not analyze the full 20 s record. Instead, it:
%   1. Converts the first column to time by t = sample_index / 5e6
%   2. Shows a short time-domain segment around one OPR pulse
%   3. Shows one full revolution using 6 OPR pulses per revolution
%
% If you want to save the figures manually later, do it from the MATLAB UI.

clc;
close all;

low_speed_dir = 'E:\试验数据\20250527\试验20250527\20250526_910';
pinlv = 5e6;
preview_file_id = 1000;
channels = 1:8;

channel_map = struct();
channel_map.capacitance = [1, 3, 6];
channel_map.opr = 4;
channel_map.other = [7, 8];
channel_map.unused = [2, 5];

opr_threshold_v = 2.0;
opr_gap_points = 1000;
opr_pulses_per_rev = 6;
short_half_window_s = 0.010;

if ~isfolder(low_speed_dir)
    error('Low-speed directory not found: %s', low_speed_dir);
end

WaveformData = struct();
SummaryRows = repmat(struct( ...
    'channel_id', NaN, ...
    'sample_count', NaN, ...
    'v_min', NaN, ...
    'v_max', NaN, ...
    'v_mean', NaN, ...
    'v_std', NaN, ...
    'group_name', ""), numel(channels), 1);

for i = 1:numel(channels)
    ch = channels(i);
    [t, v, file_name] = load_channel_raw_local(low_speed_dir, ch, preview_file_id, pinlv);
    WaveformData(ch).channel_id = ch; %#ok<SAGROW>
    WaveformData(ch).t = t;
    WaveformData(ch).v = v;
    WaveformData(ch).file_name = file_name;

    SummaryRows(i).channel_id = ch;
    SummaryRows(i).sample_count = numel(v);
    SummaryRows(i).v_min = min(v);
    SummaryRows(i).v_max = max(v);
    SummaryRows(i).v_mean = mean(v);
    SummaryRows(i).v_std = std(v);
    SummaryRows(i).group_name = string(classify_channel_local(ch, channel_map));
end

SummaryTable = struct2table(SummaryRows);
disp(SummaryTable);

t_opr = WaveformData(channel_map.opr).t;
v_opr = WaveformData(channel_map.opr).v;
opr_times = extract_opr_rising_edges_local(t_opr, v_opr, opr_threshold_v, opr_gap_points);
if numel(opr_times) < opr_pulses_per_rev + 1
    error('Not enough OPR pulses found in CH4 to define one revolution.');
end

last_valid_start = numel(opr_times) - opr_pulses_per_rev;
target_opr_idx = min(max(10, floor(last_valid_start / 3)), last_valid_start);

short_center_t = opr_times(target_opr_idx);
short_t0 = short_center_t - short_half_window_s;
short_t1 = short_center_t + short_half_window_s;

rev_t0 = opr_times(target_opr_idx);
rev_t1 = opr_times(target_opr_idx + opr_pulses_per_rev);

fprintf('Preview file id: %d\n', preview_file_id);
fprintf('Short window: [%.6f, %.6f] s\n', short_t0, short_t1);
fprintf('One-revolution window: [%.6f, %.6f] s\n', rev_t0, rev_t1);
fprintf('OPR pulses per revolution: %d\n', opr_pulses_per_rev);

fig_all_short = figure('Name', '20250527 Low-Speed All Channels - Short Window', ...
    'Color', 'w', 'Position', [60, 60, 1400, 900], 'Visible', 'on');
tiledlayout(fig_all_short, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ch = channels
    nexttile;
    hold on;
    grid on;
    box on;
    mask = WaveformData(ch).t >= short_t0 & WaveformData(ch).t <= short_t1;
    plot(WaveformData(ch).t(mask), WaveformData(ch).v(mask), 'b-', 'LineWidth', 0.8);
    title(sprintf('CH%d | %s', ch, classify_channel_local(ch, channel_map)));
    xlabel('Time (s)');
    ylabel('Voltage (V)');
end

fig_group_short = figure('Name', '20250527 Low-Speed Grouped - Short Window', ...
    'Color', 'w', 'Position', [100, 100, 1400, 900], 'Visible', 'on');
tiledlayout(fig_group_short, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_group_local(WaveformData, channel_map.capacitance, short_t0, short_t1, ...
    'Capacitance Sensors: CH1, CH3, CH6');
plot_group_local(WaveformData, channel_map.opr, short_t0, short_t1, ...
    'OPR Sensor: CH4');
plot_group_local(WaveformData, channel_map.other, short_t0, short_t1, ...
    'Other Displacement Sensors: CH7, CH8');

fig_group_rev = figure('Name', '20250527 Low-Speed Grouped - One Revolution', ...
    'Color', 'w', 'Position', [140, 140, 1400, 900], 'Visible', 'on');
tiledlayout(fig_group_rev, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_group_local(WaveformData, channel_map.capacitance, rev_t0, rev_t1, ...
    'Capacitance Sensors: CH1, CH3, CH6');
plot_group_local(WaveformData, channel_map.opr, rev_t0, rev_t1, ...
    'OPR Sensor: CH4');
plot_group_local(WaveformData, channel_map.other, rev_t0, rev_t1, ...
    'Other Displacement Sensors: CH7, CH8');

fig_unused = figure('Name', '20250527 Low-Speed Unused Channels - Short Window', ...
    'Color', 'w', 'Position', [180, 180, 1200, 360], 'Visible', 'on');
tiledlayout(fig_unused, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(channel_map.unused)
    ch = channel_map.unused(i);
    nexttile;
    hold on;
    grid on;
    box on;
    mask = WaveformData(ch).t >= short_t0 & WaveformData(ch).t <= short_t1;
    plot(WaveformData(ch).t(mask), WaveformData(ch).v(mask), 'k-', 'LineWidth', 0.8);
    title(sprintf('CH%d | likely unused / low-activity', ch));
    xlabel('Time (s)');
    ylabel('Voltage (V)');
end

drawnow;
fprintf('Figures opened:\n');
fprintf('  1. All channels in a short time window\n');
fprintf('  2. Grouped channels in a short time window\n');
fprintf('  3. Grouped channels in one revolution\n');
fprintf('  4. Unused channels in a short time window\n');

function [t, v, file_name] = load_channel_raw_local(data_dir, channel_id, target_file_id, pinlv)
file_name = fullfile(data_dir, sprintf('4-%d-%d.mat', channel_id, target_file_id));
if ~isfile(file_name)
    error('File not found: %s', file_name);
end

loaded = load(file_name);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];

t = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function t_cross = extract_opr_rising_edges_local(t, v, threshold, gap_points)
idx_above = find(v > threshold);
if isempty(idx_above)
    t_cross = [];
    return;
end

split_points = [1; find(diff(idx_above) > gap_points) + 1];
start_idx = idx_above(split_points);
t_cross = zeros(numel(start_idx), 1);

for i = 1:numel(start_idx)
    idx = start_idx(i);
    if idx <= 1
        t_cross(i) = t(idx);
        continue;
    end

    t1 = t(idx - 1);
    t2 = t(idx);
    v1 = v(idx - 1);
    v2 = v(idx);
    if abs(v2 - v1) < eps
        t_cross(i) = t2;
    else
        t_cross(i) = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
    end
end
end

function plot_group_local(WaveformData, channels, t0, t1, title_str)
nexttile;
hold on;
grid on;
box on;
colors = lines(numel(channels));
for i = 1:numel(channels)
    ch = channels(i);
    mask = WaveformData(ch).t >= t0 & WaveformData(ch).t <= t1;
    plot(WaveformData(ch).t(mask), WaveformData(ch).v(mask), '-', ...
        'Color', colors(i, :), 'LineWidth', 0.8, ...
        'DisplayName', sprintf('CH%d', ch));
end
title(title_str);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Location', 'best');
end

function group_name = classify_channel_local(ch, channel_map)
if ismember(ch, channel_map.capacitance)
    group_name = 'capacitance';
elseif ismember(ch, channel_map.opr)
    group_name = 'opr';
elseif ismember(ch, channel_map.other)
    group_name = 'other';
else
    group_name = 'unused_or_low_activity';
end
end
