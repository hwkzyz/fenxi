clear; clc; close all;

data_root = 'E:\试验数据\20250527\应变片数据20250527';
target_folders = {
    fullfile(data_root, '20250527_3150')
    fullfile(data_root, '20250527_2500_3500_time200')
    fullfile(data_root, '20250527_2500_3500_time400')
    fullfile(data_root, '20250527_2500_3500_time800')
    };

channel_tags = {'AI1-03', 'AI1-04'};

window_sec = 0.25;
hop_sec = 0.05;
fft_low_cut_hz = 20;
fft_high_cut_hz = 1200;
score_threshold_ratio = 0.72;
score_min_quantile = 0.85;
min_region_duration_sec = 0.25;
max_regions_to_plot = 3;
local_fft_num_segments = 80;

fprintf('=== 20250527 应变片无先验共振候选区诊断 ===\n');
fprintf('方法: 短时能量 + RMS + 包络 + 局部FFT峰值 + 频谱熵\n');

for iFolder = 1:numel(target_folders)
    folderpath = target_folders{iFolder};
    folder_label = get_folder_label(folderpath);

    for iChan = 1:numel(channel_tags)
        channel_tag = channel_tags{iChan};
        if ~has_channel_file(folderpath, channel_tag)
            continue;
        end

        filepath = pick_latest_channel_file(folderpath, channel_tag);
        [t, x, Fs] = load_strain_signal(filepath);

        [tc, features] = compute_short_time_features(x, Fs, window_sec, hop_sec, fft_low_cut_hz, fft_high_cut_hz);
        score = build_resonance_score(features);

        regions = find_candidate_regions(tc, score, score_threshold_ratio, score_min_quantile, min_region_duration_sec);
        regions = sort_regions_by_score(regions);
        for k = 1:numel(regions)
            idx = regions(k).index_range;
            [~, local_idx] = max(score(idx));
            peak_idx = idx(local_idx);
            regions(k).peak_freq = features.peak_freq(peak_idx);
        end

        fprintf('\n[%s - %s]\n', folder_label, channel_tag);
        fprintf('  Fs = %.1f Hz, duration = %.3f s\n', Fs, t(end) - t(1));
        fprintf('  Candidate regions found = %d\n', numel(regions));
        for k = 1:numel(regions)
            fprintf('  Region %d: [%.3f, %.3f] s, duration = %.3f s, peak score = %.3f, dominant freq ~= %.3f Hz\n', ...
                k, regions(k).t_start, regions(k).t_end, regions(k).duration, regions(k).peak_score, regions(k).peak_freq);
        end

        h = figure('Name', sprintf('STE Resonance Candidates - %s - %s', folder_label, channel_tag), ...
            'Color', 'w', 'NumberTitle', 'off');

        ax1 = subplot(4, 1, 1, 'Parent', h);
        plot(ax1, t, x, 'k');
        grid(ax1, 'on');
        xlabel(ax1, 'Time (s)');
        ylabel(ax1, 'Strain voltage');
        title(ax1, sprintf('%s - raw strain (%s)', folder_label, channel_tag), 'Interpreter', 'none');
        hold(ax1, 'on');
        mark_regions(ax1, regions, max_regions_to_plot);

        ax2 = subplot(4, 1, 2, 'Parent', h);
        plot(ax2, tc, features.ste, 'b', 'LineWidth', 1.0);
        hold(ax2, 'on');
        plot(ax2, tc, features.rms, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0);
        plot(ax2, tc, features.env_mean, 'Color', [0.47, 0.67, 0.19], 'LineWidth', 1.0);
        grid(ax2, 'on');
        xlabel(ax2, 'Time (s)');
        ylabel(ax2, 'Feature value');
        title(ax2, 'Short-time energy / RMS / envelope mean');
        legend(ax2, {'STE', 'RMS', 'Envelope mean'}, 'Location', 'best');
        mark_regions(ax2, regions, max_regions_to_plot);

        ax3 = subplot(4, 1, 3, 'Parent', h);
        yyaxis(ax3, 'left');
        plot(ax3, tc, features.peak_amp, 'Color', [0 0.45 0.74], 'LineWidth', 1.0);
        ylabel(ax3, 'FFT peak amplitude');
        yyaxis(ax3, 'right');
        plot(ax3, tc, features.peak_freq, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
        ylabel(ax3, 'Peak frequency (Hz)');
        grid(ax3, 'on');
        xlabel(ax3, 'Time (s)');
        title(ax3, 'Local FFT peak tracking');
        mark_regions(ax3, regions, max_regions_to_plot);

        ax4 = subplot(4, 1, 4, 'Parent', h);
        yyaxis(ax4, 'left');
        plot(ax4, tc, score, 'm', 'LineWidth', 1.2);
        ylabel(ax4, 'Resonance score');
        yyaxis(ax4, 'right');
        plot(ax4, tc, features.spec_entropy, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0);
        ylabel(ax4, 'Spectral entropy');
        grid(ax4, 'on');
        xlabel(ax4, 'Time (s)');
        title(ax4, 'Score and spectral entropy');
        mark_regions(ax4, regions, max_regions_to_plot);

        drawnow;

        nPlot = min(max_regions_to_plot, numel(regions));
        for k = 1:nPlot
            plot_local_3d_fft(t, x, Fs, regions(k), folder_label, channel_tag, local_fft_num_segments);
        end
    end
end

fprintf('\n建议解释方式:\n');
fprintf('1. 先看 STE/RMS/包络和综合分数最高的连续区域, 把它们当作候选共振区。\n');
fprintf('2. 再看这些候选区对应的局部三维频谱图, 只保留有清晰窄带峰脊的区域。\n');
fprintf('3. 之后再用这些时间窗去对齐 BTT 数据做标定后辨识。\n');


function tf = has_channel_file(folderpath, channel_tag)
d = dir(fullfile(folderpath, sprintf('%s_*.mat', channel_tag)));
tf = ~isempty(d);
end


function filepath = pick_latest_channel_file(folderpath, channel_tag)
pattern = sprintf('%s_*.mat', channel_tag);
d = dir(fullfile(folderpath, pattern));
if isempty(d)
    error('No file found for %s in %s', channel_tag, folderpath);
end

names = {d.name};
is_new = contains(names, '20250527');
if any(is_new)
    d = d(is_new);
end

[~, idx] = max([d.datenum]);
filepath = fullfile(folderpath, d(idx).name);
end


function [t, x, Fs] = load_strain_signal(filepath)
S = load(filepath);
t = S.Datas(:,1);
x = S.Datas(:,2);
if isfield(S, 'SampleFrequency')
    Fs = str2double(S.SampleFrequency);
else
    Fs = 1 / median(diff(t));
end
end


function [tc, features] = compute_short_time_features(x, Fs, window_sec, hop_sec, fft_low_cut_hz, fft_high_cut_hz)
x = detrend(x(:));
win = max(64, round(window_sec * Fs));
hop = max(8, round(hop_sec * Fs));
nFrames = floor((numel(x) - win) / hop) + 1;

tc = zeros(nFrames, 1);
ste = zeros(nFrames, 1);
rmsv = zeros(nFrames, 1);
envm = zeros(nFrames, 1);
peak_amp = zeros(nFrames, 1);
peak_freq = zeros(nFrames, 1);
spec_entropy = zeros(nFrames, 1);

freq = Fs * (0:floor(win/2)) / win;
freq_mask = freq >= fft_low_cut_hz & freq <= fft_high_cut_hz;

for i = 1:nFrames
    idx = (i-1)*hop + (1:win);
    seg = x(idx);
    tc(i) = ((idx(1) + idx(end)) / 2 - 1) / Fs;

    ste(i) = sum(seg.^2);
    rmsv(i) = sqrt(mean(seg.^2));
    envm(i) = mean(abs(hilbert(seg)));

    segw = seg .* hann(win);
    Y = fft(segw);
    P2 = abs(Y / win);
    P1 = P2(1:floor(win/2)+1);
    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end

    Pband = P1(freq_mask);
    fband = freq(freq_mask);
    [peak_amp(i), idx_pk] = max(Pband);
    peak_freq(i) = fband(idx_pk);

    p = Pband .^ 2;
    p = p / max(sum(p), eps);
    spec_entropy(i) = -sum(p .* log(p + eps));
end

features.ste = ste;
features.rms = rmsv;
features.env_mean = envm;
features.peak_amp = peak_amp;
features.peak_freq = peak_freq;
features.spec_entropy = spec_entropy;
end


function score = build_resonance_score(features)
z_ste = robust_standardize(features.ste);
z_rms = robust_standardize(features.rms);
z_env = robust_standardize(features.env_mean);
z_peak = robust_standardize(features.peak_amp);
z_ent = robust_standardize(features.spec_entropy);

score = z_ste + 0.8*z_rms + 0.8*z_env + 1.0*z_peak - 0.5*z_ent;
score = movmean(score, 5);
end


function z = robust_standardize(x)
x = x(:);
medv = median(x, 'omitnan');
madv = median(abs(x - medv), 'omitnan');
scale = 1.4826 * max(madv, eps);
z = (x - medv) / scale;
end


function regions = find_candidate_regions(tc, score, threshold_ratio, min_quantile, min_duration_sec)
score = score(:).';
tc = tc(:).';
if isempty(score)
    regions = struct('t_start', {}, 't_end', {}, 'duration', {}, 'peak_score', {}, 'peak_freq', {}, 'index_range', {});
    return;
end

thr = max(threshold_ratio * max(score), quantile(score, min_quantile));
mask = score >= thr;
regions = build_regions_from_mask(mask, tc, score, min_duration_sec);

if isempty(regions)
    [~, idx_max] = max(score);
    fallback_mask = false(size(mask));
    fallback_mask(max(1, idx_max-1):min(numel(mask), idx_max+1)) = true;
    regions = build_regions_from_mask(fallback_mask, tc, score, 0);
end
end


function regions = build_regions_from_mask(mask, tc, score, min_duration_sec)
d = diff([false, mask, false]);
start_idx = find(d == 1);
end_idx = find(d == -1) - 1;

regions = struct('t_start', {}, 't_end', {}, 'duration', {}, 'peak_score', {}, 'peak_freq', {}, 'index_range', {});
count = 0;
for k = 1:numel(start_idx)
    idx = start_idx(k):end_idx(k);
    t_start = tc(idx(1));
    t_end = tc(idx(end));
    duration = t_end - t_start;
    if duration < min_duration_sec
        continue;
    end
    count = count + 1;
    regions(count).t_start = t_start;
    regions(count).t_end = t_end;
    regions(count).duration = duration;
    regions(count).peak_score = max(score(idx));
    regions(count).peak_freq = NaN;
    regions(count).index_range = idx;
end
end


function regions = sort_regions_by_score(regions)
if isempty(regions)
    return;
end
[~, order] = sort([regions.peak_score], 'descend');
regions = regions(order);
end


function mark_regions(ax, regions, max_regions_to_plot)
nPlot = min(max_regions_to_plot, numel(regions));
for k = 1:nPlot
    xline(ax, regions(k).t_start, '--r', sprintf('R%d start', k));
    xline(ax, regions(k).t_end, '--r', sprintf('R%d end', k));
end
end


function plot_local_3d_fft(t, x, Fs, region, folder_label, channel_tag, num_segments)
margin_sec = max(0.2, 0.2 * max(region.duration, 0.2));
t1 = max(t(1), region.t_start - margin_sec);
t2 = min(t(end), region.t_end + margin_sec);
mask = t >= t1 & t <= t2;
t_local = t(mask);
x_local = x(mask);

if numel(x_local) < num_segments * 8
    return;
end

x_local = detrend(x_local(:));
N = numel(x_local);
segment_length = floor(N / num_segments);
if segment_length < 16
    return;
end

usable_length = segment_length * num_segments;
x_local = x_local(1:usable_length);
t_local = t_local(1:usable_length);

freq = Fs * (0:floor(segment_length/2)) / segment_length;
amplitudes = zeros(num_segments, numel(freq));
segment_t = zeros(num_segments, 1);

for i = 1:num_segments
    idx = (i-1)*segment_length + (1:segment_length);
    segment_data = x_local(idx) .* hann(segment_length);
    segment_t(i) = mean(t_local(idx));
    segment_fft = fft(segment_data);
    segment_p2 = abs(segment_fft / segment_length);
    segment_p1 = segment_p2(1:floor(segment_length/2)+1);
    if numel(segment_p1) > 2
        segment_p1(2:end-1) = 2 * segment_p1(2:end-1);
    end
    segment_p1(1) = 0;
    amplitudes(i, :) = segment_p1;
end

h = figure('Name', sprintf('Local 3D FFT - %s - %s - %.3f to %.3f s', ...
    folder_label, channel_tag, region.t_start, region.t_end), ...
    'Color', 'w', 'NumberTitle', 'off');

ax1 = subplot(2, 1, 1, 'Parent', h);
plot(ax1, t_local, x_local, 'k');
hold(ax1, 'on');
xline(ax1, region.t_start, '--r', 'Res start');
xline(ax1, region.t_end, '--r', 'Res end');
grid(ax1, 'on');
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Strain voltage');
title(ax1, sprintf('%s - %s local strain', folder_label, channel_tag), 'Interpreter', 'none');

ax2 = subplot(2, 1, 2, 'Parent', h);
surf(ax2, segment_t, freq, amplitudes.', 'EdgeColor', 'none');
ylim(ax2, [0, 1200]);
view(ax2, 30, 30);
grid(ax2, 'on');
colormap(ax2, jet);
colorbar(ax2);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Frequency (Hz)');
zlabel(ax2, 'Amplitude');
title(ax2, 'Local 3D FFT (fft1-style)');
end


function label = get_folder_label(folderpath)
[~, label] = fileparts(folderpath);
end
