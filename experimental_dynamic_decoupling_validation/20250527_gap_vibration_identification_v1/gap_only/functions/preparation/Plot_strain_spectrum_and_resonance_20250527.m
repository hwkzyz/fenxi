clear; clc; close all;

data_root = 'E:\试验数据\20250527\应变片数据20250527';

steady_folder = fullfile(data_root, '20250527_3150');
ramp_folders = {
    fullfile(data_root, '20250527_2500_3500_time200')
    fullfile(data_root, '20250527_2500_3500_time400')
    fullfile(data_root, '20250527_2500_3500_time800')
    };

steady_channel_tags = {'AI1-03', 'AI1-04'};
ramp_channel_tag = 'AI1-03';

fft_band_hz = [400, 800];
display_band_hz = [0, 1200];
ridge_search_band_hz = [450, 700];
ridge_amp_threshold = 3e5;
ridge_amp_smooth_frames = 7;
ridge_merge_gap_sec = 0.60;
ridge_min_duration_sec = 0.80;
ridge_freq_cluster_gap_hz = 18;

fprintf('=== 20250527 应变片频谱与共振区诊断 ===\n');

steady_files = cell(size(steady_channel_tags));
peak_freqs = nan(size(steady_channel_tags));
peak_amps = nan(size(steady_channel_tags));

for k = 1:numel(steady_channel_tags)
    steady_files{k} = pick_latest_channel_file(steady_folder, steady_channel_tags{k});
    [t, x, Fs] = load_strain_signal(steady_files{k});
    [f_fft, amp_fft] = compute_fft_spectrum(x, Fs);

    mask = f_fft >= fft_band_hz(1) & f_fft <= fft_band_hz(2);
    [peak_amps(k), idx_pk] = max(amp_fft(mask));
    f_band = f_fft(mask);
    peak_freqs(k) = f_band(idx_pk);

    h = figure('Name', sprintf('20250527 Steady %s FFT', steady_channel_tags{k}), ...
        'Color', 'w', 'NumberTitle', 'off');
    ax = axes('Parent', h);
    plot(ax, f_fft, amp_fft, 'LineWidth', 1.0);
    hold(ax, 'on');
    xline(ax, peak_freqs(k), '--r', sprintf('Peak %.3f Hz', peak_freqs(k)), ...
        'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
    xlim(ax, display_band_hz);
    grid(ax, 'on');
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Amplitude');
    title(ax, sprintf('Steady 3150 FFT - %s', steady_channel_tags{k}), 'Interpreter', 'none');

    fprintf('\nSteady file: %s\n', steady_files{k});
    fprintf('  Fs = %.1f Hz, duration = %.3f s\n', Fs, t(end) - t(1));
    fprintf('  Peak frequency in %.0f-%.0f Hz = %.3f Hz\n', ...
        fft_band_hz(1), fft_band_hz(2), peak_freqs(k));
end

reference_freq_hz = mean(peak_freqs, 'omitnan');

fprintf('\nReference resonance frequency from steady 3150 = %.3f Hz\n', reference_freq_hz);
fprintf('STFT ridge search band = [%.3f, %.3f] Hz\n', ridge_search_band_hz(1), ridge_search_band_hz(2));
fprintf('Ridge amplitude threshold = %.0f\n', ridge_amp_threshold);
fprintf('Ridge smoothing frames = %d, merge gap = %.2f s, min duration = %.2f s\n', ...
    ridge_amp_smooth_frames, ridge_merge_gap_sec, ridge_min_duration_sec);
fprintf('Ridge frequency clustering gap = %.1f Hz\n', ridge_freq_cluster_gap_hz);

for i = 1:numel(ramp_folders)
    ramp_file = pick_latest_channel_file(ramp_folders{i}, ramp_channel_tag);
    [t, x, Fs] = load_strain_signal(ramp_file);

    win = 8192;
    hop = 1024;
    nfft = 8192;
    [Sxx, F, T] = spectrogram(detrend(x), hann(win), win - hop, nfft, Fs, 'yaxis');
    P = abs(Sxx).^2;

    ridge_mask = F >= ridge_search_band_hz(1) & F <= ridge_search_band_hz(2);
    A_ridge = abs(Sxx(ridge_mask, :));
    F_ridge = F(ridge_mask);

    [ridge_amp, ridge_idx] = max(A_ridge, [], 1);
    ridge_amp = ridge_amp(:);
    ridge_freq = F_ridge(ridge_idx(:));
    ridge_amp_smooth = movmean(ridge_amp, ridge_amp_smooth_frames);

    [peak_amp, idx_peak] = max(ridge_amp);
    peak_time = T(idx_peak);
    peak_freq = ridge_freq(idx_peak);

    ridge_regions = build_ridge_regions(T, ridge_amp_smooth, ridge_amp, ridge_freq, ...
        ridge_amp_threshold, ridge_merge_gap_sec, ridge_min_duration_sec);
    [ridge_regions, ridge_groups] = cluster_ridge_regions_by_frequency(ridge_regions, ridge_freq_cluster_gap_hz);

    h = figure('Name', sprintf('20250527 Ramp %s', get_folder_label(ramp_folders{i})), ...
        'Color', 'w', 'NumberTitle', 'off');

    ax1 = subplot(4, 1, 1, 'Parent', h);
    plot(ax1, t, x, 'k');
    hold(ax1, 'on');
    mark_regions(ax1, ridge_regions);
    grid(ax1, 'on');
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Strain voltage');
    title(ax1, sprintf('%s - raw strain (%s)', get_folder_label(ramp_folders{i}), ramp_channel_tag), ...
        'Interpreter', 'none');

    ax2 = subplot(4, 1, 2, 'Parent', h);
    surf(ax2, T, F, 20*log10(abs(Sxx) + eps), 'EdgeColor', 'none');
    view(ax2, 45, 60);
    ylim(ax2, [0, 1200]);
    xlim(ax2, [T(1), T(end)]);
    hold(ax2, 'on');
    plot3(ax2, T, ridge_freq, max(20*log10(abs(Sxx) + eps), [], 1) + 5, 'w', 'LineWidth', 1.2);
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Frequency (Hz)');
    zlabel(ax2, 'Magnitude (dB)');
    title(ax2, '3D STFT spectrum with ridge');
    colormap(ax2, turbo);
    colorbar(ax2);

    ax3 = subplot(4, 1, 3, 'Parent', h);
    imagesc(ax3, T, F, 20*log10(abs(Sxx) + eps));
    axis(ax3, 'xy');
    ylim(ax3, [0, 1200]);
    hold(ax3, 'on');
    plot(ax3, T, ridge_freq, 'w', 'LineWidth', 1.2);
    yline(ax3, peak_freq, '--w', sprintf('Peak %.2f Hz', peak_freq));
    mark_regions(ax3, ridge_regions);
    colormap(ax3, turbo);
    colorbar(ax3);
    grid(ax3, 'on');
    xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'Frequency (Hz)');
    title(ax3, 'Spectrogram with resonance ridge');

    ax4 = subplot(4, 1, 4, 'Parent', h);
    yyaxis(ax4, 'left');
    plot(ax4, T, ridge_amp, 'Color', [0.60, 0.78, 1.00], 'LineWidth', 0.8);
    hold(ax4, 'on');
    plot(ax4, T, ridge_amp_smooth, 'b', 'LineWidth', 1.4);
    ylabel(ax4, 'Ridge amplitude');
    hold(ax4, 'on');
    yline(ax4, ridge_amp_threshold, '--r', sprintf('Threshold %.0f', ridge_amp_threshold));
    yyaxis(ax4, 'right');
    plot(ax4, T, ridge_freq, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0);
    ylabel(ax4, 'Ridge frequency (Hz)');
    xline(ax4, peak_time, '--k', 'Peak time');
    yyaxis(ax4, 'left');
    mark_regions(ax4, ridge_regions);
    grid(ax4, 'on');
    xlabel(ax4, 'Time (s)');
    title(ax4, sprintf('Ridge tracking in %.1f-%.1f Hz', ridge_search_band_hz(1), ridge_search_band_hz(2)));

    fprintf('\nRamp file: %s\n', ramp_file);
    fprintf('  Ridge peak time = %.3f s\n', peak_time);
    fprintf('  Ridge peak frequency = %.3f Hz\n', peak_freq);
    if isempty(ridge_regions)
        fprintf('  No ridge regions exceeded the threshold.\n');
    else
        for k = 1:numel(ridge_regions)
            fprintf('  Resonance region %d = [%.3f, %.3f] s, duration = %.3f s, peak ridge amp = %.3f, peak ridge freq = %.3f Hz, group = G%d (%.3f Hz)\n', ...
                k, ridge_regions(k).t_start, ridge_regions(k).t_end, ridge_regions(k).duration, ...
                ridge_regions(k).peak_amp, ridge_regions(k).peak_freq, ...
                ridge_regions(k).group_id, ridge_regions(k).group_center_freq);
        end
        fprintf('  Frequency groups:\n');
        for g = 1:numel(ridge_groups)
            fprintf('    G%d: center = %.3f Hz, n = %d, total duration = %.3f s\n', ...
                ridge_groups(g).group_id, ridge_groups(g).center_freq, ...
                numel(ridge_groups(g).region_indices), ridge_groups(g).total_duration);
        end
    end
end

fprintf('\n建议：\n');
fprintf('1. 用 20250526_910 的 BTT 低速数据做静态标定。\n');
fprintf('2. 若处理稳态 3150 工况，可直接把 3150 作为共振态辨识案例。\n');
fprintf('3. 若处理 2500-3500 变速工况，优先截取上面输出的推荐共振时间窗做辨识。\n');


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
if ~isfield(S, 'Datas') || size(S.Datas, 2) < 2
    error('File %s does not contain Datas(:,2).', filepath);
end
t = S.Datas(:,1);
x = S.Datas(:,2);

if isfield(S, 'SampleFrequency')
    Fs = str2double(S.SampleFrequency);
else
    dt = median(diff(t));
    Fs = 1 / dt;
end
end


function [f, amp] = compute_fft_spectrum(x, Fs)
x = detrend(x);
N = min(numel(x), 262144);
x = x(1:N) .* hann(N);
Y = fft(x);
P2 = abs(Y / N);
amp = P2(1:floor(N/2) + 1);
amp(2:end-1) = 2 * amp(2:end-1);
f = Fs * (0:floor(N/2)) / N;
end


function label = get_folder_label(folderpath)
[~, label] = fileparts(folderpath);
end


function regions = build_ridge_regions(T, ridge_amp_detect, ridge_amp_raw, ridge_freq, ...
    ridge_amp_threshold, merge_gap_sec, min_duration_sec)
mask = (ridge_amp_detect(:).' >= ridge_amp_threshold);

dt = median(diff(T));
gap_frames = max(0, round(merge_gap_sec / max(dt, eps)));
if gap_frames > 0
    mask = close_small_gaps(mask, gap_frames);
end

d = diff([false, mask, false]);
start_idx = find(d == 1);
end_idx = find(d == -1) - 1;

regions = struct('t_start', {}, 't_end', {}, 'duration', {}, 'peak_amp', {}, 'peak_freq', {}, ...
    'index_range', {}, 'group_id', {}, 'group_center_freq', {});
count = 0;
for k = 1:numel(start_idx)
    idx = start_idx(k):end_idx(k);
    t_start = T(idx(1));
    t_end = T(idx(end));
    duration = t_end - t_start;
    if duration < min_duration_sec
        continue;
    end
    count = count + 1;
    [peak_amp, local_idx] = max(ridge_amp_raw(idx));
    peak_freq = ridge_freq(idx(local_idx));
    regions(count).t_start = t_start;
    regions(count).t_end = t_end;
    regions(count).duration = duration;
    regions(count).peak_amp = peak_amp;
    regions(count).peak_freq = peak_freq;
    regions(count).index_range = idx;
    regions(count).group_id = NaN;
    regions(count).group_center_freq = NaN;
end
end


function [regions, groups] = cluster_ridge_regions_by_frequency(regions, gap_hz)
groups = struct('group_id', {}, 'center_freq', {}, 'region_indices', {}, 'total_duration', {});
if isempty(regions)
    return;
end

freqs = [regions.peak_freq];
[freqs_sorted, order] = sort(freqs);
raw_group_ids = zeros(size(order));
centers = [];
group_count = 0;

for i = 1:numel(freqs_sorted)
    if group_count == 0
        group_count = 1;
        centers(group_count) = freqs_sorted(i); %#ok<AGROW>
        raw_group_ids(i) = group_count;
        continue;
    end

    if abs(freqs_sorted(i) - centers(group_count)) <= gap_hz
        raw_group_ids(i) = group_count;
        idx_prev = order(raw_group_ids == group_count);
        centers(group_count) = mean([regions(idx_prev).peak_freq]);
    else
        group_count = group_count + 1;
        centers(group_count) = freqs_sorted(i); %#ok<AGROW>
        raw_group_ids(i) = group_count;
    end
end

for i = 1:numel(order)
    regions(order(i)).group_id = raw_group_ids(i);
end

for gid = 1:group_count
    idx = find([regions.group_id] == gid);
    groups(gid).group_id = gid;
    groups(gid).center_freq = mean([regions(idx).peak_freq]);
    groups(gid).region_indices = idx;
    groups(gid).total_duration = sum([regions(idx).duration]);
end

[~, sort_idx] = sort([groups.total_duration], 'descend');
groups = groups(sort_idx);

map_old_to_new = zeros(1, group_count);
for new_gid = 1:numel(groups)
    old_gid = groups(new_gid).group_id;
    map_old_to_new(old_gid) = new_gid;
    groups(new_gid).group_id = new_gid;
end

for k = 1:numel(regions)
    regions(k).group_id = map_old_to_new(regions(k).group_id);
    regions(k).group_center_freq = groups(regions(k).group_id).center_freq;
end

for g = 1:numel(groups)
    groups(g).region_indices = find([regions.group_id] == g);
    groups(g).center_freq = mean([regions(groups(g).region_indices).peak_freq]);
    groups(g).total_duration = sum([regions(groups(g).region_indices).duration]);
end
end


function mask = close_small_gaps(mask, max_gap_frames)
mask = mask(:).';
if ~any(mask)
    return;
end

d = diff([false, mask, false]);
zero_start = find(d == -1);
zero_end = find(d == 1) - 1;
for k = 1:numel(zero_start)
    gap_len = zero_end(k) - zero_start(k) + 1;
    if gap_len <= max_gap_frames
        left_ok = zero_start(k) > 1;
        right_ok = zero_end(k) < numel(mask);
        if left_ok && right_ok
            mask(zero_start(k):zero_end(k)) = true;
        end
    end
end
end


function mark_regions(ax, regions)
for k = 1:numel(regions)
    xline(ax, regions(k).t_start, '--r');
    xline(ax, regions(k).t_end, '--r');
end
end
