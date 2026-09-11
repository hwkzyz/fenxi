clear; clc; close all;

sensor_data_dir = 'E:\试验数据\20250527\试验20250527\20250526_3150';
strain_data_dir = 'E:\试验数据\20250527\应变片数据20250527\20250527_3150';
strain_channel_tag = 'AI1-03';
sensor_channel_ids = [1, 3, 6];

sensor_fs = 5e6;
feature_dt = 0.01;
sensor_bin_sec = 0.01;
strain_smooth_sec = 0.05;

coarse_step_sec = 0.05;
fine_half_range_sec = 0.30;
fine_step_sec = 0.002;
min_overlap_sec = 5.0;

fprintf('=== 20250527 传感器-应变片时间对齐 ===\n');
fprintf('Sensor channels: [%s]\n', num2str(sensor_channel_ids));
fprintf('Strain channel: %s\n', strain_channel_tag);

sensor_feature_list = cell(size(sensor_channel_ids));
sensor_time_list = cell(size(sensor_channel_ids));

for i = 1:numel(sensor_channel_ids)
    ch = sensor_channel_ids(i);
    [t_feat, feat] = build_sensor_feature(sensor_data_dir, ch, sensor_fs, sensor_bin_sec);
    sensor_time_list{i} = t_feat;
    sensor_feature_list{i} = feat;
end

[t_sensor, sensor_feature] = combine_sensor_features(sensor_time_list, sensor_feature_list, feature_dt);

strain_file = pick_latest_strain_file(strain_data_dir, strain_channel_tag);
[t_strain, v_strain, Fs_strain] = load_strain_signal(strain_file);
[t_strain_feat, strain_feature] = build_strain_feature(t_strain, v_strain, Fs_strain, feature_dt, strain_smooth_sec);

fprintf('Sensor feature span: [%.3f, %.3f] s\n', t_sensor(1), t_sensor(end));
fprintf('Strain feature span: [%.3f, %.3f] s\n', t_strain_feat(1), t_strain_feat(end));

offset_min = t_sensor(1) - t_strain_feat(end);
offset_max = t_sensor(end) - t_strain_feat(1);

[coarse_best_offset, coarse_best_score, coarse_offsets, coarse_scores] = ...
    search_best_offset(t_sensor, sensor_feature, t_strain_feat, strain_feature, ...
    offset_min, offset_max, coarse_step_sec, min_overlap_sec);

[fine_best_offset, fine_best_score, fine_offsets, fine_scores] = ...
    search_best_offset(t_sensor, sensor_feature, t_strain_feat, strain_feature, ...
    coarse_best_offset - fine_half_range_sec, coarse_best_offset + fine_half_range_sec, ...
    fine_step_sec, min_overlap_sec);

fprintf('Coarse best offset = %.4f s, score = %.4f\n', coarse_best_offset, coarse_best_score);
fprintf('Final best offset  = %.4f s, score = %.4f\n', fine_best_offset, fine_best_score);

t_strain_aligned = t_strain_feat + fine_best_offset;
strain_feature_aligned = strain_feature;

sensor_on_sensor = interp1(t_sensor, sensor_feature, t_sensor, 'linear', NaN);
strain_on_sensor = interp1(t_strain_aligned, strain_feature_aligned, t_sensor, 'linear', NaN);
valid = isfinite(sensor_on_sensor) & isfinite(strain_on_sensor);

figure('Name', '20250527 Alignment Search', 'Color', 'w', 'NumberTitle', 'off');

ax1 = subplot(3,1,1);
plot(ax1, coarse_offsets, coarse_scores, 'b', 'LineWidth', 1.2);
hold(ax1, 'on');
xline(ax1, coarse_best_offset, '--r', sprintf('Coarse %.3f s', coarse_best_offset));
grid(ax1, 'on');
xlabel(ax1, 'Offset applied to strain (s)');
ylabel(ax1, 'Correlation score');
title(ax1, 'Coarse alignment search');

ax2 = subplot(3,1,2);
plot(ax2, fine_offsets, fine_scores, 'm', 'LineWidth', 1.2);
hold(ax2, 'on');
xline(ax2, fine_best_offset, '--r', sprintf('Final %.3f s', fine_best_offset));
grid(ax2, 'on');
xlabel(ax2, 'Offset applied to strain (s)');
ylabel(ax2, 'Correlation score');
title(ax2, 'Fine alignment search');

ax3 = subplot(3,1,3);
plot(ax3, t_sensor, sensor_feature, 'b', 'LineWidth', 1.2, 'DisplayName', 'Sensor envelope feature');
hold(ax3, 'on');
plot(ax3, t_strain_aligned, strain_feature_aligned, 'Color', [0.1 0.6 0.2], 'LineWidth', 1.0, ...
    'DisplayName', 'Aligned strain envelope');
if any(valid)
    xlim(ax3, [t_sensor(find(valid,1,'first')) - 0.5, t_sensor(find(valid,1,'last')) + 0.5]);
end
grid(ax3, 'on');
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Normalized feature');
title(ax3, sprintf('Aligned features (offset %.4f s)', fine_best_offset));
legend(ax3, 'Location', 'best');

figure('Name', '20250527 Alignment Detail', 'Color', 'w', 'NumberTitle', 'off');

ax1 = subplot(2,1,1);
plot(ax1, t_sensor, sensor_feature, 'b', 'LineWidth', 1.2);
hold(ax1, 'on');
plot(ax1, t_strain_feat, strain_feature, 'Color', [0.3 0.7 0.3], 'LineWidth', 1.0);
grid(ax1, 'on');
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Normalized feature');
title(ax1, 'Before alignment');
legend(ax1, {'Sensor feature', 'Strain feature'}, 'Location', 'best');

ax2 = subplot(2,1,2);
plot(ax2, t_sensor, sensor_feature, 'b', 'LineWidth', 1.2);
hold(ax2, 'on');
plot(ax2, t_strain_aligned, strain_feature_aligned, 'Color', [0.1 0.6 0.2], 'LineWidth', 1.0);
grid(ax2, 'on');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Normalized feature');
title(ax2, sprintf('After alignment (offset %.4f s)', fine_best_offset));
legend(ax2, {'Sensor feature', 'Aligned strain feature'}, 'Location', 'best');

fprintf('\n建议下一步:\n');
fprintf('1. 先用这个 offset 把应变片时间轴平移到传感器时间轴。\n');
fprintf('2. 再把当前共振候选区时间窗同步映射到传感器/BTT 数据。\n');
fprintf('3. 若后续得到 jilublade 位移结果, 可复用同样流程再做一次更精细的位移级对齐。\n');


function [t_feat, feat] = build_sensor_feature(sensor_data_dir, channel_id, sensor_fs, sensor_bin_sec)
d = dir(fullfile(sensor_data_dir, sprintf('4-%d-*.mat', channel_id)));
if isempty(d)
    error('No raw sensor files found for channel %d.', channel_id);
end

ids = zeros(numel(d),1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    ids(i) = str2double(tok{1});
end
[~, order] = sort(ids);
d = d(order);

t_feat = [];
feat = [];
bin_pts = round(sensor_bin_sec * sensor_fs);

for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name));
    fn = fieldnames(S);
    x = S.(fn{1});
    t = x(:,1) / sensor_fs;
    v = x(:,2);
    baseline = median(v);
    v = v - baseline;

    usable = floor(numel(v) / bin_pts) * bin_pts;
    if usable < bin_pts
        continue;
    end

    v = reshape(v(1:usable), bin_pts, []);
    t = reshape(t(1:usable), bin_pts, []);

    % 用高分位绝对值提取脉冲强度，比直接 RMS 更能反映叶尖波形峰值变化
    local_feat = prctile(abs(v), 99.5, 1).';
    local_t = mean(t, 1).';

    t_feat = [t_feat; local_t]; %#ok<AGROW>
    feat = [feat; local_feat]; %#ok<AGROW>
end

[t_feat, idx] = unique(t_feat);
feat = feat(idx);
feat = normalize_feature(feat);
end


function [t_common, feat_common] = combine_sensor_features(sensor_time_list, sensor_feature_list, feature_dt)
t_start = max(cellfun(@(x) x(1), sensor_time_list));
t_end = min(cellfun(@(x) x(end), sensor_time_list));
t_common = (t_start:feature_dt:t_end).';

feat_matrix = nan(numel(t_common), numel(sensor_time_list));
for i = 1:numel(sensor_time_list)
    feat_matrix(:,i) = interp1(sensor_time_list{i}, sensor_feature_list{i}, t_common, 'linear', NaN);
end

feat_common = mean(feat_matrix, 2, 'omitnan');
feat_common = normalize_feature(feat_common);
end


function filepath = pick_latest_strain_file(strain_data_dir, channel_tag)
d = dir(fullfile(strain_data_dir, sprintf('%s_*.mat', channel_tag)));
if isempty(d)
    error('No strain files found for %s.', channel_tag);
end

names = {d.name};
is_new = contains(names, '20250527');
if any(is_new)
    d = d(is_new);
end

[~, idx] = max([d.datenum]);
filepath = fullfile(d(idx).folder, d(idx).name);
end


function [t, v, Fs] = load_strain_signal(filepath)
S = load(filepath);
t = S.Datas(:,1);
v = S.Datas(:,2);
if isfield(S, 'SampleFrequency')
    Fs = str2double(S.SampleFrequency);
else
    Fs = 1 / mean(diff(t));
end
end


function [t_feat, feat] = build_strain_feature(t, v, Fs, feature_dt, strain_smooth_sec)
v = detrend(v(:));
t = t(:);

bin_pts = round(feature_dt * Fs);
usable = floor(numel(v) / bin_pts) * bin_pts;
v = reshape(v(1:usable), bin_pts, []);
t = reshape(t(1:usable), bin_pts, []);

env = prctile(abs(v), 99, 1).';
t_feat = mean(t, 1).';

smooth_frames = max(1, round(strain_smooth_sec / feature_dt));
feat = movmean(env, smooth_frames);
feat = normalize_feature(feat);
end


function y = normalize_feature(x)
x = x(:);
x = x - median(x, 'omitnan');
scale = iqr(x);
if scale < eps
    scale = std(x, 'omitnan');
end
if scale < eps
    scale = 1;
end
y = x / scale;
end


function [best_offset, best_score, offsets, scores] = search_best_offset(t_sensor, sensor_feature, ...
    t_strain, strain_feature, offset_min, offset_max, step_sec, min_overlap_sec)
offsets = (offset_min:step_sec:offset_max).';
scores = nan(size(offsets));

for i = 1:numel(offsets)
    offset = offsets(i);
    strain_shifted = interp1(t_strain + offset, strain_feature, t_sensor, 'linear', NaN);
    valid = isfinite(sensor_feature) & isfinite(strain_shifted);
    if sum(valid) < 3
        continue;
    end
    overlap_duration = t_sensor(find(valid,1,'last')) - t_sensor(find(valid,1,'first'));
    if overlap_duration < min_overlap_sec
        continue;
    end

    a = sensor_feature(valid);
    b = strain_shifted(valid);
    a = a - mean(a, 'omitnan');
    b = b - mean(b, 'omitnan');
    denom = norm(a) * norm(b);
    if denom < eps
        continue;
    end
    scores(i) = (a' * b) / denom;
end

[best_score, idx] = max(scores);
best_offset = offsets(idx);
end
