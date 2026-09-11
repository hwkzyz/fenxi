function Plot_Step1_Sensor_Config_20250527()
%PLOT_STEP1_SENSOR_CONFIG_20250527 Visualize low-speed fingerprint config.

cfg = Get_20250527_BTT_Config();
config_path = fullfile(cfg.reference_output_dir, 'Sensor_Config_20250527.mat');
if ~isfile(config_path)
    Step1_Build_Sensor_Config_20250527();
end

loaded = load(config_path, 'Sensor_Config');
Sensor_Config = loaded.Sensor_Config;
case_data = Extract_BTT_Features_20250527(cfg.low_speed_case, cfg);

sensor_ids = Sensor_Config.Sensor_IDs;
blades_num = Sensor_Config.Blades_Num;
n_show = min(4 * blades_num, min_pulse_count_local(case_data, sensor_ids));

fig1 = figure('Name', '20250527 Step1 Fingerprint Alignment', ...
    'Color', 'w', 'Position', [80, 80, 1400, 900], 'NumberTitle', 'off');
tiledlayout(fig1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

anchor_sid = cfg.reference_sensor_id;
ref_fp = Sensor_Config.Fingerprints(anchor_sid);
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    start_idx = Sensor_Config.Target_Indices(sid);
    pulse_idx = 1:n_show;

    nexttile;
    hold on;
    grid on;
    box on;
    plot(pulse_idx, feats(pulse_idx), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, ...
        'DisplayName', 'Raw pulse feature');
    sel_idx = start_idx:(start_idx + blades_num - 1);
    plot(sel_idx, feats(sel_idx), 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'Selected fingerprint');
    if sid ~= anchor_sid
        plot(sel_idx, ref_fp(:), 's--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Reference CH1');
    end
    xlabel('Pulse index');
    ylabel('Peak feature');
    title(sprintf('CH%d fingerprint alignment', sid));
    legend('Location', 'best');
end

fig2 = figure('Name', '20250527 Step1 Standard Relative Angles', ...
    'Color', 'w', 'Position', [120, 120, 1200, 520], 'NumberTitle', 'off');
tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
angle_mat = Sensor_Config.Standard_Relative_Angles(sensor_ids, :);
imagesc(angle_mat);
axis tight;
grid on;
colorbar;
set(gca, 'XTick', 1:blades_num, 'YTick', 1:numel(sensor_ids), ...
    'YTickLabel', compose('CH%d', sensor_ids));
xlabel('Blade ID');
ylabel('Sensor');
title('Standard relative angles (deg)');

nexttile;
hold on;
grid on;
box on;
colors = lines(numel(sensor_ids));
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    plot(1:blades_num, angle_mat(i, :), 'o-', 'LineWidth', 1.5, ...
        'Color', colors(i, :), 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('CH%d', sid));
end
xlabel('Blade ID');
ylabel('Angle relative to previous OPR (deg)');
title('Per-sensor blade angle curves');
legend('Location', 'best');

figure('Name', '20250527 Step1 Low-Speed RPM', ...
    'Color', 'w', 'Position', [160, 160, 1200, 420], 'NumberTitle', 'off');
hold on;
grid on;
box on;
plot(case_data.omega_time_s, case_data.omega_rpm, '-', 'LineWidth', 1.2, ...
    'Color', [0.2 0.2 0.2]);
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('Low-speed RPM trend: %s', cfg.low_speed_case), 'Interpreter', 'none');

fprintf('Opened Step1 visualization figures for %s.\n', cfg.low_speed_case);
end


function n = min_pulse_count_local(case_data, sensor_ids)
counts = zeros(numel(sensor_ids), 1);
for i = 1:numel(sensor_ids)
    counts(i) = numel(case_data.channels(sensor_ids(i)).peak_features);
end
n = min(counts);
end
