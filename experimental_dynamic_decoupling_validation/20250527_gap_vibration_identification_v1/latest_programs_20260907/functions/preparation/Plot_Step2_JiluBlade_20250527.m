function Plot_Step2_JiluBlade_20250527(case_name)
%PLOT_STEP2_JILUBLADE_20250527 Visualize Step2 jilublade extraction result.

cfg = Get_20250527_BTT_Config();
if nargin < 1 || isempty(case_name)
    case_name = '20250526_3150';
end

output_dir = fullfile(cfg.output_root, case_name);
if ~exist(output_dir, 'dir')
    Step2_Extract_JiluBlade_20250527(case_name);
end

case_data = Extract_BTT_Features_20250527(case_name, cfg);
sensor_ids = cfg.sensor_ids;
blades_num = cfg.blades_num;
colors = lines(blades_num);

figure('Name', sprintf('20250527 Step2 RPM - %s', case_name), ...
    'Color', 'w', 'Position', [80, 80, 1200, 400], 'NumberTitle', 'off');
hold on;
grid on;
box on;
plot(case_data.omega_time_s, case_data.omega_rpm, '-', 'LineWidth', 1.2, ...
    'Color', [0.1 0.1 0.1]);
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('OPR-derived RPM trend: %s', case_name), 'Interpreter', 'none');

fig2 = figure('Name', sprintf('20250527 Step2 Blade Sequence - %s', case_name), ...
    'Color', 'w', 'Position', [120, 120, 1450, 900], 'NumberTitle', 'off');
tiledlayout(fig2, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    filepath = fullfile(output_dir, sprintf('jilublade_probe%d.mat', sid));
    if ~isfile(filepath)
        continue;
    end

    loaded = load(filepath, 'jilublade');
    jilublade = loaded.jilublade;
    n_show = min(size(jilublade, 1), 120);

    nexttile;
    hold on;
    grid on;
    box on;
    for blade_id = 1:blades_num
        mask = jilublade(1:n_show, 4) == blade_id;
        scatter(jilublade(mask, 3), jilublade(mask, 4), 24, colors(blade_id, :), ...
            'filled', 'DisplayName', sprintf('B%d', blade_id));
    end
    xlabel('Arrival time (s)');
    ylabel('Blade ID');
    ylim([0.5, blades_num + 0.5]);
    title(sprintf('CH%d first %d labeled arrivals', sid, n_show));
    legend('Location', 'eastoutside');
end

fig3 = figure('Name', sprintf('20250527 Step2 Pulse Feature - %s', case_name), ...
    'Color', 'w', 'Position', [160, 160, 1450, 900], 'NumberTitle', 'off');
tiledlayout(fig3, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    arr = case_data.channels(sid).arrival_times(:);
    if isempty(feats)
        continue;
    end

    nexttile;
    hold on;
    grid on;
    box on;
    plot(arr, feats, '-', 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], ...
        'DisplayName', 'Pulse feature');
    [pks, locs] = findpeaks(feats, 'NPeaks', min(15, numel(feats)), 'SortStr', 'descend');
    scatter(arr(locs), pks, 28, [0.85 0.33 0.10], 'filled', 'DisplayName', 'Local strong peaks');
    xlabel('Arrival time (s)');
    ylabel('Peak feature');
    title(sprintf('CH%d pulse feature trend', sid));
    legend('Location', 'best');
end

fprintf('Opened Step2 visualization figures for %s.\n', case_name);
end
