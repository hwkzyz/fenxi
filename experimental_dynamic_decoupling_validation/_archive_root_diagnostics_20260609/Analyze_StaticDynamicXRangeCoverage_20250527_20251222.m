%% Analyze_StaticDynamicXRangeCoverage_20250527_20251222
% Compare static calibration x-domain with dynamic waveform query ranges.

clear; clc; close all;

root_dir = fileparts(mfilename('fullpath'));
case_list = struct([]);
case_list(1).tag = '20250527';
case_list(1).folder = fullfile(root_dir, '20250527_low_speed_rotating_calibration');
case_list(1).template_file = fullfile(case_list(1).folder, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S136_GradientXRange030_20250527.mat');
case_list(1).result_file = fullfile(case_list(1).folder, 'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_GradientXRange030_20250527.mat');

case_list(2).tag = '20251222';
case_list(2).folder = fullfile(root_dir, '20251222_low_speed_rotating_calibration');
case_list(2).template_file = fullfile(case_list(2).folder, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
case_list(2).result_file = fullfile(case_list(2).folder, 'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_GradientXRange030_20251222.mat');

rows = table();
for ic = 1:numel(case_list)
    rows = [rows; audit_case_local(case_list(ic))]; %#ok<AGROW>
end

output_csv = fullfile(root_dir, 'StaticDynamicXRangeCoverage_Audit.csv');
writetable(rows, output_csv);

figure_dir = fullfile(root_dir, 'figures');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end
plot_coverage_summary_local(rows, fullfile(figure_dir, 'StaticDynamicXRangeCoverage_Audit.png'));

fprintf('Saved audit table: %s\n', output_csv);
fprintf('Saved audit figure: %s\n', fullfile(figure_dir, 'StaticDynamicXRangeCoverage_Audit.png'));
print_case_summary_local(rows);

function rows = audit_case_local(case_cfg)
loaded_template = load(case_cfg.template_file, 'Template');
loaded_result = load(case_cfg.result_file, 'Result');
Template = loaded_template.Template;
Result = loaded_result.Result;

rows = table();
for iw = 1:numel(Result.WindowResult)
    bundle = Result.WindowResult(iw).bundle;
    fit = Result.WindowResult(iw).Result;
    x_query = bundle.X - fit.dx_c_id - fit.u_est;

    for is = 1:numel(bundle.sensor_ids)
        sid = bundle.sensor_ids(is);
        sensor_mask = bundle.sensor_index == is;
        Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
        x_dyn = bundle.X(sensor_mask);
        xq = x_query(sensor_mask);
        weights = bundle.W(sensor_mask);

        row = table(string(case_cfg.tag), iw, sid, ...
            Tpl.x_domain(1), Tpl.x_domain(2), diff(Tpl.x_domain), ...
            min(x_dyn), max(x_dyn), prctile(x_dyn, 1), prctile(x_dyn, 99), ...
            min(xq), max(xq), prctile(xq, 1), prctile(xq, 99), ...
            min(xq) - Tpl.x_domain(1), Tpl.x_domain(2) - max(xq), ...
            prctile(xq, 1) - Tpl.x_domain(1), Tpl.x_domain(2) - prctile(xq, 99), ...
            mean(weights, 'omitnan'), nnz(xq < Tpl.x_domain(1) | xq > Tpl.x_domain(2)) / max(nnz(sensor_mask), 1), ...
            fit.EO_id, fit.A_id, fit.weighted_voltage_rmse, ...
            'VariableNames', {'CaseTag','Window','SensorID', ...
            'StaticLeft','StaticRight','StaticWidth', ...
            'DynMin','DynMax','DynP01','DynP99', ...
            'QueryMin','QueryMax','QueryP01','QueryP99', ...
            'MarginMinLeft','MarginMinRight','MarginP01Left','MarginP99Right', ...
            'MeanWeight','ClampFrac','EO','AmplitudeMM','WeightedRMSE'});
        rows = [rows; row]; %#ok<AGROW>
    end
end
end

function plot_coverage_summary_local(rows, output_png)
fig = figure('Color', 'w', 'Position', [120 120 1250 780]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
cases = unique(rows.CaseTag, 'stable');
colors = lines(numel(cases));
for ic = 1:numel(cases)
    mask_case = rows.CaseTag == cases(ic);
    sensors = unique(rows.SensorID(mask_case), 'stable');
    for is = 1:numel(sensors)
        mask = mask_case & rows.SensorID == sensors(is);
        xpos = ic + 0.12 * (is - (numel(sensors)+1)/2);
        static_left = mean(rows.StaticLeft(mask), 'omitnan');
        static_right = mean(rows.StaticRight(mask), 'omitnan');
        query_p01 = min(rows.QueryP01(mask));
        query_p99 = max(rows.QueryP99(mask));
        query_min = min(rows.QueryMin(mask));
        query_max = max(rows.QueryMax(mask));
        plot([static_left static_right], [xpos xpos], '-', 'Color', colors(ic,:), 'LineWidth', 5);
        plot([query_p01 query_p99], [xpos xpos], 'k-', 'LineWidth', 2);
        plot([query_min query_max], [xpos xpos], 'r:', 'LineWidth', 1.2);
        text(static_right + 0.05, xpos, sprintf('%s-S%d', cases(ic), sensors(is)), 'FontSize', 9);
    end
end
grid on;
xlabel('x (mm)');
ylabel('Case / sensor');
title('Static x-domain vs dynamic query x-range');
legend({'static x-domain','query P01-P99','query min-max'}, 'Location', 'best');

nexttile;
hold on;
for ic = 1:numel(cases)
    mask_case = rows.CaseTag == cases(ic);
    sensors = unique(rows.SensorID(mask_case), 'stable');
    for is = 1:numel(sensors)
        mask = mask_case & rows.SensorID == sensors(is);
        label = sprintf('%s-S%d', cases(ic), sensors(is));
        scatter(repmat(ic + 0.12 * (is - (numel(sensors)+1)/2), nnz(mask), 1), ...
            rows.ClampFrac(mask) * 100, 18, colors(ic,:), 'filled', 'DisplayName', label);
    end
end
grid on;
ylabel('Clamp fraction (%)');
xticks(1:numel(cases));
xticklabels(cases);
title('Per-window query overshoot/clamp fraction');

exportgraphics(fig, output_png, 'Resolution', 300);
end

function print_case_summary_local(rows)
cases = unique(rows.CaseTag, 'stable');
for ic = 1:numel(cases)
    M = rows(rows.CaseTag == cases(ic), :);
    fprintf('\n%s\n', cases(ic));
    fprintf('Static width mean %.3f, range [%.3f %.3f] mm\n', ...
        mean(M.StaticWidth), min(M.StaticWidth), max(M.StaticWidth));
    fprintf('Dynamic raw P01-P99 global [%.3f %.3f] mm\n', min(M.DynP01), max(M.DynP99));
    fprintf('Query P01-P99 global [%.3f %.3f] mm; query min-max [%.3f %.3f] mm\n', ...
        min(M.QueryP01), max(M.QueryP99), min(M.QueryMin), max(M.QueryMax));
    fprintf('Robust margins: left %.3f mm, right %.3f mm; mean clamp %.4f%%, max clamp %.4f%%\n', ...
        min(M.MarginP01Left), min(M.MarginP99Right), 100 * mean(M.ClampFrac), 100 * max(M.ClampFrac));
end
end
