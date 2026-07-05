function Scan_StartTime_GradientXRange030_20251222()
%% Scan_StartTime_GradientXRange030_20251222
% Check whether shifting the dynamic start time removes the first-window EO
% outlier while keeping the same ratio=0.30 static template.

clc; close all;

route_dir = fileparts(mfilename('fullpath'));
start_times = parse_numeric_list_env_local('STEP02_START_TIME_SCAN', [50.0, 50.1, 50.2, 50.3, 50.5]);

template_file = fullfile(route_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
if exist(template_file, 'file') ~= 2
    run(fullfile(route_dir, 'Step01_Main_Build_GradientXRange030_Template_20251222.m'));
end

rows = table();
for st = start_times(:).'
    tag = sprintf('Start%05dms', round(st * 1000));
    fprintf('\n=== Start-time candidate %.3f s (%s) ===\n', st, tag);

    setenv('STEP02_ANALYSIS_START_TIME', sprintf('%.8g', st));
    setenv('STEP02_DYNAMIC_SUFFIX', tag);
    run_script_in_base_local(route_dir, 'Step02_Build_Dynamic_Map_20251222.m');

    dynamic_file = fullfile(route_dir, 'output', 'dynamic_maps', ...
        sprintf('DynamicMap_B1_S123_SlidingWindows_%s_20251222.mat', tag));
    setenv('STEP03_TEMPLATE_FILE', template_file);
    setenv('STEP03_DYNAMIC_MAP_FILE', dynamic_file);
    setenv('STEP03_RESULT_SUFFIX', ['GradientXRange030_', tag]);
    setenv('STEP03_ANALYSIS_SENSORS', '1 2 3');
    setenv('STEP03_MAIN_TOP_K_EO', '3');
    setenv('STEP03D_PULSE_MODE', 'all');
    setenv('STEP03D_DOMAIN_SELECTION_MODE', 'soft');
    setenv('STEP03D_DOMAIN_SOFT_MARGIN_MM', '0.2');
    setenv('STEP03D_OVERSHOOT_PENALTY_WEIGHT', '100');
    run_script_in_base_local(route_dir, 'Step03_Main_VPTop3SynchronousWaveform_20251222.m');

    result_file = fullfile(route_dir, 'output', 'identification', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_GradientXRange030_%s_20251222.mat', tag));
    loaded = load(result_file, 'Result');
    T = loaded.Result.Trend;
    row = table();
    row.StartTimeSec = st;
    row.Tag = string(tag);
    row.Window1EO = T.EO_id(1);
    row.Window1RMSE = T.weighted_voltage_rmse(1);
    row.DominantEO = mode(T.EO_id);
    row.EOConsistency = mean(T.EO_id == row.DominantEO, 'omitnan');
    row.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
    row.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
    row.ResultFile = string(result_file);
    rows = [rows; row]; %#ok<AGROW>
end

clear_scan_env_local();
summary_file = fullfile(route_dir, 'output', 'identification', ...
    'Scan_StartTime_GradientXRange030_20251222.csv');
writetable(rows, summary_file);
disp(rows(:, {'StartTimeSec','Window1EO','DominantEO','EOConsistency','MeanRMSE','MedianRMSE'}));

figure_dir = fullfile(route_dir, 'output', 'figures', 'start_time_scan');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end
fig = figure('Name', 'Start-time scan', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 16, 8]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(rows.StartTimeSec, rows.Window1EO, 'o-', 'LineWidth', 1.2); hold on;
plot(rows.StartTimeSec, rows.DominantEO, 's--', 'LineWidth', 1.2);
ylabel('EO');
legend({'Window 1', 'Dominant'}, 'Location', 'best');
grid on; box on; set(gca, 'FontName', 'Times New Roman', 'TickDir', 'in');
nexttile;
yyaxis left;
plot(rows.StartTimeSec, rows.EOConsistency, 'o-', 'LineWidth', 1.2);
ylabel('EO consistency');
yyaxis right;
plot(rows.StartTimeSec, rows.MeanRMSE, 's-', 'LineWidth', 1.2);
ylabel('Mean RMSE (V)');
xlabel('Analysis start time (s)');
grid on; box on; set(gca, 'FontName', 'Times New Roman', 'TickDir', 'in');
exportgraphics(fig, fullfile(figure_dir, 'Scan_StartTime_GradientXRange030_20251222.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, 'Scan_StartTime_GradientXRange030_20251222.pdf'), 'ContentType', 'vector');
fprintf('Saved summary: %s\n', summary_file);
end

function vals = parse_numeric_list_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    vals = default_value;
    return;
end
vals = sscanf(raw, '%f').';
if isempty(vals)
    error('%s must contain numeric values separated by spaces.', name);
end
end

function run_script_in_base_local(route_dir, script_name)
evalin('base', sprintf("run('%s');", fullfile(route_dir, script_name)));
end

function clear_scan_env_local()
vars = {'STEP02_ANALYSIS_START_TIME','STEP02_DYNAMIC_SUFFIX','STEP03_TEMPLATE_FILE', ...
    'STEP03_DYNAMIC_MAP_FILE','STEP03_RESULT_SUFFIX','STEP03_ANALYSIS_SENSORS', ...
    'STEP03_MAIN_TOP_K_EO','STEP03D_PULSE_MODE','STEP03D_DOMAIN_SELECTION_MODE', ...
    'STEP03D_DOMAIN_SOFT_MARGIN_MM','STEP03D_OVERSHOOT_PENALTY_WEIGHT'};
for i = 1:numel(vars)
    setenv(vars{i}, '');
end
end
