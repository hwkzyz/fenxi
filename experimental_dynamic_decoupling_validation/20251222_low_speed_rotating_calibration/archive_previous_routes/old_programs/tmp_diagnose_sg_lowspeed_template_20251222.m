%% tmp_diagnose_sg_lowspeed_template_20251222
% Diagnostic comparison only. This script does not modify the identification
% pipeline; it compares the low-speed rotating template against the original
% SG calibration on the same S123/W08 dynamic data.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(route_dir, 'output', 'diagnostics', 'sg_vs_lowspeed_template');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

sg_dir = ['E:\0小论文+程序\0博士期间小论文+程序\7超高斯模型-权重-瞬态\程序\' ...
    '直叶片验证\实验验证\超高斯\20251222适配-3号传感间隙变化-改'];
sg_result_file = fullfile(sg_dir, 'Result_single_sync_20251222_B1_S123_1000_2500_3500_Start50p0s.mat');
ls_allpulse_file = fullfile(route_dir, 'output', 'identification', ...
    'Result_TemplateOnlySGLikeSingleSync_B1_S123_20251222.allpulse_balanced.bak.mat');
template_file = fullfile(route_dir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_20251222.mat');

if ~isfile(sg_result_file)
    error('Missing SG result file: %s', sg_result_file);
end
if ~isfile(ls_allpulse_file)
    error('Missing low-speed all-pulse backup result file: %s', ls_allpulse_file);
end
if ~isfile(template_file)
    error('Missing low-speed template file: %s', template_file);
end

loaded_sg = load(sg_result_file, 'Result_Struct');
SG = loaded_sg.Result_Struct;
loaded_ls = load(ls_allpulse_file, 'Result');
LS = loaded_ls.Result;
loaded_tpl = load(template_file, 'Template');
Template = loaded_tpl.Template;

window_id = 8;
sg_bundle = SG.BestWindow.bundle;
if SG.BestWindow.window_id ~= window_id
    error('Expected SG best window %d, got %d.', window_id, SG.BestWindow.window_id);
end
ls_bundle = LS.WindowResult(window_id).bundle;
ls_ct = LS.WindowResult(window_id).CandidateTable;
ls_eo14_row = ls_ct(ls_ct.EO == 14, :);
if isempty(ls_eo14_row)
    error('EO14 row not found in low-speed candidate table for W%d.', window_id);
end

fig = figure('Name', 'SG vs low-speed template diagnostic W08', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 20]);
tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'ls_static_plain_rmse_v', NaN, ...
    'sg_static_plain_rmse_v', NaN, ...
    'ls_eo14_plain_rmse_v', NaN, ...
    'sg_eo14_plain_rmse_v', NaN, ...
    'ls_eo14_bias_v', NaN, ...
    'sg_eo14_bias_v', NaN, ...
    'ls_point_count', NaN, ...
    'sg_point_count', NaN), 3, 1);

for is = 1:3
    sid = is;
    tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    sg_calib = load(fullfile(sg_dir, sprintf('SGCalib_20251222_B1_S%d.mat', sid)), 'Calib');
    cal = sg_calib.Calib.Sensor;

    m_ls = ls_bundle.sensor_index == is;
    x_abs_ls = ls_bundle.X(m_ls) + tpl.xc;
    v_ls = ls_bundle.V(m_ls);
    x_rel_ls = x_abs_ls - tpl.xc;
    v_tpl_at_ls = interp1(tpl.x_grid(:), tpl.v_grid(:), x_rel_ls, 'pchip', NaN);
    sg_static_at_ls = sg_forward_local(x_abs_ls, cal);

    u_ls14 = ls_eo14_row.A .* sin(14 .* ls_bundle.Theta(m_ls) + ls_eo14_row.phi) + ls_eo14_row.d0;
    v_ls14 = interp1(tpl.x_grid(:), tpl.v_grid(:), x_rel_ls - u_ls14, 'pchip', NaN);

    m_sg = sg_bundle.S == sid;
    x_abs_sg = sg_bundle.X(m_sg);
    v_sg = sg_bundle.V(m_sg);
    sg_static_at_sg = sg_forward_local(x_abs_sg, cal);
    ls_static_at_sg = interp1(tpl.x_grid(:), tpl.v_grid(:), x_abs_sg - tpl.xc, 'pchip', NaN);
    v_sg14 = SG.BestWindow.V_pred(m_sg);

    r_ls_static = v_ls - v_tpl_at_ls;
    r_sg_static = v_sg - sg_static_at_sg;
    r_ls14 = v_ls - v_ls14;
    r_sg14 = v_sg - v_sg14;

    rows(is).sensor_id = sid;
    rows(is).ls_static_plain_rmse_v = sqrt(mean(r_ls_static.^2, 'omitnan'));
    rows(is).sg_static_plain_rmse_v = sqrt(mean(r_sg_static.^2, 'omitnan'));
    rows(is).ls_eo14_plain_rmse_v = sqrt(mean(r_ls14.^2, 'omitnan'));
    rows(is).sg_eo14_plain_rmse_v = sqrt(mean(r_sg14.^2, 'omitnan'));
    rows(is).ls_eo14_bias_v = mean(r_ls14, 'omitnan');
    rows(is).sg_eo14_bias_v = mean(r_sg14, 'omitnan');
    rows(is).ls_point_count = nnz(m_ls);
    rows(is).sg_point_count = nnz(m_sg);

    x_plot = linspace(max([min(tpl.x_grid + tpl.xc), min(x_abs_ls), min(x_abs_sg)]), ...
        min([max(tpl.x_grid + tpl.xc), max(x_abs_ls), max(x_abs_sg)]), 500).';
    v_low_plot = interp1(tpl.x_grid(:), tpl.v_grid(:), x_plot - tpl.xc, 'pchip', NaN);
    v_sg_plot = sg_forward_local(x_plot, cal);

    nexttile;
    plot(x_abs_ls, v_ls, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 5); hold on;
    plot(x_plot, v_low_plot, 'r-', 'LineWidth', 1.2);
    plot(x_plot, v_sg_plot, 'b--', 'LineWidth', 1.2);
    xlabel('absolute x (mm)');
    ylabel('V');
    title(sprintf('S%d dynamic points + static curves', sid));
    legend({'dynamic W08', 'low-speed T_{low}', 'SGCalib'}, 'Location', 'best');
    grid on;

    nexttile;
    histogram(r_ls14, 60, 'FaceColor', [0.85 0.25 0.20], 'FaceAlpha', 0.55, 'EdgeColor', 'none'); hold on;
    histogram(r_sg14, 60, 'FaceColor', [0.15 0.35 0.85], 'FaceAlpha', 0.55, 'EdgeColor', 'none');
    xlabel('EO14 residual V_{meas}-V_{pred} (V)');
    ylabel('count');
    title(sprintf('S%d EO14 residuals', sid));
    legend({'low-speed', 'SG'}, 'Location', 'best');
    grid on;
end

summary_table = struct2table(rows);
disp(summary_table);

csv_file = fullfile(out_dir, 'SG_vs_LowSpeed_Template_W08_S123_summary.csv');
writetable(summary_table, csv_file);

fig_file = fullfile(out_dir, 'SG_vs_LowSpeed_Template_W08_S123.png');
exportgraphics(fig, fig_file, 'Resolution', 220);

fprintf('Saved summary: %s\n', csv_file);
fprintf('Saved figure:  %s\n', fig_file);

function v = sg_forward_local(x_abs, cal)
v = cal.B .* exp(-abs((x_abs - cal.xc) ./ cal.w) .^ cal.n) + cal.baseline;
end
