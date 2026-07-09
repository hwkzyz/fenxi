%% Visualize_OriginalConfig_Waveform_Crops_20250527_20251222
% Read-only diagnostic figures for the original validated configurations.
% It compares static calibration waveform crops and dynamic vibration
% waveform crops for 20250527 and 20251222.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'diagnostics_original_config_visuals');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

caseList = build_case_list_local(rootDir);

for ic = 1:numel(caseList)
    C = caseList(ic);
    fprintf('\n=== %s ===\n', C.label);
    Tbase = load_template_local(C.baseTemplateFile);
    Tgrad = load_template_local(C.gradientTemplateFile);
    Dbase = load_dynamic_local(C.baseDynamicFile);
    Dgrad = load_dynamic_local(C.gradientDynamicFile);
    Rmain = load_result_local(C.mainResultFile);

    fprintf('Main result EO14 = %d/%d\n', nnz(Rmain.Trend.EO_id == 14), height(Rmain.Trend));
    fprintf('Main DynamicMap: %s\n', Rmain.DynamicMapFile);
    fprintf('Main pulse/domain: %s / %s\n', ...
        Rmain.MethodSettings.pulse_selection_mode, Rmain.MethodSettings.domain_selection_mode);

    fig1 = plot_static_template_crop_local(C, Tbase, Tgrad);
    save_figure_local(fig1, outDir, sprintf('%s_static_calibration_crop', C.tag));

    fig2 = plot_dynamic_crop_overlay_local(C, Tgrad, Dbase, Dgrad, Rmain);
    save_figure_local(fig2, outDir, sprintf('%s_dynamic_vibration_crop_compare', C.tag));
end

fprintf('\nSaved figures under:\n  %s\n', outDir);

function caseList = build_case_list_local(rootDir)
caseList = repmat(struct(), 1, 2);

caseList(1).tag = '20250527_B1_S136';
caseList(1).label = '20250527 original validated config';
caseList(1).sensorIds = [1 3 6];
caseList(1).mainWindow = 1;
caseList(1).baseTemplateFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S136_20250527.mat');
caseList(1).gradientTemplateFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S136_GradientXRange030_20250527.mat');
caseList(1).baseDynamicFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S136_SlidingWindows_20250527.mat');
caseList(1).gradientDynamicFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S136_SlidingWindows_GradientXRange030_20250527.mat');
caseList(1).mainResultFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'identification', 'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_GradientXRange030_20250527.mat');

caseList(2).tag = '20251222_B1_S123';
caseList(2).label = '20251222 original validated config';
caseList(2).sensorIds = [1 2 3];
caseList(2).mainWindow = 1;
caseList(2).baseTemplateFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S123_20251222.mat');
caseList(2).gradientTemplateFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
caseList(2).baseDynamicFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S123_SlidingWindows_20251222.mat');
caseList(2).gradientDynamicFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S123_SlidingWindows_GradientXRange030_20251222.mat');
caseList(2).mainResultFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'identification', 'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_GradientXRange030_20251222.mat');
end

function Template = load_template_local(file)
if exist(file, 'file') ~= 2
    error('Missing template file: %s', file);
end
S = load(file, 'Template');
Template = S.Template;
end

function DynamicMap = load_dynamic_local(file)
if exist(file, 'file') ~= 2
    error('Missing DynamicMap file: %s', file);
end
S = load(file, 'DynamicMap');
DynamicMap = S.DynamicMap;
end

function Result = load_result_local(file)
if exist(file, 'file') ~= 2
    error('Missing result file: %s', file);
end
S = load(file, 'Result');
Result = S.Result;
end

function fig = plot_static_template_crop_local(C, Tbase, Tgrad)
fig = figure('Name', [C.tag, ' static crop'], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 7.5]);
tiledlayout(fig, 1, numel(C.sensorIds), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(C.sensorIds)
    sid = C.sensorIds(is);
    B = get_sensor_local(Tbase, sid);
    G = get_sensor_local(Tgrad, sid);
    xBaseInGradientFrame = B.x_grid + B.xc - G.xc;
    nexttile;
    plot(xBaseInGradientFrame, B.v_grid, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0); hold on;
    plot(G.x_grid, G.v_grid, 'r-', 'LineWidth', 1.2);
    xline(0, 'k:', 'LineWidth', 0.8);
    xline(B.xc - G.xc, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    xlim(make_crop_xlim_local(G.x_grid, 1.1));
    grid on; box on;
    title(sprintf('CH%d, \\Delta xc=%.4f mm', sid, G.xc - B.xc), 'FontWeight', 'normal');
    xlabel('x relative to gradient template center (mm)');
    if is == 1
        ylabel('Template voltage (V)');
    end
    legend({'base template', 'GradientXRange030 template', 'gradient xc', 'base xc in gradient frame'}, ...
        'Location', 'best', 'FontSize', 7);
end
end

function fig = plot_dynamic_crop_overlay_local(C, Tgrad, Dbase, Dgrad, Rmain)
iw = min(C.mainWindow, numel(Dbase.Window));
Wbase = Dbase.Window(iw);
Wgrad = Dgrad.Window(iw);
fig = figure('Name', [C.tag, ' dynamic crop'], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 12]);
tiledlayout(fig, 2, numel(C.sensorIds), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(C.sensorIds)
    sid = C.sensorIds(is);
    Tpl = get_sensor_local(Tgrad, sid);
    Db = get_dynamic_sensor_local(Wbase, sid);
    Dg = get_dynamic_sensor_local(Wgrad, sid);
    xlimCrop = make_crop_xlim_local(Tpl.x_grid, 1.0);

    nexttile;
    plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.0); hold on;
    scatter_downsample_local(Db.x_rel, Db.V, 900, [0.1 0.35 0.9], 10);
    scatter_downsample_local(Dg.x_rel, Dg.V, 900, [0.9 0.15 0.1], 10);
    xline(0, 'k:', 'LineWidth', 0.8);
    xlim(xlimCrop);
    grid on; box on;
    title(sprintf('CH%d x-V crop, window %d', sid, iw), 'FontWeight', 'normal');
    xlabel('x_{rel} (mm)');
    if is == 1
        ylabel('Voltage (V)');
    end
    legend({'Gradient template', 'base DynamicMap', 'Gradient DynamicMap'}, ...
        'Location', 'best', 'FontSize', 7);

    nexttile;
    [tb, vb] = downsample_time_local(Db.t, Db.V, 1800);
    [tg, vg] = downsample_time_local(Dg.t, Dg.V, 1800);
    plot(tb - min(tb), vb, '-', 'Color', [0.1 0.35 0.9], 'LineWidth', 0.8); hold on;
    plot(tg - min(tg), vg, '-', 'Color', [0.9 0.15 0.1], 'LineWidth', 0.8);
    grid on; box on;
    title(sprintf('CH%d vibration waveform crop', sid), 'FontWeight', 'normal');
    xlabel('Time in selected window (s)');
    if is == 1
        ylabel('Voltage (V)');
    end
end

sgtitle(sprintf('%s | main: %s / %s | EO14 %d/%d', C.label, ...
    Rmain.MethodSettings.pulse_selection_mode, Rmain.MethodSettings.domain_selection_mode, ...
    nnz(Rmain.Trend.EO_id == 14), height(Rmain.Trend)), 'Interpreter', 'none');
end

function S = get_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template missing CH%d.', sid);
end
S = Template.Sensor(idx);
end

function D = get_dynamic_sensor_local(Wmap, sid)
idx = find([Wmap.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('DynamicMap window missing CH%d.', sid);
end
D = Wmap.Sensor(idx);
end

function xlimCrop = make_crop_xlim_local(x, scale)
lo = prctile(x(:), 8);
hi = prctile(x(:), 92);
xc = 0.5 * (lo + hi);
hw = 0.5 * (hi - lo) * scale;
xlimCrop = [xc - hw, xc + hw];
end

function scatter_downsample_local(x, y, nmax, color, markerSize)
idx = find(isfinite(x) & isfinite(y));
if numel(idx) > nmax
    idx = idx(round(linspace(1, numel(idx), nmax)));
end
scatter(x(idx), y(idx), markerSize, 'MarkerFaceColor', color, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25);
end

function [td, yd] = downsample_time_local(t, y, nmax)
idx = find(isfinite(t) & isfinite(y));
if numel(idx) > nmax
    idx = idx(round(linspace(1, numel(idx), nmax)));
end
td = t(idx);
yd = y(idx);
end

function save_figure_local(fig, outDir, stem)
pngFile = fullfile(outDir, [stem, '.png']);
pdfFile = fullfile(outDir, [stem, '.pdf']);
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
fprintf('Saved %s\n', pngFile);
end
