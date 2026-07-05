%% Step05V_VisualizeLowSpeedTemplateLibrary_20241106
% Visualize the low-speed template audit saved by Step05.
% This script is intentionally separate because Step05 reads all low-speed
% raw data and may be slower when repeated only for figure tuning.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();
if exist(P.files.lowSpeedTemplateAudit, 'file') ~= 2
    error(['Missing Step05 audit artifact:\n  %s\n\n' ...
        'Run Step05_BuildLowSpeedTemplateLibrary_20241106 first.'], ...
        P.files.lowSpeedTemplateAudit);
end
if exist(P.files.lowSpeedTemplateLibrary, 'file') ~= 2
    error('Missing Step05 library artifact:\n  %s', P.files.lowSpeedTemplateLibrary);
end

loadedAudit = load(P.files.lowSpeedTemplateAudit, 'LowSpeedTemplateAudit');
LowSpeedTemplateAudit = loadedAudit.LowSpeedTemplateAudit;
loadedLibrary = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loadedLibrary.LowSpeedTemplateLibrary;

figDir = fullfile(P.view.figureDir, '05A_low_speed_template_library');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

for ib = 1:numel(LowSpeedTemplateAudit.blade)
    bladeAudit = LowSpeedTemplateAudit.blade(ib);
    bladeId = bladeAudit.blade_id;
    if ~isfinite(bladeId)
        continue;
    end

    fig = figure('Name', sprintf('Step05V B%d low-speed template audit', bladeId), ...
        'Color', 'w', 'Position', [80, 80, 1500, 900], 'NumberTitle', 'off');
    tiledlayout(fig, numel(bladeAudit.sensor), 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for is = 1:numel(bladeAudit.sensor)
        A = bladeAudit.sensor(is);
        Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, A.sensor_id);
        [xWidePlot, vWidePlot] = downsample_cloud_local(A.x_wide, A.v_wide, 25000);
        [xSelPlot, vSelPlot] = downsample_cloud_local(A.x_selected, A.v_selected, 25000);

        nexttile;
        plot(xWidePlot, vWidePlot, '.', 'Color', [0.80 0.80 0.80], 'MarkerSize', 3, ...
            'DisplayName', 'wide cloud'); hold on;
        plot(xSelPlot, vSelPlot, 'k.', 'MarkerSize', 3, 'DisplayName', 'trust cloud');
        plot(Tpl.x_grid, Tpl.v_grid, 'r-', 'LineWidth', 1.5, 'DisplayName', 'template');
        xline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'final xc');
        if isfinite(A.xc_detected_rel)
            xline(A.xc_detected_rel, 'm:', 'LineWidth', 1.0, 'DisplayName', 'detected xc');
        end
        xline(A.x_domain_rel(1), 'b-.', 'LineWidth', 1.0, 'DisplayName', 'domain');
        xline(A.x_domain_rel(2), 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(A.threshold, 'Color', [0.40 0.40 0.40], 'LineStyle', '--', 'DisplayName', 'threshold');
        xlabel('x relative to final center (mm)');
        ylabel('Voltage (V)');
        title(sprintf('B%d CH%d cloud and template', bladeId, A.sensor_id));
        grid on; box on;
        legend('Location', 'best');

        nexttile;
        yyaxis left;
        plot(A.x_profile, A.v_profile, 'k-', 'LineWidth', 1.2);
        ylabel('Smoothed profile (V)');
        yyaxis right;
        plot(A.x_profile, A.gradient_abs, 'b-', 'LineWidth', 1.2, 'DisplayName', '|dV/dx|');
        hold on;
        yline(A.gradient_threshold, 'r--', 'LineWidth', 1.0, 'DisplayName', 'gradient threshold');
        ylabel('|dV/dx| (V/mm)');
        xlabel('x relative to final center (mm)');
        title(sprintf('B%d CH%d gradient profile', bladeId, A.sensor_id));
        grid on; box on;

        nexttile;
        plot(A.x_profile, double(A.effective_mask), 'g-', 'LineWidth', 1.2, 'DisplayName', 'effective mask'); hold on;
        scatter(A.x_selected, repmat(0.5, size(A.x_selected)), 4, A.weight_selected, 'filled', ...
            'DisplayName', 'selected points / weight');
        xline(A.x_domain_rel(1), 'b-.', 'LineWidth', 1.0, 'DisplayName', 'domain');
        xline(A.x_domain_rel(2), 'b-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        ylim([-0.05, 1.25]);
        xlabel('x relative to final center (mm)');
        ylabel('Mask / weight');
        title(sprintf(['B%d CH%d trust decision\n' ...
            '%s, rev=%d, pulse=%d'], ...
            bladeId, A.sensor_id, char(A.stable_plan_mode), ...
            numel(A.selected_revolution_ids), numel(A.selected_pulse_indices)));
        grid on; box on; colorbar;
    end

    if P.view.saveFigures
        exportgraphics(fig, fullfile(figDir, sprintf( ...
            'Step05V_LowSpeedTemplateAudit_B%d_20241106.png', bladeId)), 'Resolution', 300);
    end
end

fprintf('\n=== Step05V: visualize low-speed template audit ===\n');
fprintf('Audit source: %s\n', P.files.lowSpeedTemplateAudit);
fprintf('Figure folder: %s\n', figDir);

function Tpl = load_template_sensor_local(LowSpeedTemplateLibrary, bladeId, sid)
entry = LowSpeedTemplateLibrary.entry(bladeId);
if ~entry.exists
    error('Missing template bundle for B%d.', bladeId);
end
loadedBundle = load(char(entry.template_file), 'TemplateBundle');
row = loadedBundle.TemplateBundle.SensorFiles(loadedBundle.TemplateBundle.SensorFiles.sensor_id == sid, :);
if height(row) ~= 1
    error('Template bundle for B%d does not contain CH%d.', bladeId, sid);
end
loadedSensor = load(char(row.template_file(1)), 'sensor_template');
Tpl = loadedSensor.sensor_template.Sensor;
end

function [xDs, vDs] = downsample_cloud_local(x, v, maxPoints)
keep = isfinite(x) & isfinite(v);
x = x(keep);
v = v(keep);
if numel(x) > maxPoints
    idx = round(linspace(1, numel(x), maxPoints));
else
    idx = 1:numel(x);
end
xDs = x(idx);
vDs = v(idx);
end

