%% Audit_TemplateCenter_OPRAnchored_20250527_20251222
% Audit whether template centers and dynamic-map centers follow the agreed
% OPR-anchored coordinate chain:
%   base low-speed template defines xc;
%   GradientXRange030 should inherit that xc and only change domain/weights;
%   BaseFrame is retained only as a diagnostic route.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'diagnostics_template_center_audit');
figDir = fullfile(outDir, 'figures');
if exist(outDir, 'dir') ~= 7; mkdir(outDir); end
if exist(figDir, 'dir') ~= 7; mkdir(figDir); end

cases = {
    '20250527', 'S136', [1 3 6], ...
        fullfile(rootDir, '20250527_low_speed_rotating_calibration', 'output');
    '20251222', 'S123', [1 2 3], ...
        fullfile(rootDir, '20251222_low_speed_rotating_calibration', 'output')
    };

allTemplateRows = {};
allDynamicRows = {};
for ic = 1:size(cases, 1)
    dataset = cases{ic, 1};
    sensorTag = cases{ic, 2};
    sensorIds = cases{ic, 3};
    outputRoot = cases{ic, 4};

    templateRows = audit_templates_local(outputRoot, dataset, sensorTag, sensorIds);
    dynamicRows = audit_dynamic_maps_local(outputRoot, dataset, sensorTag, sensorIds);
    allTemplateRows{ic} = templateRows; %#ok<AGROW>
    allDynamicRows{ic} = dynamicRows; %#ok<AGROW>

    plot_template_waveforms_local(outputRoot, dataset, sensorTag, sensorIds, figDir);
end

TemplateCenterAudit = vertcat(allTemplateRows{:});
DynamicCenterAudit = vertcat(allDynamicRows{:});
writetable(TemplateCenterAudit, fullfile(outDir, 'TemplateCenterAudit_OPRAnchored_20250527_20251222.csv'));
writetable(DynamicCenterAudit, fullfile(outDir, 'DynamicCenterAudit_OPRAnchored_20250527_20251222.csv'));
plot_center_bars_local(TemplateCenterAudit, DynamicCenterAudit, figDir);

fprintf('\nSaved template-center audit under:\n  %s\n', outDir);
disp(TemplateCenterAudit);
disp(DynamicCenterAudit);

%% Local functions
function rows = audit_templates_local(outputRoot, dataset, sensorTag, sensorIds)
templateDir = fullfile(outputRoot, 'templates');
templateSpecs = {
    'base', sprintf('Template_LowSpeedRotating_B1_%s_%s.mat', sensorTag, dataset)
    'opranchored', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_OPRAnchored_%s.mat', sensorTag, dataset)
    'gradient', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_%s.mat', sensorTag, dataset)
    'baseframe', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_BaseFrame_%s.mat', sensorTag, dataset)
    };

rowCells = {};
for it = 1:size(templateSpecs, 1)
    kind = templateSpecs{it, 1};
    file = fullfile(templateDir, templateSpecs{it, 2});
    if ~isfile(file)
        continue;
    end
    loaded = load(file, 'Template');
    Template = loaded.Template;
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        idx = find([Template.Sensor.sensor_id] == sid, 1);
        if isempty(idx)
            continue;
        end
        S = Template.Sensor(idx);
        xcDetected = NaN;
        if isfield(S, 'xc_detected') && ~isempty(S.xc_detected)
            xcDetected = S.xc_detected;
        end
        source = "";
        if isfield(S, 'xc_reference_source') && ~isempty(S.xc_reference_source)
            source = string(S.xc_reference_source);
        end
        dxDetectedMinusUsed = xcDetected - S.xc;
        rowCells{end+1, 1} = table(string(dataset), string(sensorTag), sid, string(kind), ... %#ok<AGROW>
            string(file), S.xc, xcDetected, dxDetectedMinusUsed, source, ...
            S.x_domain(1), S.x_domain(2), diff(S.x_domain), numel(S.x_grid), ...
            'VariableNames', {'dataset','sensor_tag','sensor_id','template_kind', ...
            'template_file','xc_used_mm','xc_detected_mm','dx_detected_minus_used_mm', ...
            'xc_reference_source','x_domain_left_mm','x_domain_right_mm', ...
            'x_domain_width_mm','x_grid_count'});
    end
end
if isempty(rowCells)
    rows = table();
else
    rows = vertcat(rowCells{:});
end
end

function rows = audit_dynamic_maps_local(outputRoot, dataset, sensorTag, sensorIds)
dynamicDir = fullfile(outputRoot, 'dynamic_maps');
mapSpecs = {
    'base', sprintf('DynamicMap_B1_%s_SlidingWindows_%s.mat', sensorTag, dataset)
    'gradient', sprintf('DynamicMap_B1_%s_SlidingWindows_GradientXRange030_%s.mat', sensorTag, dataset)
    };

rowCells = {};
for im = 1:size(mapSpecs, 1)
    kind = mapSpecs{im, 1};
    file = fullfile(dynamicDir, mapSpecs{im, 2});
    if ~isfile(file)
        continue;
    end
    loaded = load(file, 'DynamicMap');
    DynamicMap = loaded.DynamicMap;
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        xCenterMm = NaN;
        if isfield(DynamicMap, 'XCenterBySensor')
            X = DynamicMap.XCenterBySensor;
            names = X.Properties.VariableNames;
            if all(ismember({'SensorID', 'XCenterMm'}, names))
                row = X.SensorID == sid;
                if any(row)
                    xCenterMm = X.XCenterMm(find(row, 1));
                end
            elseif all(ismember({'sensor_id', 'xc_mm'}, names))
                row = X.sensor_id == sid;
                if any(row)
                    xCenterMm = X.xc_mm(find(row, 1));
                end
            else
                error('Unsupported XCenterBySensor variable names in %s.', file);
            end
        end
        rowCells{end+1, 1} = table(string(dataset), string(sensorTag), sid, string(kind), ... %#ok<AGROW>
            string(file), xCenterMm, ...
            'VariableNames', {'dataset','sensor_tag','sensor_id','dynamic_map_kind', ...
            'dynamic_map_file','x_center_mm'});
    end
end
if isempty(rowCells)
    rows = table();
else
    rows = vertcat(rowCells{:});
end
end

function plot_template_waveforms_local(outputRoot, dataset, sensorTag, sensorIds, figDir)
templateDir = fullfile(outputRoot, 'templates');
templateSpecs = {
    'base', sprintf('Template_LowSpeedRotating_B1_%s_%s.mat', sensorTag, dataset), [0.15 0.15 0.15]
    'opranchored', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_OPRAnchored_%s.mat', sensorTag, dataset), [0.00 0.50 0.25]
    'gradient', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_%s.mat', sensorTag, dataset), [0.85 0.10 0.10]
    'baseframe', sprintf('Template_LowSpeedRotating_B1_%s_GradientXRange030_BaseFrame_%s.mat', sensorTag, dataset), [0.00 0.35 0.75]
    };

fig = figure('Name', ['Template center audit ' dataset], 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 7.5]);
tiledlayout(fig, 1, numel(sensorIds), 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    nexttile;
    hold on;
    for it = 1:size(templateSpecs, 1)
        file = fullfile(templateDir, templateSpecs{it, 2});
        if ~isfile(file)
            continue;
        end
        loaded = load(file, 'Template');
        Template = loaded.Template;
        idx = find([Template.Sensor.sensor_id] == sid, 1);
        if isempty(idx)
            continue;
        end
        S = Template.Sensor(idx);
        plot(S.x_grid(:), S.v_grid(:), '-', 'Color', templateSpecs{it, 3}, ...
            'LineWidth', 1.2, 'DisplayName', templateSpecs{it, 1});
        xline(S.x_domain(1), ':', 'Color', templateSpecs{it, 3}, 'HandleVisibility', 'off');
        xline(S.x_domain(2), ':', 'Color', templateSpecs{it, 3}, 'HandleVisibility', 'off');
    end
    xline(0, 'k--', 'LineWidth', 0.9, 'DisplayName', 'x=0');
    title(sprintf('%s CH%d', dataset, sid), 'Interpreter', 'none');
    xlabel('x relative to used xc (mm)');
    ylabel('Template voltage (V)');
    box on; grid on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8);
    if is == numel(sensorIds)
        legend('Location', 'best');
    end
end
export_figure_local(fig, fullfile(figDir, sprintf('TemplateWaveformCenterAudit_%s_%s.png', dataset, sensorTag)));
end

function plot_center_bars_local(TemplateCenterAudit, DynamicCenterAudit, figDir)
fig = figure('Name', 'Template and dynamic-map center audit', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 12]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
if ~isempty(TemplateCenterAudit)
    labels = TemplateCenterAudit.dataset + " CH" + string(TemplateCenterAudit.sensor_id) + " " + TemplateCenterAudit.template_kind;
    bar(categorical(labels), TemplateCenterAudit.xc_used_mm);
    ylabel('Template xc used (mm)');
    title('Template center used for coordinate frame');
    xtickangle(35);
    grid on; box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8);
end

nexttile;
if ~isempty(DynamicCenterAudit)
    labels = DynamicCenterAudit.dataset + " CH" + string(DynamicCenterAudit.sensor_id) + " " + DynamicCenterAudit.dynamic_map_kind;
    bar(categorical(labels), DynamicCenterAudit.x_center_mm);
    ylabel('DynamicMap x center (mm)');
    title('Dynamic-map saved center');
    xtickangle(35);
    grid on; box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8);
end

export_figure_local(fig, fullfile(figDir, 'TemplateDynamicCenterAudit_OPRAnchored_20250527_20251222.png'));
end

function export_figure_local(fig, pngFile)
[folder, name] = fileparts(pngFile);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
end
