%% Step06E: Plot low-speed template and gap-library curves together
% This script is intentionally simple: it only overlays the measured
% low-speed template and the static gap-library curve in the same OPR-x
% coordinate figure. It helps visually check whether the library and
% low-speed waveform have the same x-axis position and shape.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06e_low_template_gap_overlay');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.showAllTrainGapCurves = true;

sensorOverride = strtrim(getenv('STEP06E_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06E_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
completeLibFile = fullfile(outDir, sprintf('Step06D_OPR_Anchored_GapLibrary_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
templateDir = fullfile(rootDir, '20250527_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(templateFile)
    error('Missing low-speed template file: %s', templateFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
St = load(templateFile, 'Template');
Template = St.Template;

if isfile(completeLibFile)
    C = load(completeLibFile, 'CompleteGapLibrary');
    CompleteGapLibrary = C.CompleteGapLibrary;
else
    warning('Missing Step06D complete library. Using tau=0 and middle training gap for overlay.');
    CompleteGapLibrary = make_fallback_complete_library_local(responseSurface, cfg);
end

fprintf('\n=== Step06E: low-speed template vs gap-library overlay ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Complete library: %s\n', completeLibFile);

figFile = fullfile(figDir, sprintf('Step06E_LowTemplate_GapLibrary_Overlay_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
plot_low_template_gap_overlay_local(Template, responseSurface, CompleteGapLibrary, cfg, figFile);

fprintf('Step06E complete. Figure saved and left open:\n  %s\n', figFile);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', targetBlade, exactTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20250527.mat', targetBlade)));
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)(?:_GradientXRange030)?_20250527\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    available = arrayfun(@(c) str2double(c), token{1});
    if all(ismember(sensorIds, available))
        templateFile = fullfile(files(i).folder, files(i).name);
        return;
    end
end
error('No low-speed template file found under %s for sensors %s.', templateDir, mat2str(sensorIds));
end

function CompleteGapLibrary = make_fallback_complete_library_local(responseSurface, cfg)
sensor = repmat(struct('sensorId', NaN, 'tauMm', 0, 'g0Mm', median(responseSurface.gTrainMm)), ...
    numel(cfg.analysisSensors), 1);
for i = 1:numel(cfg.analysisSensors)
    sensor(i).sensorId = cfg.analysisSensors(i);
end
CompleteGapLibrary = struct('sensor', sensor, 'xRelation', 'fallback: x_lib = x_OPR');
end

function plot_low_template_gap_overlay_local(Template, responseSurface, CompleteGapLibrary, cfg, figFile)
fig = figure('Name', 'Step06E low-template and gap-library overlay', 'Color', 'w', ...
    'Position', [50, 50, 1550, 390 * numel(cfg.analysisSensors)]);
tiledlayout(numel(cfg.analysisSensors), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    Slib = get_complete_library_sensor_local(CompleteGapLibrary, sid);
    x = Tsen.x_grid(:);
    vLow = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
    valid = isfinite(x) & isfinite(vLow);
    x = x(valid);
    vLow = vLow(valid);

    vGapRaw = eval_response_surface_local(responseSurface, Slib.g0Mm, x - Slib.tauMm);
    [vGapAffine, affine] = affine_fit_local(vLow, vGapRaw);

    nexttile; hold on; grid on; box on;
    if cfg.showAllTrainGapCurves
        for g = responseSurface.gTrainMm(:).'
            vg = eval_response_surface_local(responseSurface, g, x - Slib.tauMm);
            plot(x, vg, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6, ...
                'HandleVisibility', 'off');
        end
    end
    plot(x, vLow, 'k.', 'MarkerSize', 5, 'DisplayName', 'T low');
    plot(x, vGapRaw, 'r-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('F raw g=%.3f', Slib.g0Mm));
    title(sprintf('CH%d raw mV overlay | tau %.4f mm', sid, Slib.tauMm));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(x, vLow, 'k.', 'MarkerSize', 5, 'DisplayName', 'T low');
    plot(x, vGapAffine, 'r-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('affine F: a=%.3f b=%.1f', affine.gain, affine.offsetMv));
    title(sprintf('CH%d affine voltage overlay', sid));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(x, normalize_shape_local(vLow), 'k.', 'MarkerSize', 5, 'DisplayName', 'T low z');
    plot(x, normalize_shape_local(vGapRaw), 'r-', 'LineWidth', 1.4, 'DisplayName', 'F raw z');
    title(sprintf('CH%d normalized shape overlay', sid));
    xlabel('x_{OPR} (mm)'); ylabel('z-score'); legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function Slib = get_complete_library_sensor_local(CompleteGapLibrary, sid)
idx = find([CompleteGapLibrary.sensor.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('CompleteGapLibrary does not contain CH%d.', sid);
end
Slib = CompleteGapLibrary.sensor(idx);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function [pred, affine] = affine_fit_local(v, model)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model);
if nnz(valid) < 3
    affine = struct('gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
beta = H \ v(valid);
pred(valid) = H * beta;
affine = struct('gain', beta(1), 'offsetMv', beta(2));
end

function y = normalize_shape_local(x)
valid = isfinite(x);
y = nan(size(x));
if nnz(valid) < 3
    return;
end
sg = std(x(valid), 'omitnan');
if ~isfinite(sg) || sg < eps
    sg = 1;
end
y(valid) = (x(valid) - mean(x(valid), 'omitnan')) ./ sg;
end
