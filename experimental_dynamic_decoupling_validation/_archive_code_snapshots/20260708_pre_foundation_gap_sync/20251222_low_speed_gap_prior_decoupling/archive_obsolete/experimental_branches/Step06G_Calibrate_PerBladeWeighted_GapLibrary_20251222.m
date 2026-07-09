%% Step06G: Calibrate per-blade weighted gap library
% Experimental branch. It does not overwrite Step06F outputs.
%
% For each sensor, every per-blade static response family F_b(g,x) is fitted
% to the low-speed rotating template. The low-speed fit quality defines a
% fixed blade-family weight used later by Step07D.

clear; clc; close all;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06g_per_blade_weighted_gap_library');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 2, 3];
cfg.gBoundsMm = [];
cfg.tauBoundsMm = [-0.80, 0.80];
cfg.xScaleBounds = [0.65, 1.45];
cfg.muBounds = [-0.35, 0.35];
cfg.minFitPoints = 100;
cfg.maxLowFitPoints = 500;
cfg.weightTemperatureMv = 3;
cfg.minBladeWeight = 0.02;
cfg.maxBladeFamilies = inf;
cfg.overshootPenaltyMvPerMm = 600;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 220, 'MaxFunEvals', 750, ...
    'TolX', 1e-6, 'TolFun', 1e-6);

sensorOverride = strtrim(getenv('STEP06G_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06G_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
tempOverride = str2double(strtrim(getenv('STEP06G_WEIGHT_TEMPERATURE_MV')));
if isfinite(tempOverride) && tempOverride > 0
    cfg.weightTemperatureMv = tempOverride;
end
maxFamiliesOverride = str2double(strtrim(getenv('STEP06G_MAX_BLADE_FAMILIES')));
if isfinite(maxFamiliesOverride) && maxFamiliesOverride > 0
    cfg.maxBladeFamilies = floor(maxFamiliesOverride);
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

perBladeFile = fullfile(outDir, 'Step05G_PerBlade_Response_Surface_20251222.mat');
templateDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(perBladeFile)
    error('Run Step05G first. Missing file: %s', perBladeFile);
end
if ~isfile(templateFile)
    error('Missing low-speed template file: %s', templateFile);
end

Slib = load(perBladeFile, 'PerBladeResponseLibrary');
St = load(templateFile, 'Template');
PerBladeResponseLibrary = Slib.PerBladeResponseLibrary;
Template = St.Template;
if isempty(cfg.gBoundsMm)
    cfg.gBoundsMm = [min(PerBladeResponseLibrary.trueGapMm), max(PerBladeResponseLibrary.trueGapMm)];
end

fprintf('\n=== Step06G: per-blade weighted gap-library calibration ===\n');
fprintf('Per-blade library: %s\n', perBladeFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Sensors: %s, blade families: %s\n', mat2str(cfg.analysisSensors), mat2str(PerBladeResponseLibrary.bladeIds));

sensor = repmat(empty_sensor_local(), numel(cfg.analysisSensors), 1);
fitRows = {};
weightRows = {};

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    candidate = repmat(empty_candidate_local(), numel(PerBladeResponseLibrary.blade), 1);
    for ib = 1:numel(PerBladeResponseLibrary.blade)
        B = PerBladeResponseLibrary.blade(ib);
        if isempty(B.coeff)
            continue;
        end
        fit = fit_tilt_corrected_low_template_local(Tsen, B, cfg);
        fit.bladeId = B.bladeId;
        fit.bladeIndex = ib;
        candidate(ib) = fit;
        fitRows{end+1, 1} = table(sid, B.bladeId, ib, fit.g0Mm, fit.tauMm, ...
            fit.xScale, fit.muGapPerXMm, fit.tiltAngleDeg, fit.voltageGain, ...
            fit.voltageOffsetMv, fit.lowFitRmseMv, fit.overshootRmseMm, fit.pointCount, ...
            'VariableNames', {'sensorId','bladeId','bladeIndex','g0Mm','tauMm', ...
            'xScale','muGapPerXMm','tiltAngleDeg','voltageGain','voltageOffsetMv', ...
            'lowFitRmseMv','overshootRmseMm','pointCount'}); %#ok<SAGROW>
    end
    candidate = candidate(isfinite([candidate.lowFitRmseMv]));
    if isempty(candidate)
        error('CH%d has no valid per-blade low-speed fit.', sid);
    end
    [candidate, weights] = assign_blade_weights_local(candidate, cfg);
    for ib = 1:numel(candidate)
        candidate(ib).bladeWeight = weights(ib);
        weightRows{end+1, 1} = table(sid, candidate(ib).bladeId, candidate(ib).bladeIndex, ...
            candidate(ib).lowFitRmseMv, candidate(ib).bladeWeight, ...
            'VariableNames', {'sensorId','bladeId','bladeIndex','lowFitRmseMv','bladeWeight'}); %#ok<SAGROW>
    end
    sensor(is).sensorId = sid;
    sensor(is).candidate = candidate;
    sensor(is).bestBladeId = candidate(1).bladeId;
    sensor(is).bestBladeIndex = candidate(1).bladeIndex;
    sensor(is).bestLowFitRmseMv = candidate(1).lowFitRmseMv;
    sensor(is).weightedLowFitRmseMv = sqrt(sum(weights(:) .* [candidate.lowFitRmseMv].'.^2) / max(sum(weights), eps));
    sensor(is).activeBladeIds = [candidate.bladeId];
    sensor(is).activeBladeWeights = [candidate.bladeWeight];

    fprintf('  CH%d: best B%d RMSE %.2f mV | active %s | weights %s\n', ...
        sid, sensor(is).bestBladeId, sensor(is).bestLowFitRmseMv, ...
        mat2str(sensor(is).activeBladeIds), mat2str(sensor(is).activeBladeWeights, 3));
end

fitTable = vertcat(fitRows{:});
weightTable = vertcat(weightRows{:});

PerBladeWeightedGapLibrary = struct();
PerBladeWeightedGapLibrary.dataset = '20251222';
PerBladeWeightedGapLibrary.method = 'low_template_anchored_per_blade_weighted_gap_increment';
PerBladeWeightedGapLibrary.description = ['Each static blade family is fitted to the low-speed template. ' ...
    'Blade-family weights are fixed before high-speed identification; no per-window blade switching is allowed.'];
PerBladeWeightedGapLibrary.perBladeFile = perBladeFile;
PerBladeWeightedGapLibrary.templateFile = templateFile;
PerBladeWeightedGapLibrary.targetBlade = cfg.targetBlade;
PerBladeWeightedGapLibrary.analysisSensors = cfg.analysisSensors(:).';
PerBladeWeightedGapLibrary.PerBladeResponseLibrary = PerBladeResponseLibrary;
PerBladeWeightedGapLibrary.sensor = sensor;
PerBladeWeightedGapLibrary.cfg = cfg;
PerBladeWeightedGapLibrary.formula = ['T_low(x)+sum_b w_b*a_b*(F_b(g0_b+dg+mu_b*(x-tau_b),k_b*(x-tau_b))' ...
    '-F_b(g0_b+mu_b*(x-tau_b),k_b*(x-tau_b)))'];

matFile = fullfile(outDir, sprintf('Step06G_PerBladeWeighted_GapLibrary_20251222_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
fitCsv = fullfile(outDir, sprintf('Step06G_PerBladeWeighted_GapLibrary_Fit_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
weightCsv = fullfile(outDir, sprintf('Step06G_PerBladeWeighted_GapLibrary_Weights_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figFit = fullfile(figDir, sprintf('Step06G_PerBlade_LowTemplate_Fits_20251222_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figWeights = fullfile(figDir, sprintf('Step06G_BladeFamily_Weights_20251222_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'PerBladeWeightedGapLibrary', 'fitTable', 'weightTable', '-v7.3');
writetable(fitTable, fitCsv);
writetable(weightTable, weightCsv);
plot_low_fit_candidates_local(sensor, figFit);
plot_weight_summary_local(sensor, figWeights);

fprintf('\nStep06G complete.\n');
disp(weightTable);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, weightCsv, figWeights);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
error('No low-speed template file found under %s for %s.', templateDir, sensorTag);
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function sensor = empty_sensor_local()
sensor = struct('sensorId', NaN, 'candidate', [], 'bestBladeId', NaN, ...
    'bestBladeIndex', NaN, 'bestLowFitRmseMv', NaN, 'weightedLowFitRmseMv', NaN, ...
    'activeBladeIds', [], 'activeBladeWeights', []);
end

function fit = empty_candidate_local()
fit = struct('sensorId', NaN, 'bladeId', NaN, 'bladeIndex', NaN, ...
    'g0Mm', NaN, 'tauMm', NaN, 'xScale', NaN, 'muGapPerXMm', NaN, ...
    'tiltAngleDeg', NaN, 'voltageGain', NaN, 'voltageOffsetMv', NaN, ...
    'lowFitRmseMv', NaN, 'overshootRmseMm', NaN, 'pointCount', NaN, ...
    'bladeWeight', NaN, 'x', [], 'vLowMv', [], 'vFitMv', [], ...
    'xLib', [], 'gEff', [], 'residualMv', []);
end

function fit = fit_tilt_corrected_low_template_local(Tsen, bladeSurface, cfg)
x = Tsen.x_grid(:);
vLow = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
mask = isfinite(x) & isfinite(vLow);
x = x(mask);
vLow = vLow(mask);
if numel(x) < cfg.minFitPoints
    error('CH%d has too few low-template points.', Tsen.sensor_id);
end
if numel(x) > cfg.maxLowFitPoints
    idx = unique(round(linspace(1, numel(x), cfg.maxLowFitPoints)));
    x = x(idx);
    vLow = vLow(idx);
end

gStarts = unique([median(bladeSurface.gTrainMm), 1.04]);
tauStarts = [-0.25, 0, 0.25];
kStarts = 1.0;
muStarts = [-0.08, 0, 0.08];
best = struct('objective', inf);
for g0 = gStarts(:).'
    for tau = tauStarts
        for k = kStarts
            for mu = muStarts
                p0 = [g0, tau, k, mu];
                fun = @(p) tilt_low_template_objective_scalar_local(p, x, vLow, bladeSurface, cfg);
                pOpt = fminsearch(fun, p0, cfg.fminOptions);
                [obj, detail] = evaluate_tilt_low_template_local(pOpt, x, vLow, bladeSurface, cfg);
                if obj < best.objective
                    best = detail;
                    best.objective = obj;
                end
            end
        end
    end
end
fit = empty_candidate_local();
fit.sensorId = Tsen.sensor_id;
fit.g0Mm = best.g0Mm;
fit.tauMm = best.tauMm;
fit.xScale = best.xScale;
fit.muGapPerXMm = best.muGapPerXMm;
fit.tiltAngleDeg = atan(best.muGapPerXMm) * 180 / pi;
fit.voltageGain = best.voltageGain;
fit.voltageOffsetMv = best.voltageOffsetMv;
fit.lowFitRmseMv = best.lowFitRmseMv;
fit.overshootRmseMm = best.overshootRmseMm;
fit.pointCount = best.pointCount;
fit.x = x;
fit.vLowMv = vLow;
fit.vFitMv = best.vFitMv;
fit.xLib = best.xLib;
fit.gEff = best.gEff;
fit.residualMv = vLow - best.vFitMv;
end

function obj = tilt_low_template_objective_scalar_local(p, x, vLow, bladeSurface, cfg)
[obj, ~] = evaluate_tilt_low_template_local(p, x, vLow, bladeSurface, cfg);
end

function [obj, detail] = evaluate_tilt_low_template_local(pRaw, x, vLow, bladeSurface, cfg)
p = pRaw(:).';
[g0, tau, k, mu, penalty] = clamp_tilt_params_local(p, cfg);
[modelRaw, xLib, gEff, overshoot] = eval_raw_tilt_path_local(bladeSurface, g0, tau, k, mu, x);
[rmse, vFit, affine] = affine_rmse_local(vLow, modelRaw, ones(size(vLow)));
overshootRmse = sqrt(mean(overshoot.^2, 'omitnan'));
obj = rmse + cfg.overshootPenaltyMvPerMm * overshootRmse + penalty;
detail = struct('g0Mm', g0, 'tauMm', tau, 'xScale', k, 'muGapPerXMm', mu, ...
    'voltageGain', affine.gain, 'voltageOffsetMv', affine.offsetMv, ...
    'lowFitRmseMv', rmse, 'overshootRmseMm', overshootRmse, ...
    'pointCount', nnz(isfinite(vFit)), 'vFitMv', vFit, ...
    'xLib', xLib, 'gEff', gEff, 'overshoot', overshoot);
end

function [g0, tau, k, mu, penalty] = clamp_tilt_params_local(p, cfg)
gRaw = p(1); tauRaw = p(2); kRaw = p(3); muRaw = p(4);
g0 = min(max(gRaw, cfg.gBoundsMm(1)), cfg.gBoundsMm(2));
tau = min(max(tauRaw, cfg.tauBoundsMm(1)), cfg.tauBoundsMm(2));
k = min(max(kRaw, cfg.xScaleBounds(1)), cfg.xScaleBounds(2));
mu = min(max(muRaw, cfg.muBounds(1)), cfg.muBounds(2));
penalty = 1e5 * ((gRaw - g0)^2 + (tauRaw - tau)^2 + ...
    (kRaw - k)^2 + (muRaw - mu)^2);
end

function [model, xLib, gEff, overshoot] = eval_raw_tilt_path_local(bladeSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr(:) - tau);
gEffRaw = g0 + mu .* (xOpr(:) - tau);
xMin = min(bladeSurface.xGrid(:));
xMax = max(bladeSurface.xGrid(:));
gMin = min(bladeSurface.gTrainMm(:));
gMax = max(bladeSurface.gTrainMm(:));
xLib = min(max(xLibRaw, xMin), xMax);
gEff = min(max(gEffRaw, gMin), gMax);
overshoot = hypot(max(xMin - xLibRaw, 0) + max(xLibRaw - xMax, 0), ...
    max(gMin - gEffRaw, 0) + max(gEffRaw - gMax, 0));
model = eval_response_surface_local(bladeSurface, gEff, xLib);
end

function F = eval_response_surface_local(bladeSurface, g, xq)
B0 = interp1(bladeSurface.xGrid, bladeSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(bladeSurface.xGrid, bladeSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(bladeSurface.xGrid, bladeSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ bladeSurface.g0Mm);
end

function [rmse, pred, affine] = affine_rmse_local(v, model, w)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model) & isfinite(w);
if nnz(valid) < 3
    rmse = inf;
    affine = struct('gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
beta = (H .* sw) \ (v(valid) .* sw);
pred(valid) = H * beta;
res = v(valid) - pred(valid);
rmse = sqrt(sum(w(valid) .* res.^2) / max(sum(w(valid)), eps));
affine = struct('gain', beta(1), 'offsetMv', beta(2));
end

function [candidate, weights] = assign_blade_weights_local(candidate, cfg)
[~, order] = sort([candidate.lowFitRmseMv], 'ascend');
candidate = candidate(order);
if isfinite(cfg.maxBladeFamilies)
    candidate = candidate(1:min(numel(candidate), cfg.maxBladeFamilies));
end
rmse = [candidate.lowFitRmseMv].';
rmse0 = min(rmse);
weights = exp(-0.5 * ((rmse - rmse0) ./ cfg.weightTemperatureMv).^2);
weights(weights < cfg.minBladeWeight) = 0;
if ~any(weights > 0)
    weights(1) = 1;
end
weights = weights ./ sum(weights);
end

function plot_low_fit_candidates_local(sensor, figFile)
style = paper_style_local();
fig = figure('Name', 'Step06G per-blade low-template fits', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.5, 4.9 * numel(sensor)]);
tiledlayout(numel(sensor), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(sensor)
    S = sensor(is);
    nexttile; hold on; box on;
    C = S.candidate;
    if isempty(C)
        continue;
    end
    plot(C(1).x, C(1).vLowMv, '.', 'Color', [0.25 0.25 0.25], ...
        'MarkerSize', 4.0, 'DisplayName', 'T low');
    nShow = min(3, numel(C));
    colors = lines(nShow);
    for ib = 1:nShow
        plot(C(ib).x, C(ib).vFitMv, '-', 'Color', colors(ib, :), ...
            'LineWidth', 1.1, 'DisplayName', sprintf('B%d, w=%.2f', C(ib).bladeId, C(ib).bladeWeight));
    end
    title(sprintf('CH%d per-blade low-speed fit', S.sensorId));
    xlabel('x_{OPR} (mm)', 'Interpreter', 'tex'); ylabel('mV');
    legend('Location', 'eastoutside', 'Box', 'off');
    format_axes_local(gca, style);
end
export_paper_figure_local(fig, figFile);
end

function plot_weight_summary_local(sensor, figFile)
style = paper_style_local();
fig = figure('Name', 'Step06G blade-family weights', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.5, 9.0]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on;
for is = 1:numel(sensor)
    C = sensor(is).candidate;
    plot([C.bladeId], [C.lowFitRmseMv], '-o', 'LineWidth', 1.0, ...
        'MarkerSize', 4, 'DisplayName', sprintf('CH%d', sensor(is).sensorId));
end
xlabel('Static-library blade ID');
ylabel('Low-template fit RMSE (mV)');
title('Blade-family matching error');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
for is = 1:numel(sensor)
    C = sensor(is).candidate;
    stem([C.bladeId] + 0.06 * (is - 1), [C.bladeWeight], 'filled', ...
        'LineWidth', 1.0, 'DisplayName', sprintf('CH%d', sensor(is).sensorId));
end
xlabel('Static-library blade ID');
ylabel('Weight');
title('Fixed blade-family weights');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.tickFontSize = 8.2;
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on');
grid(ax, 'off');
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
end

function export_paper_figure_local(fig, pngFile)
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
try
    exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
catch
end
try
    print(fig, fullfile(folder, [name, '.emf']), '-dmeta', '-r300');
catch
end
end

