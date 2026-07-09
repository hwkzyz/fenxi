%% Cross-apply Step07J point selections to a no-gap direct-template model
% Purpose:
%   Compare the influence of point selection alone by fitting the same
%   low-speed direct-template/no-gap waveform model on:
%     1) Step07J direct-bundle selected points
%     2) Step07J no-direct self-selected points
%
% This script does not modify Step03/Step07J pipelines.

clear; clc;

rootDir = pwd;
rotDir = fullfile(rootDir, 'experimental_dynamic_decoupling_validation', ...
    '20250527_low_speed_rotating_calibration');
gapDir = fullfile(rootDir, 'experimental_dynamic_decoupling_validation', ...
    '20250527_low_speed_gap_prior_decoupling');
outDir = fullfile(rotDir, 'output', 'selection_effect_audit_20250527');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

directBundleFile = fullfile(gapDir, 'outputs', ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt.mat');
noDirectFile = fullfile(gapDir, 'outputs', ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt_nodirect.mat');
assert_file_local(directBundleFile);
assert_file_local(noDirectFile);

Sd = load(directBundleFile, 'Result');
Sn = load(noDirectFile, 'Result');
Rd = Sd.Result;
Rn = Sn.Result;

eoFixed = 14;
Td = fit_all_windows_local(Rd, "direct_bundle_points", eoFixed);
Tn = fit_all_windows_local(Rn, "no_direct_self_selected_points", eoFixed);
Tall = [Td; Tn];
summary = summarize_local(Tall);
pair = compare_pair_local(Tn, Td);

writetable(Tall, fullfile(outDir, 'NoGapDirectTemplate_OnStep07JSelections_Window_20250527_B1_S136.csv'));
writetable(summary, fullfile(outDir, 'NoGapDirectTemplate_OnStep07JSelections_Summary_20250527_B1_S136.csv'));
writetable(pair, fullfile(outDir, 'NoGapDirectTemplate_OnStep07JSelections_PairCompare_20250527_B1_S136.csv'));

fprintf('\n=== No-gap direct-template model on Step07J selections ===\n');
disp(summary);
fprintf('\n=== Pair comparison: no-direct selection -> direct-bundle selection ===\n');
disp(pair_summary_local(pair));
fprintf('\nSaved under:\n  %s\n', outDir);

function assert_file_local(f)
if exist(f, 'file') ~= 2
    error('Missing file: %s', f);
end
end

function T = fit_all_windows_local(Result, selectionName, eoFixed)
n = numel(Result.WindowResult);
rows = repmat(struct('selection', '', 'window_id', NaN, 'lap_start', NaN, ...
    'lap_end', NaN, 'point_count', NaN, 'mean_points_per_sensor', NaN, ...
    'EO', NaN, 'frequency_hz', NaN, 'amplitude_mm', NaN, ...
    'phase_rad', NaN, 'dx_mm', NaN, 'weighted_rmse_mV', NaN, ...
    'plain_rmse_mV', NaN, 'objective', NaN), n, 1);
for iw = 1:n
    wr = Result.WindowResult(iw);
    bundle = wr.bundle;
    fit = fit_no_gap_direct_template_local(bundle, eoFixed);
    rows(iw).selection = char(selectionName);
    rows(iw).window_id = iw;
    if isfield(wr, 'lapRange')
        rows(iw).lap_start = wr.lapRange(1);
        rows(iw).lap_end = wr.lapRange(end);
    else
        rows(iw).lap_start = iw;
        rows(iw).lap_end = iw;
    end
    rows(iw).point_count = numel(bundle.V);
    rows(iw).mean_points_per_sensor = numel(bundle.V) / max(numel(bundle.sensorIds), 1);
    rows(iw).EO = eoFixed;
    rows(iw).frequency_hz = eoFixed * bundle.rotFreqMeanHz;
    rows(iw).amplitude_mm = fit.A;
    rows(iw).phase_rad = fit.phi;
    rows(iw).dx_mm = fit.dx;
    rows(iw).weighted_rmse_mV = fit.weightedRmseMv;
    rows(iw).plain_rmse_mV = fit.plainRmseMv;
    rows(iw).objective = fit.objective;
end
T = struct2table(rows);
end

function fit = fit_no_gap_direct_template_local(bundle, eo)
seedA = 0.35;
seedPhi = 0;
seedDx = 0;
if isfield(bundle, 'VPSeedTable') %#ok<ISFLD>
end
theta0 = [seedA, seedPhi, seedDx];
opts = optimset('Display', 'off', 'MaxIter', 450, 'MaxFunEvals', 1800, ...
    'TolX', 1e-5, 'TolFun', 1e-5);
fun = @(th) objective_local(th, eo, bundle);
thetaOpt = fminsearch(fun, theta0, opts);
[objective, weightedRmse, plainRmse] = objective_local(thetaOpt, eo, bundle);
thetaOpt = bound_theta_local(thetaOpt);
fit = struct('A', thetaOpt(1), 'phi', thetaOpt(2), 'dx', thetaOpt(3), ...
    'objective', objective, 'weightedRmseMv', weightedRmse, 'plainRmseMv', plainRmse);
end

function [objective, weightedRmse, plainRmse] = objective_local(theta, eo, bundle)
theta = bound_theta_local(theta);
A = theta(1);
phi = theta(2);
dx = theta(3);
u = A .* sin(eo .* bundle.Theta(:) + phi);
vPred = nan(size(bundle.V(:)));
overshoot = zeros(size(bundle.V(:)));
for is = 1:numel(bundle.sensorIds)
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex(:) == is;
    Tpl = get_template_sensor_local(bundle.Template, sid);
    xEval = bundle.X(mask) - dx - u(mask);
    [vPred(mask), overshoot(mask)] = eval_low_template_local(Tpl, xEval);
end
valid = isfinite(vPred) & isfinite(bundle.V(:)) & isfinite(bundle.W(:));
pointCount = max(numel(bundle.V), 1);
if nnz(valid) < 8
    objective = inf;
    weightedRmse = inf;
    plainRmse = inf;
    return;
end
res = bundle.V(valid) - vPred(valid);
w = max(bundle.W(valid), 0.05);
residualObj = sum(w .* res.^2);
overshootPenalty = 100 * 1e6 * sum(w .* overshoot(valid).^2);
invalidPenalty = 1e12 * nnz(~valid);
objective = residualObj + overshootPenalty + invalidPenalty;
weightedRmse = sqrt(objective / pointCount);
plainRmse = sqrt(mean(res.^2, 'omitnan'));
end

function theta = bound_theta_local(theta)
theta = theta(:).';
theta(1) = max(min(theta(1), 0.50), 0);
theta(2) = mod(theta(2) + pi, 2*pi) - pi;
theta(3) = max(min(theta(3), 0.35), -0.35);
end

function [v, overshoot] = eval_low_template_local(Tpl, x)
xRaw = x(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xClamp = min(max(xRaw, xLo), xHi);
overshoot = max(xLo - xRaw, 0) + max(xRaw - xHi, 0);
v = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, ...
    xClamp, 'pchip', NaN);
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function summary = summarize_local(T)
selection = string(T.selection);
cases = unique(selection, 'stable');
rows = repmat(struct('selection', '', 'n_windows', NaN, 'mean_points', NaN, ...
    'median_points', NaN, 'dominant_EO', NaN, 'EO_consistency', NaN, ...
    'mean_frequency_hz', NaN, 'std_frequency_hz', NaN, ...
    'mean_amplitude_mm', NaN, 'std_amplitude_mm', NaN, ...
    'mean_weighted_rmse_mV', NaN, 'median_weighted_rmse_mV', NaN, ...
    'std_weighted_rmse_mV', NaN, 'mean_plain_rmse_mV', NaN), numel(cases), 1);
for i = 1:numel(cases)
    idx = selection == cases(i);
    Ti = T(idx, :);
    rows(i).selection = char(cases(i));
    rows(i).n_windows = height(Ti);
    rows(i).mean_points = mean(Ti.point_count, 'omitnan');
    rows(i).median_points = median(Ti.point_count, 'omitnan');
    rows(i).dominant_EO = mode(Ti.EO);
    rows(i).EO_consistency = mean(Ti.EO == rows(i).dominant_EO);
    rows(i).mean_frequency_hz = mean(Ti.frequency_hz, 'omitnan');
    rows(i).std_frequency_hz = std(Ti.frequency_hz, 'omitnan');
    rows(i).mean_amplitude_mm = mean(Ti.amplitude_mm, 'omitnan');
    rows(i).std_amplitude_mm = std(Ti.amplitude_mm, 'omitnan');
    rows(i).mean_weighted_rmse_mV = mean(Ti.weighted_rmse_mV, 'omitnan');
    rows(i).median_weighted_rmse_mV = median(Ti.weighted_rmse_mV, 'omitnan');
    rows(i).std_weighted_rmse_mV = std(Ti.weighted_rmse_mV, 'omitnan');
    rows(i).mean_plain_rmse_mV = mean(Ti.plain_rmse_mV, 'omitnan');
end
summary = struct2table(rows);
end

function C = compare_pair_local(A, B)
C = table();
C.window_id = A.window_id;
C.points_A_no_direct = A.point_count;
C.points_B_direct_bundle = B.point_count;
C.point_ratio_B_over_A = B.point_count ./ A.point_count;
C.rmse_A_no_direct_mV = A.weighted_rmse_mV;
C.rmse_B_direct_bundle_mV = B.weighted_rmse_mV;
C.rmse_delta_B_minus_A_mV = B.weighted_rmse_mV - A.weighted_rmse_mV;
C.rmse_ratio_B_over_A = B.weighted_rmse_mV ./ A.weighted_rmse_mV;
C.amp_A_no_direct_mm = A.amplitude_mm;
C.amp_B_direct_bundle_mm = B.amplitude_mm;
C.amp_delta_B_minus_A_mm = B.amplitude_mm - A.amplitude_mm;
end

function S = pair_summary_local(C)
S = table();
S.mean_points_no_direct = mean(C.points_A_no_direct, 'omitnan');
S.mean_points_direct_bundle = mean(C.points_B_direct_bundle, 'omitnan');
S.mean_point_ratio = mean(C.point_ratio_B_over_A, 'omitnan');
S.mean_rmse_no_direct_mV = mean(C.rmse_A_no_direct_mV, 'omitnan');
S.mean_rmse_direct_bundle_mV = mean(C.rmse_B_direct_bundle_mV, 'omitnan');
S.mean_rmse_delta_mV = mean(C.rmse_delta_B_minus_A_mV, 'omitnan');
S.mean_rmse_ratio = mean(C.rmse_ratio_B_over_A, 'omitnan');
S.mean_amp_delta_mm = mean(C.amp_delta_B_minus_A_mm, 'omitnan');
end
