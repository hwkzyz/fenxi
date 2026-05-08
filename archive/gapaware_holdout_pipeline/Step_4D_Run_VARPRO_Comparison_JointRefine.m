%% Step 4D: compare baseline raw-scan with final local joint refinement
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_template_library.mat'), 'templateLib', 'templateLibFull');
load(fullfile(outDir, 'stage2_simulated_data.mat'), 'truth');
load(fullfile(outDir, 'stage3_gap_estimation.mat'), 'highMap', 'staticState');
load(fullfile(outDir, 'stage4_identification_results.mat'), 'fixedResult', 'gapAwareResult', 'oracleResult');

jointCfg = struct();
jointCfg.gHalfWidth = 0.05;
jointCfg.dxHalfWidth = 0.05;
jointCfg.fHalfWidth = 20;
jointCfg.ampScaleLo = 0.50;
jointCfg.ampScaleHi = 1.50;

jointRefinedResult = run_joint_refinement(highMap, templateLib, staticState, gapAwareResult, cfg, jointCfg);

ResultD = struct();
ResultD.truth = truth;
ResultD.staticState = staticState;
ResultD.fixed = fixedResult;
ResultD.rawScan = gapAwareResult;
ResultD.jointRefined = jointRefinedResult;
ResultD.oracle = oracleResult;
ResultD.jointCfg = jointCfg;
ResultD.config = cfg;

summaryTable = make_summary_table(ResultD);
save(fullfile(outDir, 'stage4d_jointrefine_results.mat'), ...
    'ResultD', 'jointRefinedResult', 'summaryTable', '-v7.3');

figure('Name', 'Step 4D - Joint Refinement Comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 11]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

labels = categorical(summaryTable.method);

nexttile;
bar(labels, summaryTable.gap_error_mm);
ylabel('Gap error (mm)');
title('Static-gap estimation error');

nexttile;
bar(labels, summaryTable.rmse);
ylabel('RMSE');
title('Waveform fitting error');

nexttile;
bar(labels, summaryTable.mean_freq_error_Hz);
ylabel('Mean frequency error (Hz)');
title('Frequency identification error');

nexttile; hold on;
idxPlot = 1:min(1200, numel(highMap.t_v));
plot(highMap.t_v(idxPlot) * 1000, highMap.V_a(idxPlot), 'k.', 'MarkerSize', 4, 'DisplayName', 'observed');
plot(highMap.t_v(idxPlot) * 1000, gapAwareResult.VFit(idxPlot), 'b-', ...
    'DisplayName', sprintf('raw-scan %.3f', gapAwareResult.g_used));
plot(highMap.t_v(idxPlot) * 1000, jointRefinedResult.VFit(idxPlot), 'm--', ...
    'DisplayName', sprintf('joint-refined %.3f', jointRefinedResult.g_used));
plot(highMap.t_v(idxPlot) * 1000, oracleResult.VFit(idxPlot), 'g-.', ...
    'DisplayName', sprintf('oracle %.3f', oracleResult.g_used));
xlabel('Time (ms)');
ylabel('Voltage');
title('Waveform fit comparison');
legend('Location', 'best', 'Box', 'off');

fprintf('[Step 4D] Local joint-refinement comparison:\n');
disp(summaryTable);

function result = run_joint_refinement(highMap, templateLib, staticState, gapAwareResult, cfg, jointCfg)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
p0 = gapAwareResult.p(:).';
theta0 = [staticState.gHat, staticState.dx0, p0];

gLo = max(0.05, staticState.gHat - jointCfg.gHalfWidth);
gHi = staticState.gHat + jointCfg.gHalfWidth;
dxLo = staticState.dx0 - jointCfg.dxHalfWidth;
dxHi = staticState.dx0 + jointCfg.dxHalfWidth;

a1 = max(0.001, jointCfg.ampScaleLo * abs(p0(1)));
a2 = max(0.001, jointCfg.ampScaleLo * abs(p0(4)));
b1 = min(0.8, max(a1 + 1e-6, jointCfg.ampScaleHi * abs(p0(1))));
b2 = min(0.8, max(a2 + 1e-6, jointCfg.ampScaleHi * abs(p0(4))));
f1Lo = max(min(cfg.f1Grid), p0(3) - jointCfg.fHalfWidth);
f1Hi = min(max(cfg.f1Grid), p0(3) + jointCfg.fHalfWidth);
f2Lo = max(min(cfg.f2Grid), p0(6) - jointCfg.fHalfWidth);
f2Hi = min(max(cfg.f2Grid), p0(6) + jointCfg.fHalfWidth);

lb = [gLo, dxLo, a1, -pi, f1Lo, a2, -pi, f2Lo];
ub = [gHi, dxHi, b1,  pi, f1Hi, b2,  pi, f2Hi];

resFun = @(theta) residual_calc_joint(theta, t, V, x, templateLib);
thetaOpt = refine_bounded_least_squares(resFun, theta0, lb, ub);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));

r = resFun(thetaOpt);
uHat = displacement_2freq(thetaOpt(3:end), t);
VFit = V - r;
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
phi_raw = [thetaOpt(4), thetaOpt(7)];

result = struct('label', 'joint-refined', ...
    'g_used', thetaOpt(1), ...
    'dx_used', thetaOpt(2), ...
    'p', thetaOpt(3:end), ...
    'theta', thetaOpt, ...
    'A_id', A_raw(order), ...
    'f_id', f_id, ...
    'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)), ...
    't', t, 'x', x, 'V', V, ...
    'VFit', VFit, 'uHat', uHat);
end

function r = residual_calc_joint(theta, t, V, x, templateLib)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = displacement_2freq(p, t);
VSim = evaluate_template_surface(templateLib, gFit, x - dx0 - u);
r = V - VSim;
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function pOpt = refine_bounded_least_squares(resFun, p0, lb, ub)
p0 = min(max(p0, lb), ub);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, 'MaxIterations', 500);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1500, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
end

function summaryTable = make_summary_table(ResultD)
labels = ["fixed_low_gap"; "raw_scan_gap"; "joint_refined_gap"; "oracle_holdout"];
methods = {ResultD.fixed, ResultD.rawScan, ResultD.jointRefined, ResultD.oracle};
[fTrue, trueOrder] = sort(ResultD.truth.f);
ATrue = ResultD.truth.A(trueOrder);
n = numel(methods);
gUsed = zeros(n, 1);
gapError = zeros(n, 1);
fErrMean = zeros(n, 1);
AErrMean = zeros(n, 1);
rmse = zeros(n, 1);
for i = 1:n
    m = methods{i};
    gUsed(i) = m.g_used;
    gapError(i) = abs(m.g_used - ResultD.truth.g_high);
    [fSort, methodOrder] = sort(m.f_id);
    ASort = m.A_id(methodOrder);
    fErrMean(i) = mean(abs(fSort - fTrue));
    AErrMean(i) = mean(abs(ASort - ATrue));
    rmse(i) = m.rmse;
end
summaryTable = table(labels, gUsed, gapError, fErrMean, AErrMean, rmse, ...
    'VariableNames', {'method', 'gap_used_mm', 'gap_error_mm', ...
    'mean_freq_error_Hz', 'mean_amp_error_mm', 'rmse'});
end

function yPred = evaluate_template_surface(templateLib, g, x)
[curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
yPred = interp1(templateLib.xGrid, curveGrid(:), x, 'pchip', 'extrap');
end

function [curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g)
zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
SSort = templateLib.S(order, :);
FxSort = templateLib.FxMat(order, :);
zq = 1 / max(g, 1e-6);
curveGrid = interp1(zSort, SSort, zq, 'linear', 'extrap');
dGrid = interp1(zSort, FxSort, zq, 'linear', 'extrap');
end

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end
