%% 单工况主程序：手动运行入口
% 用途：
%   1) 只针对一个工况，直接运行当前默认主方法；
%   2) 手动修改参数后，快速查看识别结果；
%   3) 自动生成单工况波形拟合图，便于直观看结果。

close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

% 是否额外保存较大的 .mat 结果文件。
% 手动反复试参数时，建议保持 false。
saveMatResult = false;

% 是否保存当前图窗为 PNG。
saveFigurePng = true;

%% 1. 手动设置单工况参数
caseCfg = struct();

% 真实静态间隙，必须是静态库中已有的 gap 点之一。
caseCfg.g_true_mm = 0.8;

% 双频振动幅值，单位 mm。
caseCfg.A_true_mm = [0.25, 0.15];

% 双频振动频率，单位 Hz。
caseCfg.f_true_Hz = [500, 1300];

% 双频初相位，单位 rad。
caseCfg.phi_true_rad = [pi/4, -pi/3];

% 信噪比，单位 dB。若设为 Inf，则表示无噪声。
caseCfg.snr_dB = 10;

% 随机种子。
caseCfg.seed = 20270513;

%% 2. 载入项目上下文与静态库
cfgAna = ctx.cfgAna;
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
trustInfo = load_fixed_trust_domain(ctx.thisDir);

gTrue = caseCfg.g_true_mm;
idxTrue = find(abs(gapList - gTrue) < 1e-12, 1);
if isempty(idxTrue)
    error('g_true_mm = %.6g 不在当前静态库 gapList 中。', gTrue);
end

modelDef = get_inv_log_2_model_def();
templateTruth = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);
templateLib = make_fixed_trust_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelDef, trustInfo);
FTrue = griddedInterpolant(xCell{idxTrue}, yCell{idxTrue}, 'pchip', 'nearest');

%% 3. 生成该工况下的高速波形
cfgCase = cfgAna;
cfgCase.g_holdout = gTrue;
cfgCase.A_true = caseCfg.A_true_mm;
cfgCase.f_true = caseCfg.f_true_Hz;
cfgCase.phi_true = caseCfg.phi_true_rad;
cfgCase.noiseMode = 'noise_ratio';
cfgCase.snrDb = caseCfg.snr_dB;

if isinf(caseCfg.snr_dB)
    cfgCase.noiseRatio = 0;
else
    cfgCase.noiseRatio = calibrate_single_case_snr_ratio(cfgCase, gapList, xCell, yCell, templateTruth.domain);
end

rng(caseCfg.seed, 'twister');
dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
    cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

%% 4. 映射到空间坐标并运行当前默认主方法
highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
    templateLib.domain, cfgCase.fitActiveLevel);

% The formal main route uses an independent low-speed, no-vibration record
% to initialize the gap profile. Low and high records share the same SNR.
cfgLow = cfgCase;
cfgLow.RPM_low = min(cfgCase.RPM_high, 300);
cfgLow.NumRevs_low = 8;
dataLow = simulate_low_speed_template(FTrue, cfgLow, templateLib.domain, ...
    caseCfg.snr_dB, caseCfg.seed + 1);

timerMethod = tic;
result = run_inv_log_2_low_high_main(dataLow, highMap, templateLib, cfgCase);
result.elapsed_s = toc(timerMethod);
staticState = result.low_speed_state;

%% 5. 汇总结论并保存轻量结果
truthCase = struct('g_high', gTrue, 'A', cfgCase.A_true, 'f', cfgCase.f_true, 'phi', cfgCase.phi_true);
summaryRow = make_three_method_summary(result, truthCase);
summaryRow.g_true_mm = gTrue;
summaryRow.snr_dB = caseCfg.snr_dB;
summaryRow.static_g0_mm = staticState.gHat;
summaryRow.static_gap_error_mm = abs(staticState.gHat - gTrue);
summaryRow.static_dx0_mm = staticState.dx0;
summaryRow.dx_used_mm = result.dx_used;
summaryRow.dx_update_mm = result.dx_used - staticState.dx0;
summaryRow.trust_x_min_mm = templateLib.domain(1);
summaryRow.trust_x_max_mm = templateLib.domain(2);
summaryRow.trust_half_width_mm = 0.5 * diff(templateLib.domain);
summaryRow.seed = caseCfg.seed;
summaryRow.fast_mode_used = get_struct_field(result, 'final_fast_mode_used', false);

outDir = fullfile(ctx.outDir, 'main_single_case_manual');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

csvPath = fullfile(outDir, 'single_case_result.csv');
writetable(summaryRow, csvPath);

matPath = fullfile(outDir, 'single_case_result.mat');
if saveMatResult
    save(matPath, 'caseCfg', 'cfgCase', 'staticState', 'result', 'summaryRow', '-v7.3');
end

%% 6. 自动画出单工况波形拟合图
fig = build_single_case_figure(highMap, result, truthCase, summaryRow);

pngPath = fullfile(outDir, 'single_case_waveform_fit.png');
if saveFigurePng
    exportgraphics(fig, pngPath, 'Resolution', 200);
end

%% 7. 命令行输出
disp('单工况主程序运行完成。');
disp(summaryRow);
fprintf('结果表已保存到：\n%s\n', csvPath);
if saveFigurePng
    fprintf('波形图已保存到：\n%s\n', pngPath);
end
if saveMatResult
    fprintf('MAT 文件已保存到：\n%s\n', matPath);
end

% 手动反复试参数时，不保留过多中间量。
clear dataHigh templateTruth trustInfo xCell yCell gapList idxTrue modelDef;

function fig = build_single_case_figure(highMap, result, truthCase, summaryRow)
fig = figure('Name', '单工况波形拟合', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, highMap.x_v(:), highMap.V_a(:), 'k.', 'MarkerSize', 6, 'DisplayName', '测量波形');
hold(ax1, 'on');
[xSorted, sortIdx] = sort(highMap.x_v(:), 'ascend');
plot(ax1, xSorted, result.VFit(sortIdx), '-', ...
    'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.4, 'DisplayName', '主方法拟合');
title(ax1, '(a) 波形拟合', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, '空间坐标 x (mm)', 'FontSize', 9);
ylabel(ax1, '电压 (V)', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
residual = highMap.V_a(:) - result.VFit(:);
plot(ax2, highMap.t_v(:), residual, '.', ...
    'Color', [0.85, 0.33, 0.10], 'MarkerSize', 6, 'DisplayName', '残差');
hold(ax2, 'on');
yline(ax2, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
title(ax2, '(b) 残差', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, '时间 t (s)', 'FontSize', 9);
ylabel(ax2, '残差 (V)', 'FontSize', 9);
legend(ax2, 'Location', 'best');

annotationText = sprintf([ ...
    'g_{true} = %.3f mm, g_{id} = %.3f mm\n', ...
    'A_{true} = [%.3f, %.3f] mm\n', ...
    'f_{true} = [%.1f, %.1f] Hz\n', ...
    'f_{id} = [%.2f, %.2f] Hz\n', ...
    'gap error = %.4f mm, rmse = %.4g, fast = %d'], ...
    truthCase.g_high, result.g_used, ...
    truthCase.A(1), truthCase.A(2), ...
    truthCase.f(1), truthCase.f(2), ...
    result.f_id(1), result.f_id(2), ...
    summaryRow.gap_error_mm, result.rmse, logical(summaryRow.fast_mode_used));

annotation(fig, 'textbox', [0.17, 0.78, 0.52, 0.13], ...
    'String', annotationText, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'Interpreter', 'tex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8);

apply_inv_log_2_figure_style(ax1, 170, 110);
end

function noiseRatio = calibrate_single_case_snr_ratio(cfgCase, gapList, xCell, yCell, domain)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
dataClean = simulate_rotating_waveform_from_template(Fx, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    0, domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, 'noise_ratio');
signal = dataClean.V_clean(:) - min(dataClean.V_clean(:));
signalRms = sqrt(mean(signal.^2));
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, ...
    linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatio = signalRms ./ (templateRange .* 10^(cfgCase.snrDb / 20));
end

function val = get_struct_field(S, name, defaultVal)
if isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end
