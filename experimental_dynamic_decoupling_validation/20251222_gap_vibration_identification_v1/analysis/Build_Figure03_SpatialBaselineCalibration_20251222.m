%% Figure 03: real-waveform spatial baseline calibration
% Uses one real 20251222 low-speed pulse as the extraction data unit.
% Stable tails estimate the baseline. Only the selected above-threshold
% region is fitted and displayed with the super-Gaussian model.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(scriptDir);
addpath(packageRoot);
cfg = Setup_Paths_20251222();

sourceFile = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', 'LowSpeed_Template_SourceData_20251222.mat');
templateFile = cfg.files.lowSpeedTemplate;
rawFile = fullfile(cfg.paths.rawDatasetRoot, cfg.case.lowSpeedCase, '4-1-1000.mat');
outDir = fullfile(cfg.paths.figures, 'spatial_baseline_calibration');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

assert(isfile(sourceFile), 'Missing point-cloud source: %s', sourceFile);
assert(isfile(templateFile), 'Missing low-speed template: %s', templateFile);
assert(isfile(rawFile), 'Missing raw low-speed signal: %s', rawFile);

%% 1. Load real pulse point cloud and select a representative complete pulse
Ssrc = load(sourceFile, 'point_cloud');
Stpl = load(templateFile, 'Template');
pcIndex = find([Ssrc.point_cloud.sensor_id] == 1 & ...
    [Ssrc.point_cloud.blade_id] == 1, 1, 'first');
tplIndex = find([Stpl.Template.SensorBlade.sensor_id] == 1 & ...
    [Stpl.Template.SensorBlade.blade_id] == 1, 1, 'first');
assert(~isempty(pcIndex) && ~isempty(tplIndex), 'CH1/B1 data were not found.');

P = Ssrc.point_cloud(pcIndex);
Tpl = Stpl.Template.SensorBlade(tplIndex);
thresholdV = Tpl.threshold;
pulseIds = unique(P.pulse_index(:), 'stable');
nPulse = numel(pulseIds);
peakV = nan(nPulse, 1);
widthMm = nan(nPulse, 1);
baselineDeltaV = nan(nPulse, 1);
centerTimeS = nan(nPulse, 1);

for ip = 1:nPulse
    mask = P.pulse_index == pulseIds(ip);
    x = P.x_mm(mask);
    v = P.v(mask);
    t = P.t_s(mask);
    valid = isfinite(x) & isfinite(v) & isfinite(t);
    x = x(valid); v = v(valid); t = t(valid);
    peakV(ip) = max(v);
    centerTimeS(ip) = median(t);
    high = v > thresholdV;
    if nnz(high) >= 2
        widthMm(ip) = max(x(high)) - min(x(high));
    end
    left = x < prctile(x, 15);
    right = x > prctile(x, 85);
    baselineDeltaV(ip) = abs(median(v(left), 'omitnan') - ...
        median(v(right), 'omitnan'));
end

score = robust_square_score_local(peakV) + ...
    robust_square_score_local(widthMm) + ...
    robust_square_score_local(baselineDeltaV);
[~, representativeIndex] = min(score);
representativePulseId = pulseIds(representativeIndex);
pulseMask = P.pulse_index == representativePulseId;
xPulse = P.x_mm(pulseMask);
vPulse = P.v(pulseMask);
tPulse = P.t_s(pulseMask);
[xPulse, order] = sort(xPulse(:));
vPulse = vPulse(order);
tPulse = tPulse(order);

% Retain one complete pulse window. Its outer limits include both baselines.
xOuter = [-7.5, 7.5];
insideWindow = xPulse >= xOuter(1) & xPulse <= xOuter(2) & ...
    isfinite(vPulse) & isfinite(tPulse);
xPulse = xPulse(insideWindow);
vPulse = vPulse(insideWindow);
tPulse = tPulse(insideWindow);

calMask = vPulse > thresholdV;
assert(nnz(calMask) > 100, 'Too few above-threshold calibration samples.');
xLeft = min(xPulse(calMask));
xRight = max(xPulse(calMask));
bufferWidthMm = 0.90;
leftBaseEnd = xLeft - bufferWidthMm;
rightBaseStart = xRight + bufferWidthMm;
baseMask = xPulse <= leftBaseEnd | xPulse >= rightBaseStart;
bufferMask = ~calMask & ~baseMask;
baselineV = median(vPulse(baseMask), 'omitnan');

%% 2. Fit super-Gaussian only on above-threshold samples
xCal = xPulse(calMask);
vCal = vPulse(calMask);
[peakObserved, peakIndex] = max(vCal);
initial = [log(max(peakObserved - baselineV, 0.1)), ...
    log(max((xRight - xLeft) / 2, 0.3)), log(2.5 - 1), xCal(peakIndex)];
objective = @(q) mean((vCal - super_gaussian_local(xCal, q, baselineV)).^2);
options = optimset('Display', 'off', 'MaxIter', 1200, ...
    'MaxFunEvals', 5000, 'TolX', 1e-9, 'TolFun', 1e-10);
qFit = fminsearch(objective, initial, options);
xFit = linspace(xLeft, xRight, 900).';
vFit = super_gaussian_local(xFit, qFit, baselineV);
fitU0 = exp(qFit(1));
fitW = exp(qFit(2));
fitBeta = 1 + exp(qFit(3));
fitXc = qFit(4);
rmseCalMv = 1000 * sqrt(mean((vCal - ...
    super_gaussian_local(xCal, qFit, baselineV)).^2));

%% 3. Load the corresponding real continuous low-speed signal
Sraw = load(rawFile, 'jilu01');
rawIndex = Sraw.jilu01(:, 1);
rawVoltage = Sraw.jilu01(:, 2);
rawTimeS = (rawIndex - rawIndex(1)) / cfg.machine.sampleRateHz;
selectedCenterTimeS = median(tPulse);
topHalfSpanS = 0.036;
topMask = rawTimeS >= selectedCenterTimeS - topHalfSpanS & ...
    rawTimeS <= selectedCenterTimeS + topHalfSpanS;
topTimeMs = 1000 * (rawTimeS(topMask) - selectedCenterTimeS);
topVoltage = rawVoltage(topMask);
topPlotStep = max(1, floor(numel(topTimeMs) / 12000));
topTimeMs = topTimeMs(1:topPlotStep:end);
topVoltage = topVoltage(1:topPlotStep:end);
selectedTimeWindowMs = 1000 * ([min(tPulse), max(tPulse)] - selectedCenterTimeS);

%% 4. Draw a 16:9 presentation figure
navy = [18, 55, 116] / 255;
blue = [25, 103, 184] / 255;
lightBlue = [229, 241, 252] / 255;
orange = [232, 104, 24] / 255;
lightOrange = [255, 243, 224] / 255;
lightGray = [241, 243, 246] / 255;
midGray = [0.42, 0.44, 0.47];
dark = [0.12, 0.14, 0.18];

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [1, 1, 33.87, 19.05], 'Visible', 'off', ...
    'Renderer', 'painters');

annotation(fig, 'rectangle', [0.008, 0.01, 0.984, 0.98], ...
    'Color', navy, 'LineWidth', 1.3);
annotation(fig, 'textbox', [0.05, 0.925, 0.90, 0.055], ...
    'String', '3. 无振动工况下的空间基准波形标定', ...
    'FontName', 'Microsoft YaHei', 'FontWeight', 'bold', ...
    'FontSize', 23, 'Color', navy, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none');
annotation(fig, 'textbox', [0.08, 0.885, 0.84, 0.035], ...
    'String', '完整脉冲用于基线与边界定位；仅选定的阈值以上区域参与超高斯拟合', ...
    'FontName', 'Microsoft YaHei', 'FontSize', 11.5, ...
    'Color', [0.22, 0.25, 0.30], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none');
annotation(fig, 'line', [0.02, 0.98], [0.88, 0.88], ...
    'Color', [0.74, 0.82, 0.92], 'LineWidth', 0.8);

% Left upper panel: real continuous signal in time domain
axTop = axes(fig, 'Position', [0.055, 0.595, 0.585, 0.235]);
hold(axTop, 'on');
ylTop = [-0.55, 4.20];
patch(axTop, [selectedTimeWindowMs(1), selectedTimeWindowMs(2), ...
    selectedTimeWindowMs(2), selectedTimeWindowMs(1)], ...
    [ylTop(1), ylTop(1), ylTop(2), ylTop(2)], lightBlue, ...
    'EdgeColor', blue, 'LineStyle', '--', 'FaceAlpha', 0.52);
plot(axTop, topTimeMs, topVoltage, '-', 'Color', dark, 'LineWidth', 0.65);
yline(axTop, thresholdV, '--', 'Color', orange, 'LineWidth', 1.0);
xlim(axTop, [-1000 * topHalfSpanS, 1000 * topHalfSpanS]);
ylim(axTop, ylTop);
xlabel(axTop, '相对时间  t-tc  (ms)', 'FontName', 'Microsoft YaHei', ...
    'Interpreter', 'none');
ylabel(axTop, '电压  U  (V)', 'FontName', 'Microsoft YaHei', ...
    'Interpreter', 'none');
title(axTop, '（a）连续低速实测信号（CH1，多叶片脉冲）', ...
    'FontName', 'Microsoft YaHei', 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left', 'Interpreter', 'none');
text(axTop, mean(selectedTimeWindowMs), 4.02, '选中的同一条完整脉冲', ...
    'FontName', 'Microsoft YaHei', 'FontSize', 9, 'Color', navy, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(axTop, 34.5, thresholdV + 0.09, 'U_{th}', ...
    'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
    'FontSize', 9, 'Color', orange, 'HorizontalAlignment', 'right');
style_axes_local(axTop);
set([axTop.Title, axTop.XLabel, axTop.YLabel], ...
    'FontName', 'Microsoft YaHei', 'Interpreter', 'none');

% Left lower panel: exactly the selected real pulse in spatial coordinates
axPulse = axes(fig, 'Position', [0.055, 0.245, 0.585, 0.250]);
hold(axPulse, 'on');
yl = [-0.55, 4.15];
patch_region_local(axPulse, xOuter(1), leftBaseEnd, yl, lightGray);
patch_region_local(axPulse, leftBaseEnd, xLeft, yl, lightOrange);
patch_region_local(axPulse, xLeft, xRight, yl, lightBlue);
patch_region_local(axPulse, xRight, rightBaseStart, yl, lightOrange);
patch_region_local(axPulse, rightBaseStart, xOuter(2), yl, lightGray);

displayStep = max(1, floor(numel(xPulse) / 900));
displayIndex = 1:displayStep:numel(xPulse);
hMeasured = plot(axPulse, xPulse(displayIndex), vPulse(displayIndex), 'o', ...
    'Color', [0.40, 0.40, 0.40], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 2.6, 'LineWidth', 0.55, ...
    'DisplayName', '完整实测脉冲');
calDisplayIndex = find(calMask);
calDisplayIndex = calDisplayIndex(1:max(1, floor(numel(calDisplayIndex) / 180)):end);
hCal = plot(axPulse, xPulse(calDisplayIndex), vPulse(calDisplayIndex), 'o', ...
    'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 3.0, ...
    'LineWidth', 0.5, 'DisplayName', '阈值以上样本');
hFit = plot(axPulse, xFit, vFit, '-', 'Color', blue, 'LineWidth', 1.8, ...
    'DisplayName', '超高斯拟合');
hRef = [
    yline(axPulse, thresholdV, '--', 'Color', orange, 'LineWidth', 0.9)
    xline(axPulse, xLeft, '--', 'Color', orange, 'LineWidth', 0.9)
    xline(axPulse, xRight, '--', 'Color', orange, 'LineWidth', 0.9)
    xline(axPulse, leftBaseEnd, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.9)
    xline(axPulse, rightBaseStart, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.9)];
set(hRef, 'HandleVisibility', 'off');
xlim(axPulse, xOuter);
ylim(axPulse, yl);
xlabel(axPulse, '局部空间位置  X  (mm)', 'FontName', 'Microsoft YaHei', ...
    'Interpreter', 'none');
ylabel(axPulse, '电压  U  (V)', 'FontName', 'Microsoft YaHei', ...
    'Interpreter', 'none');
title(axPulse, '（b）同一条完整脉冲经 t→X 映射后的空间域标定', ...
    'FontName', 'Microsoft YaHei', 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'left', 'Interpreter', 'none');
legend(axPulse, [hMeasured, hCal, hFit], ...
    {'完整实测脉冲','阈值以上样本','超高斯拟合'}, ...
    'Location', 'northeast', 'Box', 'off', ...
    'FontName', 'Microsoft YaHei', 'FontSize', 7.8);
text(axPulse, mean([xOuter(1), leftBaseEnd]), -0.42, ...
    '基线区', 'FontName', 'Microsoft YaHei', 'FontSize', 8, ...
    'HorizontalAlignment', 'center');
text(axPulse, mean([leftBaseEnd, xLeft]), -0.42, ...
    '缓冲区', 'FontName', 'Microsoft YaHei', 'FontSize', 8, ...
    'HorizontalAlignment', 'center', 'Color', [0.68, 0.35, 0.08]);
text(axPulse, mean([xLeft, xRight]), -0.42, ...
    '阈值拟合区  D_{cal}', 'Interpreter', 'tex', ...
    'FontName', 'Microsoft YaHei', 'FontSize', 8.5, ...
    'HorizontalAlignment', 'center', 'Color', navy);
text(axPulse, mean([xRight, rightBaseStart]), -0.42, ...
    '缓冲区', 'FontName', 'Microsoft YaHei', 'FontSize', 8, ...
    'HorizontalAlignment', 'center', 'Color', [0.68, 0.35, 0.08]);
text(axPulse, mean([rightBaseStart, xOuter(2)]), -0.42, ...
    '基线区', 'FontName', 'Microsoft YaHei', 'FontSize', 8, ...
    'HorizontalAlignment', 'center');
text(axPulse, xLeft, thresholdV + 0.11, 'X_l', ...
    'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
    'FontSize', 9, 'Color', orange, 'HorizontalAlignment', 'right');
text(axPulse, xRight, thresholdV + 0.11, 'X_r', ...
    'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
    'FontSize', 9, 'Color', orange, 'HorizontalAlignment', 'left');
text(axPulse, xOuter(1) + 0.12, 3.95, 'X_a', ...
    'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
    'FontSize', 8.5, 'Color', midGray);
text(axPulse, xOuter(2) - 0.12, 3.95, 'X_b', ...
    'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
    'FontSize', 8.5, 'Color', midGray, 'HorizontalAlignment', 'right');
style_axes_local(axPulse);
set([axPulse.Title, axPulse.XLabel, axPulse.YLabel], ...
    'FontName', 'Microsoft YaHei', 'Interpreter', 'none');

% Full-window bracket under the pulse plot
annotation(fig, 'doublearrow', [0.073, 0.625], [0.180, 0.180], ...
    'Color', navy, 'LineWidth', 1.1, 'Head1Length', 6, 'Head2Length', 6);
annotation(fig, 'textbox', [0.255, 0.150, 0.20, 0.025], ...
    'String', '完整脉冲窗口  Dwin = [Xa, Xb]', ...
    'Interpreter', 'none', 'FontName', 'Microsoft YaHei', ...
    'FontSize', 9, 'Color', navy, 'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');

% Right workflow cards
annotation(fig, 'line', [0.665, 0.665], [0.12, 0.88], ...
    'Color', [0.74, 0.82, 0.92], 'LineWidth', 0.8);
cardX = 0.690;
cardW = 0.275;
cardH = 0.165;
cardY = [0.690, 0.505, 0.320, 0.135];
cardTitle = {
    '1  稳定脉冲筛选与完整窗口截取'
    '2  时间—空间映射与区域划分'
    '3  左右稳定尾部估计基线'
    '4  选定区域超高斯标定'};
cardBody = {
    sprintf(['阈值检测全部候选脉冲；\n按峰值、宽度和基线稳定性筛选，\n' ...
        '选取一条代表性完整脉冲。'])
    sprintf(['X = R Δθ\nDwin = Dbase ∪ Dbuf ∪ Dcal\n' ...
        'Xl、Xr仅表示阈值交点。'])
    sprintf(['Ubias = median{Umeas(Xi), Xi ∈ Dbase}\n' ...
        '基线确定后在参数优化中固定。'])
    sprintf(['仅Dcal内 Umeas > Uth 的样本\n进入普通最小二乘；超高斯曲线\n' ...
        '只在Dcal内计算和绘制。'])};
for ic = 1:4
    annotation(fig, 'rectangle', [cardX, cardY(ic), cardW, cardH], ...
        'Color', [0.53, 0.66, 0.82], 'LineWidth', 0.9, ...
        'FaceColor', [0.985, 0.990, 0.998]);
    annotation(fig, 'textbox', [cardX + 0.012, cardY(ic) + cardH - 0.048, ...
        cardW - 0.024, 0.040], 'String', cardTitle{ic}, ...
        'FontName', 'Microsoft YaHei', 'FontSize', 10.5, ...
        'FontWeight', 'bold', 'Color', navy, 'EdgeColor', 'none', ...
        'VerticalAlignment', 'middle');
    annotation(fig, 'textbox', [cardX + 0.018, cardY(ic) + 0.018, ...
        cardW - 0.036, cardH - 0.070], 'String', cardBody{ic}, ...
        'FontName', 'Microsoft YaHei', 'FontSize', 8.8, ...
        'Color', [0.20, 0.22, 0.26], 'EdgeColor', 'none', ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
        'Interpreter', 'none');
end

% Compact equation and actual fitted values
annotation(fig, 'textbox', [0.690, 0.035, 0.275, 0.082], ...
    'String', sprintf(['Upre = Ubias + U0 exp[-|(X-xc)/w|^β]\n' ...
    '实测拟合：U0 = %.3f V，w = %.3f mm，β = %.2f，xc = %.3f mm\n' ...
    '选定区域RMSE = %.1f mV'], ...
    fitU0, fitW, fitBeta, fitXc, rmseCalMv), ...
    'Interpreter', 'none', 'FontName', 'Microsoft YaHei', ...
    'FontSize', 7.8, 'Color', navy, ...
    'EdgeColor', [0.78, 0.84, 0.92], ...
    'BackgroundColor', [0.97, 0.985, 1], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

annotation(fig, 'textbox', [0.045, 0.035, 0.590, 0.065], ...
    'String', ['图中灰点和黑线均来自20251222低速实测CH1信号；橙色点仅为显示抽样，' ...
    '参数优化使用Dcal内全部样本；缓冲区和基线区不参与拟合，也不进行模型外推。'], ...
    'Interpreter', 'none', 'FontName', 'Microsoft YaHei', ...
    'FontSize', 8.7, 'Color', [0.25, 0.28, 0.32], ...
    'EdgeColor', [0.78, 0.84, 0.92], 'BackgroundColor', [0.97, 0.985, 1], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Margin', 5);

%% 5. Export
pngFile = fullfile(outDir, 'Figure03_SpatialBaselineCalibration_RealWaveform_20251222.png');
pdfFile = fullfile(outDir, 'Figure03_SpatialBaselineCalibration_RealWaveform_20251222.pdf');
emfFile = fullfile(outDir, 'Figure03_SpatialBaselineCalibration_RealWaveform_20251222.emf');
exportgraphics(fig, pngFile, 'Resolution', 220);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
print(fig, emfFile, '-dmeta', '-r300');
savefig(fig, fullfile(outDir, ...
    'Figure03_SpatialBaselineCalibration_RealWaveform_20251222.fig'));
% close(fig);

fprintf('Representative pulse: CH1/B1 pulse %d, t=%.6f s.\n', ...
    representativePulseId, selectedCenterTimeS);
fprintf('Fit: U0=%.6f V, w=%.6f mm, beta=%.6f, xc=%.6f mm.\n', ...
    fitU0, fitW, fitBeta, fitXc);
fprintf('Selected-region RMSE: %.3f mV.\n', rmseCalMv);
fprintf('Saved: %s\n', pngFile);

function score = robust_square_score_local(values)
values = values(:);
center = median(values, 'omitnan');
scale = 1.4826 * median(abs(values - center), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = max(std(values, 'omitnan'), eps);
end
score = ((values - center) ./ scale) .^ 2;
score(~isfinite(score)) = inf;
end

function voltage = super_gaussian_local(x, q, baseline)
u0 = exp(q(1));
w = exp(q(2));
beta = 1 + exp(q(3));
xc = q(4);
voltage = baseline + u0 .* exp(-abs((x - xc) ./ w) .^ beta);
end

function patch_region_local(ax, x1, x2, yLimits, color)
patch(ax, [x1, x2, x2, x1], ...
    [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.72, 'HandleVisibility', 'off');
end

function style_axes_local(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 8.3, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on', ...
    'Layer', 'top', 'XGrid', 'off', 'YGrid', 'off', ...
    'Color', 'none');
end
