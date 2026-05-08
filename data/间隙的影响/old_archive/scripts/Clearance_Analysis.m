%% H-groove Blade: Tip Clearance Analysis (Final Integrated)
%% 核心拟合算法：严格对其提供的 Reference Code
%% 绘图样式：黄金布局 + 指定Hex配色 + 黑色边框 + [改进] 线条对比度增强

clc; clear; close all;

% =========================================================================
%                             1. 基础设置与美学定义
% =========================================================================
fontName = 'Times New Roman'; 
fontSize = 9; 

set(0,'DefaultAxesFontName',fontName,'DefaultAxesFontSize',fontSize);
set(0,'DefaultTextFontName', fontName, 'DefaultTextFontSize', fontSize);
set(0,'DefaultLegendFontName', fontName, 'DefaultLegendFontSize', fontSize);
set(0,'DefaultLineLineWidth',1.2,'DefaultAxesLineWidth',1.0);
set(0,'DefaultLegendInterpreter','tex');

% --- 左图/右上图 配色 (间隙梯度色: 暖->冷) ---
c_case1 = [0.85, 0.35, 0.10]; % d0=0.5mm
c_case2 = [0.85, 0.65, 0.10]; % d0=1.0mm
c_case3 = [0.20, 0.60, 0.60]; % d0=1.5mm
c_case4 = [0.10, 0.25, 0.50]; % d0=2.0mm
colors_case = [c_case1; c_case2; c_case3; c_case4];

% --- 右下角柱状图 配色 (指定 Hex) ---
% Gauss (#E5B9B5)
c_gauss = [0.8980, 0.7255, 0.7098]; 
% BPM (#D7E3BF)
c_bpm   = [0.8431, 0.8902, 0.7490]; 
% APKM (#8DB1E2)
c_apkm  = [0.5529, 0.6941, 0.8863]; 

% --- 物理常量 ---
G.r_phys       = 5.0;  
G.delta_b_phys = 10.0; 
G.hg           = 5.0;
% G.tg 和 G.d0_fixed 将在循环中设定

% --- 文件列表 (间隙分析) ---
file_list = {
    'wave_0.5_叶片8mm_凹槽4mm_间隙0.5mm_加密.txt', 0.5;
    'wave_0.5_叶片8mm_凹槽4mm_间隙1mm_加密.txt',   1.0;
    'wave_0.5_叶片8mm_凹槽4mm_间隙1.5mm_加密.txt', 1.5;
    'wave_0.5_叶片8mm_凹槽4mm_间隙2mm_加密.txt',   2.0
};
num_files = size(file_list, 1);

% 预分配
RMSE_Matrix = zeros(num_files, 3); 
Results = struct('d0', {}, 'x_fit', {}, 'y_fit', {}, 'y_u', {}, 'resid', {}, 'x_norm', {});

% =========================================================================
%                             2. 循环处理
% =========================================================================
for k = 1:num_files
    filename = file_list{k, 1};
    current_d0 = file_list{k, 2};
    
    % 设定物理参数
    G.tg = 4.0;          
    G.d0_fixed = current_d0; 
    
    % --- 数据加载 (降采样逻辑保持) ---
        T = readtable(filename,'VariableNamingRule','preserve');
        M = table2array(T); 
        
        % 仅当间隙为 0.5mm 时进行隔点采样 (数据量大)
        if current_d0 == 0.5
            x_raw = M(1:2:end, 1); 
            y_raw = M(1:2:end, 2);
        else
            x_raw = M(:, 1); 
            y_raw = M(:, 2);
        end

    % --- 阈值筛选 ---
    y_min=min(y_raw); y_max=max(y_raw); 
    mask=(y_raw >= y_min+0.08*(y_max-y_min)); % 保持参考代码的 0.08 阈值
    x_fit=x_raw(mask); y_fit=y_raw(mask);
    
    % =====================================================================
    %              Step 1: 智能 Center 预计算 (严格按照参考代码)
    % =====================================================================
    [pks, locs] = findpeaks(y_fit, x_fit, 'MinPeakProminence', 0.01*(y_max-y_min));

    if length(locs) >= 2
        left_peak = locs(1); right_peak = locs(end);
        mask_valley = (x_fit >= left_peak) & (x_fit <= right_peak);
        [~, min_idx] = min(y_fit(mask_valley));
        x_val_segment = x_fit(mask_valley);
        center_fixed = x_val_segment(min_idx);
    else
        [~, max_idx] = max(y_fit);
        center_fixed = x_fit(max_idx);
    end
    G.center_fixed = center_fixed; 

    % =====================================================================
    %              Step 2: 模型求解 (严格按照参考代码)
    % =====================================================================
    opts = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',5e4,'FunctionTolerance',1e-9);

    % --- 1. BPM ---
    bpm_fun = @(p,x) bpm_fixed_2p(p, x, G);
    A_t = A_total(0, G.r_phys, G.delta_b_phys); A_g = A_groove(0, G.r_phys, G.tg);
    s0 = (A_t-A_g)/G.d0_fixed + A_g/(G.d0_fixed+G.hg);
    eps_est = (max(y_fit)-min(y_fit))/max(s0, 1e-9);
    pB0 = [min(y_fit), eps_est]; 
    lbB = [0, 0]; ubB = [max(y_fit), Inf];
    [pB, ~] = lsqcurvefit(bpm_fun, pB0, x_fit, y_fit, lbB, ubB, opts);
    
    % --- 2. APKM (10 Param) ---
    usahpm_fun = @(p,x) usahpm_10param_model(p, x, G);
    pU0 = zeros(1,10);
    pU0(1:2) = pB; pU0(3) = G.r_phys; pU0(4) = G.delta_b_phys; pU0(5) = G.tg;
    lbU = [0, 0, G.r_phys*0.5, G.delta_b_phys*0.5, G.tg*0.5, -Inf,-Inf,-Inf,-Inf,-Inf];
    ubU = [max(y_fit), Inf, G.r_phys*2.5, G.delta_b_phys*2, G.tg*2, Inf, Inf, Inf, Inf, Inf];
    [pU, ~] = lsqcurvefit(usahpm_fun, pU0, x_fit, y_fit, lbU, ubU, opts);
    
    % --- 3. Gauss ---
    ft_g = fittype('a*exp(-(x-b).^2/c^2)+d','independent','x');
    % StartPoint 参考代码逻辑
    [fit_gauss, ~] = fit(x_fit, y_fit, ft_g, 'StartPoint', [max(y_fit)-min(y_fit), center_fixed, 2, min(y_fit)]);
    
    % --- 记录与计算 ---
    calc_rmse = @(obs, pred) sqrt(mean((obs - pred).^2));
    rmse_g = calc_rmse(y_fit, feval(fit_gauss, x_fit));
    rmse_b = calc_rmse(y_fit, bpm_fun(pB, x_fit));
    rmse_u = calc_rmse(y_fit, usahpm_fun(pU, x_fit));
    
    RMSE_Matrix(k, :) = [rmse_g, rmse_b, rmse_u];
    
    Results(k).d0 = current_d0;
    Results(k).x_fit = x_fit;
    Results(k).y_fit = y_fit;
    Results(k).y_u = usahpm_fun(pU, x_fit);
    Results(k).resid = y_fit - Results(k).y_u;
    Results(k).x_norm = x_fit - center_fixed;
end

% =========================================================================
%                             绘图 (Final Style Dashboard)
% =========================================================================
figure('Name','Tip Clearance Analysis','Position',[150, 150, 700, 350], 'Color', 'w');

% --- 左图: 波形演变 ---
ax1 = subplot(2, 2, [1, 3]); 
set(ax1, 'Position', [0.08, 0.14, 0.44, 0.80]); 
hold on;
for k = 1:num_files
    % [改进]: 计算较浅的颜色用于Marker，深色用于Line，拉开视觉层次
    c_main = colors_case(k,:);
    c_fade = c_main + (1 - c_main) * 0.4; % 混入60%白色，使点变淡
    
    % 数据点 (Marker: 浅色, 稍小)
    plot(Results(k).x_norm, Results(k).y_fit, 'o', ...
        'MarkerSize', 4, 'MarkerEdgeColor', c_fade, ...
        'MarkerFaceColor', 'none', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
    % 拟合线 (Line: 深色, 加粗到2.0)
    [x_sort, idx_sort] = sort(Results(k).x_norm);
    plot(x_sort, Results(k).y_u(idx_sort), '-', ...
        'Color', c_main, 'LineWidth', 1.0, ...
        'DisplayName', sprintf('d0 = %.1f mm', Results(k).d0));
end
xlabel('Position relative to Center (mm)'); ylabel('Capacitance (pF)'); 
title('Effect of Tip Clearance on Waveform'); 
legend('Location','northeast', 'Box', 'off'); 
grid off; box on; 

axis tight; 
final_xlim = xlim; 
yl = ylim; ylim([yl(1)-0.02, yl(2)+0.02]);

% --- 右上: 残差对比 ---
ax2 = subplot(2, 2, 2);
set(ax2, 'Position', [0.60, 0.60, 0.36, 0.34]); 
hold on;
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
for k = 1:num_files
    % 残差线 (保持原样)
    plot(Results(k).x_norm, Results(k).resid, '.-', ...
        'Color', colors_case(k,:), 'LineWidth', 1.0, 'MarkerSize', 10);
end
title(''); xlabel('Position (mm)'); ylabel('APKM Residuals (pF)'); 
xlim(final_xlim); 
ylim([-4e-3, 4e-3]); 
grid off; box on;

% --- 右下: 全模型精度对比 ---
ax3 = subplot(2, 2, 4);
set(ax3, 'Position', [0.60, 0.14, 0.36, 0.34]); 
hold on;

b = bar(RMSE_Matrix, 0.85);

% 设置颜色 + 黑色边框 + 线宽
b(1).FaceColor = c_gauss; b(1).EdgeColor = 'k'; b(1).LineWidth = 0.5;
b(2).FaceColor = c_bpm;   b(2).EdgeColor = 'k'; b(2).LineWidth = 0.5;
b(3).FaceColor = c_apkm;  b(3).EdgeColor = 'k'; b(3).LineWidth = 0.5;

title(''); ylabel('RMSE (pF)');
set(gca, 'XTick', 1:num_files, ...
    'XTickLabel', {'0.5mm', '1.0mm', '1.5mm', '2.0mm'}); 
xlim([0.5, num_files+0.5]); 
max_rmse = max(RMSE_Matrix,[],'all');
% ylim([0, max_rmse * 1.45]); 
ylim([0, 0.08]); 

legend({'Gauss', 'BPM', 'APKM'}, 'Location', 'north', ...
       'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8, 'AutoUpdate', 'off');

% 数值标签
x_offsets = b(3).XEndPoints; 
y_heights = b(3).YEndPoints;
for k = 1:num_files
    text(x_offsets(k), y_heights(k), sprintf('%.1e', RMSE_Matrix(k,3)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 7.5, 'Color', 'k', 'FontWeight', 'bold');
end

grid off; box on;

% =========================================================================
%                             结果输出
% =========================================================================
fprintf('\n=== RMSE Comparison (Gap Analysis) ===\n');
fprintf('%-8s %-12s %-12s %-12s\n', 'Gap(d0)', 'Gauss', 'BPM', 'APKM');
for k=1:num_files
    fprintf('%.1fmm     %.2e     %.2e     %.2e\n', ...
        Results(k).d0, RMSE_Matrix(k,1), RMSE_Matrix(k,2), RMSE_Matrix(k,3));
end

% ================== 模型函数定义 (参考代码原版) ==================
function C = usahpm_10param_model(p, x, G)
    C0=p(1); epsv=p(2); r=p(3); db=p(4); tg=p(5); b = p(6:10);
    x0 = G.center_fixed; d0 = G.d0_fixed; d_groove = d0 + G.hg;
    A_t = A_total(x - x0, r, db);
    A_g = A_groove(x - x0, r, tg);
    shape = (A_t - A_g)./d0 + A_g./d_groove;
    L = max(1e-6, r + 0.5*db); z = (x - x0)./L;
    M = 1 + b(1)*z + b(2)*z.^2 + b(3)*z.^3 + b(4)*z.^4 + b(5)*z.^5;
    C = C0 + epsv * (shape .* M);
end
function C = bpm_fixed_2p(p, x, G)
    C0 = p(1); epsv = p(2);
    x0 = G.center_fixed; d0 = G.d0_fixed; d_groove = d0 + G.hg;
    A_t = A_total(x - x0, G.r_phys, G.delta_b_phys);
    A_g = A_groove(x - x0, G.r_phys, G.tg);
    shape = (A_t - A_g)./d0 + A_g./d_groove;
    C = C0 + epsv * shape;
end
function S = A_total(x, r, delta_b), lo=x-delta_b/2; hi=x+delta_b/2; S=S_strip(lo,hi,r); end
function S = A_groove(x, r, tg), lo=x-tg/2; hi=x+tg/2; S=S_strip(lo,hi,r); end
function S = S_strip(lo, hi, r)
lo_c = max(lo, -r); hi_c = min(hi, r);
S = max(0, Fint(hi_c, r) - Fint(lo_c, r));
end
function A = Fint(t, r)
tc = max(-r, min(r, t)); A = tc.*sqrt(max(0,r.^2-tc.^2)) + r.^2.*asin(tc./r);
end