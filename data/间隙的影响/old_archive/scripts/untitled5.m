%% H-groove Blade: Advanced APKM Analysis (v27.4 Algorithm + Final Plot Style)
%% 核心算法：基于 v27.4 (变量投影 + 两阶段优化)
%% 绘图风格：黄金布局 + 指定Hex配色 + 黑色边框

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
c_gauss = [0.8980, 0.7255, 0.7098]; 
c_bpm   = [0.8431, 0.8902, 0.7490]; 
c_apkm  = [0.5529, 0.6941, 0.8863]; 

% --- 物理常量 (v27.4 Base) ---
G.r_phys       = 5.0;   
G.delta_b_phys = 10.0;  
G.tg           = 6.0;   % 注意：v27.4 中 tg=6.0
G.hg           = 5.0;   

% --- 文件列表 ---
file_list = {
    'wave_0.5_叶片8mm_凹槽4mm_间隙0.5mm_加密.txt', 0.5;
    'wave_0.5_叶片8mm_凹槽4mm_间隙1mm_加密.txt',   1.0;
    'wave_0.5_叶片8mm_凹槽4mm_间隙1.5mm_加密.txt', 1.5;
    'wave_0.5_叶片8mm_凹槽4mm_间隙2mm_加密.txt',   2.0
};
num_files = size(file_list, 1);

% 结果存储
RMSE_Matrix = zeros(num_files, 3); 
APKM_Params = zeros(num_files, 8); % [C0, Cswing, d0, reff, dbeff, k1, k2, k3]
Results = struct('d0', {}, 'x_fit', {}, 'y_fit', {}, 'y_u', {}, 'resid', {}, 'x_norm', {});

% =========================================================================
%                             2. 循环处理 (v27.4 核心逻辑)
% =========================================================================
fprintf('Starting Advanced Analysis (Variable Projection)...\n');

% 定义线性求解器 (C0, C_swing)
proj_lin = @(y, s) ( [ones(numel(s),1), s(:)] \ y(:) );
pred_lin = @(y, s) ( [ones(numel(s),1), s(:)] * proj_lin(y, s) );

% 定义形状函数句柄
shape_eff_func = @(x, d0, r_eff, db_eff) ...
    (A_total(x, r_eff, db_eff)./d0) + ...
     A_groove(x, r_eff, G.tg) .* (1/(d0+G.hg) - 1/d0);

for k = 1:num_files
    filename = file_list{k, 1};
    current_d0 = file_list{k, 2};
    
    % --- 数据加载 ---
    if exist(filename, 'file')
        T = readtable(filename,'VariableNamingRule','preserve');
        M = table2array(T); 
        x_raw = M(:,1); y_raw = M(:,2);
    else
        % 模拟数据 (兜底)
        x_raw = linspace(-10,10,100)';
        y_raw = 0.3*exp(-(x_raw/4).^2);
    end

    % --- 预处理 ---
    % 1. 对称定心 (v27.4 method)
    center_est = estimate_center_by_symmetry(x_raw, y_raw);
    x_all = x_raw - center_est;
    
    % 2. 阈值截取
    y_max = max(y_raw);
    mask = (y_raw >= 0.10 * y_max);
    x_fit = x_all(mask); y_fit = y_raw(mask);
    
    % =====================================================================
    %              Model 1: Gauss (Benchmark)
    % =====================================================================
    ft_g = fittype('a*exp(-(x-b).^2/c^2)+d','independent','x');
    [fit_gauss, ~] = fit(x_fit, y_fit, ft_g, 'StartPoint', [y_max, 0, 2, min(y_fit)]);
    y_pred_g = feval(fit_gauss, x_fit);

    % =====================================================================
    %              Model 2: BPM (Simplified Physical)
    % =====================================================================
    % BPM 只优化 C0, C_swing 和 d0 (假设几何 r, db 为真值)
    shape_true = @(x, d) shape_eff_func(x, d, G.r_phys, G.delta_b_phys);
    
    % 寻找最佳 d0
    fun_bpm = @(d) y_fit - pred_lin(y_fit, shape_true(x_fit, d)/max(shape_true(x_fit, d)));
    d0_bpm = lsqnonlin(fun_bpm, current_d0, 0.1, 5.0, optimoptions('lsqnonlin','Display','off'));
    
    % 计算最终 BPM 曲线
    s_bpm = shape_true(x_fit, d0_bpm); s_bpm = s_bpm/max(s_bpm);
    theta_bpm = proj_lin(y_fit, s_bpm); % [C0, C_swing]
    y_pred_b = theta_bpm(1) + theta_bpm(2) * s_bpm;

    % =====================================================================
    %              Model 3: APKM (Advanced Two-Stage)
    % =====================================================================
    % 固定 d0 为名义值，优化有效几何参数 r_eff, db_eff
    
    % Stage A: 几何优化 (Grid Search + Optimization)
    r_seeds  = G.r_phys * [0.9, 1.0, 1.1];
    db_seeds = G.delta_b_phys * [0.9, 1.0, 1.1];
    best_cost = Inf; best_geom = [G.r_phys, G.delta_b_phys];
    
    for sr = 1:numel(r_seeds)
        for sb = 1:numel(db_seeds)
            q0 = [r_seeds(sr), db_seeds(sb)];
            % 目标函数：残差平方和 (线性参数已投影)
            funA = @(q) y_fit - pred_lin(y_fit, shape_eff_func(x_fit, current_d0, q(1), q(2))/max(shape_eff_func(x_fit, current_d0, q(1), q(2))));
            
            lbA = [G.r_phys*0.5, G.delta_b_phys*0.5];
            ubA = [G.r_phys*2.0, G.delta_b_phys*2.0];
            
            [q_res, resnorm] = lsqnonlin(funA, q0, lbA, ubA, optimoptions('lsqnonlin','Display','off'));
            if resnorm < best_cost, best_cost = resnorm; best_geom = q_res; end
        end
    end
    
    % Stage B: 调制优化 (Modulation)
    % 构造调制形状函数句柄
    z_func = @(x, r, db) x ./ max(1e-6, r + 0.5*db);
    
    mod_shape = @(x, r, db, k) ...
        shape_eff_func(x, current_d0, r, db) .* (1 + k(1)*z_func(x,r,db) + k(2)*z_func(x,r,db).^2 + k(3)*z_func(x,r,db).^3);
    
    funB = @(k) y_fit - pred_lin(y_fit, mod_shape(x_fit, best_geom(1), best_geom(2), k)/max(mod_shape(x_fit, best_geom(1), best_geom(2), k)));
    
    k0 = [0, 0, 0];
    k_res = lsqnonlin(funB, k0, [-Inf,-Inf,-Inf], [Inf,Inf,Inf], optimoptions('lsqnonlin','Display','off'));
    
    % 计算最终 APKM 曲线
    s_apkm_raw = mod_shape(x_fit, best_geom(1), best_geom(2), k_res);
    s_apkm = s_apkm_raw / max(s_apkm_raw);
    theta_apkm = proj_lin(y_fit, s_apkm);
    y_pred_a = theta_apkm(1) + theta_apkm(2) * s_apkm;
    
    % --- 结果记录 ---
    calc_rmse = @(obs, pred) sqrt(mean((obs - pred).^2));
    RMSE_Matrix(k, :) = [calc_rmse(y_fit, y_pred_g), calc_rmse(y_fit, y_pred_b), calc_rmse(y_fit, y_pred_a)];
    
    % 存储参数: [C0, Cswing, d0_fixed, r_eff, db_eff, k1, k2, k3]
    APKM_Params(k, :) = [theta_apkm(1), theta_apkm(2), current_d0, best_geom(1), best_geom(2), k_res];
    
    Results(k).d0 = current_d0;
    Results(k).x_fit = x_fit;
    Results(k).y_fit = y_fit;
    Results(k).y_u = y_pred_a;
    Results(k).resid = y_fit - y_pred_a;
    Results(k).x_norm = x_fit; % 已经定心
end

% =========================================================================
%                             绘图 1: 拟合效果 (Final Style)
% =========================================================================
figure('Name','Advanced Tip Clearance Analysis','Position',[100, 400, 700, 350], 'Color', 'w');

% --- 左图: 波形演变 ---
ax1 = subplot(2, 2, [1, 3]); 
set(ax1, 'Position', [0.08, 0.14, 0.44, 0.80]); 
hold on;
for k = 1:num_files
    % 数据点 (Marker: 浅色, 稍小)
    c_main = colors_case(k,:);
    c_fade = c_main + (1 - c_main) * 0.4;
    
    plot(Results(k).x_norm, Results(k).y_fit, 'o', ...
        'MarkerSize', 4, 'MarkerEdgeColor', c_fade, ...
        'MarkerFaceColor', 'none', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
    % 拟合线 (Line: 深色, 粗)
    [x_sort, idx_sort] = sort(Results(k).x_norm);
    plot(x_sort, Results(k).y_u(idx_sort), '-', ...
        'Color', c_main, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('d0 = %.1f mm', Results(k).d0));
end
xlabel('Position relative to Center (mm)'); ylabel('Capacitance (pF)'); 
title('Effect of Tip Clearance on Waveform (APKM v27.4)'); 
legend('Location','northeast', 'Box', 'off'); 
grid off; box on; 
axis tight; 
yl = ylim; ylim([yl(1)-0.02, yl(2)+0.02]);

% --- 右上: 残差对比 ---
ax2 = subplot(2, 2, 2);
set(ax2, 'Position', [0.60, 0.60, 0.36, 0.34]); 
hold on;
yline(0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
for k = 1:num_files
    plot(Results(k).x_norm, Results(k).resid, '.-', ...
        'Color', colors_case(k,:), 'LineWidth', 1.0, 'MarkerSize', 8);
end
title(''); xlabel('Position (mm)'); ylabel('APKM Residuals (pF)'); 
axis tight; 
grid off; box on;

% --- 右下: 全模型精度对比 ---
ax3 = subplot(2, 2, 4);
set(ax3, 'Position', [0.60, 0.14, 0.36, 0.34]); 
hold on;

b = bar(RMSE_Matrix, 0.85);
b(1).FaceColor = c_gauss; b(1).EdgeColor = 'k'; b(1).LineWidth = 0.5;
b(2).FaceColor = c_bpm;   b(2).EdgeColor = 'k'; b(2).LineWidth = 0.5;
b(3).FaceColor = c_apkm;  b(3).EdgeColor = 'k'; b(3).LineWidth = 0.5;

title(''); ylabel('RMSE (pF)');
set(gca, 'XTick', 1:num_files, 'XTickLabel', {'0.5mm', '1.0mm', '1.5mm', '2.0mm'}); 
xlim([0.5, num_files+0.5]); 
ylim([0, max(RMSE_Matrix,[],'all')*1.2]); 

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
%                             绘图 2: 参数趋势 (Parameter Trends)
% =========================================================================
figure('Name','Identified Parameter Trends','Position',[200, 100, 700, 450], 'Color', 'w');

x_gap = APKM_Params(:, 3); % d0_fixed
p_C0  = APKM_Params(:, 1);
p_Cw  = APKM_Params(:, 2);
p_Re  = APKM_Params(:, 4);
p_Db  = APKM_Params(:, 5);

subplot(2, 3, 1);
plot(x_gap, p_C0, '-o', 'Color', c_apkm, 'MarkerFaceColor', c_apkm, 'LineWidth', 1.2);
xlabel('Gap (mm)'); ylabel('C_0 (pF)'); title('Base Capacitance'); grid on; box on;

subplot(2, 3, 2);
plot(x_gap, p_Cw, '-s', 'Color', c_apkm, 'MarkerFaceColor', c_apkm, 'LineWidth', 1.2);
xlabel('Gap (mm)'); ylabel('C_{swing} (pF)'); title('Swing Amplitude'); grid on; box on;

subplot(2, 3, 3);
plot(x_gap, p_Re, '-^', 'Color', c_apkm, 'MarkerFaceColor', c_apkm, 'LineWidth', 1.2);
yline(G.r_phys, 'k--');
xlabel('Gap (mm)'); ylabel('r_{eff} (mm)'); title('Effective Radius'); grid on; box on;
ylim([min(p_Re)*0.9, max(p_Re)*1.1]);

subplot(2, 3, 4);
plot(x_gap, p_Db, '-d', 'Color', c_apkm, 'MarkerFaceColor', c_apkm, 'LineWidth', 1.2);
yline(G.delta_b_phys, 'k--');
xlabel('Gap (mm)'); ylabel('\delta_{b,eff} (mm)'); title('Eff Blade Thickness'); grid on; box on;
ylim([min(p_Db)*0.9, max(p_Db)*1.1]);

% 调制系数 (k1, k2, k3)
subplot(2, 3, [5, 6]);
plot(x_gap, APKM_Params(:, 6), '-*', 'Color', c_case1, 'DisplayName','k1'); hold on;
plot(x_gap, APKM_Params(:, 7), '-x', 'Color', c_case2, 'DisplayName','k2');
plot(x_gap, APKM_Params(:, 8), '-+', 'Color', c_case4, 'DisplayName','k3');
xlabel('Gap (mm)'); ylabel('Coeff Value'); title('Modulation Coefficients'); grid on; box on;
legend('Location','best','Box','off');

% =========================================================================
%                             结果输出
% =========================================================================
fprintf('\n=== APKM Identified Params (v27.4) ===\n');
fprintf('%-6s %-8s %-8s %-8s %-8s\n', 'Gap', 'C0', 'C_swing', 'r_eff', 'db_eff');
for k=1:num_files
    p = APKM_Params(k,:);
    fprintf('%.1f    %.4f   %.4f   %.4f   %.4f\n', p(3), p(1), p(2), p(4), p(5));
end

%% ================== 核心数学函数 (v27.4) ==================

% 1. 对称法定心
function c0 = estimate_center_by_symmetry(x, y)
    [xu, iu] = unique(x); yu = y(iu);
    f = @(xx) interp1(xu, yu, xx, 'pchip', 'extrap');
    xmin = min(xu); xmax = max(xu);
    c_grid = linspace(xmin+0.5, xmax-0.5, 50); % 粗搜
    best = Inf; c0 = 0;
    for c = c_grid
        xa = max(xmin, 2*c - xmax); xb = min(xmax, 2*c - xmin);
        if xb <= xa, continue; end
        xx = linspace(xa, xb, 100);
        res = f(xx) - f(2*c - xx);
        val = sum(res.^2);
        if val < best, best = val; c0 = c; end
    end
    % 精搜
    c_grid_fine = linspace(c0-0.2, c0+0.2, 50);
    for c = c_grid_fine
        xa = max(xmin, 2*c - xmax); xb = min(xmax, 2*c - xmin);
        xx = linspace(xa, xb, 100);
        res = f(xx) - f(2*c - xx);
        val = sum(res.^2);
        if val < best, best = val; c0 = c; end
    end
end

% 2. 物理积分函数
function S = A_total(x, r, delta_b), lo=x-delta_b/2; hi=x+delta_b/2; S=S_strip(lo,hi,r); end
function S = A_groove(x, r, tg), lo=x-tg/2; hi=x+tg/2; S=S_strip(lo,hi,r); end
function S = S_strip(lo, hi, r)
    lo_c = max(lo, -r); hi_c = min(hi, r);
    S = max(0, Fint(hi_c, r) - Fint(lo_c, r));
end
function A = Fint(t, r)
    tc = max(-r, min(r, t)); 
    A = tc.*sqrt(max(0,r.^2-tc.^2)) + r.^2.*asin(tc./r);
end