# 留一间隙盲测：Gap-Aware VARPRO 程序框架

## 1. 核心目标

本程序用于验证一个更接近实际高速工况的问题：

低速模板库由若干已知间隙无振动波形构成，但高速工况的静态间隙未知，并且该间隙不直接参与模板库构造。程序需要先从高速波形中估计等效静态间隙状态，再利用估计模板辨识振动参数。

因此，验证设计必须严格区分三类间隙：

```matlab
G_train     % 用于构造静态模板库的间隙集合
g_low       % 用于生成低速旋转参考波形的间隙
g_holdout   % 留出不用、作为高速未知间隙真值的间隙
```

其中 `g_holdout` 只允许用于生成高速测试数据和最终误差评估，不能进入模板库，也不能被识别程序直接读取。

---

## 2. 推荐第一版数据划分

`直叶片2mm_不同间隙.txt` 中默认间隙为：

```matlab
gapList = (0.2:0.2:1.4)';
```

第一版推荐：

```matlab
g_low = 0.8;       % 低速参考间隙
g_holdout = 1.0;   % 高速未知间隙，留一盲测
G_train = gapList(gapList ~= g_holdout);
```

理由：

- `1.0 mm` 位于 `0.8 mm` 和 `1.2 mm` 之间，属于插值验证，不是外推验证
- `0.8 mm` 作为低速参考，和高速 `1.0 mm` 不同，可以检验原固定模板法的误差
- 模板库仍有足够间隙点用于构造 `F(x,g)`

第二版可以再增加难度：

```matlab
g_low = 0.6;
g_holdout = 1.2;
```

第三版再做批量留一：

```matlab
for each g_holdout in gapList(2:end-1)
    leave-one-gap-out test
end
```

---

## 3. 新增主脚本

建议新增主脚本：

```matlab
Program_GapAware_VARPRO_Holdout.m
```

该脚本不是替代 `VARPRO_Lightning_Optimization_best.m`，而是在其前面增加“留一间隙模板库 + 未知间隙估计”层，后半段仍然继承变量投影识别思想。

主脚本建议结构：

```matlab
%% Phase A. 读取直叶片多间隙静态波形
%% Phase B. 留一间隙，构造训练模板库
%% Phase C. 生成低速参考旋转波形
%% Phase D. 生成高速未知间隙 + 振动波形
%% Phase E. OPR + alpha_k 高速空间映射
%% Phase F. 估计高速未知静态状态 gHat
%% Phase G. 在 gHat 模板上进行 VARPRO 频率搜索
%% Phase H. 非线性整波形精修
%% Phase I. 固定低速模板法 vs Gap-Aware 法对比
```

---

## 4. Phase A：读取静态多间隙波形

输入文件：

```matlab
dataFile = fullfile(rootDir, '间隙的影响', '直叶片2mm_不同间隙.txt');
```

输出：

```matlab
gapList
xCell
yCell
```

建议函数：

```matlab
[gapList, xCell, yCell] = load_stacked_gap_curves(dataFile);
```

该函数可以直接复用当前多个脚本里的 `load_stacked_curves` 逻辑。文件中每一段曲线对应一个间隙，默认识别为：

```matlab
gapList = (0.2:0.2:1.4)';
```

---

## 5. Phase B：构造留一模板库

关键原则：

`g_holdout` 必须从训练集中删除。

```matlab
holdoutMask = abs(gapList - g_holdout) < 1e-12;
trainIdx = ~holdoutMask;
```

构造共同空间网格：

```matlab
xMin = max(cellfun(@min, xCell(trainIdx)));
xMax = min(cellfun(@max, xCell(trainIdx)));
xGrid = linspace(xMin, xMax, 2001)';
```

构造训练响应面：

```matlab
S_train = zeros(nnz(trainIdx), numel(xGrid));
for ii = 1:nnz(trainIdx)
    ig = trainIds(ii);
    S_train(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
```

构造模板插值器：

```matlab
F = griddedInterpolant({G_train, xGrid}, S_train, 'linear', 'nearest');
```

构造导数模板：

```matlab
dSdx = gradient(S_train, mean(diff(xGrid)));
Fx = griddedInterpolant({G_train, xGrid}, dSdx, 'linear', 'nearest');
```

输出结构：

```matlab
templateLib.gapTrain = G_train;
templateLib.xGrid = xGrid;
templateLib.F = F;
templateLib.Fx = Fx;
templateLib.g_holdout_truth = g_holdout;   % 仅供最终评估
templateLib.g_low = g_low;
```

---

## 6. Phase C：生成低速参考旋转波形

低速参考波形来自 `g_low`。

如果 `g_low` 属于训练集，可直接用训练模板生成：

```matlab
F_low = @(x) templateLib.F(g_low, x);
```

如果想严格使用原始波形，则从 `xCell/yCell` 取对应曲线：

```matlab
idxLow = find(abs(gapList - g_low) < 1e-12, 1);
F_low_true = griddedInterpolant(xCell{idxLow}, yCell{idxLow}, 'pchip', 'nearest');
```

生成低速数据：

```matlab
Data_Low = simulate_rotating_waveform_from_template( ...
    F_low_true, RPM_low, NumRevs_low, fs, R_tip, alpha_k, theta_OPR, ...
    zeros(size(t)), noiseSetting);
```

第一版单等效传感器即可：

```matlab
NumSensors = 1;
alpha_k = alpha_0;
```

后续再扩展到多传感器：

```matlab
alpha_k = Sim.alpha;
```

低速数据用途：

- 提供 OPR
- 提供几何映射验证
- 模拟真实低速参考工况

注意：在该验证程序中，低速数据不再负责构造唯一模板，因为模板库已经来自多间隙静态数据。

---

## 7. Phase D：生成高速未知间隙振动波形

高速真实间隙来自 `g_holdout`，但算法不能知道它。

取留出的原始静态波形：

```matlab
idxHigh = find(abs(gapList - g_holdout) < 1e-12, 1);
F_high_true = griddedInterpolant(xCell{idxHigh}, yCell{idxHigh}, 'pchip', 'nearest');
```

定义振动：

```matlab
A_true = [0.30, 0.20];       % mm
f_true = [500, 1300];        % Hz
phi_true = [pi/4, -pi/3];

u_t = A_true(1) * sin(2*pi*f_true(1)*t + phi_true(1)) + ...
      A_true(2) * sin(2*pi*f_true(2)*t + phi_true(2));
```

生成高速波形：

```matlab
Data_High = simulate_rotating_waveform_from_template( ...
    F_high_true, RPM_high, NumRevs_high, fs, R_tip, alpha_k, theta_OPR, ...
    u_t, noiseSetting);
```

保存真值：

```matlab
truth.g_high = g_holdout;
truth.A = A_true;
truth.f = f_true;
truth.phi = phi_true;
truth.u_t = u_t;
```

识别阶段只允许使用：

```matlab
Data_High.t
Data_High.V_cap
Data_High.V_OPR
templateLib.F
templateLib.Fx
```

---

## 8. Phase E：高速空间映射

这一层保留 `VARPRO_Lightning_Optimization_best.m` 的核心几何：

```matlab
T_e = T_opr_h(m) + alpha_k(k) / Omega_m;
x_local = Omega_m * R_tip * (t - T_e);
```

输出：

```matlab
highMap.t_v
highMap.V_a
highMap.x_v
highMap.S_v
highMap.rev_v
```

关键修改：

此阶段不要直接计算：

```matlab
V_exp = V_fp{k}(x_v);
K_exp = K1_fp{k}(x_v);
DV = V_a - V_exp;
```

因为高速静态间隙未知。应该等 `gHat` 估计完成后，再计算：

```matlab
V_hat = templateLib.F(gHat, x_v);
K1_hat = templateLib.Fx(gHat, x_v);
DV = V_a - V_hat;
```

---

## 9. Phase F：估计高速未知间隙

第一版使用模板匹配估计：

```matlab
gapQueryGrid = linspace(min(templateLib.gapTrain), max(templateLib.gapTrain), 201);
dxGrid = linspace(-0.8, 0.8, 161);
```

代价函数：

```matlab
J_static(g) = min_dx sum_i w_i * (V_a(i) - F(g, x_v(i)-dx))^2
```

伪代码：

```matlab
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    for id = 1:numel(dxGrid)
        dx = dxGrid(id);
        V_pred = templateLib.F(g, highMap.x_v - dx);
        res = highMap.V_a - V_pred;
        J2D(ig, id) = mean(w .* res.^2);
    end
    [J_static(ig), bestDxIdx] = min(J2D(ig, :));
    dxStatic(ig) = dxGrid(bestDxIdx);
end

[~, idxBest] = min(J_static);
staticState.gHat = gapQueryGrid(idxBest);
staticState.dx0 = dxStatic(idxBest);
staticState.costCurve = J_static;
```

这里 `dx0` 是辅助对齐量，只用于防止小幅振动污染静态间隙估计，不作为最终振动位移输出。

第二版可加入等电平宽度辅助：

```matlab
widthObs = extract_equal_level_widths(highMap.V_a, highMap.x_v, levelList);
gByWidth = match_width_signature(widthObs, templateLib.widthSignature);
```

然后用 `gByWidth` 限制 `gapQueryGrid`：

```matlab
gapQueryGrid = linspace(gByWidth-0.2, gByWidth+0.2, 101);
```

---

## 10. Phase G：Gap-Aware VARPRO 输入构造

固定 `gHat`：

```matlab
V_hat = templateLib.F(staticState.gHat, highMap.x_v);
K1_hat = templateLib.Fx(staticState.gHat, highMap.x_v);
DV = highMap.V_a - V_hat;
```

构造变量投影字典：

```matlab
for i = 1:length(f_grid)
    fi = f_grid(i);
    Dict(:, 2*i-1) = -K1_hat .* sin(2*pi*fi * highMap.t_v);
    Dict(:, 2*i) = -K1_hat .* cos(2*pi*fi * highMap.t_v);
end
```

后续的频率峰值搜索、候选频率组合和线性系数估计可以复用 `VARPRO_Lightning_Optimization_best.m` 中的 Phase 1 和 Phase 1.5。

---

## 11. Phase H：非线性整波形精修

新增残差函数：

```matlab
function r = residual_calc_6D_gapaware(p, t, V, x, templateLib, gHat)
    u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
        p(4) * sin(2*pi*p(6)*t + p(5));

    V_sim = templateLib.F(gHat, x - u);
    r = V - V_sim;
end
```

优化入口：

```matlab
[p_f, cost_f] = lsqnonlin( ...
    @(p) residual_calc_6D_gapaware(p, highMap.t_v, highMap.V_a, ...
    highMap.x_v, templateLib, staticState.gHat), ...
    p0, lb, ub, opts_lsq);
```

第一版不建议同时优化 `g` 和 `p`。原因是 `g` 和低频/静态位移可能互相抢解释权，导致验证不清楚。

第二版可做小范围联合精修：

```matlab
q = [p, g];
g in [gHat - 0.05, gHat + 0.05]
```

但应该放在独立对比项里。

---

## 12. 必须设置的对照组

为了证明间隙影响确实被消除，同一组高速盲测数据至少跑三种方法。

### 方法 1：固定低速模板

```matlab
V_sim = F(g_low, x - u(t;p));
```

这代表原程序隐含假设：

```matlab
g_high = g_low
```

### 方法 2：Gap-Aware 两阶段法

```matlab
gHat = estimate_static_gap(...)
V_sim = F(gHat, x - u(t;p));
```

这是推荐主方法。

### 方法 3：真值模板上限

```matlab
V_sim = F_high_true(x - u(t;p));
```

这不是实际算法，只作为性能上限，用于判断误差来自模板库插值还是来自 VARPRO 识别。

---

## 13. 输出指标

建议输出：

```matlab
Result.truth.g_high
Result.truth.A
Result.truth.f
Result.truth.phi

Result.fixed.g_used
Result.fixed.A_id
Result.fixed.f_id
Result.fixed.phi_id
Result.fixed.rmse

Result.gapAware.gHat
Result.gapAware.g_error
Result.gapAware.A_id
Result.gapAware.f_id
Result.gapAware.phi_id
Result.gapAware.rmse

Result.oracle.A_id
Result.oracle.f_id
Result.oracle.rmse
```

核心评价：

```matlab
gap_error = abs(gHat - g_holdout);
freq_error = abs(sort(f_id) - sort(f_true));
amp_error = abs(sort(A_id) - sort(A_true));
rmse_ratio = Result.gapAware.rmse / Result.fixed.rmse;
```

---

## 14. 推荐图表

第一版至少画 4 张图：

1. 模板库训练间隙与留出间隙示意
2. `J_static(g)` 代价曲线，标出 `gHat` 和 `g_holdout`
3. 固定模板法与 Gap-Aware 法的高速波形拟合对比
4. 真实振动 `u(t)` 与识别振动 `uHat(t)` 对比

图名建议：

```matlab
fig1_template_holdout.png
fig2_static_gap_cost.png
fig3_waveform_fit_comparison.png
fig4_vibration_reconstruction.png
```

---

## 15. 最小可运行版本的范围

第一版为了快速闭环，可以只做：

- 单等效传感器
- `alpha_k = alpha_0`
- `a_k = 1`
- `b_k = 0`
- 高速静态间隙恒定
- 双频振动
- `g_holdout` 只取一个内部间隙

等单传感器跑通后，再扩展：

- 多传感器
- 多个 `g_holdout` 批量留一
- 等电平宽度辅助间隙估计
- `gHat` 小范围联合精修
- 实验数据接入

---

## 16. 一句话程序逻辑

从 `直叶片2mm_不同间隙.txt` 中留出一条未知间隙波形作为高速真值，其余间隙构造 `F(x,g)` 模板库；用另一条训练间隙波形模拟低速参考；高速数据只通过 OPR 和 `alpha_k` 映射到空间坐标，然后先估计未知静态间隙 `gHat`，再固定 `gHat` 使用原 VARPRO 流程识别振动参数，并与固定低速模板法进行对比。
