# GapAware_Holdout 方法求解思路说明

## 1. 问题定义

当前要解决的问题不是普通的振动参数辨识，而是：

> 当高速工况下的静态间隙状态与低速参考工况不同，且该高速静态间隙未知时，如何先消除静态间隙影响，再识别动态振动参数。

高速观测模型写成：

\[
V^H(t)=F(x(t)-u(t),g_H)+\varepsilon(t)
\]

其中：

- \(F(x,g)\) 为静态间隙为 \(g\) 时的无振动空间模板
- \(x(t)\) 为无振动参考空间坐标
- \(u(t)\) 为动态振动位移
- \(g_H\) 为高速真实静态间隙
- \(\varepsilon(t)\) 为噪声与模型误差

传统固定模板法默认：

\[
g_H=g_L
\]

即默认高速静态间隙与低速参考间隙相同。这样一来，静态间隙变化引起的模板失配会被误判为振动，进而造成频率、幅值和整波形拟合结果偏差。

因此，本方法的目标是实现：

\[
\text{静态间隙影响} \quad \text{与} \quad \text{动态振动影响}
\]

的分层解耦。

---

## 2. 总体思路

当前 `gapaware_holdout_pipeline` 的主线程已经切换为：

```text
Step 3: raw-scan 估计静态 gap 初值与静态平移初值
Step 4: 基于初值做 joint_wide_free 联合精修
Step 5: 输出最终 gap-aware 结果与对比图
Step 6: 分析 holdout 位置和模板库稀疏度影响
```

对应的求解逻辑可以写成：

\[
\text{模板库 }F(x,g)
\rightarrow
(\hat g_0,\hat d_{x,0})
\rightarrow
\hat p_0
\rightarrow
(\hat g,\hat d_x,\hat p)
\]

这里：

- \((\hat g_0,\hat d_{x,0})\) 是 `Step 3` 的静态初值
- \(\hat p_0\) 是在 \(\hat g_0\) 上通过 VARPRO 得到的双频振动初值
- \((\hat g,\hat d_x,\hat p)\) 是最终 `joint_wide_free` 联合精修结果

也就是说，当前主线程已经不是“只估计 \(g\) 再固定它识别振动”，而是：

> 先用 raw-scan 获得可靠静态初值，再在小邻域内把间隙、静态平移和振动参数一起联合优化。

---

## 3. 盲测验证问题的构造

总入口：

[Run_All_GapAware_Holdout.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Run_All_GapAware_Holdout.m)

当前不是直接上实验高速数据，而是先做一个严格的“留一间隙盲测”验证。

使用静态数据文件：

`直叶片2mm_不同间隙.txt`

将其中的多条静态无振动波形划分为三类：

- `G_train`：用于构造静态模板库的训练间隙集合
- `g_low`：低速参考工况对应的静态间隙
- `g_holdout`：故意留出、不参与建库、作为高速未知真实间隙的目标值

默认设置为：

```matlab
g_low = 0.8;
g_holdout = 1.0;
```

这意味着：

- `1.0 mm` 不进入训练模板库
- 高速波形真实来自 `1.0 mm`
- 算法只能通过扫描和插值去估计该未知间隙

这样就把问题严格变成了一个“高速静态间隙未知”的盲测识别任务。

---

## 4. 求解步骤

## 4.1 Step 1：构造多间隙静态模板库

文件：

[Step_1_Build_Holdout_TemplateLibrary.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_1_Build_Holdout_TemplateLibrary.m)

该步骤从训练间隙集合中构造：

\[
F(x,g), \qquad F_x(x,g)
\]

其中：

- \(F(x,g)\) 为静态模板面
- \(F_x(x,g)\) 为模板对空间坐标的一阶导数

程序上，这一步输出的不是单一模板，而是一张可插值的二维静态模板面。

---

## 4.2 Step 2：构造低速参考和高速未知间隙数据

文件：

[Step_2_Simulate_Holdout_Data.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_2_Simulate_Holdout_Data.m)

这一步做两件事：

1. 用 `g_low` 生成低速无振动参考数据
2. 用 `g_holdout` 生成高速含双频振动数据

高速观测模型为：

\[
V^H(t)=F(x(t)-u(t),g_{holdout})+\varepsilon(t)
\]

当前主线程是双频振动辨识，振动位移模型为：

\[
u(t)=A_1\sin(2\pi f_1 t+\phi_1)+A_2\sin(2\pi f_2 t+\phi_2)
\]

---

## 4.3 Step 3：空间映射与静态初值 raw-scan

文件：

[Step_3_Map_And_Estimate_StaticGap.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_3_Map_And_Estimate_StaticGap.m)

首先利用 OPR 和传感器安装相位 \(\alpha_k\)，把高速时间波形映射到空间域：

\[
T_e=T_{OPR,m}+\frac{\alpha_k}{\Omega_m}
\]

\[
x=\Omega_m R_{tip}(t-T_e)
\]

从而得到高速空间采样点：

\[
\{x_i,V_i^H\}
\]

随后执行静态初值扫描：

\[
(\hat g_0,\hat d_{x,0})
=
\arg\min_{g,d_x}
\sum_i
\left[
V_i^H-F(x_i-d_x,g)
\right]^2
\]

这里：

- \(\hat g_0\) 是高速静态间隙初值
- \(\hat d_{x,0}\) 是静态平移初值

当前程序已经允许 gap 搜索在训练区间外侧小范围外扩，因此边界 holdout 工况不再被硬性截断在训练端点。

这一步的定位不是得到最终结果，而是给后续联合精修提供一个可靠静态初值。

---

## 4.4 Step 4 第一层：固定 \(\hat g_0\) 做 VARPRO 初值辨识

文件：

[Step_4_Run_VARPRO_Comparison.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_4_Run_VARPRO_Comparison.m)

在 `joint_wide_free` 主线程中，程序不会直接从 `Step 3` 跳到最终结果，而是先在 \(\hat g_0\) 上做一次标准 gap-aware 振动识别，得到振动初值：

\[
\hat p_0=[A_1,\phi_1,f_1,A_2,\phi_2,f_2]
\]

具体过程是：

1. 构造静态模板和导数模板

\[
V_0=F(x,\hat g_0),\qquad K_1=F_x(x,\hat g_0)
\]

2. 构造线性化残差

\[
\Delta V=V^H-V_0
\]

3. 在双频网格上做 VARPRO 粗搜索
4. 得到双频振动参数的非线性精修初值

这一步得到的 `rawScanResult` 是联合精修前的动态初值解。

---

## 4.5 Step 4 第二层：joint_wide_free 联合精修

这一步是当前主线程最重要的改动。

当前正式的 `gapAwareResult` 不再等于 `rawScanResult`，而是使用：

\[
(\hat g,\hat d_x,\hat p)
=
\arg\min_{g,d_x,p}
\sum_i
\left[
V_i^H-F(x_i-d_x-u(t_i;p),g)
\right]^2
\]

但该优化不是全局无约束进行，而是以上一步的初值为中心，在一个小邻域内进行局部联合精修：

\[
g \in [\hat g_0-0.05,\ \hat g_0+0.05]\ \text{mm}
\]

\[
d_x \in [\hat d_{x,0}-0.05,\ \hat d_{x,0}+0.05]\ \text{mm}
\]

\[
f_1 \in [f_{1,0}-20,\ f_{1,0}+20]\ \text{Hz}
\]

\[
f_2 \in [f_{2,0}-20,\ f_{2,0}+20]\ \text{Hz}
\]

振幅在初值附近做比例约束，两个相位保持在：

\[
\phi_1,\phi_2 \in [-\pi,\pi]
\]

这就是当前所说的 `joint_wide_free`：

- `joint`：同时优化 \(g,d_x,p\)
- `wide`：相对不是特别窄的局部盒约束
- `free`：不额外给 \(g,d_x\) 加软惩罚项

它的作用不是推翻前面的分层流程，而是：

> 以 raw-scan 和 VARPRO 已得到的可靠初值为起点，在局部邻域内释放静态与动态参数的耦合自由度，从而进一步降低残差、改善振动参数识别。

---

## 4.6 Step 5：结果可视化与总表输出

文件：

[Step_5_Make_Figures.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_5_Make_Figures.m)

当前绘图层已经明确区分三类结果：

- `fixed low gap`
- `raw-scan initial`
- `gap-aware (joint refined)`
- `oracle`

因此，程序现在不仅能看到最终 `joint_wide_free` 结果，还能直观看到：

- raw-scan 初值结果
- 联合精修后的提升幅度

---

## 4.7 Step 6：模板库覆盖与稀疏度分析

文件：

[Step_6_Analyze_Library_Influence.m](E:/0小论文+程序/0博士期间小论文+程序/间隙和振动解耦/gapaware_holdout_pipeline/Step_6_Analyze_Library_Influence.m)

当前 `Step 6` 也已经切换为联合精修主线程。

也就是说，在 holdout 位置分析和稀疏库分析中：

- `gHat_raw_mm` 记录 `Step 3` 的静态初值
- `g_used_final_mm` 记录 `joint_wide_free` 最终使用的 gap
- 后续频率、幅值和 RMSE 统计都基于联合精修后的 `gapAwareResult`

这样可以同时分析：

1. raw-scan 初值本身有多准
2. 联合精修后最终 gap 会不会进一步调整
3. 模板库稀疏度和边界 holdout 对最终结果的影响

---

## 5. 当前主线程与其他方法的关系

当前 `Step 4` 中保留的对比方法是：

### 5.1 fixed_low_gap

直接使用低速参考 gap：

\[
g=g_{low}
\]

不做间隙修正。它代表传统固定模板法，是错误基线。

### 5.2 gap_aware

这已经不再是早期版本里“固定 \(\hat g_0\) 后只优化振动参数”的方法。

当前 `gap_aware` 的正式含义是：

> 先用 raw-scan 得到 \((\hat g_0,\hat d_{x,0})\)，再用 `joint_wide_free` 在小邻域内联合精修 \(g,d_x,p\)。

### 5.3 oracle_holdout

直接把真实高速 gap 告诉算法，只作为理论上限参考。

---

## 6. 为什么主线程改成 joint_wide_free

前一阶段的多工况、多信噪比对比已经表明：

- `raw_scan_only` 在 gap 误差上更稳
- `joint_wide_free` 在平均频率误差、平均幅值误差和平均 RMSE 上更优

因此，如果当前目标是：

> 让最终振动参数辨识和整波形拟合结果更优

那么把主线程切换为 `joint_wide_free` 是合理的。

当前推荐的理解方式是：

- `Step 3 raw-scan`：负责给出可信静态初值
- `Step 4 joint_wide_free`：负责把静态初值与振动初值一起联合打磨到更优

也就是说，`raw_scan_only` 现在是初始化层，`joint_wide_free` 才是正式输出层。

---

## 7. 一句话总结

当前 `GapAware_Holdout` 主线程的完整求解思路是：

> 先用多间隙静态模板库建立 \(F(x,g)\)，再将高速波形映射到空间域，通过 raw-scan 获得静态间隙和静态平移初值 \((\hat g_0,\hat d_{x,0})\)，随后在 \(\hat g_0\) 上通过 VARPRO 获得双频振动初值 \(\hat p_0\)，最后采用 `joint_wide_free` 在小邻域内对 \((g,d_x,p)\) 进行联合非线性精修，从而实现静态间隙变化与动态振动参数的分层识别和最终联合优化。
