# R4 盲态未知间隙动态辨识：独立改进设计

## 1. 目的与边界

本目录是 R4 冻结程序之外的独立改进支线。它回答的问题是：在低速无振动、高速有振动，且当前低速与高速绝对间隙都不提供给反演器时，能否仅由两段波形估计

\[
g_L,\quad g_H,\quad \Delta g=g_H-g_L,
\]

并在高速波形中分离振动。

`Main_04_TransferredSurfaceDynamicValidation.m` 保留为历史转移面接口验证。当前 `Main_06_NoiselessDynamicGate.m` 已改为盲 EO 动态 gate；它不再接收真值 EO 对、竞争 EO 或真值附近振动初值，正式盲态未知间隙结果以 `Main_09_BlindGapOnlyVPDynamicValidation.m` 和 Main06 的盲 gate 输出为准。

## 2. 受控真值与盲态输入必须分离

真值间隙不随机取样，也不由估计器选取。每个案例从已有 FE 间隙网格中预先指定：

\[
g_L^{\rm true},g_H^{\rm true}\in\{0.5,0.7,0.9,1.1,1.3,1.5\}\ {\rm mm}.
\]

首批案例应只选取校准范围内部的组合，例如 `(0.7,0.9)`、`(0.9,1.1)`、`(1.1,0.9)` 和 `(1.3,1.1)` mm；每一组合同时覆盖正、负和零附近的间隙变化。倾角也从已有的 `0.5:0.5:3.5 deg` FE 面中预先指定，并实行 leave-one-tilt-out。

波形生成允许读取真值：

\[
y_L(x)=V_{\rm FE}(x;z^{\rm true},g_L^{\rm true}),
\]

\[
y_H(t)=V_{\rm FE}(x(t)-u(t);z^{\rm true},g_H^{\rm true}),
\]

其中低速没有振动，`u(t)` 仅在高速生成阶段出现。随后应构造一个没有 `truth` 字段的 `observation`：

```matlab
observation = struct('lowWave', lowWave, 'highWave', highWave, ...
    'sampling', sampling, 'eoPair', prescribedEoPair);
```

估计函数只能接收 `observation`、留出目标面后的 R4 流形以及预先固定的物理边界。`truth` 只能传给最终的评价函数。特别禁止将 `truth.gLow` 或 `truth.gHigh` 用于初值、搜索中心、边界收缩、候选筛选、模型缓存或振动初值。

## 3. 低速：未知间隙的流形定位

R4 原有的已知 anchor 公式应替换为

\[
(\hat z_k,\hat g_{L,k},\hat d_{L,k})
\in\operatorname{localmin}_{z,g,d}
J_L(z,g,d),
\]

\[
J_L(z,g,d)=\frac{1}{N_L}\sum_i\left[y_L(x_i)-V_{\rm R4}(x_i-d;z,g)\right]^2.
\]

其中 `z` 是 PC1 系数流形坐标，`g` 是低速绝对间隙，`d` 是低速空间配准偏移。第一版仿真可令真值 `d_L=0`，但估计器仍应搜索有限的、预先规定的范围，例如 `[-0.10,0.10] mm`，以避免把同域恒等配准误写成算法前提。

现有 `Main_04_BoundaryIdentifiability.m` 的 `J(z,g)` 已实现了核心的二维未知 anchor-gap profile；新程序应复用其数学形式，但不复用其仅针对端点面、仅生成图、且不含偏移和高速联合排序的脚本结构。

推荐的数值顺序：

1. 用 leave-one-tilt-out 训练面拟合 gap-law 与 PC1 系数流形；
2. 在整个预设 `z x g x d` 网格计算 `J_L`；
3. 每个 `z-g` 邻域只保留一个局部谷值，防止相邻网格点重复占用候选名额；
4. 对前 `K=5` 个、且在 `z/g` 上彼此分离的谷值做有界连续最小二乘精修；
5. 保留全部精修候选以及其低速残差，不在低速阶段强制只选一个解。

搜索网格必须由 FE 可支持的范围预先确定，例如 `gSearch=0.5:0.01:1.5 mm`；不是以当前真值为中心的 `gTrue +/- h`。

## 4. 高速：在每个低速候选上拟合未知间隙与振动

第一版固定 EO 阶次，以隔离“未知间隙”而不是“频率搜索”带来的困难。对每个低速候选 `k`，求解

\[
\min_{g_H,d_H,a_1,b_1,a_2,b_2}J_{H,k},
\]

\[
u(\theta)=a_1\sin(EO_1\theta)+b_1\cos(EO_1\theta)
+a_2\sin(EO_2\theta)+b_2\cos(EO_2\theta),
\]

\[
J_{H,k}=\frac{1}{N_H}\sum_j\left[y_H(t_j)-V_{\rm R4}(x_j-d_H-u(\theta_j);\hat z_k,g_H)\right]^2.
\]

约束使用训练库全局范围：`gH in [min(gapTrain), max(gapTrain)]`，`dH in [-0.20,0.20] mm`，振幅使用与 R4 既有动态程序一致但不围绕真值的绝对上界。初值应来自预定义的全局多起点或低速候选的解析零振动投影，不能使用 `gHigh`、真实振幅或真实相位。

高速阶段暂不允许重新训练 PC1，也不允许把高波形直接当作另一个静态 anchor。`z` 由低速无振动波形确定；高速只检验该静态几何状态能否支持间隙与振动解耦。后续若需要处理低速到高速的几何状态变化，应把它明确定义为新的物理参数，而不能静默地让 `z` 漂移。

## 5. 联合候选选择与可辨识性结论

对所有低速候选的高速精修结果，用量纲一致的波形平方残差联合排序：

\[
J_{\rm joint,k}=N_LJ_{L,k}+N_HJ_{H,k}.
\]

输出首选解，以及第二优解的相对目标间隔：

\[
m=\frac{J_{\rm joint,2}-J_{\rm joint,1}}{\max(J_{\rm joint,1},\epsilon)}.
\]

不应只因最优解接近仿真真值就宣布成功。若第二优的间隙与 `z` 参数实质不同而 `m` 小于预先注册阈值，则标记为 `ambiguous_joint`。低速候选本身不分离时标记为 `ambiguous_low_speed`。

此外，在最终候选处数值计算低速灵敏度的归一化相关系数：

\[
\rho_{zg}=\operatorname{corr}\left(\frac{\partial V}{\partial z},\frac{\partial V}{\partial g}\right).
\]

`abs(rho_zg)` 接近 1 时，低速单波形没有足够信息将流形形状与间隙分开；应报告此限制，而非通过真值初始化掩盖它。

## 6. 输出和验收顺序

新脚本将只写入：

```text
blind_unknown_gap/output/<run_id>/
  case_results.csv
  low_profile_candidates.csv
  joint_profile_candidates.csv
  profile_surfaces.mat
  identifiability_summary.csv
  run_configuration.mat
```

`case_results.csv` 至少含有 `gLow_true/gLow_hat`、`gHigh_true/gHigh_hat`、`deltaG_true/deltaG_hat`、`z_hat`、低/高速 RMSE、联合目标间隔、`rho_zg` 和 `identifiability_status`。真值列只在此评价产物中出现。

实施与验证按以下层级推进：

1. 无噪声、固定 EO、零空间偏移；确认真值工况没有传入估计接口；
2. 同一真值表加入受控空间偏移，验证 `d_L/d_H` 不把间隙吸收；
3. 加入多个 SNR 重复；
4. 仅在前三步通过后，再增加未知 EO、单频/双频模型选择。

最小 smoke test 应覆盖至少一个增隙和一个减隙案例，并明确检查：估计器函数签名和输入结构中无 `truth`、无 `anchorGap`、无 `gLowTrue/gHighTrue` 字段。该检查能发现“真值泄漏”；若失败，结果不写入正式汇总。
