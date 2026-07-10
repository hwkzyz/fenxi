# 静态 dg 与残余 tilt 耦合的理论结论

## 1. 先纠正一个判断

如果 CH6 只有等效静态偏移 `dg6`，`gap_tilt` 的正确退化关系应当是：

```text
gap_tilt_residual(dmu = 0) == gap_fixed_tilt / 当前 gap_only
```

它不一定等于当前 `fixed`。因为当前 `fixed` 是“不允许高速 dg、不允许高速 dmu”的模型；
只有当 `fixed` 所用静态模板本身已经包含真实 `dg6` 时，它才会和纯 `dg6` 情况一致。

因此更准确的判据是：

```text
纯 dg6:
gap_only 应恢复 A_true；
gap_tilt 如果允许 dmu，也应该选择 dmu -> 0，并退回 gap_only。
```

## 2. 为什么当前 gap_tilt 会把 A 拉低

当前残余 tilt 的有效间隙可以写成：

```text
g_eff = g0 + dg + (mu0 + dmu) * (x_eval - tau)
x_eval = X - dx - A sin(EO theta + phi) - eta
```

所以 `dmu` 对电压的修正近似包含：

```text
dV_mu ~ F_g * dmu * (X - dx - tau - eta)
        - F_g * dmu * A sin(EO theta + phi)
```

第二项就是 EO 同频项。也就是说，`dmu` 不是一个纯静态参数；一旦它乘在含有
`u = A sin(...)` 的 `x_eval` 上，它就天然可以产生与振动同频的修正。优化器于是会在
`A / phi / dx / dg / dmu` 之间交易误差：

```text
减小 A + 调整 dmu + 调整 dg
```

这会降低电压残差，但得到的振动幅值不再是独立可解释的主振动幅值。

## 3. 仿真诊断结果说明什么

新建研究文件夹中的 Step04 已经验证：

```text
纯 dg6 clean 情况下，真实解存在：
A = 0.37109 mm, dg6 = 0.05000 mm, dmu6 = 0
```

从 `gap_only` 的零 `dmu` 点启动 `gap_tilt`，也能回到近似同一个解：

```text
A ~= 0.37116-0.37127 mm, dg6 ~= 0.049 mm, dmu6 ~= 0
```

但独立多初值的非线性 `gap_tilt` 会落入较差局部盆地：

```text
A ~= 0.343-0.357 mm
dg6 被压小
dmu6 出现非零补偿
```

所以问题不是“纯 dg6 物理上必须导致 gap_tilt 低估”，而是：

```text
残余 dmu 自由度让目标函数出现不良局部极小值和参数不可辨识通道。
```

## 4. 耦合的解决原则

主结果不要让高速自由 `dmu` 参与振动幅值辨识。主线模型应改成：

```text
gap_fixed_tilt:
  允许 dg_s
  固定 mu_s = 低速标定 corr.muGapPerXMm
  不允许高速自由 dmu_s
```

当前代码中它基本就是 `gap_only`。因此论文主结果建议报告 `gap_only`，并在命名上改为
`gap_fixed_tilt`，避免被误解成“没有考虑 tilt”。

当前 `gap_tilt` 应改名为：

```text
gap_tilt_residual
```

只作为诊断/对照，不作为主振动幅值。

## 5. 如果必须保留 residual dmu

应当只让 `dmu` 解释主模型解释不了的残差。线性化写法：

```text
J_q = [J_A, J_phi, J_dx, J_dg]
J_mu_perp = (I - P_Jq) J_mu
```

其中 `P_Jq` 是加权投影到主参数子空间的投影算子。然后只用
`J_mu_perp` 拟合残余 `dmu`：

```text
min || W^(1/2) (r - J_mu_perp dmu) ||^2
```

Step05 验证了这个思路：CH6 的原始 `dmu` 列约 `79.5%` 可投影到
`[A, phi, d0, dg_s]` 子空间；正交化后，`dmu` 仍可降低残差，但不会一阶改变
主振动幅值。

## 6. 间隙方向和倾斜方向怎么解耦

只看静态部分时，一阶模型可以写成：

```text
J_g  = dV/dg
J_mu = dV/dg * (x - tau)
```

如果 `tau` 不是当前加权观测点集的中心，`J_mu` 里会包含一部分 `J_g`。这时
优化器调 `dmu`，等价于同时调了一部分 `dg`。解耦方法是先做加权 Gram-Schmidt：

```text
alpha = <J_g, J_mu>_W / <J_g, J_g>_W
J_mu_gap_perp = J_mu - alpha * J_g
```

如果 `J_mu = J_g * (x - tau)` 严格成立，这等价于把倾斜支点移动到：

```text
tau_star = tau + alpha
```

也就是选择满足下面条件的支点：

```text
<J_g, J_g * (x - tau_star)>_W = 0
```

写成显式形式就是：

```text
tau_star = sum(W * (dV/dg)^2 * x) / sum(W * (dV/dg)^2)
```

所以更稳的静态参数化不是：

```text
dg + dmu * (x - tau_old)
```

而是：

```text
dg_star + dmu * (x - tau_star)
```

这里 `dg_star` 表示在正交支点 `tau_star` 处的等效间隙偏移。若需要换回旧支点：

```text
dg_old = dg_star - dmu * (tau_star - tau_old)
```

但辨识时应使用 `dg_star`，因为它和 `dmu` 一阶正交。

这只解决 `dg` 与 `dmu` 的静态互相吞噬。若要保护振动幅值，还要进一步把
`J_mu` 对 `[J_A, J_phi, J_dx/d0, J_g]` 整体正交化：

```text
J_mu_safe = (I - P_[J_A,J_phi,J_dx,J_g]) J_mu
```

因此推荐实现顺序是：

```text
第一层：gap_fixed_tilt 主辨识，估计 [A, phi, dx, dg_star]
第二层：只在残差中拟合 J_mu_safe * dmu，作为诊断
```

更重要的是，第二层不能回改第一层的 `A`。残余倾斜可以解释“剩余波形结构”，但不能
重新定义振动幅值。否则 `dmu` 又会通过同频项把振动吃掉。

实际实现时可以采用两阶段：

```text
Stage 1:
  nonlinear fit q = [A, phi, dx/d0, dg_s], dmu_s = 0
  report A from this stage

Stage 2:
  compute residual r = V_obs - V_pred(q)
  compute J_mu at q
  J_mu_safe = (I - P_q) J_mu
  solve dmu_res = argmin ||W^(1/2)(r - J_mu_safe dmu)||^2
  report dmu_res only as diagnostic
```

如果为了画图需要重构波形，可以用：

```text
V_plot = V_pred(q) + J_mu_safe dmu_res
```

但论文表中的振动幅值仍取 Stage 1 的 `A`。

这里有一个必须承认的辨识边界：

```text
J_mu = P_q J_mu + (I - P_q) J_mu
```

其中 `P_q J_mu` 与 `[A, phi, dx/d0, dg_s]` 不可区分。只靠同一段高速电压波形，
无法判断这一部分到底来自真实 `dmu`，还是来自振动幅值/相位/间隙偏移的改变。能够
安全估计的只有：

```text
(I - P_q) J_mu
```

所以 residual dmu 的安全估计不是完整物理倾斜估计，而是“主模型不能解释的倾斜型残差”。
若真实倾斜的大部分落在 `P_q J_mu` 中，单靠当前高速辨识无法无偏恢复 `A`。这时必须
引入额外信息：

```text
低速 tilt 标定先验；
多窗口/多工况共享同一个 dmu_s；
剔除或降权疑似有 tilt 误差的传感器；
应变片或有限元传递系数给出的独立 A 参考；
更多传感器或不同周向/轴向位置以改变 J_mu 与 J_A 的相关性。
```

因此论文主线不能声称“高速同时自由辨识 A、dg、dmu 后完全解耦”。更稳妥的表述是：

```text
本文用低速标定锁定 tilt 方向，仅在高速中估计 gap offset 与振动；
残余 tilt 通过正交投影作为诊断项，用于识别未建模静态波形误差，
但不参与主振动幅值定义。
```

## 7. 推荐实验/仿真结论写法

可以把方法结论写成：

```text
低速标定给出 muGapPerXMm，是物理先验；
高速辨识只估计 dg_s 和振动参数；
残余 dmu_s 与振动幅值存在结构耦合，因此不进入主辨识；
残余 dmu_s 只作为诊断项，并需经过正交化、正则化和显著性检验。
```

这样理论上自洽，也能解释为什么当前 `gap_tilt` 幅值偏低。
