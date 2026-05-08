# 场机理标定展开模型验证结果

## 1. 验证目的

完整理论文档 `THEORY_COMPLETE_FIELD_CALIBRATED_TEMPLATE_CN.md` 给出了主理论模型：

```math
\widehat F(g,x)
=
B_0(x)+B_1(x)/g+B_2(x)\log g
```

该文档只讨论理论模型和辨识机理。本文件只记录与该模型相关的程序验证结果。

## 2. 静态模板库 leave-one-out 结果

对比模型包括：

```text
baseline_1_over_g: [1, 1/g]
inv_log_2: [1, 1/g, log(g)]
power_law_p025_order2: [1, g^{-0.25}, g^{-1.25}]
inv_poly_3: [1, 1/g, 1/g^2, 1/g^3]
cap_edge_1: [1, 1/g, log(g), 1/g^2]
```

主要结果如下：

```text
inv_log_2:
mean template RMSE       ≈ 0.00121
mean derivative RMSE     ≈ 0.00626
boundary template RMSE   ≈ 0.00279
boundary derivative RMSE ≈ 0.01285

power_law_p025_order2:
mean template RMSE       ≈ 0.00143
mean derivative RMSE     ≈ 0.00677
boundary template RMSE   ≈ 0.00346
boundary derivative RMSE ≈ 0.01464

baseline_1_over_g:
mean template RMSE       ≈ 0.00564
mean derivative RMSE     ≈ 0.00637
boundary template RMSE   ≈ 0.01343
boundary derivative RMSE ≈ 0.01266

inv_poly_3:
mean template RMSE       ≈ 0.00759
mean derivative RMSE     ≈ 0.05664
boundary template RMSE   ≈ 0.02507
boundary derivative RMSE ≈ 0.18428
```

结论：

```text
[1, 1/g, log(g)] 在静态模板值、导数和边界外推综合表现上最好；
[1, 1/g, 1/g^2, 1/g^3] 虽有物理来源，但边界外推明显失稳。
```

## 3. 高阶反比项的验证结论

`1/g^2` 和 `1/g^3` 来源于：

```math
\frac{1}{g+s(x)}
=
\frac{1}{g}
-
\frac{s(x)}{g^2}
+
\frac{s^2(x)}{g^3}
\cdots
```

但程序结果表明，当前模板数量有限时，高阶反比项会放大小间隙边界外推误差。

例如 `g=0.2 mm` 边界留一时：

```text
inv_poly_3:
template RMSE   ≈ 0.04973
derivative RMSE ≈ 0.36439

inv_log_2:
template RMSE   ≈ 0.00430
derivative RMSE ≈ 0.02169
```

因此高阶反比项目前只适合作为物理候选项和对照，不适合作为主模板重构模型。

## 4. 高速闭环 smoke 结果

在 `g=0.2 mm, 10 dB` 小间隙边界工况下：

```text
inv_log_2:
gap error       ≈ 0.00256 mm
frequency error ≈ 0.166 Hz
amplitude error ≈ 0.00087 mm

power_law_p025_order2:
gap error       ≈ 0.00262 mm
frequency error ≈ 0.177 Hz
amplitude error ≈ 0.00286 mm

baseline_1_over_g:
gap error       ≈ 0.01571 mm
frequency error ≈ 24.875 Hz
amplitude error ≈ 0.02997 mm
```

结论：

```text
inv_log_2 与 power_law_p025_order2 都能显著改善小间隙边界识别；
inv_log_2 在该工况下略优。
```

## 5. 高速闭环 quick 结果

代表性 gap 工况下：

```text
inv_log_2:
mean gap error       ≈ 0.00818 mm
mean frequency error ≈ 0.137 Hz
mean amplitude error ≈ 0.00215 mm
mean template RMSE   ≈ 0.00208

power_law_p025_order2:
mean gap error       ≈ 0.00680 mm
mean frequency error ≈ 0.147 Hz
mean amplitude error ≈ 0.00332 mm
mean template RMSE   ≈ 0.00262

baseline_1_over_g:
mean gap error       ≈ 0.01557 mm
mean frequency error ≈ 41.645 Hz
```

结论：

```text
power_law_p025_order2 的平均 gap 误差略小；
inv_log_2 的频率误差、幅值误差、波形 RMSE 和模板 RMSE 略优；
结合理论可解释性，inv_log_2 更适合作为主理论模型。
```

## 6. 当前推荐

当前建议模型定位为：

```text
主理论模型: inv_log_2 = [1, 1/g, log(g)]
强数值对照: power_law_p025_order2
传统基线: baseline_1_over_g
高阶反比候选: inv_poly_3，仅作对照，不作为主模型
```
