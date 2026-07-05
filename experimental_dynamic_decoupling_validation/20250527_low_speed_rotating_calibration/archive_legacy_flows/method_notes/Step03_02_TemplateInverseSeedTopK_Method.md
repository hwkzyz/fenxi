# Step03_02 模板反查种子 top-K 方法说明

对应程序：

```text
Step03_02_Run_TemplateInverseSeedTopK_20250527.m
```

这是一个实验加速版本，不是当前主程序。

## 方法定位

`Step03_02` 的目标是在保持最终波形判据不变的前提下，减少进入完整波形优化的 EO 数量。

它的思路是：

```text
先用模板电压反查得到位移观测 -> 快速给 EO 排名 -> 只保留 top-K EO 做完整波形优化
```

## 输入数据

与 `Step03_01` 相同：

```text
Template_LowSpeedRotating_B1_<sensor_tag>_20250527.mat
DynamicMap_B1_<sensor_tag>_SlidingWindows_20250527.mat
```

## 种子构造

对每个动态采样点，程序先用当前传感器的低速模板反查电压对应的静态坐标：

```text
x_stat = T_s^{-1}(V)
```

然后构造位移观测：

```text
Y_obs = x - x_stat
```

这一步把电压波形暂时转换成一个近似位移序列。

## EO 快速评分

对每个 EO，拟合：

```text
Y_obs ~= A * sin(EO * theta + phi) + d0
```

等价地，用线性最小二乘求：

```text
Y_obs ~= b_s * sin(EO * theta) + b_c * cos(EO * theta) + d0
```

然后转换为：

```text
A = sqrt(b_s^2 + b_c^2)
phi = atan2(b_c, b_s)
```

得到每个 EO 的初值和快速评分。

## top-K 选择

程序按快速评分排序，只保留 top-K EO 候选。默认 top-K 由环境变量控制：

```text
STEP03_02_TOP_K_EO
```

若不设置，则使用脚本默认值。

## 最终优化

被保留的 EO 仍然进入和 `Step03_01` 相同的完整低速模板波形优化：

```text
V_pred = T_s(x - A sin(EO theta + phi) - d0)
```

最终 EO 仍由 post-refinement 的 `weighted_voltage_rmse` 决定。

## 优点

- 比 `Step03_01` 快。
- 最终判据仍是完整波形 RMSE。

## 风险

- 模板电压反查可能有多值或局部歧义。
- 若快速评分阶段把正确 EO 排除，最终优化无法恢复。

因此该方法只作为加速对照，不作为当前主程序。
