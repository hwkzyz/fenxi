# Step03_03 一阶模板 VP top-K 方法说明

对应程序：

```text
Step03_03_Run_FirstOrderVPTopK_20250527.m
```

这是一个实验加速版本，不是当前主程序。

## 方法定位

`Step03_03` 使用一阶模板变量投影快速筛选 EO 候选，然后只对 top-K EO 做完整波形优化。

它比 `Step03_02` 更接近理论中的变量投影思想，因为它直接在电压域使用模板导数 `T_s'(x)`。

## 一阶线性化

低速模板前向模型为：

```text
V = T_s(x - u)
```

在小位移近似下：

```text
T_s(x - u) ~= T_s(x) - T_s'(x) * u
```

因此：

```text
V - T_s(x) ~= -T_s'(x) * u
```

又因为：

```text
u = d0 + b_s * sin(EO * theta) + b_c * cos(EO * theta)
```

所以对每个 EO，可构造线性系统：

```text
V - T_s(x) ~= -T_s'(x) * [d0 + b_s sin(EO theta) + b_c cos(EO theta)]
```

## 变量投影打分

对每个 EO，程序用加权最小二乘直接求：

```text
d0, b_s, b_c
```

然后转换为：

```text
A = sqrt(b_s^2 + b_c^2)
phi = atan2(b_c, b_s)
```

该结果只用于快速打分和初值生成，不作为最终结果。

## top-K 选择

程序按 VP 阶段的评分保留 top-K EO。默认 top-K 由环境变量控制：

```text
STEP03_03_TOP_K_EO
```

## 最终优化

保留下来的 EO 进入完整模板前向波形优化：

```text
V_pred = T_s(x - A sin(EO theta + phi) - d0)
```

最终 EO 仍按 post-refinement 的 `weighted_voltage_rmse` 选择。

## 优点

- 比模板反查方式更接近前向模型。
- VP 阶段速度快。
- 适合做加速候选生成器。

## 风险

- top-K 过小可能提前丢掉正确 EO。
- 一阶近似只在位移不太大、模板导数稳定时可靠。

因此该方法适合作为 `Step03_01` 的加速对照。
