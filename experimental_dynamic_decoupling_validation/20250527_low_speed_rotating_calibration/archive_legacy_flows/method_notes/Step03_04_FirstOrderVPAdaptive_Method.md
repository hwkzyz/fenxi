# Step03_04 一阶模板 VP top-3 同步方法说明

对应程序：

```text
Step03_04_Run_FirstOrderVPAdaptive_20250527.m
```

这是一个实验加速版本，不是当前主程序。

## 方法定位

`Step03_04` 当前默认设置为严格 VP top-3，并采用同步振动约束：

```text
所有 EO 候选 -> 一阶模板 VP 打分 -> 只保留 top-3 EO
-> 对每个 EO 优化 A,phi,dx_c -> 按 weighted RMSE 选最终结果
```

最终前向模型为：

```text
V_pred = T_low(x - dx_c - A*sin(EO*theta + phi))
```

其中：

- `f = EO * rot_freq_mean`，不作为独立连续变量优化。
- `dx_c` 是空间对齐/常值位移等效项，不再解释为独立振动偏置。
- 当前实验只考虑同步主频；连续频率搜索留给异步振动场景。

## 基本流程

每个窗口内：

1. 生成所有整数 EO 候选。
2. 用一阶模板 VP 给所有 EO 打分。
3. 保留 VP 排名前 `3` 的 EO。
4. 对这 3 个 EO 做完整低速模板波形优化，优化变量为 `A,phi,dx_c`。
5. 按 post-refinement `weighted_voltage_rmse` 选择最终 EO。

## top-K 设置

top-K 数量由环境变量控制：

```text
STEP03_04_TOP_K_EO
```

当前默认是 `3`。

## 可选回退设置

程序仍保留可选回退入口，便于以后做稳健性实验，但默认全部关闭：

```text
STEP03_04_ALL_EO_WARMUP_WINDOWS
STEP03_04_GAP_RATIO_FALLBACK
STEP03_04_LINEAR_GAP_RATIO_FALLBACK
```

当前默认不会触发全 EO。

## 优点

- 比 `Step03_01` 快。
- 候选数量固定为每窗 3 个。
- 对同步振动更物理：频率由 EO 和转频确定，不给错误 EO 额外连续频率自由度。
- `dx_c` 的物理含义与理论中的空间对齐项一致。

## 风险

- 如果 VP 阶段把正确 EO 排除，最终完整波形优化无法恢复。
- 若后续要处理异步振动，需要另建连续频率版本，而不应混入当前同步版本。

当前主程序仍以 `Step03_01` 为准。
