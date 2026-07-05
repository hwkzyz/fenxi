# 当前 Step03 使用说明

当前有效范围：只运行 `S13` 和 `S136`。

## 当前主程序

`Step03_01_Run_DirectLowSpeedWaveform_20250527.m`

这是当前文件夹中的主程序。它采用直接低速模板波形优化，不使用 VP、top-K 候选筛选、EO 锁定或旧超高斯频率先验。

主判据为：

```text
所有整数 EO 候选 -> 每个 EO 都做完整波形优化 -> 选择加权波形 RMSE 最小的 EO
```

## 批处理入口

运行当前主程序：

```matlab
run('Step03_Run_Current_S13_S136_20250527.m')
```

该批处理会依次运行：

```text
S13  = sensors [1, 3]
S136 = sensors [1, 3, 6]
```

## 手动运行方式

运行 `S13`：

```matlab
setenv('STEP03_ANALYSIS_SENSORS','1 3');
run('Step03_01_Run_DirectLowSpeedWaveform_20250527.m')
```

运行 `S136`：

```matlab
setenv('STEP03_ANALYSIS_SENSORS','1 3 6');
run('Step03_01_Run_DirectLowSpeedWaveform_20250527.m')
```

## 程序版本顺序

1. `Step03_01_Run_DirectLowSpeedWaveform_20250527.m`
   当前主程序。无加速。每个 EO 候选都完整优化，最终按优化后的加权波形 RMSE 选 EO。

2. `Step03_02_Run_TemplateInverseSeedTopK_20250527.m`
   实验加速版本。先用模板电压反查构造近似位移种子，再保留 top-K EO 进入完整波形优化。

3. `Step03_03_Run_FirstOrderVPTopK_20250527.m`
   实验加速版本。先用一阶模板 VP 给 EO 候选打分，再保留 top-K EO 进入完整波形优化。

4. `Step03_04_Run_FirstOrderVPAdaptive_20250527.m`
   实验加速版本。默认使用一阶模板 VP 严格 top-3，每个窗口只保留 3 个 EO；同步频率由 `EO * rot_freq_mean` 给出，最终阶段优化 `A,phi,dx_c`。该版本不做 EO 邻域扩展，也不默认 warmup 或回退到全 EO。

## 详细方法说明文件

详细说明放在：

```text
method_notes/
```

建议阅读顺序：

1. `method_notes/Step03_00_Method_Structure_Overview.md`
2. `method_notes/Step03_01_DirectLowSpeedWaveform_Method.md`
3. `method_notes/Step03_02_TemplateInverseSeedTopK_Method.md`
4. `method_notes/Step03_03_FirstOrderVPTopK_Method.md`
5. `method_notes/Step03_04_FirstOrderVPAdaptive_Method.md`

## 当前汇总结果

总说明文件：

```text
method_notes/Step03_00_Method_Structure_Overview.md
```

该文件已经汇总四个方法在 `S13` 和 `S136` 上的运行结果、速度和振动参数辨识结果。

结果 CSV 位于：

```text
output/identification/Step03_Method_RunSummary_20250601.csv
output/identification/Step03_Method_ResultSummary_20250601.csv
output/identification/Step03_Method_WindowParameters_20250601.csv
```

当前结论：`Step03_01`、`Step03_02`、`Step03_03`、`Step03_04` 在 `S13` 和 `S136` 的 18 个窗口上都辨识为 `EO14`。`Step03_01` 仍是当前主程序；其他三个文件是加速对照方法。

## 已归档的旧先验路线

以下脚本不是当前主程序：

- `archive_step03_experimental/Step03_Run_Template_Matching_Identification_20250527.m`
- `archive_step03_experimental/Step03_Visualize_Identification_Results_20250527.m`

原因：这些脚本包含早期窗口 `GlobalEO` / `EOLock` 工作流，并保留旧超高斯对照常数，不适合作为当前无先验基线。
