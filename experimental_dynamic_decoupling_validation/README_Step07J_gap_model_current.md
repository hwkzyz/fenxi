# Step07J 当前 gap-aware 模型路线说明

更新日期：2026-06-09

## 1. 当前主程序

当前 Step07J 主程序为：

- `20250527_low_speed_gap_prior_decoupling/Step07J_NestedStaticWarp_VPFullWave_20250527.m`
- `20251222_low_speed_gap_prior_decoupling/Step07J_NestedStaticWarp_VPFullWave_20251222.m`

主程序默认只运行主模型 `gap_tilt`，结果后缀为：

```text
_main_gaptilt
```

方法对比不再混在主程序输出里，而是通过 Step09 单独生成：

- `20250527_low_speed_gap_prior_decoupling/Step09_Compare_Step07J_Methods_20250527.m`
- `20251222_low_speed_gap_prior_decoupling/Step09_Compare_Step07J_Methods_20251222.m`

Step09 会生成统一对比结果：

```text
_method_compare
```

其中只比较：

- `direct/fixed`
- `gap_only`
- `gap_tilt`

`gap_tilt_shift` 已从当前对比中移除，不再作为有效方法使用。

## 2. 当前结果文件口径

主结果：

```text
20250527_low_speed_gap_prior_decoupling/outputs/Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt.mat
20251222_low_speed_gap_prior_decoupling/outputs/Step07J_NestedStaticWarp_VPFullWave_20251222_B1_S123_main_gaptilt.mat
```

方法对比结果：

```text
20250527_low_speed_gap_prior_decoupling/outputs/Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_method_compare.mat
20251222_low_speed_gap_prior_decoupling/outputs/Step07J_NestedStaticWarp_VPFullWave_20251222_B1_S123_method_compare.mat
```

旧的 `phase_safe_prevref_auto_all`、无后缀 Step07J、`OPRCenterStd_Dx035_NoEta` 等结果已经归档，不再作为当前程序默认输入。

## 3. 可视化脚本

当前波形重构流程图：

- `20250527_low_speed_gap_prior_decoupling/Step11_Visualize_Step07J_ReconPipeline_20250527.m`
- `20251222_low_speed_gap_prior_decoupling/Step11_Visualize_Step07J_ReconPipeline_20251222.m`

这两个脚本优先读取 `_main_gaptilt`，若不存在则回退读取 `_method_compare`。若二者都不存在，会提示先运行 Step00 或 Step09。

方法分图/叠加图：

- `20250527_low_speed_gap_prior_decoupling/Step10_Visualize_Step07J_MethodSeparatePipeline_20250527.m`
- `20251222_low_speed_gap_prior_decoupling/Step10_Visualize_Step07J_MethodSeparatePipeline_20251222.m`
- `20250527_low_speed_gap_prior_decoupling/Step12_Visualize_Step07J_MethodWaveformOverlay_20250527.m`

这些脚本需要 `_method_compare`，用于比较 `direct/fixed`、`gap_only` 和 `gap_tilt`。

20251222 的应变片/TRAC 对比：

- `20251222_low_speed_gap_prior_decoupling/Step12_Visualize_Step07J_TRAC_StrainSpectrum_20251222.m`
- `20251222_low_speed_gap_prior_decoupling/Step13_Visualize_Step07J_vs_StrainFFT_20251222.m`

Step12 优先读取 `_main_gaptilt`，Step13 读取 `_method_compare`。

## 4. 推荐运行顺序

每个日期文件夹中，主流程为：

```text
Step00_Run_Main_LowSpeedGapPrior_*.m
```

如果只需要当前主模型，运行 Step00 即可。若需要方法对比和多方法可视化，再运行：

```text
Step09_Compare_Step07J_Methods_*.m
Step10_Visualize_Step07J_MethodSeparatePipeline_*.m
Step11_Visualize_Step07J_ReconPipeline_*.m
```

20251222 若需要应变片验证，再运行：

```text
Step12_Visualize_Step07J_TRAC_StrainSpectrum_20251222.m
Step13_Visualize_Step07J_vs_StrainFFT_20251222.m
```

## 5. 方法判断

当前统一主模型为 `gap_tilt`：

- 20250527 中 `gap_only` 容易收敛到低 EO，`gap_tilt` 对倾斜传感器引起的等效间隙随横坐标变化更稳健；
- 20251222 中 `gap_only` 和 `gap_tilt` 都能稳定 EO14，但使用 `gap_tilt` 作为统一主模型不会破坏稳定性；
- `gap_tilt_shift` 物理解释过柔，容易吸收相位/位置误差，不再进入当前结果、图表或论文结论。

旧结果仍可在归档目录中用于方法演化复核，但当前程序和说明均以 `_main_gaptilt` 与 `_method_compare` 为准。
