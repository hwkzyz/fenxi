# R5 试验辨识方法与完整程序冻结说明

冻结日期：2026-09-06  
冻结对象：`latest_programs_20260905` 目录  
统一生产入口：`Run_R5_Main_20251222.m`

本文件冻结当前 R5 试验方法、程序边界、输入输出合同和运行顺序。后续任何算法、参数或文件接口修改，都必须建立新的版本目录或新的冻结标识，不得静默覆盖本目录的正式定义。

## 一、冻结的方法定义

### 1. 低速前端

对目标叶片采用目标叶片排除的 PC1 响应面：

```text
target-excluded PC1 waveform family
        -> B2 registration
        -> low-speed anchor localization
        -> target-specific fixed sidecar
```

PC1 坐标只用于低速阶段定位目标叶片响应面。高速阶段固定该定位结果，不把 PC1 坐标作为动态状态，不引入动态 `z_w`。

### 2. 高速正式模型

每个高速窗口独立辨识：

```text
EO, frequency, A, phase, dx_c, dg_1, dg_2, dg_3
```

其中 `dg_s` 是每个活动传感器独立的动态间隙变化；`eta_s=0`；不启用旧的 `delta_mu`、`delta_tau` 或动态 PC1 状态。EO 为离散候选，频率在 EO 对应频率附近进行连续小范围细化。

正式方法名为 `R5_PC1_AnchorSurface_SensorWiseDg_FullWave`（中文：目标排除 PC1 锚定响应面—逐传感器动态间隙全波形联合辨识方法）。代码中的 `gap_only` 仅是内部模型分支名，不是正式方法名称。其零间隙限制 `dg_s=0` 必须逐窗口严格复现 `fixed` 正演模型，作为 Zero-gap equivalence 门禁。

### 3. Foundation 的身份

Foundation 只生成公共 raw-window bundle，并提供历史 no-gap 对比结果。Foundation 拟合得到的 `EO、A、phase、dx、eta、VPred` 不会被 Main06 继承；Main06 重新完成候选筛选和全波形辨识。

```text
Foundation bundle != Foundation identification result
```

## 二、正式程序顺序

统一生产入口：

```matlab
setenv('BLADE_RESONANCE_REGION_ID','1');  % B1 / R01
Run_R5_Main_20251222()
```

或：

```matlab
setenv('BLADE_RESONANCE_REGION_ID','4');  % B5 / R04
Run_R5_Main_20251222()
```

完整 Main 链为：

```text
Main01_Extract_OPR_Blade_Timing_20251222
    -> Main02_Calculate_FullTime_BTT_Displacement_20251222
    -> Main03_Prepare_Strain_Resonance_Evidence_20251222
    -> Main04_Detect_BTT_STE_Resonance_20251222
    -> Main05_Foundation_NoGap_VPTopK_20251222
    -> Main05_R5_PC1_FrontEnd_20251222
    -> Main06_GapAware_FullWave_Identification_20251222
    -> R5_PreMain07_HardGate + R5_Audit_SensorWiseDg
    -> R5_Compare_Proposed_vs_Foundation
    -> Main07_Validate_Identification_With_Strain_20251222  [显式开启]
    -> Main08_Visualize_Strain_BTT_Waveforms_20251222       [可选]
```

`Run_R5_Main_20251222` 统一执行 Foundation、R5 前端、Main06、审计和对齐比较；Main07 默认不运行。正式入口拒绝遗留的 `STEP05_*`、`STEP07J_*` 数值或路径覆盖，只允许区域选择、Main07 开关和绘图显示开关。Main07/Main08 不再按时间戳选择“最新文件”，必须读取 HardGate 已检查的精确结果文件。B1/R01 的应变验证属于目标叶片直接交叉验证；B5/R04 的应变验证属于同工况独立跨模态动态验证，结果中必须明确标注 `strain_reference_blade=1`。应变选频沿用 20241106 的独立链：完整窗口固定物理频带 FFT→全局单频细化→逐窗参考频率 ±3 Hz 拟合；不允许以 BTT 频率作为应变搜索中心。

## 三、完整程序清单

### 顶层入口和检查

```text
Check_Latest_Program_20260905.m
Setup_Paths_20251222.m
Config_20251222.m
Run_R5_Main_20251222.m
Run_Latest_20260905_ASCII.cmd
```

### 主程序

```text
Main01_Extract_OPR_Blade_Timing_20251222.m
Main02_Calculate_FullTime_BTT_Displacement_20251222.m
Main03_Prepare_Strain_Resonance_Evidence_20251222.m
Main04_Detect_BTT_STE_Resonance_20251222.m
Main05_Foundation_NoGap_VPTopK_20251222.m
Main05_R5_PC1_FrontEnd_20251222.m
Main06_GapAware_FullWave_Identification_20251222.m
R5_Compare_Proposed_vs_Foundation.m
Main07_Validate_Identification_With_Strain_20251222.m
Main08_Visualize_Strain_BTT_Waveforms_20251222.m
```

历史名称对应关系：

```text
Foundation Main05  ≈ 历史 Main09 Foundation
R5/Main06          ≈ 历史 Main10 GapAware
Main07             ≈ 历史 Main12 strain validation
```

历史名称仅用于追溯，正式运行不得从旧目录调用。

### R5 前端子程序

```text
analysis/r5_anchor_guided_experimental_20260903/
├─ R5_Build_R4_ExperimentalSurfaceFamily.m
├─ R5_Calibrate_B2_AnchorRegistration.m
├─ R5_Localize_R4_Surface_FromLowSpeed.m
└─ R5_Export_Main10_CompatibleLibrary.m
```

### 结果审计程序

```text
R5_Check_Simulation_Experiment_DynamicContract.m
R5_PreMain07_HardGate.m
R5_Audit_SensorWiseDg.m
R5_Compare_Proposed_vs_Foundation.m
```

`R5_Audit_SensorWiseDg` 按每个窗口、所选 EO 和传感器对应的实际
`[deltaGapLowerMm, deltaGapUpperMm]` 判断边界命中；全局 `±deltaGapLimitMm`
仅是外层安全边界，不再作为正式 boundary-hit 判据。

## 四、区域和窗口合同

正式区域由 `BLADE_RESONANCE_REGION_ID` 唯一决定：

| region | target blade | start time | target laps | expected windows | nominal EO |
|---|---:|---:|---:|---:|---:|
| R01 | B1 | 50.2 s | 20 | 18 | EO14 |
| R04 | B5 | 31.8712 s | 20 | 18 | EO18 |

程序拒绝与所选区域不一致的 `STEP05_ANALYSIS_START_TIME`、`STEP05_TARGET_LAPS`、`STEP05_WINDOW_LAPS` 和 `STEP05_SLIDING_STEP_LAPS`。切换 B1/B5 时只修改 `BLADE_RESONANCE_REGION_ID`。

完整运行必须满足：

```text
Foundation coverage = 18/18
Main06 windows      = 18
```

## 五、正式结果门禁

Main07 前必须检查目标叶片、resonance region、起始时间、target laps、窗口覆盖、EO/frequency、`dg_s` 边界、RMSE、Zero-gap equivalence、sidecar manifest、LOBO 排除目标和 sensor calibration 合同。

`PASS_WITH_FLAGS` 表示来源和数值完整性通过，但存在 EO 分支、频率离散或其他科学质量标记；不能写成“全窗口稳定”。

## 六、冻结边界

当前正式方法不包含动态 PC1 `z_w`、共享动态间隙、`delta_mu_s`、`delta_tau_s`、Foundation 参数继承、应变 EO 强制锁定、跨窗口人工 EO 修正或旧目录静默回退。任何扩展都必须建立新的消融或冻结版本。

## 七、当前结果解释

B1 的早期 R5 PC1-Anchor SensorWiseDynamicGap Full-Wave 结果与当前正式入口在相同输入下逐窗口一致。R04 已冻结为 31.8712 s；在该正式配置下 B5 已完成 18/18 且为 18/18 EO18。该起点变更属于区域窗口合同更新，不是高速辨识器的后验 EO 强制。

本次正式重算的 Main06 结果为：平均频率 632.12 Hz、频率标准差 1.78 Hz、平均 RMSE 35.597 mV、Zero-gap equivalence 最大值 0 mV；R5 硬门禁与 sensor-wise `dg_s` 审计均为 PASS。正式结果文件为 `results/gap_aware/Main_GapAware_VPTopK_FullWave_20251222_B5_S123_r04_gapaware.mat`，对应趋势表为同名 `.csv`。Main07 已按统一时间坐标重新运行。

当前正式对齐比较如下。Foundation 是“不考虑动态间隙变化、直接使用低速波形模板”的既有思路；`gap_only` 不出现在正式方法名称中，仅保留为 MAT 内部字段。

| case | method | EO dominant fraction | mean frequency | frequency SD | mean RMSE |
|---|---|---:|---:|---:|---:|
| B1/R01 | Foundation direct template | 13/18 | 570.98 Hz | 19.18 Hz | 62.07 mV |
| B1/R01 | R5_PC1_AnchorSurface_SensorWiseDg_FullWave | 18/18 EO14 | 582.55 Hz | 0.49 Hz | 33.93 mV |
| B5/R04 | Foundation direct template | 13/18 | 626.68 Hz | 188.04 Hz | 58.48 mV |
| B5/R04 | R5_PC1_AnchorSurface_SensorWiseDg_FullWave | 18/18 EO18 | 632.12 Hz | 1.78 Hz | 35.60 mV |

B1 和 B5 的平均逐窗口 RMSE 降幅分别为 45.29% 和 39.14%。该表只陈述对齐计算结果；间隙参数必要性的论证还必须结合 `dg_s` 局部边界、连续性与 Zero-gap equivalence，不只依据 RMSE。

## 八、2026-09-07 最终生产冻结补充

应变和 BTT 的统一时间坐标固定为：

```text
t_BTT = t_strain_raw + strain_time_offset_total_sec - best_tau_sec
best_tau_sec = 30.045883997 s
```

Main04 的原始应变分段 FFT、Main07 和三维频谱共同使用该式。修正后的约 630 Hz 应变区域为 31.545–33.012 s（BTT 坐标），峰值约 31.911 s，与 R04 正式起点 31.8712 s 一致。R04 是正式试验方案编号；修正后的 Main04 自动检测表内该区域编号为 3，两者不得混用。

为避免 Windows MAX_PATH 导致 MATLAB v7.3/HDF5 读写失败，正式 Foundation 输出目录固定为 `step05_r01` 和 `step05_r04`；AdaptiveSG 模板使用等内容短文件名 `Template_AdaptiveSG.mat`，定位结果使用 `loc.mat`。正式结果不得再引用旧 `Z:` 盘路径。

B1、B5 已从统一入口重新生成。两者均满足 18/18 window coverage、Zero-gap equivalence = 0 mV、实际逐传感器 `dg_s` 边界无命中，并且 HardGate 与 SensorWiseDg audit 均为 PASS。Main07 也已在统一时间坐标下重新运行。`Check_Latest_Program_20260905` 现同时检查正式结果非空、Foundation 来源可读取及两套 HardGate，不再只检查源码文件是否存在。
