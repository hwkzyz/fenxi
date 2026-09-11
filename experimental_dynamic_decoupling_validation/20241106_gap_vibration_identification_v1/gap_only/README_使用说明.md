# 20241106 间隙—振动辨识正式程序包

本目录按照实际物理流程排列。根目录中的 `Main01`—`Main13` 都是可以直接打开查看的真实程序，不是调用隐藏大程序的包装入口。

## 正式流程

```text
BTT与应变选区
  Main01  BTT位移和物理叶片编号
  Main02  应变频谱与共振区域

低速与间隙标定
  Main03  静态间隙响应面
  Main04  叶片偏置/倾斜修正响应面
  Main05  全叶片低速模板库
  Main06  CH5/CH7间隙标定库

振动辨识
  Main07  Foundation固定间隙辨识
  Main08  GapAware逐窗口动态映射
  Main09  GapAware全波形辨识
  Main10  两种正式方法比较

应变物理验证
  Main11  应变—BTT全局时间对齐
  Main12  应变频率、幅值与TRAC验证
  Main13  参数趋势及实测—预测波形可视化
```

## 常用运行顺序

第一次从公共原始数据和标定数据完整生成时，依次运行 `Main01`—`Main13`。如果 `inputs` 中的正式标定库已经确认有效，日常辨识只需要运行：

```matlab
Main01_Build_BTT_Displacement_20241106
Main07_Foundation_FixedGap_Identification_20241106
Main08_Build_GapAware_DynamicMap_20241106
Main09_GapAware_FullWave_Identification_20241106
Main10_Compare_FixedGap_GapAware_20241106
```

应变验证另外运行：

```matlab
Main02_Prepare_Strain_Resonance_Evidence_20241106
Main11_Estimate_Strain_BTT_TimeAlignment_20241106
Main12_Validate_Identification_With_Strain_20241106(86.6)
Main13_Visualize_Strain_BTT_Waveforms_20241106(86.6)
```

`Main11` 默认需要配置文件中 `cfg.strain.alignmentCasesSec` 指定的多个时间点结果同时存在。各时间点结果文件名包含时间标签，不再互相覆盖。

## 正式辨识约束

- 20圈输入，3圈滑动窗口，步长1圈，共18个窗口；
- 搜索频率范围 `[300,1000] Hz`；
- VP Top-K默认3，仅用于候选生成；
- 每个Top-K候选独立进行连续频率细化；
- 每个滑动窗口独立辨识，不跨窗口共享频率、相位或EO；
- 不设置参考阶次，不允许固定或强制EO；
- 理论偏置 `eta=0`，不加载非零固定eta先验；
- 最终全波形候选排序使用普通、不加权电压RMSE；
- CH5、CH7允许辨识窗口内间隙修正，CH2保持直接模板通道。

## 文件夹用途

- `functions`：配置适配、原始数据重建和数学辅助函数，不放正式主流程；
- `diagnostics`：OPR、窗口选取等可选诊断程序；
- `inputs`：随程序包保存的标定输入和准备数据；
- `results`：标定、辨识、比较和应变验证结果；
- `docs`：冻结基线和补充说明；
- `tools`：可选批处理工具；
- `archive`：旧包装版和过渡程序，仅用于追溯。

运行前可以执行：

```matlab
Check_Program_20241106('all')
```

## V1 freeze notes (2026-07-18)

- Main05 now rebuilds the formal B4/S2357 template bank only from the bundled Step04 template. It does not read the historical rotating-calibration folder.
- The bundled low-speed template records one common lap range, 1--30, for CH2/CH5/CH7. Check_Program verifies this provenance.
- Foundation Main07 and GapAware Main09 use the same per-window continuous-frequency refinement setting (search range 300--1000 Hz, local half-width +/-2 Hz).
- Foundation candidate ranking uses ordinary unweighted voltage residual RMSE. VP/gradient weights are initialization diagnostics only.
- Foundation results for 86.6 s were regenerated after the change, and Main10 was rerun to refresh the fixed-gap versus GapAware comparison.
- The bundled Step04 source currently contains the validated target blade B4. Main05 therefore refuses to fabricate missing B1--B3/B5--B6 entries.
- Main03 and Main04 now read the bundled static multi-gap waveform library and blade-offset table under `inputs/calibration`; neither step reads a sibling analysis/program folder.
- Main06 calibrates the blade IDs actually present in the formal common-window template bank (currently B4).
- Main09 prefers the package-local rebuilt template/gap banks in `results/calibration`, with the frozen banks in `inputs/calibration/gap` used only as fallback.
- The refreshed 86.6 s comparison verifies that the GapAware `fixed` branch exactly reproduces Foundation in all 18 windows. The `gap_only` branch identifies EO12 in 18/18 windows, with mean frequency 631.93 Hz, mean amplitude 0.36523 mm and mean voltage RMSE 38.293 mV.
- Step04 was fully rerun on 2026-07-18. The formal B4 CH2/CH3/CH5/CH7 templates are all rebuilt from laps 1--30, all record `common_window_from_package_config`, and all pass with `quality_status=good`.
- Target-specific runtime banks now use truthful names: `LowSpeedTemplateBank_20241106_B4_S2357.mat` and `GapCalibrationBank_20241106_B4_S57.mat`. The older `B1toB6` filenames remain only for frozen bundled fallback compatibility.

该检查会验证13个正式程序、标定输入、准备数据、`eta=0`、普通RMSE目标，以及根目录中不存在旧 `Run_*` 包装入口和隐藏的 `src` 主程序目录。
