# 20250527 低速间隙先验解耦程序说明

本文件夹按文件名排序即可看到当前工作顺序：

```text
Step00 : 主流程总入口
Step01-Step08J : 当前主程序，Step07J 主流程默认 gap_tilt + causal_vp_state dx_ref
Step09 : 额外方法对比，会单独生成 fixed/gap_only/gap_tilt 对比结果
Step10-Step12 : 波形重构可视化
Step90 : 诊断分析程序
```

优先使用下面两个入口：

```matlab
Step00_Run_Main_LowSpeedGapPrior_20250527
Step09_Compare_Step07J_Methods_20250527
```

## 主流程程序

`Step00_Run_Main_LowSpeedGapPrior_20250527.m` 是主流程入口，内部按顺序调用：

```matlab
Step01_Extract_OPR_Blade_Timing_20250527
Step02_Calc_FullTime_BTT_Displacement_20250527
Step03_Locate_Vibration_By_Strain_OPR_20250527
Step04_Detect_Resonance_By_BTT_STE_20250527
Step05_Build_Response_Surface_20250527
Step05I_Learn_OffsetTilt_Shared_Response_Surface_20250527
Step06_Build_HighMap_20250527
Step06I_Calibrate_OffsetTiltShared_GapLibrary_20250527
Step06K_Publish_GapCalibrationLibrary_20250527
Step07J_NestedStaticWarp_VPFullWave_20250527
Step08J_CrossValidate_StaticWarp_20250527
```

`ProjectionFlow_Config_20250527.m` 是当前主流程配置中心，集中给出默认目标叶片、传感器组、起始时间、目标叶片通过次数、窗口长度、滑动步长和 `Delta g_proj` 参数。默认值与当前已验证结果保持一致：

```text
target blade       : B1
analysis sensors   : S136
start time         : 1.5 s
target blade pass  : 20
window/step        : 3 / 1 blade passes
OPR convention     : 6 pulses per rotor revolution
```

对固定目标叶片而言，每转只出现一次目标叶片通过事件，因此 `targetBladePasses=20` 对应连续 20 转的目标叶片波形。程序仍按 `jilublade(:,4)==targetBlade` 且时间晚于起始时间的规则选波形，以保持现有结果不变。

如果只想从中间继续，可以在 `Step00_Run_Main_LowSpeedGapPrior_20250527.m` 顶部修改：

```matlab
runFromStep = 1;
runToStep = 8;
```

## 方法对比程序

`Step09_Compare_Step07J_Methods_20250527.m` 是当前唯一推荐的方法对比入口。主程序的 `Step07J` 默认采用 `gap_tilt`，并以 `causal_vp_state` 产生在线 `dx_ref`；`Step09` 会用 `STEP07J_RUN_MODE=comparison` 单独生成/读取对比辨识结果：

```text
outputs/Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_method_compare.mat
```

对比方法只包含：

```text
direct_low_template_main
fixed
gap_only
gap_tilt
```

`gap_tilt_shift` 不再进入当前对比，因为它额外引入传感器位置平移自由度，训练误差可能更低，但物理约束更弱，不作为论文主方法。

`Step09` 会输出方法汇总表、窗口对比表和趋势图，并可自动调用 `Step10_Visualize_Step07J_MethodSeparatePipeline_20250527.m` 生成分方法波形图。

## 当前主方法

当前主流程按“先标定库、后识别取波形”的方式组织。`Step06I` 在保留原组合库输出的同时，会把每个 blade-sensor 的标定项发布到：

```text
outputs/gap_calibration_library/
```

文件形如：

```text
GapCalib_20250527_B1_S1.mat
GapCalib_20250527_B1_S3.mat
GapCalib_20250527_B1_S6.mat
```

并生成对应索引 CSV。`Step07J` 默认仍读取原来的组合 `CorrectedGapLibrary`，因此数值结果不因这个库化发布层改变；该库化层用于明确论文流程：所有 blade-sensor 标定先完成，识别阶段只根据配置选择目标波形并调用对应标定项。

如果已经存在旧的 `Step06I` 组合标定结果、但还没有库化索引，可单独运行：

```matlab
Step06K_Publish_GapCalibrationLibrary_20250527
```

该脚本只读取已有 `CorrectedGapLibrary` 并发布 `GapCalib_*.mat`，不会重新标定，也不会改变 `Step07J` 的默认输入。

当前统一的论文主方法是带投影搜索尺度的 `gap_tilt`，并默认使用 `causal_vp_state` 产生在线 `dx_ref`。`Step07J` 先在 VP/direct/previous-window 给出的候选 EO 初值处计算

```text
Delta g_proj = epsilon_perp / S_g
Delta g_loc  = alpha * Delta g_proj
```

再用 `Delta g_loc` 作为每个传感器间隙增量 `dg_s` 的局部搜索半宽。`STEP07J_DELTA_G_LIMIT_MM = 0.25 mm` 现在只作为投影尺度的上限和退化兜底；只有显式设置 `STEP07J_USE_DELTA_G_PROJ=0` 时，才回到固定 `0.25 mm` 边界，用于消融对比或调试。

主流程中的物理模型仍是 `gap_tilt`：

```text
fixed    : 不引入高速间隙修正
gap_only : 每个传感器只识别一个间隙增量 dg_s
gap_tilt : 在 dg_s 基础上增加小的 x 相关斜率增量 dmu_s
```

20250527 中有一个传感器存在倾斜安装影响，`gap_only` 只允许整体间隙增量，不能充分描述倾斜造成的波形变化，因此容易把 waveform objective 推向错误 EO。`gap_tilt` 增加小斜率修正后能更稳定地回到 EO14。

## 可视化程序

```matlab
Step10_Visualize_Step07J_MethodSeparatePipeline_20250527
Step11_Visualize_Step07J_ReconPipeline_20250527
Step12_Visualize_Step07J_MethodWaveformOverlay_20250527
```

推荐优先看 `Step10`。它默认读取 `method_compare` 对比结果，为 `fixed`、`gap_only`、`gap_tilt` 分别出图，不把三种方法挤在同一个窗口里。

图中曲线含义：

```text
low-speed template       : 低速标定得到的模板波形
static gap library       : 低速模板结合间隙库后的静态间隙拟合波形
high-speed data          : 扣除识别振动位移后的高速采样点
method final             : 当前方法识别参数下的最终连续重构波形
prediction samples       : 当前方法在高速采样点处的预测值
```

## 诊断程序

`Step90_Analyze_CH1_HighStatic_VibrationAware_20250527.m` 是诊断分析程序，不属于当前主流程。

`archive_previous_versions/` 和 `archive_obsolete/` 中的程序只用于回溯，不作为当前结果入口。
