# 最新试验程序（冻结版）

这是当前唯一推荐的试验程序目录。旧 Main05--Main09、Main11 和
`LegacyFixed025` 不复制到这里；它们按当前方法视为错误入口，仍只在上级
`old_program/legacy_pre_pc1_20260905` 中留作历史追溯。

统一生产入口是 `Run_R5_Main_20251222.m`。它按当前区域依次执行 Foundation、R5
PC1/anchor 前端和 `R5_PC1_AnchorSurface_SensorWiseDg_FullWave` 方法，并在结束时自动执行来源门禁、基于实际局部范围的 sensor-wise
`dg` 审计，以及与 Foundation 直接低速模板方法的逐窗口对齐比较；Main07 默认不运行。

代码中的 `gap_only` 只是内部数据结构分支名，不是最新方法名称。B1/R01
可以在硬门禁后进入应变验证。B1/R01 是目标叶片与应变参考一致的直接交叉验证；
B5/R04 使用同一试验工况下的应变片作为独立跨模态动态证据，结果中明确标注
应变参考叶片和验证范围，不把它包装成 B5 叶片的直接应变真值。
应变验证沿用旧版的独立选频链：先在预注册物理频带内用完整 18 窗时段确定
应变全局频率，再以该频率 ±3 Hz 逐窗拟合；BTT 辨识频率只在最终比较阶段使用，
不得作为应变搜索中心。

完整冻结说明、全部程序清单和详细运行顺序见 `R5_METHOD_FREEZE_20260906.md`。

## 顺序

1. `Main01_Extract_OPR_Blade_Timing_20251222.m`
2. `Main02_Calculate_FullTime_BTT_Displacement_20251222.m`
3. `Main03_Prepare_Strain_Resonance_Evidence_20251222.m`
4. `Main04_Detect_BTT_STE_Resonance_20251222.m`
5. `Main05_Foundation_NoGap_VPTopK_20251222.m`
6. `Main05_R5_PC1_FrontEnd_20251222.m`
7. `Main06_GapAware_FullWave_Identification_20251222.m`
8. `R5_PreMain07_HardGate.m`、`R5_Audit_SensorWiseDg.m`、`R5_Compare_Proposed_vs_Foundation.m`
9. `Main07_Validate_Identification_With_Strain_20251222.m`（默认不运行）
10. `Main08_Visualize_Strain_BTT_Waveforms_20251222.m`（可选）

## 程序调用关系

主程序本身不互相 `run`；每个 Main 只读取前一步已经保存的结果，并调用本目录内的
MATLAB 函数。实际依赖关系如下：

```text
Setup_Paths_20251222
  └─ Config_20251222

Main01
  └─ Get_20251222_BTT_Config
      └─ Step2_Extract_JiluBlade_20251222
          └─ functions/preparation 下的 Step1、配置和数据读取子程序

Main02
  └─ Get_20251222_BTT_Config
      └─ Step3_Displacement_Calculation_20251222
          └─ functions/preparation 下的位移、应变同步子程序

Main03
  └─ Get_20251222_BTT_Config
      └─ Step3_Align_StrainSpectrum_With_RPM_20251222
          └─ functions/preparation 下的频谱/转速对齐子程序

Main04
  └─ 读取 Main03 的 Step3_Spectrum_RPM_Alignment 结果
  └─ 读取 Main02 的 Step3_Result 结果
  └─ 执行本文件末尾的 resonance/STFT 本地函数

Main05
  ├─ R5_Build_R4_ExperimentalSurfaceFamily
  ├─ R5_Calibrate_B2_AnchorRegistration
  ├─ R5_Localize_R4_Surface_FromLowSpeed
  └─ R5_Export_Main10_CompatibleLibrary
      └─ analysis/r5_anchor_guided_experimental_20260903 下的 R5 子程序

Main06
  ├─ 读取 Main05 生成的 R5 sidecar + R5_Formal_Manifest
  ├─ 读取 Main01--04 形成的试验/模板结果
  └─ 调用 functions/foundation、calibration、gap_aware、utilities 子程序

Main07
  └─ 读取 Main06 的 gap-aware Result
  └─ 读取 Main04 的 Step04_BTT_STE_Resonance_Regions
  └─ 调用本文件末尾的应变参考和绘图本地函数

Main08（可选）
  └─ 独立读取 Main06/Main04 结果并生成另一组应变-BTT 波形图

## R5 冻结动态合同与审计

正式 Main06 的高速动态参数为 EO、A、phase、dx 和每个活动传感器独立的
`dg_s`。PC1 坐标只在 Main05 低速 anchor 阶段定位目标叶片响应面，不能在高速
窗口中作为动态状态自由变化；正式模式也不启用旧的 delta-mu/delta-tau 分支。

Main06 生成新结果后，可显式调用：

```matlab
R5_Check_Simulation_Experiment_DynamicContract()
R5_Audit_SensorWiseDg('<新生成的 Main06 MAT 结果>')
R5_PreMain07_HardGate('<新生成的 Main06 MAT 结果>')
```

这些程序只读取显式指定的新结果文件，不会自动搜索或引用历史结果。科学质量问题
（例如 `dg_s` 接近边界、窗口间跳变或 RMSE 偏大）作为标记保留；只有 provenance、
窗口完整性、NaN、求解失败和支持域错误等硬错误阻止 Main07。
```

因此推荐执行顺序是：

```text
Main01 → Main02 → Main03 → Main04 → Main05 → Main06 → Main07
                                                     └→ Main08（可选）
```

Main01--Main04 是数据准备链；Foundation Main05 复现原有 no-gap Step05 公共 bundle；
R5 Main05 是冻结的 R5 前端；Main06 是正式动态辨识；
Main07 是应变一致性验证；Main08 只是补充制图，不参与 Main06 的辨识，也不能替代
Main07。若 Main05 未成功生成对应 blade 的 sidecar 和 manifest，Main06 会停止而不会
静默使用旧库。

Foundation Main05 完全沿用原有 no-gap Step05 目标函数，只负责生成 Main06 所需的公共
raw-window bundle。R5 Main05 内部调用冻结的 R5 组件：PC1 波形面、B2 registration、低速 anchor 定位和
Main10 sidecar 导出。Main06 是原 Main10，算法和目标函数不变；正式模式下 R5
sidecar 缺失、目标 blade 不匹配、LOBO/manifest 不匹配都会直接报错。历史旧库不属于
本冻结目录，应在 `old_program` 中单独复现。

## 调用边界

所有可执行 MATLAB 源码都在本目录内。Main01--Main04 只调用本目录
`functions/preparation`、`functions/foundation`、`functions/calibration`、
`functions/gap_aware` 和 `functions/utilities` 下的子程序；Main05 只调用本目录
`analysis/r5_anchor_guided_experimental_20260903` 下的 4 个 R5 子程序；Main06--08
只调用本目录下的函数和本目录结果。没有调用上级目录、`old_program` 或其它版本的
Main 程序。`analysis`、`functions`、`inputs` 都是实体目录，不是 junction/link。

这里的“独立”指程序和标定/模板输入已随包复制；原始试验数据仍按 `Config_20251222.m`
中的路径从试验数据盘读取，不把原始数据重复复制进程序包。

## 不移动程序时的 MATLAB 7.3 路径兼容处理

如果当前工程必须保留在中文目录，不要直接让 `matlab -batch` 读取该目录下的
MATLAB 7.3 文件。先在 Windows 命令行建立临时 ASCII 盘符：

也可以直接双击本目录中的 `Run_Latest_20260905_ASCII.cmd`。该入口会自动建立并
清理 `Z:` 映射，执行 Check 和统一的 `Run_R5_Main_20251222`。需要运行 Main07
时，必须先人工检查门禁结果，再设置 `R5_RUN_MAIN07=1` 后重新调用统一入口。

```bat
subst Z: "E:\0小论文+程序\0博士期间小论文+程序\间隙和振动解耦"
```

然后在 MATLAB 中始终通过 `Z:` 进入本程序包并运行主程序：

```matlab
cd('Z:\experimental_dynamic_decoupling_validation\20251222_gap_vibration_identification_v1\latest_programs_20260905')
Check_Latest_Program_20260905()
Main05_R5_PC1_FrontEnd_20251222
Main06_GapAware_FullWave_Identification_20251222
```

程序本身不需要移动；`Z:` 只是 Windows 对原目录的 ASCII 别名。正式重跑时应始终使用
同一个盘符生成和读取 sidecar/Manifest，避免把 `Z:` 生成的 Manifest 与中文路径结果混用。
实验数据目录若同样出现 HDF5 读取错误，也用同样方式映射，例如：

```bat
subst Y: "E:\试验数据"
```

再将 `Config_20251222.m` 中的原始数据根目录改成 `Y:\20251222\...`。任务结束后可解除映射：

```bat
subst Z: /d
subst Y: /d
```

## 完整性

本目录包含入口程序、配置、检查程序、全部 MATLAB 子程序、R5 helper 和完整输入数据。
`functions`、`inputs`、`analysis` 均为本目录内的实体目录，不是目录链接；删除上级旧目录
不会影响本目录的程序依赖。

运行前：

```matlab
cd('<本目录绝对路径>')
Check_Latest_Program_20260905()
setenv('BLADE_RESONANCE_REGION_ID','4')  % R04/B5, start=31.8712 s
```

正式区域运行时，不要再设置 `STEP05_ANALYSIS_START_TIME`、
`STEP05_TARGET_LAPS`、`STEP05_WINDOW_LAPS` 或
`STEP05_SLIDING_STEP_LAPS` 去覆盖区域配置；程序现在会拒绝与所选区域不一致的覆盖值，
以防止把 B1 的窗口误接到 B5 或把不完整 bundle 当作正式结果。需要切换 B1/B5 时，
只修改 `BLADE_RESONANCE_REGION_ID`。

Main05 的结果位于本目录 `analysis/r5_anchor_guided_experimental_20260903/results/r4_frontend_only`，
其中包含 `R5_Formal_Manifest.mat`；Main06/07 的正式结果写入本目录 `results`。如需
强制指定库，可设置 `STEP07J_CORRECTED_LIB_FILE`，但正式结果必须使用 manifest 中的
R5 sidecar。
