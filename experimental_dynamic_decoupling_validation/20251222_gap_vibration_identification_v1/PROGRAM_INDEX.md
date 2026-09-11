# 20251222 V1 程序索引（冻结版）

以下才是当前试验的有效程序链。旧 Main05--Main09、Main11 和
`LegacyFixed025` 按当前方法判定为错误入口，仅保留在
`old_program/legacy_pre_pc1_20260905` 供追溯，禁止作为正式结果来源。

| 顺序 | 有效程序 | 作用 | 主要输出 |
|---|---|---|---|
| 01 | `Main01_Extract_OPR_Blade_Timing_20251222.m` | 提取 OPR、转速及叶片到达时刻 | `results/prepared` |
| 02 | `Main02_Calculate_FullTime_BTT_Displacement_20251222.m` | 计算全时段 BTT 位移并同步应变 | `results/prepared` |
| 03 | `Main03_Prepare_Strain_Resonance_Evidence_20251222.m` | 生成应变侧共振证据 | `results/prepared` |
| 04 | `Main04_Detect_BTT_STE_Resonance_20251222.m` | 检测 BTT-STE 共振区域 | `results/prepared` |
| 05 | `analysis/r5_anchor_guided_experimental_20260903/R5_Run_R4_FrontEnd_Only_20251222.m` | 执行 PC1 波形面、B2 registration、低速 anchor 定位及 Main10 sidecar 导出 | `results/r5_frontend` |
| 06 | `Main10_GapAware_FullWave_Identification_20251222.m` | 使用 sidecar 的原始 Main10 动态辨识 | `results/gap_aware` |
| 07 | `Main12_Validate_Identification_With_Strain_20251222.m` | 独立应变频率/幅值验证 | `results/strain_validation` |
| 08 | `Main13_Visualize_Strain_BTT_Waveforms_20251222.m` | 可选的最终图形输出 | `results/figures` |

第 05 步入口内部固定调用：
`R5_Build_R4_ExperimentalSurfaceFamily`、`R5_Calibrate_B2_AnchorRegistration`、
`R5_Localize_R4_Surface_FromLowSpeed` 和 `R5_Export_Main10_CompatibleLibrary`。
`R45_*` 仅是统一验证、消融、敏感性和报告程序，不属于生产主链。

共振区域通过环境变量 `BLADE_RESONANCE_REGION_ID` 选择：R1=`1`、R4=`4`、
R7=`7`。程序据此自动选择时间、叶片和文件标签。
