# R4/R5 冻结方法与程序顺序

## 1. 冻结的方法

R4 仿真和 R5 试验采用同一套前端方法：

1. 用 gap-law 建立静态波形面
   \[
   V(x,g)=B_0(x)+B_1(x)/g+B_2(x)\log(g/g_{ref}).
   \]
2. 将每个叶片的系数拼成向量 `q`。
3. 对系数向量做 PC1 排序，得到一维系数流形 `q(z)`。
4. 用分段线性插值生成候选波形面。
5. 使用低速 anchor 对 `z` 和有效 gap 做 profile 定位。
6. 输出目标叶片 sidecar 波形库，交给原 Main10 做动态识别。

R4 的 registration 是同域恒等映射；R5 使用冻结的 B2 registration 处理传感器/坐标域差异。除该适配外，R4/R5 不再分叉出两套算法。

## 2. 正式主程序

### R4 仿真

目录：`main_inv_log_2_method/最新程序/33_tilt_calibration_transfer`

执行顺序：

1. `Build_TiltGap_WaveformSurfaceData.m`
2. `Run_01_StaticCalibrationTransferAudit.m`
3. `Run_02_NoiselessDynamicGate.m`
4. `Summarize_TiltCalibrationTransfer.m`

该链用于选择和冻结 PC1 方法，不调用真实试验 Main10。

### R5 试验

目录：`experimental_dynamic_decoupling_validation/20251222_gap_vibration_identification_v1`

正式计算顺序：

1. `Main01_Extract_OPR_Blade_Timing_20251222.m`
2. `Main02_Calculate_FullTime_BTT_Displacement_20251222.m`
3. `Main03_Prepare_Strain_Resonance_Evidence_20251222.m`
4. `Main04_Detect_BTT_STE_Resonance_20251222.m`
5. `analysis/r5_anchor_guided_experimental_20260903/R5_Run_R4_FrontEnd_Only_20251222.m`
6. `Main10_GapAware_FullWave_Identification_20251222.m`
7. `Main12_Validate_Identification_With_Strain_20251222.m`

其中 R5 front-end sidecar 和 Main10 是本方法的核心计算节点；Main12 是独立验证，不参与反向调参。

## 3. R5 验证程序

目录：`analysis/r5_anchor_guided_experimental_20260903`

这些程序用于审计和复现，不改变正式主程序定义：

- `R45_Validate_UnifiedTransfer.m`：统一 R4/R5 静态留出验证；
- `R45_Prepare_ExperimentalComparison.m`：Oracle/CommonB2/R5-LOBO 前端对照；
- `R45_Run_DynamicBatch.ps1`：批量运行 18 窗动态识别；
- `R45_Summarize_Dynamics.m`：输入契约和动态结果汇总；
- `R45_IndependentStrain.m`：独立应变频率对照；
- `R45_Check_RidgeResponseSpread.m`：ridge 响应面敏感性；
- `R45_Plot_Closure.m`、`R45_Write_ClosureReport.m`：图表和最终报告。

## 4. 程序边界

- Main05/Main07/Main10 保持为正式生产链，不在验证程序中复制或改写。
- R4 负责方法选择和可控条件验证；R5 负责真实部署、动态稳定性和独立验证。
- `R5-LOBO` 是泛化压力测试，不是生产运行模式。
- 超出训练面覆盖范围时应标记低覆盖/高不确定性，不能宣称唯一绝对 gap。
- 结果输出写入各自 `results` 目录；历史结果不覆盖。
- 冻结前的 R5 方案和消融脚本副本位于 `analysis/r5_anchor_guided_experimental_20260903/archive_legacy_20260905`，仅用于历史追溯。

## 5. 当前冻结结论

PC1 排序系数流形是当前已比较方法中的主方法。该结论限定于现有数据、评价指标和候选方法，不等同于对所有可能算法的全局最优证明。
