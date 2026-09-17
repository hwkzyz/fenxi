# R3 Completion Audit

## 1. Final verdict

**C. R3 EVIDENCE COMPLETE, WITH IDENTIFIED LIMITATION.**

R3 的正式非零振动证据、固定法错误归因机制、adaptive 的 8 条残余失败机制，以及跨间隙空载误报机制均已有可追溯数据和交叉检查。尚未解决的是跨间隙空载误报控制，而不是 R3 主科学问题的证据缺失。

该 verdict 不等于可以直接提交当前图表：主结果文件中的旧 `state_matched` 是 absolute-coordinate 版本，不能再作为论文控制组。修正后的控制组已完成并保存为独立结果，论文汇总/绘图必须切换至该结果后才可进入制图。

## 2. Protocol consistency

| 检查项 | 结论 | 证据 |
|---|---|---|
| 校准/留出状态 | PASS | `R3_Protocol.m`: calibration `[0.2 0.3 0.4 0.6 0.8 1.1 1.3 1.4 1.5]`，held-out `[0.5 0.7 0.9 1.0 1.2]`，参考 0.8 mm。 |
| 冻结 C0 与 gRef | PASS | `Run_02_FreezeReferenceCalibration.m` 冻结单一 C0；明细 `g_low_model_mm=0.8111225572`。 |
| EO 真值/未知反演 | PASS | `R3_Protocol.m` 定义 `eoTrue=10`；`r3_run_three_methods_global_final.m` 未将 EO=10 固定给 adaptive。 |
| 幅值、相位、重复数 | PASS | 7 个幅值（含 A=0），75 个预定义相位-高种子配对；所有正式 cell 均为 n=75。 |
| 噪声冻结 | PASS | `Run_02_FreezeReferenceCalibration.m` 生成 reference-equivalent `sigmaV`，`Run_05_Main15dB_Worker.m` 对每一正式记录复用它。 |
| 同一 noisy record 比较三方法 | PASS | `r3_fit_rows.m` 对一个 `r3_generate_high_case` 的 `highMap` 连续运行三种方法。 |
| 非 R3 路线混入 | PASS | R3 正式运行/汇总代码没有 Step26 或 synthetic-spectrum 调用；仅 README 提到旧研究。 |

## 3. Data integrity

- 正式明细 `output/05_main_15db/main_detail_merged.csv` 有 **9450** 行，恰等于 `6 gaps x 7 amplitudes x 75 repeats x 3 methods`。
- 以 `(method, g_truth_mm, A_true_mm, high_seed)` 计的唯一键也是 **9450**；126 个 method x gap x amplitude cell 均为 n=75；没有非有限的 `g_est_mm`、`A_est_mm` 或 `rmse_V`。
- `main_result.mat` 与 `main_detail_merged.csv` 的 Detail 完全一致；MAT 和 CSV 的 Summary/Paired 数值字段一致（表的字符串读入类型不同）。
- `Merge_05_Main15dB_Workers.m` 自带行数、重复 key、每 cell 三方法 triplet 和 method-set 审计。因此正式矩阵完整性为 **PASS**。

求解器 `exitflag` 中存在 0/2/3，但不构成未完成记录：所有行均有有限拟合量，且主端点直接由复算后的 EO、振幅、相位成功条件定义。它应作为补充数据披露，不应被误称为“所有优化都以 exitflag=1 结束”。

## 4. Comparator correctness

- **Fixed**：`r3_run_three_methods_global_final.m` 将 `stateFixed.gHat=C0.pathCal.g0`，并令 gap search width 为零。
- **Adaptive**：同一文件使用 `C0.templateModel`，在 protocol 的 0.4--1.3 mm 数值区间联合搜索；未读取 held-out raw curve 作为 inverse template。
- **旧 State-matched（deprecated）**：同一文件使用 `stateMatched.gHat=gTruth`。这把物理 absolute gap 错当成内部 reference-anchored query coordinate，因而旧 formal `state_matched` 结果不能作为最终控制组。
- **修正 State-matched（论文应使用）**：`Run_09_StateMatchedOffReferenceDiagnostic.m` 第 21 行采用
  `qDelta=C0.pathCal.g0+(g-P.gReferenceMm)`，即“冻结低速锚点 + 真正运行间隙相对参考态的变化”。`Run_16_CorrectedStateMatchedBenchmark.m` 汇总该 `SM_delta` 结果，并在参考态继承 fixed 的参考等效阈值。

修正控制组使用相同 calibration-only 模板和求解器，不使用 held-out raw truth waveform 进行反演；只给定正确的模型内部状态坐标。因此 comparator correctness 对 **Fixed/Adaptive 为 PASS**，对 **最终 corrected State-matched 为 PASS**，但当前 `main_result.mat` 的旧 State-matched 字段为 **deprecated**。

## 5. Non-zero vibration evidence

正式非零记录每方法 n=2700：

| 方法 | P_vib | 直接证据 |
|---|---:|---|
| Fixed | 16.67% | `main_summary.csv`；五个 held-out gap 的 pooled P_vib 均为 0%。 |
| Adaptive | 99.70% (2692/2700) | 单 cell 为 96.0--100%，各 gap pooled 为 99--100%。 |
| Corrected State-matched | 98.63% (2663/2700) | `output/16_corrected_state_matched_benchmark/corrected_state_matched_*`。 |

`main_paired_summary.csv` 对 fixed/adaptive 使用同一 `high_seed` 配对：五个 held-out gap 中只有 adaptive rescue，没有 adaptive-only penalty；0.8 mm 匹配控制两者一致。此结果足以支持：**对非零同步振动，显式运行间隙自适应显著减轻固定参考间隙失配对 EO、振幅与相位恢复的污染。**

Fixed failure decomposition 已有 `output/10_fixed_failure_attribution`：2700 条非零记录中 450 条无失败模式，2027 条为 `EO + amplitude + phase`，195 条为 `EO + amplitude`，27 条为 `EO + phase`，1 条仅 EO；计数完整闭合为 2700。它支持“失配被吸收进振动参数”的机制性解释。

**审计注记**：`Run_10c_FixedFailureAttribution.m` 以 10 deg 而非 protocol 的 20 deg 作为 `phase_fail` 分类门槛。它不改变主 P_vib（由 `r3_pack_result.m` 的 protocol 20 deg 定义），但论文若报告该分解中的 phase 子类比例，必须以 20 deg 重新生成该附属表。

## 6. Adaptive residual failures

正式 adaptive 非零失败恰为 8 条，均为 EO=16，位于 g=0.9 mm（A=0.20/0.25）和 g=1.0 mm（A=0.20）。`output/12_eo_competition_delta_coordinate/eo_competition_delta_summary.csv` 显示 8/8 为 `joint_search_rank_inversion_confirmed`：在 corrected true-state coordinate 下 EO10 优于 EO16，而自由联合 gap search 后 EO16 优于 EO10。

`output/07_phase_noise_separation_diagnostic/phase_noise_separation_summary.csv` 的三条代表波形显示：固定失败噪声扫相位时 EO10 成功率为 96.0--98.7%，固定异常相位扫噪声时仅为 0--18.7%。因此可判为“残余失败机制已充分诊断”为 phase--finite-sampling--joint-search rank inversion；不能宣称唯一物理机理已被证明。

## 7. Null-state analysis

空载问题已完成诊断，但未解决：

- 正式 adaptive 空载 P_FP 为 **272/450=60.44%**；0.5、0.9、1.2 mm 均为 100%。
- 修正 state-matched 在统一 reference-equivalent fixed threshold 下 P_FP 为 **83.78%**（`output/16_corrected_state_matched_benchmark`），说明问题不只来自 free-gap search。
- `output/13_adaptive_null_distribution_audit` 显示空载幅值、EO 和内部坐标随 gap 系统漂移，并非少数 outlier。
- `output/14_adaptive_null_cost_profile` 对 zero-vibration、vibration 与 true/near-true-state candidate 的残差进行了直接比较。
- `output/17b_state_threshold_transfer_audit` 给出 corrected state-matched null 的 gap-specific q95/q99；相对 fixed threshold，q95 比例为 0.42 (0.8 mm)、2.68 (0.7 mm)、3.10 (0.9 mm)、5.53 (1.0 mm)、9.39 (1.2 mm)、10.68 (0.5 mm)。

故下述表述已有充分证据：**空载判决统计量显著依赖运行间隙，参考状态冻结的单一阈值不能保证跨间隙误报受控。** 机制包含 (A) response-family/state-dependent residual 与 threshold transfer，及 (B) free-gap search 的模型竞争/过拟合；后者不是唯一来源。

## 8. Unresolved issues

### A. 会影响主 claim

无。非零同步振动主 claim 不依赖“P_FP 已受控”，且 Fixed/Adaptive 的正式比较完整、配对且机制闭合。

### B. limitation / supplementary

1. 跨间隙 null detection 尚未解决，论文必须将其作为 limitation，不能宣称完整可靠检测域。
2. adaptive 的 8/2700 EO=16 残余失败必须披露为有限相位/采样下的 rank-inversion 风险。
3. optimizer exitflag 分布应在补充材料如实说明。

### C. 已解释但需防止误用

1. 旧 absolute-coordinate State-matched 已被 corrected benchmark 替代。
2. 旧 `R3_ANALYSIS_SUMMARY.md`、`R3_STATE_MATCHED_COORDINATE_AUDIT*.md` 的部分 state-matched 数字不是最终论文口径；以 `R3_ANALYSIS_CORRECTED_FINAL.md` 和 output 16/17b 为准。
3. `Plot_R3_Fig3.m` 与 `Plot_R3_Supplementary.m` 默认读取 `main_result.mat`，仍会画旧 State-matched；不可直接用于最终论文图。

## 9. What NOT to do next

- 不要为降低 P_FP 而重跑 R3 主矩阵；这属于后续检测器方法，而非本 R3 的完成条件。
- 不要将 75 个配对重复误称为 75 种噪声，或继续扩展为 75 x 75 全因子试验来追求当前主结论。
- 不要再对 8 条 EO16 失败做大规模 Monte Carlo；其 rank inversion 与 phase/noise 敏感性已足以作为残余机制诊断。
- 不要使用旧 formal State-matched 或其单独校准的 0.0312257 mm 阈值支撑横向 null 比较。

## 10. Remaining actions before writing

1. 生成一份只包含 Fixed、Adaptive、**corrected** State-matched 的统一论文汇总表；明确旧 state-matched 已弃用。
2. 更新 Fig. 3 和 supplementary 绘图输入，使其读取 corrected State-matched benchmark，而不是 `main_result.mat` 的旧字段。
3. 以 protocol 的 20 deg 门槛重算 Fixed failure composition 附属表。
4. 在图注/方法中明确 75 为均匀相位--seed 配对设计，并报告 Wilson interval。
5. 在讨论中分开写非零恢复主结论与空载判决 limitation。

从科学证据角度，R3 已经可以停止继续分析并进入论文制图/写作阶段；在制图前必须完成上述数据口径同步，且不能把跨间隙空载误报描述为已经解决。
