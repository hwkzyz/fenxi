# R3 论文数据与制图交接说明

## 冻结状态

R3 的科学状态冻结为：**R3 EVIDENCE COMPLETE, WITH IDENTIFIED LIMITATION**。

论文数据的唯一入口为 `output/18_paper_data_sync/r3_paper_result.mat`。其中只包含 `fixed`、`adaptive` 和 `state_matched_corrected`；旧正式矩阵中的 absolute-coordinate `state_matched` 不得用于正文、图表或补充材料。

## 主文建议

Fig. 3 建议保留四个面板：

- (a) Fixed 非零振动恢复域；
- (b) Adaptive 非零振动恢复域；
- (c) corrected State-matched 非零振动恢复域；
- (d) Fixed→Adaptive 的逐记录 rescue probability。

正文主数字：Fixed pooled `P_vib=16.67%`，且五个 held-out gap 均为 0%；Adaptive 为 `2692/2700=99.70%`，单 cell 为 96–100%；corrected State-matched 为 `2663/2700=98.63%`。五个 held-out gap 的 paired comparison 只有 Adaptive rescue，没有相对 Fixed 的 penalty。

建议正文表述：

> 对于非零同步振动，显式运行间隙自适应显著减轻了固定参考间隙失配对 EO、振幅和相位恢复的污染。在五个留出运行间隙中，固定参考法未恢复任何非零振动记录，而自适应法恢复了 2692/2700 条记录；剩余 8 条失败均属于 EO=16 的联合搜索排名反转。

## Fig. 3 图注建议

> **Fig. 3 | Recovery of non-zero synchronous vibration under held-out clearance states.** (a–c) Probability of successful EO, amplitude and phase recovery for the fixed-reference, clearance-adaptive and corrected state-matched methods, respectively. (d) Paired probability that the adaptive method succeeds when the fixed-reference method fails. Each cell contains 75 uniformly spaced phase values, each paired one-to-one with a deterministic high-speed noise seed. The dashed contour denotes the predeclared 0.90 success criterion. The corrected state-matched control uses the internal coordinate `q_delta=g0+(g_truth-g_ref)` and does not use held-out raw response curves in the inverse model.

## Supplementary 建议

Supplementary 主图放两类证据：

- 各 gap 的空载 `P_FP`，明确 corrected State-matched 使用统一 reference-equivalent threshold；
- 各 gap 的平均 `P_joint`，用于说明非零恢复与间隙恢复同时满足容差。

建议图注：

> **Supplementary Fig. Sx | State dependence of null discrimination and joint recovery.** (a) False-positive probability at zero vibration for the fixed-reference, clearance-adaptive and corrected state-matched methods. (b) Mean joint-success probability over the six non-zero vibration amplitudes. Each condition contains 75 uniformly spaced phase–seed pairs. The persistent gap dependence of the null statistic shows that a threshold frozen at the reference state does not control false positives across operating clearances.

其余补充表依次提供：完整 gap×amplitude 矩阵及 Wilson 区间、20° Fixed failure composition、8 条 EO16 全部案例、phase/noise separation、corrected State-matched 完整表、gap-specific null q95/q99，以及正式矩阵完整性审计。

## Discussion 边界

非零振动恢复和空载检测必须分开写：

> The clearance-adaptive formulation substantially reduced the contamination of synchronous-vibration parameters caused by a frozen reference-clearance assumption. This improvement concerns parameter recovery conditional on non-zero vibration. It does not establish clearance-invariant null discrimination: the null amplitude statistic drifted systematically with operating clearance, and a threshold calibrated at the reference state did not control false positives across held-out states. State-conditioned residual correction or threshold calibration is therefore required for a deployable detector.

不能使用“普适扩展可靠域”“所有工况可靠”或“跨间隙误报已解决”等表述。

## 可追溯文件

- 统一结果：`output/18_paper_data_sync/r3_paper_result.mat`
- 统一汇总：`output/18_paper_data_sync/r3_paper_summary.csv`
- 配对汇总：`output/18_paper_data_sync/r3_paper_paired_summary.csv`
- 20° Fixed 分解：`output/18_paper_data_sync/fixed_failure_composition_20deg_*.csv`
- 主图脚本：`Plot_R3_Fig3_Final.m`
- 补充图脚本：`Plot_R3_Supplementary_Final.m`
- 最终导出：`output/figures_final/`
