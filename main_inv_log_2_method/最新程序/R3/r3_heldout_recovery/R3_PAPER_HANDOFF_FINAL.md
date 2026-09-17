# R3 论文图表最终口径

R3 分析保持冻结。主图使用 `Plot_R3_Fig3_PaperFinal.m`：

1. Fixed `P_vib`；
2. Adaptive `P_vib`；
3. Paired adaptive rescue `P_rescue`；
4. Fixed failure-mode composition（按正式 20° 相位门槛）。

主标题为 **Synchronous-vibration recovery under clearance mismatch**。颜色只表示已测试离散工况的经验比例，不表示对连续 `(Delta g, A)` 区域进行了插值或可靠域外推。Panel (c) 的 `Delta g=0` 为零 rescue，含义是 Fixed 在匹配状态本来已成功，不是 Adaptive 失败。

Supplementary 使用 `Plot_R3_Supplementary_Paper.m`：左图为三种方法的空载 PFP，右图只比较 Adaptive 的 `P_vib` 与 `P_joint`，避免把 Fixed/State-matched 的人为已知状态与 Adaptive 的联合间隙识别混写成同一个 joint endpoint。

`state-matched benchmark` 是诊断控制，不称为 upper bound 或 ceiling；`corrected` 只在 Methods 中说明内部坐标公式 `q_delta=g0+(g_truth-g_ref)`，不作为正式方法名。

最终导出位于 `output/figures_final/`：

- `Fig3_R3_recovery_paper_final.*`
- `FigS_R3_null_and_joint_paper.*`
