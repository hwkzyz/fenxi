# R3 六 Panel 主图最终方案

主图脚本：`Plot_R3_Fig3_SixPanel_Release.m`。

Panel (a) 展示冻结 calibration states、held-out operating states 和 reference state；(b) 展示同一正式 paired record 的 noisy waveform，并标注 Formal Fixed/Adaptive 估计值；(c,d) 为 Fixed/Adaptive 的离散 `P_vib` 矩阵；(e) 为 `P_rescue=P(C_A=1,C_F=0)`；(f) 为按正式 20° 阈值得到的 Fixed failure composition。

主标题：**Synchronous-vibration recovery under matched and held-out clearance states**。

图注必须说明：颜色仅表示已测试离散工况的经验比例，不表示连续区域插值；`Delta g=0` 的 zero rescue 表示 Fixed 在匹配状态没有失败；代表性记录的 waveform 只用于展示同一 noisy record 下的参数归因差异，不替代全矩阵统计。

最终文件：`output/figures_final/Fig3_R3_six_panel_evidence_release.*`。Supplementary 继续使用 `Plot_R3_Supplementary_Paper.m`，右图仅展示 Adaptive 的 `P_vib` 与 `P_joint`。
