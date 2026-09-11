# gap_only 冻结说明（2026-09-09）

本目录冻结为 20251222 现有历史 V1 gap_only 程序包。正式 GapAware 主程序已按原有冻结方法运行复现；程序源文件、输入和已有结果不再修改。

正式主链：`Main01–Main04` → `Main06_GapAware_FullWave_Identification_20251222` → `Main07_Validate_Identification_With_Strain_20251222` → `Main08_Visualize_Strain_BTT_Waveforms_20251222`。

说明：本目录同时保留旧版固定间隙对照程序；它们不是正式 GapAware 主链入口。

复核记录（2026-09-09）：`Check_Program_20251222('all',false)` 通过；根目录仅保留正式 `Main05_Foundation_NoGap_VPTopK_20251222.m`，其余两个历史 Main05 位于 `old_program/legacy_pre_pc1_20260905`。Main01 已实际运行并成功复用/刷新 OPR 中心结果，Main06 正式 R5 结果、Main07 应变验证和 Main08 波形可视化均已运行通过。
