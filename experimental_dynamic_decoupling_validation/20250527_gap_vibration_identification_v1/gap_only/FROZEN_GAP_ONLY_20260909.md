# gap_only 冻结说明（2026-09-09）

本目录冻结为 20250527 历史 V1 gap_only 程序包。程序源文件、输入和已有结果不再修改；旧运行包装器已移入 `old_program/run_wrappers_20260909`。

已通过 `Check_Program_20250527('all')`。

复核记录（2026-09-09）：13 个 Main 文件均通过 MATLAB `checkcode`；`Check_Program_20250527('all',false)` 通过；Main01 实际运行成功，既有固定间隙、GapAware 对比和应变验证结果保持不变。

补充复核：Main02、Main03、Main04、Main07 已实际运行成功。Main05/Main06 运行时发现历史输入目录曾随程序移动而缺失，已从同一历史 V1 的 `20241106` 输入包恢复到本包 `inputs/calibration`，不改变计算代码或模型。

路径复核：20250527 程序要求的历史公共标定库位于版本目录同级位置，已恢复 `static_gap_waveform_library_20260524_2000Hz` 和 `reference_blade_gap_analysis/results/analysis_04_blade_offset_summary_matched.csv`。随后 Main05、Main06 均已实际运行成功。
