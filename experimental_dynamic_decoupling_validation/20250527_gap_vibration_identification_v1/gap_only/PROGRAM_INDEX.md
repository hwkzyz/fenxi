# 20250527 V1正式程序索引

| 顺序 | 程序 | 作用 |
|---|---|---|
| 01 | `Main01_Extract_OPR_Blade_Timing_20250527.m` | 提取OPR、转速和叶片通过时序 |
| 02 | `Main02_Calculate_FullTime_BTT_Displacement_20250527.m` | 计算全时段BTT位移 |
| 03 | `Main03_Prepare_Strain_Resonance_Evidence_20250527.m` | 准备应变频谱与共振证据 |
| 04 | `Main04_Detect_BTT_STE_Resonance_20250527.m` | 使用BTT-STE检测共振区 |
| 05 | `Main05_Build_Static_Gap_Response_Surface_20250527.m` | 建立静态间隙响应面 |
| 06 | `Main06_Correct_OffsetTilt_Response_Surface_20250527.m` | 学习叶片偏置和倾斜修正 |
| 07 | `Main07_Calibrate_LowSpeed_GapLibrary_20250527.m` | 标定低速传感器间隙库 |
| 08 | `Main08_Build_GapAware_DynamicMap_20250527.m` | 建立高速共同窗口和动态映射 |
| 09 | `Main09_Foundation_FixedGap_Identification_20250527.m` | 固定间隙Foundation辨识 |
| 10 | `Main10_GapAware_FullWave_Identification_20250527.m` | 考虑间隙变化的全波形辨识 |
| 11 | `Main11_Compare_FixedGap_GapAware_20250527.m` | 对比两种正式辨识结果 |
| 12 | `Main12_Validate_Identification_With_Strain_20250527.m` | 应变频率、幅值和TRAC验证 |
| 13 | `Main13_Visualize_Strain_BTT_Waveforms_20250527.m` | 最终参数和波形可视化 |

根目录中的Main程序均包含真实计算流程。诊断、消融和旧版本分别位于`analysis`和`archive`，不参与正式计算。
