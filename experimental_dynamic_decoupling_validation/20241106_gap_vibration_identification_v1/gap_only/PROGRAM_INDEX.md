# 正式程序索引

| 顺序 | 程序 | 主要输入 | 主要输出 | 日常辨识必需 |
|---|---|---|---|---|
| 01 | `Main01_Build_BTT_Displacement_20241106.m` | 原始BTT、低速参考、物理叶片标签 | `Step02_BTT_Displacement_20241106.mat` | 是 |
| 02 | `Main02_Prepare_Strain_Resonance_Evidence_20241106.m` | AI1应变、OPR | 共振区域、应变频谱证据 | 应变验证需要 |
| 03 | `Main03_Build_Static_Gap_Response_Surface_20241106.m` | 多间隙静态波形库 | 静态间隙响应面 | 重新标定时 |
| 04 | `Main04_Correct_OffsetTilt_Response_Surface_20241106.m` | 静态响应面、叶片偏置 | 偏置/倾斜修正响应面 | 重新标定时 |
| 05 | `Main05_Build_LowSpeed_TemplateBank_20241106.m` | 低速单叶片模板 | 全叶片模板库 | 重新标定时 |
| 06 | `Main06_Calibrate_GapLibrary_20241106.m` | 修正响应面、模板库 | CH5/CH7间隙标定库 | 重新标定时 |
| 07 | `Main07_Foundation_FixedGap_Identification_20241106.m` | 公共高速波形、低速模板 | 固定间隙逐窗口结果 | 是 |
| 08 | `Main08_Build_GapAware_DynamicMap_20241106.m` | BTT位移、模板库、间隙库 | 当前时间点DynamicMap | 是 |
| 09 | `Main09_GapAware_FullWave_Identification_20241106.m` | Foundation结果、DynamicMap | GapAware逐窗口结果 | 是 |
| 10 | `Main10_Compare_FixedGap_GapAware_20241106.m` | 两种正式结果 | 逐窗口比较表和图 | 是 |
| 11 | `Main11_Estimate_Strain_BTT_TimeAlignment_20241106.m` | 多时间点GapAware结果、应变证据 | 全局时间对齐 | 应变验证需要 |
| 12 | `Main12_Validate_Identification_With_Strain_20241106.m` | GapAware结果、对齐应变 | 频率、幅值、TRAC验证 | 应变验证需要 |
| 13 | `Main13_Visualize_Strain_BTT_Waveforms_20241106.m` | 辨识和应变验证结果 | 最终趋势及波形图 | 论文作图需要 |

根目录只有正式配置、检查程序、13个真实主程序和两份说明文件。可选诊断、辅助函数和旧版本分别位于 `diagnostics`、`functions` 和 `archive`。
