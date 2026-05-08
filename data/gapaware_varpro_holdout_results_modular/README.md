# gapaware_varpro_holdout_results_modular

旧 GapAware/VARPRO 流水线保存结果目录。

常用文件：

- `stage0_config.mat`: 公共配置。
- `stage1_template_library.mat`: 静态间隙模板库。
- `stage2_simulated_data.mat`: 单场景仿真真值和数据。
- `stage3_gap_estimation.mat`: 单场景高速映射数据和初始静态间隙估计。
- `stage4d_jointrefine_results.mat`: 单场景旧方法对比结果。
- `stage6g_compare_methods_across_scenarios.mat`: 多场景/多噪声旧方法对比结果。

`../GapVib_Ortho_Verify` 会读取这些文件作为输入和旧基线，但不会在新方法验证时重新计算旧方法。
