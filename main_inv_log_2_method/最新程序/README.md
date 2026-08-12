# main_inv_log_2_method 说明

这是当前“反比-对数双基函数间隙-振动解耦方法”的主程序目录。

## 当前默认设置

当前默认主方法入口：

```matlab
run_inv_log_2_main_method.m
```

当前默认实际实现：

```matlab
run_low_high_gap_fast_vp_method.m  % Route 30
```

当前默认静态间隙库：

```text
data/间隙的影响/直叶片2mm_统一间隙库0.2_0.1_1.5mm.txt
```

该库包含：

```matlab
0.2:0.1:1.5  % mm
```

旧库路径仍保留在 `cfg.legacyDataFile` 中，只作为对照或历史分析使用。

## 手动跑单个工况

推荐入口：

```matlab
main_single_case_manual/Run_Main_Single_Case.m
```

当前默认单工况使用：

```matlab
g_true_mm = 0.8;
A_true_mm = [0.25, 0.15];
f_true_Hz = [500, 1300];
snr_dB = 10;
```

`0.8 mm` 位于统一库中部，比边界点更适合作为默认演示工况。

## 文件夹作用

- `01_trust_domain_calibration`：基于当前默认静态库重新标定可信空间。
- `02_fixed_trust_domain_pipeline`：固定可信空间下的主流程分析。
- `03_dx_true_injection_pipeline`：显式注入空间偏移 `dx_true` 的仿真分析。
- `04_vibration_amplitude_boundary`：振动幅值边界、失效机制和效率分析。
- `05_gap_recovery_analysis`：统一间隙库下的间隙恢复、留一验证、库密度敏感性和可视化分析。
- `11_simulation_completion`：围绕论文三点创新补齐的主对照仿真，包括联合波形模型混淆对照、变量投影消融和自适应间隙局部搜索范围验证。
- `main_single_case_manual`：只跑一个工况的手动入口。
- `local_func`：主方法依赖的底层函数。
- `trust_domain_result`：当前固定可信空间结果。

## 当前状态

主程序已经切换到统一间隙库 `0.2:0.1:1.5 mm`。如果后续更新 COMSOL 静态库，需要先重新生成统一库，再运行：

```matlab
01_trust_domain_calibration/Step00_Calibrate_TrustDomain.m
```

然后再跑主程序或 `05_gap_recovery_analysis`。
