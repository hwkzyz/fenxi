# main_single_case_manual 说明

本文件夹只保留一个最直接的手动运行入口：

```matlab
Run_Main_Single_Case.m
```

## 用途

用于手动修改一个工况，然后直接运行当前默认主方法。

适合：

- 快速检查某个间隙和振动组合；
- 手动改参数后看识别结果；
- 生成单工况波形拟合图；
- 避免批量分析脚本运行时间过长。

## 当前默认工况

```matlab
caseCfg.g_true_mm = 0.8;
caseCfg.A_true_mm = [0.25, 0.15];
caseCfg.f_true_Hz = [500, 1300];
caseCfg.phi_true_rad = [pi/4, -pi/3];
caseCfg.snr_dB = 10;
```

`0.8 mm` 是统一静态库的中部点，适合作为默认演示工况。若要检查边界性能，可以手动改成 `0.3 mm`、`1.4 mm` 或 `1.5 mm`。

## 当前默认静态库

脚本通过 `load_inv_log_2_project_context.m` 读取当前默认静态库：

```text
data/间隙的影响/直叶片2mm_统一间隙库0.2_0.1_1.5mm.txt
```

即：

```matlab
0.2:0.1:1.5  % mm
```

## 输出

默认保存：

```text
results/main_inv_log_2_method/main_single_case_manual/single_case_result.csv
results/main_inv_log_2_method/main_single_case_manual/single_case_waveform_fit.png
```

默认不保存大的 `.mat` 文件。若确实需要保存完整中间变量，可把脚本里的：

```matlab
saveMatResult = false;
```

改为：

```matlab
saveMatResult = true;
```
