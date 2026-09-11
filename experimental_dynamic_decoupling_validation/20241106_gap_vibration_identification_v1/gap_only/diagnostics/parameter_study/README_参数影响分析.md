# 参数影响分析

本目录只负责参数敏感性分析，不改变正式 `Main01`—`Main13` 的算法。
每个案例使用单独的结果目录，并保存配置、环境参数、运行状态和汇总指标，
因此不同参数结果不会互相覆盖。

## 1. 查看可用参数案例

```matlab
cd diagnostics/parameter_study
Run_ParameterStudy_20241106
```

当前预置了以下单变量分析组：

- `vp_top_k`：Top-K取2、3、5；
- `frequency_refine`：连续频率细化范围取正负0.5、2、5 Hz；
- `gap_limit`：间隙增量边界取正负0.05、0.10、0.25 mm；
- `window_size`：3圈和5圈滑动窗口；
- `vp_weight_power`：VP梯度权重指数取1、2、3。

正式约束 `eta=0` 和最终普通电压RMSE不作为可调因素。

## 2. 预览而不运行

```matlab
Run_ParameterStudy_20241106('vp_top_k', false)
```

## 3. 执行一个参数组

```matlab
Run_ParameterStudy_20241106('vp_top_k', true)
```

也可以只执行一个案例：

```matlab
Run_ParameterStudy_20241106('window_5lap', true)
```

每个案例都会依次运行正式程序07、08、09。运行时间可能较长，
但这种方式保证改变窗口、频率范围等上游参数时不会错误复用旧结果。

## 4. 比较结果

```matlab
Compare_ParameterStudy_20241106([], 'vp_top_k')
```

汇总指标包括：

- EO12窗口数和主导EO；
- 平均频率与平均幅值；
- 普通电压RMSE；
- 平均间隙修正量；
- EO候选裕量；
- 连续频率是否触及细化边界；
- 相对正式默认案例的幅值、频率和RMSE变化率；
- 每个案例计算时间和失败信息。

## 5. 增加新案例

只编辑 `Config_ParameterStudy_20241106.m`。对于根配置已有的参数，
写入案例的 `override`；对于Main09保留的诊断环境参数，写入
`environment`。不要直接在正式主程序中建立消融开关。

结果保存在：

```text
results/study/Cxx/
```

短编号目录用于避免Windows长路径限制；案例名称记录在清单文件中。
每个目录中的 `StudyOverride_20241106.mat` 和
`StudyManifest_20241106.mat` 用于记录完整设置和计算状态。
