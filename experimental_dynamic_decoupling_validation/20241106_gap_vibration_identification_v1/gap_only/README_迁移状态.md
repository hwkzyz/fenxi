# 20241106 新理论迁移状态

当前目录用于将 `20251222/latest_programs_20260905` 的 R5 方法迁移到 20241106 数据。

## 已复制

- R5 surface family / B2 registration / low-speed localization / sidecar export
- R5 hard gate、sensor-wise dg audit、Foundation 对比程序

## 本试验固定数据合同

- target blade: B4
- analysis sensors: S2/S5/S7
- gap-method sensors: S5/S7（S2 保留在分析链，但不进入 gap-aware 反演）

新增 `analysis/R5_Audit_SensorRoles_20241106`，用于在正式运行前确认 S2 没有被误送入 gap-aware 定位，而 S5/S7 仍属于总体分析传感器。

新增 `analysis/R5_CompareSharedIndependentLatent_20241106`，只比较 S5/S7 的共享与独立 latent；S2 不进入该诊断。由于当前低速模板仍缺少非 B4 锚点，该程序必须在锚点条件满足后使用，不能绕过前置条件。

当前 `cfg.r5.latentMode='sensor_conditioned'`，只对 S5/S7 分别保留 `z_{s,b}`；不会因同一 B4 而强制共享响应坐标。
- dynamic case: `3000_3150`
- analysis start: 86.6 s
- low-speed common laps: 1--30
- windows: 20 target blade passes, 3-pass window, 1-pass step, 18 windows
- frequency search: 300--1000 Hz
- low-speed template: `inputs/calibration/foundation/step04_low_speed_template/Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S2357_20241106.mat`

## 锚点适配与中文路径

原始 gap bank 内部其实已经包含 `GapCalibrationBank.responseSurface`，因此不需要从 CSV 反推波形。新增 `Build_ResponseSurface_Adapter_20241106`，用于将该内部响应面展开为 R5 所需的顶层 `responseSurface` MAT。CSV 仅保留用于旧结果对照，不能作为新理论响应面来源。

该适配器已经在英文临时路径中实际执行成功，生成的 `Step05_Response_Surface_20241106.mat` 已复制到本目录的 `inputs/calibration`。响应面尺寸为 801 个空间采样点 × 6 个叶片 × 16 个间隙等级，满足 R5 surface-family 的输入字段要求。

现已确认并使用 `LowSpeedTemplateBank_20241106_B1toB6_S2357.mat` 中的 B1--B6 低速模板和 `GapCalibrationBank_20241106_B1toB6_S57.mat` 中的已知间隙标定。新增 `R5_Build_LowSpeedTemplateBankAdapter_20241106`，将 S5/S7 的 mV 模板转换为 R5 所需的基线去除 V 模板；B2 的 S5、S7 锚点间隙分别使用标定库中的 1.28603 mm 和 1.41665 mm，不再错误地共用 1.0 mm。

因此当前 20241106 路线是“已有标定库辅助的 target-excluded R5”，不是把 B4 自身当作锚点。锚点来源、传感器独立间隙和适配文件均记录在 `r5_preparation_manifest.mat`。

MATLAB 对中文长路径下的 v7.3 HDF5 MAT 直接读取不稳定。所有 R5 入口先通过 `R5_StageMatForMatlab` 将输入复制到短 ASCII 临时路径再读取，源文件不移动、不覆盖，输出中的路径仍保留原始来源。暂存文件名采用每次唯一的 ASCII 名称，两个工况可以并行运行而不会互相覆盖。

## 已验证的动态结果（20260908）

完整执行 `Run_R5_Preparation_20241106` 和动态识别后，18/18 个窗口通过输入与有限性检查，EO 均为 12。R5 动态 RMSE 中位数为 63.346 mV，Foundation 对照为 43.654 mV，中位差为 +19.539 mV；S5/S7 的独立 `dg_s` 均值约为 0.0689/0.0709 mm。按 1e-5 数值容差检查，18 个窗口中有 16 个振幅接近上界 0.5 mm，说明当前振动幅值参数存在明显饱和，结果不能解释为 R5 已经优于 Foundation。

本试验的动态验证使用归档的 Foundation 结果作为对照，因为当前迁移目录没有新的 Foundation 输出；该来源路径已写入动态结果结构中。下一步应优先检查 B4 的目标模板与响应族一致性，并在解除振幅饱和后重新进行窗口识别。

已完成振幅上界敏感性诊断：将上界从 0.50 mm 放宽到 0.80/1.00 mm 后，振幅中位数约为 0.5675 mm，R5 RMSE 中位数仅由 63.346 mV 降至 62.883 mV，仍比 Foundation 高约 19.23 mV。因此 0.50 mm 上界会暴露饱和，但不是主要误差来源；不能简单放宽上界来宣称方法改进。诊断结果见 `results/r5_amplitude_bound_audit.csv`。
## 动态接口冻结（20260907）

正式动态接口已改为 sensor-conditioned reference-anchored 形式：静态阶段冻结 `z_{s,b}`、`g^L_{s,b}`、registration、gain 和低速模板；动态窗口只优化 EO、A、phase、dx 以及各活动传感器独立的 `dg_s`。不使用 `Δmu_s` 或 `Δtau_s`。可由 `R5_Build_SensorConditionedDynamicModels_20241106` 生成模型；当前 gap-aware 传感器仅为 S5/S7，S2 不进入 gap 增量反演。

20260909 对非锚定目标叶片的低速定位增加了冻结的 `target_x_offset_mm` 静态合同参数（±0.30 mm 扫描），动态阶段仍保留独立 `dg_s`，不增加动态 `Δtau_s`。随后复核动态查询坐标后，20241106 最终采用原低速模板查询坐标（目标偏移只保留为定位诊断量），重新识别中位 RMSE 约 62.16 mV、均值约 60.90 mV，较原约 63.35 mV 有小幅改善；S7 仍需继续复核。
## Main01--Main13 程序归档

旧程序 Main01--Main13 已放入本目录，作为 legacy/reference 复现链；它们可用于原始 BTT、应变、旧响应面和旧 GapAware 结果的生成与对照，但不替代 R5 正式入口。旧程序中的动态倾斜/平移自由度不得带入 R5 正式识别。
