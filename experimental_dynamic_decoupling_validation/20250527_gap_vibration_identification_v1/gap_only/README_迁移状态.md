# 20250527 新理论迁移状态

当前目录用于将 `20251222/latest_programs_20260905` 的 R5 方法迁移到 20250527 数据。

## 已复制

- R5 surface family / B2 registration / low-speed localization / sidecar export
- R5 hard gate、sensor-wise dg audit、Foundation 对比程序

## 本试验固定数据合同

- target blade: B1
- analysis/gap sensors: S136
- dynamic case: `20250526_2500-3500_t400`
- low-speed case: `20250526_910`（沿用旧程序正式配置）
- analysis start: 1.5 s
- windows: 20 target blade passes, 3-pass window, 1-pass step, 18 windows
- frequency search: 300--1000 Hz
- static response surface: `inputs/calibration/Step05_Response_Surface_20250527.mat`
- gap library: `inputs/calibration/Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136.mat`

## 当前可执行入口

- `Check_R5_Migration_20250527`
- `Run_R5_Preparation_20250527`
- `Run_R5_SensorConditionedDynamic_20250527`
- `R5_Audit_SensorConditionedDynamic_20250527`

## 当前迁移状态与正式边界

动态前端已经完成参数化，窗口数据由本试验的动态包读取，不会误用 20251222 动态结果。由于当前模型是 sensor-conditioned，Main10 的单一共享响应面导出不适用，预处理清单将其明确标记为 `not_applicable_sensor_conditioned`；这不是动态链路失败。

共享 latent 的 12% 离散度只作为诊断信息，不能作为强行压缩三个传感器 `z_{s,b}` 的理由，也不再阻塞独立 latent 的动态识别。

## 已验证的动态结果（20260908）

完整执行 `Run_R5_Preparation_20250527` 和动态识别后，18/18 个窗口通过输入与有限性检查，EO 均为 14，独立 `dg_s` 均未触及边界。R5 动态 RMSE 中位数为 66.754 mV，Foundation 对照为 62.954 mV，中位差为 +4.231 mV；独立 `dg_s` 的均值约为 S1=0.0200 mm、S3=0.0400 mm、S6=0.0183 mm。

上述结果说明新链路已经可复现运行，但当前 R5 尚未优于 Foundation，不能据此宣称性能提升。下一步应在同一 sensor-conditioned 接口下复核 S3 横坐标/基线配准、S6 电压增益及有效模板，并增加重复低速窗口验证。

B2 配准参数也支持这一判断：S3 的 `tau_mm=-0.0671 mm`、`voltage_gain=0.8780` 明显偏离 S1/S6；S6 的 `voltage_offset=-175.67 mV` 且 B2 锚点 RMSE 为 54.42 mV。后续应先做传感器级校准复核，再决定是否进入动态识别。

新增 `analysis/R5_CompareSharedIndependentLatent_20250527`。它将每个传感器独立最优 latent 与同一叶片共享 latent 的 RMSE、gap 和惩罚量逐项输出。该程序是诊断程序，不把 latent 点跨度直接当作失败条件；正式门控前应结合低速重复性噪声、profile 近优区间和 leave-one-sensor-out 稳定性解释结果。

当前 `cfg.r5.latentMode='sensor_conditioned'`：低速定位保留每个传感器自己的 `z_{s,b}`，共享 `z_b` 只作为诊断模式（设为 `shared` 才启用）。不会把不同传感器的响应坐标强行压成一个值。

新增 `analysis/R5_Audit_TransferRegistration_20250527`，用于并排比较 R5 B2 配准和旧 P2 registration 的 `tau/xScale/gain/RMSE`。旧 P2 参数只作为对照，不会自动替换到 R5 response family 中。

审计结果显示：R5 与旧 P2 的 `xScale` 存在系统差异（R5 约 0.996--1.005，旧 P2 约 1.023--1.037），`tau` 和 gain 也有传感器相关差异。旧 P2 的 RMSE 略低，但它是在旧生产 response surface 和旧 transfer 模型上得到的，不能据此判定 R5 配准错误；下一步必须在同一 R5 response family 上做留出区间/重复性验证。
## 动态接口冻结（20260907）

正式动态接口已改为 sensor-conditioned reference-anchored 形式：静态阶段冻结 `z_{s,b}`、`g^L_{s,b}`、registration、gain 和低速模板；动态窗口只优化 EO、A、phase、dx 以及各活动传感器独立的 `dg_s`。不使用 `Δmu_s` 或 `Δtau_s`。可由 `R5_Build_SensorConditionedDynamicModels_20250527` 生成模型；该模型把 `x-dx-u(t)` 先送入冻结 registration，再用 `F_{s,b}(x_F,g^L+dg_s)-F_{s,b}(x_F,g^L)` 计算增量。

## 目标叶片横坐标合同修正（20260909）

低速定位增加了一个仅针对非锚定目标叶片的静态目标波形横坐标偏移 `target_x_offset_mm`，范围为 ±0.30 mm；它不是动态自由度，也不是 `Δtau`。该偏移用于表达目标叶片低速模板与 B2 锚定注册坐标之间的固定合同差异，并在 sidecar 中冻结。动态窗口仍只优化 `EO/A/phase/dx/dg_s`。

重新执行准备与动态识别后，20250527 的 R5 中位 RMSE 由约 66.75 mV 降至约 53.22 mV（18/18 窗口有效）；这说明原先主要误差确实包含目标叶片横坐标合同不一致。该改进不是把旧程序的动态 `Δtau` 带回 R5，而是把低速阶段识别出的固定目标注册偏移冻结在静态定位合同中；动态模型仍只使用冻结注册和独立 `dg_s`。
## Main01--Main13 程序归档

旧程序 Main01--Main13 已放入本目录，作为 legacy/reference 复现链；它们可用于原始 BTT、应变、旧响应面和旧 GapAware 结果的生成与对照，但不替代 R5 正式入口。旧程序中的动态倾斜/平移自由度不得带入 R5 正式识别。
