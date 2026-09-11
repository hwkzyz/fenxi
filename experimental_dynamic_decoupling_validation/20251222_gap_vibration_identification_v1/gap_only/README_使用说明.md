# 20251222 间隙-振动辨识 V1

这是依据 `20250527_gap_vibration_identification_v1` 的内容框架整理出的
20251222 独立版本。本目录拥有自己的主程序、函数、
标定输入、Foundation 准备数据和结果目录，不调用旧的 20251222 程序文件夹。

正式路线已按 `20251222_V1_FROZEN_20260717` 及 `FROZEN_METHOD_R4R5_20260905.md` 固定。具体不可混用的模型边界、
变量和目标函数见 `FROZEN_METHOD_20260717.md`。后续诊断只能写入
`diagnostics`，不得直接改变 Main09/Main10 的正式定义。

## 正式方法

- 低速模板固定使用第 2–21 圈，三个传感器共用同一圈段。
- 六个 OPR 事件采用实测非等间隔角度，不再假设机械等间隔 60°。
- Foundation 和 GapAware 是两种独立方法，候选 EO 不相互混用。
- 两种方法均全 EO 扫描，VP 只负责加权初值与 Top-3 筛选。
- 最终 EO 由各自 Top-3 的完整正演电压拟合决定，最终波形残差不加权。
- `eta_s=0`；GapAware 主模型为 `gap_only`，每个窗口独立辨识。

## 直接运行

1. 先阅读 `FROZEN_METHOD_R4R5_20260905.md`，确认 R4/R5 统一方法和程序顺序。
2. 在 MATLAB 中进入本目录。
3. 运行 `Check_Program_20251222('all')` 检查程序和随附输入。
4. 设置共振区域，例如 R4：

   `setenv('BLADE_RESONANCE_REGION_ID','4')`

5. 按冻结文档运行 Main01–Main04 的原始数据准备链。
6. 运行 `analysis/r5_anchor_guided_experimental_20260903/R5_Run_R4_FrontEnd_Only_20251222.m`，
   完成 B2 registration、PC1 叶片面定位和 Main10 兼容 sidecar 导出。
7. 使用已有动态映射运行 `Main10_GapAware_FullWave_Identification_20251222.m`。
8. 运行 `Main12_Validate_Identification_With_Strain_20251222.m`。

原来的 Main05--Main08 是冻结 PC1 方法之前的共享面/offset-tilt/path 路线，已移入
`old_program/legacy_pre_pc1_20260905`，不再作为当前方法入口。`Main09` 和 `Main11`
属于固定间隙对照，不是 R4/R5 PC1 主链的必要步骤。

同一目录中的 `LegacyFixed025_20251222.m` 也已移入该旧程序目录；这些文件按当前
冻结方法视为错误入口，仅用于历史追溯，不能用于生成最新结果。

默认 R1 使用 B1/50.2 s，R7 使用 B1/39.2 s，R4 使用 B5/31.9 s。

本文件夹的完整试验程序说明见 `程序说明_试验冻结方法_20260905.md`；仿真程序说明位于 R4 仿真目录，不与本试验目录混用。
所有正式结果只写入本目录的 `results`，不会覆盖旧程序及旧结果。

## 目录职责

- `functions`：本版本私有函数，按准备、Foundation、标定、GapAware 和工具分类。
- `inputs/calibration`：OPR 槽角、响应面及 B1/B5 低速间隙库。
- `inputs/prepared/foundation`：已验证的 Step01–Step04 Foundation 输入。
- `results`：新运行结果；固定间隙与考虑间隙严格分开。
- `analysis`、`diagnostics`：后续消融和诊断程序，不放入正式主流程。
- `archive`：仅存放本版本未来淘汰的代码快照。
