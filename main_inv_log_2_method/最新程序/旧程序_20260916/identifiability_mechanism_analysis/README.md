# 有限时长完整波形联合可辨识性机制分析

本目录是冻结联合辨识程序之上的独立理论—仿真分析区。它不修改 `low_high_gap_fast_bridge`、根目录主入口或 `local_func` 中的冻结求解器。

## 科学问题

1. 有限 passage 内的相位演化提供了多少不重复的动态信息？
2. 这些信息在排除间隙、公共偏移和其他干扰方向后还剩余多少？
3. 为什么连续参数在给定 EO 时可恢复，不同 EO 结构仍可能发生选择翻转？

## 首轮最小闭环

```text
P0 冻结证据基线
  ↓
P1 r2 / exact–frozen 相位信息几何
  ↓
P2 nuisance-conditioned local recoverability
  ↓
P3 constrained model distance / bias / residual margin
```

首轮不做 TOA/peak 基线、标定协方差、第二物理几何、全局最坏情形或新的 Funnel 调参。

## 目录约定

```text
identifiability_mechanism_analysis/
├─ README.md
├─ FROZEN_BASELINE.md
├─ ANALYSIS_PROTOCOL.md
├─ CLAIM_EVIDENCE_LEDGER.md
├─ Run_00_Check_Frozen_Baseline.m
├─ Run_01_Verify_Phase_Diversity_Identity.m
├─ config/
├─ local_func/
└─ output/
```

## 运行顺序

```matlab
Run_00_Check_Frozen_Baseline
Run_01_Verify_Phase_Diversity_Identity
Run_02_ExactFrozen_ConditionalGeometry
Run_03_Window_Aperture_Causal_Scan
Run_04_X01_Evidence_Synthesis
Run_05_Formal20_Predictive_Association
```

`Run_00` 只检查依赖和冻结文件；`Run_01` 只验证相位集中度解析恒等式，不运行联合辨识器。

`Run_02`–`Run_05` 依次完成 exact/frozen 条件几何、窗口孔径扫描、X01 冻结证据综合及正式20工况预测关联。它们只在 Route 31 写入新结果。

## 文件边界

- 冻结主程序只读取，不在本目录中建立修改副本。
- 新仿真结果统一写入 `output/<analysis_id>/`。
- 每个结果必须与预先规定的问题、指标和停止条件对应。
- 首轮分析不改 SG、PCHIP、Funnel V2、Top-24、Top-3 或正式成功阈值。
