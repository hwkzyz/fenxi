# 旧 V1 归档审计（20250527）

## 结论

`archive_historical/20250527_gap_vibration_identification_v1` 确实包含完整的旧 V1 间隙变化流程，但该目录不是“未经修改的最初版本快照”。目录中的 `Main10_GapAware_FullWave_Identification_20250527.m` 已包含 2026 年统一理论说明和多个后续诊断输出。因此，比较时应冻结“程序合同 + 结果文件”，不能把整个目录不加区分地称为单一原始 V1。

## 可确认的旧 V1 主流程

1. `Main05_Build_Static_Gap_Response_Surface_20250527.m`：静态响应面。
2. `Main06_Correct_OffsetTilt_Response_Surface_20250527.m`：偏置/倾斜修正。
3. `Main07_Calibrate_LowSpeed_GapLibrary_20250527.m`：低速间隙库。
4. `Main08_Build_GapAware_DynamicMap_20250527.m`：高速动态映射。
5. `Main10_GapAware_FullWave_Identification_20250527.m`：全波形、间隙增量和振动联合辨识。

旧路线保留了每个传感器的低速倾斜/间隙斜率 `muGapPerXMm`，并在动态阶段至少支持 `gap_only`（固定低速斜率，仅识别 `dg_s`）。`gap_tilt` 只应作为敏感性比较，不能直接作为无泄漏正式结果。

## 结果证据

归档结果 `results/gap_aware/Step07J_NestedStaticWarp_VPFullWave_Summary_20250527_B1_S136_oneway_joint_prod_20260831.csv` 给出：EO=14、18 个窗口，平均波形 RMSE=59.03 mV，中位数=58.16 mV。

同目录 `..._formal_recheck.csv` 给出平均 RMSE=59.85 mV；`..._nob_latest_20260830.csv` 为 81.52 mV。由此可见，旧路线性能高度依赖正确的静态修正/传感器坐标合同，不能只比较文件名或固定间隙基础分支。

`analysis/production_calibration_gate_20260831/outputs/ProductionHighSpeedComparison_20250527.csv` 进一步显示，采用完整 transfer 后 18/18 窗口 EO=14，平均 RMSE=59.12 mV；其中 CH3 的平均间隙增量约 0.0282 mm，CH1/CH6 约为 0。这与“各传感器独立间隙增量”的物理设定一致。

## 与当前 R5 比较时的关键判据

- 必须固定同一高速原始数据、窗口、OPR 槽映射和普通电压 RMSE。
- 必须分别报告低速模板误差、静态 registration 误差、迁移重建误差和动态拟合误差。
- 旧程序的 `muGapPerXMm` 是低速标定结果，可以作为冻结的外部校准；不能把高速拟合得到的 `dmu` 或 `eta` 直接写回正式 R5。
- 20250527 的 S6 倾斜效应应作为传感器独立 observation/registration 处理，不能通过共享 latent 或动态 `dTau` 强行消除。

因此，下一步正确工作不是复制旧 `Main10`，而是用上述归档结果冻结一条“旧 V1 参考合同”，逐项替换 R5 的低速模板、registration、静态斜率和动态增量接口，定位 RMSE 增量来自哪一层后再改动对应程序。
