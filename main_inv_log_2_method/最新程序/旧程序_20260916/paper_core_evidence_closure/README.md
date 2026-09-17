# Route 32：论文核心证据闭环

本目录只补充决定论文主张强度的三项分析，不修改 Route 30 冻结主程序，不调整 Funnel、模板或正式仿真参数。

## 执行顺序

1. `Run_01_EnergyMatchedPhaseAperture.m`
2. `Run_02_ConditionalCRLBMonteCarlo.m`
3. `Run_03_X01ModelDistanceBiasProjection.m`
4. `Run_04_SummarizePaperValue.m`

## 可否证门槛

- **Gate A**：在相同样本、相同局部波形灵敏度、并匹配条件 Jacobian 总能量后，增大 passage 内相位孔径仍应提高最弱信息方向；否则“相位多样性”降级为特定 frozen 反事实。
- **Gate B**：固定正确 EO、模型正确时，条件 CRLB 与 Monte Carlo 标准差应在高/中 SNR 下保持同量级和一致趋势；否则 FIM 只作局部解释指标。
- **Gate C**：X01 的完整重优化残差裕量翻转必须保留；偏差投影只作局部诊断，不宣称普适阈值。

所有输出写入 `output/`。既有 Route 30/31 结果只读复用。
