# R2 主程序说明

R2 验证固定参考状态造成的运动参数偏置，以及引入状态坐标后的 profiled recovery。

## 正式盲辨识入口

论文 R2 的唯一统一盲辨识入口在上一级目录：

```matlab
Run_R2_BlindSingleFrequency
```

它调用本目录的 `Run_R2_UnifiedBlindBackend`。本目录中的
`Run_R2_ProfiledMechanism` 是机制剖析入口，不是 R2 盲辨识主结果入口。

## 机制审计顺序

1. `Run_R2_ProfiledMechanism.m`：机制结果。
2. `Run_R2_CanonicalCurvatureAudit.m`：曲率与 Schur-complement 审计。
3. `Run_R2_LeaveOneStateOutProfile.m`：留一状态验证。
4. `Run_R2_MultistartAudit.m`：多初值一致性。
5. `Run_R2_NoiseProfile15dB.m`：15 dB 噪声验证。

`Run_R2_IdentityAudit.m` 为目标函数恒等性补充审计。其余补充程序和重复版本位于 `归档_补充审计`。

R2 不负责 R3 的全域成功率，也不负责 R4 的 PC1 响应面迁移；相关结论以各自目录的冻结文档为准。

本目录当前只保留 R2 统一盲后端实现。旧机制审计脚本位于
`旧程序_20260916`，不作为主程序调用。
