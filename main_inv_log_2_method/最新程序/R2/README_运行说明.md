# R2 运行说明

从 `Main_R2.m` 开始。主程序按“参数配置 → 调用正式入口 → 显示结果”的顺序组织，案例数只需修改文件顶部的 `nCases`。正式入口再按“路径与数据 → 仿真配置 → 逐案例生成观测和盲辨识 → 保存结果”执行；数值优化细节保留在 `local_func`，便于复用和单独测试。

```matlab
cd(fileparts(which('Main_R2.m')));
Main_R2
```

输出保存在 `r2_profiled_recovery/output/unified_blind_backend/`。
