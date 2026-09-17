# R3 运行说明

从 `Main_R3.m` 开始。R3 只验证同一响应面族内的留出间隙点，不改变响应面或倾角。主程序只保留案例数和结果显示，留出点生成、观测构造和盲辨识在正式入口中按顺序执行。

```matlab
cd(fileparts(which('Main_R3.m')));
Main_R3
```

输出保存在 `r3_heldout_recovery/output/unified_blind_backend/`。
