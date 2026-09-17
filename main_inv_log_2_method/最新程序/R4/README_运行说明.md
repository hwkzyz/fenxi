# R4 运行说明

从 `Main_R4.m` 开始。R4 用多个倾角响应面建立 PC1 流形并进行盲迁移验证。网格数和低速候选数集中在主程序顶部，响应面构造、候选搜索和结果写入由唯一正式入口连续完成。

```matlab
cd(fileparts(which('Main_R4.m')));
Main_R4
```

输出目录由 `Main_10_UnifiedBlindBackendValidation` 的 `opts.outputDir` 控制。
