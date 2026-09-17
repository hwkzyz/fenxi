# 同步振动整数阶次快速搜索（2026-08-09）

## 实现

正式入口 `run_inv_log_2_main_method` 新增显式频率结构路由：

```matlab
cfg.route30FrequencyStructureMode = "dual_sync_sync";
fit = run_inv_log_2_main_method(highMap, templateLib, cfg, lowState);
```

可用模式：

- `general`：通用单/双频辨识；
- `single_sync`：整数 EO 单同步；
- `single_async`：连续频率单异步；
- `dual_sync_sync`：整数 EO 双同步；
- `dual_sync_async`：一个整数 EO 加一个连续异步频率。

同步基优先使用实测轴角 `theta_v`，没有轴角时才使用名义恒转速角度。双同步候选由搜索带内全部整数 EO 自动组合，不需要真实 EO。

若激励阶次来自外部物理信息而不是辨识真值，还可进一步设置：

```matlab
cfg.singleSyncCandidateEO = [14 24];       % 单同步允许的已知 EO 集合
cfg.dualSyncCandidateEOPairs = [14 24];    % 已知双同步 EO 对
cfg.syncAsyncCandidateEO = 14;             % 混合振动中已知同步 EO
```

仿真 benchmark 没有设置这些字段，始终枚举完整 EO 范围，满足“不用真实频率或真实 EO 生成候选”的约束。

## 候选规模

当前转频为 50 Hz，300–1500 Hz 对应 EO `6:30`：

- 整数 EO 数：25；
- 双同步无序 EO 对：300；
- 通用 5 Hz 无序频率对：28,203；
- 单个间隙上的候选对减少：`94.0×`。

快速结构化预算采用 15 个间隙锚点，每个间隙评估全部 300 个 EO 对，共 4,500 次廉价粗筛；每个间隙保留40个，候选库共600个；随后只做400个完整波形迭代重排和12个非线性精修。

## 固定种子结果

固定种子 `20260809`，双同步 `700+1200 Hz`，间隙 `0.8 -> 0.5 mm`，振幅 `0.25+0.15 mm`，测试 `25/15/5 dB`。整数 EO 路径只收到搜索范围和 RPM，没有收到真实 EO。

| 方法 | 正确率 | 平均时间 | 最大频率误差 |
|---|---:|---:|---:|
| 通用快速双频 | 3/3 | 16.891 s | 0.0555 Hz |
| 整数 EO 搜索 | 3/3 | 4.730 s | 0 Hz |

整数 EO 搜索的端到端平均加速为 `3.57×`。其平均阶段耗时为：

- EO-间隙粗筛：1.878 s；
- 完整波形迭代重排：2.069 s；
- Top-K 非线性精修：0.758 s。

## 建议

1. 已知振动为同步结构时，正式程序应显式设置 `route30FrequencyStructureMode`，不要先运行通用双频搜索。
2. 只知道“同步”而不知道具体 EO 时，搜索完整整数 EO 范围；这已经获得约3.6倍端到端加速。
3. 只有当 EO 来自叶片数、齿数、控制命令或已知外部激励时，才设置固定 EO 字段；不得在仿真验证中注入真实 EO。
4. 同步假设不可靠时继续使用 `general`。混合同步-异步使用 `dual_sync_async`，不能把异步分量强制量化到整数 EO。
5. 结构化路径的低预算通过当前三组 SNR 回归，但仍应保留通用路径用于结构假设错误或置信度不足的工况。

结果文件：

- `output/synchronous_order_fast_benchmark/case_metrics.csv`
- `output/synchronous_order_fast_benchmark/summary.csv`

