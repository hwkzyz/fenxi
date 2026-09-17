# 冻结基线

## 程序链

- 低速标定：`../calibrate_inv_log_2_low_speed.m`
- 结构化辨识路由：`../run_inv_log_2_structured_main.m`
- 高速低速模板增量生成：`../local_func/simulate_highspeed_from_low_increment.m`
- 低速 adaptive SG：`../local_func/build_low_speed_templates_adaptive.m`
- 双同步 Funnel V2 求解器：`../local_func/run_dual_sync_voltage_vp_funnel.m`
- 正式 support-aware 冻结验证：`../low_high_gap_fast_bridge/Run_DualSyncSupportAwarePostFreezeValidation.m`

### 可复现性注记

`最新程序/local_func/make_fixed_trust_template_library.m` 仍调用未包含在“最新程序”根目录中的 `make_response_template_library.m`。Route 31 已将该原样构库函数放入自身 `local_func` 作为显式分析依赖。这不改变算法，但表明正式发布前还需对“最新程序”做一次自包含性整理。

## 冻结方法

- 完整 DAQ 时间记录加噪，然后进行 OPR 空间映射和可信窗口截取。
- 低速模板使用 support-aware adaptive SG，窗宽由整转分组交叉验证选择，空间查询使用 PCHIP。
- 高速前向模型保留低速实测模板，静态响应面仅提供高低速间隙增量。
- 四类结构：单同步、单异步、双同步、同步—异步。
- 双同步 Funnel V2：所有 EO 对先获得一次 joint-GN，Top-24 再获得一次 joint-GN，Top-3 进行完整非线性精修。
- Funnel V2 解决候选召回，不保证近退化候选的最终正确选择。

## 已有正式证据

| DAQ SNR | EO 正确 | 全参数成功 | Recall@3 |
|---:|---:|---:|---:|
| 25 dB | 198/200 | 198/200 | 200/200 |
| 15 dB | 192/200 | 190/200 | 200/200 |
| 5 dB | 163/200 | 148/200 | 199/200 |

X01 支持“近共线 EO 子空间与确定性低速模板偏差共同导致零 DAQ 噪声下的模型选择翻转”。X02 支持“弱分量在有限 SNR 下的实用模型选择与连续参数恢复退化”。

## 尚未得到的结论

- passage 内相位演化在 nuisance conditioning 后必然增加独立信息；
- exact waveform 必然增大所有 EO 候选间的模型距离；
- 局部条件信息可以预测盲测案例的成功/失败；
- 当前结论可以跨叶尖几何或传感器响应泛化。
