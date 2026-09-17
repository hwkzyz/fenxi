# low_increment 一致性修正

本次只修正正式仿真前向模型与文档定义不一致的问题，未改动旧程序、旧结果或论文图片。

正式模型固定为

\[
V_s(x)=T_s^{low}(x)+\kappa_s\{R_s(x,g_{0,s}+\Delta g_s)-R_s(x,g_{0,s})\}.
\]

其中低速波形 `T_low` 和标定间隙 `g0` 在低速阶段冻结，高速估计器优化的变量是相对间隙增量 `delta_g`，不是绝对高速间隙。`make_low_increment_path_model` 统一封装该接口；增量查询使用 `eval_path_increment_template`，而 `eval_gap_template` 在固定路径模型上仍保留绝对间隙查询兼容性。

R4 `Main_08`、`Main_10` 的真值生成、高速 VP 筛选和全波精修已统一调用该接口；结果输出仍将 `gHigh = gLow + delta_g` 转换为绝对间隙，便于误差评价。真实间隙只保留在波形生成和最终评价中。

`run_inv_log_2_structured_main` 在 `route30ForwardModel="low_increment"` 时，如果输入不是已完成低速标定的固定路径模型，将直接报错，不再静默回退到绝对响应面。公共 Route-30 求解器内部仍可用绝对查询间隙 `g=g0+delta_g` 访问固定路径评估器；这只是内部坐标表示，实际正演仍是上式的增量模型。`absolute` 仍可用于显式历史/消融配置。

一致性审计脚本为 `low_high_gap_fast_bridge/Run_LowIncrementConsistencyAudit.m`。缓存查询与逐点精确公式之间采用 `1e-4 V` 数值容差；本次审计实测最大差为 `1.64e-6 V`。本次未运行大规模 R2/R3/R4 重算；应先运行该审计和小规模无噪声单频测试，再进行论文数据重跑。
