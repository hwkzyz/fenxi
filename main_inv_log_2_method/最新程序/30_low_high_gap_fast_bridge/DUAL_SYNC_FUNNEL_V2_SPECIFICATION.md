# 双同步 Funnel V2 正式方案说明

## 1. 适用范围

本方案用于已知振动结构为双同步振动的高速叶尖传感器完整波形辨识。输出两个 EO、两个振幅与相位、低速基准间隙、相对间隙变化和高速拟合残差。

前向模型采用低速模板增量形式：

\[
\widehat V_H(t)=T_L(q(t))+R_{stat}(g_0+\Delta g,q(t))-R_{stat}(g_0,q(t)).
\]

其中 \(q(t)=x(t)-d-u(t)\)，双同步位移由两个 EO 正弦分量组成。

## 2. 固定计算流程

正式求解器不再使用 Top-40 或 Top-80 截断作为主结构：

\[
1496\ \mathrm{VP}
\rightarrow272\ (\mathrm{exact0+fixed\text{-}g\ GN1})
\rightarrow136\times1\ (\mathrm{joint\ GN})
\rightarrow24\times1\ (\mathrm{joint\ GN})
\rightarrow3\ (\mathrm{full\ refinement}).
\]

1. VP 对全部 EO--间隙候选作线性化筛选。
2. 每个 EO 对保留两个间隙状态，进行 exact replay 和一次固定间隙 GN1。
3. 对全部 136 个离散 EO 对各做一次允许 \(g,d,A,\phi\) 联合更新的 joint-GN。此步之后才允许删减 EO 对。
4. 按 joint-GN1 完整电压残差保留 \(K_2=24\) 个候选，再做一次 joint-GN。
5. 最终 3 个候选进入有界完整非线性最小二乘精修。

默认 \(K_2=24\) 是安全裕量。代表性 15 dB 验证中 \(K_2=8\) 已足够，但不应仅据该小范围结果缩小正式预算。

## 3. 噪声与数据处理协议

噪声定义在传感器采集层，而不是在最终可信空间窗口上定义：

\[
s_{DAQ}=\operatorname{RMS}(V^0-\overline{V^0}),\qquad
\sigma=s_{DAQ}10^{-SNR_{DAQ}/20}.
\]

对完整无噪声时间记录加入独立高斯白噪声后，才进行 OPR 空间映射、可信区截取和参数辨识。

- 低速：20 个完整转子周期；
- 高速：8 个完整转子周期；
- 低速、高速记录分别按照各自完整 DAQ 波形的 AC RMS 设置噪声；
- 不对每个随机种子强制缩放为精确 SNR；
- 记录实际 DAQ SNR、噪声标准差及高速振动增量 SNR。

低速模板通过 20 圈空间分箱均值和低速数据内的自适应 SG 交叉验证建立；一套低速标定可复用到同一测量设置下的多个高速记录。

## 4. 正式入口

```matlab
% 默认：15 dB、10 个高速噪声种子、K2=24
[Detail,Summary] = Run_DualSyncFunnelV2Generalization;

% 完整 20-seed 验证
[Detail,Summary] = Run_DualSyncFunnelV2Generalization(15,1:20);

% 指定工况或 K2 的诊断比较
[Detail,Summary] = Run_DualSyncFunnelV2Generalization(...
    15,1:20,{"C06_amp_15_10","C14_amp_18_12"},8);
```

结果写入：

```text
output/dual_sync_funnel_v2_generalization_daq_snr/
```

历史 V1 / Top-80 结果保留在独立目录，仅用于消融比较。

## 5. 当前验证结论

在完整 DAQ 25 dB、20 seeds 下，16 个常规工况和 `X03_corner_equal`、`X04_corner_wide` 均为 20/20 EO 正确及参数成功。`X02_corner_lowamp` 为 19/20。

在完整 DAQ 15 dB、20 seeds 下，所有工况的真 EO 均进入最终 Top-3；弱振幅 `(0.15,0.10) mm` 与 `(0.18,0.12) mm` 均实现 20/20 正确恢复。第一轮全 EO joint-GN 消除了 V1 中真 EO 在固定间隙 GN1 排名 45--58 而被截断的问题。

## 6. 边界与报告规则

`(10,12)` 近邻 EO、`X01_corner_near` 和 `X02_corner_lowamp` 在 15 dB 下可出现完整非线性残差竞争：真 EO 已进入 Top-3，但另一组合的 full-refinement 残差更低。这不是 Funnel 候选筛选漏检。

因此正式结果必须分别报告：

- `true_eo_rank_gn1`：固定间隙 GN1 排名；
- `true_eo_rank_joint1`：第一次 joint-GN 后排名；
- `recall_at_3`：真 EO 是否进入最终完整精修集；
- 最终 EO、振幅和间隙误差；
- 第一、第二竞争解的残差裕量或置信标记。

论文主结论应以常规工况和明确的适用范围表述。近邻 EO 的竞争解属于可辨识性边界，应单列分析，不能把它表述为候选漏检或通过增加 K2 可以解决的问题。
