function R45_Write_ClosureReport()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','unified_closure_20260905');
S=readtable(fullfile(out,'UnifiedTransfer_PerAnchor.csv'),'TextType','string');
D=readtable(fullfile(out,'Dynamics_AllWindows.csv'),'TextType','string');
R=readtable(fullfile(out,'Ridge_PerWindowSpans.csv'));
P=readtable(fullfile(out,'IndependentStrain_PerWindow.csv'));
E=readtable(fullfile(out,'Ridge_FixedTrajectoryResponseSpread.csv'));
C=readtable(fullfile(out,'Dynamics_AllEOCandidateCosts.csv'),'TextType','string');
assert(height(D)==288 && height(R)==36 && height(P)==144);
file=fullfile(out,'R45_Final_Closure_Report_20260905.md');
fid=fopen(file,'w','n','UTF-8'); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# R4/R5 统一方法最终验证报告（2026-09-05）\n\n');
fprintf(fid,'本报告汇总本次补齐的验证，所有数值从保存的结果表自动生成。Main05、Main07、Main10 源码保持不变。旧结果保留；本报告对隐藏 gap、物理一致性和 LOBO 范围的解释优先于此前讨论。\n\n');
fprintf(fid,'## 执行完成情况\n\n');
fprintf(fid,'- 统一静态留出：R4 7 个面 × 6 个 anchor，实验静态 6 个叶片 × 16 个 anchor；各自测试 known/hidden gap 和 Common/Nearest/PC1，共 %d 行。\n',height(S));
fprintf(fid,'- 动态：两叶片各 5 个 ridge 候选及 Oracle/CommonB2/R5LOBO 三组，共 16 次、288 个窗口结果；ridge 部分为 180 个窗口结果。\n');
fprintf(fid,'- 对每个动态窗口核对时间、原始 X/T/V、传感器索引、圈号与 control 完全一致；sensor 修正保持一致，零间隙增量基线等价误差小于 1e-9 mV。\n');
fprintf(fid,'- 独立应变：两通道、两叶片时间段、18 窗、两组固定频带，共 144 行。选频不使用 BTT 频率。\n\n');
fprintf(fid,'## 统一静态验证\n\n');
fprintf(fid,'| 数据 | anchor gap | 方法 | 整面中位 NRMSE (%%) | 最大 NRMSE (%%) | 中位 gap 增量误差 (mV) |\n|---|---|---|---:|---:|---:|\n');
for data=["R4","R5_static"]
 for hidden=[false true]
  for method=["Common","Nearest","PC1"]
   a=S(S.dataset==data & S.hidden_gap==hidden & S.method==method,:);
   modes=["known","hidden"];
   fprintf(fid,'| %s | %s | %s | %.6f | %.6f | %.6f |\n',data,modes(1+hidden),method,median(a.held_surface_nrmse_pct),max(a.held_surface_nrmse_pct),median(a.gap_increment_rmse_mv));
  end
 end
end
fprintf(fid,'\n误差只在未用于定位的 gap 上计算；每个 anchor 对应整组留出 gap，因此此表不是旧审计逐波形指标的更新值。各 fold 在完整 x 网格上重拟合 PCA，目标系数未参与训练；有效评价 support 沿用原校准域。\n\n');
fprintf(fid,'R4 在不外推、隐藏倾角及 anchor gap 的条件下仍支持高精度整面迁移。实验静态数据隐藏 gap 后误差明显上升，因此不能把 R4 的精度直接推广为真实数据的普遍唯一整面恢复。PC1 的中位表现支持保留该候选，但并不证明全局最优或每个案例都优于简单方法。\n\n');
fprintf(fid,'![统一静态验证](UnifiedTransfer_KnownVsHidden.png)\n\n图：每个点对应同一目标面和 anchor 在 known/hidden gap 条件下的留出整面误差；虚线表示两者相等。对数坐标，不同点并非独立实验重复。\n\n');
fprintf(fid,'## 固定后端的 ridge 敏感性\n\n');
fprintf(fid,'| 叶片 | 同窗口最大频率跨度 (Hz) | 最大幅值跨度 (mm) | 最大 dx 跨度 (mm) | 最大 RMSE 跨度 (mV) | EO 改变窗口数 |\n|---|---:|---:|---:|---:|---:|\n');
for blade=[1 5]
 a=R(R.blade==blade,:);
 fprintf(fid,'| B%d | %.9f | %.9f | %.9f | %.6f | %d |\n',blade,max(a.frequency_span_hz),max(a.amplitude_span_mm),max(a.dx_span_mm),max(a.rmse_span_mv),sum(a.eo_count>1));
end
fprintf(fid,'\n| 叶片 | 固定 control 轨迹下候选间最大增量电压 RMS 差 (mV) |\n|---|---:|\n');
for blade=[1 5], a=E(E.blade==blade,:); fprintf(fid,'| B%d | %.6f |\n',blade,max(a.max_pair_increment_difference_rms_mv)); end
fprintf(fid,'\n该增量分析固定 control 的位移、幅值、相位和逐传感器 dg，按 Main10 相同坐标映射及 clipping 评价五个面，确认变化确实进入响应模型。动态重新优化后仍可能发生幅值/dg 补偿。具体 g0、dg、g0+dg 均保存在 Dynamics_AllWindows.csv；这些值不是非 B2 的绝对物理间隙真值。+1 mV 是敏感性范围，五点不构成连续流形证明或统计置信区间。\n\n');
fprintf(fid,'![动态敏感性](Ridge_DynamicSpans.png)\n\n图：每个窗口内五个候选响应面得到的最大值减最小值；两叶片分别绘制。\n\n');
fprintf(fid,'## Oracle / CommonB2 / 前端 LOBO\n\n');
fprintf(fid,'| 叶片 | 方法 | 平均频率 (Hz) | 平均幅值 (mm) | 平均电压 RMSE (mV) | 相对 control EO 不同窗口数 |\n|---|---|---:|---:|---:|---:|\n');
for blade=[1 5]
 for method=["oracle","commonb2","r5lobo"]
  a=D(D.blade==blade & D.candidate==method,:);
  fprintf(fid,'| B%d | %s | %.9f | %.9f | %.6f | %d |\n',blade,method,mean(a.frequency_hz),mean(a.amplitude_mm),mean(a.rmse_mv),sum(~a.eo_matches_control));
 end
end
for blade=[1 5]
 a=D(D.blade==blade & D.candidate=="oracle",:); b=D(D.blade==blade & D.candidate=="commonb2",:); c=D(D.blade==blade & D.candidate=="r5lobo",:);
 fprintf(fid,'\nB%d：LOBO 相对 CommonB2 的平均 RMSE 差 = %.6f mV；逐窗更低 %d/18；预期 Oracle ≤ LOBO < CommonB2 的逐窗排序成立 %d/18。\n\n',blade,mean(c.rmse_mv-b.rmse_mv),sum(c.rmse_mv<b.rmse_mv),sum(a.rmse_mv<=c.rmse_mv & c.rmse_mv<b.rmse_mv));
end
fprintf(fid,'三组使用相同冻结传感器修正。该修正历史上来自全叶片共享静态库，因此本试验只隔离前端系数替换收益，不是完全无目标数据依赖的端到端 LOBO。R5LOBO 的 PCA、候选构建和定位未使用目标多 gap 面；目标面仅进入 Oracle 分支及事后评价。详见 Protocol_and_DataUse.md。\n\n');
fprintf(fid,'## 独立应变与 B5 第 18 窗\n\n');
fprintf(fid,'主表使用 300–1000 Hz 固定频带；520–640 Hz 为预先存在的物理频带复核。沿用保存的时间对齐，以相同中心的三圈长度取样。FFT 定位后由应变自身进行正弦拟合，没有以 BTT 频率为搜索中心。\n\n');
fprintf(fid,'| 叶片 | 通道 | 历史 R5-BTT 减应变平均频差 (Hz) | 最大绝对频差 (Hz) | 最近整数阶一致窗口数 |\n|---|---|---:|---:|---:|\n');
for blade=[1 5]
 for channel=[1 3]
  a=P(P.blade==blade & P.channel==channel & P.band_low_hz==300,:);
  fprintf(fid,'| B%d | AI1-%02d | %.6f | %.6f | %d/18 |\n',blade,channel,mean(a.btt_minus_strain_hz),max(abs(a.btt_minus_strain_hz)),sum(a.btt_eo==a.strain_nearest_order));
 end
end
fprintf(fid,'\nB5 历史结果第 18 窗选择 EO14、约 494 Hz，而两路应变约 635 Hz、对应 EO18。此前约 32.6 Hz 的频率标准差不能解释为单纯运行状态扫频；18 窗 R5/control EO 相同也不等于 18 窗物理阶次都正确。\n\n');
fprintf(fid,'| 新三组对照 | 叶片 | 相对 AI1-03 中位绝对频差 (Hz) | 最大绝对频差 (Hz) | 阶次一致窗口数 |\n|---|---|---:|---:|---:|\n');
for blade=[1 5]
 p=P(P.blade==blade & P.channel==3 & P.band_low_hz==300,:);
 for method=["oracle","commonb2","r5lobo"]
  a=D(D.blade==blade & D.candidate==method,:); assert(isequal(a.window,p.window));
  e=abs(a.frequency_hz-p.strain_frequency_hz);
  fprintf(fid,'| %s | B%d | %.6f | %.6f | %d/18 |\n',method,blade,median(e),max(e),sum(a.eo==p.strain_nearest_order));
 end
end
fprintf(fid,'\n');
fprintf(fid,'| B5 第18窗候选面 | 选中 EO | 选中频率 (Hz) | EO14 完整代价 RMSE (mV) | EO18 完整代价 RMSE (mV) |\n|---|---:|---:|---:|---:|\n');
for method=unique(D.candidate,'stable').'
 a=D(D.blade==5 & D.window==18 & D.candidate==method,:); c=C(C.blade==5 & C.window==18 & C.surface_candidate==method,:);
 c14=c.plainRmseMv(c.EO==14); c18=c.plainRmseMv(c.EO==18);
 if isempty(c14), c14=NaN; end; if isempty(c18), c18=NaN; end
 fprintf(fid,'| %s | %d | %.6f | %.6f | %.6f |\n',method,a.eo,a.frequency_hz,c14(1),c18(1));
end
fprintf(fid,'\nEO18 已在历史候选集中，历史完整目标函数仍偏好 EO14。该异常不能归因于没有搜索 EO18。表中结果用于判定静态面替换能否改变这一选择；未用应变反向修改或挑选 Main10 输出。\n\n');
fprintf(fid,'![独立应变频率](IndependentStrain_Frequency.png)\n\n图：历史已审计的 R5 结果与两路独立选频应变逐窗比较，保留全部窗口。应变幅值仅报告原始单位；未重新验证应变片对应叶片及模态位移换算，因此不宣称绝对位移幅值和绝对间隙已获物理验证。\n\n');
fprintf(fid,'## 结论与尚未解决的科学问题\n\n');
fprintf(fid,'1. 本次规定的数值分析已执行完毕，完成不代表全部预期假设成立。\n');
fprintf(fid,'2. 建议保留统一 PC1 排序系数流形作为当前候选。算法统一为同一 gap-law/PC1/interpolation/profile 规则；R4 同域配准为恒等，实验使用 B2 冻结配准。\n');
fprintf(fid,'3. 实验 hidden-gap 整面恢复存在明显失败案例。可用响应面与唯一正确物理状态必须分开表述；不能仅凭 anchor RMSE 或传感器 gate 宣称整面正确。\n');
fprintf(fid,'4. 完全严格的端到端 LOBO 仍需排除目标数据对历史共享库、传感器修正及支持选择的影响。这会改变当前固定后端对照的前提，本次没有把这种额外验证冒充完成。\n');
fprintf(fid,'5. B5 末窗的物理错阶问题以独立应变及候选代价表为依据保留。当前工作不能宣称所有窗口物理识别已闭环。绝对幅值需独立模态换算依据；非 B2 绝对 gap 仍需独立真值。\n');
fprintf(fid,'6. 最优结论限定于已测试方法、任务、数据与评价量。静态中位误差更优不构成端到端全局最优证明，也不能代替上述失败边界。\n\n');
fprintf(fid,'## 文件\n\n');
fprintf(fid,'- UnifiedTransfer_PerAnchor.csv / UnifiedTransfer_Summary.csv：统一静态及增量验证。\n- Dynamics_AllWindows.csv / Ridge_PerWindowSpans.csv：288 窗结果与同窗跨度。\n- Dynamics_AllEOCandidateCosts.csv：所有候选阶次代价，含异常窗。\n- IndependentStrain_PerWindow.csv：独立应变逐窗数据。\n- Ridge_FixedTrajectoryResponseSpread.csv：固定轨迹下观测模型变化。\n- Protocol_and_DataUse.md：数据使用边界、历史差异和复现说明。\n');
fprintf('Report saved: %s\n',file);
end
