# 预先分析协议

## P1：动态相位信息

H1：在总灵敏度能量相当时，passage 内相位分散使单谐波正余弦 Gram 矩阵的最小特征值增大。

H2：对至少一类近退化 EO 对，exact passage 比 frozen passage 提供更大的候选子空间主角或受约束模型距离。

\[
r_2=\frac{|\sum_n w_n e^{i2\psi_n}|}{\sum_nw_n},\qquad
\lambda_\pm=\frac{W_0}{2}(1\pm r_2).
\]

同时报告每传感器和联合的 `r2`、`lambda_min_gram`、`condition_gram`、主角与模型距离。

对照分为：固定权重的纯相位 exact/frozen、允许 `F_x(x-u(t))` 变化的完整物理对照，以及 exact 数据由 frozen 模型拟合的模型失配对照。

**Gate 1：**如果 exact 仅改变未条件化 `r2`，但条件信息、EO 主角和模型距离基本不变，则不再将 trajectory-segment 作为主要科学主张。

## P2：局部条件可恢复性

第一轮只在正确 EO、正确前向模型和无噪/高 SNR 下分析连续参数：

\[
F_{p|\nu}=\widetilde J_p^TQ_\nu\widetilde J_p.
\]

目标坐标优先使用间隙与谐波正余弦系数，不在弱振幅附近直接对相位做 CRLB 解释。

**Gate 2：**若缩放后的最小条件信息不能预测正确 EO 条件下的经验误差趋势，则先检查参数尺度、非线性、边界和优化收敛，不进入盲测分类。

## P3：离散模型分离与偏差竞争

\[
d_{r\rightarrow c}(\theta_r^*)=
\inf_{\theta_c\in\Theta_c^{adm}}
\|F_r(\theta_r^*)-F_c(\theta_c)\|_W.
\]

受约束参数域必须声明振幅下界、支撑域、分量排序、间隙/偏移范围及 nuisance 是否重新优化。

\[
\beta_{rc}=|\langle W^{1/2}b,h_{rc}\rangle|,
\qquad M_{rc}=R_c-R_r.
\]

`beta_rc` 只作为局部偏差方向诊断，最终翻转以完整重优化后的 `M_rc` 判断。

**Gate 3：**若 `beta_rc / d_rc` 只能事后解释 X01，不能预测未参与指标构造的近退化 EO 对，则将它限定为案例诊断，不提升为一般失效定律。

