# `inv_log_2` 间隙-振动解耦方法

## 1. 问题定义

本文针对叶尖定时（blade tip timing, BTT）观测中的间隙识别与叶尖振动重构问题，提出 `inv_log_2` 间隙-振动解耦方法。设静态间隙为 \(g\)，空间坐标为 \(x\)，静态无振动波形为 \(F_g(x)\)。在高速工况下，叶尖在时刻 \(t\) 的实际空间位置为

\[
x(t)=x_0(t)-u(t),
\]

其中 \(x_0(t)\) 为名义扫掠坐标，\(u(t)\) 为叶尖径向振动位移。相应的高速观测模型写为

\[
V(t)=F_g(x_0(t)-u(t))+\varepsilon(t),
\]

其中 \(\varepsilon(t)\) 为噪声及未建模误差。由此，观测信号同时包含静态间隙效应和振动扰动效应，二者在高速数据中呈耦合状态。本文的目标是在有限静态标定库基础上重构连续静态模板族，并据此实现间隙初值估计、振动投影修正与联合精修。

## 2. 静态模板族重构

静态标定库由若干离散间隙下的无振动波形构成：

\[
\mathcal D_s=\{F_{g_i}(x),\;i=1,\dots,N\}.
\]

由于高速工况的真实间隙 \(g^\ast\) 一般不与标定间隙 \(g_i\) 重合，故需由离散库构造连续模板族

\[
\widehat F(g,x).
\]

在固定空间位置 \(x\) 上，可将静态波形视为间隙变量的平滑函数 \(v_x(g)=F_g(x)\)。在研究区间内，采用低维基函数展开对其进行参数化：

\[
\widehat F(g,x)\approx B_0(x)+B_1(x)\phi_1(g)+B_2(x)\phi_2(g).
\]

本文采用的间隙基函数为

\[
\phi(g)=\left[1,\frac{1}{g},\log g\right].
\]

因此，连续模板族写为

\[
\widehat F(g,x)=B_0(x)+B_1(x)\frac{1}{g}+B_2(x)\log g.
\]

其中，常数项用于吸收基线和偏置项，\(1/g\) 表示理想平行场主导变化，\(\log g\) 表示有限尺寸和边缘场引起的弱非线性修正。

对静态标定库中的每个空间采样点 \(x_j\)，令

\[
y_j=
\begin{bmatrix}
F_{g_1}(x_j)\\
F_{g_2}(x_j)\\
\vdots\\
F_{g_N}(x_j)
\end{bmatrix},
\qquad
\Phi=
\begin{bmatrix}
1 & 1/g_1 & \log g_1\\
1 & 1/g_2 & \log g_2\\
\vdots & \vdots & \vdots\\
1 & 1/g_N & \log g_N
\end{bmatrix},
\]

则系数向量 \(b_j=[B_0(x_j),B_1(x_j),B_2(x_j)]^T\) 通过最小二乘估计：

\[
b_j^\ast=\arg\min_b\|\Phi b-y_j\|^2.
\]

由此得到随空间位置变化的系数函数 \(B_0(x)\)、\(B_1(x)\) 和 \(B_2(x)\)，并据此构造连续模板族 \(\widehat F(g,x)\)。

## 3. 可信空间确定

并非整条波形都适于间隙辨识。叶尖通过波形的两端通常接近基线，局部斜率和曲率较小，间隙变化在该区域内的可分性不足。为此，需要先从静态库中确定可信空间 \(\Omega_x=[x_L,x_R]\)。

可信空间的确定基于三项指标：

1. 区间内响应活动应足够显著；
2. 不同间隙在该区间内应具有足够区分度；
3. 留一间隙验证中的模板误差和导数误差应较小。

对候选窗口 \(\Omega_x\)，定义模板误差、导数误差和间隙敏感度项，并构造综合评分

\[
S(\Omega_x)=E_F(\Omega_x)+\lambda_dE_{F'}(\Omega_x)+\lambda_sR_{\mathrm{sens}}(\Omega_x).
\]

其中 \(E_F\) 为模板重构误差，\(E_{F'}\) 为导数误差，\(R_{\mathrm{sens}}\) 为间隙敏感度不足的惩罚项。最终选择在误差可接受条件下尽可能宽的窗口，以兼顾信息保留与数值稳定性。

## 4. 高速观测模型

对高速旋转工况，第 \(m\) 圈、第 \(k\) 个叶片的叶尖通过时刻记为 \(T_{m,k}\)。若转子半径为 \(R\)，角速度为 \(\Omega\)，则名义扫掠坐标为

\[
x_0(t)=R\Omega(t-T_{m,k}).
\]

考虑叶尖径向振动位移 \(u(t)\) 后，实际采样位置变为 \(x_0(t)-u(t)\)，从而高速观测可写为

\[
V_{\mathrm{hs}}(t)=\widehat F(g,x_0(t)-u(t))+\varepsilon(t).
\]

该模型表明，振动首先通过空间坐标扰动进入观测过程，其次才体现为波形形状变化。因此，后续修正应优先针对坐标扰动建模。

## 5. 静态间隙初值估计

在可信空间 \(\Omega_x\) 内，首先忽略振动项，仅以静态模板匹配高速波形主轮廓。引入空间平移参数 \(d_x\)，构造静态目标函数

\[
J_{\mathrm{static}}(g,d_x)=\frac{1}{M}\sum_{j=1}^{M}\bigl[V(x_j)-\widehat F(g,x_j-d_x)\bigr]^2,\quad x_j\in\Omega_x.
\]

其中 \(d_x\) 用于吸收名义坐标与实际对齐之间的整体偏差。通过先对 \(d_x\) 最小化，再对 \(g\) 搜索，可避免将对齐误差误判为间隙变化。由于该目标函数通常为非凸形式，采用二维网格或分阶段搜索比直接局部优化更稳健。

记静态初值为 \(\hat g_0\)。该值仅用于为后续频率识别和振动修正提供合适的起点，并不代表最终间隙结果。

## 6. 振动投影修正

静态初值仅能解释波形主轮廓，不能解释由振动引起的细节偏移。若振动幅值相对于波形空间尺度较小，则可对模板作一阶展开：

\[
\widehat F(g,x_0(t)-u(t))\approx \widehat F(g,x_0(t))-\widehat F_x'(g,x_0(t))u(t),
\]

其中

\[
\widehat F_x'(g,x)=\frac{\partial \widehat F(g,x)}{\partial x}.
\]

若位移可表示为少量谐波叠加

\[
u(t)=\sum_{q=1}^{Q}A_q\sin(\omega_q t+\varphi_q),
\]

则一阶振动项属于由 \(\widehat F_x'\) 及其正弦、余弦调制项张成的有限维子空间。由此，定义投影残差

\[
r_\perp(t;g)=\bigl(I-\mathcal P_{\mathcal B}\bigr)\bigl[V(t)-\widehat F(g,x_0(t))\bigr],
\]

其中 \(\mathcal B\) 为振动解释子空间，\(\mathcal P_{\mathcal B}\) 为到该子空间的最小二乘投影算子。

该步骤的作用在于分离能够由振动模型解释的残差成分，从而减弱振动对间隙搜索的污染。投影修正并非一次性的后处理，而是对每个候选间隙均需执行的评价环节。

## 7. 联合辨识与变量投影

在静态初值和投影修正基础上，进一步联合求解间隙参数与振动参数。设振动参数记为 \(\theta_u\)，则联合模型为

\[
V(t)\approx \widehat F(g,x_0(t)-u(t;\theta_u)).
\]

若直接同时搜索 \(g\) 与 \(\theta_u\)，则问题维数较高且耦合显著。为降低搜索复杂度，采用变量投影思想：外层搜索间隙，内层对振动引入的线性或准线性子问题进行最小二乘消元。这样可将非线性部分主要保留在间隙维度上，而幅值、相位与频率参数在每个候选间隙下解析或半解析求解。

在实现上，频率列表可先由固定间隙下的变量投影识别得到，或由静态残差频谱初筛得到，再带入间隙修正模型中继续精化。由此，振动投影与联合辨识形成闭环：先稳定间隙，再提升振动参数识别精度。

## 8. 位移重构与结果解释

当最终获得 \(\hat g\) 与 \(\hat\theta_u\) 后，可重构叶尖径向位移时程

\[
\hat u(t)=u(t;\hat\theta_u).
\]

该结果可用于两个层面的验证：

1. 时间域上，检查振动的幅值、相位和频率是否合理；
2. 空间域上，代回 \(\widehat F(\hat g,x_0(t)-\hat u(t))\) 检查是否能够重构高速 BTT 观测。

因此，位移重构是间隙-振动联合模型的最终输出，而非独立附加步骤。若重构波形在可信空间内与观测数据保持一致，则表明当前模型下估计的间隙和振动参数具有内部一致性。

值得注意的是，间隙估计更准确并不必然意味着后续振动辨识更优。原因在于，间隙目标函数主要约束静态模板主轮廓，而振动重构还依赖模板导数、频率候选、谐波分离以及最终非线性联合优化。两类目标函数所强调的物理量并不完全相同。

## 9. `inv_log_2` 基函数模型

`inv_log_2` 的核心在于将间隙变量上的连续响应面表示为低阶可验证模型。与高阶多项式相比，\(1/g\) 和 \(\log g\) 更适合刻画小间隙段的快速变化与较大间隙段的缓变过渡；与纯样条插值相比，该形式更便于进行快速搜索和导数分析。

本文将 `inv_log_2` 解释为“反比-对数二项响应面”方法，即以 \(1/g\) 和 \(\log g\) 为核心的双基函数间隙模型，并辅以常数项吸收背景与偏置。其适用性以静态库上的留一验证、模板误差和导数误差为判据，而非由先验物理定律直接给出。

## 10. 方法边界

该方法适用于以下条件：

1. BTT 信号可稳定获取，且同一测点下的静态库已建立；
2. 间隙响应在研究区间内近似单调且平滑，模板族可连续化；
3. 振动幅值相对静态模板空间变化不大，小扰动一阶近似仍然适用；
4. 主要振动可由少量谐波或少数频率成分描述。

当振动幅值过大、模态过多、波形失真严重或间隙响应出现强非线性拐点时，单一 `inv_log_2` 模板可能不足以覆盖全部情形，需要扩展间隙基函数、引入分段模板或增加振动子空间维数。

## 11. 相关文献依据

与本文相关的代表性文献包括：

1. Heath, S., and Imregun, M. A survey of blade tip-timing measurement techniques for turbomachinery vibration. *Journal of Engineering for Gas Turbines and Power*, 1998.
2. Heath, S. A new technique for identifying synchronous resonances using tip-timing. *Journal of Engineering for Gas Turbines and Power*, 2000.
3. Carrington, I., et al. A comparison of blade tip timing data analysis methods. *Proceedings of the Institution of Mechanical Engineers, Part G*, 2001.
4. Dimitriadis, G., Carrington, I., Wright, J., and Cooper, J. Blade-tip timing measurement of synchronous vibrations of rotating bladed assemblies. *Mechanical Systems and Signal Processing*, 2002, DOI: 10.1006/mssp.2002.1489.
5. Lin, et al. Sparse reconstruction of blade tip-timing signals for multi-mode blade vibration monitoring. *Mechanical Systems and Signal Processing*, 2016, DOI: 10.1016/j.ymssp.2016.03.020.
6. An improved multiple signal classification for nonuniform sampling in blade tip timing. *IEEE Transactions on Instrumentation and Measurement*, 2020, DOI: 10.1109/TIM.2020.2980912.
7. Golub, G. H., & Pereyra, V. Separable nonlinear least squares: the variable projection method and its applications. *Inverse Problems*, 2003.

上述文献表明，BTT 的核心任务是由离散叶尖观测反推叶片振动。本文在此基础上进一步引入静态间隙库、可信空间筛选和振动投影修正，以处理间隙与振动强耦合的辨识场景。
