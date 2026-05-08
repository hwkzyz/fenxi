# 场效应传感器响应模型文献整理

## 1. 结论先行

目前检索到的文献并不支持把叶片经过电容或电涡流传感器时的输出电压简单写成一个固定解析式，例如：

```math
V(g,x)=a+b/g
```

更常见、也更严谨的表达链条是：

```text
电场/磁场作用机理
→ 电容 C 或线圈阻抗 Z / 电流 I 的变化
→ 调理电路输出电压/频率/幅值
→ 通过标定曲线或拟合模型反算 tip clearance
```

因此，对当前工作的启发是：

```text
1/g 可以作为理想平板电容的原理近似；
但真实波形 F_g(x) 应被看成“标定得到的传感器响应模板”，
而不是由单一解析式直接决定。
```

## 2. 电容式传感器

### 2.1 理想电容表达

综述文献通常用平板电容公式解释电容法测间隙的基本原理：

```math
C=\frac{\varepsilon_0\varepsilon_r A}{d}
```

其中 `A` 是探头有效面积，`d` 是探头到叶尖的距离。这个公式说明电容与距离在理想条件下呈反比关系，因此 `1/g` 插值有物理来源。

但同一类综述也强调，实际系统需要先得到传感器特性曲线，常用查表、最小二乘、Lagrange 插值或样条插值等方式把输出信号映射到间隙。

代表文献：

- Yu et al., 2020, *A review of blade tip clearance-measuring technologies for gas turbine engines*, Measurement and Control, DOI: `10.1177/0020294019877514`.

### 2.2 实际电容-间隙关系通常作为标定函数

Fabian et al. 的微型燃气轮机电容传感器工作非常关键。该文献并没有把最终间隙估计建立在单一 `1/d` 解析式上，而是明确把探头电容看成间隙的非线性函数：

```math
\hat C_e=f_c(c)
```

间隙估计通过反标定函数得到：

```math
\hat c=f_c^{-1}(\hat C_e)
```

这对我们非常重要。它说明，即使电容法有 `1/d` 的基本物理来源，实际工程中仍然通过校准函数把电容或电压输出转换为间隙。

代表文献：

- Fabian, Prinz, Brasseur, 2005, *Capacitive Sensor for Active Tip Clearance Control in a Palm-Sized Gas Turbine Generator*, IEEE Transactions on Instrumentation and Measurement, DOI: `10.1109/TIM.2005.847233`.

### 2.3 电路输出与电容的关系也要单独处理

Sarma and Barranger 的 NASA/IEEE 工作讨论了电容式叶尖间隙测量的调理电路。文献指出，叶尖与探头形成的电容是运放反馈元件之一；在特定时间约束下：

```text
斜坡输入时，输出电压与静态叶尖间隙电容成比例；
直流输入时，输出与间隙电容导数成比例，再积分恢复动态电容。
```

这说明实际测得的电压不是直接等于 `1/g`，而是经过具体电路结构后才与电容或电容导数产生关系。

代表文献：

- Sarma and Barranger, 1992, *Capacitance-type blade-tip clearance measurement system using a dual amplifier with ramp/DC inputs and integration*, IEEE Transactions on Instrumentation and Measurement, DOI: `10.1109/19.177341`.
- Barranger, 1987, *Low-cost FM oscillator for capacitance type of blade tip clearance measurement system*, NASA-TP-2746.

### 2.4 近年电容法更常把波形建成脉冲模型

近年的电容同步测量文献常常不直接写电场解析式，而是把叶尖经过探头时的电压信号建成高斯脉冲或对称函数。间隙与脉冲幅值通过预标定关系对应，TOA 则由峰值时间或对称约束估计。

这与我们当前的“静态模板波形库”思路接近：重点不是假定 `V=1/g`，而是利用标定波形或标定幅值关系来重构/反算间隙。

近期线索：

- Li et al., 2026, *Capacitive sensing-based method for blade tip clearance and timing measurement through symmetric inverse function fitting in aero-engines*, Aerospace Science and Technology, DOI: `10.1016/j.ast.2025.111430`.
- Cao et al., 2025, *A High-Precision Measurement Method for Blade Tip Vibration and Clearance Based on Capacitive Sensor*, IEEE Transactions on Instrumentation and Measurement, DOI: `10.1109/TIM.2025.3561389`.

这两篇属于较新的文献线索，正式写论文时建议以 DOI 页面和出版社页面再次核对出版状态；主理论支撑仍应优先使用综述、NASA/IEEE 电容电路和已稳定发表的标定类文献。

## 3. 电涡流/电感式传感器

### 3.1 基本链条是电磁场到线圈阻抗

电涡流或电感式传感器通常用如下机理描述：

```text
线圈激励产生交变磁场；
金属叶尖经过磁场后感生涡流；
涡流产生反向磁场；
线圈等效阻抗、感抗或输出电压发生变化；
通过标定曲线得到间隙。
```

综述文献会把等效阻抗写成与距离、导磁率、电阻率、线圈和目标几何、激励频率有关的函数：

```math
Z=Z(d,\mu,\rho,\text{geometry},\omega)
```

因此，电涡流响应本身比电容理想平板模型更难简化为某个固定 `1/g` 或 `g^{-p}` 形式。

代表文献：

- Yu et al., 2020, *A review of blade tip clearance-measuring technologies for gas turbine engines*, Measurement and Control, DOI: `10.1177/0020294019877514`.

### 3.2 更严格的电涡流模型通常走 Maxwell 方程或有限元

Jamia et al. 对电涡流叶尖定时传感器输出进行了二维/三维有限元仿真。他们从 Maxwell 方程、材料关系、运动导体中的 Ohm 定律出发，建立磁准静态模型，再用 COMSOL 模拟传感器输出。

文献中的关键结论是：

```text
电涡流传感器输出同时受间隙、传感器位置、叶片几何和转速影响；
因此真实输出不应只看成单一距离函数。
```

这直接支持我们对 `F_g(x)=H(g+s(x))` 的降级判断：对于场效应传感器，`x` 不只是等效距离平移项，还可能改变场耦合状态本身。

代表文献：

- Jamia et al., 2018, *Simulating eddy current sensor outputs for blade tip timing*, Advances in Mechanical Engineering, DOI: `10.1177/1687814017748020`.

### 3.3 高温电感/电涡流传感器仍强调标定

Liu/Zhao/Lyu 等关于高温电感式传感器的研究给出了等效电路、电压分压关系和线圈磁场计算，同时通过实验标定评价温度和间隙影响。该类文献的共同点是：

```text
理论模型解释工作原理；
具体间隙估计仍依赖实验标定曲线或拟合修正。
```

代表文献：

- Liu et al., 2019, *Experimental Investigation of Inductive Sensor Characteristic for Blade Tip Clearance Measurement at High Temperature*, Sensors, DOI: `10.3390/s19173694`.
- Zhao et al., 2019, *Experimental Investigation of High Temperature-Resistant Inductive Sensor for Blade Tip Clearance Measurement*, Sensors, DOI: `10.3390/s19010061`.
- Wu et al., 2019, *Eddy Current Sensor System for Blade Tip Clearance Measurement Based on a Speed Adjustment Model*, Sensors, DOI: `10.3390/s19040761`.

## 4. 对当前理论的修正建议

基于上述文献，当前理论应避免这样的强表述：

```text
传感器输出电压服从某个确定的 1/g、g^{-p} 或 exp(-g/lambda) 形式。
```

更稳妥的表达是：

```text
电容式传感器在理想平板极限下给出 C ∝ 1/d 的物理来源；
电涡流/电感式传感器的输出则由电磁场与线圈阻抗决定；
但实际叶尖经过探头时的电压波形受到探头几何、叶尖形状、材料、电路和安装状态共同影响。
因此工程上通常通过标定曲线、查表、插值、拟合或有限元辅助模型建立输出与间隙之间的映射。
```

这意味着我们的静态模板库方法是合理的：

```math
V(t)=F_g(x(t)-u(t))
```

其中 `F_g(x)` 不应被解释为从某个简单解析式直接推出，而应解释为：

```text
在给定传感器、叶尖几何、安装状态和电路条件下，
由静态标定得到的全波形响应模板。
```

## 5. 对当前二阶幂律模型的定位

当前最优的

```math
F_g(x)\approx A_0(x)+A_1(x)g^{-p}+A_2(x)g^{-(p+1)}
```

不应写成“由电场理论严格推出”。更准确的定位是：

```text
受电容理想反比响应和场效应传感器非线性标定思想启发的低维等效响应展开。
```

它的合理性来自两层证据：

1. 物理上，电容法存在 `C ∝ 1/d` 的理想极限，电涡流法也具有明显非线性距离响应；
2. 数据上，该模型在静态模板留一重构和高速闭环识别中优于当前 `1/g` 基线。

因此论文中建议使用如下表述：

```text
考虑到场效应传感器输出通常需要通过标定曲线建立与间隙的非线性映射，
本文不预设简单的解析电压-间隙关系，
而是将静态波形族视为由标定模板定义的响应流形。
在该流形上，采用物理启发的低维幂律展开来提高有限模板库下的连续间隙重构精度。
```

## 6. 可引用文献清单

1. Yu, B., Ke, H., Shen, E., Zhang, T. (2020). *A review of blade tip clearance-measuring technologies for gas turbine engines*. Measurement and Control. DOI: `10.1177/0020294019877514`.
2. Sarma, G. R., Barranger, J. P. (1992). *Capacitance-type blade-tip clearance measurement system using a dual amplifier with ramp/DC inputs and integration*. IEEE Transactions on Instrumentation and Measurement. DOI: `10.1109/19.177341`.
3. Barranger, J. P. (1987). *Low-cost FM oscillator for capacitance type of blade tip clearance measurement system*. NASA-TP-2746.
4. Fabian, T., Prinz, F. B., Brasseur, G. (2005). *Capacitive Sensor for Active Tip Clearance Control in a Palm-Sized Gas Turbine Generator*. IEEE Transactions on Instrumentation and Measurement. DOI: `10.1109/TIM.2005.847233`.
5. Salinas, S., Castillo, A., Bloxham, M., Paniagua, G. (2023). *Comprehensive capacitance sensor calibration for high-speed fluid-machinery tip clearance characterization*. Measurement. DOI: `10.1016/j.measurement.2023.113117`.
6. Jamia, N., Friswell, M. I., El-Borgi, S., Fernandes, R. (2018). *Simulating eddy current sensor outputs for blade tip timing*. Advances in Mechanical Engineering. DOI: `10.1177/1687814017748020`.
7. Liu, Z., Zhao, Z., Lyu, Y., Zhao, L. (2019). *Experimental Investigation of Inductive Sensor Characteristic for Blade Tip Clearance Measurement at High Temperature*. Sensors. DOI: `10.3390/s19173694`.
8. Zhao, Z., Liu, Z., Lyu, Y., Gao, Y. (2019). *Experimental Investigation of High Temperature-Resistant Inductive Sensor for Blade Tip Clearance Measurement*. Sensors. DOI: `10.3390/s19010061`.
9. Wu, J., Wen, B., Zhou, Y., Zhang, Q., Ding, S. et al. (2019). *Eddy Current Sensor System for Blade Tip Clearance Measurement Based on a Speed Adjustment Model*. Sensors. DOI: `10.3390/s19040761`.
10. Li, F., Liu, H., Duan, F., Guo, G., Teng, G., Zhou, X., Zhou, Q., Xiao, F. (2025/2026). *Capacitive sensing-based method for blade tip clearance and timing measurement through symmetric inverse function fitting in aero-engines*. Aerospace Science and Technology. DOI: `10.1016/j.ast.2025.111430`. 近期线索，引用前需再次核对出版状态。
11. Cao, H., Chen, Y., Zhang, Y., Yuan, M., Li, H. (2025). *A High-Precision Measurement Method for Blade Tip Vibration and Clearance Based on Capacitive Sensor*. IEEE Transactions on Instrumentation and Measurement. DOI: `10.1109/TIM.2025.3561389`. 近期线索，引用前需再次核对出版状态。
