# local_func 说明

## 这个文件夹是做什么的

这个文件夹存放的是**底层函数库**。

这里不是给你直接点开就跑分析的地方，而是给主程序、固定可信空间流程、边界分析流程提供共用函数。

## 这里面的函数大致分成什么类型

主要有这几类：

1. **模板构建类**
   - 例如 `build_gap_template_library.m`
   - 用于构造静态响应面和相关预计算量。

2. **静态间隙初值估计类**
   - 例如 `estimate_highspeed_static_gap_raw.m`
   - 用于从高速波形中先估计静态间隙初值。

3. **变量投影与频率识别类**
   - 例如 `fit_vp_main.m`
   - 用于双频或单频振动参数的候选识别。

4. **联合优化细化类**
   - 例如 `refine_projected_joint.m`
   - 用于在初值附近做最终有界非线性优化。

5. **波形仿真与映射类**
   - 例如 `simulate_rotating_waveform_from_template.m`
   - 用于生成高速波形，或者做空间映射。

6. **汇总与辅助工具类**
   - 例如 `aggregate_multi_case_table.m`
   - 用于汇总结果、打包输出、生成表格。

## 什么时候需要看这里

只有下面这些情况才建议直接进这个文件夹看：

- 你要改底层算法；
- 你要诊断具体是哪一步算错了；
- 你要继续优化运行效率；
- 你要追某个脚本内部到底调用了哪些步骤。

## 平时不要把这里当入口

如果你只是想运行程序，不建议从这里开始。

你平时应该优先用这些入口：

- 单工况手动运行：
  `main_single_case_manual/Run_Main_Single_Case.m`

- 常规多工况主分析：
  `02_fixed_trust_domain_pipeline`

- 边界与失效分析：
  `04_vibration_amplitude_boundary`
