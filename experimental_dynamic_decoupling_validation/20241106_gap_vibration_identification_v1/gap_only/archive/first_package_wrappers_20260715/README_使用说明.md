# 20241106 间隙—振动辨识独立程序包

本目录是新的独立版本。原来的
`20241106_low_speed_gap_prior_decoupling` 和
`20241106_btt_data_foundation` 没有覆盖、没有改名，继续作为历史基线保留。

## 最简单的使用方式

通常只需要修改根目录的 `Config_20241106.m`，然后选择一个入口运行：

1. `Run_01_Prepare_Data_20241106`：准备并检查公共输入；
2. `Run_02_Foundation_NoGap_20241106`：只运行不考虑高速间隙变化的方法；
3. `Run_03_GapAware_20241106`：运行考虑高速间隙变化的方法；缺少前序结果时会自动补算；
4. `Run_04_Compare_Results_20241106`：生成两种方法的逐窗口对比；缺少前序结果时会自动补算。

四个入口均可从MATLAB命令行直接调用，不需要进入`src`目录，也不需要手动设置环境变量。

## 文件夹结构

```text
20241106_gap_vibration_identification_v1
├─ Config_20241106.m                 唯一正式配置
├─ Check_Program_20241106.m          完整性和独立性检查
├─ Run_01_Prepare_Data_20241106.m    数据准备入口
├─ Run_02_Foundation_NoGap_20241106.m
├─ Run_03_GapAware_20241106.m
├─ Run_04_Compare_Results_20241106.m
├─ inputs
│  ├─ calibration                    低速模板与间隙标定库
│  └─ prepared                       随包保存的公共准备数据
├─ src
│  ├─ foundation                     不考虑间隙的方法核心
│  ├─ preparation                    数据准备与DynamicMap构造
│  ├─ gap_aware                      考虑间隙的方法核心
│  └─ utilities
└─ results
   ├─ prepared                       当前工况的文件夹内输入
   ├─ foundation_no_gap              不考虑间隙结果
   ├─ gap_aware                      考虑间隙结果
   ├─ comparison                     两种方法对比
   └─ figures
```

## 两种辨识方法的边界

### Foundation：不考虑高速间隙变化

正演模型中不释放`dg`，作为固定间隙对比方法独立运行。该方法可以出现混合EO，不能把它提前选出的EO锁死到GapAware方法中。

### GapAware：考虑高速间隙变化

每个3圈窗口独立执行：

```text
加权伪位移VP扫描
→ 保留Top-K EO候选
→ 每个候选分别进行间隙感知正演
→ 连续频率细化
→ 按最终普通电压RMSE选择EO、幅值、相位、dx和dg
```

正式设置为：

- `eta=0`；
- 不强制EO12；
- 不跨窗口共享EO、频率或相位；
- VP权重只用于伪位移初始化和候选筛选；
- 最终电压波形残差不加权；
- 普通RMSE模式下不再使用加权固定模型误差下限；
- 只有配置中的CH5、CH7释放间隙增量。

## 当前回归结果

默认配置为B4、CH2/CH5/CH7、86.6 s、20圈输入、3圈窗口、1圈步长，共18个窗口。

- Foundation：EO12为14/18窗口，其余出现EO6、EO13和EO17；
- GapAware：18/18窗口均为EO12；
- GapAware平均频率：631.93 Hz；
- GapAware平均幅值：0.357970 mm；
- GapAware平均普通电压RMSE：38.381832 mV。

这些数值是新目录独立运行得到的回归检查点，不是对所有时间区域的物理结论。

## 数据说明

代码、标定库和已准备输入均位于新目录，不读取旧辨识目录的程序或结果。高速原始波形仍通过公共试验数据路径读取；如原始数据位置改变，只修改`Config_20241106.m`中的`cfg.paths.rawDatasetRoot`。

运行完整检查：

```matlab
Check_Program_20241106('all', true)
```
