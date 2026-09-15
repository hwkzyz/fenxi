# R5 正式试验结果清单

生成日期：2026-09-07  
正式方法：`R5_PC1_AnchorSurface_SensorWiseDg_FullWave`  
统一检查：`Check_Latest_Program_20260905` = PASS

## 1. 正式结果文件

| case | 正式辨识结果 | Foundation 公共窗口包 | 应变验证 |
|---|---|---|---|
| B1 / R01 | `results/gap_aware/Main_GapAware_VPTopK_FullWave_20251222_B1_S123_r01_gapaware.mat` | `inputs/prepared/foundation/step05_r01/1000_2500_3500/Result_Step05_NoEtaVPTopKDirectTemplate_B1_S123_20251222.mat` | `results/strain_validation/StrainBTT_B1_StrainReferenceValidation.csv` |
| B5 / R04 | `results/gap_aware/Main_GapAware_VPTopK_FullWave_20251222_B5_S123_r04_gapaware.mat` | `inputs/prepared/foundation/step05_r04/1000_2500_3500/Result_Step05_NoEtaVPTopKDirectTemplate_B5_S123_20251222.mat` | `results/strain_validation/StrainBTT_B5_StrainReferenceValidation.csv` |

两套结果均为 18/18 窗口，HardGate = PASS，SensorWiseDg audit = pass，Zero-gap equivalence = 0 mV。

## 2. Foundation 对比

| case | method | EO dominant fraction | mean frequency (Hz) | frequency SD (Hz) | mean RMSE (mV) |
|---|---|---:|---:|---:|---:|
| B1 / R01 | Foundation direct low-speed template | 13/18 | 570.981 | 19.181 | 62.070 |
| B1 / R01 | R5 formal method | 18/18 EO14 | 582.554 | 0.488 | 33.927 |
| B5 / R04 | Foundation direct low-speed template | 13/18 | 626.684 | 188.037 | 58.477 |
| B5 / R04 | R5 formal method | 18/18 EO18 | 632.121 | 1.782 | 35.597 |

对应的正式对齐比较表位于 `results/gap_aware/comparison/`。

## 3. 应变验证口径

时间转换固定为：

```text
t_BTT = t_strain_raw + strain_time_offset_total_sec - 30.045883997 s
```

| case | validation scope | BTT mean frequency (Hz) | strain frequency range (Hz) | mean BTT-strain difference (Hz) |
|---|---|---:|---:|---:|
| B1 / R01 | 同目标叶片直接交叉验证 | 582.5545 | AI1-03: 582.4806-582.7806 | -0.0650 |
| B5 / R04 | B1 应变片的同工况跨模态验证 | 632.1209 | AI1-03: 629.8139-634.1139 | 0.1737 |

B5 的 AI1-03 幅值趋势 Pearson 相关系数为 0.9516；AI1-01 为 0.7698。B1 振幅变化较小，不能将其接近零的趋势相关系数解释为方法失效，应以频率一致性和绝对幅值差为主。

## 4. 可用于论文的结论边界

- Foundation 表示既有的“直接低速波形模板、无动态间隙变化”思路。
- R5 方法在两组正式工况中均降低平均 RMSE，并将 EO 和频率稳定为单一支路。
- B1 可作为同叶片应变交叉验证；B5 只能表述为同工况、独立应变通道支持的动态一致性验证，不能称为 B5 位移真值。
- 动态间隙参数的必要性应结合 RMSE、局部边界无命中、参数连续性和 Zero-gap equivalence 共同论证。

## 5. 不纳入正式结果的内容

- `shift1`、`nested_compare`、`r04fix` 等历史诊断文件。
- 旧 `Z:` 路径记录或零字节 MAT 文件。
- `gap_only`：仅为 MAT 内部模型分支字段，不是正式方法名称。
