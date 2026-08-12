# trust_domain_result 说明

## 这个文件夹是做什么的

这个文件夹存放的是**可信空间标定结果文件**。

它不是脚本目录，也不是运行入口目录，而是 `01_trust_domain_calibration` 计算后的结果缓存目录。

## 这里面的文件代表什么

- `inv_log_2_trust_domain.mat`
  - 保存可信空间的核心结果，供后续固定可信空间流程读取。

- `inv_log_2_trust_domain_summary.csv`
  - 保存可信空间筛选的简要汇总结果。

- `inv_log_2_trust_domain_metrics.csv`
  - 保存可信空间相关的误差指标。

## 谁会用到这里

主要是这些目录会读取这里的结果：

- `02_fixed_trust_domain_pipeline`
- `03_dx_true_injection_pipeline`
- 单工况主程序中与固定可信空间相关的部分

## 什么时候需要关心这里

你只在下面两种情况下需要关心这里：

1. 想确认当前可信空间结果有没有成功生成；
2. 想看 `01_trust_domain_calibration` 输出了哪些文件。

平时不需要直接在这里运行任何程序。
