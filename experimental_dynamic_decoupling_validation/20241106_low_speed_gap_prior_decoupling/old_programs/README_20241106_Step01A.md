# 20241106 Step01A 脉冲提取说明

脚本：

```text
Step01A_Extract_OPR_Probe_From_ChannelFiles_20241106.m
```

这个脚本用于从 `4-1-1000.mat`、`4-2-1000.mat` 等通道文件中重新检测 OPR 和各叶尖传感器脉冲。它是对原始 `maichongcompute.m` 的简化重写，便于调试和阅读。

## 默认行为

默认只打印诊断并打开图，不覆盖已有文件。

```matlab
saveRecomputedFiles = false;
```

如果需要保存，会保存为带后缀的新文件，不覆盖旧文件：

```text
jiluOPR_recomputed_step01a.mat
jilublade_probe2_recomputed_step01a.mat
...
```

## 通道定义

当前按 20241106 目录里的 `maichongcompute.m` 设置：

```matlab
oprChannel = 1;
probeChannels = [2 3 4 5 6 7];
```

注意：这里的 OPR 通道是 20241106 这批数据自己的通道定义，不沿用 20250527/20251222 的配置。

## 阈值

顶部集中设置：

```matlab
thresholdByChannel(1) = 2.0;   % OPR
thresholdByChannel(2) = -0.2;
thresholdByChannel(3) = -0.2;
thresholdByChannel(4) = 2.0;
thresholdByChannel(5) = 1.0;
thresholdByChannel(6) = 2.0;
thresholdByChannel(7) = 1.0;
```

如果某个通道漏检或误检，优先调对应阈值。

## 快速调试

先只读取少数 block：

```matlab
maxBlocksToRead = 1;
```

全量处理时改回：

```matlab
maxBlocksToRead = [];
```

也可以不改脚本，直接用环境变量临时指定：

```matlab
setenv('STEP01A_MAX_BLOCKS','1')
Step01A_Extract_OPR_Probe_From_ChannelFiles_20241106
setenv('STEP01A_MAX_BLOCKS','')
```

## 只处理指定工况、通道或窗口

脚本顶部有这些选择项：

```matlab
selectedCaseIndex = [];        % []=低速和高速都处理；1=只低速；2=只高速
selectedProbeChannels = [2 3 4 5 6 7];
selectedBlockRange = [];       % []=全部 block；[36000 36000]=只处理 4-x-36000.mat
selectedTimeRangeSec = [];     % []=不按时间截取；[80 90]=只处理 80-90 s
```

也可以用环境变量临时覆盖：

```matlab
setenv('STEP01A_CASE_INDEX','2')          % 只处理高速 3000_3150
setenv('STEP01A_PROBE_CHANNELS','5 7')    % 只处理 probe 5 和 7
setenv('STEP01A_BLOCK_RANGE','36000 36000')
setenv('STEP01A_TIME_RANGE_SEC','70 72')
Step01A_Extract_OPR_Probe_From_ChannelFiles_20241106
```

清空环境变量：

```matlab
setenv('STEP01A_CASE_INDEX','')
setenv('STEP01A_PROBE_CHANNELS','')
setenv('STEP01A_BLOCK_RANGE','')
setenv('STEP01A_TIME_RANGE_SEC','')
```

注意：如果只处理 `4-x-36000.mat` 这种中间 block，脚本不会执行开头 20000 点裁剪；裁剪只对真实第一块 `4-x-1000.mat` 生效。

## 图怎么看

OPR 图：

- 第一幅：由重新检测 OPR 得到的转速曲线。
- 第二幅：新检测结果和已有 `jiluOPR.mat` 的时间差。
- 第三幅：若保存了示例脉冲，会显示局部脉冲波形和到达时刻。

Probe 图：

- 第一幅：相邻叶片到达时间间隔，正常应随转速变化平滑。
- 第二幅：新检测结果和已有 `jilublade_probe*.mat` 的时间差；若没有旧文件，则显示示例脉冲。

## 当前注意

这个脚本中的 probe 到达时刻使用“阈值以上面积一半”的简单定义。原 `maichongcompute.m` 对部分电容通道还做了多项式/LM 拟合，因此两者不一定完全一致。该脚本首先用于检查通道、阈值和脉冲分段是否正确；正式替代旧结果前，需要逐通道比较时间差和波形。
