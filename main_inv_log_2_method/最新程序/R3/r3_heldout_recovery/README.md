# R3: held-out clearance-state recovery

This folder contains the COMSOL-based R3 analysis. It is independent of the
older Step26 hybrid experimental-model study, which remains unchanged.

## Scientific roles

- Nine calibration states build one response family:
  `0.2, 0.3, 0.4, 0.6, 0.8, 1.1, 1.3, 1.4, 1.5 mm`.
- Five disjoint raw COMSOL states generate held-out operating truth:
  `0.5, 0.7, 0.9, 1.0, 1.2 mm`.
- The raw held-out curves never enter fixed, adaptive, or state-matched
  inversion. All three use the same calibration-only family.
- The `0.8 mm` calibration node is an additional matched control.
- The main analysis is conditional on one frozen noisy low-speed calibration.

## Run order

R3 论文的统一盲辨识主入口在上一级目录：

```matlab
Run_R3_BlindSingleFrequency
```

它调用本目录的 `Run_R3_UnifiedBlindBackend`，输出写入
`output/unified_blind_backend`。下面的 `Run_00`--`Run_05` 是 R3 的冻结标定、
阈值和正式 Monte Carlo 主分析链，回答的是留出状态成功率问题；两者输出
不能混为同一张表。

1. `Run_00_AuditGapLibrary`
2. `Run_01_StateMatchedGate`
3. `Run_02_FreezeReferenceCalibration`
4. `Run_03_SmokeTest`
5. `Run_04_CalibrateDetectionThreshold`
6. `Run_05_Main15dB`
7. `Plot_R3_Fig3`
8. `Plot_R3_Supplementary`

`Run_05_Main15dB` refuses to start until the method-specific null thresholds
contain all 500 predeclared realizations. Nonzero cells start with 25 paired
realizations and expand to 75 only when the 95% Wilson interval crosses the
predeclared `P_vib = 0.90` criterion. Null cells use 75 realizations.

The formal plotting entry points are `Plot_R3_Fig3_PaperFinal.m` and
`Plot_R3_Supplementary_Paper.m`. Earlier plot variants and diagnostic runners
are preserved under `归档_历史诊断与出图` and are not part of the frozen route.

## Interpretation

`P_vib` is vibration-recovery success for `A >= 0.10 mm`. `P_FP` is a separate
null statistic for `A = 0`. `P_joint` additionally requires the gap estimate to
meet its tolerance and must not replace `P_vib`. The state-matched result is a
calibrated diagnostic benchmark, not an oracle using held-out raw responses.

The smoke test is engineering evidence only. Formal claims and figures must be
drawn from `output/05_main_15db/main_result.mat` after the complete threshold
calibration and Monte Carlo run.

The current blind-backend implementation is `Run_R3_UnifiedBlindBackend.m`.
Historical data-building and aggregation scripts are under
`旧程序_20260916` and are not part of the default run chain.
