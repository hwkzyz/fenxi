# Inv Log 2 Main

## Layout

- `01_trust_domain_calibration`
- `02_fixed_trust_domain_pipeline`
- `trust_domain_result`

## Run Order

```matlab
cd('E:\0小论文+程序\0博士期间小论文+程序\间隙和振动解耦\GapVib_Ortho_Verify\04_inv_log_2_main')
run('01_trust_domain_calibration\Step00_Calibrate_TrustDomain.m')
run('01_trust_domain_calibration\Fig01_TrustDomain_StaticLibrary.m')
run('01_trust_domain_calibration\Fig02_TrustDomain_WindowScan.m')
run('01_trust_domain_calibration\Fig03_TrustDomain_ReconstructionCompare.m')
run('02_fixed_trust_domain_pipeline\Step01_KeyCases_FixedTrust.m')
run('02_fixed_trust_domain_pipeline\Step02_ClosedLoop_FixedTrust.m')
run('02_fixed_trust_domain_pipeline\Step03_Displacement_FixedTrust.m')
run('02_fixed_trust_domain_pipeline\Fig01_KeyCaseRobustness.m')
run('02_fixed_trust_domain_pipeline\Fig02_ClosedLoopTrends.m')
run('02_fixed_trust_domain_pipeline\Fig03_DisplacementVisualization.m')
run('02_fixed_trust_domain_pipeline\Fig04_HighSpeedWindowing.m')
run('02_fixed_trust_domain_pipeline\Fig05_GapCostCurves.m')
run('02_fixed_trust_domain_pipeline\Fig06_WaveformFit.m')
run('02_fixed_trust_domain_pipeline\Fig07_PerformanceCompare.m')
```
