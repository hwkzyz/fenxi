# Trust-Domain Calibration

Run this folder when the static gap library or the response model changes.
It scans candidate spatial windows, performs leave-one-gap validation, and
saves the selected window for the fixed-trust pipeline.

```matlab
Step00_Calibrate_TrustDomain
Plot01_TrustDomain_StaticLibrary
Plot02_TrustDomain_WindowScan
Plot03_TrustDomain_ReconstructionCompare
```

The saved result is written to:

```text
../trust_domain_result/inv_log_2_trust_domain.mat
../trust_domain_result/inv_log_2_trust_domain_summary.csv
../trust_domain_result/inv_log_2_trust_domain_metrics.csv
```
