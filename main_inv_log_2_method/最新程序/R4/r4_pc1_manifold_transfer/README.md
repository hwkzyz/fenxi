# Tilt calibration-transfer analysis

The R4 synthetic route is frozen as the method-selection and controlled-validation
route for the unified R4/R5 PC1 coefficient-manifold method. The current PC1
simulation chain is:

1. `Main_01_BuildSimulationSurface`
2. `Main_02_PC1SurfaceValidation`
3. `Main_03_OpenSetTransfer`
4. `Main_04_BoundaryIdentifiability` (boundary diagnostic)
5. `Main_05_PC1Figures` (figures after analyses)
6. `Main_06_NoiselessDynamicGate` (latest panel-f dynamic consequence data)
7. `Main_07_SummarizeDynamicGate` (dynamic CSV aggregation and gate summary)

The `Main_*` files are complete standalone programs. They do not call the old
`Analyze_*` files. Historical `Analyze_*`, `Plot_*`, `Run_*`, and `Summarize_*`
files are kept under `旧程序_20260905`.

`Main_06_NoiselessDynamicGate` and `Main_07_SummarizeDynamicGate` now run the
blind-EO dynamic gate over the C1/C2/FE-oracle geometry models. C1/C2/FE-oracle
remain geometry-transfer ablations; they are not EO priors. Run them after the
static chain when regenerating the data consumed by
`paper_figures/fig04/Plot_R4_Fig4_GeometryMismatch_v12.m`.

The older `Run_01_StaticCalibrationTransferAudit`, `Run_02_NoiselessDynamicGate`,
and `Summarize_TiltCalibrationTransfer` files remain archived historical copies.

The corresponding R5 experimental order and the boundary between production
programs and audit scripts are documented in
`experimental_dynamic_decoupling_validation/20251222_gap_vibration_identification_v1/FROZEN_METHOD_R4R5_20260905.md`.

For the complete local simulation procedure, inputs, outputs, and file roles,
see `程序说明_仿真冻结方法_20260905.md`.

This directory contains the frozen R4 simulation chain without modifying the
Route 30/31 solver. Historical audit files are archived under `旧程序_20260905`.
