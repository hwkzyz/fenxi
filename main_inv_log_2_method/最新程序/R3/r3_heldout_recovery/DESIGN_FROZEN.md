# COMSOL-based R3 analysis

This analysis uses nine finite-element clearance states to build one calibrated
response family and five disjoint raw states to generate held-out operating
truth. The matched 0.8 mm state is a calibration-node sanity control. The main
Monte Carlo is conditional on one predeclared noisy low-speed calibration.

Run order:

1. `Run_00_AuditGapLibrary`
2. `Run_01_StateMatchedGate`
3. `Run_02_FreezeReferenceCalibration`
4. `Run_03_SmokeTest`
5. `Run_04_CalibrateDetectionThreshold`
6. `Run_05_Main15dB`

The smoke test is engineering evidence only. The formal main script refuses
to start unless the method-specific null thresholds contain 500 independent
high-speed realizations. Raw held-out curves are used only by the forward
generator; fixed, adaptive and calibrated state-matched fits all use the same
calibration-only family.
